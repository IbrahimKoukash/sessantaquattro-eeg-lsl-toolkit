'''
cap_checks.py: cross-channel checks and the full cap assessment.

  neighbour_scores: does each channel track its nearest neighbours?
                    Also decides the cap-level "is there a head?" state.
  identical_groups: channels that are copies of each other (not on the
                    scalp, or bridged).
  reference_check:  mains reaching every channel through the shared
                    reference / ground electrode.
  assess_cap:       all of the above combined with the own-signal scores.
'''

import numpy as np
from scipy import signal

from . import runtime as rt
from . import montage
from .config import (trust_min, head_min_trusted, corr_reach, corr_k,
                     head_r_min, corr_z_bad, corr_z_ok, corr_abs_bad,
                     corr_abs_ok, corr_r_bad, corr_r_good,
                     copy_dc_tol_uv, copy_rms_uv, copy_min_channels,
                     amp_zero_uv, ref_floor_uv, ref_group_span, ref_similar,
                     ref_min_channels, ref_coherence, line_freq)
from .filters import display_filter
from .quality import own_signal_scores


# Cap-level result of the last neighbour_scores() call.
'''
state "ok": enough well-contacted channels, head structure present
state "nohead": enough well-contacted channels, but no spatial
                EEG structure among them (empty cap, no reference)
state "insufficient": too few well-contacted channels to judge
'''
head_stat = {"r": None, "n_trusted": 0, "state": "insufficient"}


def neighbour_scores(eeg_window, trusted=None):
    '''
    eeg_window: (samples, n_eeg) raw microvolts.
    trusted: boolean mask of well-contacted channels, derived from each
             channel's own-signal score if not given.
    Returns (scores 0-100, r per channel). Entries are NaN where a channel
    has too few well-contacted neighbours to be judged, "not assessed",
    never "failed".

    Only well-contacted channels are used as references.

    Correlations are computed after re-referencing to the average of the
    well-contacted channels, which removes whatever the single common
    reference electrode puts on every channel and leaves only genuine
    spatial structure: on a real head neighbours agree (r = +0.54 to +0.79
    measured); on an empty cap nothing is left but noise (r = 0).
    '''
    n = eeg_window.shape[1]
    nan = np.full(n, np.nan)
    if eeg_window.shape[0] < 50 or n < 3:
        head_stat.update(r=None, n_trusted=0, state="insufficient")
        return nan, nan.copy()
    if trusted is None:
        trusted = own_signal_scores(eeg_window)[0] >= trust_min
    trusted = np.asarray(trusted, dtype=bool)
    k_tr = int(trusted.sum())
    if k_tr < head_min_trusted:
        head_stat.update(r=None, n_trusted=k_tr, state="insufficient")
        return nan, nan.copy()

    x = display_filter(eeg_window)
    x = x - x[:, trusted].mean(axis=1, keepdims=True) # average reference
    sd = x.std(axis=0)
    live = sd > 1e-9
    c = np.full((n, n), np.nan)
    if live.sum() >= 2:
        idx = np.flatnonzero(live)
        c[np.ix_(idx, idx)] = np.corrcoef(x[:, live].T)

    order = montage._neighbour_order
    r = np.full(n, np.nan)
    for i in range(n):
        if order is not None and order.shape[0] == n:
            near = order[i][:corr_reach]
        else:
            near = [j for j in range(n) if j != i]
        cand = [j for j in near if trusted[j]][:corr_k]
        if len(cand) < 3:
            continue # not enough evidence -> NaN
        v = np.nan_to_num(c[i, cand], nan=0.0)
        # Upper quartile, so one poor neighbour cannot condemn a good one.
        r[i] = np.quantile(v, 0.75)

    rt_ = r[trusted & np.isfinite(r)]
    if rt_.size < head_min_trusted // 2:
        head_stat.update(r=None, n_trusted=k_tr, state="insufficient")
        return nan, r
    med = float(np.median(rt_))
    mad = 1.4826 * float(np.median(np.abs(rt_ - med)))
    head_stat.update(r=med, n_trusted=k_tr,
                     state="ok" if med >= head_r_min else "nohead")
    z = (r - med) / max(mad, 0.03)

    # Condemned only if both an outlier among its peers and poorly
    # correlated in absolute terms (good channels on a real head sat at r >= +0.61; planted defects at r <= +0.28)
    relative = 100.0 * np.clip((z - corr_z_bad) / (corr_z_ok - corr_z_bad),
                               0, 1)
    absolute = 100.0 * np.clip(
        (r - corr_abs_bad) / (corr_abs_ok - corr_abs_bad), 0, 1)
    score = np.maximum(relative, absolute)
    # Clearly anti-correlated with its neighbours: inverted/mis-referenced.
    score = np.where(r < -0.3, 0.0, score)
    score = np.where(np.isfinite(r), score, np.nan)
    return score, r


def identical_groups(eeg_window, fs_=None):
    '''
    Groups of >= copy_min_channels channels that are copies of each other:
    offsets within copy_dc_tol_uv and RMS difference above 1 Hz below
    copy_rms_uv. Returns a list of index arrays (empty if none).
    '''
    fs_ = rt.fs if fs_ is None else fs_
    w = np.asarray(eeg_window, dtype=np.float64)
    n = w.shape[1]
    if w.shape[0] < fs_ or n < copy_min_channels:
        return []
    dc = w.mean(axis=0)
    hp = signal.butter(2, 1.0, btype="highpass", fs=fs_, output="sos")
    x = signal.sosfiltfilt(hp, signal.detrend(w, axis=0), axis=0)
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n - 1):
        near = np.flatnonzero(np.abs(dc[i + 1:] - dc[i]) < copy_dc_tol_uv) + i + 1
        if near.size:
            d = np.sqrt(np.mean((x[:, near] - x[:, [i]]) ** 2, axis=0))
            j = near[d < copy_rms_uv]
            adj[i, j] = adj[j, i] = True
    # connected groups
    seen = np.zeros(n, dtype=bool)
    groups = []
    for i in range(n):
        if seen[i] or not adj[i].any():
            continue
        stack, comp = [i], []
        seen[i] = True
        while stack:
            k = stack.pop()
            comp.append(k)
            for j in np.flatnonzero(adj[k] & ~seen):
                seen[j] = True
                stack.append(j)
        if len(comp) >= copy_min_channels:
            groups.append(np.array(sorted(comp)))
    taken = set(int(i) for g in groups for i in g)
    out = []
    for g in groups:
        g = list(g)
        level = float(np.median(dc[g]))
        g += [int(j) for j in range(n)
              if j not in taken and abs(dc[j] - level) < copy_dc_tol_uv]
        taken.update(g)
        out.append(np.array(sorted(g)))
    return out


def reference_check(eeg_window, parts, fs_=None):
    '''
    Is mains reaching every channel through the shared reference/ground?

    Returns {"state": "ok"|"ref", "channels": [indices], "mains_uv",
    "floor_uv", "coherence"}. "ref" only when (1) no live channel on the
    cap has mains below ref_floor_uv, one good channel proves the
    reference is fine, since its mains would be on every channel, and
    (2) at least ref_min_channels EEG-like channels near that floor carry
    mains of similar size with near-identical waveforms.
    '''
    fs_ = rt.fs if fs_ is None else fs_
    out = {"state": "ok", "channels": [], "mains_uv": None,
           "floor_uv": None, "coherence": None}
    mains = np.array([p["mains_uv"] for p in parts])
    rsd = np.array([p["rsd"] for p in parts])
    amp = np.array([p["amp"] for p in parts])
    live = (rsd >= amp_zero_uv[0]) & ~np.array([bool(p.get("copy"))
                                                 for p in parts])
    if not live.any() or eeg_window.shape[0] < fs_:
        return out
    floor = float(mains[live].min())
    out["floor_uv"] = floor
    if floor < ref_floor_uv:
        return out
    # EEG-like on its own signal: plausible amplitude (the display notch has already taken the mains out of that measure)
    grp = np.flatnonzero(live & (amp >= 60) & (mains <= ref_group_span * floor))
    if grp.size < ref_min_channels:
        return out
    m = mains[grp]
    if np.percentile(m, 90) / max(np.percentile(m, 10), 1e-9) > ref_similar:
        return out
    hi = min(line_freq + 1.0, fs_ / 2 - 1)
    sos_ = signal.butter(2, [line_freq - 1.0, hi], btype="bandpass",
                         fs=fs_, output="sos")
    seg = eeg_window[:, grp].astype(np.float64)
    nb = signal.sosfiltfilt(sos_, seg - seg.mean(axis=0), axis=0)
    c = np.corrcoef(nb.T)
    np.fill_diagonal(c, np.nan)
    # members: channels whose mains waveform matches the rest of the group
    each = np.nanmedian(c, axis=1)
    keep = each >= ref_coherence
    if keep.sum() < ref_min_channels:
        return out
    sub = c[np.ix_(keep, keep)]
    coh = float(np.nanmedian(sub))
    out.update(channels=[int(i) for i in grp[keep]],
               mains_uv=float(np.median(m[keep])), coherence=coh, state="ref")
    return out


def assess_cap(eeg_window, fs_=None):
    '''
    Full quality assessment of one window. Returns (overall[n], parts[n], r[n], cap) where cap is a copy of
    head_stat plus cap["ref"], the result of reference_check().
    '''
    fs_ = rt.fs if fs_ is None else fs_
    own, parts = own_signal_scores(eeg_window, fs_)
    copies = identical_groups(eeg_window, fs_)
    is_copy = np.zeros(len(own), dtype=bool)
    for g in copies:
        is_copy[g] = True
    for i in np.flatnonzero(is_copy):
        parts[i]["copy"] = True
    own = np.where(is_copy, 0.0, own)
    trusted = own >= trust_min
    cs, rr = neighbour_scores(eeg_window, trusted)
    overall = own.copy()
    for i in range(len(own)):
        if np.isfinite(cs[i]):
            parts[i]["corr"] = float(cs[i])
            overall[i] = min(overall[i], cs[i])

    '''
    No head under an otherwise clean-looking cap: every channel "looks
    fine" on its own, so scale down every channel that has not itself
    shown real neighbour agreement. A channel that has (r >= corr_abs_ok)
    has proved there is structure around it and keeps its score.
    '''
    if head_stat["state"] == "nohead":
        gate = float(np.clip((head_stat["r"] - corr_r_bad)
                             / (corr_r_good - corr_r_bad), 0, 1))
        proven = np.isfinite(rr) & (rr >= corr_abs_ok)
        overall = np.where(proven, overall, overall * gate)
    overall = np.where(is_copy, 0.0, overall)
    rr = np.where(is_copy, np.nan, rr)
    cap = dict(head_stat)
    cap["copies"] = [[int(i) for i in g] for g in copies]
    cap["ref"] = reference_check(eeg_window, parts, fs_)
    return overall, parts, rr, cap
