'''
quality.py: per-channel quality, judged on each channel's OWN signal.

Sub-scores (0-100, overall = weakest link):
  line: mains, as a ratio to the EEG and as an absolute amplitude
  amp:  robust sd of the 1-45 Hz signal (dead / railed / plausible)
  hf:   broadband power above the EEG band, mains excluded
  corr: optional neighbour-correlation score from cap_checks
'''

import numpy as np
from scipy import signal

from . import runtime as rt
from .config import (eeg_band, line_freq, line_ratio_ok, line_ratio_bad,
                     mains_abs_ok, mains_abs_bad, good_rms_uv, amp_zero_uv,
                     hf_ratio_ok, hf_ratio_bad)
from .filters import display_filter

try:
    _trapz = np.trapezoid # numpy >= 2.0
except AttributeError:
    _trapz = np.trapz # numpy < 2.0


def _bandpower(f, pxx, band):
    mask = (f >= band[0]) & (f <= band[1])
    if not mask.any():
        return 0.0
    return _trapz(pxx[mask], f[mask])


def channel_quality(x, fs_, corr_score=None):
    '''x: one channel's raw window in microvolts. -> (overall_pct, parts).'''
    x = np.asarray(x, dtype=np.float64)
    x = x - np.mean(x)

    '''
    Linear detrend per segment: the front end is DC-coupled, and drift
    leaking into the lowest bins would otherwise inflate "in-band"
    power and flatter a drifting channel.
    '''
    f, pxx = signal.welch(x, fs=fs_, nperseg=min(len(x), fs_),
                          detrend="linear")
    total = _bandpower(f, pxx, eeg_band)

    # Mains, scored on a LOG scale of amplitude ratio (mains / EEG).
    line_pow = _bandpower(f, pxx, (line_freq - 2, line_freq + 2))
    line_amp_ratio = float(np.sqrt(line_pow / (total + 1e-12)))
    if line_amp_ratio <= line_ratio_ok:
        ratio_score = 100.0
    else:
        ratio_score = 100.0 * float(np.clip(
            1 - np.log10(line_amp_ratio / line_ratio_ok)
            / np.log10(line_ratio_bad / line_ratio_ok), 0, 1))
    # on the absolute mains amplitude, so a channel whose broadband junk is as large as its mains cannot hide behind a small ratio.
    mains_uv = float(np.sqrt(line_pow))
    if mains_uv <= mains_abs_ok:
        abs_score = 100.0
    else:
        abs_score = 100.0 * float(np.clip(
            1 - np.log10(mains_uv / mains_abs_ok)
            / np.log10(mains_abs_bad / mains_abs_ok), 0, 1))
    line_score = min(ratio_score, abs_score)

    # amplitude: dead, railed, or plausible. Measured in-band, so DC drift from the DC-coupled front end doesn't dominate.
    xf = display_filter(x[:, None])[:, 0]
    rms = float(np.sqrt(np.mean(xf ** 2)))
    # robust sd: sparse large events (blinks on a well-gelled frontal channel) barely move it; continuous junk (a dry electrode) does.
    rsd = 1.4826 * float(np.median(np.abs(xf - np.median(xf))))
    peak = float(np.max(np.abs(x)))
    lo, hi = good_rms_uv
    zlo, zhi = amp_zero_uv
    if peak >= rt.rail_uv:
        amp_score = 0.0
    elif rsd < lo:
        amp_score = 100.0 * float(np.clip((rsd - zlo) / (lo - zlo), 0, 1))
    elif rsd > hi:
        amp_score = 100.0 * float(np.clip(
            1 - np.log10(rsd / hi) / np.log10(zhi / hi), 0, 1))
    else:
        amp_score = 100.0

    # Broadband noise above the EEG band, with the mains line and its harmonic cut out. 60 Hz sits inside 40-100 Hz.
    hf_hi = min(fs_ / 2 - 1, 100)
    hf = _bandpower(f, pxx, (eeg_band[1], hf_hi))
    # +/-5 Hz, not +/-3: at this 1 Hz resolution the mains skirt spills well past +/-3 Hz on a real recording.
    for h in (1, 2):
        f0 = line_freq * h
        if eeg_band[1] < f0 < hf_hi:
            hf -= _bandpower(f, pxx, (f0 - 5, f0 + 5))
    hf = max(hf, 0.0)
    hf_ratio = hf / (total + 1e-12)
    # in-band power above 40 Hz (normal muscle tone + noise floor).
    hf_score = 100.0 * float(np.clip(
        (hf_ratio_bad - hf_ratio) / (hf_ratio_bad - hf_ratio_ok), 0, 1))

    parts = {"line": line_score, "amp": amp_score, "hf": hf_score,
             "rms": rms, "rsd": rsd, "mains_uv": mains_uv}
    subs = [line_score, amp_score, hf_score]
    # corr is optional: None or NaN means "not assessed" (too few well-contacted neighbours), which must not count as a failure.
    if corr_score is not None and np.isfinite(corr_score):
        parts["corr"] = float(corr_score)
        subs.append(float(corr_score))
    overall = min(subs) # weakest link
    return overall, parts


def own_signal_scores(eeg_window, fs_=None):
    '''Per-channel (overall, parts) from each channel's OWN signal only.'''
    fs_ = rt.fs if fs_ is None else fs_
    res = [channel_quality(eeg_window[:, i], fs_, None)
           for i in range(eeg_window.shape[1])]
    return np.array([r_[0] for r_ in res]), [r_[1] for r_ in res]


def quality_color(pct):
    if pct >= 70:
        return "tab:green"
    if pct >= 40:
        return "tab:orange"
    return "tab:red"
