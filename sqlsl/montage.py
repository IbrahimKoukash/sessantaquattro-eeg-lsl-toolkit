'''
montage.py: electrode coordinates and the nearest-neighbour map.
'''

import numpy as np

from . import runtime as rt
from .config import CHANNEL_LABELS, corr_k


def electrode_positions(labels):
    '''
    Look up 3-D coordinates by label so downstream tools get a montage
    for free. Silently returns {} if mne isn't installed.
    '''
    try:
        import mne
        pos = None
        for nm in ("colin27_1005", "standard_1005"):
            try:
                pos = mne.channels.make_standard_montage(
                    nm).get_positions()["ch_pos"]
                break
            except Exception:
                continue
        if pos is None:
            return {}
    except Exception:
        return {}
    lowered = {k.lower(): k for k in pos}
    out = {}
    for lb in labels:
        key = lowered.get(lb.lower())
        if key is not None:
            out[lb] = tuple(float(c) for c in pos[key])
    return out


# Neighbour correlation

_neighbours = None # (n_eeg, k) indices of nearest electrodes
_neighbour_order = None # (n_eeg, n_eeg-1) all others, nearest first


def build_neighbours():
    '''Nearest electrodes by 3-D position. Falls back to all-channels.'''
    global _neighbours, _neighbour_order
    _neighbours = None
    _neighbour_order = None
    if CHANNEL_LABELS is None:
        return
    n_eeg = rt.n_eeg
    locs = electrode_positions(rt.eeg_labels)
    if len(locs) < n_eeg:
        return
    p = np.array([locs[l] for l in rt.eeg_labels], dtype=float)
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    order = np.argsort(d, axis=1)[:, :n_eeg - 1] # self sorts last (inf)
    _neighbour_order = order
    _neighbours = order[:, :min(corr_k, n_eeg - 1)]
