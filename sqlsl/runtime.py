'''
runtime.py: shared, mutable stream state.

Channel layout (nch, n_eeg, labels), the ring buffers the acquisition
thread writes and the viewer reads, the lock, and the stop flag.

Other modules must read these as attributes (runtime.nch,
runtime.eeg_labels, ...) rather than "from runtime import nch":
reconfigure() rebinds them when the device turns out to send a
different channel count, and a copied name would go stale.
'''

import threading
import numpy as np

from .config import (CHANNEL_LABELS, ACCESSORY_NAMES, n_accessory,
                     window_sec)
from .protocol import create_bin_command, lsb_microvolts


# Shared state
command, nch, fs, bps = create_bin_command()


def eeg_count(total):
    '''How many of total channels are cap electrodes.'''
    if CHANNEL_LABELS is not None and len(CHANNEL_LABELS) <= total - n_accessory:
        return len(CHANNEL_LABELS)
    return total - n_accessory


def trailing_labels(total, n_e):
    '''Names for the non-EEG channels: any AUX inputs, then Buffer/Ramp.'''
    extra = total - n_e
    names = (["AUX{}".format(i + 1) for i in range(max(0, extra - 2))]
             + ACCESSORY_NAMES)
    return names[-extra:] if extra else []


n_eeg = eeg_count(nch)
win = window_sec * fs
to_uv = lsb_microvolts(bps)

# Full-scale, used to detect railing (same value either resolution)
full_scale_uv = (2 ** (8 * bps - 1)) * to_uv
rail_uv = 0.95 * full_scale_uv

if CHANNEL_LABELS is None:
    eeg_labels = ["EEG{}".format(i + 1) for i in range(n_eeg)]
else:
    if len(CHANNEL_LABELS) != n_eeg:
        raise ValueError(
            "CHANNEL_LABELS has {} entries, expected {}".format(
                len(CHANNEL_LABELS), n_eeg))
    eeg_labels = list(CHANNEL_LABELS)

aux_labels = trailing_labels(nch, n_eeg)

ring = np.zeros((win, nch), dtype=np.float32)
ring_w = {"i": 0} # circular write index
ring_lock = threading.Lock()
stop_event = threading.Event()
sample_counter = {"n": 0}


# Display-filtered EEG, written in lock-step with ring (same index).
# Filtered continuously in the acquisition thread with carried state (see filters.disp_filter_block())
disp_ring = np.zeros((win, n_eeg), dtype=np.float32)
disp_zi = {"z": None} # carried filter state for the live display


def _circ_write(buf, block, w):
    m = block.shape[0]
    if m >= win:
        buf[:] = block[-win:]
        return
    end = w + m
    if end <= win:
        buf[w:end] = block
    else:
        k = win - w
        buf[w:] = block[:k]
        buf[:end - win] = block[k:]


def ring_write(block, dblock=None):
    '''block: (m, nch) raw; dblock: (m, n_eeg) filtered. Hold ring_lock.'''
    m = block.shape[0]
    w = ring_w["i"]
    _circ_write(ring, block, w)
    if dblock is not None:
        _circ_write(disp_ring, dblock, w)
    ring_w["i"] = 0 if m >= win else (w + m) % win


def _unroll(buf):
    w = ring_w["i"]
    if w == 0:
        return buf.copy()
    return np.concatenate((buf[w:], buf[:w]), axis=0)


def ring_snapshot():
    '''Oldest-to-newest copy of the raw window. Caller holds ring_lock.'''
    return _unroll(ring)


def disp_snapshot():
    '''Oldest-to-newest copy of the filtered EEG. Caller holds ring_lock.'''
    return _unroll(disp_ring)


def reconfigure(new_nch):
    '''Rebuild the globals that depend on the channel count.'''
    global nch, n_eeg, ring, disp_ring, aux_labels, eeg_labels
    nch = new_nch
    n_eeg = eeg_count(nch)
    ring = np.zeros((win, nch), dtype=np.float32)
    disp_ring = np.zeros((win, n_eeg), dtype=np.float32)
    disp_zi["z"] = None
    ring_w["i"] = 0
    aux_labels = trailing_labels(nch, n_eeg)
    if CHANNEL_LABELS is not None and len(CHANNEL_LABELS) == n_eeg:
        eeg_labels = list(CHANNEL_LABELS)
    else:
        if CHANNEL_LABELS is not None:
            print("WARNING: CHANNEL_LABELS has {} entries but {} EEG "
                  "channels were detected; falling back to EEG1....EEG{}."
                  .format(len(CHANNEL_LABELS), n_eeg, n_eeg))
        eeg_labels = ["EEG{}".format(i + 1) for i in range(n_eeg)]
    print("layout: {} EEG + {}".format(n_eeg, ", ".join(aux_labels)))
