'''
filters.py: display filters (applied to the plot and the quality
analysis only, never to the LSL stream).
'''

import numpy as np
from scipy import signal

from . import runtime as rt
from .config import disp_hp_hz, disp_lp_hz, line_freq


def build_display_filters():
    fs = rt.fs
    sos = [signal.butter(2, disp_hp_hz, btype="highpass", fs=fs, output="sos")]
    nyq = fs / 2.0
    if disp_lp_hz < nyq - 1:
        sos.append(signal.butter(4, disp_lp_hz, btype="lowpass",
                                 fs=fs, output="sos"))
    for h in (1, 2): # mains and first harmonic
        f0 = line_freq * h
        if f0 < nyq - 2:
            b, a = signal.iirnotch(f0, Q=30.0, fs=fs)
            sos.append(signal.tf2sos(b, a))
    return np.vstack(sos)


disp_sos = build_display_filters()


def disp_filter_block(block_eeg):
    '''
    Causal filter for the live display, with state carried between blocks.

    Re-filtering the 5-s window every frame makes the filter ring at both
    ends of the window, with a DC-coupled amplifier and hundreds of mV of drifting offset, those
    edge transients swamp the plot. A filter that has been running
    continuously has no window edges at all. It is initialised to the
    steady state of the first sample, so there is no start-up step
    either. (Causal means a few ms of phase delay, irrelevant for viewing)
    '''
    disp_zi = rt.disp_zi # carried state lives in runtime (reset by reconfigure)
    x = np.asarray(block_eeg, dtype=np.float64)
    if disp_zi["z"] is None or disp_zi["z"].shape[-1] != x.shape[1]:
        zi0 = signal.sosfilt_zi(disp_sos)                # (sections, 2)
        disp_zi["z"] = zi0[:, :, None] * x[0][None, None, :]
    y, disp_zi["z"] = signal.sosfilt(disp_sos, x, axis=0, zi=disp_zi["z"])
    return y.astype(np.float32)


def display_filter(x):
    '''
    Zero-phase filter for analysis of a finished window (quality scores,
    correlation). Linear detrend first: removing the offset and the bulk
    of the drift before filtering is what keeps the window edges from
    ringing into the RMS and correlation estimates.
    '''
    x = signal.detrend(np.asarray(x, dtype=np.float64), axis=0, type="linear")
    return signal.sosfiltfilt(disp_sos, x, axis=0)
