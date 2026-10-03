'''
config.py: every user-tunable constant in one place.

Electrode labels, network settings, display settings and all the
quality-score thresholds. Nothing here changes at runtime.
'''


# Electrode Labels
'''
The 64 labels below are transcribed from the montage diagram, read
top-to-bottom / left-to-right. All 64 resolve against MNE's
standard_1005, so the SET is confirmed.

The ORDER below was read out of an OTBioLab+ .otb+ recording of this
cap, so it is the device's real pin order rather than an inference.
'''

CHANNEL_LABELS = [
    # Verified against an OTBioLab+ recording of this cap (ADEEGCap64SE /
    # ELSCH064EEG).
    "Fp1", "Fpz", "Fp2", "AF7", "AF3", "AFz", "AF4", "AF8",
    "F9", "F7", "F3", "F1", "Fz", "F2", "F4", "F8",
    "F10", "FT9", "FT7", "FC5", "FC3", "FC1", "FCz", "FC2",
    "FC4", "FC6", "FT8", "FT10", "T9", "T7", "C5", "C3",
    "C1", "Cz", "C2", "C4", "C6", "T8", "T10", "TP9",
    "TP7", "CP5", "CP3", "CP1", "CPz", "CP2", "CP4", "CP6",
    "TP8", "TP10", "P9", "P7", "P3", "P1", "Pz", "P2",
    "P4", "P8", "P10", "PO7", "POz", "PO8", "O1", "O2",
]

# Configuration

host, port = "0.0.0.0", 45454
window_sec = 5

# 16 or 24. Same 286 nV per count either way; the extra byte extends
# the range (+/- 2.4 V vs +/- 9.4 mV), not the resolution. 24-bit is the
# safe default, at 16-bit a few mV of electrode offset can clip.
RESOLUTION_BITS = 24

line_freq = 60.0 # confirmed from a real recording: 60 Hz
disp_hp_hz = 1.0 # display high-pass; device is DC coupled
disp_lp_hz = 45.0 # display low-pass

channels_per_page = 16 # how many lanes on screen at once
samples_per_read = 16 # samples pulled from the socket per recv block

# The last 2 channels are accessory (buffer + free-running sample
# counter), not biopotential. They are Unsigned, sign-extending them
# corrupts the counter halfway through its range and would also mangle
# the SyncMini flag on bit 15 in 16-bit mode.
n_accessory = 2

# The trailing channels are OTBioLab+'s AdapterControl pair
# inputs. RampChannel is the free-running sample counter the dropout monitor relies on.
ACCESSORY_NAMES = ["BufferChannel", "RampChannel"]

'''
The Ramp channel is a 16-bit counter even in 24-bit mode: a real
recording shows it climbing to 65535 and wrapping to 0. Computing the
wrap modulo 2^24 instead would report a phantom "dropout" of 16.7 M
samples every 131 s.
'''
RAMP_MOD = 1 << 16

'''
Quality: each channel is first judged on its own signal, with no
reference to the other channels, so a cap with only two gelled
electrodes shows those two as good instead of as outliers of a bad
majority. Cross-channel checks come second, and only when there are
enough well-contacted channels to compare against.
'''
'''
Amplitude window for the robust sd (1.4826 x MAD) of the 1-45 Hz signal.
Robust, so a gelled frontal channel's blinks don't count against it.
Real good contact on this cap measured 15-22 uV; poor contact 100-600.
'''
good_rms_uv = (3.0, 50.0) # full marks inside this range
amp_zero_uv = (1.0, 150.0) # zero at or beyond these (linear below, log above)
eeg_band = (0.5, 40.0) # physiological band of interest, Hz
'''
Absolute mains amplitude (uV at the mains frequency). The ratio-to-EEG
test alone lets a channel whose junk is also large look fine. Measured
on this cap: good contact 5-20 uV, poor contact 1,000-15,000 uV.
'''
mains_abs_ok = 30.0 # at or below -> 100
mains_abs_bad = 300.0 # at or above -> 0 (log scale between)

# Identical-copy check
copy_dc_tol_uv = 2000.0 # offsets within 2 mV
copy_rms_uv = 2.5 # RMS difference (above 1 Hz) below this
copy_min_channels = 3 # group size that counts

# Reference / Ground check
ref_floor_uv = 75.0 # lowest mains on the cap must be at least this
ref_group_span = 3.0 # channels within this x the floor form the group
ref_similar = 2.5 # 90th/10th percentile mains ratio within the group
ref_min_channels = 3 # group size needed
ref_coherence = 0.9 # median correlation of their mains waveforms
ref_history = 3 # banner follows the majority of the last 3 checks

# A channel whose own-signal score reaches this is "well contacted" and
# may be used as a reference point for its neighbours.
trust_min = 60.0
'''
Cross-channel checks (neighbour correlation, head check) need at least
this many well-contacted channels; below it they are reported as "not
assessed" instead of being allowed to fail good channels.
'''
head_min_trusted = 8

# Neighbours are only sought among this many nearest electrodes, so a
# channel is never compared with a far-away one just because it's good.
corr_reach = 12

# Neighbour-correlation sub-score.

'''
Volume conduction makes nearby electrodes see overlapping activity, so
a well-attached channel tracks its neighbours. A bad one does not.
This is the one check that catches artifact sitting inside the EEG
band, which the ratio-based line/hf scores cannot see.

Per channel, r = the 75th percentile of its correlations with the
corr_k nearest electrodes. The score is then how far r falls below the
other channels. Two absolutes still apply: r < 0 scores zero (an inverted or mis-referenced
channel), and if the whole cap is decorrelated every score is scaled
down, so uniformly bad cannot read as uniformly fine.
'''

line_ratio_ok = 1.0 # mains <= EEG amplitude -> 100
line_ratio_bad = 10.0 # mains >= 10x EEG amplitude -> 0

# High-frequency sub-score: power 40-100 Hz (mains excluded) as a
# fraction of 0.5-40 Hz power.
hf_ratio_ok = 0.35 # at or below -> 100
hf_ratio_bad = 1.5 # at or above -> 0

corr_k = 8 # neighbours compared against
corr_r_bad = 0.20 # group median r here or below -> scores -> 0
corr_r_good = 0.45 # group median r here or above -> no scaling

'''
Cap-level "is there a head?" check, on the same average-referenced
neighbour correlation. Measured: a real head gave +0.62 to +0.92 in
every 5-s window, simulated empty caps gave +0.03 to +0.18. The
decision uses the median of the last few evaluations (10 s), because
a single window can dip.
'''

head_r_min = 0.35
head_history = 5
corr_abs_bad = 0.30 # own neighbour r at or below -> 0
corr_abs_ok = 0.55 # (if also an outlier), at or above -> always fine
corr_z_bad = -4.0 # this many MADs below the group -> 0
corr_z_ok = -2.0 # this far below the group is still fine
