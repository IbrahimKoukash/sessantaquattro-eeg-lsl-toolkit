'''
sq_lsl_viewer.py: Sessantaquattro (OT Bioelettronica, 64 ch) -> LSL,
                     with a live viewer and per-channel quality scores

Acquisition thread: TCP from the device -> layout check -> decode ->
                     LSL outlet (raw, microvolts) + continuous display
                     filter -> ring buffers
Main thread: matplotlib viewer, 16 channels per page, every trace
                     confined to its own lane, quality % per channel

Documentation: README.md
'''

import socket
import threading
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Button
from scipy import signal
from pylsl import StreamInfo, StreamOutlet, local_clock

try:
    _trapz = np.trapezoid # numpy >= 2.0
except AttributeError:
    _trapz = np.trapz # numpy < 2.0

# Left/right arrows are bound by default to the navigation toolbar's
# back/forward history, which can fight the paging controls. Release
# them before any figure is created.
for _km, _keys in (("keymap.back", ["left", "c", "backspace"]),
                   ("keymap.forward", ["right", "v"])):
    try:
        matplotlib.rcParams[_km] = [k for k in matplotlib.rcParams[_km]
                                    if k not in _keys]
    except Exception:
        pass


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


# Protocol helpers

def integer_to_bytes(command):
    return int(command).to_bytes(2, byteorder="big")


def create_bin_command(start=1, imp=0):
    '''
    Build the 2-byte Sessantaquattro command word.

    Field meanings follow OT Bioelettronica's reference script
    (OTB-Matlab):
      start (GO): bit 0 / 1 = send settings and start data transfer
      rec: bit 1 / 1 = record to the device's SD card
      trig: bits 2-3 / 0 = transfer/SD controlled remotely,
                       3 = SD recording from the pushbutton
      ext: (EXTEN) bits 4-5 / INPUT RANGE: 0 standard, 1 x2, 2 x4, 3 x8.
                              Keep 0: the 0.286 uV/count factor assumes it.
      hpf:bit 6  / 0 = DC coupled, 1 = hardware high-pass on
      hres: bit 7 / 0 = 16-bit, 1 = 24-bit samples
      mode: bits 8-10 / 0 monopolar, 1 bipolar, 2 differential,
                        3 accelerometers, 6 impedance check, 7 test
      nch: bits 11-12 / 0/1/2/3 = 8/16/32/64 channels
      fsamp: bits 13-14 /  0/1/2 = 500/1000/2000 Hz (mode != 3)
    imp=1 selects impedance-check mode (6). The viewer does not decode
    impedance data; it exists only so a caller can request it.
    '''
    rec = 0
    trig = 0 # data transfer controlled remotely by this command
    ext = 0 # standard input range (other values change the scaling)
    hpf = 0 # DC coupled
    hres = 1 if RESOLUTION_BITS == 24 else 0
    mode = 0 if imp == 0 else 6 # 0 monopolar; 6 impedance check
    nch = 3 # 64 biopotential inputs
    fsamp = 0 # 500 Hz
    getset = 0
    command = (start + rec * 2 + trig * 4 + ext * 16 + hpf * 64
               + hres * 128 + mode * 256 + nch * 2048
               + fsamp * 8192 + getset * 32768)
    bps = 3 if hres else 2

    # Total channels in the stream: 64 EEG + 2 AUX + Buffer + Ramp = 68.
    # OT Bioelettronica's reference script uses NumChan = 68 for 64-channel
    total = 68
    return integer_to_bytes(command), total, 500, bps


def lsb_microvolts(bps):
    '''
    286 nV per bit, for both resolutions.

    OT Bioelettronica's reference script: it applies ConvFact = 0.000286 mV to 24-bit
    and 16-bit data alike (OTB-Matlab).
    '''
    return 0.286


def disconnect_from_sq(conn):
    if conn is None:
        return
    stop, _, _, _ = create_bin_command(start=0)
    try:
        conn.send(stop)
        conn.shutdown(socket.SHUT_RDWR)
    except OSError:
        pass
    conn.close()

# Shared state
command, nch, fs, bps = create_bin_command()

# The trailing channels are OTBioLab+'s AdapterControl pair
# inputs. RampChannel is the free-running sample counter the dropout monitor relies on.
ACCESSORY_NAMES = ["BufferChannel", "RampChannel"]


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
# Filtered continuously in the acquisition thread with carried state (see disp_filter_block())
disp_ring = np.zeros((win, n_eeg), dtype=np.float32)


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

# Decoding

def read_block(conn, nsamp):
    '''
    Read exactly nsamp samples worth of bytes.

    Returns (buf, t) where t is taken the instant the last byte lands,
    before any decoding. Decoding is numpy work that can be preempted by
    the GUI thread, and any delay between arrival and local_clock()
    inflates the timestamp for the whole chunk.
    '''
    size = nsamp * nch * bps
    buf = bytearray(size)
    view = memoryview(buf)
    got = 0
    while got < size:
        n = conn.recv_into(view[got:], size - got)
        if n == 0:
            raise ConnectionError("Connection closed by device")
        got += n
    return buf, local_clock()


def decode_block(buf, nsamp):
    '''
    Big-endian -> (samples, channels) float32, EEG in microvolts.

    EEG and AUX are two's complement; the trailing accessory channels
    are unsigned. Also returns the raw integer array, so the sample
    counter can be read without microvolt scaling applied.
    '''
    a = np.frombuffer(buf, dtype=np.uint8).reshape(nsamp, nch, bps)
    a = a.astype(np.int64)
    if bps == 3:
        v = (a[:, :, 0] << 16) | (a[:, :, 1] << 8) | a[:, :, 2]
    else:
        v = (a[:, :, 0] << 8) | a[:, :, 1]

    n_signed = nch - n_accessory
    half = 1 << (8 * bps - 1)
    full = 1 << (8 * bps)
    sv = v[:, :n_signed]
    v[:, :n_signed] = np.where(sv >= half, sv - full, sv)

    out = v.astype(np.float32)
    out[:, :n_eeg] *= to_uv # AUX and accessory stay raw
    return out, v


# Channel-count / alignment probe
'''
The number of channels per sample decides how the byte stream is cut
up. Get it wrong and the cut slips by a fixed amount every sample, so
each "channel" is really a rotating mixture of physical inputs, and
with flat or near-identical channels that is completely invisible.

The Ramp channel is a free-running counter, so the correct layout is
the one where some column advances by exactly +1 every sample. That is
a strong signature: a wrong channel count almost never produces one.
'''

def detect_channel_count(conn, candidates=(66, 68, 70, 72), nsamp=48):
    '''Probe the live stream. Returns (nch, ramp_col, leftover_bytes).'''
    need = max(candidates) * bps * (nsamp + 2)
    buf = bytearray()
    while len(buf) < need:
        chunk = conn.recv(need - len(buf))
        if not chunk:
            raise ConnectionError("Connection closed while probing")
        buf += chunk

    raw = np.frombuffer(bytes(buf), dtype=np.uint8).astype(np.int64)
    mod = 1 << (8 * bps)
    hits = []

    for c in candidates:
        per = c * bps
        for phase in range(per): # in case the stream is offset
            usable = (len(raw) - phase) // per
            if usable < 8:
                continue
            a = raw[phase:phase + usable * per].reshape(usable, c, bps)
            if bps == 3:
                v = (a[:, :, 0] << 16) | (a[:, :, 1] << 8) | a[:, :, 2]
            else:
                v = (a[:, :, 0] << 8) | a[:, :, 1]
            # The Ramp is a 16-bit counter whatever the resolution, so a
            # step of +1 modulo 65536 counts as a clean increment.
            d = np.diff(v, axis=0) % RAMP_MOD
            ramp = np.flatnonzero(np.all(d == 1, axis=0))
            if ramp.size:
                hits.append((c, int(ramp[0]), phase))
                break

    if not hits:
        return None, None, buf
    # Prefer the smallest layout that explains the data; a larger one can
    # sometimes alias onto a real ramp by coincidence.
    hits.sort(key=lambda h: h[0])
    c, col, phase = hits[0]
    if len(hits) > 1:
        print("  note: more than one layout fits ({}); taking the smallest."
              .format(", ".join(str(h[0]) for h in hits)))
    '''
    The probe has already taken these bytes off the socket, so they must
    all be handed back, dropping a partial sample here would shift
    every later read and defeat the point of probing at all. Top the
    buffer up to a whole number of samples instead.
    '''
    per = c * bps
    rem = (len(buf) - phase) % per
    if rem:
        want = per - rem
        while want > 0:
            chunk = conn.recv(want)
            if not chunk:
                raise ConnectionError("Connection closed while realigning")
            buf += chunk
            want -= len(chunk)
    return c, col, bytes(buf)[phase:]


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



# Sample-counter / dropout monitor
'''
The device sends a free-running counter on one accessory channel.
It is the ONLY way to detect samples lost on the wireless link:
TCP delivers a continuous byte stream with no gap, so without this
check a dropout silently shifts every later sample earlier in time
relative to your other LSL streams.
'''

counter = {"col": None, "last": None, "lost": 0, "events": 0,
           "checked": False}
'''
The Ramp channel is a 16-bit counter even in 24-bit mode: a real
recording shows it climbing to 65535 and wrapping to 0. Computing the
wrap modulo 2^24 instead would report a phantom "dropout" of 16.7 M
samples every 131 s.
'''
RAMP_MOD = 1 << 16
counter_mod = RAMP_MOD


def find_counter_column(v):
    '''Pick the accessory column that steps by exactly +1 per sample.'''
    if v.shape[0] < 4:
        return None
    for c in range(nch - n_accessory, nch):
        d = np.diff(v[:, c]) % counter_mod
        if np.all(d == 1):
            return c
    return None


def check_counter(v):
    '''Returns number of samples lost immediately before this block.'''
    if not counter["checked"]:
        counter["col"] = find_counter_column(v)
        counter["checked"] = True
        if counter["col"] is None:
            print("Warning: no sample counter found on the accessory "
                  "channels. Dropout detection is OFF.")
        else:
            print("Sample counter on channel index {} "
                  "(dropout detection active).".format(counter["col"]))

    c = counter["col"]
    if c is None:
        return 0

    col = v[:, c]
    lost = 0

    if counter["last"] is not None:
        gap = int((col[0] - counter["last"] - 1) % counter_mod)
        if gap:
            lost += gap

    d = np.diff(col) % counter_mod
    bad = d != 1
    if bad.any():
        lost += int((d[bad] - 1).sum())

    counter["last"] = int(col[-1])
    if lost:
        counter["lost"] += lost
        counter["events"] += 1
        print("DROPOUT: {} sample(s) lost on the link "
              "({:.1f} ms). Total {} in {} event(s)."
              .format(lost, 1000.0 * lost / fs,
                      counter["lost"], counter["events"]))
    return lost


# LSL outlet

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


def make_outlet():
    info = StreamInfo(
        name="Sessantaquattro",
        type="EEG",
        channel_count=nch,
        nominal_srate=fs,
        channel_format="float32",
        source_id="sessantaquattro_64ch")

    locs = electrode_positions(eeg_labels)
    chns = info.desc().append_child("channels")

    for lb in eeg_labels:
        c = chns.append_child("channel")
        c.append_child_value("label", lb)
        c.append_child_value("unit", "microvolts")
        c.append_child_value("type", "EEG")
        if lb in locs:
            x, y, z = locs[lb]
            loc = c.append_child("location")
            loc.append_child_value("X", "{:.6f}".format(x))
            loc.append_child_value("Y", "{:.6f}".format(y))
            loc.append_child_value("Z", "{:.6f}".format(z))

    for lb in aux_labels:
        c = chns.append_child("channel")
        c.append_child_value("label", lb)
        c.append_child_value("unit", "raw")
        # Not AUX inputs -- OTBioLab+ calls these AdapterControl.
        c.append_child_value("type", "Misc")

    acq = info.desc().append_child("acquisition")
    acq.append_child_value("manufacturer", "OT Bioelettronica")
    acq.append_child_value("model", "Sessantaquattro")
    acq.append_child_value("resolution_bits", str(8 * bps))
    acq.append_child_value("lsb_microvolts", "{:.8f}".format(to_uv))
    if not locs:
        print("Note: no electrode coordinates written "
              "(mne not installed, or labels not in standard_1005).")

    return StreamOutlet(info, chunk_size=samples_per_read, max_buffered=360)


# Display filters (applied to the plot only, never to the LSL stream)

def build_display_filters():
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
disp_zi = {"z": None} # carried filter state for the live display


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

# Quality

def _bandpower(f, pxx, band):
    mask = (f >= band[0]) & (f <= band[1])
    if not mask.any():
        return 0.0
    return _trapz(pxx[mask], f[mask])

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
    locs = electrode_positions(eeg_labels)
    if len(locs) < n_eeg:
        return
    p = np.array([locs[l] for l in eeg_labels], dtype=float)
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    order = np.argsort(d, axis=1)[:, :n_eeg - 1] # self sorts last (inf)
    _neighbour_order = order
    _neighbours = order[:, :min(corr_k, n_eeg - 1)]


# Cap-level result of the last neighbour_scores() call.
'''
state "ok": enough well-contacted channels, head structure present
state "nohead": enough well-contacted channels, but no spatial
                EEG structure among them (empty cap, no reference)
state "insufficient": too few well-contacted channels to judge
'''
head_stat = {"r": None, "n_trusted": 0, "state": "insufficient"}


def own_signal_scores(eeg_window, fs_=None):
    '''Per-channel (overall, parts) from each channel's OWN signal only.'''
    fs_ = fs if fs_ is None else fs_
    res = [channel_quality(eeg_window[:, i], fs_, None)
           for i in range(eeg_window.shape[1])]
    return np.array([r_[0] for r_ in res]), [r_[1] for r_ in res]


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

    r = np.full(n, np.nan)
    for i in range(n):
        if _neighbour_order is not None and _neighbour_order.shape[0] == n:
            near = _neighbour_order[i][:corr_reach]
        else:
            near = [j for j in range(n) if j != i]
        cand = [j for j in near if trusted[j]][:corr_k]
        if len(cand) < 3:
            continue # not enough evidence -> NaN
        v = np.nan_to_num(c[i, cand], nan=0.0)
        # Upper quartile, so one poor neighbour cannot condemn a good one.
        r[i] = np.quantile(v, 0.75)

    rt = r[trusted & np.isfinite(r)]
    if rt.size < head_min_trusted // 2:
        head_stat.update(r=None, n_trusted=k_tr, state="insufficient")
        return nan, r
    med = float(np.median(rt))
    mad = 1.4826 * float(np.median(np.abs(rt - med)))
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
    fs_ = fs if fs_ is None else fs_
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
    fs_ = fs if fs_ is None else fs_
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
    fs_ = fs if fs_ is None else fs_
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
    if peak >= rail_uv:
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


def quality_color(pct):
    if pct >= 70:
        return "tab:green"
    if pct >= 40:
        return "tab:orange"
    return "tab:red"


def print_quality(snapshot):
    ov, parts, rr, hs = assess_cap(snapshot[:, :n_eeg])
    print("\n--- EEG quality (%) ---")
    print("{:>6}  {:>6}  {:<20} {:>5} {:>5} {:>5} {:>5} {:>8} {:>8} {:>6}"
          .format("ch", "score", "", "line", "amp", "hf", "corr",
                  "rsd_uV", "60Hz_uV", "r"))
    for ch in range(n_eeg):
        p = parts[ch]
        print("{:>6}  {:5.1f}%  {:<20} {:5.0f} {:5.0f} {:5.0f} {:>5} "
              "{:8.1f} {:8.1f} {:>6}"
              .format(eeg_labels[ch], ov[ch],
                      "COPY" if p.get("copy") else "#" * int(ov[ch] / 5),
                      p["line"], p["amp"], p["hf"],
                      "{:.0f}".format(p["corr"]) if "corr" in p else "-",
                      p["rsd"], p["mains_uv"],
                      "{:+.2f}".format(rr[ch]) if np.isfinite(rr[ch])
                      else "-"))
    good = [eeg_labels[c] for c in range(n_eeg) if ov[c] >= 70]
    bad = [eeg_labels[c] for c in np.argsort(ov) if ov[c] < 40]
    print("\n{} good (>=70%): {}".format(len(good), ", ".join(good) or "none"))
    print("{} bad  (<40%): {}".format(len(bad), ", ".join(bad[:20]) +
                                        (" ..." if len(bad) > 20 else "")))
    print("cap-level: {} ({} well-contacted{})".format(
        hs["state"], hs["n_trusted"],
        "" if hs["r"] is None else ", spatial r {:+.2f}".format(hs["r"])))
    for g in hs.get("copies", []):
        print("IDENTICAL COPIES: {} channels have the same offset and the same "
              "signal (not touching the scalp, or bridged): {}".format(
                  len(g), ", ".join(eeg_labels[c] for c in g)))
    ref = hs.get("ref", {})
    if ref.get("state") == "ref":
        print("REFERENCE/GROUND: {} channels look like EEG but carry the same "
              "~{:.0f} uV of {:.0f} Hz (waveforms {:.2f} alike): {}".format(
                  len(ref["channels"]), ref["mains_uv"], line_freq,
                  ref["coherence"],
                  ", ".join(eeg_labels[c] for c in ref["channels"])))
        print("  No channel on the cap is below {:.0f} uV of mains, so it is "
              "coming in through the shared reference or ground electrode."
              .format(ref["floor_uv"]))
        print("  Check those electrodes first; re-gelling the channels "
              "above will not help.")
    print("  '-' in corr/r = not assessed: fewer than 3 well-contacted "
          "neighbours to compare with.")

# Acquisition thread

def acquisition(conn, outlet, primed=b""):
    try:
        # Samples already pulled off the socket by the layout probe.
        if primed:
            m = len(primed) // (nch * bps)
            if m:
                block, raw = decode_block(primed[:m * nch * bps], m)
                check_counter(raw)
                outlet.push_chunk(block.tolist(), local_clock())
                dblock = disp_filter_block(block[:, :n_eeg])
                with ring_lock:
                    ring_write(block, dblock)
                    sample_counter["n"] += m

        while not stop_event.is_set():
            buf, t_recv = read_block(conn, samples_per_read)
            block, raw = decode_block(buf, samples_per_read)
            check_counter(raw)

            '''
            One timestamp for the chunk = the arrival time of its last
            sample; pylsl back-fills the rest from the nominal rate.
            Residual per-chunk jitter is removed downstream by LSL
            dejittering (see the note in the header).
            '''
            outlet.push_chunk(block.tolist(), t_recv)

            # Display filter runs here, continuously, not per frame.
            dblock = disp_filter_block(block[:, :n_eeg])
            with ring_lock:
                ring_write(block, dblock)
                sample_counter["n"] += block.shape[0]
    except (ConnectionError, OSError) as e:
        print("Acquisition stopped:", e)
    finally:
        stop_event.set()

# Viewer

def main():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind((host, port))
    sock.listen(1)
    print("Waiting for Sessantaquattro connection...")
    conn, addr = sock.accept()
    print("Connected from:", addr)
    conn.send(command)
    print("Start command sent.")
    print("(If no data flows, press the Sessantaquattro pushbutton.)")

    # Verify the stream layout before a single sample is trusted.
    print("Probing stream layout (expecting {} channels)...".format(nch))
    found, ramp_col, primed = detect_channel_count(conn)
    if found is None:
        print("\n  WARNING: could not find the Ramp counter in any of the")
        print("  layouts tried, so the channel count could not be confirmed.")
        print("  Continuing with {} channels, but if the traces look wrong"
              .format(nch))
        print("  this is the first thing to suspect.")
        primed = b""
    else:
        if found != nch:
            print("\n  Device is sending {} channels, not {}. Adjusting."
                  .format(found, nch))
            reconfigure(found)
        else:
            print("  Confirmed: {} channels, Ramp counter on index {}."
                  .format(found, ramp_col))
        counter["col"] = ramp_col
        counter["checked"] = True

    build_neighbours()
    if _neighbours is not None:
        print("  neighbour map: {} nearest electrodes per channel "
              "(e.g. {} -> {})".format(
                  _neighbours.shape[1], eeg_labels[0],
                  ", ".join(eeg_labels[j] for j in _neighbours[0])))
    else:
        print("  neighbour map: no electrode coordinates; correlation "
              "will use all channels instead of nearest neighbours.")

    outlet = make_outlet()
    print("LSL outlet 'Sessantaquattro' created "
          "({} ch, {} EEG, {} Hz, {}-bit, {:.6f} uV/bit)."
          .format(nch, n_eeg, fs, 8 * bps, to_uv))
    print("Streaming to LSL + viewer...")

    th = threading.Thread(target=acquisition,
                          args=(conn, outlet, primed), daemon=True)
    th.start()

    n_pages = int(np.ceil(n_eeg / channels_per_page))
    state = {"page": 0, "gain": 1.0, "auto": True, "step": 50.0,
             "keys_seen": 0, "clip": True}

    fig, ax = plt.subplots(figsize=(11, 8))
    t = np.linspace(-window_sec, 0, win)
    offsets = np.arange(channels_per_page) * 1.0 # rescaled each frame

    lines = []
    qtexts = []
    base_colors = []
    for i in range(channels_per_page):
        (ln,) = ax.plot(t, np.zeros(win), linewidth=0.7)
        lines.append(ln)
        base_colors.append(ln.get_color())
        txt = ax.text(1.005, i, "", transform=ax.get_yaxis_transform(),
                      va="center", ha="left", fontsize=8, clip_on=False)
        qtexts.append(txt)

    ax.set_xlabel("seconds")
    ax.set_xlim(-window_sec, 0)
    ax.set_yticks(offsets)

    # scale bar
    scale_line, = ax.plot([], [], color="k", linewidth=2.5,
                          solid_capstyle="butt")
    scale_text = ax.text(0, 0, "", fontsize=9, va="bottom", ha="left",
                         bbox=dict(facecolor="white", edgecolor="none",
                                   alpha=0.75, pad=1.0))

    def page_channels():
        lo = state["page"] * channels_per_page
        return list(range(lo, min(lo + channels_per_page, n_eeg)))

    def change_page(delta):
        state["page"] = int(np.clip(state["page"] + delta, 0, n_pages - 1))

    def change_gain(factor):
        state["gain"] *= factor

    def dump_quality(_evt=None):
        with ring_lock:
            snap = ring_snapshot()
        print_quality(snap)

    def on_key(event):
        state["keys_seen"] += 1
        k = event.key
        if k in ("down", "pagedown", "right", "n", "j", " "):
            change_page(+1)
        elif k in ("up", "pageup", "left", "p", "k"):
            change_page(-1)
        elif k in ("+", "="):
            change_gain(1.3)
        elif k in ("-", "_"):
            change_gain(1 / 1.3)
        elif k == "a":
            state["auto"] = not state["auto"]
        elif k == "0":
            state["gain"] = 1.0
        elif k == "q":
            dump_quality()
        elif k == "c":
            state["clip"] = not state["clip"]
        elif k == "home":
            state["page"] = 0
        elif k == "end":
            state["page"] = n_pages - 1

    def on_scroll(event):
        change_page(+1 if event.button == "down" else -1)

    fig.canvas.mpl_connect("key_press_event", on_key)
    fig.canvas.mpl_connect("scroll_event", on_scroll)
    
    '''
    On-screen controls. Keyboard focus is unreliable across backends
    and absent entirely in inline/notebook ones, so every action has a
    button too. Buttons must stay referenced or they stop responding.
    '''
    fig.subplots_adjust(left=0.10, right=0.90, top=0.90, bottom=0.14)

    def mkbutton(x, w, label, cb):
        bax = fig.add_axes([x, 0.035, w, 0.055])
        b = Button(bax, label)
        b.on_clicked(cb)
        return b

    buttons = [
        mkbutton(0.10, 0.11, "< Prev", lambda e: change_page(-1)),
        mkbutton(0.22, 0.11, "Next >", lambda e: change_page(+1)),
        mkbutton(0.36, 0.08, "Zoom +", lambda e: change_gain(1.3)),
        mkbutton(0.45, 0.08, "Zoom -", lambda e: change_gain(1 / 1.3)),
        mkbutton(0.56, 0.10, "Auto",
                 lambda e: state.update(auto=not state["auto"])),
        mkbutton(0.68, 0.10, "Reset",
                 lambda e: state.update(gain=1.0)),
        mkbutton(0.80, 0.10, "Quality", dump_quality),
    ]
    state["_buttons"] = buttons # keep them alive

    backend = matplotlib.get_backend()
    print("matplotlib backend: {}".format(backend))
    if backend.lower() in ("agg", "template") or "inline" in backend.lower():
        print("  WARNING: this backend is non-interactive. The window will")
        print("  not accept clicks or key presses. Run the script from a")
        print("  terminal (python3 sq_lsl_viewer.py), or select a GUI")
        print("  backend such as TkAgg / QtAgg before importing pyplot.")
    else:
        print("  Paging: click < Prev / Next >, scroll the wheel, or use the")
        print("  arrow keys after clicking once ON THE PLOT to give it focus.")
        print("  Traces are confined to their own lane; noisy ones turn grey")
        print("  and are marked CLIP. Press 'c' to toggle confinement.")

    last_qual = {"n": -10 ** 9, "vals": np.zeros(n_eeg)}
    head = {"hist": [], "state": None, "n_trusted": 0}
    refst = {"hist": [], "on": False, "info": None}
    copyst = {"chans": set(), "n_shown": 0}
    head_banner = fig.text(0.5, 0.925, "", ha="center", va="center",
                           fontsize=12, fontweight="bold", color="white",
                           bbox=dict(facecolor="#c53030", edgecolor="none",
                                     pad=5), visible=False)

    def update(_frame):
        with ring_lock:
            snap = ring_snapshot()
            dsnap = disp_snapshot()
            total = sample_counter["n"]

        chans = page_channels()

        # quality, every 2 s, for all channels. Computed first, because the display scale and the colours below
        # both depend on which channels are well contacted.
        if total - last_qual["n"] >= 2 * fs and total >= fs:
            ov, _, _, hs = assess_cap(snap[:, :n_eeg])
            last_qual["vals"][:] = ov
            last_qual["n"] = total
            if hs["state"] in ("ok", "nohead"):
                head["hist"] = (head["hist"] + [hs["r"]])[-head_history:]
                st_ = ("ok" if float(np.median(head["hist"])) >= head_r_min
                       else "nohead")
            else:
                head["hist"] = []
                st_ = "insufficient"
            if st_ != head["state"]:
                print("HEAD CHECK: {}".format({
                    "ok": "EEG structure present (r = {:+.2f})",
                    "nohead": "NO HEAD SIGNAL (r = {:+.2f})",
                    "insufficient": "only {} well-contacted channels - "
                                    "need {} to check the cap as a whole"
                }[st_].format(*((float(np.median(head["hist"])),)
                                if st_ != "insufficient" else
                                (hs["n_trusted"], head_min_trusted)))))
            head["state"] = st_
            head["n_trusted"] = hs["n_trusted"]
            cp = set(c for g in hs.get("copies", []) for c in g)
            if len(cp) != copyst["n_shown"]:
                if cp:
                    names = [eeg_labels[c] for c in sorted(cp)]
                    print("IDENTICAL COPIES: {} channels have the same offset "
                          "and signal -- not touching the scalp, or bridged "
                          "({}{})".format(len(names), ", ".join(names[:12]),
                                          ", ..." if len(names) > 12 else ""))
                else:
                    print("IDENTICAL COPIES: cleared")
                copyst["n_shown"] = len(cp)
            copyst["chans"] = cp
            ref = hs.get("ref", {"state": "ok"})
            refst["hist"] = (refst["hist"] + [ref["state"] == "ref"])[-ref_history:]
            on = 2 * sum(refst["hist"]) > len(refst["hist"]) # majority
            if ref["state"] == "ref":
                refst["info"] = ref
            if on != refst["on"]:
                if on:
                    names = [eeg_labels[c] for c in refst["info"]["channels"]]
                    print("REFERENCE/GROUND: {} channels look like EEG but share "
                          "the same ~{:.0f} uV of mains -- check the "
                          "reference and ground electrodes ({}{})".format(
                              len(names), refst["info"]["mains_uv"],
                              ", ".join(names[:12]),
                              ", ..." if len(names) > 12 else ""))
                else:
                    print("REFERENCE/GROUND: cleared")
            refst["on"] = on
        qual = last_qual["vals"]

        disp = dsnap[:, chans].astype(np.float64)
        disp = disp - np.median(disp, axis=0)

        '''
        Display scale from the WELL-CONTACTED channels, anywhere on the
        cap, so good EEG fills its lane and bad channels clip. Scaling to
        the median channel instead meant that on a mostly-dry cap the
        two gelled channels were squashed to flat lines by 585 uV lanes.
        Falls back to this page's median channel if nothing is good yet.
        7 robust sd per lane leaves a healthy channel +/-3.3 sd of room.
        '''
        if state["auto"]:
            good = np.flatnonzero(qual >= 70)
            if good.size < 4 and refst["on"] and refst["info"]:
                # reference problem: the EEG-carrying channels set the scale
                good = np.array(refst["info"]["channels"])
            elif good.size < 4:
                '''
                few good channels: widen to the fair ones too, so one or
                two small-amplitude channels can't shrink the lanes and
                clip the real EEG next to them. Copies are never used.
                '''
                good = np.array([c for c in np.flatnonzero(qual >= 30)
                                 if c not in copyst["chans"]], dtype=int)
            if good.size >= 2:
                g = dsnap[:, good].astype(np.float64)
                g = g - np.median(g, axis=0)
                sd = 1.4826 * float(np.median(np.median(np.abs(g), axis=0)))
            else:
                sd = 1.4826 * float(np.median(np.median(np.abs(disp),
                                                        axis=0)))
            target = max(7.0 * sd, 1.0) # uV, floor at 1 uV
            state["step"] = 0.85 * state["step"] + 0.15 * target
        step = state["step"] / max(state["gain"], 1e-6)
        half = 0.47 * step # lane half-height

        off = np.arange(len(chans)) * step
        clip_frac = np.zeros(len(chans))
        for i in range(channels_per_page):
            if i < len(chans):
                y = disp[:, i]
                if state["clip"]:
                    # Confine the trace to its own lane.
                    clip_frac[i] = float(np.mean(np.abs(y) > half))
                    y = np.clip(y, -half, half)
                lines[i].set_ydata(y + off[i])
                '''
                Grey means BAD CONTACT, not merely "big". A well-gelled
                frontal channel clips during blinks and must keep its
                colour; a dry channel clips because it is junk.
                '''
                bad = qual[chans[i]] < 40 and clip_frac[i] > 0.10
                if refst["on"] and refst["info"] and \
                        chans[i] in refst["info"]["channels"]:
                    bad = False # EEG is there; the mains is the reference's
                if chans[i] in copyst["chans"]:
                    bad = True # not a real electrode signal
                lines[i].set_color("#a3a8ae" if bad else base_colors[i])
                lines[i].set_linewidth(0.6 if bad else 0.8)
                lines[i].set_visible(True)
            else:
                lines[i].set_visible(False)
                qtexts[i].set_text("")

        ax.set_yticks(off)
        ax.set_yticklabels([eeg_labels[c] for c in chans], fontsize=9)
        ax.set_ylim(-step, off[-1] + step)

        ref_chans = (set(refst["info"]["channels"])
                     if refst["on"] and refst["info"] else set())
        if copyst["chans"]:
            head_banner.set_text(
                "{} CHANNELS ARE IDENTICAL COPIES  -  same offset and signal: "
                "not touching the scalp, or bridged".format(len(copyst["chans"])))
            head_banner.get_bbox_patch().set_facecolor("#2d3748")
            head_banner.set_visible(True)
        elif refst["on"]:
            head_banner.set_text(
                "CHECK REFERENCE / GROUND  -  {} channels look like EEG but share "
                "the same {:.0f} uV of {:.0f} Hz".format(
                    len(ref_chans), refst["info"]["mains_uv"], line_freq))
            head_banner.get_bbox_patch().set_facecolor("#6b46c1")
            head_banner.set_visible(True)
        elif head["state"] == "nohead":
            head_banner.set_text(
                "NO HEAD SIGNAL  -  channels show no spatial EEG structure "
                "(r = {:+.2f})".format(float(np.median(head["hist"]))))
            head_banner.get_bbox_patch().set_facecolor("#c53030")
            head_banner.set_visible(True)
        elif head["state"] == "insufficient" and total >= win:
            head_banner.set_text(
                "{} of {} channels well contacted  -  cap-level checks "
                "start at {}".format(int((qual >= trust_min).sum()), n_eeg,
                                     head_min_trusted))
            head_banner.get_bbox_patch().set_facecolor("#b7791f")
            head_banner.set_visible(True)
        else:
            head_banner.set_visible(False)

        for i, ch in enumerate(chans):
            pct = qual[ch]
            qtexts[i].set_position((1.005, off[i]))
            if ch in copyst["chans"]:
                qtexts[i].set_text("  0%  COPY")
                qtexts[i].set_color("#4a5568")
                continue
            if ch in ref_chans:
                qtexts[i].set_text("{:3.0f}%  REF".format(pct))
                qtexts[i].set_color("#6b46c1")
                continue
            qtexts[i].set_text("{:3.0f}%{}".format(
                pct, "  CLIP" if (pct < 40 and clip_frac[i] > 0.10) else ""))
            qtexts[i].set_color(quality_color(pct))

        # scale bar
        x0 = -window_sec * 0.985
        y0 = off[0] - step * 0.95
        scale_line.set_data([x0, x0], [y0, y0 + step])
        scale_text.set_position((x0 + window_sec * 0.012, y0))
        scale_text.set_text("{:.0f} uV".format(step))

        drop = ("  |  LOST {} ({} ev)".format(counter["lost"],
                                             counter["events"])
                if counter["lost"] else "")
        fig.suptitle(
            "Sessantaquattro live  |  page {}/{}  |  ch {}-{}  |  "
            "{} uV/div  |  {}{}{}"
            .format(state["page"] + 1, n_pages, chans[0] + 1, chans[-1] + 1,
                    int(round(step)),
                    "AUTO" if state["auto"] else "MANUAL", drop,
                    "" if not head["hist"] else
                    "  |  head r {:+.2f}".format(
                        float(np.median(head["hist"])))),
            fontsize=10)

        if stop_event.is_set():
            plt.close(fig)
        return lines

    ani = animation.FuncAnimation(fig, update, interval=50, blit=False,
                                  cache_frame_data=False)
    try:
        plt.show()
    except KeyboardInterrupt:
        pass
    finally:
        stop_event.set()
        th.join(timeout=1.0)
        disconnect_from_sq(conn)
        sock.close()
        print("Disconnected.")


if __name__ == "__main__":
    main()
