'''
decoding.py: reading and decoding the raw TCP byte stream, plus the
channel-count / alignment probe.
'''

import numpy as np
from pylsl import local_clock

from . import runtime as rt
from .config import n_accessory, RAMP_MOD


def read_block(conn, nsamp):
    '''
    Read exactly nsamp samples worth of bytes.

    Returns (buf, t) where t is taken the instant the last byte lands,
    before any decoding. Decoding is numpy work that can be preempted by
    the GUI thread, and any delay between arrival and local_clock()
    inflates the timestamp for the whole chunk.
    '''
    size = nsamp * rt.nch * rt.bps
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
    nch, bps = rt.nch, rt.bps
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
    out[:, :rt.n_eeg] *= rt.to_uv # AUX and accessory stay raw
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
    bps = rt.bps
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
