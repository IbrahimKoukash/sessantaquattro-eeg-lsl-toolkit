'''
counter.py: sample-counter / dropout monitor.

The device sends a free-running counter on one accessory channel.
It is the ONLY way to detect samples lost on the wireless link:
TCP delivers a continuous byte stream with no gap, so without this
check a dropout silently shifts every later sample earlier in time
relative to your other LSL streams.
'''

import numpy as np

from . import runtime as rt
from .config import n_accessory, RAMP_MOD


counter = {"col": None, "last": None, "lost": 0, "events": 0,
           "checked": False}
counter_mod = RAMP_MOD # see config.RAMP_MOD: 16-bit even in 24-bit mode


def find_counter_column(v):
    '''Pick the accessory column that steps by exactly +1 per sample.'''
    if v.shape[0] < 4:
        return None
    for c in range(rt.nch - n_accessory, rt.nch):
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
              .format(lost, 1000.0 * lost / rt.fs,
                      counter["lost"], counter["events"]))
    return lost
