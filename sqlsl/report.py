'''
report.py: the printed quality table (the "Quality" button / 'q' key).
'''

import numpy as np

from . import runtime as rt
from .config import line_freq
from .cap_checks import assess_cap


def print_quality(snapshot):
    n_eeg, eeg_labels = rt.n_eeg, rt.eeg_labels
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
