'''
viewer.py: the live matplotlib viewer (main thread).

16 channels per page, every trace confined to its own lane, quality %
per channel, cap-level banners (copies / reference / no head).
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Button

from . import runtime as rt
from .config import (window_sec, channels_per_page, line_freq, trust_min,
                     head_min_trusted, head_r_min, head_history, ref_history)
from .counter import counter
from .cap_checks import assess_cap
from .quality import quality_color
from .report import print_quality


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


def run_viewer():
    '''
    Build the figure and run the animation until the window is closed or
    the acquisition thread stops. Blocks in plt.show(); the caller does
    the cleanup.
    '''
    # Read after any reconfigure(), so these match the detected layout.
    n_eeg = rt.n_eeg
    win = rt.win
    fs = rt.fs

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
        with rt.ring_lock:
            snap = rt.ring_snapshot()
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
        print("  terminal (python3 main.py), or select a GUI")
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
        eeg_labels = rt.eeg_labels
        with rt.ring_lock:
            snap = rt.ring_snapshot()
            dsnap = rt.disp_snapshot()
            total = rt.sample_counter["n"]

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

        if rt.stop_event.is_set():
            plt.close(fig)
        return lines

    ani = animation.FuncAnimation(fig, update, interval=50, blit=False,
                                  cache_frame_data=False)
    state["_ani"] = ani # keep the animation alive
    plt.show()
