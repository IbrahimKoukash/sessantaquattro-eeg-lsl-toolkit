# sessantaquattro_eeg_lsl_toolkit

Stream a 64-channel **OT Bioelettronica Sessantaquattro** to [Lab Streaming Layer](https://labstreaminglayer.org) (LSL), with a live viewer that scores the contact quality of every electrode.

- **Raw stream to LSL** in microvolts, with real electrode labels and 3-D positions in the stream metadata, so a recording opens in MNE or EEGLAB with its montage attached.
- **A live viewer** showing 16 channels per page. Each trace stays in its own row however noisy it is. Badly contacted channels are greyed out, and each channel has a quality percentage.
- **Contact quality you can watch while gelling.** Each electrode is judged on its own signal, so channels turn green one at a time as you apply gel, even when most of the cap is still dry.
- **Data integrity checks.** It confirms the stream layout on connect, detects samples lost on the wireless link, and warns when there is no head under the cap.

---

## Contents
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Using the viewer](#using-the-viewer)
- [How quality is scored](#how-quality-is-scored)
- [The LSL stream](#the-lsl-stream)
- [Configuration](#configuration)
- [How it works](#how-it-works)
- [Troubleshooting](#troubleshooting)
- [Companion tools](#companion-tools)
- [Known limitations](#known-limitations)
- [References](#references)

---

## Requirements
| | Version | Needed for |
|---|---|---|
| Python | 3.8+ | |
| numpy | 1.21+ | |
| scipy | 1.8+ | filtering, spectra |
| matplotlib | 3.5+ | the viewer (needs an interactive backend, e.g. TkAgg or QtAgg) |
| pylsl | 1.16+ | the LSL outlet |
| mne | 1.0+ | *optional*: electrode positions for the stream metadata and the neighbour-based quality checks |

```bash
pip install -r requirements.txt
```

`pylsl` bundles the LSL library on Windows and macOS. On Linux, if `pylsl` cannot find `liblsl`, install it with `conda install -c conda-forge liblsl`.

Hardware: a Sessantaquattro with a 64-channel EEG cap connected in **monopolar** mode. The cap used during development was ELSCH064EEG on the ADEEGCap64SE adapter.

---

## Quick start

1. **Put the computer and the Sessantaquattro on the same Wi-Fi network.** The computer acts as the server: it listens on **TCP port 45454** and the device connects to it (this is also how OT Bioelettronica's reference code works). On Windows, allow Python through the firewall for incoming connections the first time it asks.

2. **Start the viewer from a terminal:**
   ```bash
   python sq_lsl_viewer.py
   ```
   Run it from a terminal, not a notebook, Spyder or VS Code's interactive window. Those use a non-interactive plot backend, so the window won't respond to clicks or keys, and the viewer warns about this at startup.

3. **Turn on / connect the Sessantaquattro.** A healthy start looks like:
   ```
   Waiting for Sessantaquattro connection...
   Connected from: ('192.168.1.2', 50231)
   Start command sent.
   Probing stream layout (expecting 68 channels)...
     Confirmed: 68 channels, Ramp counter on index 67.
     neighbour map: 8 nearest electrodes per channel (e.g. Fp1 -> AF3, AF7, Fpz, AFz, F7, Fp2, F3, AF4)
   LSL outlet 'Sessantaquattro' created (68 ch, 64 EEG, 500 Hz, 24-bit, 0.286000 uV/bit).
   Streaming to LSL + viewer...
   HEAD CHECK: EEG structure present (r = +0.55)
   ```

4. **Record** with [LabRecorder](https://github.com/labstreaminglayer/App-LabRecorder): select the `Sessantaquattro` stream (type `EEG`) and start recording. The viewer only displays; recording is LabRecorder's job.

Closing the viewer window sends the stop command to the device and closes the connection.

---

## Using the viewer

### Controls
Every action has an on-screen button. Keys work too, once you've clicked on the plot to give it focus.

| Button | Keys | Action |
|---|---|---|
| `< Prev` / `Next >` | `←` `→` `↑` `↓`, `PgUp` `PgDn`, `p`/`n`, `k`/`j`, `Space`, mouse wheel | Previous / next page of 16 channels |
| | `Home` / `End` | First / last page |
| `Zoom +` / `Zoom -` | `+` / `-` | Scale traces up / down (×1.3 per step) |
| `Auto` | `a` | Toggle automatic scaling |
| `Reset` | `0` | Reset zoom to 1× |
| `Quality` | `q` | Print a full per-channel quality table to the terminal |
| | `c` | Toggle keeping each trace inside its own row (default on) |

### Reading the screen

| Element | Meaning |
|---|---|
| **Title bar** | Page, channel range, µV per row, AUTO/MANUAL scaling, lost-sample count (only shown if samples were lost), and `head r` (see [cap-level checks](#cap-level-checks)) |
| **Scale bar** (bottom left) | The height of one row, in µV |
| **Percentage** beside each channel | Contact quality: **green ≥ 70**, **orange 40–69**, **red < 40**. Updated every 2 s for all 64 channels |
| **Grey trace + `CLIP`** | The channel keeps hitting the edge of its row **and** scores below 40%: a bad contact. A good channel that briefly clips (a blink on a gelled frontal electrode) keeps its colour |
| **Amber banner** | Fewer than 8 channels are well contacted, so cap-level checks are not running yet (normal while gelling) |
| **Red banner** | `NO HEAD SIGNAL`: enough channels look clean on their own, but together they show none of the spatial structure a head produces (an empty cap, or a disconnected reference/ground) |

The display is filtered 1–45 Hz with notches at the mains frequency and its first harmonic, for viewing only. **The LSL stream is never filtered.** The display scale follows the well-contacted channels, so good EEG fills its row and poor channels clip against the edges.

<img width="875" height="625" alt="image" src="https://github.com/user-attachments/assets/3376c855-d7f0-4dd0-b88e-00706173e43d" />

*A real session with only F7, F9 and AF7 gelled. Those three score 99–100% while the 61 dry channels are greyed out, and the amber banner says cap-level checks are waiting for 8 good channels.*

---

## How quality is scored

Each channel's percentage is the **lowest** of its sub-scores, so one serious problem can't be hidden by good results elsewhere. Every channel is first judged on **its own signal**, without reference to the other channels. Cross-channel checks are added only when enough channels are well contacted to compare against.

### Per-channel sub-scores (5-second window)

| Sub-score | Measures | Full marks | Zero |
|---|---|---|---|
| **amp** | Typical amplitude: robust SD (1.4826 × MAD) of the 1–45 Hz signal. Blinks barely move it; continuous noise does | 3–50 µV | ≤ 1 µV or ≥ 150 µV |
| **line** | Mains pickup, both **absolute** and **relative to the channel's EEG** (the worse of the two) | ≤ 30 µV and ≤ 1× EEG | ≥ 300 µV or ≥ 10× EEG |
| **hf** | 40–100 Hz power, with the mains band removed, relative to 0.5–40 Hz power (muscle, broadband noise) | ≤ 0.35 | ≥ 1.5 |
| **corr** | Agreement with neighbouring electrodes (below). Only computed when possible | see below | |

Any sample within 5% of full scale scores 0 (clipping at the amplifier).

