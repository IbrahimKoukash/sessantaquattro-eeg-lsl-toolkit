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

