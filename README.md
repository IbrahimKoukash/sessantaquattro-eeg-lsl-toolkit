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




