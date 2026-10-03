'''
main.py: Sessantaquattro (OT Bioelettronica, 64 ch) -> LSL,
         with a live viewer and per-channel quality scores

Acquisition thread: TCP from the device -> layout check -> decode ->
                     LSL outlet (raw, microvolts) + continuous display
                     filter -> ring buffers
Main thread: matplotlib viewer, 16 channels per page, every trace
                     confined to its own lane, quality % per channel

Run: python3 main.py
Settings live in sqlsl/config.py. Documentation: README.md
'''

import socket
import threading

from sqlsl import runtime as rt
from sqlsl import montage
from sqlsl.config import host, port
from sqlsl.protocol import disconnect_from_sq
from sqlsl.decoding import detect_channel_count
from sqlsl.counter import counter
from sqlsl.outlet import make_outlet
from sqlsl.acquisition import acquisition
from sqlsl.viewer import run_viewer


def main():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind((host, port))
    sock.listen(1)
    print("Waiting for Sessantaquattro connection...")
    conn, addr = sock.accept()
    print("Connected from:", addr)
    conn.send(rt.command)
    print("Start command sent.")
    print("(If no data flows, press the Sessantaquattro pushbutton.)")

    # Verify the stream layout before a single sample is trusted.
    print("Probing stream layout (expecting {} channels)...".format(rt.nch))
    found, ramp_col, primed = detect_channel_count(conn)
    if found is None:
        print("\n  WARNING: could not find the Ramp counter in any of the")
        print("  layouts tried, so the channel count could not be confirmed.")
        print("  Continuing with {} channels, but if the traces look wrong"
              .format(rt.nch))
        print("  this is the first thing to suspect.")
        primed = b""
    else:
        if found != rt.nch:
            print("\n  Device is sending {} channels, not {}. Adjusting."
                  .format(found, rt.nch))
            rt.reconfigure(found)
        else:
            print("  Confirmed: {} channels, Ramp counter on index {}."
                  .format(found, ramp_col))
        counter["col"] = ramp_col
        counter["checked"] = True

    montage.build_neighbours()
    if montage._neighbours is not None:
        print("  neighbour map: {} nearest electrodes per channel "
              "(e.g. {} -> {})".format(
                  montage._neighbours.shape[1], rt.eeg_labels[0],
                  ", ".join(rt.eeg_labels[j] for j in montage._neighbours[0])))
    else:
        print("  neighbour map: no electrode coordinates; correlation "
              "will use all channels instead of nearest neighbours.")

    outlet = make_outlet()
    print("LSL outlet 'Sessantaquattro' created "
          "({} ch, {} EEG, {} Hz, {}-bit, {:.6f} uV/bit)."
          .format(rt.nch, rt.n_eeg, rt.fs, 8 * rt.bps, rt.to_uv))
    print("Streaming to LSL + viewer...")

    th = threading.Thread(target=acquisition,
                          args=(conn, outlet, primed), daemon=True)
    th.start()

    try:
        run_viewer()
    except KeyboardInterrupt:
        pass
    finally:
        rt.stop_event.set()
        th.join(timeout=1.0)
        disconnect_from_sq(conn)
        sock.close()
        print("Disconnected.")


if __name__ == "__main__":
    main()
