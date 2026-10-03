'''
acquisition.py: the acquisition thread.

TCP from the device -> decode -> dropout check -> LSL outlet (raw,
microvolts) + continuous display filter -> ring buffers.
'''

from pylsl import local_clock

from . import runtime as rt
from .config import samples_per_read
from .decoding import read_block, decode_block
from .counter import check_counter
from .filters import disp_filter_block


def acquisition(conn, outlet, primed=b""):
    try:
        # Samples already pulled off the socket by the layout probe.
        if primed:
            m = len(primed) // (rt.nch * rt.bps)
            if m:
                block, raw = decode_block(primed[:m * rt.nch * rt.bps], m)
                check_counter(raw)
                outlet.push_chunk(block.tolist(), local_clock())
                dblock = disp_filter_block(block[:, :rt.n_eeg])
                with rt.ring_lock:
                    rt.ring_write(block, dblock)
                    rt.sample_counter["n"] += m

        while not rt.stop_event.is_set():
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
            dblock = disp_filter_block(block[:, :rt.n_eeg])
            with rt.ring_lock:
                rt.ring_write(block, dblock)
                rt.sample_counter["n"] += block.shape[0]
    except (ConnectionError, OSError) as e:
        print("Acquisition stopped:", e)
    finally:
        rt.stop_event.set()
