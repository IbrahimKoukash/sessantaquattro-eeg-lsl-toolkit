'''
outlet.py: the LSL outlet, with channel labels, units and electrode
locations in its metadata.
'''

from pylsl import StreamInfo, StreamOutlet

from . import runtime as rt
from .config import samples_per_read
from .montage import electrode_positions


def make_outlet():
    info = StreamInfo(
        name="Sessantaquattro",
        type="EEG",
        channel_count=rt.nch,
        nominal_srate=rt.fs,
        channel_format="float32",
        source_id="sessantaquattro_64ch")

    locs = electrode_positions(rt.eeg_labels)
    chns = info.desc().append_child("channels")

    for lb in rt.eeg_labels:
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

    for lb in rt.aux_labels:
        c = chns.append_child("channel")
        c.append_child_value("label", lb)
        c.append_child_value("unit", "raw")
        # Not AUX inputs -- OTBioLab+ calls these AdapterControl.
        c.append_child_value("type", "Misc")

    acq = info.desc().append_child("acquisition")
    acq.append_child_value("manufacturer", "OT Bioelettronica")
    acq.append_child_value("model", "Sessantaquattro")
    acq.append_child_value("resolution_bits", str(8 * rt.bps))
    acq.append_child_value("lsb_microvolts", "{:.8f}".format(rt.to_uv))
    if not locs:
        print("Note: no electrode coordinates written "
              "(mne not installed, or labels not in standard_1005).")

    return StreamOutlet(info, chunk_size=samples_per_read, max_buffered=360)
