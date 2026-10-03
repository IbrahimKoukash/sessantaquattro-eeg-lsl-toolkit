'''
protocol.py: Sessantaquattro command word and connection helpers.
'''

import socket

from .config import RESOLUTION_BITS


def integer_to_bytes(command):
    return int(command).to_bytes(2, byteorder="big")


def create_bin_command(start=1, imp=0):
    '''
    Build the 2-byte Sessantaquattro command word.

    Field meanings follow OT Bioelettronica's reference script
    (OTB-Matlab):
      start (GO): bit 0 / 1 = send settings and start data transfer
      rec: bit 1 / 1 = record to the device's SD card
      trig: bits 2-3 / 0 = transfer/SD controlled remotely,
                       3 = SD recording from the pushbutton
      ext: (EXTEN) bits 4-5 / INPUT RANGE: 0 standard, 1 x2, 2 x4, 3 x8.
                              Keep 0: the 0.286 uV/count factor assumes it.
      hpf:bit 6  / 0 = DC coupled, 1 = hardware high-pass on
      hres: bit 7 / 0 = 16-bit, 1 = 24-bit samples
      mode: bits 8-10 / 0 monopolar, 1 bipolar, 2 differential,
                        3 accelerometers, 6 impedance check, 7 test
      nch: bits 11-12 / 0/1/2/3 = 8/16/32/64 channels
      fsamp: bits 13-14 /  0/1/2 = 500/1000/2000 Hz (mode != 3)
    imp=1 selects impedance-check mode (6). The viewer does not decode
    impedance data; it exists only so a caller can request it.
    '''
    rec = 0
    trig = 0 # data transfer controlled remotely by this command
    ext = 0 # standard input range (other values change the scaling)
    hpf = 0 # DC coupled
    hres = 1 if RESOLUTION_BITS == 24 else 0
    mode = 0 if imp == 0 else 6 # 0 monopolar; 6 impedance check
    nch = 3 # 64 biopotential inputs
    fsamp = 0 # 500 Hz
    getset = 0
    command = (start + rec * 2 + trig * 4 + ext * 16 + hpf * 64
               + hres * 128 + mode * 256 + nch * 2048
               + fsamp * 8192 + getset * 32768)
    bps = 3 if hres else 2

    # Total channels in the stream: 64 EEG + 2 AUX + Buffer + Ramp = 68.
    # OT Bioelettronica's reference script uses NumChan = 68 for 64-channel
    total = 68
    return integer_to_bytes(command), total, 500, bps


def lsb_microvolts(bps):
    '''
    286 nV per bit, for both resolutions.

    OT Bioelettronica's reference script: it applies ConvFact = 0.000286 mV to 24-bit
    and 16-bit data alike (OTB-Matlab).
    '''
    return 0.286


def disconnect_from_sq(conn):
    if conn is None:
        return
    stop, _, _, _ = create_bin_command(start=0)
    try:
        conn.send(stop)
        conn.shutdown(socket.SHUT_RDWR)
    except OSError:
        pass
    conn.close()
