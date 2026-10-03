'''
sqlsl: Sessantaquattro (OT Bioelettronica, 64 ch) -> LSL library.

Modules
  config: all user-tunable constants (labels, network, thresholds)
  protocol: command word, LSB scale, disconnect
  runtime: shared mutable state: layout, ring buffers, lock, stop flag
  decoding: socket reads, byte decoding, channel-count probe
  counter: Ramp sample counter / dropout monitor
  filters: display filters (causal live + zero-phase analysis)
  montage: electrode coordinates and nearest-neighbour map
  outlet: LSL outlet with channel metadata
  quality: per-channel own-signal quality scores
  cap_checks: neighbour correlation, copies, reference check, assess_cap
  report: printed quality table
  acquisition: acquisition thread
  viewer: matplotlib live viewer
'''
