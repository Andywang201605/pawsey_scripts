#!/usr/bin/env python3
"""
Concatenate SIGPROC filterbank files ALONG FREQUENCY into one filterbank,
processing in time chunks so that memory use stays small.

Features
--------
* Files may start at slightly different times (tstart): all files are trimmed
  to the common overlapping time range (offset rounded to whole samples).
* Files may cover different frequency ranges; channel grids must have the same
  channel width and be aligned to the same grid (within tolerance).
* Per-file edge-channel flagging (and arbitrary channel ranges). Flagged
  channels are replaced by a fill value (median / mean / zero / constant).
* Gaps between sub-bands are filled with the same fill value, so the output
  has a uniform frequency axis.
* Overlapping channels: the unflagged one wins. If two files both have
  unflagged data in the same channel, an error is raised (flag one of them).
* Reads/writes in chunks; supports 8/16/32-bit data, nifs == 1.

Usage
-----
    python fil_concat.py a.fil b.fil c.fil -o out.fil \
        --flag-default 4,4 \
        --flag 1:10,2 \
        --flag-range 2:100-120 \
        --fill median --chunk-mb 256

    --flag-default LO,HI   flag first LO and last HI channels (as stored in the
                           file) of every file without its own --flag
    --flag IDX:LO,HI       same, for file IDX (0-based, order on command line)
    --flag-range IDX:A-B   flag channels A..B inclusive (stored order), repeatable

As a library:
    from fil_concat import concat_filterbanks
    concat_filterbanks(paths, "out.fil", flags={0: dict(lo=5, hi=5)})
"""
import argparse
import struct
import sys
from collections import OrderedDict

import numpy as np

INT_KEYS = {"machine_id", "telescope_id", "data_type", "nchans", "nbits",
            "nifs", "scan_number", "barycentric", "pulsarcentric",
            "nbeams", "ibeam", "nsamples"}
DOUBLE_KEYS = {"az_start", "za_start", "src_raj", "src_dej", "tstart",
               "tsamp", "fch1", "foff", "refdm", "period", "fchannel",
               "FREQUENCY_START", "FREQUENCY_END"}
STR_KEYS = {"source_name", "rawdatafile"}
CHAR_KEYS = {"signed"}


# --------------------------------------------------------------------------
# Header I/O
# --------------------------------------------------------------------------
def _read_str(f):
    n = struct.unpack("<i", f.read(4))[0]
    if not 0 < n < 256:
        raise ValueError("Not a valid SIGPROC header (bad string length)")
    return f.read(n).decode("ascii")


def _write_str(f, s):
    b = s.encode("ascii")
    f.write(struct.pack("<i", len(b)) + b)


def read_header(path):
    """Return (header OrderedDict, header size in bytes)."""
    hdr = OrderedDict()
    with open(path, "rb") as f:
        if _read_str(f) != "HEADER_START":
            raise ValueError(f"{path}: missing HEADER_START")
        while True:
            key = _read_str(f)
            if key == "HEADER_END":
                break
            if key in STR_KEYS:
                hdr[key] = _read_str(f)
            elif key in INT_KEYS:
                hdr[key] = struct.unpack("<i", f.read(4))[0]
            elif key in DOUBLE_KEYS:
                hdr[key] = struct.unpack("<d", f.read(8))[0]
            elif key in CHAR_KEYS:
                hdr[key] = struct.unpack("<b", f.read(1))[0]
            else:
                raise ValueError(f"{path}: unknown header key '{key}'")
        size = f.tell()
    return hdr, size


def write_header(f, hdr):
    _write_str(f, "HEADER_START")
    for k, v in hdr.items():
        _write_str(f, k)
        if k in STR_KEYS:
            _write_str(f, v)
        elif k in INT_KEYS:
            f.write(struct.pack("<i", int(v)))
        elif k in DOUBLE_KEYS:
            f.write(struct.pack("<d", float(v)))
        elif k in CHAR_KEYS:
            f.write(struct.pack("<b", int(v)))
    _write_str(f, "HEADER_END")


def _dtype(hdr):
    nbits = hdr["nbits"]
    if nbits == 8:
        return np.dtype(np.int8 if hdr.get("signed", 0) else np.uint8)
    if nbits == 16:
        return np.dtype("<u2")
    if nbits == 32:
        return np.dtype("<f4")
    raise NotImplementedError(f"nbits={nbits} not supported (use 8/16/32)")


# --------------------------------------------------------------------------
# Main routine
# --------------------------------------------------------------------------
class _Band:
    """One input file, with orientation matched to the output frequency axis."""

    def __init__(self, path, idx, flag):
        self.path, self.idx = path, idx
        self.hdr, self.hsize = read_header(path)
        h = self.hdr
        if h.get("nifs", 1) != 1:
            raise NotImplementedError(f"{path}: only nifs=1 is supported")
        self.dtype = _dtype(h)
        self.nchan = h["nchans"]
        self.rowbytes = self.nchan * self.dtype.itemsize
        import os
        self.nsamp = (os.path.getsize(path) - self.hsize) // self.rowbytes
        self.flag = flag or {}

    def stored_flag_mask(self):
        m = np.zeros(self.nchan, bool)
        lo, hi = self.flag.get("lo", 0), self.flag.get("hi", 0)
        if lo:
            m[:lo] = True
        if hi:
            m[self.nchan - hi:] = True
        for a, b in self.flag.get("ranges", []):
            m[max(a, 0):b + 1] = True
        return m


def concat_filterbanks(paths, output, flags=None, default_flag=None,
                       chunk_mb=256, fill="median", verbose=True):
    """
    paths        : list of input .fil paths
    output       : output path
    flags        : {file_index: {"lo": int, "hi": int, "ranges": [(a, b), ...]}}
                   channel numbers refer to the channel order stored in that file
    default_flag : same dict form, used for files not in `flags`
    fill         : 'median' | 'mean' | 'zero' | number
    """
    flags = flags or {}
    bands = [_Band(p, i, flags.get(i, default_flag)) for i, p in enumerate(paths)]
    log = print if verbose else (lambda *a, **k: None)

    # ---- consistency checks --------------------------------------------
    ref = bands[0]
    for b in bands:
        if b.hdr["nbits"] != ref.hdr["nbits"]:
            raise ValueError(f"nbits differ: {b.path}")
        if abs(b.hdr["tsamp"] / ref.hdr["tsamp"] - 1) > 1e-6:
            raise ValueError(f"tsamp differs: {b.path}")
        if abs(abs(b.hdr["foff"]) / abs(ref.hdr["foff"]) - 1) > 1e-4:
            raise ValueError(f"channel width differs: {b.path}")
    tsamp = ref.hdr["tsamp"]
    dtype = ref.dtype
    foff_out = ref.hdr["foff"]            # output keeps first file's direction
    sign = np.sign(foff_out)

    # ---- frequency layout ----------------------------------------------
    for b in bands:
        b.flip = np.sign(b.hdr["foff"]) != sign
        if b.flip:   # view channels in output frequency order
            b.f_start = b.hdr["fch1"] + (b.nchan - 1) * b.hdr["foff"]
        else:
            b.f_start = b.hdr["fch1"]
    f0 = (min if sign > 0 else max)(b.f_start for b in bands)
    for b in bands:
        k = (b.f_start - f0) / foff_out
        b.k0 = int(round(k))
        if abs(k - b.k0) > 0.1:
            raise ValueError(f"{b.path}: channel grid misaligned with others "
                             f"by {abs(k - b.k0):.2f} channels")
    nch_out = max(b.k0 + b.nchan for b in bands)

    claims = np.zeros(nch_out, int)
    for b in bands:
        m = b.stored_flag_mask()
        b.flagged = m[::-1] if b.flip else m
        claims[b.k0:b.k0 + b.nchan] += ~b.flagged
    if (claims > 1).any():
        bad = np.where(claims > 1)[0]
        raise ValueError(f"{len(bad)} output channels (e.g. {bad[:5].tolist()}) "
                         "have unflagged data in more than one file: bands "
                         "overlap; flag the overlapping edge channels")
    bad_mask = claims == 0
    for b in bands:  # contiguous runs of usable channels (in oriented order)
        pad = np.concatenate([[0], (~b.flagged).astype(np.int8), [0]])
        d = np.diff(pad)
        b.runs = list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))

    # ---- time alignment ------------------------------------------------
    t0 = max(b.hdr["tstart"] for b in bands)
    for b in bands:
        exact = (t0 - b.hdr["tstart"]) * 86400.0 / tsamp
        b.skip = int(round(exact))
        log(f"[{b.idx}] {b.path}: nchan={b.nchan} nsamp={b.nsamp} "
            f"skip={b.skip} samples (residual {b.skip - exact:+.3f} samp) "
            f"flagged={int(b.flagged.sum())}")
    nout = min(b.nsamp - b.skip for b in bands)
    if nout <= 0:
        raise ValueError("files have no common time range")

    # ---- output header -------------------------------------------------
    hdr = OrderedDict(ref.hdr)
    hdr.update(nchans=nch_out, nifs=1, fch1=f0, foff=foff_out, tstart=t0)
    if "nsamples" in hdr:
        hdr["nsamples"] = nout
    for k in ("FREQUENCY_START", "FREQUENCY_END", "fchannel"):
        hdr.pop(k, None)
    log(f"Output: {nch_out} channels, {nout} samples, fch1={f0}, foff={foff_out}, "
        f"{int(bad_mask.sum())} filled channels")

    chunk = max(1, int(chunk_mb * 2**20 // (nch_out * dtype.itemsize)))
    handles = [open(b.path, "rb") for b in bands]
    try:
        with open(output, "wb") as fo:
            write_header(fo, hdr)
            for start in range(0, nout, chunk):
                n = min(chunk, nout - start)
                out = np.empty((n, nch_out), dtype)
                for b, fh in zip(bands, handles):
                    fh.seek(b.hsize + (b.skip + start) * b.rowbytes)
                    data = np.fromfile(fh, dtype, count=n * b.nchan)
                    if data.size != n * b.nchan:
                        raise IOError(f"short read from {b.path}")
                    data = data.reshape(n, b.nchan)
                    if b.flip:
                        data = data[:, ::-1]
                    for s, e in b.runs:
                        out[:, b.k0 + s:b.k0 + e] = data[:, s:e]
                if bad_mask.any():
                    if fill == "zero":
                        val = 0
                    elif isinstance(fill, (int, float)):
                        val = fill
                    else:
                        sub = out[::max(1, n // 64)][:, ~bad_mask]
                        val = np.median(sub) if fill == "median" else sub.mean()
                    out[:, bad_mask] = np.asarray(val).astype(dtype)
                out.tofile(fo)
                log(f"  wrote samples {start}-{start + n} / {nout}")
    finally:
        for fh in handles:
            fh.close()
    return output


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def _parse_lohi(s):
    lo, hi = (int(x) for x in s.split(","))
    return dict(lo=lo, hi=hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--flag-default", type=_parse_lohi, metavar="LO,HI")
    ap.add_argument("--flag", action="append", default=[], metavar="IDX:LO,HI")
    ap.add_argument("--flag-range", action="append", default=[], metavar="IDX:A-B")
    ap.add_argument("--fill", default="median",
                    help="median | mean | zero | <number>  (default median)")
    ap.add_argument("--chunk-mb", type=float, default=256)
    a = ap.parse_args()

    flags = {}
    for s in a.flag:
        i, v = s.split(":")
        flags.setdefault(int(i), {}).update(_parse_lohi(v))
    for s in a.flag_range:
        i, v = s.split(":")
        lo, hi = (int(x) for x in v.split("-"))
        flags.setdefault(int(i), {}).setdefault("ranges", []).append((lo, hi))
    # files with only --flag-range still inherit default lo/hi
    if a.flag_default:
        for i, d in flags.items():
            for k, v in a.flag_default.items():
                d.setdefault(k, v)

    try:
        fill = float(a.fill)
    except ValueError:
        fill = a.fill
    concat_filterbanks(a.files, a.output, flags=flags,
                       default_flag=a.flag_default, chunk_mb=a.chunk_mb,
                       fill=fill)


if __name__ == "__main__":
    sys.exit(main())