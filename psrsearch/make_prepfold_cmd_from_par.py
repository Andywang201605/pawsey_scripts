#!/usr/bin/env python3
"""
Build a `prepfold` command (with -p / -pd / -pdd) for a binary pulsar (ELL1)
from a .par file.

For a binary pulsar the *apparent* spin period seen at the barycenter is
modulated by the orbit, so over a short fold the period is approximated by
a Taylor series about the observation start t0:

    P(t0 + dt) ~= p + pd*dt + 0.5*pdd*dt^2

where p, pd, pdd are derived from the par file as follows.

Model (ELL1 Roemer delay, first order in EPS1/EPS2):
    Phi(t) = 2*pi*(t - TASC)/PB
    D(t)   = A1*[sin(Phi) + 0.5*(EPS2*sin(2Phi) - EPS1*cos(2Phi))]
    F_int  = F0 + F1*(t-PEPOCH) + 0.5*F2*(t-PEPOCH)^2
    P(t)   = (1 + dD/dt) / F_int(t)

and pd, pdd are the exact first/second time derivatives of P(t) at t0.
(Shapiro delay, orbital-phase-dependent terms of order (v/c)^2 and PBDOT
are neglected.)

Usage example:
    python make_prepfold_cmd.py J1646-44.par --tstart 61263.1234 \
        --tobs 1200 --file obs.fits
"""

import argparse
import sys
import numpy as np

from astropy.io import fits

DAY_S = 86400.0
SECPERDAY = 86400.0


# ----------------------------------------------------------------- par file
def read_par(path):
    par = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            par[parts[0].upper()] = parts[1]
    return par


def getf(par, key, default=0.0):
    return float(par[key].replace("D", "E")) if key in par else default


# ------------------------------------------------------------- orbital model
class BinaryModel:
    def __init__(self, par):
        if par.get("BINARY", "").upper() not in ("ELL1", "ELL1H", "ELL1K"):
            print(f"WARNING: BINARY={par.get('BINARY')}; this script only "
                  f"implements the ELL1 model.", file=sys.stderr)
        self.F0 = getf(par, "F0")
        self.F1 = getf(par, "F1")
        self.F2 = getf(par, "F2")
        self.PEPOCH = getf(par, "PEPOCH")
        self.PB = getf(par, "PB") * DAY_S          # s
        self.A1 = getf(par, "A1")                  # light-s
        self.TASC = getf(par, "TASC")              # MJD
        self.E1 = getf(par, "EPS1")
        self.E2 = getf(par, "EPS2")
        self.w = 2 * np.pi / self.PB

    # derivatives of the Roemer delay D(t): returns D', D'', D'''
    def delay_derivs(self, t_mjd):
        t = np.asarray(t_mjd, dtype=float)
        Phi = self.w * (t - self.TASC) * DAY_S
        s1, c1 = np.sin(Phi), np.cos(Phi)
        s2, c2 = np.sin(2 * Phi), np.cos(2 * Phi)
        x, w, e1, e2 = self.A1, self.w, self.E1, self.E2
        d1 = x * w * (c1 + e2 * c2 + e1 * s2)
        d2 = x * w**2 * (-s1 - 2 * e2 * s2 + 2 * e1 * c2)
        d3 = x * w**3 * (-c1 - 4 * e2 * c2 - 4 * e1 * s2)
        return d1, d2, d3

    def spin_freq(self, t_mjd):
        dt = (np.asarray(t_mjd, dtype=float) - self.PEPOCH) * DAY_S
        return self.F0 + self.F1 * dt + 0.5 * self.F2 * dt**2

    def period(self, t_mjd):
        d1, _, _ = self.delay_derivs(t_mjd)
        return (1.0 + d1) / self.spin_freq(t_mjd)        

    def p_pd_pdd(self, t_mjd):
        """Apparent period and its 1st/2nd derivatives at t_mjd."""
        d1, d2, d3 = self.delay_derivs(t_mjd)
        dt = (t_mjd - self.PEPOCH) * DAY_S
        F = self.F0 + self.F1 * dt + 0.5 * self.F2 * dt**2
        Fd = self.F1 + self.F2 * dt
        Fdd = self.F2
        # h = 1/F ; h' = -F'/F^2 ; h'' = 2F'^2/F^3 - F''/F^2
        h = 1.0 / F
        h1 = -Fd / F**2
        h2 = 2 * Fd**2 / F**3 - Fdd / F**2
        # g = 1 + D'  ; g' = D'' ; g'' = D'''
        g, g1, g2 = 1.0 + d1, d2, d3
        p = g * h
        pd = g1 * h + g * h1
        pdd = g2 * h + 2 * g1 * h1 + g * h2
        return float(p), float(pd), float(pdd)


class PsrFits:
    def __init__(self, filename):
        self._load_file(filename)

    def _load_file(self, filename):
        hdul = fits.open(filename)
        self.bandwidth = hdul[0].header["OBSBW"]        # total bandwidth in MHz
        self.ctrfreq = hdul[0].header["OBSFREQ"]        # center frequency in MHz  
        self.nchan = hdul["SUBINT"].header["NCHAN"]     # number of channels 
        self.npol = hdul["SUBINT"].header["NPOL"]       # number of polarizations

        self.nsubint = hdul["SUBINT"].header["NAXIS2"]  # number of subintegrations
        self.nsblk = hdul["SUBINT"].header["NSBLK"]     # samples per row
        self.nout = self.nsubint * self.nsblk           # spectra per file
        self.tsamp = hdul["SUBINT"].header["TBIN"]      # sample time in seconds

        self.chanwidth = float(self.bandwidth / self.nchan)
        f0 = float(self.ctrfreq) - float(self.bandwidth) / 2.0 + self.chanwidth / 2.0
        self.freqs = np.arange(f0, f0 + self.bandwidth, self.chanwidth)  # frequency array in MHz

        ### load starting mjd...
        imjd = int(hdul[0].header["STT_IMJD"])
        smjd = float(hdul[0].header["STT_SMJD"])
        offs = float(hdul[0].header.get("STT_OFFS", 0.0))

        extra_days, sec_of_day = divmod(smjd + offs, SECPERDAY)
        self.start_imjd = imjd + int(extra_days)       # integer MJD
        self.start_fmjd = sec_of_day / SECPERDAY       # fractional day
        self.tobs = self.nout * self.tsamp             # length in seconds

        hdul.close()

    def __str__(self):
        return f"""Total Bandwidth: {self.bandwidth} MHz
Center Frequency: {self.ctrfreq} MHz
Number of Channels: {self.nchan}
Number of Polarizations: {self.npol}
Sample Time: {self.tsamp} s
Start MJD: {self.start_imjd} + {self.start_fmjd}
Total Observed Time: {self.tobs} s
Spectra per file: {self.nout}
"""

# ------------------------------------------------------------ time handling
def utc_to_bary_tdb_mjd(mjd_utc, ra, dec, lat=None, lon=None, height=None, site=None):
    """Convert a topocentric UTC MJD to barycentric TDB MJD (astropy)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord, EarthLocation
    from astropy.time import Time
    if site is not None:
        site = EarthLocation.of_site(site)
    elif lat is not None and lon is not None and height is not None:
        site = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=height * u.m)
    else:
        raise ValueError("Either site or lat/lon/height must be provided")
    coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame="icrs")
    t = Time(mjd_utc, format="mjd", scale="utc", location=site)
    ltt = t.light_travel_time(coord, kind="barycentric", ephemeris="builtin")
    return float((t.tdb + ltt).mjd)


# -------------------------------------------------------------- validation
def taylor_phase_error(model, t0_mjd, tobs, p, pd, pdd, n=20001):
    """Max |phase difference| (cycles) between exact orbit and prepfold's
    (f, fd, fdd) polynomial over [t0, t0+tobs]."""
    t = np.linspace(0.0, tobs, n)
    F_exact = 1.0 / model.period(t0_mjd + t / DAY_S)
    phi_exact = np.concatenate(
        [[0.0], np.cumsum(0.5 * (F_exact[1:] + F_exact[:-1]) * np.diff(t))])
    f = 1.0 / p
    fd = -pd / p**2
    fdd = (2 * pd**2 - p * pdd) / p**3
    phi_poly = f * t + fd * t**2 / 2 + fdd * t**3 / 6
    return np.max(np.abs(phi_exact - phi_poly))



def max_valid_tobs(model, t0_mjd, p, pd, pdd, tol=0.05):
    """Longest span (s) over which the polynomial stays within tol cycles."""
    lo, hi = 1.0, model.PB
    if taylor_phase_error(model, t0_mjd, hi, p, pd, pdd, n=4001) < tol:
        return hi
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if taylor_phase_error(model, t0_mjd, mid, p, pd, pdd, n=4001) < tol:
            lo = mid
        else:
            hi = mid
    return lo


def make_prepfold_cmd(args, nfold=1, ifold=1):
    par = read_par(args.parfile)
    model = BinaryModel(par)

    ### get start time...
    if args.tstart is None:
        ### load data from fits file...
        psrfits = PsrFits(args.file)
        t0 = psrfits.start_imjd + psrfits.start_fmjd
        filelen = psrfits.tobs
    else:
        t0 = args.tstart
        filelen = args.tobs

    if nfold > 1:
        # t0 = t0 + (ifold - 1) * filelen / nfold / DAY_S
        t0 = t0 + (ifold - 1 / 2) * filelen / nfold / DAY_S
    else:
        t0 = t0 + 0.5 * filelen / DAY_S

    if args.site:
        if "," in args.site:
            lat, lon, h = (float(v) for v in args.site.split(","))
            t0_bary = utc_to_bary_tdb_mjd(t0, par["RAJ"], par["DECJ"], lat, lon, h)
        else:
            t0_bary = utc_to_bary_tdb_mjd(t0, par["RAJ"], par["DECJ"], site=args.site)

        print(f"# topocentric UTC MJD {t0:.9f} -> barycentric TDB MJD "
              f"{t0_bary:.9f} (shift {(t0_bary - t0) * DAY_S:+.2f} s)",
              file=sys.stderr)
        t0 = t0_bary

    p, pd, pdd = model.p_pd_pdd(t0)
    phase = ((t0 - model.TASC) * DAY_S / model.PB) % 1.0
    p_int = 1.0 / model.spin_freq(t0)

    print(f"# orbital phase at start (from TASC): {phase:.4f}", file=sys.stderr)
    print(f"# intrinsic period: {p_int:.15f} s", file=sys.stderr)
    print(f"# apparent period : {p:.15f} s  (x{p / p_int:.8f})", file=sys.stderr)
    print(f"# pd  = {pd:.6e} s/s", file=sys.stderr)
    print(f"# pdd = {pdd:.6e} s/s^2", file=sys.stderr)

    tmax = max_valid_tobs(model, t0, p, pd, pdd)
    print(f"# (p,pd,pdd) stays within 0.05 cycles for ~{tmax:.0f} s "
          f"({tmax / model.PB:.3f} PB) from this start time", file=sys.stderr)

    if args.tobs:
        tobs = filelen / nfold
        err = taylor_phase_error(model, t0, tobs, p, pd, pdd)
        frac = tobs / model.PB
        print(f"# tobs/PB = {frac:.3f}; max phase error of (p,pd,pdd) "
              f"polynomial over tobs = {err:.4f} cycles", file=sys.stderr)
        if err > 0.05:
            print("# WARNING: polynomial is a poor approximation over this "
                  "span. Fold shorter segments (-start/-end), or use "
                  "prepfold's binary options / -par instead.", file=sys.stderr)

    dm = args.dm if args.dm is not None else getf(par, "DM")
    out = args.out or par.get("PSRJ", "psr")
    if nfold > 1:
        out += f"_{nfold}c{ifold}"

    cmd = [
        "prepfold",
        f"-p {p:.15f}",
        f"-pd {pd:.8e}",
        # f"-pdd {pdd:.8e}",
        f"-dm {dm:g}",
    ]
    if args.nsub:
        cmd.append(f"-nsub {args.nsub}")
    if args.npart:
        cmd.append(f"-npart {args.npart}")
    if args.bins:
        cmd.append(f"-n {args.bins}")
    cmd += ["-psrfits -noscales -nooffsets -noweights", "-noxwin", f"-o {out}"]
    if args.extra:
        cmd.append(args.extra)
    if nfold > 1:
        startfrac = (ifold - 1) / nfold
        endfrac = ifold / nfold
        cmd += [f"-start {startfrac:.2f}", f"-end {endfrac:.2f}"]
    cmd.append(args.file)

    print(" ".join(cmd))
    return " ".join(cmd)

def make_sbatch_cmds(cmdlst):
    with open("prepfold_cmds.txt", "w") as f:
        for cmd in cmdlst:
            f.write(cmd + "\n")
    
    cmd = f"""#!/bin/bash
#SBATCH --account=pawsey0272
#SBATCH --output=./logs/%x_%A_%a.out
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --job-name=prepfold
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4

module use /software/projects/pawsey0272/setonix/2025.08/modules/zen3/gcc/14.2.0
module load singularity/4.1.0-psr

export psrsrc="singularity exec /software/projects/pawsey0272/zwan4817/container/psr-search.sif"
export psraly="singularity exec /software/projects/pawsey0272/zwan4817/container/psr-analysis.sif" 

CMD=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" prepfold_cmds.txt)
echo "Running: $CMD"
eval $psrsrc $CMD
"""
    with open("prepfold.sbatch", "w") as f:
        f.write(cmd)


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-par", "--parfile", type=str, required=True, help="par file")
    ap.add_argument("--tstart", type=float, required=False, default=None,
                    help="MJD of the start of the data to fold")
    ap.add_argument("--tobs", type=float, default=None,
                    help="observation length in seconds (used for the "
                         "Taylor-expansion validity check)")
    ap.add_argument("--site", default=None, metavar="LAT,LON,HEIGHT_M/SITE_NAME",
                    help="if given, --tstart is treated as a topocentric UTC "
                         "MJD and converted to barycentric TDB, e.g. "
                         "-32.9984,148.2635,414.8 (Parkes) or "
                         "Parkes (site name)")
    ap.add_argument("--file", required=True, help="data file to fold")
    ap.add_argument("--out", default="blinky", help="output basename (-o)")
    ap.add_argument("--nfold", default=4, help="number of folds to perform")
    ap.add_argument("--dm", type=float, default=None,
                    help="override DM from par file")
    ap.add_argument("--nsub", type=int, default=64)
    ap.add_argument("--npart", type=int, default=64)
    ap.add_argument("--bins", type=int, default=32, help="-n profile bins")
    ap.add_argument("--extra", default="-npfact 4 -ndmfact 2 -pstep 1 -pdstep 1 -dmstep 3", help="extra prepfold options")
    args = ap.parse_args()

    cmdlst = []
    if args.nfold > 1:
        for i in range(args.nfold):
            cmdlst.append(make_prepfold_cmd(args, nfold=args.nfold, ifold=i + 1))

    make_sbatch_cmds(cmdlst)


if __name__ == "__main__":
    main()
    # psrfitspath = "/scratch/pawsey0272/zwan4817/J1646/260803/add/J1646.260803.sf"
    # psrfits = PsrFits(psrfitspath)
    # print(psrfits)