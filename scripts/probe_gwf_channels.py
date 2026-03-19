"""
Probe a single GWF file: stage it from HPSS and list all channels
matching a pattern. Useful for checking channel naming conventions.

Usage:
    python probe_gwf_channels.py --gps-day 1364774400 --pattern HWS
    python probe_gwf_channels.py --gps-day 1364774400 --pattern INF_NI
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--gps-day", type=int, required=True,
                   help="GPS start of day (integer, multiple of 86400)")
    p.add_argument("--pattern", default="HWS",
                   help="Substring to filter channel names (case-insensitive)")
    p.add_argument("--hpss-base",
                   default="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend")
    p.add_argument("--stage-dir", default=None)
    return p.parse_args()


def main():
    args = parse_args()
    gps = args.gps_day
    year = 2000 + (gps - 630720013) // 31557600  # rough GPS→year
    # More accurate: use datetime
    import datetime
    t0 = datetime.datetime(1980, 1, 6) + datetime.timedelta(seconds=gps - 18)
    year = t0.year

    remote = f"{args.hpss_base}/{year}/V-trend-{gps}-86400.gwf"
    stage_dir = Path(args.stage_dir) if args.stage_dir else Path(tempfile.mkdtemp(prefix="probe_gwf_"))
    gwf = stage_dir / f"V-trend-{gps}-86400.gwf"

    print(f"Staging {remote} → {gwf}", flush=True)
    ret = subprocess.run(["rfcp", remote, str(gwf)], capture_output=True)
    if ret.returncode != 0 or not gwf.exists():
        sys.exit(f"rfcp failed: {ret.stderr.decode()}")
    print(f"Staged OK ({gwf.stat().st_size/1e6:.1f} MB)", flush=True)

    # List channels using FrChannels (part of FrameLib/framel)
    # fallback: gwpy
    pat = args.pattern.upper()
    try:
        result = subprocess.run(
            ["FrChannels", str(gwf)],
            capture_output=True, text=True
        )
        channels = [l.split()[0] for l in result.stdout.splitlines() if pat in l.upper()]
        if channels:
            print(f"\nChannels matching '{pat}':")
            for ch in sorted(channels):
                print(f"  {ch}")
            print(f"\nTotal: {len(channels)}")
            return
        elif result.returncode == 0:
            all_ch = result.stdout.splitlines()
            print(f"No channels matching '{pat}' (total channels: {len(all_ch)})")
            return
    except FileNotFoundError:
        print("FrChannels not found, trying gwpy...", flush=True)

    # gwpy fallback
    try:
        from gwpy.io.gwf import iter_channel_names
        channels = [ch for ch in iter_channel_names(str(gwf)) if pat in ch.upper()]
        print(f"\nChannels matching '{pat}':")
        for ch in sorted(channels):
            print(f"  {ch}")
        print(f"\nTotal: {len(channels)}")
    except Exception as e:
        sys.exit(f"gwpy channel listing failed: {e}")


if __name__ == "__main__":
    main()
