import subprocess, warnings
warnings.filterwarnings("ignore")
gwf = "/tmp/check_sr_ty_test.gwf"
r = subprocess.run(["rfcp", "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend/2026/V-trend-1456790400-86400.gwf", gwf], capture_output=True)
if r.returncode != 0:
    print("rfcp failed:", r.stderr.decode()[:200])
else:
    from gwpy.io.gwf import iter_channel_names
    channels = [c for c in iter_channel_names(gwf) if "ASC_SR_T" in c]
    for c in sorted(channels):
        print(c)
