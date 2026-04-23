from gwpy.io.gwf import get_channel_names
import re
channels = get_channel_names('/tmp/nh_ch_probe/V-trend-1448928000-86400.gwf')
ni = [c for c in channels if '_NI_' in c or c.startswith('V1:NI_')]
print(f'Total NI: {len(ni)}')
te  = [c for c in ni if re.search(r'TE[0-9]|TEMP|therm|BODY', c, re.I)]
pwr = [c for c in ni if re.search(r'RH|CURR|VOLT|PWR|HEAT|_I_|_P_|INF_NI', c, re.I)]
print('\n=== ALL NI ===')
for c in sorted(ni): print(c)
print('\n=== TEMPERATURE ===')
for c in sorted(te): print(c)
print('\n=== POWER/CURRENT/RH ===')
for c in sorted(pwr): print(c)
