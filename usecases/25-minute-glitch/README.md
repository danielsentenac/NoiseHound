# 25-Minute Glitches — NOISEHOUND Investigation

## Introduction

The 25-minute glitches are a persistent family of broadband noise transients in the Virgo strain channel (Hrec_hoft): SNR ~400, peak frequency 40–50 Hz, recurrence interval ~25–32 min (seasonal). First documented before the O4b run start, they survived every hardware intervention through early 2026 and their origin remained unidentified until this investigation. Around March 2026 the BruCo online monitoring was disabled for this family as the glitches appear to have disappeared.

An independent cross-check is provided by the EXCAVATor tool (Swinkels, Nikhef), which ran on a 10-hour epoch in April 2023 (GPS 1366750818–1366786818, 26 triggers, 11 807 channels on 1 Hz trend). EXCAVATor independently confirms that LSC/DARM channels dominate (BPC_Y_COUPLING rank 1, IMC_LINE rank 2), fully consistent with the NOISEHOUND Jan 2025 ranking. Two additional pointers emerge from EXCAVATor not present in the NOISEHOUND top 20: (a) V1:ENV_CEB_UPS_CURR_R_mean (CEB UPS current, rank 15) — first per-event evidence for a CEB electrical coupling, consistent with the logbook suspicion; (b) V1:SBE_SNEB_GEO_GRWE_raw_mean (Suspended North End Bench geophone W-E, rank 5) and GRNS (rank 20) — seismic/mechanical activity at the North End bench correlates with the glitches. Both channels have been added to the slow-sensor rate-correlation campaign (Step 3), together with the SWEB counterparts as controls. The EXCAVATor report is in Appendix C.

This document reports the NOISEHOUND investigation. Section 2 summarises prior knowledge from the Virgo logbook. Section 3 outlines the investigation roadmap (Steps 1–6). Section 4 gives the detailed results for each step. The full logbook history is in Appendix A; trigger CSV files in Appendix B; the EXCAVATor cross-check in Appendix C.

---

## 2 — Background: prior knowledge from the Virgo logbook

### Physical characteristics

Recurring broadband glitches in Hrec_hoft with SNR ~400, peak frequency 40–50 Hz, recurrence interval ~1600 ± 60 s (~26 min). The interval varies seasonally: ~23 min in winter, up to ~32 min in summer. By late 2025 the interval had drifted to ~30 min and the family is sometimes called "30-minute glitches". The waveform is broadband and approximately symmetric, consistent with a step-like strain discontinuity. A residual 50 Hz tail appeared from mini-ER December 2023 onwards, attributed to feedforward subtraction.

### Established conclusions

Best correlator: V1:INF_NI_BOTTOM_TE1 (NI tower bottom temperature), Pearson r = −0.72 with glitch rate (logbook #67414, direnzo, Aug 2025) — significantly stronger than any mirror temperature channel. The seasonal variation of the recurrence period tracks ambient temperature and this channel. Source location is suspected near the North Input (NI) tower; the CEB area and UPS mains line have also been suspected.

What is ruled out: NE mirror replacement, SWEB controls, SDB2_B1 channels (UPV analysis). V1:Sc_WI_FF50HZ_P_ERR is flagged 'Danger' in the ER16 channel safety study (#66923) — its correlation with the glitches is non-causal (feedforward reacting to the loud strain transient).

### Data quality impact

Glitches at SNR ~400 cause measurable drops in the BNS range when the PSD is estimated with a mean-based Welch method. Median-based estimators (e.g. gwistat) are robust against them.

See Appendix A for the full chronological logbook history.

---

## 3 — Investigation roadmap

### Step 1 — Per-event NOISEHOUND ranking on 1 Hz trend data  [completed]

Applied the NOISEHOUND ranking pipeline to a 3-hour epoch on 2025-01-01 (GPS 1419724818–1419735618, 7 glitches, z-score 28–44). Extracted ±5 s windows from 9161 1 Hz trend channels; ranked by z-score, hit fraction, and cross-correlation lag. Result: LSC/DARM channels dominate (IMC line, PSTAB0, SIB1, lag +1.5–2 s). Top non-DARM channel Sc_WI_FF50HZ_P_ERR is non-causal (positive lag, 'Danger' flag). → Detailed results: Section 4.1.

### Step 2 — Causality analysis on the 3-hour epoch  [completed]

Three methods applied to the top 40 ranked channels: cross-correlation lag sign, Granger causality (VAR F-test), and Transfer Entropy (Schreiber 2000). All 40 channels show positive lag and TE direction ←hrec: strain drives auxiliary in every case. The ±5 s window is blind to the actual causal mechanism, which operates on thermal timescales. → Detailed results: Section 4.2.

### Step 3 — Glitch rate correlation with slow sensors over full O4b  [ongoing]

Pearson and Spearman correlation of the hourly glitch rate against 30 slow channels over the full O4b dataset (Apr 2023 – Jan 2026, 24 212 triggers, 4 824 one-hour bins). A cross-correlation lag scan over ±7 days detects any thermal lead time. Channels were selected from three sources: (a) the Virgo logbook history of this glitch family, (b) physical reasoning about which slow sensors could modulate the thermal state of the input mirrors, and (c) the EXCAVATor per-event ranking (Appendix C). → Detailed results: Section 4.3.

The 30 channels and their selection rationale:

| Channel (V1:...) | Description | Selection rationale / source |
|---|---|---|
| INF_NI_BOTTOM_TE1 | NI tower bottom temperature | Logbook #67414 (direnzo, Aug 2025): best known rate correlator, Pearson r=−0.72. Logbook #68210: remains best predictor after Sep 2025 CEB HVAC failure. |
| INF_WI_BOTTOM_TE1 | WI tower bottom temperature | WI counterpart; comparison with NI to localise the thermal source. |
| INF_NI_MIR_COIL_UL_TE | NI mirror coil upper-left temperature | INF thermometry close to NI suspension; tracks local heat load on mirror coil. |
| INF_WI_MIR_COIL_DR_TE | WI mirror coil down-right temperature | WI counterpart. |
| TCS_HWS_NI_TE1_mean | NI HWS bench thermistor 1 | Temperature of the NI Hartmann Wavefront Sensor optics bench; tracks thermal state of the NI TCS system. |
| TCS_HWS_NI_TE2_mean | NI HWS bench thermistor 2 | Second NI HWS thermistor for redundancy. |
| TCS_HWS_WI_TE1_mean | WI HWS bench thermistor 1 | WI counterpart. |
| TCS_HWS_WI_TE2_mean | WI HWS bench thermistor 2 | WI counterpart. |
| TCS_HWS_NE_TE1_mean | NE HWS bench thermistor (control) | North End bench, geographically remote from NI; serves as control channel. |
| ENV_TCS_CO2_NI_TE | NI CO2 bench ambient temperature | Ambient thermal environment of the NI CO2 laser bench. |
| TCS_NI_TE_CO2Laser | NI CO2 laser body temperature | Direct CO2 laser temperature; thermally coupled to output power stability. |
| ENV_TCS_CO2_WI_TE | WI CO2 bench ambient temperature | WI counterpart. |
| TCS_WI_TE_CO2Laser | WI CO2 laser body temperature | WI counterpart. |
| TCS_NI_CO2_PWRLAS_mean | NI CO2 laser output power | Power proxy for NI TCS actuation strength; complements the temperature channels. |
| LSC_Etalon_NI_RH_SET_mean | NI ring heater setpoint [W] | Ring heater is the primary slow thermal actuator on the NI input mirror; its setpoint tracks the required thermal correction. |
| LSC_Etalon_NI_RH_OUT_mean | NI ring heater output [W] | Actual delivered power; may differ from setpoint due to control saturation. |
| LSC_Etalon_NI_RH_IN_mean | NI ring heater input [W] | Input power monitor. |
| LSC_Etalon_NI_RH_ERR_mean | NI ring heater control error | Loop residual; tracks deviation of mirror temperature from thermal target. |
| LSC_Etalon_WI_RH_SET_mean | WI ring heater setpoint [W] | WI counterpart. |
| LSC_Etalon_WI_RH_OUT_mean | WI ring heater output [W] | WI counterpart. |
| LSC_Etalon_WI_RH_ERR_mean | WI ring heater control error | WI counterpart. |
| ENV_NEB_UPS_VOLT_R_mean | NEB UPS voltage phase R [V] | Mains voltage at NEB; logbook suspicion of UPS mains coupling to glitch rate. |
| ENV_CEB_UPS_VOLT_R_mean | CEB UPS voltage phase R [V] | CEB area mains voltage; CEB area suspected as source location in several logbook entries. |
| ENV_WEB_UPS_VOLT_R_mean | WEB UPS voltage phase R [V] | WEB mains voltage; control channel. |
| ENV_MCB_IPS_CURR_T_mean | MCB IPS current phase T [A] | Main control building current; proxy for overall mains electrical load. |
| ENV_CEB_UPS_CURR_R_mean | CEB UPS current phase R [A] | EXCAVATor per-event rank 15: first quantitative evidence of CEB electrical coupling. Added after EXCAVATor cross-check (Appendix C). |
| SBE_SNEB_GEO_GRWE_raw_mean | SNEB geophone W-E [counts] | EXCAVATor per-event rank 5: Suspended North End Bench geophone W-E; seismic/mechanical activity at SNEB correlates per-event with glitches. Added after EXCAVATor cross-check. |
| SBE_SNEB_GEO_GRNS_raw_mean | SNEB geophone N-S [counts] | EXCAVATor per-event rank 20: N-S component of SNEB seismic sensor. |
| SBE_SWEB_GEO_GRWE_raw_mean | SWEB geophone W-E [counts] | West End Bench counterpart; control channel (remote from NI). |
| SBE_SWEB_GEO_GRNS_raw_mean | SWEB geophone N-S [counts] | SWEB N-S control channel. |

### Step 4 — Lag refinement with high-frequency frames  [planned]

Once Step 3 identifies the leading channel(s) and approximate lag, the lag will be refined to sub-second precision using raw Virgo frames from HPSS. → Results: Section 4.4.

### Step 5 — Convergent Cross Mapping  [planned, if needed]

If Granger/TE remain inconclusive due to nonlinearity, CCM (Sugihara et al. 2012, Science 338:496) will be applied. → Results: Section 4.5.

### Step 6 — Control epoch: 2026 glitch-free period  [planned]

Re-run the rate-correlation pipeline on Jan–Mar 2026 (glitch-free period) to cross-validate the thermal correlation and characterise the disappearance. → Results: Section 4.6.

---

## 4 — Detailed results

### 4.1 — Step 1: per-event NOISEHOUND ranking

#### Run parameters

- Epoch: 2025-01-01, GPS 1419724818–1419735618 (3 hours)
- Glitches: 7, z-score range 28–44
- Channels ranked: 9161 (1 Hz trend)
- Window: ±5 s around each glitch peak

#### Top channels

![Overlay of top auxiliary channels around 25-minute glitch events](overlay_aux_channels.png)

| Channel | rank_score | median_z | hit_frac | lag (s) |
|---|---|---|---|---|
| V1:LSC_DARM_IMC_LINE_mag_100Hz_mean | 175.1 | 207.8 | 1.0 | +1.5 |
| V1:LSC_DARM_PSTAB0_COUPLING_100Hz_mean | 44.5 | 82.8 | 1.0 | +1.5 |
| V1:LSC_DARM_IMC_LINE_Q_100Hz_mean | 39.5 | 40.8 | 0.875 | +2.0 |
| V1:LSC_DARM_PSTAB0_I_FS_mean | 37.9 | 53.1 | 1.0 | +2.0 |
| V1:LSC_DARM_PSTAB0_I_100Hz_mean | 37.7 | 51.7 | 1.0 | +2.0 |
| V1:LSC_DARM_SIB1_LINE_mag_100Hz_mean | 31.8 | 44.9 | 1.0 | +1.5 |
| V1:Sc_WI_FF50HZ_P_ERR_mean | 22.8 | 36.2 | 1.0 | +1.0 |
| V1:LSC_DARM_SIB1_LINE_Q_100Hz_mean | 21.3 | 31.3 | 0.875 | +2.0 |
| V1:LSC_DARM_IMC_LINE_I_100Hz_mean | 20.6 | 27.3 | 1.0 | +1.0 |
| V1:Sc_WI_FF50HZ_PHASE_mean | 14.7 | 16.6 | 0.875 | +2.0 |
| V1:LSC_DARM_PSTAB0_Q_100Hz_mean | 11.1 | 32.3 | 1.0 | +2.0 |
| V1:LSC_DARM_PSTAB0_Q_FS_mean | 10.8 | 31.1 | 1.0 | +2.0 |
| V1:ASC_SR_TY_DCP_mag_B1_mag_10Hz_FS_mean | 8.3 | 8.2 | 0.625 | +5.0 |
| V1:LSC_DARM_BPC_TY_COUPLING_100Hz_mean | 8.1 | 18.9 | 0.875 | +1.5 |
| V1:SQB1_Cam1_FitWaistY_mean | 8.0 | 126.4 | 0.875 | +2.0 |
| V1:LSC_DARM_BPC_Y_COUPLING_100Hz_mean | 7.4 | 27.1 | 1.0 | +1.5 |
| V1:Sc_WI_FF50HZ_G_ERR_mean | 6.0 | 17.7 | 0.625 | +1.0 |
| V1:LSC_DARM_SIB1_LINE_I_100Hz_mean | 5.2 | 16.0 | 0.625 | +1.5 |
| V1:ASC_DIFFp_TY_DCP_mag_B1_mag_10Hz_FS_mean | 4.8 | 7.8 | 0.625 | +5.0 |
| V1:LSC_DARM_BPC_TX_COUPLING_100Hz_mean | 3.6 | 11.1 | 0.875 | +1.5 |

#### Top subsystems

| Subsystem | max rank_score |
|---|---|
| LSC | 175.1 |
| ASC | 8.3 |
| SQB1 (squeezing bench camera) | 8.0 |
| FCEB | 2.2 |
| SDB1 / SDB2 | 1.6 |
| SPRB | 1.5 |
| SNEB | 1.5 |
| EQB1 | 0.9 |
| PSL | 0.8 |

#### Key findings

All top-ranked channels belong to LSC/DARM control loops, with positive lags of +1.5–2 s (strain leads auxiliary). The top non-DARM channel (Sc_WI_FF50HZ_P_ERR) carries a 'Danger' flag and is confirmed non-causal (feedforward reacting to the loud strain transient). No slow environmental or thermal channel ranks in the top 40 — consistent with a causal mechanism operating on timescales much longer than the ±5 s ranking window.

### 4.2 — Step 2: causality analysis

All 40 top-ranked channels were tested with three methods over the same 3-hour epoch (GPS 1419724818–1419735618):

- **Cross-correlation lag**: positive lag for all 40 channels (strain leads auxiliary by 1–5 s). A negative lag would be required for any causal candidate; none found.

- **Granger causality** (VAR model, lag order selected by AIC, F-test): the direction hrec→aux is significant at p < 0.001 for all channels tested, including ASC_SR_TY (p = 0.000). The reverse direction (aux→hrec) is not significant.

- **Transfer entropy** (Schreiber 2000, histogram MI estimator, 5-bin): TE(hrec→aux) > TE(aux→hrec) for all 40 channels. The net information flow is from strain to auxiliary in every case. No channel shows a causal signature within the ±5 s event window.

**Overall conclusion**: the per-event ranking approach correctly identifies channels that couple to the glitch, but cannot identify its origin because the causal mechanism operates on timescales much longer than the ±5 s window.

### 4.3 — Step 3: glitch rate correlation  [results pending]

SLURM array job running on CC-IN2P3. This section will be updated with the correlation table, scatter plots, and lag-scan results once the job completes. Expected key output: Pearson r and optimal lag for INF_NI_BOTTOM_TE1 and the 29 other slow channels over the full O4b baseline (Apr 2023 – Jan 2026, 24 212 triggers, 4 824 one-hour bins).

### 4.4 — Step 4: lag refinement  [pending Step 3]

To be completed after Step 3 identifies the dominant channel and approximate thermal lead time.

### 4.5 — Step 5: Convergent Cross Mapping  [pending Step 3]

To be performed if Granger/TE results in Step 3 are inconclusive.

### 4.6 — Step 6: 2026 control epoch  [planned]

To be run after Step 3 completes, using the same pipeline on Jan–Mar 2026 data.

---

## Appendix A — Virgo logbook entries (25-minute glitches)

Entries retrieved from logbook.virgo-gw.eu, search keyword '25-minute', listed chronologically.

| Entry | Date | Author(s) | Summary |
|---|---|---|---|
| #59763 | 12 Apr 2023 | direnzo | First report of periodic loud glitches in Hrec_hoft and LSC_DARM during Apr 7–11 long locks. Omicron characterisation: median spacing 28m 32s (range 26m 40s – 29m 20s), peak frequency ~47.6 Hz. Brute-force correlation with 10k rms trend channels inconclusive; best Pearson r=38% on LSC_DARM_PSTAB0_COUPLING_100Hz_rms. |
| #59766 | 12 Apr 2023 | andrew.lundgren | Notes similarity to LLO chiller glitch (electrical transient at power line frequency from AC chiller cycling on/off), consistent with a ~50 Hz short glitch. |
| #59791 | 14 Apr 2023 | direnzo | Histogram of trigger spacings; brute-force correlation with rms and derivative trend channels. Best correlator: LSC_DARM_PSTAB0_COUPLING_100Hz_rms (Pearson r=38%), inconclusive. Top-100 channel list attached to logbook entry. |
| #59826 | 17 Apr 2023 | robinet | Omicron running online on LSC_DARM; glitches clearly visible in VIM. Peak frequency ~70 Hz, jumped suddenly to 85 Hz around 05:00 UTC — possible ON/OFF switch event. |
| #61837 | 2023 (pre-O4b) | — | Confirmed presence before run start; interval ~26 min |
| #66292 | Mar 2025 | Paoletti | Omicron trigger CSV files published; SNR∈[100,500], freq∈[30,50] Hz, sep>15 min |
| #66628 | Apr 2025 | narnaud | Post-intervention DQ overview: 25-min glitches confirmed still present |
| #66923 | May 2025 | bersanetti | Sc_WI_FF50HZ_* flagged as glitchy/Danger in ER16 channel safety study; correlation likely non-causal |
| #66972 | Jun 2025 | direnzo | Despite other glitch families disappearing after maintenance, 25-min remain (SNR~400, 40-50 Hz) |
| #66975 | Jun 2025 | direnzo, mwas | Dec 2024 high-rate period survey; 25-min glitches are the constant background family |
| #67414 | Aug 2025 | direnzo | Pearson r=−0.72 between glitch rate and INF_NI_BOTTOM_TE1; stronger than mirror temperature |
| #67431 | Aug 2025 | direnzo | Full O4b glitch population survey; 25-min visible as continuous 40-50 Hz band throughout run |
| — | Aug 2025 | CEB team | Interval varies 23-32 min seasonally tracking ambient temperature; source suspected near NI tower, possibly UPS mains |
| #68110 | Nov 2025 | bersanetti | DeepExtractor ML waveform clustering on Jun 11 data; confirms multiple glitch families co-existing |
| #68210 | Nov 2025 | direnzo | After HVAC failure in CEB (Sep 2025): etalon disruption; NI tower bottom temp remains best predictor; now called ~30-min glitches |
| #68511 | 17 Jan 2026 | narnaud | BruCo check for 25-minute glitches disabled: those glitches seem to have disappeared. |

## Appendix B — Trigger CSV files (direnzo)

Omicron triggers for the 25-minute family published by M. Di Renzo on the Virgo scientists portal. Selection: SNR∈[100,500], frequency∈[30,50] Hz, inter-trigger separation >15 min. Access requires EGO SSO authentication. Files are mirrored in `NOISEHOUND/data/` on CCA (CC-IN2P3).

Base URL: https://scientists.virgo-gw.eu/DataAnalysis/DetCharDev/users/direnzo/glitches/25minute/

- **Extended catalogue** — 24 212 triggers, Apr 2023 – Jan 2026: `full_25min_glitches_ER16-O4b.csv`
- **O4b catalogue** — 19 933 triggers, Mar 2024 – Mar 2026: `25min_glitches_ER16-O4b.csv`
- **mini-ER Dec 2023** — 2 106 triggers, Dec 2023 – Feb 2024: `25min_glitches_mER_Dec23.csv`

## Appendix C — EXCAVATor cross-check report

EXCAVATor (Swinkels, Nikhef, 2013) is an independent Virgo tool that correlates glitch occurrence with the value of auxiliary channels using two statistics: gain of detection probability (Eff/DT-based) and the Kolmogorov-Smirnov test. It operates on 1 Hz trend data, similarly to NOISEHOUND, but uses a multi-glitch epoch rather than per-event windows.

Report URL: https://scientists.virgo-gw.eu/DataAnalysis/Excavator/test/half_hour_glitch/

### Run parameters

- Epoch: GPS 1366750818–1366786818 (2023-04-23, 10 hours)
- Triggers: 26 (Omicron, V1:Hrec_hoft_16384Hz, f<100 Hz, SNR>20)
- Channels analysed: 11 807 (1 Hz trend mean)
- Blacklisted channels: 68 370

### Top results — gain ranking (top 5 + CEB UPS)

| Rank | Channel | gain | Eff | DT |
|---|---|---|---|---|
| 1 | V1:LSC_DARM_BPC_Y_COUPLING_100Hz_mean | 12.58 | 0.923 | 0.016 |
| 2 | V1:LSC_DARM_IMC_LINE_mag_100Hz_mean | 11.91 | 0.923 | 0.043 |
| 3 | V1:LSC_DARM_BPC_X_COUPLING_100Hz_mean | 9.75 | 0.923 | 0.134 |
| 4 | V1:LSC_DARM_BPC_TY_COUPLING_100Hz_mean | 4.63 | 0.885 | 0.269 |
| 5 ★ | V1:SBE_SNEB_GEO_GRWE_raw_mean | 4.61 | 0.962 | 0.579 |
| 15 ★ | V1:ENV_CEB_UPS_CURR_R_mean | 2.41 | 0.962 | 0.696 |

★ Non-DARM channels of particular interest added to Step 3 channel list.

### Comparison with NOISEHOUND results

Both tools independently identify the same LSC/DARM family as the dominant per-event correlators across different epochs (Apr 2023 vs Jan 2025) and different algorithms. BPC_Y_COUPLING tops EXCAVATor while IMC_LINE tops NOISEHOUND — both are DARM channels from the same control loop; the ordering difference reflects the different ranking metrics (gain vs z-score).

Two channels appear in EXCAVATor but not in the NOISEHOUND top 20: ENV_CEB_UPS_CURR_R (CEB electrical supply — first per-event evidence for a mains coupling) and SBE_SNEB_GEO_GRWE/GRNS (North End Bench geophones — mechanical activity at the NE suspension). Both have been added to the Step 3 channel list (`rate_correlation_direct.py`), together with SWEB equivalents as arm-asymmetry controls.
