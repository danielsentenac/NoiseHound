# 25-minute glitch — Investigation notes
## Interval 1426891202–1426897935 (2025-03-24 22:39 – 2025-03-25 00:31 UTC)

---

## Observation

The glitch is clearly visible in **V1:Hrec_hoft** and **V1:LSC_DARM_CORR** as a sharp
spike at 40 Hz (fstart–fend: 37.2–44.0 Hz, Q≈5.5, SNR≈490).
No clear response is seen in the AUX channels corresponding to mirror correction signals.

---

## Hypotheses (in order of likelihood)

### 1. Glitch injected downstream of the control loops *(most likely)*

The AUX correction channels shown are mostly **angular** (ASC) with control bandwidth
of only ~1–5 Hz. A 40 Hz transient is completely invisible to them — they simply
cannot react. The longitudinal DARM loop has higher bandwidth and DARM_CORR *does*
show a spike, consistent with this picture.

**Next step:** confirm the ASC bandwidth and check whether any longitudinal correction
channel (other than DARM_CORR) shows a response.

---

### 2. Glitch is a sensing / electronic artifact

If the noise enters *after* the error point — in the photodetection chain,
demodulation, or signal reconstruction (calibration pipeline) — the control system
never sees it and never sends corrections to the mirrors. The ITF is physically
undisturbed but Hrec is corrupted. This is the hallmark of a **fake strain injection
by electronics**.

**Next step:** compare the raw DARM error signal (before calibration) with the
reconstructed Hrec. If the glitch is absent in the raw error signal but present in
Hrec, the injection point is in the calibration/reconstruction pipeline.

---

### 3. Pendulum mechanical filtering

At 40 Hz the mirror pendulum response goes as ~(f₀/f)² ≈ (0.6/40)² ≈ 10⁻⁴.
Even if the loop pushes, the actual displacement is suppressed by ~10,000×.
However, the correction *force* signal should still appear in the actuator channels
even if the physical motion is tiny — so this alone does not explain the silence
in the correction channels.

**Status:** unlikely as the sole explanation; may contribute.

---

## Conclusion (current)

The most consistent interpretation is that this glitch is **not a physical mirror
displacement** but an artifact injected into the sensing or signal reconstruction
chain. The ITF is intact — the glitch fakes a strain event. This rules out seismic
or acoustic coupling to the mirrors as the primary source and is important context
for veto/classification.
