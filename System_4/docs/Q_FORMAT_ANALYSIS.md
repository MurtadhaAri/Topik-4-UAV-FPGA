# Q_FORMAT_ANALYSIS — Analisis Bit-Width, Range, dan Rekomendasi Q-format

> **Status dokumen.** Semua angka di bawah berasal dari `results/csv/*.csv` dan `results/reports/analysis_report.md`
> (dibangkitkan oleh `build/bin/run_analysis`, deterministik, dataset **sintetis**, 60 s @500 Hz). Parameter yang bukan bagian
> kesepakatan tim ditandai **PROPOSED**. Pengukuran dilakukan pada *study pipeline* System 4 (bukan ES-EKF tim) →
> **PENDING INTEGRATION**. Reproduksi: `make analysis`.

## 1. Notasi dan rumus
`Qm.n`: *m* bit integer **termasuk tanda**, *n* bit pecahan, word *W = m+n*. Untuk two's complement:

| Besaran | Rumus |
|---|---|
| Resolusi = LSB | 2⁻ⁿ |
| Maks. galat rounding | ½·2⁻ⁿ |
| Maks. galat truncation (floor) | 2⁻ⁿ |
| Nilai min | −2^(m−1) |
| Nilai maks | 2^(m−1) − 2⁻ⁿ |

Perbandingan format (`results/csv/qformat_table.csv`; sign bits = 1 untuk semua):

| Format | W | Integer bits (tanpa tanda) | Frac | Resolusi | Min | Max | Galat rounding maks |
|---|---|---|---|---|---|---|---|
| Q1.15 | 16 | 0 | 15 | 3,05·10⁻⁵ | −1 | 0,999969 | 1,53·10⁻⁵ |
| Q2.14 | 16 | 1 | 14 | 6,10·10⁻⁵ | −2 | 1,99994 | 3,05·10⁻⁵ |
| Q7.9 | 16 | 6 | 9 | 1,95·10⁻³ | −64 | 63,998 | 9,77·10⁻⁴ |
| Q1.23 | 24 | 0 | 23 | 1,19·10⁻⁷ | −1 | 0,99999988 | 5,96·10⁻⁸ |
| Q2.22 | 24 | 1 | 22 | 2,38·10⁻⁷ | −2 | 1,9999998 | 1,19·10⁻⁷ |
| Q7.17 | 24 | 6 | 17 | 7,63·10⁻⁶ | −64 | 63,99999 | 3,81·10⁻⁶ |
| Q2.30 | 32 | 1 | 30 | 9,31·10⁻¹⁰ | −2 | 1,999999999 | 4,66·10⁻¹⁰ |
| Q1.31 | 32 | 0 | 31 | 4,66·10⁻¹⁰ | −1 | 0,9999999995 | 2,33·10⁻¹⁰ |
| Q16.16 | 32 | 15 | 16 | 1,53·10⁻⁵ | −32768 | 32767,99998 | 7,63·10⁻⁶ |

## 2. Quaternion (Q1.15 vs Q1.23 vs Q2.30)
| Format | +1,0 | Galat pada q = [1,0,0,0] | Resolusi sudut (≈ 2·LSB) |
|---|---|---|---|
| Q1.15 | **saturasi** → 0,999969 | 3,05·10⁻⁵ | 6,1·10⁻⁵ rad (0,0035°) |
| Q1.23 | **saturasi** → 0,99999988 | 1,19·10⁻⁷ | 2,4·10⁻⁷ rad |
| Q2.30 | **tepat** (1<<30) | 0 | 1,9·10⁻⁹ rad |

Verifikasi: `test_quant_stress` (saturasi tanpa wrap-around), `test_adapter` (qw=+1 tepat di Q2.30).
**Rekomendasi: Q2.30** — sesuai ICD Role 4 (“integrasi q wajib 32-bit”), +1 tepat, dan memberi 1 bit headroom untuk
renormalisasi (|q|² sedikit di atas 1). Q1.15/Q1.23 hanya cukup bila quaternion disimpan *hanya* untuk keluaran, bukan
integrasi: resolusi 6,1·10⁻⁵ rad × 500 Hz akan menelan inkremen ½ωΔt < LSB (underflow).
*Catatan cakupan:* di pipeline, q selalu Q2.30; Q1.15/Q1.23 diuji sebagai representasi (bukan sebagai jalur integrasi penuh).

## 3. Normalisasi variabel (PROPOSED)
`x_norm = x_phys / s`, `x_phys = x_norm · s`. Skala dipilih **pangkat dua** terkecil ≥ batas ICD → normalisasi = *bit-shift*.
Biaya: range ternormalisasi tidak penuh [−1,1] untuk posisi/kecepatan (0,625) dan gyro bias (0,8).

| Variable | Physical Range (ICD / terukur) | Scale Factor | Normalized Range | Proposed Format |
|---|---|---|---|---|
| δp (error posisi) | ±10 m | 16 (2⁴) | ±0,625 | dx: **Q4.x** (lihat §5) |
| δv | ±10 m/s | 16 (2⁴) | ±0,625 | idem |
| δθ | ±0,5 rad | 0,5 (2⁻¹) | ±1,0 | idem |
| δb_a | ±2 m/s² | 2 (2¹) | ±1,0 | idem |
| δb_g | ±0,2 rad/s | 0,25 (2⁻²) | ±0,8 | idem |
| Quaternion q | [−1,1] (unit) | 1 (tidak dinormalisasi) | [−1,1] | **Q2.30** |
| p (nominal) | terukur maks. 1113 m (highrange); ICD: “ribuan m” | — (satuan fisik) | — | **Q16.16** (±32768 m) — lihat §6 |
| v (nominal) | terukur maks. 22,3 m/s | — | — | **Q8.24** (±128 m/s) |
| bias nominal | |b| ≪ 8 | — | — | **Q4.28** (±8) |
| U (Uₙ off-diag) | terukur **maks. 30,0** (static/constvel 60 s) | s_j/s_i | tak terbatas (lihat §4) | **Q7.x** (±64) |
| D (Dₙ) | terukur maks. 0,206; min (float) 7,2·10⁻⁶ | 1/s² | [0, 0,21] | **Q1.x** |
| ΔΦₙ (Φ−I ternormalisasi) | ≪ 1 (O(10⁻²) per langkah kov. 0,02 s) | s_j/s_i | (−1,1) | **Q1.x** |
| Register kerja W (Thornton) | maks. |Uₙ| ≈ 30 | — | **Q8.24** (int32) |
| α, g, k (Bierman) | ≪ 4 | — | — | **Q3.29** |
| IMU a, ω, dt | a ≤ 128 m/s², ω ≤ 128 rad/s | — | — | **Q8.24**, dt **Q2.30** |

Faktor skala *terpisah per kelompok* (bukan satu skala global) — Position/Velocity: 16, Attitude: 0,5, Accel bias: 2,
Gyro bias: 0,25, Covariance: turunan (1/s², s_j/s_i), Quaternion: 1.

## 4. Covariance U dan D — analisis rentang berbasis pengujian
`results/csv/range_analysis.csv` (float reference, sampling 0,1 s):

| Skenario | maks. \|Uₙ\| off-diag | maks. Dₙ | puncak \|dxₙ\| sebelum injeksi | bit integer (incl. tanda) U / D / dx |
|---|---|---|---|---|
| static | 29,99 | 0,206 | 0 | 6 / 1 / 1 |
| constvel | 29,45 | 0,200 | 0,124 | 6 / 1 / 1 |
| dynamic | 8,02 | 0,200 | 0,188 | 5 / 1 / 1 |
| highrange | 6,55 | 0,200 | **3,52** | 4 / 1 / **3** |

Interpretasi matematis:
* **D** adalah *variansi* (satuan², bukan satuan state). Setelah normalisasi Dₙ = D/s² ≤ 0,21 → Q1.x cukup; yang menentukan
  adalah **LSB**: Dₙ minimum 7,2·10⁻⁶ (float) → LSB harus ≪ 7·10⁻⁶ agar galat relatif ≤ 1 %: Q1.15 (3,05·10⁻⁵) *lebih besar dari
  D minimum* (underflow), Q1.23 (1,2·10⁻⁷, 1,7 %) marginal, Q1.31 (4,7·10⁻¹⁰) lapang.
* **U** *tidak terbatas*: yaw tak teramati (tanpa magnetometer) sehingga Uₙ(δθ_z, δb_gz) ≈ t·s_bg/s_θ = 0,5·t → ≈ 30 pada 60 s. Q7.x (±64)
  menampung ≈ 128 s tanpa koreksi yaw. **Magnetometer (PENDING) akan membatasi nilai ini** — jika tersedia, Q5.x cukup (dynamic: 8,0).
  Diagonal U = 1,0 → butuh ≥ 2 bit integer.
* Catatan koreksi: awalnya D diberi 2 bit integer karena P₈₈ (yaw) tumbuh ke ≈ 1,6; analisis rentang menunjukkan itu berasal dari
  Uₙ²·D (U ≈ 30), bukan D. D dikembalikan ke Q1.x (satu bit pecahan lebih banyak). Penyebab asli galat besar saat itu adalah saturasi U di ±16.
* **Process noise di bawah LSB (perhitungan dari parameter PROPOSED, bukan simulasi):** per langkah kovariansi (Δt_c = 0,02 s)
  Qd_n(δb_g) = 10⁻¹⁰·0,02/0,0625 = 3,2·10⁻¹¹ < LSB Q1.31 (4,7·10⁻¹⁰); Qd_n(δb_a) = 5·10⁻⁹ < LSB Q1.23. Pada pipeline, noise
  ini dijumlahkan di akumulator int64 tetapi hilang saat D dibulatkan ke storage. Mitigasi (belum diimplementasikan):
  stochastic rounding / carry residual / D berformat floating-blok.

## 5. Hasil pengukuran: 16 vs 24 vs 32 bit
Konfigurasi (storage word length jalur error; jalur nominal identik: p Q16.16, v Q8.24, q Q2.30, bias Q4.28):

| cfg | dx | U | D | Φ |
|---|---|---|---|---|
| W16 | Q4.12 | Q7.9 | Q1.15 | Q1.15 |
| W24 | Q4.20 | Q7.17 | Q1.23 | Q1.23 |
| W32 | Q4.28 | Q7.25 | Q1.31 | Q1.31 |

RMSE fixed vs float32 (60 s; `summary_rmse.csv`, kind `fixed_vs_float`):

| Skenario | Besaran | W16 | W24 | W32 |
|---|---|---|---|---|
| constvel | posisi [m] | 6,41·10⁻² | 3,17·10⁻³ | 1,01·10⁻³ |
| | kecepatan [m/s] | 2,78·10⁻² | 2,06·10⁻³ | 1,18·10⁻⁴ |
| | attitude [rad] | 3,48·10⁻² | 1,56·10⁻³ | 3,14·10⁻⁵ |
| | P_norm diag | 2,15·10⁻² | 4,35·10⁻⁴ | 1,48·10⁻⁵ |
| dynamic | posisi [m] | 6,45·10⁻² | 3,17·10⁻³ | 2,12·10⁻⁴ |
| | kecepatan [m/s] | 2,90·10⁻² | 2,19·10⁻³ | 4,16·10⁻⁴ |
| | attitude [rad] | 1,95·10⁻² | 5,11·10⁻⁴ | 8,05·10⁻⁴ |
| | P_norm diag | 6,18·10⁻³ | 2,73·10⁻⁵ | 9,58·10⁻⁵ |
| highrange | posisi [m] | 1,17e-1 | 9,57e-3 | 4,78e-4 |
| | kecepatan [m/s] | 1,57e-1 | 9,28e-3 | 2,39e-4 |
| | attitude [rad] | 2,48e-1 | 6,46e-3 | 4,24e-4 |
| | P_norm diag | 3,64e-3 | 2,36e-5 | 7,47e-7 |

Bacaan yang *sahih*: W16 jauh lebih buruk (10–100×) di semua skenario; W24→W32 memperbaiki posisi 3–20×.
Bacaan yang **tidak sahih**: urutan W24 vs W32 pada `dynamic` untuk *attitude* dan *P_norm diag* — selisihnya berada pada
level galat float32 reference sendiri (P_ii float32 vs dense double: 5,3·10⁻⁵…7,8·10⁻⁵ relatif). Dan pada `dynamic` W32 attitude (8,05·10⁻⁴) tidak
lebih baik daripada W24 (5,11·10⁻⁴) — penyebab **tidak teridentifikasi tuntas**; hipotesis “presisi register kerja W” sudah diuji (WFRAC 24→26)
dan **tidak terbukti**. Galat terhadap *truth* identik antara float dan W24/W32 (mis. dynamic posisi 0,134 m untuk keduanya; W16 0,154 m).

### Ablation per kelompok variabel (60 s; `ablation.csv`) — posisi [m] / attitude [rad]
| Varian (dx/U/D/Φ) | constvel | dynamic | highrange |
|---|---|---|---|
| W16 base (16/16/16/16) | 6,40e-2 / 3,44e-2 | 6,45e-2 / 1,95e-2 | 1,17e-1 / 2,48e-1 |
| W24 base (24/24/24/24) | 3,17e-3 / 1,56e-3 | 3,18e-3 / 5,09e-4 | 9,57e-3 / 6,46e-3 |
| W32 base (32/32/32/32) | 1,01e-3 / 3,74e-5 | 2,12e-4 / 8,05e-4 | 4,78e-4 / 4,24e-4 |
| W24 + D32 (24/24/32/24) | 1,01e-3 / 2,01e-3 | 2,20e-4 / 7,52e-4 | 5,16e-4 / 8,72e-4 |
| W24 + D32 + U32 (24/32/32/24) | 1,01e-3 / 1,26e-4 | 2,20e-4 / 8,01e-4 | 4,88e-4 / 4,49e-4 |
| W32 with D24 (32/32/24/32) | 3,17e-3 / 5,11e-4 | 3,18e-3 / 4,58e-4 | 9,64e-3 / 6,13e-3 |
| W32 with U24 (32/24/32/32) | 1,00e-3 / 1,91e-3 | 2,10e-4 / 7,57e-4 | 5,21e-4 / 9,07e-4 |
| W32 with dx24 (24/32/32/32) | 1,01e-3 / 6,71e-5 | 2,19e-4 / 8,12e-4 | 4,85e-4 / 4,21e-4 |
| W32 with Phi24 (32/32/32/24) | 1,01e-3 / 1,48e-4 | 2,13e-4 / 8,05e-4 | 4,81e-4 / 4,27e-4 |
| W16 + D32 (16/16/32/16) | 5,46e-3 / 4,31e-2 | 8,27e-3 / 1,95e-2 | 1,10e-1 / 8,07e-2 |
| mixed-A dx16 U24 D32 Phi16 (16/24/32/16) | 4,82e-3 / 2,05e-2 | 4,32e-3 / 3,79e-3 | 6,93e-3 / 3,25e-3 |
| mixed-B dx24 U24 D32 Phi24 (24/24/32/24) | 1,01e-3 / 2,01e-3 | 2,20e-4 / 7,52e-4 | 5,16e-4 / 8,72e-4 |

**Temuan utama: lebar word D adalah tuas akurasi dominan** (posisi 3–15× lebih baik; dikunci di `test_models`: “D32 pada W24 ≥ 2× lebih baik”,
terukur 7,3×/7,6×/5,5× pada 20 s). dx dan Φ hampir tidak berpengaruh bila diturunkan ke 24 bit; U memengaruhi terutama *attitude*.

## 6. Rekomendasi (berbasis data di atas; keputusan akhir ada pada tim)
| Kelompok | Rekomendasi | Alasan matematis / bukti |
|---|---|---|
| Quaternion (nominal) | **Q2.30** | +1 tepat; resolusi sudut 1,9·10⁻⁹ rad; ICD Role 4 |
| Posisi (nominal) | **Q16.16** hanya bila ditambah *residual carry*; bila tidak, Q12.20 (±2048 m) | Terukur: v=0,01 m/s → galat 23,7 % di Q16.16 vs 0,14 % di Q12.20/Q8.24. Q12.20 tak cukup untuk “ribuan meter” ICD. *Residual carry belum diimplementasikan/diuji* |
| Kecepatan (nominal) | **Q8.24** | maks. terukur 22,3 m/s; ±128 headroom |
| Bias nominal | **Q4.28** | ±8 ≫ batas ICD |
| dx (error state) | **24 bit, Q4.20** | puncak \|dxₙ\| 3,52 → 3 bit integer wajib. *Eksperimen selama pengembangan* dengan dx Q2.x (±2): highrange mencatat 4–5 saturasi dan galat attitude W24/W32 ≈ 0,20 rad (20–40× lebih besar daripada dengan Q4.x pada saat itu, ketika masih memakai normalisasi Taylor). Hasil itu tidak dapat direproduksi dengan konfigurasi default saat ini; untuk mengulang, ubah field `dx` di `k_cfg` (`src/model/ekf_fixed_adapter.c`) ke Q2.x. Diturunkan ke 24 bit tidak mengubah hasil (ablation) |
| U | **24 bit (Q7.17)** minimum; **32 bit (Q7.25)** bila attitude < 10⁻³ rad diperlukan | U dengan Q5.x (±16) saturasi pada skenario 60 s (terukur selama pengembangan: galat varians yaw −63 % pada semua word length; tidak lagi direproduksi oleh default Q7.x); U24 → attitude 2,0·10⁻³ rad (constvel) vs 3,1·10⁻⁵ rad (W32) |
| **D** | **32 bit (Q1.31)** | tuas dominan (ablation); Q1.15 di bawah D minimum; 16-bit tidak layak |
| Φ (ΔΦ) | **24 bit (Q1.23)** | W32→Φ24 tidak mengubah posisi; hanya selisih kecil attitude constvel (3,1·10⁻⁵ → 1,4·10⁻⁴ rad) |
| Kandidat keseluruhan | **mixed-B: dx Q4.20 / U Q7.17 / D Q1.31 / Φ Q1.23** (+ jalur nominal di atas) | posisi setara W32 di ketiga skenario; attitude lebih kasar di constvel (2,0·10⁻³ vs 3,1·10⁻⁵ rad). Pilih W32 penuh bila margin attitude dibutuhkan |
| **16-bit seragam** | **tidak direkomendasikan** | RMSE posisi 6–12 cm vs ≤ 1 cm pada W24/W32; Dₙ min < LSB |

## 7. Batasan yang harus dibaca bersama rekomendasi
1. Dataset sintetis, satu *seed*, parameter noise IMU dan Q **PROPOSED**; belum ada data terbang.
2. Pipeline studi bukan ES-EKF tim; tanpa magnetometer; kovariansi diprediksi tiap 10 sampel IMU (50 Hz) — **PROPOSED**.
3. **Normalisasi quaternion.** Aturan order-1 Taylor dari ICD gagal untuk injeksi besar: pada `highrange` \|\|q\|−1\| = 0,29/0,44/0,44 (W16/W24/W32). Model B kini memakai *inverse-sqrt Newton-Raphson* (4 iterasi, tanpa sqrt/pembagian): \|\|q\|−1\| ≤ 5,6·10⁻⁴ (W16) dan ≈ 4·10⁻⁸ (W24/W32) pada highrange, ≈ 10⁻⁹ di skenario lain. Aturan Taylor tersedia via `-DS4_QNORM_TAYLOR`. Angka highrange di §5 sudah memakai Newton (attitude W32: 4,9·10⁻³ → 4,2·10⁻⁴ rad).
4. Ambang penerimaan (mis. “1 mm posisi”) untuk skenario dinamis **bukan** bagian ICD; ambang di test (`test_models`) adalah PROPOSED.
