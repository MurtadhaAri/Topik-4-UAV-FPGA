# NUMERICAL_ERROR_ANALYSIS — Metodologi dan Hasil Analisis Galat Numerik

> Semua angka berasal dari `results/` (deterministik; `make analysis`). Dataset **sintetis**. Pipeline yang diukur adalah
> *study pipeline* System 4 (**PENDING INTEGRATION** dengan ES-EKF tim). Parameter noise/Q ber-label PROPOSED.

## 1. Definisi metrik (`src/analysis/s4_metrics.c`)
Untuk pasangan (referensi *rᵢ*, uji *tᵢ*), *i = 1…N* (N = jumlah sampel × komponen):

| Metrik | Rumus |
|---|---|
| RMSE | √( (1/N) Σ (rᵢ − tᵢ)² ) |
| Maximum Absolute Error | max \|rᵢ − tᵢ\| |
| Mean Absolute Error | (1/N) Σ \|rᵢ − tᵢ\| |
| Relative Error | (1/N′) Σ \|rᵢ − tᵢ\|/\|rᵢ\|, hanya \|rᵢ\| > 10⁻⁹ |
| Quantization Error | galat *round-trip* float→fixed→float pada state referensi (kind `quantization_roundtrip`) |
| Overflow Count | nilai pra-clip tidak muat int32 *container* |
| Saturation Count | hasil di-clip ke word target (16/24/32 bit) |
| Underflow Count | operan non-nol menghasilkan 0 pasti (mis. inkremen posisi < ½ LSB) |
| Sudut attitude | 2·atan2(\|r_v\|, \|r_w\|) dari q_ref⁻¹⊗q_uji → **q dan −q setara** (galat 0) |
| Komponen quaternion | dihitung setelah menyamakan hemisfer tanda (`s4_quat_sign`) |

Kelompok dievaluasi terpisah: position, velocity, attitude, accel bias, gyro bias, quaternion, U (Uₙ off-diag), D (per grup state), P = UDUᵀ
(diag/off-diag). U dan D dibandingkan di domain **ternormalisasi** (U tunggal sering *ill-conditioned*: Uᵢⱼ = Pᵢⱼ/Dⱼ besar bila Dⱼ kecil;
P = UDUᵀ lebih bermakna).

## 2. Metodologi
* Model A (float32) dan Model B (fixed) dijalankan *lock-step* pada frame `sensor_frame_t` yang sama; perbandingan tiap 0,1 s (600 titik/skenario).
* Skenario (`src/analysis/s4_dataset.c`): **static** (IMU+GNSS ideal, tanpa noise; kasus sanity ICD), **constvel** (v=[10,5,−0,5] m/s), **dynamic**
  (rate gyro & akselerasi sinusoidal), **highrange** (error awal ≈95 % batas ICD: 9,5 m, 9,5 m/s, ≈0,485 rad, bias 1,9 m/s²/0,19 rad/s; filter mulai dari nol).
  constvel/dynamic memakai *warm start* kecepatan (GNSS-style); highrange sengaja *cold start*. Noise GNSS/baro mengikuti R ICD; noise IMU PROPOSED
  (a: 0,05 m/s², ω: 0,003 rad/s). IMU 500 Hz, GNSS 5 Hz, baro 20 Hz. Prediksi kovariansi tiap 10 sampel.
* Test 5 (quantization stress) dan Test 6 (arithmetic unit test) berdiri sendiri di `tests/`.

## 3. Akurasi *reference* itu sendiri (penting untuk menafsirkan hasil)
| Pemeriksaan | Hasil |
|---|---|
| Float32 UDU vs *dense double* P-EKF independen (galat relatif terburuk diag P, akhir 60 s) | static 7,8·10⁻⁵; constvel 7,1·10⁻⁵; dynamic 5,3·10⁻⁵; highrange 7,6·10⁻⁵ |
| Sensitivitas float32 terhadap perturbasi input ≈ 1 ε float32 (RMSE) — posisi / kecepatan / attitude | constvel 3,5·10⁻⁶ m / 5,1·10⁻⁷ / 1,2·10⁻⁶ rad; dynamic 7,3·10⁻⁶ / 1,9·10⁻⁶ / 1,2·10⁻⁶; highrange 1,5·10⁻⁵ / 7,3·10⁻⁶ / 1,1·10⁻⁶ |
| Sensitivitas P_norm diag terhadap perturbasi yang sama | 5·10⁻⁹…4·10⁻⁸ (diukur selama pengembangan, constvel/dynamic/highrange) |
Artinya selisih fixed-vs-float pada level ≈ 10⁻⁴ (relatif, kovariansi) tidak dapat dihubungkan ke fixed-point; pada nominal state *floor*-nya jauh lebih rendah (10⁻⁶–10⁻⁵).
Verifikasi Thornton/Bierman float juga dijalankan di `test_models` (< 10⁻³ relatif; terukur 3,9·10⁻⁵ static / 2,4·10⁻⁵ dynamic pada 30 s).

## 4. Hasil utama
Tabel RMSE W16/W24/W32 dan ablation per kelompok ada di `docs/Q_FORMAT_ANALYSIS.md` §5 (sumber: `summary_rmse.csv`, `ablation.csv`).
Ringkasan terukur (60 s, fixed vs float):
* **W16**: RMSE posisi 6,4–11,6 cm, attitude 2–28·10⁻² rad → tidak layak.
* **W24**: posisi 3,2·10⁻³ (constvel/dynamic) – 9,3·10⁻³ m (highrange).
* **W32**: posisi 1,0·10⁻³ / 2,1·10⁻⁴ / 1,4·10⁻³ m (constvel/dynamic/highrange).
* Galat terhadap **truth**: float, W24, W32 praktis identik (constvel posisi 0,109 m; dynamic 0,134 m; highrange ≈ 0,40 m) → pada skenario realistis
  kuantisasi ≥ 24 bit tidak menurunkan akurasi estimasi; W16 menurunkan (constvel 0,136 m, dynamic 0,154 m).
* Skenario **static** (ICD sanity): drift posisi 0,000000 mm dan Dₙ > 0 pada W16/W24/W32, 0 klem positif-definit (`static_criteria.csv`).
  *Ini pemeriksaan pipeline studi, bukan kriteria penerimaan ES-EKF penuh.*

## 5. Precision degradation, overflow, underflow
| Skenario | Saturasi | Overflow | Underflow (W16 / W24 / W32) | Klem D≤0 |
|---|---|---|---|---|
| static, constvel | 0 | 0 | 0 / 0 / 0 | 0 |
| dynamic | 0 | 0 | 72 / 96 / 70 (semua di inkremen posisi nominal) | 0 |
| highrange | 0 | 0 | 1 / 1 / 1 | 0 |
* Dengan dx Q4.x, U Q7.x, D Q1.x **tidak ada saturasi/overflow** pada 4 skenario. Dengan format sebelumnya saturasi *terdeteksi* (lihat Q_FORMAT_ANALYSIS §6, §4) — penghitung berfungsi sebagai alat desain.
* **Underflow** ada pada integrasi posisi nominal (Q16.16 = 15,3 µm): `pos_integration.csv`

| Format p | v = 0,01 m/s | v = 0,1 m/s | v = 1 m/s |
|---|---|---|---|
| Q16.16 | galat 23,7 % (0,142 m dari 0,600 m) | 0,82 % | 0,055 % |
| Q12.20 | 0,14 % | 0,14 % | 0,0073 % |
| Q8.24 | 0,14 % | 0,013 % | 0,0013 % |
  Penyebab: pembulatan inkremen tetap yang konsisten bias (mis. 20,97 LSB → 21). Bukan acak, sehingga *berakumulasi*.
* **Normalisasi quaternion** (`quat_norm.csv`): max \|\|q\|−1\| dengan Newton inverse-sqrt = 0 (static), ≈ 1·10⁻⁹ (constvel, dynamic), 5,6·10⁻⁴/4,3·10⁻⁸/4,1·10⁻⁸ (highrange, W16/W24/W32). Dengan aturan Taylor orde-1 ICD (`-DS4_QNORM_TAYLOR`) nilai highrange 0,29/0,44/0,44 (≈ 3·10⁻⁶ di dynamic) — Taylor hanya valid bila \|q\| ≈ 1; injeksi \|δθ\| > ~0,5 rad melanggarnya.

## 6. Modul aritmetika
**CORDIC** (30 iterasi, 40 001 sudut di [−π,π] di `test_cordic`; tabel `cordic_error.csv` memakai 20 001 sudut): RMS 3,2·10⁻⁹, maks. 1,4·10⁻⁸ (≈ 14,7 LSB Q2.30), maks. \|sin²+cos²−1\| 2,8·10⁻⁸.

| Iterasi | 8 | 12 | 16 | 20 | 24 | 28 | 30 |
|---|---|---|---|---|---|---|---|
| maks. galat sin | 7,8·10⁻³ | 4,9·10⁻⁴ | 3,0·10⁻⁵ | 1,9·10⁻⁶ | 1,3·10⁻⁷ | 1,7·10⁻⁸ | 1,4·10⁻⁸ |
Galat turun ≈ 4× per 2 iterasi (2 bit/iterasi) hingga lantai pembulatan ≈ 15 LSB (galat akumulasi ≈ N·½ LSB).

**Divider** (`divider_compare.csv`, 100 000 pasangan acak Q16.16):

| Mode | Beda vs baseline | Galat maks. vs eksak | Iterasi rata-rata / maks. | Estimasi siklus rata-rata | x/0 |
|---|---|---|---|---|---|
| baseline (int64 `/`) | 0 | 1,00 LSB | — (perangkat lunak) | tidak berlaku | saturasi ±max + flag |
| restoring | **0** (bit-exact) | 1,00 LSB | 38,0 / 39 | 38,0 | idem |
| non-restoring | **0** (bit-exact) | 1,00 LSB | 38,0 / 39 | 39,0 (+1 koreksi) | idem |
Estimasi siklus = jumlah iterasi radix-2 (1 bit/siklus), **bukan** hasil sintesis; frekuensi clock/latensi aktual PENDING. Iterasi = lebar bit dividen
(dividen = |a|<<frac); perangkat keras memakai batas tetap lebar maksimum.

## 7. Menjalankan ulang / memverifikasi
`./scripts/run_all_tests.sh --analysis` → semua tabel di `results/reports/analysis_report.md` dibangkitkan ulang. Pipeline divider-agnostik: `test_models`
memeriksa bahwa baseline/restoring/non-restoring menghasilkan RMSE **identik bit-demi-bit** pada skenario dynamic.

## 8. Keterbatasan (wajib dibaca)
1. Data sintetis, satu *seed*, tanpa magnetometer, tanpa gating outlier; bukan ES-EKF tim.
2. Float32 reference memiliki galat sendiri 5·10⁻⁵–8·10⁻⁵ relatif pada P.
3. Beberapa perilaku **belum dijelaskan tuntas**: W32 tidak lebih baik dari W24 pada attitude `dynamic` (8,05·10⁻⁴ vs 5,11·10⁻⁴ rad).
4. Tabel RMSE highrange memakai normalisasi Newton; angka lama berbasis Taylor (attitude W32 4,9·10⁻³ rad) hanya tersisa sebagai catatan.
5. Divisi dilakukan per pasangan (i,j) (bukan satu reciprocal + perkalian); implementasi hardware dengan reciprocal akan mengubah pembulatan — PENDING.
6. Lebar register kerja (int32 Q8.24/Q3.29, akumulator int64) belum dipetakan ke DSP48 — PENDING.
