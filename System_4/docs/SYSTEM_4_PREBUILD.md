# SYSTEM_4_PREBUILD — Arsitektur & Prebuild Fixed-Point / Bit-Accurate Model

## 1. Ruang lingkup
**Dalam lingkup:** analisis kuantisasi, normalisasi, aritmetika fixed-point bit-accurate (saturasi/rounding terdefinisi), CORDIC, divider hardware-oriented,
adapter float↔fixed kompatibel ICD, dataset sintetis, pipeline RMSE, dokumentasi.
**Di luar lingkup:** algoritma ES-EKF resmi (Orang 1–3), tuning Q/R, dataset terbang nyata, ekspor `.mem`, RTL, sintesis, magnetometer, Git/remote.

### Hasil Phase 1 — inspeksi
* Dokumen ditemukan di `System_4/`: `Kesepakatan Teknis.docx`, `Kesepakatan_Teknis_ICD_ES-EKF.docx`; 6 paper PDF di folder induk (tidak diparse ulang di sini).
* **Repositori `micelele1/Topik-4-UAV-FPGA` tidak ada secara lokal dan `System_4` bukan repositori Git.** Tidak ada source/header anggota lain → tidak ada yang dimodifikasi.
  Tidak ada commit/push/branch yang dilakukan. Branch `System_4` belum dibuat.
* `ekf_types.h` hanya ada di dalam dokumen → ditranskripsi ke `third_party/icd/` (lihat INTERFACE_CONTRACT).
* Toolchain: gcc, make, cmake, python3 (tanpa numpy/matplotlib → plot dibuat SVG stdlib).

### Pemisahan sumber informasi
| Kategori | Isi |
|---|---|
| **Kesepakatan tim (ICD)** | struktur state & indeks, `U_IDX`, NED/FRD/Hamilton, g = 9,80665, batas error (10 m, 10 m/s, 0,5 rad, 2 m/s², 0,2 rad/s), D₀, R, aturan kode (statis, 1-D, tanpa inversi/sqrt di loop kov.), Bierman/Thornton, format `sensor_frame_t`, kriteria static (≤1 mm, D>0), saran Q16.16/Q2.30 untuk nominal |
| **Referensi tambahan** | PDF penelitian & playlist HLS (konsep hardware-friendly: loop tetap, aritmetika sederhana, validasi model sebelum RTL). *Tidak* diklaim membahas ES-EKF fixed-point |
| **PROPOSED oleh System 4** | skala pangkat-dua, seluruh Q-format error path, Q_c, noise IMU, laju sensor, desimasi kovariansi (10), format register kerja, iterasi CORDIC (30), ambang test, skenario sintetis |

## 2. Arsitektur
```
        sensor_frame_t (ICD, float)                      parameter PROPOSED (s4_params_t)
                 │                                                   │
   ┌─────────────┴──────────────┐                                    │
   │                            │ sensor_frame_to_fixed()            │
   ▼                            ▼                                    ▼
 MODEL A  float32 reference    MODEL B  fixed-point (int32/int64, saturating, rounding eksplisit)
 (ekf_s4_ref_float.c)          (ekf_s4_ref_fixed.c)
  nominal ─ ΦU ─ Thornton       nominal(Q16.16/Q8.24/Q2.30) ─ ΔΦU ─ Thornton(W=Q8.24) ─ Bierman(Q3.29) ─ inject/reset
  Bierman ─ inject/reset        dx/U/D/ΔΦ disimpan WL ∈ {16,24,32}, domain ternormalisasi
        │                            │              ▲ CORDIC (sin/cos)   ▲ divider (baseline/restoring/non-restoring)
        └──────────┬─────────────────┘              │                     │ (ekf_fixed_trig.c, ekf_fixed_divider.c; aritmetika: ekf_fixed_math.c)
                   ▼
        s4_compare.c  (lock-step, metrik: RMSE/MAE/max/rel, quantization round-trip, counters, vs truth, sensitivity floor)
                   ▼
        tools/run_analysis.c ──► results/csv/*.csv ──► scripts/plot_results.py ──► results/plots/*.svg
                   └────────────► results/reports/analysis_report.md
 s4_dense_ref.c : P-EKF dense *double* independen, hanya untuk menguji Model A
```
Keputusan arsitektur: (1) *storage vs compute*: dx/U/D/ΔΦ disimpan WL bit, register kerja 32-bit + akumulator int64 (pembulatan tunggal di akhir dot-product);
(2) ΔΦ = Φ−I disimpan sehingga identitas tidak membuang presisi dan +1 tak perlu representable di Φ; (3) kovariansi di domain ternormalisasi (similaritas, §Theory 7);
(4) `fxp_stats` per kelompok variabel agar saturasi/overflow/underflow dapat ditelusuri ke sumbernya.

## 3. Pembagian modul
| Modul | File | Isi |
|---|---|---|
| Aritmetika | `ekf_fixed_math.{h,c}` | add/sub/mul/div sat., rounding, shift, konversi, statistik |
| Divider | `ekf_fixed_divider.{h,c}` | baseline, restoring, non-restoring (SRT **tidak**) |
| Trigonometri | `ekf_fixed_trig.{h,c}` | CORDIC rotation mode, tabel atan, reduksi kuadran |
| Tipe & konfigurasi | `ekf_fixed_{types,config}.h` | `es_ekf_fixed_t`, `ekf_qcfg_t`, konstanta |
| Adapter | `ekf_fixed_adapter.{h,c}` | float↔fixed, normalisasi, tabel konfigurasi W16/W24/W32 |
| Model studi | `ekf_s4_ref_{float,fixed}.c`, `ekf_s4_models.h` | Model A dan B |
| Dataset | `s4_dataset.{h,c}` | 4 skenario sintetis deterministik |
| Konfigurasi terbaik + template DUT | `include/system4.h`, `src/model/system4_dut.c`, `tests/proof/system4_proof.c` | `s4_best_cfg()`, ambang penerimaan (PROPOSED), `s4_dut_t`, template TODO, bukti satu-file |
| Analisis | `s4_{metrics,compare,dense_ref}.{h,c}`, `tools/run_analysis.c` | metrik, lock-step, ablation, laporan |

## 4. File yang digunakan / dihasilkan
Input: 2 dokumen docx. Header ICD (transkripsi): `third_party/icd/*.h`. Output: `results/csv/*` (qformat_table, cordic_error, divider_compare, quant_hist,
pos_integration, range_analysis, quat_norm, ablation, summary_rmse, summary_counters, static_criteria, 12 `*_timeseries.csv`), `results/plots/*.svg` (12), `results/reports/analysis_report.md`.

## 5. Aturan kode (diperiksa otomatis)
C99, `-Wall -Wextra -Wpedantic` tanpa warning; tanpa `malloc/calloc/realloc/free` di seluruh repo; tanpa `float/double`/libm di inti fixed
(`ekf_fixed_{math,divider,trig}.c`, `ekf_s4_ref_fixed.c`) di luar blok `BEGIN/END-CONVERSION` (konversi I/O & parameter awal); array 1-D statis; tanpa rekursi;
loop bound tetap pada aritmetika inti. Pemeriksaan: `scripts/check_core_rules.py` (diuji negatif: mendeteksi pelanggaran yang disisipkan).
*Catatan jujur:* `divu_*` memakai jumlah iterasi = lebar bit dividen (variabel per panggilan, ≤ 64); perangkat keras akan memakai batas tetap. Pemeriksaan batas indeks array
dilakukan lewat struktur loop konstan, tidak ada `assert` runtime.

## 6. Hubungan dengan anggota 1–3
Lihat INTERFACE_CONTRACT §6. Ringkas: Model B mengikuti persamaan ICD dan menjadi *drop-in comparator*; setiap fungsi bertanda `_fx` memiliki padanan nama di ICD.

## 7. Deliverables
Library aritmetika, CORDIC, divider, adapter, 2 model + comparator, dataset sintetis, tool RMSE + plot + laporan, 6 test suite (464 check), 6 dokumen, Makefile + CMake, `run_all_tests.sh`.

## 8. Status modul (definisi: IMPLEMENTED=ada; COMPILED=bebas warning; TESTED=ada test lulus; VERIFIED=dibandingkan ke referensi independen dengan galat terukur; PENDING INTEGRATION=bergantung modul tim)
| Modul | IMPL | COMPILED | TESTED | VERIFIED | Catatan |
|---|---|---|---|---|---|
| `ekf_fixed_math` | ✔ | ✔ | ✔ (89 check) | ✔ vs `double` (property test 100k) | |
| `ekf_fixed_divider` | ✔ | ✔ | ✔ (62) | ✔ bit-exact vs `/` native, ≤ 1 LSB vs eksak | SRT tidak ada; latensi hardware nyata PENDING |
| `ekf_fixed_trig` (CORDIC) | ✔ | ✔ | ✔ (26) | ✔ vs libm: RMS 3,2·10⁻⁹, maks 1,4·10⁻⁸ | tidak dipakai pipeline (q memakai aproksimasi ICD); CORDIC tersedia untuk sudut besar |
| `ekf_fixed_adapter` | ✔ | ✔ | ✔ (49) | — | **PENDING INTEGRATION** (header tim asli) |
| Quant stress (Test 5) | ✔ | ✔ | ✔ (105) | ✔ galat maks. 0,5000 LSB | |
| Model A (float32) | ✔ | ✔ | ✔ | ✔ vs dense double: 5·10⁻⁵–8·10⁻⁵ | **bukan ES-EKF tim** |
| Model B (fixed) | ✔ | ✔ | ✔ (133 bersama Model A) | **Parsial**: hanya relatif terhadap Model A float32 | PENDING INTEGRATION; tanpa magnetometer |
| Dataset sintetis | ✔ | ✔ | ✔ | — | adapter data terbang **tidak dibuat** (tidak ada data) |
| Analisis RMSE / CSV / plot | ✔ | ✔ | ✔ | — | |
| Ekspor `.mem` golden vector | ✘ | | | | Role 5 |

## 9. Kriteria penerimaan (§12 brief)
| # | Kriteria | Status |
|---|---|---|
| 1 | Direktori & file utama | ✔ IMPLEMENTED |
| 2 | Header fixed-point | ✔ |
| 3 | Library aritmetika dikompilasi | ✔ COMPILED (gcc + CMake) |
| 4 | Saturation & rounding diuji | ✔ TESTED |
| 5 | Q-format dapat dikonfigurasi | ✔ (`FXP_DEFAULT_FRAC`, fungsi `*_q`, `ekf_qcfg_t`) |
| 6 | CORDIC dikompilasi & diuji | ✔ VERIFIED vs `sin/cos` |
| 7 | Divider dikompilasi & diuji | ✔ VERIFIED |
| 8 | Antarmuka pembanding float vs fixed | ✔ `s4_compare` — PENDING INTEGRATION dengan EKF tim |
| 9 | RMSE analysis dapat dijalankan | ✔ `build/bin/run_analysis` |
| 10 | Ekspor CSV | ✔ |
| 11 | Seluruh unit test lulus | ✔ 6/6 suite, 464 check (ctest 6/6) |
| 12 | Dokumentasi arsitektur & interface | ✔ |
| 13 | Tidak merusak kode anggota lain | ✔ (tidak ada kode anggota lain yang tersedia; tidak ada Git) |
| 14 | Tanpa angka hasil penelitian fiktif | ✔ semua angka dibangkitkan/ditautkan ke CSV; parameter PROPOSED diberi label |
**Peringatan:** skenario *static* lulus pemeriksaan ICD (drift ≤ 1 mm, D>0) **hanya pada study pipeline**; ini bukan bukti bahwa ES-EKF penuh memenuhi kriteria penerimaan.

## 10. Pekerjaan lanjutan yang direkomendasikan
1. Integrasi dengan modul Orang 1–3 begitu tersedia (ganti fungsi `_fx`); gunakan `s4_compare` sebagai harness.
2. Magnetometer (bounding U yaw) dan outlier gating.
3. ~~Modul inverse-sqrt untuk normalisasi quaternion~~ — selesai (Newton, 4 iterasi); pemetaan ke hardware PENDING.
4. *Residual carry* / akumulator posisi lebar untuk jalur nominal (temuan underflow Q16.16) dan stochastic rounding untuk Qd bias.
5. Telusuri anomali W32 vs W24 pada attitude `dynamic`; jalankan banyak seed noise; data terbang nyata.
6. Reciprocal+multiply pada Bierman/Thornton, pemetaan lebar ke DSP48, ekspor `stimulus_in.mem`/`expected_out.mem`, kandidat SRT.
