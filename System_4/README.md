# System 4 — Fixed-Point & Bit-Accurate C Model
**FPGA Sensor Fusion UAV — 15-DOF Quaternion ES-EKF** · Peran: *Orang 4 — Analis Kuantisasi Fixed-Point & Bit-Accurate Model*

## 1. Tujuan
System 4 **tidak membuat algoritma EKF baru**. System 4 menjawab pertanyaan: *berapa bit, format Q berapa, dan skala normalisasi
apa* yang dibutuhkan agar ES-EKF milik Orang 1–3 (Thornton UDUᵀ, Bierman, injeksi/reset error) tetap akurat saat
dipindahkan ke aritmetika fixed-point untuk FPGA. Deliverable: library aritmetika bit-accurate, CORDIC dan divider
berorientasi hardware, adapter float↔fixed yang kompatibel ICD, serta pipeline analisis RMSE yang *dapat diulang*.

## 2. Hubungan dengan sistem utama
| Anggota | Modul | Posisi di System 4 |
|---|---|---|
| 1 | `ekf_predict_nominal()`, Φₖ | Dimodelkan oleh *study pipeline* (`src/model/ekf_s4_ref_*.c`) mengikuti persamaan ICD — **PENDING INTEGRATION** |
| 2 | Measurement update, `ekf_inject_and_reset()` | Idem (GNSS pos+vel, barometer; **magnetometer belum**) |
| 3 | Thornton MWGS, Bierman | Idem |
| 5 | Dataset, testbench | System 4 membuat dataset *sintetis*; ekspor `.mem` **belum** dibuat |

> Study pipeline (Model A float32 / Model B fixed) hanyalah wahana ukur kuantisasi. Ia bukan ES-EKF tim dan harus
> diganti per-fungsi saat modul Orang 1–3 tersedia (lihat `docs/INTERFACE_CONTRACT.md`).

## 3. Prasyarat & instalasi
`gcc` (C99), `make`, `python3` (hanya stdlib, untuk plot SVG), opsional `cmake ≥ 3.10`. Tidak ada dependensi lain.
Header ICD (`ekf_types.h`, `ekf_sensor_packet.h`) ditranskripsi dari dokumen ke `third_party/icd/` karena repositori tim
tidak tersedia secara lokal; letakkan versi resmi tim lebih awal di include path dan salinan ini otomatis dilewati.

## 4. Compile
```bash
make            # library, 6 unit test, system4_proof, run_analysis   (-std=c99 -Wall -Wextra -Wpedantic)
```
Alternatif CMake:
```bash
cmake -S . -B build-cmake && cmake --build build-cmake && (cd build-cmake && ctest)
```

## 5. Menjalankan test
```bash
./scripts/run_all_tests.sh              # build + cek aturan kode + semua unit test (exit code ≠ 0 bila gagal)
./scripts/run_all_tests.sh --analysis   # + analisis 60 s dan plot
```
Analisis saja (≈ 3 detik):
```bash
make analysis        # build/bin/run_analysis results 60  &&  python3 scripts/plot_results.py results
```
Output: `results/csv/` (CSV), `results/plots/` (SVG), `results/reports/analysis_report.md` (dibangkitkan otomatis;
semua angka hasil eksekusi, tidak ada yang diketik manual).

## 6. Struktur folder
```
System_4/
├── include/
│   ├── system4.h               ← HEADER UTAMA: konfigurasi terbaik, ambang penerimaan, template DUT
│   ├── core/                   aritmetika & tipe: ekf_fixed_{math,divider,trig,config,types}.h
│   ├── model/                  adapter ICD + model studi: ekf_fixed_adapter.h, ekf_s4_models.h
│   └── analysis/               dataset sintetis, metrik, komparator: s4_{dataset,metrics,compare,dense_ref}.h
├── src/{core,model,analysis}/  implementasi (inti fixed-point: core/* dan model/ekf_s4_ref_fixed.c)
│   └── model/system4_dut.c     ← TEMPLATE ALGORITMA (S4_DUT_TEMPLATE) + binding S4_DUT_STUDY
├── tests/
│   ├── unit/                   6 unit test (math, divider, cordic, adapter, quant_stress, models)
│   └── proof/system4_proof.c   ← SATU FILE bukti dummy-test konfigurasi terbaik
├── tools/run_analysis.c        analisis penuh → CSV + plot + laporan
├── scripts/                    run_all_tests.sh, check_core_rules.py, plot_results.py
├── third_party/icd/            salinan header ICD (ekf_types.h, ekf_sensor_packet.h)
├── docs/                       5 dokumen (lihat docs/README.md)
├── results/{csv,plots,reports,proof}/
├── build/{obj,bin}/            (dihasilkan, di-.gitignore)
├── Makefile, CMakeLists.txt
└── Kesepakatan*.docx           dokumen referensi (tidak diubah)
```

### Bukti konfigurasi terbaik + template (satu file)
```bash
make proof        # = build/bin/system4_proof results  → results/proof/best_config_proof.{md,csv}
```
`tests/proof/system4_proof.c` (A) menguji aritmetika dummy, (B) menjalankan semua kandidat 16/24/32-bit pada 4 skenario dan memilih
yang termurah yang lolos semua ambang di `system4.h`, (C) memastikan hasilnya sama dengan `s4_best_cfg()`, (D) membuktikan template
TODO **ditolak** (harness bisa gagal). Cara memakai sebagai template: isi tahap `tpl_*` di `src/model/system4_dut.c`
(atau buat `s4_dut_t` sendiri), lalu jalankan `make proof`.

## 7. Status singkat
Lihat tabel status per modul di `docs/SYSTEM_4_PREBUILD.md` §8 dan keterbatasan di `docs/NUMERICAL_ERROR_ANALYSIS.md` §8.
Git: **tidak ada** operasi Git yang dilakukan (repositori tidak tersedia lokal; folder ini bukan repositori Git).
