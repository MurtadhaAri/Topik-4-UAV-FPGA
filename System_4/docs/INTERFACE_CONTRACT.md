# INTERFACE_CONTRACT — Kontrak Integrasi System 4

Sumber kontrak: `Kesepakatan_Teknis_ICD_ES-EKF.docx` (ICD) dan `Kesepakatan Teknis.docx`. Bagian yang tidak ada di ICD ditandai **PROPOSED**.

## 1. Header yang digunakan
| Header | Asal | Status |
|---|---|---|
| `third_party/icd/ekf_types.h` | Transkripsi ICD Bagian I butir 5 (guard `EKF_TYPES_H`) | **Salinan** — repositori tim tidak tersedia lokal; ganti dengan versi resmi (lebih awal di include path) |
| `third_party/icd/ekf_sensor_packet.h` | Transkripsi ICD Role 5 butir 1 (`sensor_type_e`, `sensor_frame_t`) | Salinan, idem |
| `include/core/ekf_fixed_config.h` | System 4 | `NOM_DIM`, `ERR_DIM`, `U_UPPER_SIZE`, `U_IDX` **tidak didefinisikan ulang** — diimpor dari `ekf_types.h` |
| `include/core/ekf_fixed_types.h` | System 4 | Cermin fixed dari struct ICD + `ekf_qcfg_t` (konfigurasi Q) + `sensor_frame_fixed_t` |
| `include/core/ekf_fixed_{math,divider,trig}.h`, `include/model/ekf_fixed_adapter.h`, `ekf_s4_models.h`, `include/system4.h` | System 4 | — |

**Penyimpangan dari ICD (disengaja, harus disetujui tim):** ICD berkata `typedef float ekf_scalar_t; // Akan diubah oleh Orang 4`.
System 4 **tidak mengubah** `ekf_scalar_t`; fixed-point disediakan sebagai tipe paralel (`es_ekf_fixed_t`) sehingga header tim tidak rusak.
Gravitasi mengikuti ICD **9,80665 m/s²** (brief tugas menyebut 9,81; ICD yang dipakai).

## 2. Kontrak data (tidak boleh diubah)
```c
#define NOM_DIM 16     // x_nom: p[0..2] v[3..5] q[6..9] ba[10..12] bg[13..15]
#define ERR_DIM 15     // dx   : dp[0..2] dv[3..5] dth[6..8] dba[9..11] dbg[12..14]
#define U_UPPER_SIZE 120
#define U_IDX(r,c) (((c)*((c)+1))/2 + (r))   // r <= c, packed upper-triangular, column-major, DIAGONAL DISIMPAN (=1,0)
```
Dikunci oleh `tests/unit/test_adapter.c`: `U_IDX` bijektif ke 0..119, `U_IDX(14,14)=119`, `offsetof(q)=6`, `offsetof(ba)=10`, `offsetof(bg)=13`.
Konvensi: NED, FRD, Hamilton scalar-first `[qw,qx,qy,qz]`, q: body→navigation, perturbasi lokal/kanan, satuan SI.

## 3. Format fixed-point
Jalur **nominal** (int32, `ekf_qcfg_t`): p `Q16.16`, v `Q8.24`, q `Q2.30`, ba/bg `Q4.28`.
Jalur **error/kovariansi** (storage W16/W24/W32; domain **ternormalisasi**): dx `Q4.x`, U `Q7.x`, D `Q1.x`, ΔΦ `Q1.x` (`src/model/ekf_fixed_adapter.c`, `k_cfg`).
Normalisasi: `dxₙ = dx/sᵢ`, `Uₙ(i,j) = U(i,j)·sⱼ/sᵢ`, `Dₙ(i) = D(i)/sᵢ²`, `sᵢ = 2^s4_scale_log2[i]` = {16,16,16, 16,16,16, 0.5,0.5,0.5, 2,2,2, 0.25,0.25,0.25} **(PROPOSED)**.
Paket sensor fixed (`sensor_frame_fixed_t`): `dt` Q2.30; `imu_acc`, `imu_gyr` Q8.24; `meas[0..2]` format p; `meas[3..5]` format v. Paket float tetap `sensor_frame_t` (ICD).
Konversi unit (derajat/ms → SI) **hanya** di parser sensor (ICD 1.4); System 4 menerima SI.

## 4. Fungsi adapter
| Fungsi | Keterangan |
|---|---|
| `ekf_float_to_fixed(const es_ekf_t*, es_ekf_fixed_t*, const ekf_qcfg_t*)` | normalisasi + kuantisasi (rounding half-away, saturasi) |
| `ekf_fixed_to_float(const es_ekf_fixed_t*, es_ekf_t*, const ekf_qcfg_t*)` | de-kuantisasi + de-normalisasi |
| `sensor_frame_to_fixed(const sensor_frame_t*, sensor_frame_fixed_t*, const ekf_qcfg_t*)` | |
| `s4_float_Un/Dn`, `s4_fixed_Un/Dn` | akses terukur ternormalisasi untuk analisis |
Perbedaan dari brief: tiap konverter menerima argumen `cfg` tambahan karena Q-format dapat dikonfigurasi. Adapter memakai `double` (I/O dan analisis saja).
Test: `test_adapter` — galat round-trip ≤ ½ LSB (+ε float32) untuk p, v, q, dx, U, D di W16/W24/W32; diagonal U tepat 1,0.
**Cakupan konversi:** state nominal, dx, U, D, paket sensor. **Tidak** mencakup operasi EKF (prediksi, update) milik Orang 1–3.

## 5. Primitif yang disediakan System 4 untuk Orang 1–3
```c
int32_t fxp_add_sat/sub_sat/mul/mul_round/div/shift_right/saturate(...);  float fxp_to_float(int32_t); int32_t fxp_from_float(float);
int32_t fxp_mul_q(a,b,fa,fb,fr); fxp_mul_round_q(...); fxp_mul_round_wl(...,wl); fxp_div_q(a,b,fa,fb,fr); fxp_sat_wl(v,wl);
int64_t fxp_rshift_round(int64_t v,int sh);   // round-half-away; sh<0 = shift kiri saturating
void    cordic_sincos(int32_t angle_q3_29, int32_t *sin_q2_30, int32_t *cos_q2_30);
int64_t divs64_round(int64_t n,int64_t d,div_mode_t mode,div_info_t*); int32_t div_q(..., mode, info);
fxp_stats_*  // penghitung saturation/overflow/underflow/div0 per kelompok variabel
```
Semantik rounding: `fxp_mul`/`fxp_shift_right` = **floor**; `*_round` = **half away from zero**; `fxp_div` = **truncate toward zero**. `MUL_FXP(a,b)` = `fxp_mul` dengan `FXP_DEFAULT_FRAC` (Q16.16, dapat di-override `-DFXP_DEFAULT_FRAC=n`).

## 6. Kontrak per anggota (apa yang System 4 asumsikan / sediakan)
| Anggota | Fungsi tim | Yang diasumsikan | Pengganti di System 4 (yang harus diganti) |
|---|---|---|---|
| 1 | `ekf_predict_nominal()` | ICD Role 1 butir 2: Euler order-1, q Taylor-normalize | `nominal_step()` di `ekf_s4_ref_fixed.c` / blok nominal `s4f_step()` |
| 1 | Φₖ | ICD Role 1 butir 3 (blok 3×3) | `build_phi_fx()` / `build_phi()` — **catatan:** System 4 menyimpan **ΔΦ = Φ−I** (identitas tidak disimpan) |
| 2 | model h(x), H, sequential scalar update | H = eᵢ untuk GNSS pos/vel/baro (ICD) | `bierman_fx(m,i0,y_n,R_n)` dipanggil per kanal |
| 2 | `ekf_inject_and_reset()` | ICD Role 2 butir 3 | `inject_reset_fx()` / `inject_reset()` |
| 3 | `bierman_scalar_update()` | ICD Role 3 butir 2 | `bierman_fx` / `bierman` |
| 3 | Thornton MWGS | ICD Role 3 butir 3, dasar P = ΦUDUᵀΦᵀ + Q_d diagonal | `thornton_fx` / `thornton` |
| 5 | dataset & vektor uji RTL | `sensor_frame_t` | `s4_dataset.c` (sintetis); ekspor `stimulus_in.mem` / `expected_out.mem` **belum ada** |
**Yang System 4 butuhkan dari tim:** (a) `ekf_types.h` resmi; (b) bila modul tim memakai `float`, mereka cukup memanggil adapter untuk uji kuantisasi;
(c) keputusan Q-format akhir (rekomendasi: Q_FORMAT_ANALYSIS §6); (d) nilai Q_c, noise IMU, laju prediksi kovariansi (kini **PROPOSED**); (e) magnetometer.
**Yang System 4 tidak lakukan:** mengubah ICD atau kode anggota lain; membuat EKF baru.

## 6b. Template DUT (`include/system4.h`)
Anggota tim dapat menghubungkan implementasi fixed-point ES-EKF miliknya dengan mengisi `s4_dut_t` (`init`, `step`, `export_state`, opsional `pd_clamps`).
Urutan tahap kanonik: nominal → time update kovariansi → update skalar per kanal → injeksi/reset. Kerangka siap-isi: `S4_DUT_TEMPLATE`
(`src/model/system4_dut.c`, tahap `tpl_*` bertanda TODO). `tests/proof/system4_proof.c` menjalankan DUT lewat harness yang sama dan menolak template kosong.
**Batasan:** komparator `s4_compare` (analisis penuh) masih terikat ke pipeline studi; hanya harness proof yang generik terhadap DUT.

## 7. Pertanyaan terbuka untuk tim
1. Laju prediksi kovariansi: tiap sampel IMU atau didesimasi? (pipeline studi: tiap 10; berpengaruh besar pada underflow Q_d).
2. Normalisasi quaternion: ICD mengizinkan Taylor orde-1 *atau* modul inverse-sqrt; terukur Taylor gagal untuk injeksi besar; Model B kini memakai Newton inverse-sqrt (NUMERICAL_ERROR_ANALYSIS §5).
3. Apakah ΔΦ (tanpa identitas) dapat dijadikan konvensi bersama untuk jalur fixed?
4. Ambang penerimaan akurasi fixed-vs-float (ICD hanya mendefinisikan kriteria static 1 mm dan D>0).
