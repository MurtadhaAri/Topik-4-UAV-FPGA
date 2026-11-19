# FIXED_POINT_THEORY — Dasar Teori Fixed-Point yang Dipakai System 4

## 1. Floating-point vs fixed-point
| | Floating-point (float32) | Fixed-point (int32 + Q-format) |
|---|---|---|
| Representasi | tanda · mantissa 24 bit · eksponen 8 bit | integer bertanda; titik biner *implisit* pada posisi *F* |
| Resolusi | relatif (≈ 6·10⁻⁸ × nilai) | **absolut** konstan 2⁻ᶠ |
| Rentang dinamis | ≈ 10⁻³⁸…10³⁸ | 2^(W−1) · 2⁻ᶠ, sempit |
| Biaya FPGA | FPU/DSP besar, latensi tinggi | DSP48 + shifter, latensi deterministik |
| Risiko | rounding relatif | **overflow, underflow, saturasi**, bias pembulatan |
Konsekuensi: pada fixed-point, *skala variabel adalah bagian dari desain*; itulah alasan normalisasi (Q_FORMAT_ANALYSIS §3).

## 2. Two's complement dan notasi Q
Word *W* bit bertanda: nilai integer *n* ∈ [−2^(W−1), 2^(W−1)−1]. **Notasi System 4:** `Qm.n` = *m* bit integer
**termasuk bit tanda**, *n* bit pecahan, W = m+n (Q1.15 = 16 bit, Q1.23 = 24 bit, Q2.30 = 32 bit, Q16.16 = 32 bit).
Nilai nyata: x = n·2⁻ⁿ.  Properti (dihitung otomatis di `results/csv/qformat_table.csv`):

* Resolusi (LSB) = 2⁻ⁿ  ·  maks. galat rounding = ½·2⁻ⁿ  ·  maks. galat truncation = 2⁻ⁿ
* Min = −2^(m−1)  ·  Max = 2^(m−1) − 2⁻ⁿ  (asimetris: **+2^(m−1) tidak dapat direpresentasikan**)

Akibat penting: pada Q1.15/Q1.23, +1,0 *tidak representable* → disaturasi ke 1−LSB (3,05·10⁻⁵ dan 1,19·10⁻⁷). Di Q2.30
+1,0 = 1<<30 tepat (test: `fxp_from_double_q(1.0,32,30)==1<<30`). Diagonal U bernilai 1,0 sehingga **U butuh ≥ 2 bit integer**.

## 3. Quantization
Nilai kontinu *x* → bilangan terdekat berkelipatan LSB. Galat *e = x̂ − x*. Untuk rounding-to-nearest dan sinyal “aktif”
(melintasi banyak LSB) *e* ≈ seragam pada [−½,½]·LSB: RMS = LSB/√12. Test `test_quant_stress` mengukur galat maks. 0,5000 LSB
pada semua format; `results/plots/plot6_quant_hist_*.svg` memperlihatkan histogramnya.
Nilai < ½ LSB menjadi 0 (**underflow**); System 4 menghitung kejadiannya (`fxp_counters_t.underflow`).

## 4. Rounding
| Metode | Fungsi | Perilaku | Catatan |
|---|---|---|---|
| Floor (arithmetic shift) | `fxp_mul`, `fxp_shift_right` | → −∞ | bias rata-rata −½ LSB per operasi (test: −0,5 LSB → −1, +0,5 LSB → 0) |
| Round half away from zero | `*_round`, `fxp_rshift_round` | simetris | tanpa bias tanda; dipakai di pipeline studi |
| Truncate toward zero | `fxp_div` (C99 `/`) | → 0 | galat < 1 LSB |
Bias pembulatan *berakumulasi*: integrasi `p += v·Δt` dengan inkremen tetap menghasilkan galat sistematis, bukan acak.
Terukur (`results/csv/pos_integration.csv`, v=0,01 m/s, 60 s @500 Hz, inkremen 2·10⁻⁵ m/langkah): Q16.16 (LSB 1,5·10⁻⁵)
galat relatif **23,7 %**; Q12.20/Q8.24 galat 0,14 % (bias konstan karena membulatkan 20,97 LSB → 21).

## 5. Saturation
Hasil di luar [min,max] di-*clip*, bukan *wrap-around*. Wrap-around pada filter menyebabkan lonjakan tanda (mis. +1,0 → −1,0).
`fxp_sat_wl(v, wl)` men-clip ke word 16/24/32 bit; `fxp_saturate` ke int32. Dua penghitung dibedakan:
`sat` (hasil di luar word target) dan `ovf` (nilai pra-clip bahkan tidak muat di int32 *container*).

## 6. Bit growth
* Perkalian *a(Qm₁.n₁)·b(Qm₂.n₂)* → integer W₁+W₂ bit dengan *n₁+n₂* pecahan. Karena itu semua perkalian memakai
  intermediate **int64**, lalu digeser `(a*b) >> (fa+fb−fr)` ke format hasil `fr` (format operan boleh berbeda).
* Penjumlahan *N* suku: +⌈log₂N⌉ bit. Akumulator dot-product Thornton (30 suku) di-*int64*; pembulatan dilakukan **sekali** di akhir.
* Pembagian: pembilang digeser kiri `fr+fb−fa` bit sebelum membagi.
* Pembatas: |a·b| ≤ 2⁶² sehingga int64 tak pernah overflow untuk satu perkalian; penjumlahan banyak suku memakai `fxp_add64_sat`.

## 7. Normalisasi sebagai transformasi similaritas
Dengan *S = diag(sᵢ)*, sᵢ = 2^k (pangkat dua → normalisasi = *shift*): dxₙ = S⁻¹dx. Agar P = UDUᵀ tetap berbentuk UDUᵀ dengan U
unit-upper-triangular:  **Uₙ(i,j) = U(i,j)·sⱼ/sᵢ**, **Dₙ(i) = D(i)/sᵢ²**  (cek: Uₙ Dₙ Uₙᵀ = S⁻¹ U D Uᵀ S⁻¹). Diagonal Uₙ tetap 1.
Untuk H = eᵢ (GNSS, baro) Hₙ = eᵢ, y_n = y/sᵢ, R_n = R/sᵢ² — sehingga baris Jacobian tidak berubah.

## 8. CORDIC & divider (ringkas)
* **CORDIC rotation mode**: x' = x − d·(y≫i), y' = y + d·(x≫i), z' = z − d·atan(2⁻ⁱ), d = sign(z). Gain Π√(1+2⁻²ⁱ)→1,6468;
  dikompensasi dengan x₀ = K = 0,60725293… (Q2.30 = 652032874). Rentang konvergensi |z| ≤ 1,7433 rad ⊃ [−π/2, π/2] → reduksi kuadran.
  Galat ≈ 2⁻ᴺ (≈ 2 bit/iterasi, terukur di `cordic_error.csv`).
* **Divider**: restoring (shift, trial-subtract, restore) dan non-restoring (sisa boleh negatif, +1 langkah koreksi); keduanya
  1 bit/iterasi (radix-2). **SRT tidak diimplementasikan.**
