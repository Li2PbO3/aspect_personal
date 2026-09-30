# 热力学平衡应该用哪个压力？—— `melt_thermodynamic_equilibrium` 的压力场选择

> 结论先行：**可以，而且本工作区里 ASPECT 的标准熔体模型早就这么做了**
> （`melt_simple.cc` 一直用 `get_adiabatic_conditions().pressure()`）。
> 已经实现成一个运行期开关 `Pressure for thermodynamic equilibrium`
> （`fluid pressure` | `adiabatic pressure` | `solid pressure`，本分支默认
> `adiabatic pressure`），并用你 9 月 23 日那版 `case6_05_uniphase` 做了 A/B 对照。

---

## 1. 现状：平衡计算用的是哪个压力

`MeltThermodynamicEquilibrium::evaluate()` / `melt_fractions()` 里，
传给 `solve_eq_melt_fraction()` 的压力是

```cpp
// 改动前（三处：.cc 的 melt_fractions / evaluate 主循环 / MeltOutputs 填充）
const double pressure_for_material
  = (fluid_pressure_input != nullptr && ... && std::isfinite(...))
    ? fluid_pressure_input->fluid_pressure[q]     // ← 优先用 p_f（解向量里的 fluid pressure）
    : in.pressure[q];                             // ← 回退到固体压力 p
```

`FluidPressureInputs` 是在上一步 commit `dbda05b35`（"now TEQ calc is based on p_f"）里加的，
直接从解向量取 `fluid pressure` 分量。所以当前行为是 **p_f**。

这个选择在物理上是有道理的：两相流里固相和液相各自有压力，
`p_f = p_s − p_c`，而相平衡（熔融曲线）原则上由**熔体那一侧的压力**控制。
ASPECT 的 melt 体系里 `p_f` 是主变量，`p_s` 是后处理重构出来的
（`melt.cc:1331-1398`：`p_s = p_c/(1-φ) + p_f`）。

但代价就是你观察到的：`p_f` 里含有

* **动压**（对流/ Stokes 流动引起，本算例很小），以及
* **压实压力 `p_c`**（双相流特有，由 compaction 本构决定），

而 `p_c` 是耦合系统的一个自由度，它的空间结构（孔隙波/压实波）会被
`T_m(P) = T_m0 + A·P + B·P²` 直接映射到固相线温度上，于是平衡熔体分数、
`c_s/c_l`、`debug_K` 都继承了这个波结构。对"熔体如何产生、迁移"这类科学问题，
这不是我们想要的信号。

---

## 2. 量级：`p_c` 到底有多大（你现有算例的实测）

从 `carbon_melting/data/output-melt_TEQ_carbon_case6_05`（15000 单元 × 150 km 柱，
2D 但 x 方向只有 0/5/10 m 三个节点，实际上是 **1D 柱**）

| 量 | 值 |
|---|---|
| `p_f` 范围 | 3.060×10⁹ – 8.063×10⁹ Pa |
| `p_c` 最大值 | **3.94×10⁶ Pa**（≈ 3.9 MPa） |
| `p_c` RMS | 4.7×10⁵ Pa |
| `p_c/p_f` 中位数 | ~1.7×10⁻⁸ |
| `p_c` 的垂直波长 | ~6.8–7.9 km |

把它换算成固相线温度扰动 `ΔT_m = A·p_c`（`A` = 45/112/38 K/GPa）：

| 组分 | RMS ΔT_m | 最大 ΔT_m |
|---|---|---|
| dunite | 0.020 K | 0.15 K |
| morb | 0.050 K | 0.38 K |
| cmorb | 0.017 K | 0.13 K |

也就是说，**压实压力对平衡温度的影响 < 0.4 K**，而同一算例里温度的侧向
异常是 1.5–40 K。所以从"回答科学问题"的角度，这个信号确实可以忽略——
你的判断是对的。

---

## 3. 改动：`Pressure for thermodynamic equilibrium`

新增运行期参数（`melt_thermodynamic_equilibrium.cc` 的 `declare_parameters` /
`parse_parameters`），三处压力取值统一走新函数

```cpp
double MeltThermodynamicEquilibrium<dim>::
equilibrium_pressure (const MaterialModel::MaterialModelInputs<dim> &in,
                      const unsigned int q,
                      const FluidPressureInputs<dim> *fluid_pressure_input) const
```

* `fluid pressure` —— 保持原行为（p_f）。
* **`adiabatic pressure`（本分支默认）** ——
  `this->get_adiabatic_conditions().pressure(in.position[q])`，
  即 `AdiabaticConditions::ComputeProfile` 从 **Surface pressure** 出发、
  用侧向平均参考密度 `∫ρ_ref g dz` 积分出来的**静岩压力剖面**。
  它只是深度的光滑函数，逐点、且不随流动变，完全不含压实波。
  在参考剖面还没建好（第一个时间步里材料模型被用来构造它自己）时，
  自动回退到 `in.pressure[q]`，避免自引用。
* `solid pressure` —— 用 `in.pressure[q]`（Stokes 的 `p` 变量）。

> ⚠️ 只改**热力学平衡**用的压力。密度、剪切/压实粘度、熔体密度里对压力的
> 依赖（`pressure_temperature_viscosity_factor`、`melt_compressibility`）
> 仍然按原样使用 p_f。若也想统一，需要再单独讨论。

---

## 4. A/B 验证（你 9/23 那版 `case6_05_uniphase`，12 步 / 200 yr）

* 两个 deck 只差 `Output directory` 和这一个新参数，其余逐字相同；
* 同一二进制 `build_TEQ_bulk_release/aspect`（HEAD `390046b3b`），各 8 MPI ≈ 90 s。

| 场 | max&#124;Δ&#124; | RMS Δ | 尺度 |
|---|---|---|---|
| **T** | **0.000** | 0 | 1.59×10³ K |
| **density** | **0.000** | 0 | 3.4×10³ |
| **viscosity** | **0.000** | 0 | 9.9×10²³ |
| **porosity** | 6.1×10⁻⁷ | 7.9×10⁻⁸ | 3.0×10⁻³（相对 RMS 2×10⁻⁴） |
| permeability | 8.2×10⁻²⁰ | 9.3×10⁻²¹ | 2.7×10⁻¹⁶ |
| p_c | 5.3×10³ | 2.3×10² | 3.9×10⁶ |
| p_f | 1.2×10⁴ | 5.3×10² | 8.1×10⁹ |
| debug_1（K₀） | 6.9×10⁻³ | 8.7×10⁻⁴ | 2.0×10² |

读数：**换成静岩压力后，温度场逐位不变，孔隙度只动了 ~10⁻⁷（相对 2×10⁻⁴）**，
压力场那点差别纯粹是 Stokes 解对 10⁻⁷ 量级 RHS 扰动的舍入响应。
`debug_K` 几乎逐位相同 —— 平衡内核本身就是对压力不敏感的。

这符合预期：`p_c` 只有几 MPa，而这里 `T` 距固相线有上百 K。
在**真正有熔融**的算例里，效应会正比于 `A·p_c`（最多 ~0.4 K），
相对温度异常仍是小量。

---

## 5. 顺带核对：平衡内核与参考实现一致

对着你工作区里的 `melt_equilibrium_comparison/aspect_kernel_reference.cpp`
（据称是当前内核的逐行转写，用 `g++ -O2` 编成 `teq_ab/akr`）
在 100 km 处逐点核对（`p_f = 6.395463×10⁹ Pa`、`T = 1298.95 °C`、
`c_bulk = 0.7 / 0.299333 / 6.66665×10⁻⁴`）：

| 量 | 算例存储值 | 参考内核 | 差 |
|---|---|---|---|
| `debug_1` = K₀ | 14.3543 | 14.35433 | ~10⁻⁵ |
| `debug_2` = K₁ | 7.7522 | 7.75224 | ~10⁻⁵ |
| `debug_3` = K₂ | 6.1756×10⁻⁴ | 6.1760×10⁻⁴ | ~10⁻⁶ |
| `c_l`(dunite) | 0.048771 | 0.048771 | 逐位 |
| `c_l`(morb) | 0.038616 | 0.038616 | 逐位 |
| `c_l`(cmorb) | 0.912613 | 0.912612 | ~10⁻⁶ |
| `φ` | 1.130092×10⁻⁴ | 1.169728×10⁻⁴ | 3.4 % |

结论：**平衡内核本身是对的**——平衡常数与液相浓度与参考实现吻合到 ~10⁻⁵。
唯一 3.4 % 的差异在 `φ`，量级 ~4×10⁻⁶ 绝对、且与 A/B 两组之间的差
（~10⁻⁷）完全不同尺度，最可能来自温度场以 float32 存储、而 `φ` 对 `T`
敏感（`∂lnφ/∂T` 在固相线附近可达 ~14 K⁻¹，故 0.015 K 的温度差就足以解释）。

顺带两点观察（供参考，与本次改动无关）：

* 液相场满足 `c_l = c_bulk/(f+(1−f)K)`，`Σ_cl` 接近 1 只是因为
  cmorb 高度不相容（`K₂ ~ 6×10⁻⁴`），其液相质量分数本来就接近 1；
  这不是归一化错误。
* `porosity` 与 `T` 从第 0 步到第 12 步几乎逐位不变，速度 ~10⁻⁷ m/yr，
  即这个短算例实际处于"温度场冻结"的状态。这也解释了 A/B 差异只有 10⁻⁷ 量级。

---

## 6. 建议

1. **采用静岩压力做平衡**（默认已改）。它把热力学状态变成
   `(depth, T, c) → (φ, c_s, c_l)` 的确定函数，压实波不再进入平衡，
   你也不用再为"为什么 K 里有波"写解释。量级上损失 < 0.4 K 固相线位移。
2. 若想把这件事写进论文/文档，可用的表述：
   * 两相流的相平衡原则上由熔体压力控制（`p_f = p_s − p_c`）；
   * 但 `p_c/p_lith = O(ξ C / ρ g H) ≲ 10⁻³`，对固相线温度的影响
     < 0.4 K，远小于热异常；
   * 因此平衡计算采用**参考静岩压力**（等同 `melt_simple` 的既有做法），
     以把压实波这一双相流内模从热力学闭合中剔除。
3. 注意适用边界：若将来做**动压很大**的算例（例如强烈对流、`p_dyn` 达到
   10⁸–10⁹ Pa），`adiabatic pressure` 与 `p_f` 的差别就不再只是压实压力，
   平衡会明显不同。届时可切回 `fluid pressure`，或考虑第三个选项。
4. `porosity` 与参考内核的 3.4 %（绝对 ~4×10⁻⁶）差异建议顺手确认一下
   （见 §5）：怀疑是 float32 存储的 `T` 与 `φ(T)` 的高灵敏度所致，
   但如果要发表绝对熔体分数，值得把它压到 10⁻⁷ 量级。

## 附：复现方式

```bash
# 两个 deck（只差输出目录与该参数）
teq_ab/fluid.prm
teq_ab/adiabatic.prm
# 运行（各 8 MPI，约 90 s）
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
mpirun -np 8 --oversubscribe ./build_TEQ_bulk_release/aspect teq_ab/fluid.prm
mpirun -np 8 --oversubscribe ./build_TEQ_bulk_release/aspect teq_ab/adiabatic.prm
# 定量对照
python3 teq_ab/compare.py
```

`teq_ab/ab_summary.json`、`teq_ab/ab_summary2.json` 是上面表格里数字的原始输出；
`teq_ab/akr` 是编译好的参考内核。

## 附：独立验证（2026-09-29）

见 **`TEQ_PRESSURE_VERIFICATION.md`**（配套算例
`region5_carbon/case6_practical_0724/verify_teq_pressure/`）。要点：

* 三个口径用的压力被逐点确证（内核重写 vs ASPECT step 0，吻合 1e-10）：
  `fluid` → `p_f`，`adiabatic` → `Surface pressure + ρ_ref g z`。
* §4 的 A/B 结论成立，但要注意它**测的只是压力口径之差**：这里重新测到
  `φ(fluid) − φ(adiabatic)` max **6.08e-7**、rms 7.9e-8，与 `p_c` 的相关系数
  **+0.987**；相对波幅只有 **2.6e-5**。
* 旧口径造成的外部初始条件不一致 = max **6.3e-7**；新口径把这一项压到 ≤7.7e-8。
* 但 step 0 的总失配（两种口径都）是 **~2.1e-6**，主导项是 `.prm` 的 cMORB
  熔融曲线 `A`、`B` 被舍入成 38.0e-9 / −4.00e-18（与 IC 生成器的标定值不一致），
  不是压力口径。想让初始条件真正落在平衡态上，这一项也要改。
* `p_f − P_lith` 的 max 是 **9.25 MPa**（不只 `p_c` 的 3.94 MPa，还有动压），
  对应 morb 熔点位移最大 **+0.53 K**，比 §2 的"< 0.4 K"略大。
* `solid pressure` 在 t=0 不可用（第一次材料模型求值时 `in.pressure` 还不是
  静岩剖面，见验证报告 §6）。
