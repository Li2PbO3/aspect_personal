# `Pressure for thermodynamic equilibrium` 改进的验证报告

> 结论先行：
>
> 1. **改进已生效。** 运行期开关实现正确：`fluid pressure` → 用 `p_f`，
>    `adiabatic pressure` → 用 `P_ad = Surface pressure + ∫ρ_ref g dz`，
>    两者与一份独立的平衡内核重写在 **1e-10** 量级上逐点吻合。
> 2. **改进前（用 `p_f`）确实让外部生成的初始条件对 ASPECT 不平衡，但量级很小**：
>    平衡压力偏离初始条件的静岩压力 rms 1.11 MPa、最大 9.25 MPa，
>    由此产生的孔隙度不一致最大 **6.1e-7**（rms 7.9e-8），只有 3e-3 孤立波幅值的
>    **2.6e-5**，且与压实压力 `p_c` 的相关系数 **+0.987**——确实带着压实波的印迹。
>    改成 `adiabatic pressure` 后这一项降到 ≤7.7e-8。
> 3. **但这一项不是初始不平衡的主要来源。** step 0 的
>    `φ_ASPECT − φ_IC` 对两种压力口径都是 **~2.1e-6**（rms 1.4e-6），
>    主导项是 **.prm 的 cMORB 熔融曲线系数与 IC 生成器不一致**（见 §5），
>    压力口径只占其中约 1/4。想让初始条件真正落在平衡态上，这一项也得一起改。

验证算例（新建）：`region5_carbon/case6_practical_0724/verify_teq_pressure/`
（三个 deck + 运行脚本 + 分析脚本 + 图）。原始数据：
`analysis_summary.json`；图：`analysis.png`。

---

## 1. 被验证的改动与二进制

* 改动位于 `aspect_personal` 工作区（**尚未提交**，HEAD = `390046b3b`）：
  * `include/aspect/material_model/melt_thermodynamic_equilibrium.h`
  * `source/material_model/melt_thermodynamic_equilibrium.cc`
* 新增运行期参数 `Pressure for thermodynamic equilibrium`
  ∈ {`fluid pressure`, `adiabatic pressure`, `solid pressure`}，默认
  `adiabatic pressure`；三处压力取值统一走新函数 `equilibrium_pressure()`。
* `fluid pressure` 分支与改动前的内联表达式**逐字等价**，因此**同一个二进制**
  就能同时复现"改动前"和"改动后"的行为，A/B 对照不需要换二进制、不存在编译差异。
* 使用的二进制：`aspect_personal/build_TEQ_bulk_release/aspect`
  （37662832 字节，2026-09-29 22:25；`find source include -newer` 无更新，
  即源码与二进制一致）。

---

## 2. 验证算例

基例取自 `region5_carbon/case6_practical_0724` 的
`melt_TEQ_carbon_case6_05_uniphase.prm`（15000×1 单元 1D 柱，100–250 km，
200 ppm CO₂ 背景 + 225 km 处孤立波，`φ_peak = 3e-3`）。

新建 `verify_teq_pressure/`，三个 deck **只差 `Output directory` 和
`Pressure for thermodynamic equilibrium`**，其余逐字相同；都读同一个外部
python 生成器写出的静岩初始条件
（`initial_composition_for_carbon_case06_uniphase.txt`）。各跑 45 yr / 3 步。
验证只看 **step 0**：ASPECT 在 step 0 会把 prescribed 场
（`porosity`、`*_liquid`、`melting_rate`、`debug_*`）用**自己的平衡内核**重算一遍，
所以 step 0 的输出就是"ASPECT 眼中的初始平衡态"，与 IC 文件之差就是初始不平衡。

```bash
cd /home/zhiqianli/workspace/region5_carbon/case6_practical_0724
bash verify_teq_pressure/run.sh 8      # 三个算例，各 ~20 s
python3 verify_teq_pressure/analyze.py # 表格 + analysis_summary.json + analysis.png
```

---

## 3. 结果

### 3.1 开关是否生效（内核重写 vs ASPECT step 0，去掉上下 1 km 边界节点）

| 口径 | ASPECT vs `内核@应该用的压力`（max） |
|---|---|
| `fluid pressure` | **1.36e-10**（用 `p_f`） |
| `adiabatic pressure` | **1.16e-10**（用 `Surface pressure + ρ_ref g z`） |
| `solid pressure` | 2.38e-3 ✗（见 §6） |

吻合到 1e-10 说明内核重写与 ASPECT 完全一致，也**确证了每个口径究竟用了哪个压力**。

### 3.2 step 0 的初始不平衡 `φ_ASPECT − φ_IC` 及其分解

| 量 | max | rms |
|---|---|---|
| **内核差** `PRM@P_py − PY@P_py` | **2.158e-6** | **1.409e-6** |
| `P_ad − P_py`（静岩廓线的 0.72 MPa 舍入差） | 7.74e-8 | 4.03e-8 |
| **`p_f − P_py`**（= 旧行为引入的不一致） | **6.35e-7** | 8.66e-8 |
| `p_f − P_ad`（新旧行为之差） | **6.04e-7** | 7.90e-8 |
| `p_s − P_ad` | 3.46e-7 | 4.53e-8 |
| ── 合计：`fluid` | 2.452e-6 | 1.377e-6 |
| ── 合计：`adiabatic` | 2.123e-6 | 1.372e-6 |

分解关系（可逐项验证）：

```
φ_adiabatic − φ_IC = [内核差]  + [P_ad − P_py]
                   = -2.1e-6   + ±7.7e-8
φ_fluid     − φ_IC = [内核差]  + [p_f  − P_py]
                   = -2.1e-6   + 6.3e-7
```

即：**压力口径的贡献（~6e-7）被内核差的贡献（~2.1e-6）盖住了**。

### 3.3 压力差与压实波印迹

| 量 | rms | max |
|---|---|---|
| `p_f − P_ad`（旧平衡实际看到的扰动） | 1.115 MPa | 9.246 MPa |
| `p_s − P_ad` | 0.641 MPa | 5.300 MPa |
| `p_c` | 0.475 MPa | 3.938 MPa |
| `φ(fluid) − φ(adiabatic)` | 7.95e-8 | 6.08e-7 |

* `corr( φ(fluid)−φ(adiabatic), p_c ) = +0.987` —— **旧口径确实把压实波结构
  直接搬进了平衡熔融分数**，这正是 `TEQ_PRESSURE_CHOICE.md` 想消除的信号。
* 换成 `adiabatic` 后，平衡压力只剩一个**只随深度变化的常数偏移**
  `P_ad − P_py = −0.72 MPa`，与 `p_c` 无关。
* 相对量级：rms 压力效应 = 波幅 `3e-3` 的 **2.6e-5**。

### 3.4 灵敏度与外推

* `dφ/dP`（在 IC 状态下中心差分）= max **1.08e-13 /Pa**，rms 5.6e-14 /Pa。
* `|Δp| = 9.25 MPa` 对应的熔融曲线位移：
  dunite `+0.118 K`、morb `+0.533 K`、cmorb `−0.245 K`
  （比 `TEQ_PRESSURE_CHOICE.md` 里"< 0.4 K"略大：那里只算了 `p_c`，
  没算动压；`p_f − P_ad` 最大 9.25 MPa 而 `p_c` 最大只有 3.94 MPa）。
* 本算例是弱流动（v ~ 1e-4 m/yr），所以动压小。若将来动压达到
  `1e8 / 1e9 Pa`，同一状态下的平衡熔融分数会偏 `~1e-5 / ~1e-4`
  （波幅的 0.3% / 3%）——那时压力口径才会变成主要矛盾。

---

## 4. 逐条回答

| 问题 | 回答 |
|---|---|
| 改进是否生效？ | **是**。开关可用、默认已切到 `adiabatic`、每个口径用的压力都被逐点确证到 1e-10。 |
| 改进前用 `p_f` 是否让 IC 对 ASPECT 不平衡？ | **是**：旧口径在 IC 处把平衡算在 `P_ad + (p_c + 动压)` 上，与生成 IC 用的静岩压力差 rms 1.11 / max 9.25 MPa，产生 ≤6.1e-7 的孔隙度不一致，且 `corr = +0.987` 地带着压实波印迹。 |
| 这个不平衡严重吗？ | 在本算例里**很小**（波幅的 2.6e-5），既不主导 step 0 的失配，也不产生发散的初始瞬变（3 步内 `φ(fluid) − φ(adiabatic)` 基本恒定在 6.1e-7）。 |
| 改成 `adiabatic` 后 IC 就平衡了吗？ | **没有完全平衡**。残留 ~2.1e-6（rms 1.4e-6），主导项是熔融曲线系数不一致（§5），不是压力。 |

---

## 5. 残留失配的主因：cMORB 系数与 `Surface pressure` 的舍入

IC 生成器 `make_uniform_phase_ic_case6.py` 用的是
`calibrating_June2026.K_class_list_3_components_recalibrated_June2026`：

| 参数 | IC 生成器 | `.prm` | 差 |
|---|---|---|---|
| cMORB `T_m0` | 617.27680612617 | 617.28 | 3.2e-3 K |
| cMORB `A` | 3.799177670844276e-8 | 38.0e-9 | 8.2e-12 |
| cMORB `B` | −4.001653781719964e-18 | −4.00e-18 | 1.7e-21 |
| dunite / morb | 与 .prm 一致 | — | 0 |
| `Surface pressure` | 静岩 P(100 km)=3.060720e9 | 3.06e9 | 7.2e5 Pa |

用 python 内核逐项打开验证（`kernel_attrib.py`，230 km 处 `φ`）：

```
prm(A=38.0e-9, B=-4.00e-18) : 1.04210e-4   ← 与 ASPECT 一致
只把 T_m0 换成生成器值        : 1.04255e-4   Δ=+4.5e-8   （可忽略）
只把 A    换成生成器值        : 1.05066e-4   Δ=+8.6e-7
只把 B    换成生成器值        : 1.05275e-4   Δ=+1.1e-6
三个都换（= 生成器）          : 1.06175e-4   Δ=+2.0e-6
```

全剖面上的 max 失配随系数组合的变化（`attrib_kernel.py`）：

| cMORB 系数 | max&#124;φ − φ_IC&#124; |
|---|---|
| `.prm`（A=38.0e-9, B=−4.00e-18） | 2.158e-6 |
| 只改 `T_m0` | 2.113e-6 |
| 只改 `A` | 1.232e-6 |
| 只改 `B` | 9.686e-7 |
| 改 `A`、`B` | 5.63e-8 |
| 三个都改（= 生成器） | 5.0e-13 |

所以 **2.1e-6 的内核差几乎全部来自 `.prm` 把 cMORB 的 `A`、`B` 舍入成
38.0e-9 / −4.00e-18**：单独改 `A` 剩 1.2e-6，单独改 `B` 剩 9.7e-7，
两个都改就只剩 5.6e-8（那是 `T_m0` 的 0.0032 K 造成的）。

建议（可选，与本改动无关）：把 `.prm` 写成

```
set Surface pressure = 3.06072e9
set Melting point for each component at surface   = 1780, 1000, 617.27680612617
set Melting curve coefficient A for each component = 45.0e-9, 112.0e-9, 3.799177670844276e-8
set Melting curve coefficient B for each component = -2.00e-18, -3.37e-18, -4.001653781719964e-18
```

这样 `adiabatic` 口径下的 step-0 失配会从 ~2.1e-6 掉到 ~1e-11（只剩投影/float32
误差；只改 `A`、`B` 是 5.6e-8，再把 `T_m0` 和 `Surface pressure` 对齐就基本归零），
初始条件真正落在 ASPECT 的平衡态上。

---

## 6. 附注：`solid pressure` 在 t=0 不可用

`solid pressure` 走 `in.pressure[q]`（Stokes 压力）。实测它在 t=0 给出

* step-0 `φ_peak = 6.1e-4`（应为 3.0e-3），`φ_ASPECT − φ_IC` 最大 **2.4e-3**；
* 第一次 Stokes 解出 `p ∈ [−1.77e10, 4.33e10] Pa`、RMS 速度 0.229 m/yr
  （`fluid`/`adiabatic` 是 1.06e-4 m/yr，差 2000 倍）；
* 到 step 1 才恢复正常（`p ∈ [3.06e9, 8.06e9]`，`φ_peak = 2.95e-3`）。

原因是在 t=0 的**第一次**材料模型求值时 `in.pressure` 还不是静岩参考剖面
（而 `p_f` 分量此时已经是静岩剖面），平衡被算在了一个非物理的压力上。
结论：**用外部静岩初始条件时不要选 `solid pressure`**；默认的
`adiabatic pressure` 是正确选择（也与 `melt_simple` 一直以来的做法一致）。

---

## 附：本报告用到的脚本与产物

| 文件 | 说明 |
|---|---|
| `verify_teq_pressure/run.sh` | 跑三个 deck |
| `verify_teq_pressure/analyze.py` | 内核重写 + 分解 + 灵敏度 + 出图 |
| `verify_teq_pressure/attrib_kernel.py` | 把残留失配归因到 cMORB 的哪个系数（§5） |

---

## 附：P0 之后的补充（2026-09-29）

本报告的所有数字都是**打 P0 补丁之前**测的。P0 =
`compute_current_constraints()` 不再给 `prescribed_field` 类型的场施加
Dirichlet 边界条件（详见 `INTERNAL_EQUILIBRIUM_IC_DESIGN.md` §5.10/§9.1）。
它带来的变化是局部的：

* 内部节点（`y ∈ (1 km, 149 km)`）的所有场，以及三个 PDE 场
  `dunite/morb/cmorb` 的**全部**节点，step 0 **逐位不变**；
* 上下边界节点上的 prescribed 场从"初始条件的值"改为"材料模型的平衡值"，
  `max|porosity − Φ(P_ad,T,c̄)|` 在边界上从 **2.133e-06** 降到 **2.106e-12**。

所以本报告 §3 的分解（内核差 2.16e-6 / 压力项 6.0e-7 等）与 §5 的结论
**不受影响**；`analyze.py` 里"去掉上下 1 km 边界节点"的做法现在只是
多余而非必要。验收算例见
`region5_carbon/case6_practical_0724/verify_prescribed_bc/`。
| `verify_teq_pressure/analysis_summary.json` | 本报告所有数字 |
| `verify_teq_pressure/analysis.png` | 四联图 |
| `verify_teq_pressure/out_{fluid,adiabatic,solid}/` | 三个算例输出 |
| `aspect_personal/teq_ab/` | 上一版 12 步 A/B（结论一致：6.09e-7） |
