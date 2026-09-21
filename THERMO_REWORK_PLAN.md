# ASPECT 热力学平衡材料模型 —— 改造指示文档

> **本文档是给"在 `aspect_personal` 工作区新开的对话"的任务书。**
> 新对话没有此前讨论的上下文，所有必要信息都在本文档里。
> 生成时间：2026-09-16。工作区：`/home/zhiqianli/workspace/aspect_personal`。
>
> ⚠️ **关于行号**：本文档引用的 `melt_thermodynamic_equilibrium.cc/.h` 行号是按
> 2026-09-16 的**工作树内容**（`git HEAD dbda05b35` + 使用者自己未提交的两处改动）
> 核对的。动手前请**用引用的代码片段做行号复核**（`grep -n`），
> 不要盲信数字——尤其是 `.h` 的改动会整体平移 `.cc` 以外的行号。
> 相对固定的锚点是**代码片段本身**与**函数名**。

---

## 0. 一句话目标

对 `melt_thermodynamic_equilibrium` 材料模型做四类改造：
**(A) 修一个被参数掩盖的物理 bug**、**(B) 消除 `p_f` 时间异步**、
**(C) 把平衡求解器提速约 24 倍**、**(D) 清理死代码与隐患**。
全部改动都必须做到**结果可验证**（数值对照 + 逐位/量级判据），不得凭感觉改。

---

## 1. 环境与证据

### 1.1 构建

| 项 | 值 |
|---|---|
| 源码根 | `/home/zhiqianli/workspace/aspect_personal` |
| git HEAD | `dbda05b35`（"now TEQ calc is based on p_f"） |
| 工作树现有改动 | `M include/aspect/material_model/melt_thermodynamic_equilibrium.h`、`MM source/material_model/melt_thermodynamic_equilibrium.cc`（**使用者自己的改动，不要动**） |
| 构建目录 | `build_TEQ_bulk_release/`（**Makefile，不是 ninja**；Release） |
| 编译宏 | `ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS=ON` |
| 重编译 | `cd build_TEQ_bulk_release && make -j12`（只改一个 .cc 约 1–2 min） |
| 运行 | `export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 && mpirun -np 20 --oversubscribe ./build_TEQ_bulk_release/aspect <deck>.prm` |

⚠️ `build/`、`build_TEQ_bulk/` 是另外两个构建目录（**没有**开这个宏），不要混用。

### 1.2 证据文档（**动手前先读**）

| 路径 | 内容 |
|---|---|
| `/home/zhiqianli/workspace/melt_equilibrium_comparison/REPORT.md` | 主报告。§9 精度控制、§10 求解效率、§11 `p_f` 异步量化、§12 Picard 实测、§13 机制澄清 |
| `/home/zhiqianli/workspace/melt_equilibrium_comparison/proposed_fix_fill_prescribed.diff` | 已测试过的 P0-2 补丁 |
| `/home/zhiqianli/workspace/melt_equilibrium_comparison/bench_aspect/` | 缩比算例与全部测量记录（含每个 deck 的运行方法） |
| `/home/zhiqianli/workspace/melt_equilibrium_comparison/solver_variants.cpp` | 5 个求解器变体的参考实现（可直接编译对跑） |
| `/home/zhiqianli/workspace/melt_equilibrium_comparison/aspect_kernel_reference.cpp` | 当前热力学内核的**逐行转写**（改造后的正确性基准） |

### 1.3 用于验证的真实算例

| 算例 | 路径 | 特点 |
|---|---|---|
| `case6_05` | `carbon_melting/data/output-melt_TEQ_carbon_case6_05` | 15000 单元 / 150 km 柱 / 6524 步 / A,B≠0（**有压力依赖、有孔隙波**）——主力验证算例 |
| `mr_48_F_100` | `carbon_melting/data/excerpted_output-melt_TEQ_carbon_test_mr_48_F_100` | A=B=0，T 恒定（**∂φ/∂P ≡ 0**）——对照算例，任何 `p_f` 相关改动在此应当**无影响** |

缩比 deck 在 `bench_aspect/`，单个 2–7 min，**优先用它们做 A/B，不要直接跑 6524 步的全量**。

---

## 2. 改动清单

优先级：**P0 = 正确性，必须先做；P1 = 效率大头；P2 = 效率补充；P3 = 结构清理。**

### P0-1 🔴 压实粘度被剪切粘度覆盖（被参数掩盖的物理 bug）

**位置**：`source/material_model/melt_thermodynamic_equilibrium.cc:1408`（ON 分支），
`:2864`（OFF 分支，同样的错）

**现状**：

```cpp
1388:  melt_out->compaction_viscosities[i] = xi_0 * std::exp(- alpha_phi * porosity);
1390-1406:  ... 再乘温度/压力依赖 ...
1408:  // finally, clip viscosity into the range from 1e15 to 1e25
       melt_out->compaction_viscosities[i] = std::max(std::min(out.viscosities[i], 1e25), 1e15);
                                                                      ^^^^^^^^^^^^^^^^^
                                                                      写成了【剪切】粘度
```

注释说 "clip viscosity"，但右边写的是 `out.viscosities[i]`（**剪切**粘度，已在 `:1351` 夹取过），
于是 1388–1406 辛苦算出的压实粘度被**整个丢掉**。
后果：

* `Reference bulk viscosity`（ξ₀）**完全失效**；
* `Thermal bulk viscosity exponent` **完全失效**；
* `compaction_viscosities ≡ viscosities`，而压实长度 δ = √(k·ξ/η_f) 被系统性算错。

**为什么现在没被发现**：两个算例的 `.prm` 都显式设成
`Reference shear viscosity = 1e20` 且 `Reference bulk viscosity = 1e20`，
`Exponential melt weakening factor = 0`，所以 ξ₀·exp(0) ≡ η₀，两种写法数值相同。

**但**：ASPECT 自身的默认值是 **剪切 `5e20`、体粘 `1e22`**（`declare_parameters`，
`.cc:1440-1449`），相差 20 倍。任何**没有显式设这两行**的 deck，ξ 都会偏小 20 倍
⟹ 压实长度偏小 √20 ≈ **4.5 倍**——这对孔隙波的宽度与传播速度是决定性的。

**改法**（一行）：

```cpp
melt_out->compaction_viscosities[i] =
  std::max(std::min(melt_out->compaction_viscosities[i], 1e25), 1e15);
```

**验收**：
1. **现有算例必须逐位/1-ULP 不变**（因 ξ₀=η₀、α_φ=0）。若跑出差异，说明改错了。
2. 造一个 `Reference bulk viscosity = 1e22`（≠1e20）的 deck，确认 `compaction_viscosity`
   输出场变成 1e22（而不是 1e20）。**当前代码下它不会变——这是 bug 的直接证据。**
3. 回顾 `REPORT.md` §11 的孔隙波位置/幅度，确认默认参数下未变。

**风险**：低（默认参数下无影响）。但**需要使用者确认** ξ₀≠η₀ 时的物理意图（见 §6）。

---

### P0-2 🔴 消除 `p_f` 时间异步（末尾补填）

**问题**：`Simulator::fill_prescribed_fields()` 排在 Stokes 求解**之前**，
所以存下来的 `porosity` 用的是**上一时间步**的流体压力 `p_f`，
而后处理假定 `porosity = φ^eq(p_f, T, c̄)`（同一时间步）。
实测残差：`case6_05` 为 8.4×10⁻⁹（相对 2.8×10⁻⁶），是后处理 Γ 的主要误差源。

**机制**（细节见 `REPORT.md` §13）：`p_f` 是**独立求解变量**
（`source/simulator/melt.cc:1696-1701` 插入，与 `p_c` 同 block），
**平流不写它**，只有 `assemble_and_solve_stokes()` 才更新。
填充语句永远排在 Stokes 前 ⟹ 存储的 φ 永远比存储的 p_f 晚一次 Stokes 求解。
只跑一遍（= 现状）时，"一次 Stokes 求解"就等于"一整个时间步"。

**改法**：在 Stokes 之后**再补一次** `fill_prescribed_fields()`。两处，共 9 行：

```cpp
// source/simulator/solver_schemes.cc:817 附近
 void Simulator<dim>::solve_single_advection_single_stokes ()
 {
   assemble_and_solve_temperature();
   assemble_and_solve_composition();
   fill_prescribed_fields();
   assemble_and_solve_stokes();
+  // 用刚算出的 p_f 再刷一次 prescribed 场，
+  // 使存储的 porosity 与存储的 p_f 属于同一个状态
+  fill_prescribed_fields();
   ...
 }

// source/simulator/solver_schemes.cc:1221 附近
// （solve_iterated_advection_and_stokes 的 while 之后）
   while (nonlinear_solver_control.check(...) == SolverControl::iterate);
+  fill_prescribed_fields();
   signals.post_nonlinear_solver(nonlinear_solver_control);
```

现成补丁：`/home/zhiqianli/workspace/melt_equilibrium_comparison/proposed_fix_fill_prescribed.diff`
（`git apply` 即可，已测试过）。

**验收**（判据是"ε_now ≪ ε_lag"）：

设第 n 步
`ε_now = max| φ^eq(p_f(n), T(n), c̄(n)) − φ_stored(n) |`、
`ε_lag = max| φ^eq(p_f(n−1), T(n), c̄(n)) − φ_stored(n) |`。
（`melt_equilibrium_comparison/picard_final.py` 里有现成实现，可搬。）

| 算例 | 改前 ε_now | 改后 ε_now 期望 |
|---|---|---|
| `case6_05` 正常步 | 8.2–8.5×10⁻⁹ | **≤ 0.7×10⁻⁹**（float32 存储底噪） |
| `mr_48_F_100` | ~0.6×10⁻⁹ | **不变**（∂φ/∂P ≡ 0，应无影响） |

并确认状态场**逐位/1-ULP 不变**：`T` 应逐位相同；`p_f`、`p`、`dunite` 差 ≤1 ULP
（此前实测：`T` 差 0.000e+00；`porosity` 系统性改变 8.2×10⁻⁹，正是要修的那一项）。

**明确不要做的事**：不要用"开启 Picard 迭代"代替本改动。
实测 `iterated Advection and Stokes, tol=1e-5` **在"第一次迭代即达标"的时间步上会静默失效**
（`REPORT.md` §12.2 步 7：ε_now 回到 4.3×10⁻⁹，与不改一样），
而末尾补填**与迭代次数无关**。若同时想要 Picard 的物理自洽，两者叠加即可
（`REPORT.md` §12.5 方案 B′：代价 2.4×）。

**风险**：低。已逐位验证物理中性。

---

### P1-1 🟠 二分法缓存 `f(a)/f(b)/f(c)`（结果逐位不变）

**位置**：`include/aspect/material_model/melt_thermodynamic_equilibrium.h:416-460`（ON 分支）

**现状**：循环体内 `f(a)`、`f(b)` 各算 1 遍，`f(c)` 在「等于 0」「与 `f(a)` 比符号」
「与 tolerance 比」三处**各算 1 遍** ⟹ **每次迭代 6 次函数求值**；
又因参数是 `std::function`，编译器无法做公共子表达式消除。
实测：一次熔体分数求解要 ~252 次求值、~1281 次 `exp()` 调用。

**改法**：

```cpp
double
bisection (const std::function<double(const double)> &f,
           const double lower_bound, const double upper_bound,
           const unsigned int max_iter  = 1000,
           const double       tolerance = 1e-10) const
{
  double a = lower_bound, b = upper_bound;
  const double fa0 = f(a);              // ← 移出循环
  const double fb0 = f(b);
  AssertThrow(fa0 * fb0 <= 0,
              ExcMessage("bisection: the function must have opposite signs at the bounds. "
                         "a=" + std::to_string(a) + " b=" + std::to_string(b)));
  double fa = fa0, fb = fb0;
  double c = 0.0, fc = 0.0;
  for (unsigned int i = 0; i < max_iter; ++i)
    {
      c  = 0.5 * (a + b);
      fc = f(c);                         // ← 每次迭代只求 1 次
      if (fc == 0.0) break;
      if (fc * fa < 0) { b = c; fb = fc; }
      else             { a = c; fa = fc; }
      if (std::fabs(fc) < tolerance) break;
    }
  (void)fb;
  return c;
}
```

中点序列与原实现**完全一致** ⟹ 结果逐位相同。

**验收**：`bench_aspect/scale_1600.prm` 改造前后 `porosity`/`p_f`/`T` **逐位相同**；
`Interpolate prescribed composition` 一节耗时降到约 1/2。
（参考 `REPORT.md` §10.1 表 B 行：exp 调用 1281→228。）

---

### P1-2 🟠 用端点符号判据取代固/液相线两次二分（exp 少 427 倍）

**位置**：`source/material_model/melt_thermodynamic_equilibrium.cc:499-571`

**依据**：`.cc:550-556` 的原注释说"这个方法还没实现，因为要先确认 `f_eq_equation`
是否单调"。**可以证明它严格单调增**：

```
F(f) = Σ cᵢ(Kᵢ−1)/(f+(1−f)Kᵢ)
dF/df = Σ cᵢ(Kᵢ−1)²/(f+(1−f)Kᵢ)² > 0
F(0) = σ − Σ cᵢ/Kᵢ,   F(1) = Σ cᵢKᵢ − σ,   σ = Σ cᵢ
```

且 `find_solidus` 解的是 `Σ cᵢ/Kᵢ = 1`、`find_liquidus` 解的是 `Σ cᵢKᵢ = 1`，
正是 `F(0)`、`F(1)` 的零点条件。因为 `Σcᵢ/Kᵢ` 随 T 增、`ΣcᵢKᵢ` 随 T 减：

```
T < T_solidus  ⟺  F(0) > σ − 1
T > T_liquidus ⟺  F(1) < 1 − σ
```

**改法**：删掉 `.cc:499-571` 的两次 `find_solidus`/`find_liquidus` 调用，替换为

```cpp
double sigma = 0.0, F0 = 0.0, F1 = 0.0;
for (unsigned int i = 0; i < n_components; ++i)
  {
    sigma += bulk_concentrations[i];
    F0    += bulk_concentrations[i] * (1.0 - 1.0 / equilibrium_constants[i]);
    F1    += bulk_concentrations[i] * (equilibrium_constants[i] - 1.0);
  }
if (F0 - (sigma - 1.0) > 0.0) return 0.0;   // 全固
if (F1 - (1.0 - sigma) < 0.0) return 1.0;   // 全液
// 否则 [0,1] 天然构成括号（F 单调增，F(0)<0<F(1)），直接对 f 求根
```

`equilibrium_constant` 已有 ±700 指数夹取，保证 `1.0/Kᵢ` 不除零。

⚠️ `find_solidus`/`find_liquidus` 仍在头文件里声明。确认没有其它调用点后再删
（或保留但标注"仅供诊断"）。`solidus`/`liquidus` 目前**只用于错误信息**，
`evaluate()` 不读它们——要保留诊断就在出错时才懒计算。

**验收**：
- `REPORT.md` §3.1 的 101 520 点网格上，改造前后 `f` 差 ≤ 1e-9；
- `bench_aspect` 的 `porosity` 相对差 ≤ 1e-12；
- `exp()` 调用从 ~228 降到 **3**（可用 `solver_variants.cpp` 的计数方式复核）；
- 熔体分数精度**不应变差**（`REPORT.md` §10.1 表 E 行：保护牛顿下中位误差为 0）。

---

### P1-3 🟠 去掉每次平衡求解的 400+ 次堆分配

**位置**：
- `include/.../melt_thermodynamic_equilibrium.h:392-394` + `.cc:383-385`：
  `solve_eq_melt_fraction(..., std::vector<double> bulk_concentrations)` —— **按值**收 vector
- `.cc:499-502`：进函数后又拷 4 份（`_melting_points` / `_bulk_concentrations` /
  `_latent_heats` / `_tuning_parameters`）
- `.cc:275`（`find_solidus`）与 `.cc:336`（`find_liquidus`）：lambda **每次求值**都
  `std::vector<double> eq_consts(n)` —— 一次固相线二分 ≈210 次求值 = **210 次 malloc**

**改法**：
1. 签名改成 `const std::vector<double> &bulk_concentrations`（头文件 + .cc 同步）；
2. 删掉 `.cc:499-502` 的 4 份拷贝，直接用原变量（它们本来就没被修改）；
3. `eq_consts` 移到 lambda **外面**声明一次复用；`n_components` 只有 3，
   用 `std::array<double, N>` 或固定长度栈数组更好；
4. 可选：把 `bisection` 的 `std::function` 参数改成模板 `template <typename F>`，
   才能真正内联（这会改头文件接口，注意编译面）。

**验收**：结果逐位不变；`case6_05` 的 `Interpolate prescribed composition`
一节耗时降到约 1/2（`REPORT.md` §10.1 表 A→B：18.3→8.5 µs/点）。

---

### P2-1 🟡 扰动求解在未用潜热时是死代码（2×）

**位置**：`.cc:1106-1175`

**现状**：对每个点**再解一次** `T + 1e-3` 的平衡（含它自己的两次固/液相线二分），
结果只写进 `ComponentPhaseExchangeOutputs`（`.cc:1214-1219`）。
而 `create_additional_named_outputs()`（`.cc:1880-1899`）**从不创建**该输出对象
——它只由 `heating_model/component_phase_exchange_heating.cc:193` 创建，
即**只有启用 `component phase exchange heating` 时才存在**。
两个算例的 `Heating model` 段都是注释掉的 ⟹ `phase_exchange_out == nullptr` ⟹ **白算一半**。

**改法（最小，一行守卫）**：

```cpp
if (phase_exchange_out != nullptr && enable_equilibrium_calculation)
  { ... 原 .cc:1100-1175 ... }
```

**改法（更好，将来真要用潜热时）**：改用**解析导数**，不要差分。
记 `Dᵢ = f + (1−f)Kᵢ`、`Kᵢ' = dKᵢ/dT = −Kᵢ(Lᵢ/rᵢ)/T²`：

```
∂F/∂f =  Σ cᵢ(Kᵢ−1)²/Dᵢ²                (> 0)
∂F/∂T = −(1/T²) Σ cᵢKᵢ(Lᵢ/rᵢ)/Dᵢ²
df/dT = −(∂F/∂T)/(∂F/∂f)
      = (1/T²)·[Σ cᵢKᵢ(Lᵢ/rᵢ)/Dᵢ²] / [Σ cᵢ(Kᵢ−1)²/Dᵢ²]

dc_l,ᵢ/dT = [−cᵢKᵢ'(1−f) + cᵢ(Kᵢ−1)·df/dT] / Dᵢ²
dc_s,ᵢ/dT = [ cᵢKᵢ'f      + cᵢKᵢ(Kᵢ−1)·df/dT] / Dᵢ²
```

不需要额外 `exp`、不需要第二次求解，且比 `1e-3` 差分**更准**（无 O(ΔT²) 截断误差）。

**验收**：`.prm` 未启用该加热模型时，`case6_05` 总墙钟降到约 **1/2**；
`porosity`/`p_f`/`T` **逐位不变**（因为结果本就被丢弃）。

---

### P2-2 🟡 `evaluate()` 里重复计算 T_m 与 K

**位置**：`.cc:964-988`（求解器内部已算过一遍，这里又算 3 次 `exp`）、
`.cc:1244-1290`（`debug_1..3` 再各算一次 T_m 与 K，共 3 次 `exp`）

**改法**：让 `solve_eq_melt_fraction` 通过出参返回 `melting_points` 与
`equilibrium_constants`（或把它们提到 `evaluate()` 里算一次传进去）。
`debug_*` 字段的值就是 K，可直接复用。

**验收**：结果逐位不变。P1-2 完成后这一项占比会明显上升，值得顺手做。
`Fill debug fields = false` 时（`mr_48_F` 算例）应完全跳过该块。

---

### P2-3 🟡 8 个 prescribed 场 = 8 次完整 `evaluate()`，只有 1 次必要（8×）

**位置**：`source/simulator/solver_schemes.cc` 的 `fill_prescribed_fields()`

**现状**：算例里 11 个组分场有 **8 个**是 `prescribed field`
（`porosity, dunite_liquid, morb_liquid, cmorb_liquid, melting_rate, debug_1..3`）。
`fill_prescribed_fields()` 对**每个** prescribed 场各调一次
`interpolate_material_output_into_advection_field()`，每次都是一遍**完整网格**
+ `material_model->evaluate()`。而 `evaluate()` 无论被谁调用都会解平衡
（守卫 `in.requests_property(...)` 在 `MaterialModelInputs` 默认
`requested_properties = all_properties` 下恒为真，见 `interface.h:263/281/326`），
且**一次 `evaluate()` 已经把全部组分场的 prescribed 值都填好了**
（`.cc:881` / `:940` / `:1085`）。

**证据**：`case6_05` 的 `log.txt` 里 `Interpolate prescribed composition` 调用数
**52192 = 8 × 6524 时间步**，恰好等于 prescribed 字段数；该节占总墙钟 **29%**。
网格缩放测试证明该节成本 **100% 在逐点上**（无固定开销）。

**改法**：把 `fill_prescribed_fields()` 改成**一遍网格循环填所有 prescribed 场**。
本算例所有组分场都是 `Composition polynomial degree = 2`，支撑点完全重合，一遍即可。
实现要点：先收集所有 prescribed 场的 index（若阶数不一致则按阶数分组，分几遍），
在单元循环里**只调一次** `material_model->evaluate(in, out)`，
再把 `prescribed_field_out->prescribed_field_outputs[i][c]` 一次性写入
`distributed_vector` 的各个 block。

**验收**：`Interpolate prescribed composition` 调用数从 8N 降到 N（N = 时间步数）；
该节耗时降到约 1/8；**`porosity`/`p_f`/`T` 逐位不变**。

> ⚠️ 这是**唯一需要改 ASPECT core（而非材料模型）** 的改动，改动面最大。
> 请放在 P0/P1 全部通过之后再做，并单独 commit 便于回退。

---

### P3 结构清理（低风险，随手做）

| # | 位置 | 问题 | 改法 |
|---|---|---|---|
| P3-1 | `.cc:1684` + `.h` | `enable_chemical_reaction_rate` 的解析被注释掉，成员**无默认初始化**（该头文件所有成员都没有），ON 分支从未读取 ⟹ 未初始化 bool。`.prm` 里 `Enable chemical reaction rate = true` 被静默忽略 | 取消注释真正解析，或从 `.h` 删除该成员与 `.cc:1508` 的 `declare_entry`，或至少 `= false` 默认初始化 |
| P3-2 | `.cc:314`、`:375`（ON）、`:2106`、`:2167`（OFF） | `AssertThrow(false, ExcMessage(""))` —— **空错误信息**，出错时无法定位 | 填上有用信息（压力、浓度、括号端点；可参考 `.cc:521-548` 已有的写法） |
| P3-3 | `.cc:36` | `const bool local_debug = true;` 所有使用都被注释 ⟹ 死变量（unused 警告） | 删除 |
| P3-4 | `.cc:1904-3323` | OFF 分支是**死代码且编译不过**（`MeltOutputs::concentrations_in_phases` 只在宏开启时声明）。另有 3 处与 ON/Python 不一致：无压力阈值外推；越界直接 `AssertThrow` 而非夹取；体成分从分相场重构 | 整段删除，或至少在文件头写明"此分支不可编译、不要维护" |
| P3-5 | `.cc:1375`、`:2832` | `* Utilities::fixed_power<0>(1.0-(porosity+background_porosity))` —— `fixed_power<0>` 恒等于 1.0。注释说明是刻意砍掉 `(1−φ)²`，但代码形态极易误读 | 直接删掉该因子，把理由留成注释 |
| P3-6 | `.h:114-141`、`.cc:1856-1858` | `CompositionFieldGradientsInputs::fill` 在 FEValues 缺 `update_gradients` 时**静默 return**，留下零梯度；守卫是"size == n_evaluation_points"，极易忽略 | 加 `AssertThrow` 或显式标记"未填充"，避免静默零梯度 |
| P3-7 | `.cc:102-189` | `match_component_index_to_composition_index` 不对称：bulk 缺失 → `AssertThrow`；**solid 缺失 → 静默容忍**；liquid 缺失 → `AssertThrow` | 明确意图后统一；若 solid 可选是有意设计，加注释说明 |

---

## 3. 建议的实施顺序

```
第 0 步  建立基线：跑 bench_aspect/scale_1600.prm，存下 porosity/p_f/T 的
        逐位指纹（后续都用 np.array_equal 对照），以及 log.txt 的计时表
第 1 步  P0-1 压实粘度          ← 一行，先做，确认默认参数下逐位不变
第 2 步  P0-2 末尾补填          ← 两处 9 行，用 ε_now/ε_lag 判据验收
第 3 步  P1-3 去堆分配          ← 逐位不变
第 4 步  P1-1 二分缓存          ← 逐位不变
第 5 步  P1-2 端点符号判据      ← 允许 ≤1e-9 的差异
第 6 步  P2-1 扰动求解守卫      ← 逐位不变，2×
第 7 步  P2-2 复用 T_m/K
第 8 步  P3-* 清理
第 9 步  P2-3 合并 prescribed 填充（改 core，单独 commit）
```

**每一步都要单独 commit**，commit message 写清本文档编号（如 `P1-1 bisection cache`），
便于出现回归时二分定位。

---

## 4. 总体验收

改造完成后应能复现以下数字（来源：`REPORT.md` §10.1、§12）：

| 指标 | 改造前 | 改造后期望 |
|---|---|---|
| 单次平衡求解 `exp()` 调用 | 1281 | **3** |
| 单次平衡求解目标函数求值 | 678 | ≤ 50 |
| `case6_05` ε_now（φ 与 p_f 自洽性） | 8.4×10⁻⁹ | **≤ 0.7×10⁻⁹** |
| `mr_48_F_100` ε_now | 0.6×10⁻⁹ | 不变 |
| `Interpolate prescribed composition` 占总墙钟 | 29–33% | **≤ 5%** |
| 平衡内核单点耗时 | ~18 µs | **≤ 1 µs** |
| 熔体分数 f 对 `aspect_kernel_reference.cpp` 的偏差 | — | ≤ 1e-9 |

---

## 5. 禁区

1. **不要改物理参数与 `.prm` 的默认值**，除非本文档明确指出（P0-1 验收需要临时造 deck）。
2. **不要为了让测试通过而放宽判据**。若某改动达不到期望，先回报，不要自行降低标准。
3. **不要动使用者现有的两处未提交改动**
   （`include/.../melt_thermodynamic_equilibrium.h`、`source/.../melt_thermodynamic_equilibrium.cc`
   的 `M`/`MM` 状态）。改之前先确认 `git status`，或先 `git stash`。
4. **不要用"开启 Picard 迭代"替代 P0-2**（`REPORT.md` §12.2 已证其会静默失效）。
5. **不要删除或覆盖 `/home/zhiqianli/workspace/melt_equilibrium_comparison/` 下的任何东西**
   —— 那是证据与基准，只读引用。
6. 改 `source/simulator/` 之前先确认 `git status` 干净；回退用
   `git checkout -- <file>`（**不要** `git checkout .`，会连带回退使用者自己的改动）。

---

## 6. 需要使用者拍板的问题

1. **P0-1**：ξ₀ ≠ η₀、α_φ ≠ 0 时，压实粘度应当用哪个表达式？
   现表达式为 `ξ = ξ₀·exp(−α_φ·φ)·(温度依赖)`。默认参数下三种写法等价，
   无法从现有结果反推意图。
2. **P2-3**：是否接受改动 ASPECT core（`solver_schemes.cc`）？
   若不想动 core，这一项（8× 收益）只能放弃。
3. **P3-4**：OFF 分支直接删除，还是保留但标注"不可编译"？
