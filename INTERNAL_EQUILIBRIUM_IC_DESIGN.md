# ASPECT 内置「平衡初始条件」模块 —— 设计方案

> 状态：**设计稿，尚未动手改代码**。本文只回答"怎么做、为什么这么做、代价多大"，
> 所有结论都已核到源码行号或实测数字。
>
> 结论先行：
>
> 1. **最大的一条简化**：在这个模型里**真正需要初始条件的只有 `dunite / morb / cmorb`
>    三个场**。`porosity`、`*_liquid`、`melting_rate`、`debug_*` 全是 `prescribed field`，
>    step 0 的第一次 `fill_prescribed_fields()` 发生在**第一次 Stokes 之前**，所以它们
>    的初值从来不被使用（只有边界节点例外，见 §5.10）。于是"我只指定几个场、让 ASPECT
>    算出其他场"这件事，本质上是**给 c̄ 一个由平衡反演出来的初值**，其余全是免费附赠的。
> 2. **模块形态**：一个新的 `InitialComposition` 插件 `thermodynamic equilibrium`，
>    它调用**材料模型自己的内核**（新增一个 ~30 行的 API），在网格节点上把
>    `(P, T, 目标) → c̄ / φ / c_s / c_l` 算出来。三种用法（前向 / 反演熔融分数 /
>    指定相成分）统一成"**一维混合族 + 一个标量**"。
> 3. **精度**：我核对了 `source/simulator/helper_functions.cc:1994-2157`——
>    prescribed 场的**求值点就是单元的 support points（= 节点）**，
>    不是求积点。所以只要 IC 也在节点上、用同一个内核，
>    `φ_stored(node) ≡ Φ(P,T,c̄)` 到二分容差（~1e-11），
>    你列的四个误差源里 **参数舍入 / ASCII 插值 / float32 三个直接归零**。
> 4. **附带一个独立小改动**：`compute_current_constraints()` 里**跳过 `prescribed_field`
>    类型的场**，不再给它们施加 Dirichlet 边界条件。你在
>    `melt_rate_pure_convection_diag` 分支的 `6f6830314` 上已经为 `*_pure_convection`
>    做过同样的事（后缀硬编码），这里把它推广成一条通用规则。

---

## 0. 版本前提

* 分支/工作树：`aspect_personal`，HEAD = `390046b3b`，`Pressure for thermodynamic equilibrium`
  的改动**尚未提交**（工作区修改）。
* 编译：`ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS=ON`；`CMakeLists.txt:47` 是
  `FILE(GLOB_RECURSE TARGET_SRC ...)`，所以**新增源文件不需要改 CMake**。
* 下面引用的行号都以当前工作树为准。

---

## 1. 现状盘点：ASPECT 侧到底发生了什么

### 1.1 场的分类（case6 / region5 的 11 个场）

| 场 | method | 谁决定它的值 | 需要初值吗 |
|---|---|---|---|
| `dunite`, `morb`, `cmorb` | `field` | **平流 PDE**（固相速度，装配里含 `divergence_u` 与液相分离项，无源项） | **需要**，这是唯一的自由数据 |
| `porosity` | `prescribed field` | 材料模型每步 `φ = Φ(P,T,c̄)`（`.cc:845-847`） | 不需要 |
| `dunite_liquid` … | `prescribed field` | 材料模型每步填 | 不需要 |
| `melting_rate` | `prescribed field` | `.cc:1067-1078`，且 `timestep_number()>0` 才填 | 不需要 |
| `debug_1..3` | `prescribed field` | `Fill debug fields=true` 时填 K₀,K₁,K₂ | 不需要 |

### 1.2 prescribed 场在哪里、以**什么点**求值（关键）

`Simulator::fill_prescribed_fields()`（`solver_schemes.cc:372-416`）
→ `interpolate_material_outputs_into_advection_fields()`（`helper_functions.cc:1994-2157`）：

```cpp
const Quadrature<dim> quadrature(
    dof_handler.get_fe().base_element(...).get_unit_support_points());   // ← 支撑点 = 节点
FEValues<dim> fe_values (*mapping, dof_handler.get_fe(), quadrature, ...);
...
in.reinit(fe_values, cell, introspection, solution);
material_model->evaluate(in, out);
...
distributed_vector(local_dof_indices[dof_idx])
    = prescribed_field_out->prescribed_field_outputs[j][adv_field.compositional_variable];
...
if (adv_field.is_temperature() || ... != density)
  current_constraints.distribute (distributed_vector);   // ← 2148-2150：Dirichlet 覆盖
```

三条结论：

1. **材料模型是在节点上被调用的**（支撑点），不是单元求积点。所以
   `φ_stored(node) = Φ(P(node), T(node), c̄(node))` **逐点成立**。
   → 在节点上生成 IC 是**精确**的，不存在节点/求积点插值不一致。
   （文档里说的"求积点 vs 节点 ε"是 `melting rate instantaneous` 插件在 **cell midpoint**
   上求值跟节点存储值比，属于诊断口径问题，不是 IC 问题。）
2. `current_constraints.distribute()` 会把 Dirichlet 边界值**盖在** prescribed 场上面
   → 这正是我们上一轮实测到的"边界节点 φ 等于 IC 值、内部等于平衡值"。
3. 调用时机（`solve_single_advection_single_stokes`，`solver_schemes.cc:817-827`）：
   `assemble_and_solve_temperature → assemble_and_solve_composition → fill_prescribed_fields
   → assemble_and_solve_stokes → fill_prescribed_fields`。
   **step 0 时 dt=0，前两个 solve 是空转，所以第一次平衡重写早于第一次 Stokes。**
   IC 里的 φ/c_l 只会被用在：step-0 之前的 postprocess，以及（在没修 §5.10 之前）
   被 Dirichlet 锁住的边界节点。

### 1.3 初始条件设置的顺序与可用的 hook

`Simulator::run()`（`core.cc:1997-2002`）：

```cpp
signals.pre_set_initial_state (triangulation);
set_initial_temperature_and_compositional_fields ();
compute_initial_pressure_field ();
signals.post_set_initial_state (*this);
```

* `set_initial_temperature_and_compositional_fields()`（`initial_conditions.cc:42-252`）：
  对 `n = 0..n_fields`，用 `VectorTools::interpolate` 把
  `initial_temperature_manager->initial_temperature(p)` /
  `initial_composition_manager->initial_composition(p, n-1)` 插到支撑点上；
  最后 `compute_current_constraints(); current_constraints.distribute(initial_solution);`
  → **边界节点立刻被 BC 覆盖**。
* `compute_initial_pressure_field()`（`initial_conditions.cc:439+`）：把
  `adiabatic_conditions->pressure(p)` 插到 **`fluid pressure`** 变量上。
  → **t=0 的 `p_f` 就等于绝热（静岩参考）剖面**，与新的默认压力口径完全一致。
* `adiabatic_conditions` 在 IC 之前就由 `ComputeProfile::initialize()`
  （`compute_profile.cc:58-192`，内部会调用材料模型）建好了。
* 插件可用的东西（`simulator_access.h`）：
  `get_initial_temperature_manager()`(649)、`get_adiabatic_conditions()`(552)、
  `get_material_model()`(526)、`get_signals()`、`get_dof_handler()`、`get_mapping()`、
  `get_solution()`。
  ⚠️ 两个 IC manager 在 **step 0 解完之后**被 `reset()`（`core.cc:2042-2043`，
  在 `start_timestep()+solve_timestep()` 之后）。所以
  `post_set_initial_state`（在 2002 行、主循环之前）里**还能**问温度 IC，
  但插件如果想把状态留到更晚用，必须自己缓存。
* `signals.post_set_initial_state` 的签名是 **`void(const SimulatorAccess<dim>&)`**
  （`simulator_signals.h:92`）⟹ 只能**读**、不能改 `solution`。
  这正是我把汇总报告挂在它上面、而不把生成逻辑挂上去的原因；
  若将来要做"能改 T 也能改 c̄"的后处理式模块（§4 方案 C），需要另加 hook。
* 插件范式（照抄 `initial_composition/ascii_data.h`）：
  ```cpp
  template <int dim> class X : public Interface<dim>, public aspect::SimulatorAccess<dim>
  ```
  manager 会 `initialize_simulator()` → `parse_parameters()` → `initialize()`。
* `Manager::initial_composition()`（`initial_composition/interface.cc:159-172`）把列表里的
  插件按 `List of model operators`（默认 `add`）**累加**，初值 0。
  → 用这个新插件时应当让它**单独**出现在 `List of model names` 里。

---

## 2. 现状盘点：你现在是怎么造 IC 的（R1–R10）

| 编号 | 位置 | 指定 | 求出 | 新模块要覆盖吗 |
|---|---|---|---|---|
| R1 | `make_initial_conditions.py:27` `calculate_equilibrium_array` | P, T, c̄ | φ, c_s, c_l | ✅ 前向模式 |
| R2 | `make_uniform_phase_ic.py` | 波形 f(x)、常数 K | c̄=(1−f)c_s,bg+f c_l,bg | ✅ tie-line 族 |
| R3 | `make_uniform_phase_ic_case6.py` | 深度相关 P/T、f_bg(z)、波形、200 ppm 背景 | 同上（逐深度） | ✅ **主要目标** |
| R4 | `make_initial_conditions.py:1722` | f_target(z)、方向 (0.7, 0.3−c, c) | brentq 反解 c_cmorb | ✅ 反演模式 |
| R5 | `make_initial_conditions.py:1318` | 孤立波 f，P/T 常数 | 闭式反解 c̄ | ✅ 反演模式（闭式是特例） |
| R6 | `make_initial_conditions.py:390` `calculate_T_from_f` | f, c̄ | **T** | ⚠️ v1.1（自由 DOF = T） |
| R7 | `make_initial_conditions.py:767` `isoenthalpy_melting_profile` | 等焓 H₀ | (T, c_an) 2×2 | ❌ v2 |
| R8 | `make_initial_conditions.py:504` | c̄ 剖面 | φ | ✅ 前向模式 |
| R9 | `make_initial_conditions.py:246/327` | f（指定 T 用 degC 传参） | 闭式 c̄（2/3 组分） | ✅ 反演模式的特例 |
| R10 | `make_initial_conditions.py:583` | 沿一条 δc 方向扫 c̄ | φ 曲线 | ✅ 反演模式 |

**顺带确认的三处坑**（设计要避开）：

1. **同名 K 列表有两份、内容不同。** `calibrating_June2026.K_class_list_3_components_recalibrated_June2026`
   = `617.27680612617 / 3.799177670844276e-8 / −4.001653781719964e-18`（`..._uniphase.txt` 用它）；
   `component_data.K_class_list_3_components_recalibrated_June2026` =
   `617 / 38.0e-9 / −4.00e-18`（= 现在 `.prm` 的值，`initial_composition_for_carbon_case06.txt` 用它）。
   两份差 ⟹ φ 差 **2.04e-6**（上一轮已实测）。
2. `calc_compaction_length` 算出 **683.13 m**，但两个 case6 脚本都硬编码 `wave_delta = 53.45`。
3. 闭式反解 R9 在硬编码 T 下给出非物理 `c_an ≈ 26.5`（无适用域检查，`main()` 里已注释掉）。

**ASCII 格式**（如果还要保留"导出对照"能力就得遵守）：
`# POINTS: 2 <n_y>` 必须有；行序 x 最快；列序 = `Names of fields`；
y 从盒底算、升序；`repr(float(...))` 全精度；否则 ASPECT 要么 `AssertThrow` 要么**静默退出**。

---

## 3. 把需求形式化：到底有几个自由度

每个节点上的未知标量：`P, T, c̄(3 个但 Σ=1 ⇒ 2 个自由), φ, c_s(2 个自由), c_l(2 个自由)` = **8**。
平衡关系：`(φ, c_s, c_l) = Φ(P, T, c̄)` ⇒ **5** 个约束（φ 1 个 + c_s 2 个 + c_l 2 个）。

⟹ **剩 3 个自由度**。

* 固定 **P**（= 绝热剖面，也就是 t=0 的 `p_f`）→ 剩 2；
* 固定 **T**（= 你现有的 `Initial temperature model`，完全自由）→ 剩 **1**。

所以"指定某几个场、其余由平衡决定"这件事，本质上永远是：

> **在一条一维族上选一个点。**

统一形式：

```
c̄(s; x) = (1 − s)·A(x) + s·B(x)          Σ A = Σ B = 1  ⇒  Σ c̄ = 1 恒成立
Φ(x, s) ≡ Φ(P(x), T(x), c̄(s;x))          s ↦ Φ 单调增（下面给证明/条件）
```

三种用法 = 这个标量 `s` 分别取：

| 用法 | 你给什么 | 怎么定 s | 备注 |
|---|---|---|---|
| **前向** | c̄ 本身（uniform / 剖面 / 函数） | 不需要 s，直接用 c̄ | 你的 R1/R8/R10；其实**今天就已可用**（见 §8） |
| **反演熔融分数** | φ*(x) 目标 | 解 `Φ(x,s) = φ*(x)` | 你的 R2/R3/R4/R5/R9/R10，**核心新增能力** |
| **指定相成分** | φ*, c_s*, c_l* | `c̄ = (1−φ*)c_s* + φ*c_l*`（质量平衡，闭式） | 过定 ⇒ 必须**校验**是否落在模型 tie-line 上 |

---

## 4. 设计空间与选型

| 方案 | 做法 | 改动面 | 自由度 | 评价 |
|---|---|---|---|---|
| **A. 新 `InitialComposition` 插件** | 新类 + 材料模型小 API | 中（1 新插件 + ~30 行 API） | 写 c̄/φ/c_s/c_l 全覆盖 | ✅ **推荐** |
| B. 材料模型加开关 | `set Fill initial equilibrium = true`，simulator 在 IC 后调用一遍 | 小 | 只能改材料模型知道的场；规格参数要挤进材料模型 subsection | 可行但耦合差，材料模型不该知道"孤立波" |
| C. 新插件类型 / `post_set_initial_state` 后处理 | 新 hook（该 signal 的签名是 `const SimulatorAccess&`，**不能改状态**） | 大 | 最大（可同时改 T 和 c̄） | 留作 v2；只有"自由 DOF = T"才真正需要 |
| D. 只修边界条件，不做反演 | 零新代码 | 极小 | 只能前向 | 作为 **Phase 0**，立刻可做 |

选 **A**，理由是：
* `initial_composition(p, n_comp)` 的签名虽然只能给一个场，但插件知道 `n_comp`，
  而且 T/P/材料模型都能通过 `SimulatorAccess` 拿到，信息是够的；
* `Manager` 的 operator 机制允许"这个插件单独负责所有场"；
* 不侵入材料模型的热路径；
* 与你已有的 `melting rate instantaneous` 插件是同一套"把外部 python 逻辑搬进内核"的思路。

---

## 5. 推荐设计

### 5.1 数据流

```
                      ┌──────────────────────────────────────────────┐
Initial temperature ──┤ initial_temperature_manager->initial_temperature(p) │──┐
  model (随便你用什么) └──────────────────────────────────────────────┘  │
                                                                        ▼
AdiabaticConditions ──► P(x) ──────────────────────────────►  ┌──────────────────────┐
  (Surface pressure + ∫ρ_ref g dz)                            │  ThermodynamicEquil.  │
                                                              │  InitialComposition   │
Background bulk  ──► c̄_bg(x) ──► Φ(P,T,c̄_bg) ──► A,B(x) ───► │  1-D family + 1 scalar│
  (uniform/剖面/函数)                                          │  solve s              │
                                                              └──────────┬───────────┘
Target (φ* 或 c_s*,c_l*) ──────────────────────────────────────────────┘
                                                                        │
                              ┌─────────────────────────────────────────┴──────────┐
                              ▼                    ▼              ▼                ▼
                        c̄ (bulk fields)      φ (porosity)    c_l (*_liquid)   c_s (若存在 *_solid)
```

### 5.2 材料模型要加的 API（约 30 行）

`MaterialModel::Interface`（`include/aspect/material_model/interface.h`）加两个**有默认实现**的虚函数，
`MeltThermodynamicEquilibrium` 覆写：

```cpp
/** 一个"固-液热力学平衡"模型的一个化学组分。 */
struct EquilibriumComponent
{
  std::string name;                  // "dunite"
  std::string solid_field_suffix;    // "_solid"（可配置！见 .cc:1637）
  std::string liquid_field_suffix;   // "_liquid"
};

/** 模型不做平衡时返回空。 */
virtual std::vector<EquilibriumComponent> get_equilibrium_components () const { return {}; }

/**
 * 直接求平衡态。pressure 用 Pa，temperature 用 K，bulk_composition 长度 = 组分数。
 * 返回 false 表示模型不支持。
 */
virtual bool
evaluate_equilibrium_state (const double pressure,
                            const double temperature,
                            const std::vector<double> &bulk_composition,
                            double &melt_fraction,
                            std::vector<double> &solid_composition,
                            std::vector<double> &liquid_composition) const { return false; }

/** 当前用于平衡的压力口径名（用于一致性自检）。 */
virtual std::string get_equilibrium_pressure_name () const { return ""; }
```

`MeltThermodynamicEquilibrium` 侧只要把已有的
`solve_eq_melt_fraction()` + `calculate_concentration_solid/liquid()` 包一层，
**不引入任何新物理**。

> 为什么不让插件直接调 `material_model->evaluate()`？
> 也行（`evaluate` 是 public，`create_additional_named_outputs` 也是 public，
> `MeltOutputs` 可以手工塞进 `additional_outputs`），但每次调用都要构造
> `MaterialModelInputs/Outputs`、算一遍密度/粘度/反应项，而反演要调用**几十次**；
> 而且 `timestep_number()==0` 的守卫、`in.strain_rate` 等约束要小心。
> 显式 API 更快、更清楚，还能被 `melting rate instantaneous` 那类后处理复用。
> 若你想零改动材料模型，这条是 fallback（写在 §10 的选项里）。

### 5.3 插件参数（完整清单）

```
subsection Initial composition model
  set List of model names = thermodynamic equilibrium

  subsection Thermodynamic equilibrium
    # ── 1. 我要指定什么 ────────────────────────────────────────
    set Specification = melt fraction        # bulk composition | melt fraction | phase compositions

    # ── 2. 压力和温度从哪来 ────────────────────────────────────
    set Pressure source = adiabatic          # adiabatic | depth profile | function
    # 温度没有选项：永远读 Initial temperature model（见 §10.4）

    # ── 3. 一维混合族（Specification = bulk composition 时不用）──
    set Mixing family = tie line
      #   tie line            A = c_s(P,T,c̄_bg), B = c_l(P,T,c̄_bg)      ← 你的 R2/R3/R4
      #   background to liquid A = c̄_bg,          B = c_l(P,T,c̄_bg)
      #   background to solid  A = c̄_bg,          B = c_s(P,T,c̄_bg)
      #   background to component <i>  A = c̄_bg,  B = e_i            ← "加 CO2" 方向
      #   component <i> to <j>         A = e_i,   B = e_j
      #   vector              A = c̄_bg,          B = normalize(c̄_bg + v)   v 由参数给
    set Family vector = 0, 0, 1              # 仅 'vector' 用

    # ── 4. 背景全岩成分 c̄_bg(x) ────────────────────────────────
    subsection Background bulk composition
      set Model = uniform                    # uniform | depth profile | function
      set Uniform composition = 0.7, 0.2993333333333333, 0.0006666666666667
      subsection Depth profile
        set Data directory = ./
        set Data file name = background_composition.txt   # 列名 = 组分名
      end
      subsection Function
        set Function expression = 0.7; 0.3-c; c
      end
    end

    # ── 5. 目标（Specification = melt fraction / phase compositions）──
    subsection Target melt fraction
      set Model = background                 # background | uniform | depth profile | function
      # 'background' = 让 φ* = Φ(P,T,c̄_bg)（即"只要背景平衡"，s = 背景混合分数）
      subsection Solitary wave
        set Add solitary wave        = true
        set Peak position            = 25000   # ASPECT 原生 y [m]（自下而上；= 225 km 深）
        set Peak porosity            = 3.0e-3  # 绝对峰值孔隙度 φ_peak
        set Truncation threshold     = 1.0e-4  # (φ_r − 1) < tol 之外不嵌入
        set Derive compaction length = true    # δ 由材料模型导出（见 §5.5.2）
        set Compaction length        = 1.0     # m；仅在上面为 false 时使用
      end
    end

    # ── 6. 反演与自检 ─────────────────────────────────────────
    set Solving tolerance            = 1e-14   # IC 独立容差，比时间步的 1e-10 紧（§10.5）
    set Consistency check tolerance  = 1e-8    # Specification = phase compositions 用
    set Report initial equilibrium   = true    # 通过 post_set_initial_state 打印统计
    set Dump generated state to file = false   # 写 ASCII 供与 python 老流程逐点对照
  end
end
```

> 失败行为固定为**报错**（你的选择），所以没有 `Unreachable target` 之类的开关：
> 目标不可达、族成员 `Σc̄≠1`、相成分不自洽，三者都在 `AssertThrow` 里带上下文抛出。

### 5.4 三种模式的算法

**(a) `Specification = bulk composition`（前向）**

```
c̄(x) = normalize(model(x))
(φ, c_s, c_l) = Φ(P, T, c̄)
```
没有 s，没有迭代。这一模式的价值是：**把 c̄ 的书写交给 ASPECT 的
`uniform/剖面/函数` 机制，而 φ/c_l 由内核算出**——即"指定全岩成分，其余自然算出"。

**(b) `Specification = melt fraction`（反演）**

```
1. (φ_bg, c_s,bg, c_l,bg) = Φ(P, T, c̄_bg)          # 背景
2. φ*(x) = TargetMeltFraction(x, φ_bg)              # 可选叠加孤立波
3. 若 s 的解析值可用（family = tie line 且 φ_bg > 0）：
       s = φ* ,  c̄ = (1−s)c_s,bg + s c_l,bg
   否则：在 [0,1] 上二分 s 使 Φ(P,T,c̄(s)) = φ*
4. (φ, c_s, c_l) = Φ(P, T, c̄)                        # 回代，得到真正写出的值
5. 自检 |Φ − φ*| ≤ Solving tolerance
```

**为什么 tie-line 族 s ≡ φ 是精确的（可以跳过二分）**：
`c_s,i = K_i c_l,i` 时

```
Σ_i c̄_i (K_i−1)/(f+(1−f)K_i)
  = Σ_i [(1−f)c_s,i + f c_l,i](K_i−1)/(f+(1−f)K_i)
  = Σ_i c_l,i (K_i−1) = Σ c_s − Σ c_l          （因为 (1−f)K_i + f = f+(1−f)K_i）
```

而平衡条件正是 `Σc_s = Σc_l`（`f>0` 时机器精度成立）⟹ `f = s` 精确是根。
这就是你 python 里 `c̄=(1−f)c_s,bg+f c_l,bg` 回代误差只有 `1.9e-19` 的原因，
也是 `ZERO_GAMMA_CALIBRATION.md` 里 Γ≡0 的代数基础。

> **前提**：`Σc_s,bg = Σc_l,bg` 只在 `φ_bg > 0` 时成立（`f=0` 时 `Σc_l ≠ 1`，
> 你自己在文档里写过：P=7.23 GPa、T=1330.67 °C、c̄_cmorb=3e-4 时 Σc_l = 0.664）。
> 所以实现里**解析路径要带 `φ_bg > 0` 守卫**，否则走二分。
> 更严重的是下面 §5.4.2：全固区里 `c_l,bg` 本身不是合法组成，整个族都不成立。

**单调性与括号**：`dΦ/ds > 0` 沿一条合理方向成立；`tie line` 族还有端点恒等式
`Φ(s=0)=0, Φ(s=1)=1`，所以 `[0,1]` 天然是括号。其余族要显式检查
`Φ(0) ≤ φ* ≤ Φ(1)`，不成立就按 `Unreachable target` 处理（默认 `error`，
并打印位置、(P,T)、φ*(x)、`Φ(0)`、`Φ(1)`）。

**(c) `Specification = phase compositions`（质量平衡 + 校验）**

```
c̄ = (1−φ*)c_s* + φ*c_l*          # 闭式，Σ 自动为 1（若输入自洽）
(φ', c_s', c_l') = Φ(P,T,c̄)       # 回代
自检 max(|φ'−φ*|, |c_s'−c_s*|, |c_l'−c_l*|) ≤ Consistency check tolerance
```
不一致时的选项（`Resolution`）：
* `error`（默认）——直接报错并打印残差；
* `keep melt fraction`——丢掉 c_s*/c_l*，把 c̄ 当"只保证 φ*"的输入（退化成模式 b）；
* `keep liquid`——保留 φ* 与 c_l*，用模式 b 解 c̄（因为 c_l 定了、tie-line 定了，c̄ 也就定了）。

后两个选项让用户在"手写的相成分不在 tie-line 上"时有个明确的降级路径。

#### 5.4.1 这套算法已经用数值原型对着你的权威 IC 文件验证过

`prototype_equilibrium_ic.py`（同目录）把上面的算法照着实现了一遍，
输入用 `initial_composition_for_carbon_case06_uniphase.txt` 的 T 与它自己的 φ*，
输出与文件里的 c̄ 逐点比较（14799 个内部节点）：

| 检查 | 结果 |
|---|---|
| 解析路径 `c̄ = (1−φ*)c_s,bg + φ*c_l,bg` vs 文件 c̄ | **max 4.5e-13**（逐分量） |
| 回代 `Φ(P,T,c̄) − φ*` | **max 2.2e-18**（机器精度） |
| 通用二分路径求出的 s vs 解析的 s | **max 2.2e-18** |
| 二分路径重建的 c̄ vs 文件 c̄ | **max 4.5e-13** |

波峰（225 km）：`φ* = 3.000e-03`，`s = 0.003000000000`，
`c̄ = (0.69811622, 0.29855198, 0.00333180)`，
背景 `φ_bg = 1.073e-04`、`c_s,bg = (0.70006989, 0.29936232, 5.6778e-04)`、
`c_l,bg = (0.04884452, 0.02924757, 0.92190791)`。
内部节点里 3554 个是全固（φ*=0，走 s=0 的退化点）。

即：**新算法能逐位复现你现在的 IC**，而且通用二分路径与 tie-line 解析路径同样精确——
所以其它混合族（加 CO2、指定端元…）可以放心用同一个二分器。

#### 5.4.2 全固区：tie-line 族在那里**不是"方向反了"，而是根本不成立**

（这是我上一版设计稿里写错的一句话，这里更正，并给出实测数字。
复现脚本：`prototype_mixing_families.py`。）

`tie line` 族的两个锚点是背景态的相成分 `A = c_s,bg`、`B = c_l,bg`。
由 `c_s = K·c_l`、`c_l = c̄/(f+(1−f)K)`：

* **`φ_bg > 0` 时**：`Σc_s,bg = Σc_l,bg = 1`（平衡条件的定义），
  所以**族里每个成员都 Σ=1**，是一个合法的、真实的 tie-line 段，
  而且 `Φ(s) ≡ s` 精确成立。
* **`φ_bg = 0` 时**：求解器走"全固"分支，`f=0` ⟹
  `c_s,bg = c̄_bg`（Σ=1 ✓），但 `c_l,bg = c̄_bg/K` 是一个**假想的液体**，
  **连 Σ=1 都不满足**。于是 `c̄(s) = (1−s)c_s,bg + s c_l,bg` 的
  `Σc̄(s) = 1 − s·(1 − Σc_l,bg) < 1`。
  **ASPECT 的平衡内核不做归一化**（`.cc:295-311` 只把负浓度夹到 0），
  所以这些成员不是合法的全岩成分，整条"族"没有意义。

实测（200 ppm CO₂ 背景，case6 的 P/T 剖面）：

| 深度 | φ_bg | `Σc_l,bg` | `Σc̄(s)`，s = 0 / 0.25 / 0.5 / 0.75 / 1 |
|---|---|---|---|
| 110 km | 0 | **0.1565** | 1.000 / 0.789 / 0.578 / 0.367 / **0.156** |
| 130 km | 0 | **0.6701** | 1.000 / 0.918 / 0.835 / 0.753 / **0.670** |
| 225 km | 1.073e-4 | 1.000000 | 1.000 / 1.000 / 1.000 / 1.000 / 1.000 |

（110 km 处 `K = (87.7, 29.4, 4.82e-3)`；cmorb 的 K≪1 所以 `c_l,bg(cmorb)=0.138`，
而 Σ 只有 0.156。）

**这正好解释了你 python 脚本里那条 mask**
`mask = (shape > 0) & (f_bg > 0)`：全固区必须排除，不是因为方向不对，
而是因为**那里没有 tie-line 可言**。

**那全固区想造熔融怎么办？** 换一个方向——富集。同一深度实测：

| family | Φ(s)，s = 0 / 0.25 / 0.5 / 0.75 / 1 | 能否到 φ*=1e-3 |
|---|---|---|
| `tie line` | 0 / 0.0407 / 0.118 / 0.285 / 1.0（成员 Σ≠1，无意义） | — |
| `bg → pure cmorb` | 0 / 0.251 / 0.507 / 0.763 / 1.0 | ✅ s≈0.004 |
| `bg → pure morb` | 0 / 0 / 0 / 0 / 0 | ❌ 永远不熔 |

（`bg → pure morb` 不熔是物理的：加难熔的 morb 是把固相线**抬高**。）

**所以实现规则定为**（你 2026-09-29 的澄清："给定的初始成分场在某个区域就在固相线下方，
那就不用考虑液体的事情了"）：

1. **"全固"的判据不是错误，而是一种正常状态**：某点 `c_l,bg` 的 Σ≠1，
   等价于该点 `φ_bg = 0`。这时**液相浓度本来就没有物理意义**，
   不需要报错、也不需要特殊处理——照材料模型的返回值写出去即可
   （它随后会被 `fill_prescribed_fields()` 重算，值本来就无关紧要）。
2. **`Mixing family = tie line` 自动只在 `φ_bg > 0` 的点上使用**。
   `φ_bg = 0` 的点一律保持背景：`c̄ = c̄_bg`、`φ = 0`。
   这条规则逐位复刻你 python 脚本里的 `mask = (shape>0) & (f_bg>0)`，
   不需要额外参数、也不需要用户操心——所以参数 `If background is all solid`
   已被删除。
3. **唯一真正会报错的情形**：在 `φ_bg = 0` 的点上却要求 `φ* > tol`
   （目标在所选族上不可达）。这时按 §10.3 报错，信息里给该点 (P,T)、
   `φ_bg`、`φ*`，并提示改用 `background to component <最不相容的组分>`。
   你自己的 case6 波形永远不会触发它——波窗（225 ± 10 km）内 `φ_bg > 0`，
   窗外 `φ* = φ_bg = 0`。

> 换句话说：**"目标熔融分数为 0" 是合法的、不需要反演**；
> 只有"要在固相线下方硬造出熔融"才是错误输入。

### 5.5 目标场 / 背景场的四个 model

都只需要 `double f(const Point<dim>&)`：

| model | 实现 | 用途 |
|---|---|---|
| `uniform` | 常数（标量或向量） | 背景、快速试验 |
| `depth profile` | `Utilities::AsciiDataProfile<dim>`（`structured_data.h:658`，支持列名、一维插值） | 你现在的 `.txt` 剖面 |
| `function` | `Functions::ParsedFunction<dim>` | 解析式 |
| `solitary wave` | 内置，见下 | 你的 case6 波形 |

`Solitary wave` 的内核逐字移植 `classes_and_functions.py:153-206`：

```
implicit_function(φ_r) = | sqrt(A+1/2) · [ 2·sqrt(A−φ_r) − (1/sqrt(A−1))·ln((sqrt(A−1)−sqrt(A−φ_r))/(sqrt(A−1)+sqrt(A−φ_r))) ] |
φ_r(x)                 = brentq(implicit_function(φ_r) − |x − x_peak|/δ, 1, A)
```

#### 5.5.1 无量纲约定：`φ_r ∈ [1, A]`（背景 = 1，峰 = A）

按孤立波理论的习惯，**无量纲形状就是 `φ_r` 本身**，尾部趋于 1、峰处等于相对振幅 `A`。
r4 版设计稿里那个 `shape = (φ_r − 1)/(A − 1) ∈ [0,1]` 是多余的归一化，已删除。
于是嵌入公式是**乘法**而不是加法：

```
φ*(x) = f_bg(x) · φ_r(x)                    ← 用当地背景孔隙度缩放当地无量纲振幅
```

#### 5.5.2 三个参数全部**导出**，用户只需给三个输入

（2026-09-30 你的方案，已采纳。）

| 用户指定 | 含义 |
|---|---|
| `Peak position` | 波峰位置（ASPECT 原生 `y`） |
| `Peak porosity` | 波峰的**绝对**孔隙度 `φ_peak` |
| `Truncation threshold` | 截断阈值 `tol`：`φ_r − 1 < tol` 之外不嵌入 |

导出过程：

```
1. φ_0 = Φ(P(x_peak), T(x_peak), c̄_bg)          # 峰处的背景平衡孔隙度
                                                 #   —— 背景成分已知，所以它是确定的
2. 用材料模型在该状态、取 φ = φ_0 求值，读回
       k = melt_out.permeabilities[0]
       ξ = melt_out.compaction_viscosities[0]
       η = out.viscosities[0]
       μ_f = melt_out.fluid_viscosities[0]
   δ = sqrt( k · (ξ + 4η/3) / μ_f )              # 压实长度，由模型自己的本构给出
3. A = φ_peak / φ_0                              # 相对振幅
4. 窗口半径 r = δ · implicit_function(1 + tol)   # **直接求值**，不需要搜索
   （implicit_function(1+tol) 对 tol>0 有限；tol→0 时 →∞）
5. 窗内 φ*(x) = f_bg(x)·φ_r(|x−x_peak|/δ)，窗外 φ*(x) = f_bg(x)
```

★ 注意第 2 步**不引入任何新参数**：`k(φ)`、`ξ(φ)`、`η`、`μ_f` 全部由材料模型在
"峰处状态 + `φ = φ_0`"下求值给出，所以压实长度自动跟随模型的本构
（包括 `Reference permeability`、`Reference bulk viscosity`、
`Viscosity activation energy` 的 E+PV 因子、`Background porosity` 的偏置等），
再也不会出现"`.prm` 里改了渗透率、波形却纹丝不动"的情况。
计算出的 δ、A、窗口半径、FWHM 会在 IC 阶段打印出来供核对。

**这个构造比旧版多了三条好性质**：

1. **全固区自动归零**：`φ*(x) = f_bg(x)·φ_r(x)`，而 `f_bg = 0` 处
   `φ* = 0` ⟹ **不再需要 `mask = (shape>0) & (f_bg>0)`**，也不需要
   "`φ_bg=0` 就报错"的规则（§5.4.2 那条自动满足）。
2. **窗外连续回到背景**：截断处的跳变是 `f_bg·tol`（本算例 `1.05e-4 × 1e-4 ≈ 1e-8`），
   可以忽略。
3. **A 有物理含义**：旧版 python 里 `background_porosity=1e-3` 与 `amplitude=20`
   相乘是 0.02，而实际峰值是 3e-3 —— 两个参数都是摆设
   （调用方随即 `/1e-3` 把它除回去，峰值另由 `PHI_PEAK` 设定）。
   新版里 `A` 就是"峰处放大倍数"，`φ_0·A = φ_peak` 精确成立。

**必须加的守卫**（报错，不做隐式夹逼）：

* `φ_0 ≤ 0`（峰位落在固相线以下）⟹ 报错，提示"把峰位放进熔化区，或提高背景不相容组分含量"；
* `A ≤ 1`（`φ_peak ≤ φ_0`，波不存在）⟹ 报错；
* `r` 超出模型域或超出熔化区 ⟹ **只警告不报错**（因为 `f_bg=0` 处 `φ*=0`，
  构造本身是安全的），但要打印 `r` 与熔化区范围供判断。

#### 5.5.3 ⚠️ 代价：波宽对背景成分**非常敏感**（`δ ∝ φ_0^{3/2}`）

因为 `k ∝ φ_0³` 而当前 `.prm` 里 `Exponential melt weakening factor = 0`
（`ξ = ξ_0·exp(−α_φ φ)` 退化成常数 `ξ_0`），所以 `δ ∝ φ_0^{3/2}`。
实测（case6 的峰位 225 km，模型常数 `ξ = η = 2.1424e20`、`μ_f = 5e-3`、`K_ref = 1e-8`）：

| 背景 CO₂ | `φ_0 = f_bg(peak)` | `δ` [m] | `A = φ_peak/φ_0` | 窗口半径 @tol=1e-4 [km] | FWHM [m] |
|---|---|---|---|---|---|
| 200 ppm（case6） | 1.0545e-04 | **34.2** | 28.45 | 2.42 | 2854 |
| 300 ppm | 4.672e-04 | 319 | 6.42 | 8.35 | 6805 |
| 500 ppm | 1.191e-03 | 1299 | 2.52 | 25.7 | 14329 |
| 1000 ppm | 3.000e-03 | 5195 | **1.000** | —（A=1，波不存在） | — |
| 2000 ppm | 6.617e-03 | 17019 | **0.453** | —（A<1，非法） | — |

即：背景 CO₂ 从 200 提到 500 ppm，波宽涨 **38 倍**（物理上自洽，但形态不再可控）；
1000 ppm 时 `φ_peak` 恰好等于背景孔隙度，`A = 1`，波消失。

**结论**：默认按你的方案**全部导出**；但保留一个逃生口
`set Derive compaction length = false` + `set Compaction length = <m>`，
供"要一条宽度受控的波"的对照实验使用。
（`A` 不需要逃生口——它由绝对峰值孔隙度唯一确定，这正是你要的语义。）

### 5.6 压力与温度来源，以及"必须一致"的检查

* **P 默认 = `adiabatic_conditions->pressure(p)`**。这不是随便挑的：
  `compute_initial_pressure_field()` 就是把这条剖面写进 `p_f` 的，
  而 `Pressure for thermodynamic equilibrium` 的默认也是它 ⟹ **t=0 时 IC 与模型看到的 P 是同一个数**。
* 提供 `depth profile` / `function` 供你做实验（比如"故意用错压力"看看差多少），
  但这时要打印一行警告。
* **一致性检查**：拿 `get_equilibrium_pressure_name()`，如果材料模型的口径不是
  `adiabatic pressure`，`initialize()` 里 `pcout` 警告（`solid pressure` 在 t=0 更是
  直接不可用，上一轮已实测：φ_peak 6.1e-4 vs 3.0e-3、速度差 2000 倍）。
* **T 直接问温度 IC**：`initial_temperature_manager->initial_temperature(p)`。
  这保证"IC 用的是你真正写进温度场的那条 T"，不需要在 `.prm` 里再抄一遍。
  温度暂不内置（你的决定）——所以 T 的来源和写法完全不变，
  `ascii data` / `function` / 将来的 `PracticalTemperatureProfile` 插件都可以。

#### 5.6.1 ⚠️ 循环依赖：绝热参考剖面建好之前，插件会被调用一次

case6/region5 的 `.prm` 写着
`set Composition reference profile = initial composition`，于是
`ComputeProfile::initialize()`（`compute_profile.cc:158`）会对每个参考点调用
`initial_composition_manager->initial_composition(p, c)` —— 也就是调用本插件。
这次调用发生在 `adiabatic_conditions->initialize()`（`core.cc:372`）**内部**，
早于 `set_initial_temperature_and_compositional_fields()`（`core.cc:1999`），
此刻 `adiabatic_conditions->is_initialized() == false`。

（`ComputeProfile::update()` 只在 `Use surface condition function = true` 时才重建，
所以本算例这条剖面**只建一次**，循环也只发生一次。）

**规则**：

* `is_initialized() == false` ⟹ 返回**背景成分 c̄_bg**（它不需要 P），派生场返回 0；
* `is_initialized() == true` ⟹ 走完整的反演流程。

对本算例无影响：密度是 `3400·(1−αΔT)·(1+βΔP)`，而
`Thermal expansion coefficient = 0`、`Solid compressibility = 0`，
参考密度恒为 3400，与参考成分无关。
（若将来密度依赖成分，用背景成分当参考剖面本来也是合理选择。）

**想彻底避开**：把那一行换成 `set Composition reference profile = function`
+ 一个常数 `Function expression`，参考剖面就与 IC 插件完全解耦。
* ⚠️ 两个 IC manager 在 step 0 解完之后被 `reset()`，
  所以插件要在 `post_set_initial_state`（主循环之前）里就把 T 用完；
  更晚需要 T 的话得自己缓存。

### 5.7 `.prm` 示例：把 case6_05_uniphase 完全搬进内核

```prm
subsection Compositional fields
  set Number of fields = 11
  set Names of fields = porosity, dunite, morb, cmorb, dunite_liquid, morb_liquid, cmorb_liquid, melting_rate, debug_1, debug_2, debug_3
  set Compositional field methods = prescribed field, field, field, field, prescribed field, prescribed field, prescribed field, prescribed field, prescribed field, prescribed field, prescribed field
end

# 温度仍然完全自由 —— 这里继续用它原来的外部文件（也可以换成 function/adiabatic）
subsection Initial temperature model
  set List of model names = ascii data
  subsection Ascii data model
    set Data directory = ./
    set Data file name = initial_temperature_for_carbon_case06_uniphase.txt
  end
end

# 全岩成分 + 熔融结构：不再需要任何 composition 的 txt
subsection Initial composition model
  set List of model names = thermodynamic equilibrium
  subsection Thermodynamic equilibrium
    set Specification  = melt fraction
    set Mixing family  = tie line
    set Pressure source = adiabatic

    subsection Background bulk composition
      set Model = uniform
      set Uniform composition = 0.7, 0.2993333333333333, 0.0006666666666667   # 200 ppm CO2
    end

    subsection Target melt fraction
      set Model = background
      subsection Solitary wave
        set Add solitary wave        = true
        set Peak position            = 25000   # ASPECT 原生 y [m]（= 225 km 深）
        set Peak porosity            = 3.0e-3
        set Truncation threshold     = 1.0e-4
        set Derive compaction length = true
      end
    end
  end
end
```

配合 §5.10 的边界条件修复（**已完成**），这份 `.prm` 不再需要
`initial_composition_for_carbon_case06_uniphase.txt`（4.5 MB）；
`initial_temperature_...txt` 继续用（温度暂不内置，见 §10.4）。
更完整的"实现后怎么用"说明见 `INTERNAL_IC_USAGE.md`。

### 5.8 自检与报告

* **每个点**在回代后算 `|Φ(c̄) − φ*|`，累计 `max`、Σ、计数；反演失败/超差按
  `Unreachable target` 处理。
* 在 `signals.post_set_initial_state`（`simulator_signals.h:92`）里打印一次汇总：

```
[TEQ initial condition] specification = melt fraction, family = tie line, P = adiabatic
  nodes                       : 45000   (interior 44402)
  points with phi* = 0        : 40102
  points solved analytically  : 4898    (tie-line family, phi_bg > 0)
  points solved by bisection  : 0
  points unreachable          : 0
  max |Phi(c_bar) - phi*|     : 3.1e-15
  max |Sigma c_bar - 1|       : 4.4e-16
  max |phi_stored - phi*|     : 8.9e-12      <-- 与 IC 的最终对照
  cMORB at wave peak          : 3.330132787e-3  (999.0 ppm CO2)
```

  最后一行 `phi_stored` 需要读一下 `solution`（signal 里有 `SimulatorAccess`），
  这是"IC 到底落在平衡态上没有"的**唯一权威判据**，也是替代
  `TEQ_PRESSURE_VERIFICATION.md` 那套手工对照的自动化版本。
* `Dump generated state to file = true` 时写一个 ASCII（同一套 `# POINTS:` 格式），
  可以直接和你现在的 python 产物逐点 `diff`，作为迁移期的回归测试。

### 5.9 需要改动的文件

| 文件 | 改动 | 估计行数 |
|---|---|---|
| `include/aspect/material_model/interface.h` | `EquilibriumComponent` + 3 个默认虚函数 | +30 |
| `include/aspect/material_model/melt_thermodynamic_equilibrium.h` | 3 个 override 声明 | +8 |
| `source/material_model/melt_thermodynamic_equilibrium.cc` | 3 个 override 实现（包一层已有函数） | +35 |
| `include/aspect/initial_composition/thermodynamic_equilibrium.h` | **新** | ~160 |
| `source/initial_composition/thermodynamic_equilibrium.cc` | **新**（含 uniform/profile/function/solitary wave 四个小 model + 二分） | ~800 |
| `source/simulator/core.cc` | §5.10 的 BC 跳过 | ~20 |
| `source/simulator/parameters.cc` | （可选）显式排除列表 | ~10 |
| `doc/` / `tests/` | 参数文档、一个小回归算例 | — |

新文件在 `source/initial_composition/` 下会被 `FILE(GLOB_RECURSE)` 自动纳入编译。

设计验证脚本（都已跑通）：
* `prototype_equilibrium_ic.py` —— 反演算法逐位复现现有 IC（§5.4.1）
* `prototype_mixing_families.py` —— 各族锚点在全固/熔融区的实测行为（§5.4.2）

### 5.10 独立小改动（**P0，必须先做**）：prescribed 场不该有 Dirichlet 边界条件

> 你 2026-09-29 的强调："务必让 prescribed 场不受边条件影响，否则边界点
> 可能一直保持初始值（比如 0）而不遵从平衡条件。" 这一条因此是硬性要求，
> 不是可选项，而且要作为**验收判据**：修好之后 `porosity`（含 `*_liquid`、
> `melting_rate`、`debug_*`）在**包括上下边界节点在内的全网格**上，
> 都必须等于材料模型的平衡值。

**位置**：`Simulator<dim>::compute_current_constraints()`，
`core.cc:717-744` 的

```cpp
for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
  for (const auto p : boundary_composition_manager.get_fixed_composition_boundary_indicators())
    VectorTools::interpolate_boundary_values (...);
```

**改法**：

```cpp
// 纯 prescribed 场不是任何 PDE 的解，它的值完全由材料模型给出；
// 施加 Dirichlet 只会（a）在边界制造非平衡值，（b）在
// interpolate_material_outputs_into_advection_fields() 的
// current_constraints.distribute() 里把重写结果再盖掉。
if (parameters.compositional_field_methods[c]
      == Parameters<dim>::AdvectionFieldMethod::prescribed_field)
  continue;
```
（温度场如果 method 是 `prescribed field` 同理。）

**为什么正确**：`prescribed_field` 的定义就是"值由材料模型给"；
`prescribed_field_with_diffusion` 仍然解一个扩散方程，保留 BC 有意义。

**影响**：case6/region5 的 8 个 prescribed 场在上下边界的节点值会从
"IC 文件的值"变成"平衡值"。上一轮实测这个差在 250 km 处是 ~1.5e-6（相对 1.5%），
更重要的是它让 `porosity` 场在**全网格**上（含边界）都严格等于平衡解，
这是 `melting rate instantaneous` 那类诊断的前提。

**先例**：分支 `melt_rate_pure_convection_diag` 的 `6f6830314` 已经用
"字段名以 `_pure_convection` 结尾就跳过"的硬编码做过同样的事；
这里把它换成按 method 判断的通用规则，并保留一个可选的显式排除列表。

---

## 6. 精度预算

| 误差源 | 现在（外部 ASCII 流程） | 内建后 | 依据 |
|---|---|---|---|
| 热力学参数舍入（python 精确值 vs `.prm`） | **2.16e-6**（step-0 `φ_ASPECT−φ_IC` max） | **0**（同一个内核、同一份参数，没有第二份数据） | 本报告 §0 实测 |
| ASCII 文件插值（100 m 文件 vs 0.625 m 节点） | 折点/曲率 **11.21×** | **0**（没有文件） | `make_solitary_wave_ic_fine.py` |
| float32 存储/往返 | 后处理链上 ε ~2e-10（相对 1.07e-6 的 0.11 倍） | **0**（IC 全程 float64，不经过 vtu） | `MELTING_RATE_NUMERICS.md` |
| 节点 vs 求积点 | 诊断口径问题（ε/Δt） | **不变**，但这不是 IC 的误差：prescribed 场本来就在节点上求值，动力学在求积点看插值后的 φ，是离散格式本身的性质 | `helper_functions.cc:2021` |
| `Surface pressure` 舍入（3.06e9 vs 3.06072e9） | 7.7e-8 in φ | **0**（IC 用 `adiabatic_conditions->pressure()`，与模型同一个数） | 本报告 §0 |
| 二分容差 | python `brentq` xtol 尾部 5e-13 → 经 1/Δt 放大过 3.45% | 用与材料模型**同一个** `Equilibrium solving tolerance`（默认 1e-10，可调 1e-14） | `MELTING_RATE_NUMERICS §12` |
| 目标值本身（`φ_peak`、背景 ppm 等） | 人手写 | 不变 | — |

**期望的最终 step-0 失配：~1e-11**（只剩浮点求和/投影），
与你文档里"要发表绝对熔体分数时压到 1e-7"的目标差两个数量级。

---

## 7. 边界情况与失败模式（实现时必须处理）

| 情况 | 处理 |
|---|---|
| `φ_bg = 0`（全固区） | **不是错误状态**：该点没有熔体，液相浓度无物理意义，照写即可。`tie line` 族在该点自动退化为"保持背景"（`c̄=c̄_bg`, `φ=0`），见 §5.4.2。 |
| 族成员不合法（`Σc̄(s) ≠ 1`） | 只在 `φ_bg = 0` 时出现，此时该族**不被使用**（自动走"保持背景"）。在 `φ_bg > 0` 的点上族成员必然 Σ=1；实现里仍加一条 `|Σc̄−1| ≤ 1e-12` 的断言作为防御。**不做隐式归一化**（ASPECT 内核也不归一化）。 |
| `φ_bg = 0` 的点上却要求 `φ* > tol` | **报错**，信息给出该点 (P,T)、`φ_bg`、`φ*`，并提示改用 `background to component <最不相容的组分>`（§5.4.2 实测：110 km 处 `bg → pure cmorb` 在 s≈0.004 就到 φ*=1e-3，而 `bg → pure morb` 永远不熔）。 |
| 目标 φ* 在所选族上不可达（`φ* ∉ [Φ(0), Φ(1)]`） | **报错**（你的选择），信息里给位置、(P,T)、φ*、`Φ(0)`、`Φ(1)`、所用 family，以及二分是否找到括号。**不做隐式夹逼**。 |
| 用户在手写相成分时分量数 ≠ 模型分量数 | 启动时 `AssertThrow`。 |
| 组分场缺失（`*_solid` 可选、`*_liquid` 必需） | 用 `get_equilibrium_components()` 返回的后缀查 `introspection`；`*_liquid` 缺 → 报错；`*_solid` 缺 → 跳过写入。 |
| 用户同时列了别的 IC 插件 | operator 是 `add`，会串味。文档里明确"这个插件应当单独出现"；`initialize()` 里检查 `List of model names` 长度并警告。 |
| 每个场调用一次 ⟹ 每节点算 11 次 | 加一个**单点缓存**（`mutable Point<dim> last_p; mutable State last_state;`）。`VectorTools::interpolate` 对同一个场遍历整网格、对 11 个场各遍历一遍，遍历顺序相同 ⟹ 缓存命中率 ≈ 10/11。 |
| 与剖面文件的坐标约定 | **不做任何深度换算**：ASPECT 的原生坐标就是 `y` 自下而上增大，所有 `depth profile` / `function` model 一律直接吃 `Point<dim>`（和 `AsciiDataProfile`、`ParsedFunction` 的既有约定一致）。"深度"只在**绘图**时换算，不进代码。（你 2026-09-29 的决定。） |
| MPI | 纯逐点计算，无通信；统计量用 `Utilities::MPI::max` 归约。 |
| `AssertThrow` 在 IC 里 | ASPECT 的 IC 插值有 try/catch（`initial_conditions.cc:120-181`），抛异常会 MPI_Abort。所以**先收集再统一报**，错误信息里带坐标。 |

---

## 8. 兼容性与迁移

* **完全增量**：不动任何现有插件与默认行为（除了 §5.10 那个边界条件修复，
  它会让 prescribed 场的边界节点值变成平衡值——这是修正，但会改变结果，
  建议单独一个 commit + 一个 A/B 算例）。
* **Phase 0 可以立刻用**：只做 §5.10，然后你的**前向模式（R1/R8/R10）
  今天就已经是内建的了**——把 `porosity` / `*_liquid` 的 IC 列随便写 0 就行，
  因为它们在第一次 Stokes 之前就被重算了。你只需要保留 c̄ 的 ASCII 文件。
  这已经能消掉"参数舍入"和"边界锁死平衡态"两个问题。
* **老 `.prm` 不动**：`Initial composition model = ascii data` 继续工作。
* **两条路径并行**（你的决定）：外部 ASCII 与内建插件都是一等公民。
  代价是必须维护"两条路径给同一套物理"，做法见 **§11**：
  ① `.prm` 写未舍入的精确值；② 修掉 python 侧两份同名 K 列表；
  ③ 用 `Dump generated state to file` 做逐点回归。
* **回归测试**：`Dump generated state to file = true` 的输出与你现在的
  `initial_composition_for_carbon_case06_uniphase.txt` 做逐点比较，
  预期差 ~1e-11（而不是现在的 2.16e-6）。这正好是一个可自动化的验收判据。

---

## 9. 分阶段实施

| 阶段 | 内容 | 产出 | 风险 |
|---|---|---|---|
| **P0** ⭐ | §5.10 边界条件跳过（按 `prescribed_field` method 自动跳过；温度场同理） | **✅ 已完成 2026-09-29**，见下面 | 低；会改变结果（边界节点不再被 IC 钉住） |
| **P1** | 材料模型 3 个 API（`get_equilibrium_components` / `evaluate_equilibrium_state`） | **✅ 已完成**：`interface.{h,cc}` + TEQ 覆写 | 低，纯包装 |
| **P2** | 插件骨架 + `Specification = bulk composition`（前向） | **✅ 已完成** | 低 |
| **P3** | `uniform` / `depth profile` 背景 model、`Specification = melt fraction`、`Mixing family = tie line`（默认）+ `background to liquid/solid/component`、`solitary wave`（δ/A/窗口全部导出）、自检报告 | **✅ 已完成**（偏差见 §9.2） | 中 |
| **P4** | `Specification = phase compositions`、背景的 `function` model、`dump to file`、压力源选项 | 未做 | 中 |
| **P5** | （暂缓，见 §10.4）注册为 `initial temperature`，加 `Free variable = temperature`（覆盖 R6） | 逆温度 | 中；两个 manager 的 `declare_parameters` 会撞名，需要 idempotent 处理 |
| **P6** | （可选，建议延后）R7 等焓 2×2 | — | 高 |

P0 建议单独一个 commit：它改的是 `source/simulator/core.cc` 的通用逻辑，
与插件无关，做完立刻可以验证（跑一次 case6_05，检查边界节点的 `porosity`
是否等于平衡值，而不是 IC 值）。

### 9.1 P0 已完成并验收（2026-09-29）

**改动**：`source/simulator/core.cc` 的 `compute_current_constraints()`——
组分场循环里跳过 `compositional_field_methods[c] == prescribed_field` 的场，
温度场同理跳过 `temperature_method == prescribed_field`。
`prescribed_field_with_diffusion` 不跳过（那种方法仍解扩散方程）。
二进制已重建：`build_TEQ_bulk_release/aspect`（17 s 增量编译）。

**验收算例**：`region5_carbon/case6_practical_0724/verify_prescribed_bc/`
（deck 与 `verify_teq_pressure/adiabatic.prm` **只差输出目录**；
基线 = 打补丁前那个 deck 的输出）。

| 检查 | 结果 |
|---|---|
| step-0 内部（`y ∈ (1 km, 149 km)`）所有场 | **逐位不变**（`max|Δ| = 0`） |
| step-0 三个 PDE 场 `dunite/morb/cmorb`（含边界） | **逐位不变** |
| step-0 边界节点 `porosity` vs 内核预测 `Φ(P_ad,T,c̄)` | pre-P0 **2.133e-06** → post-P0 **2.106e-12** |
| 最直观的例子：`debug_1`（K₀）在顶边界 | 补丁前一直是 IC 里的 `0`，补丁后 = `200.7` |
| 第 3 步传播范围 | `T`、`dunite`、`morb` 逐位不变；`cmorb` 内部 2.3e-10；`p_f` 相对 9e-7 |

细节与完整表格见该目录的 `README.md` 与 `summary.json`。
**注意**：`TEQ_PRESSURE_VERIFICATION.md` 里"边界节点 φ 等于 IC 值"的观测
是补丁前的历史现象；P0 之后边界节点改为服从平衡条件，
其中间/内部的数字不受影响（内部逐位不变）。

---

### 9.2 P1–P3 实现记录（2026-09-30）

**改动文件**（都在工作区未提交）：

| 文件 | 内容 |
|---|---|
| `include/aspect/material_model/interface.h` | `EquilibriumComponent` + 2 个默认虚函数 |
| `source/material_model/interface.cc` | 两个默认实现（返回空 / false） |
| `include/aspect/material_model/melt_thermodynamic_equilibrium.h` | 2 个 override 声明 |
| `source/material_model/melt_thermodynamic_equilibrium.cc` | 2 个 override 实现（薄包装 `solve_eq_melt_fraction` + `calculate_concentration_*`） |
| `include/aspect/initial_composition/thermodynamic_equilibrium.h` | **新**（~190 行） |
| `source/initial_composition/thermodynamic_equilibrium.cc` | **新**（~700 行） |

**验收**：`region5_carbon/case6_practical_0724/verify_internal_ic/`

| 检查 | 结果 |
|---|---|
| 生成场 vs 独立重算的平衡映射 | **5.3e-10**（float32 输出限制；插件内存残差 **7.2e-13**） |
| 生成场 vs 解析孤立波目标 | 1.2e-08 |
| **窗外** 与外部 ASCII 初始条件逐点比较 | **逐位相同（0.0）** |
| 波峰 | φ = 3.000000e-03，A = 28.45，δ（导出）= 34.2384 m，窗口 ±2419.42 m |
| 运行时间 | 3 步（45 yr）含 IC 生成共 16 s（8 MPI） |

**与设计稿的偏差**（都是缩减，不是行为改变）：

* `Peak position` 用 **ASPECT 原生竖直坐标**（2D 盒子里是 y，自下而上），
  不是 `geometry_model->depth()`（后者在盒子里是"从盒顶算起"，225 km 会写成 125000，容易混）。
* 背景成分只实现了 `uniform | depth profile`；`function` 推迟到 P4
  （deal.II 的 `ParsedFunction::declare_parameters` 需要在声明期就知道分量数，
  而分量数只能从材料模型拿到）。
* `Pressure source` 只有 `adiabatic`；`Specification = phase compositions`、
  `Dump generated state to file`、`Family vector` 一律推迟到 P4。
* 目标熔融分数的 `function` model 已实现（标量，声明期分量数=1，无问题）。

**两个必须记住的实现细节**（都已在代码里处理并有注释）：

1. **绝热参考剖面建好之前插件会被调用一次**：`Composition reference profile =
   initial composition` 使 `ComputeProfile::initialize()`（`core.cc:372`）提前调用
   初始成分插件，此时 `adiabatic_conditions->is_initialized() == false`，
   也不能安全访问温度 IC。插件直接返回**背景成分**（不需要 P）。
2. **initial temperature manager 在 step 0 之后被释放**（`core.cc:2042` 的
   `reset()`），但 `Boundary composition model = initial composition` 让
   `BoundaryComposition` 插件继续持有并调用初始成分 manager，于是 **step 1 的
   `compute_current_constraints()` 会再次调用本插件** —— 若此时再用
   `get_initial_temperature_manager()` 就会 SIGSEGV。按 ASPECT 文档的标准做法，
   在 `initialize()` 里把 `get_initial_temperature_manager_pointer()` 存进
   自己的 `shared_ptr`。

**性能备注**：单点缓存只在同一个场的遍历内有效（`VectorTools::interpolate`
对每个场各遍历一遍全网格），所以实际做了 `场数 × 节点数 ≈ 1.49e6` 次状态求值。
tie-line 族走解析路径，每次状态求值 = 2 次平衡求解；总耗时仍可接受（16 s 含 3 步）。

---

### 9.3 P4 的具体内容（待做）

P0–P3 已经把主路径跑通。P4 是把它从"够用"补到"完整"，五项，按价值排序：

**P4.1 `Specification = phase compositions`（模式 C，你最初三条诉求的第三条）—— ✅ 已完成 2026-09-30**

* 参数：`Phase composition models`（c_s* 与 c_l* 各一个 `uniform` / `depth profile` model，
  每个化学组分一列/一个值）+ 复用已有的 `Target melt fraction model` 当 φ*；
  `Consistency check tolerance`（默认 1e-8）、
  `If inconsistent = error | keep melt fraction | keep liquid`。
* 算法：`c̄ = (1−φ*)c_s* + φ*c_l*` → 断言 `|Σc̄ − 1| ≤ 1e-10` →
  回代 `(φ', c_s', c_l') = Φ(P,T,c̄)` → 与输入比。
* **必须在 `φ = 0` 处屏蔽 c_l 的比对**：全固区 `Σc_l ≠ 1` 是正常的，
  模型返回的 `c_l' = c̄/K` 与用户给的 c_l* 本来就不该相等
  （对应你说的"固相线下方就不用考虑液体的事情了"）。
* 降级路径的意义：手写的 `(c_s,c_l,φ)` 一般不在模型的 tie-line 上
  （尤其 .prm 的 K 系数被舍入时）。`keep melt fraction` 退化成模式 B；
  `keep liquid` 保留 φ* 与 c_l*（tie-line 由此确定）再解 c̄。
* 代价：~250 行 + 2 个 field model + 一组参数。

**P4.2 `Dump generated state to file = true`（迁移期回归对照）—— ✅ 已完成 2026-09-30**

* 在 `post_set_initial_state` 里（和 report 同一处）把生成的状态写成一份 ASCII：
  `# POINTS: 2 n_y`，列顺序 = `.prm` 的 `Names of fields`，
  y 用 ASPECT 原生坐标，x 列沿用外部生成器的 0/1 约定（本算例沿 x 均匀）。
* 价值：可以 (a) 直接与 `initial_composition_for_carbon_case06_uniphase.txt`
  逐点 diff，(b) 反过来用 `ascii data` 读回去做端到端回归。
* 代价：~80 行（IC 只有 3 万节点，收集到 rank 0 再排序写出即可）。

**P4.3 背景成分的 `function` model（便利性）—— 未做**

* 障碍：deal.II 的 `ParsedFunction::declare_parameters(prm, n_components)` 要求在
  **声明期**就知道分量数，而分量数只能从材料模型拿（`get_equilibrium_components()`），
  且 `declare_parameters` 是静态的、早于 `parse_parameters`。
* 解法（推荐）：自己把 `Function expression` 声明成 `Patterns::Anything`，
  parse 期按 `;` 拆成 n_components 个表达式，逐个构造 `Functions::ParsedFunction<dim>(1)`。
  （备选：声明期固定一个大分量数、parse 期 assert；或让用户再写一遍分量数——都与材料模型重复。）
* 代价：~60 行。

**P4.4 `Pressure source = depth profile | function`（逃生口）—— 未做**

* 只用于对照实验（"故意用错压力看差多少"）。默认仍是 `adiabatic`；
  非 adiabatic 时打印一行警告，因为那时 IC 与 `compute_initial_pressure_field()`
  写进 `p_f` 的剖面不再一致。
* 代价：~50 行。

**P4.5 其余 mixing family —— 未做**

* `vector`（`A = c̄_bg`, `B = normalize(c̄_bg + v)`，v 由参数给）与
  `component i to j`（`A = e_i`, `B = e_j`）。
* 代价：~40 行；价值低（`background to component` 已覆盖"富集造熔融"的物理需求）。

**明确不做 / 单独讨论**：

* `Peak position is depth` 开关：现在用原生竖直坐标（225 km 写 `25000`），
  加个 flag 只是把它写成 `125000`（`geometry_model->depth()` 在盒子里是从盒顶算起），
  反而更容易混 —— 不加。
* **P5（暂缓）**：把同一个类也注册为 `initial temperature`，实现"自由 DOF = 温度"
  （覆盖 R6 逆温度）。难点是两个 manager 的 `declare_parameters` 会撞名，
  需要 idempotent 处理。
* **P6（暂缓）**：R7 等焓 2×2。

---

### 9.4 P4.1 / P4.2 实现与验收记录（2026-09-30）

**新增/修改的代码**：

* `material_model/interface.{h,cc}`：`evaluate_equilibrium_state()` 增加一个
  `const double tolerance = -1.0` 参数（负值 = 用材料模型自己的
  `Equilibrium solving tolerance`）。
* `melt_thermodynamic_equilibrium.{h,cc}`：`solve_eq_melt_fraction()` 增加
  `const double tolerance = -1.0`（带默认值，所有既有调用点不变），
  并在二分处用 `(tolerance > 0 ? tolerance : equilibrium_tolerance)`；
  `evaluate_equilibrium_state()` 把它透传下去。
  **这修掉了一个真实缺陷**：插件此前解析了 `Solving tolerance` 却从未使用。
* `initial_composition/thermodynamic_equilibrium.{h,cc}`：
  * 模式 C：`Specification = phase compositions`，配 `Consistency check tolerance`
    与 `If inconsistent = error | keep melt fraction | keep liquid`；
    固/液相成分各有 `uniform | depth profile` 两个 model。
    **液相成分只在生成熔融分数 > 0 处参与比对**（全固区液相无定义）。
  * `Dump generated state to file` + `Output file name`：在
    `post_set_initial_state` 里写一份 `# POINTS: nx ny` 的 ASCII，
    列顺序 = `Names of fields`，可直接与外部产物 diff、也能被 `ascii data` 读回去。
  * 单点缓存从"一条"改成 `std::map`（同时服务 dump）。

**验收**（`region5_carbon/case6_practical_0724/verify_internal_ic/`）：

| 检查 | 结果 |
|---|---|
| 模式 C（把背景的精确相成分以深度剖面喂回去）vs 模式 B | **逐位相同** |
| 模式 C 自检 `max|相成分 − 输入|`（`Solving tolerance = 1e-14`） | **1.6e-14** |
| 模式 C 的 `max|Φ(c̄) − φ*|` | **6.7e-17** |
| 不一致输入 + `If inconsistent = error` | **exit 1**，报错给出残差 `0.00425`、位置、P、T |
| 不一致输入 + `If inconsistent = keep melt fraction` | 成功，输出与模式 B 逐位相同，`max|Φ(c̄)−φ*| = 6.7e-17` |
| dump 文件 | `# POINTS: 3 30001`，`max|φ_dump − φ_vtu| = 1.9e-10`（float32） |

**过程中发现并修掉的三个缺陷**（都由新功能暴露出来）：

1. **`Solving tolerance` 是死参数** —— 见上，已接通。
2. **缓存在绝热剖面构建阶段被污染**：`compute_state()` 在剖面未就绪时返回
   "背景成分 + φ=0" 的临时态，而这条路径也进了缓存；`ComputeProfile` 的 2000 个
   参考点里有 **(x=5, y=0)** 和 **(x=5, y=150000)** 恰好与网格节点重合，
   于是这两个节点的 prescribed 场被永久钉在 φ=0。
   *是 dump 功能把它暴露出来的*（dump 里多出 1998 行、且 (x=5,y=0) 的 φ 为 0）。
   修法：剖面未就绪时**不写缓存**。
3. **残差记录时机错了**：`max|Φ(c̄)−φ*|` 在模式 C 的降级回退**之前**记录，
   于是回退后的结果配着回退前的大残差（5.6e-5）。已移到一致性检查之后。

**尚未做**：P4.3（背景 `function`）、P4.4（压力源选项）、P4.5（`vector` /
`component i to j` 族）；P5（逆温度）、P6（等焓）。

---

## 10. 已拍板的设计决定（2026-09-29）

1. **外部 ASCII 与内建插件并行共存**，两者都是一等公民。
   → 因此必须保证两条路径产出**同一套物理**（见 §11）。插件的
   `Dump generated state to file = true` 就是为两条路径的逐点对照准备的。
2. **P 默认 = `adiabatic_conditions->pressure()`**，即从 `Surface pressure` 出发、
   用材料模型的参考密度（本算例恒为 3400 kg/m³）向下积分的静岩剖面。
   这正是 `compute_initial_pressure_field()` 在 t=0 写进 `p_f` 的那条剖面，
   所以 IC 与模型在 t=0 看到的是**同一个压力**。
   仍保留 `Pressure source = depth profile | function` 供对照实验用（会打印警告）。
3. **反演失败一律报错**（不可达、族不合法、相成分不自洽都算），
   不做隐式夹逼、不自动换方向。错误信息必须带：位置、(P,T)、φ*、所用 family、
   `Φ(0)`、`Φ(1)`、以及一条"该怎么办"的提示。
4. **温度暂不内置**。`Initial temperature model` 继续用现有机制
   （`ascii data` / `function` / 将来你自己的 `PracticalTemperatureProfile` 插件）；
   插件通过 `initial_temperature_manager->initial_temperature(p)` 读它，
   保证与真正写进温度场的那条 T 一致。
5. **平衡容差收紧**：IC 走一个独立的、默认更紧的容差
   （新参数 `Solving tolerance`，默认 `1e-14`），
   但**不改变**材料模型的 `Equilibrium solving tolerance` 默认值（仍 1e-10），
   以免影响正常时间步的性能与既有结果。
   > 实测依据：1e-10→1e-14 让 `max|∇c_l|`(cMORB) 降 31×（1.81e-11→5.91e-13），
   > 而 Γ 只变 1.12×；另外 1e-16 与 1e-14 逐位相同，所以 1e-14 已经到平台。
6. **全固区的 tie-line 问题**：见 §5.4.2——不是"方向反了"，而是
   该族在全固区不构成合法的一维族（`c_l,bg` 自己 Σ≠1）。
   处理规则见下面第 7 条。
7. **全固区不需要反演，也不需要报错**（你 2026-09-29 的澄清）：
   "给定的初始成分场在某个区域就在固相线下方，那就不用考虑液体的事情了"。
   → `tie line` 族自动只在 `φ_bg > 0` 的点上使用；`φ_bg = 0` 的点保持背景
   （`c̄ = c̄_bg`, `φ = 0`），液相浓度照材料模型的返回值写（无物理意义、随后即被重算）。
   这条规则**逐位复刻** python 的 `mask = (shape>0) & (f_bg>0)`，
   因此参数 `If background is all solid` 已删除。
   **只有"在固相线下方硬要求 φ*>0"才报错**（§10.3）。
8. **坐标系：不做任何深度换算**。ASPECT 原生 `y` 自下而上增大，
   插件的 `depth profile` / `function` model 直接使用 `Point<dim>`，
   与 `AsciiDataProfile`/`ParsedFunction` 的既有约定一致；
   深度只在**绘图**时换算。
9. **`Mixing family` 默认 = `tie line`**（对应你现在的两个 case6 脚本）。
10. **暂时不改动任何已存在的算例**（`.prm`、IC 文件一律不动）。
    §11.1 的精确值只用于**新建**的 deck；老 deck 保持原样也能从插件受益，
    因为插件的 IC 与它所配的那份 `.prm` 天生自洽（`617` 也好、精确值也好）。
11. **prescribed 场的边界条件必须去掉**（你的强调）：
    "否则边界点可能一直保持初始值（比如 0）而不遵从平衡条件"。
    §5.10 因此从"附带小改动"升级为 **P0、必须先做**，且是硬性验收项：
    修好之后 `porosity`（以及 `*_liquid`/`melting_rate`/`debug_*`）
    在**含边界在内的全网格**上都等于平衡值。

### 10.1 设计阶段无遗留问题

原来的三个小问题都已由你回答（坐标系 = 原生 `y`；family 默认 = `tie line`；
不动已存在的算例）。可以进入实施。

---

## 11. 精确值与"数据自洽"（对照你"可以采用精确值"的决定）

先说结论：**内建 IC 之后，"参数舍入不一致"这个问题会从根上消失**——
因为 `.prm` 成为热力学参数的**唯一**来源，IC 与时间步用同一个内核、同一份参数，
无论你写 `617` 还是 `617.27680612617`，IC 都严格落在 ASPECT 自己的平衡态上。
精确值只在两种情况下才重要：

* 你想让**内建路径**与**外部 ASCII 路径**给出同一个结果（你要求两条路并存）；
* 你想让结果与 `calibrating_June2026` 那条标定曲线严格对应（发表时的可追溯性）。

### 11.1 精确值：用于**新建** deck，不动已存在的算例

按你"暂时不要修改已经存在的 case"的要求，下面这些只写进**新** deck；
`melt_TEQ_carbon_case6_05/10/20/40.prm` 与现有 IC 文件一律保持原样。
（老 deck 用 `617` 也**不影响**插件的正确性：插件的 IC 和它所配的那份 `.prm`
用的是同一份参数，天生自洽。精确值只影响"与标定曲线的可追溯性"
和"与外部 ASCII 路径的逐点可比性"。）

```prm
set Surface pressure = 3.06072e9

set Melting point for each component at surface      = 1780, 1000, 617.27680612617
set Melting curve coefficient A for each component   = 45.0e-9, 112.0e-9, 3.799177670844276e-8
set Melting curve coefficient B for each component   = -2.00e-18, -3.37e-18, -4.001653781719964e-18
set Melting curve pressure threshold for each component = 6e9, 6e9, 4.4e9
```

* cMORB 三个系数取 `calibrating_June2026.K_class_list_3_components_recalibrated_June2026`
  的未舍入值（那是标定的产物，`617`/`38.0e-9`/`-4.00e-18` 只是它的手抄舍入）。
* `Surface pressure = 3.06072e9` 是 python 的
  `calc_simple_P_profile(2700, 3400, h_moho=40e3)` 在盒顶（100 km）给的值
  （`2700·9.81·40e3 + 3400·9.81·60e3 = 3.06072e9`，精确）。
* 实测效果（`TEQ_PRESSURE_VERIFICATION.md`）：这几行对齐后，
  **外部 ASCII** IC 的 step-0 失配从 **2.16e-6 → ~1e-11**。
* 注意：内建插件**不需要**这些精确值也能让 IC 落在平衡态上——
  它用 `adiabatic_conditions->pressure()` 和同一份 `.prm` 参数，
  舍入与否都与运行时间步完全一致。

### 11.2 python 侧有一个必须先修的陷阱：**两份同名、内容不同的 K 列表**

| 模块 | cMORB 的 (T_m0, A, B) | 谁在用 |
|---|---|---|
| `calibrating_June2026.K_class_list_3_components_recalibrated_June2026` | `617.27680612617 / 3.799177670844276e-8 / −4.001653781719964e-18` | `make_uniform_phase_ic_case6.py`、`equilibrium_postprocessor.py`、`melting_process.py`、`single_point_tests.py` |
| `component_data.K_class_list_3_components_recalibrated_June2026`（**同名**） | `617 / 38.0e-9 / −4.00e-18` | `make_initial_conditions.py:15`（它的 `import` 是从 `component_data` 来的），因而 `calc_c_from_T_P_and_f_3_components_case6()` 与旧 `initial_composition_for_carbon_case06.txt` 用的是这一份 |

磁盘上的两个 IC 文件正好各用一份（已用文件数值双向佐证），
所以"IC 与 `.prm` 是不是同一条曲线"取决于**哪个脚本生成的**——
这就是 `make_uniform_phase_ic_case6.py` docstring 里那句警告的真实来源。

**建议**：把 `component_data` 那一份改成 `calibrating_June2026` 的别名（或直接改
`make_initial_conditions.py` 的 import），让 python 侧也只有一个来源。
这样"外部 ASCII 路径"和"内建路径"才能保证同样的物理。

### 11.3 每条路径的"自洽"定义（回归判据）

| 路径 | 自洽的含义 | 判据 |
|---|---|---|
| 内建 IC | `Φ(P_ad, T, c̄) = φ_stored` 到二分容差 | 插件 `Report initial equilibrium` 打印的 `max|phi_stored − phi*| ≤ 1e-11` |
| 外部 ASCII + `.prm` | python 的 K/P/T 与 `.prm` 逐位一致 | 用 python 在 IC 的 (P,T,c̄) 上正向重解，`max|Δφ| ≤ 1e-11`（现在是 2.16e-6） |
| 两条路径互相对照 | 同一份 `.prm` 下两条路给出同一个 c̄ | `Dump generated state to file` 的输出与 `initial_composition_...txt` 逐点比对，`max|Δc̄| ≤ 1e-12` |

---
