# 内置初始条件模块：实现之后怎么写 IC —— 以 `case6_05_uniphase` 为例

> 本文描述的是**实现之后**的用法（P1–P3 还没做）。
> 设计细节见 `INTERNAL_EQUILIBRIUM_IC_DESIGN.md`，已完成的 P0 见其 §9.1。

---

## 0. 一句话对比

| | 今天 | 实现后 |
|---|---|---|
| 全岩成分 c̄ | 跑 `python run/make_uniform_phase_ic_case6.py` → `initial_composition_for_carbon_case06_uniphase.txt`（4.5 MB，30002 行） | **不需要** |
| 温度 T | 同一个脚本 → `initial_temperature_for_carbon_case06_uniphase.txt`（0.87 MB） | **不变**（温度暂不内置） |
| `.prm` 的 composition 段 | `set Model name = ascii data` + `Data file name` | `set List of model names = thermodynamic equilibrium` + 一族参数 |
| φ / c_l / c_s | 由 python 写好、再由 ASPECT 覆盖 | ASPECT 自己算，**与运行时间步用同一内核、同一份参数** |
| 边界节点 | P0 之前被 IC 值钉住 | 自动与内部平衡一致 |

工作流从「改 python → 跑 python → 拷 txt → 跑 ASPECT」变成「改 `.prm` → 跑 ASPECT」。

---

## 1. 实现后的 `.prm`（把 case6_05_uniphase 完整搬进内核）

````prm
##################### 场定义：不变 ########################
subsection Compositional fields
  set Number of fields = 11
  set Names of fields = porosity, dunite, morb, cmorb, dunite_liquid, morb_liquid, cmorb_liquid, melting_rate, debug_1, debug_2, debug_3
  set Compositional field methods = prescribed field, field, field, field, prescribed field, prescribed field, prescribed field, prescribed field, prescribed field, prescribed field, prescribed field
end

##################### 温度：不变（暂时不内置）#############
subsection Initial temperature model
  set List of model names = ascii data
  subsection Ascii data model
    set Data directory = ./
    set Data file name = initial_temperature_for_carbon_case06_uniphase.txt
  end
end

##################### 初始成分：新插件 #####################
subsection Initial composition model
  set List of model names = thermodynamic equilibrium

  subsection Thermodynamic equilibrium
    # ── 我要指定什么 ────────────────────────────────────────
    set Specification  = melt fraction
    set Mixing family  = tie line
    set Pressure source = adiabatic
    set Solving tolerance = 1e-14

    # ── 背景全岩成分 c̄_bg（200 ppm CO2）──────────────────────
    subsection Background bulk composition
      set Model = uniform
      set Uniform composition = 0.7, 0.2993333333333333, 0.0006666666666667
    end

    # ── 我想要什么熔融结构 ──────────────────────────────────
    subsection Target melt fraction
      set Model = background            # 基准 = 背景平衡值 f_bg(z)
      subsection Solitary wave
        set Add solitary wave        = true
        set Peak position            = 25000   # ASPECT 原生 y [m]（= 225 km 深）
        set Peak porosity            = 3.0e-3  # 绝对峰值孔隙度
        set Truncation threshold     = 1.0e-4  # (phi_r - 1) < tol 之外不嵌入
        set Derive compaction length = true    # δ、A、窗口全部导出
      end
    end
  end
end

##################### 绝热参考剖面：不变 ###################
subsection Adiabatic conditions model
  set Model name = compute profile
  subsection Compute profile
    set Composition reference profile = initial composition
    set Number of points = 2000
    set Use surface condition function = false
  end
end

##################### 边界条件：不变 #######################
subsection Boundary composition model
  set Fixed composition boundary indicators = top, bottom
  set List of model names = initial composition
end
````

**注意**：`Boundary composition model = initial composition` 现在会去问同一个插件，
所以 c̄ 的 Dirichlet 值自动等于内部算出的 c̄ 在边界点上的值——
"IC 与边界值不一致"这个问题从根上消失。
（`porosity`/`*_liquid`/`melting_rate`/`debug_*` 是 prescribed 场，
P0 之后完全不受边界条件影响。）

---

## 2. 逐条对照：python 变量 → `.prm` 参数

`make_uniform_phase_ic_case6.py` 与上面 `.prm` 的对应关系：

| python | 值（case6_05） | `.prm` |
|---|---|---|
| `DY` / `N_POINTS` | 10 m / 15001 | **消失**（直接在网格节点上算，不需要文件分辨率） |
| `C_CO2_PPM_BG = 200`、`C_CO2_IN_CMORB = 0.30` | → `c_cmorb,bg = 200e-6/0.30 = 6.6667e-4` | `Uniform composition = 0.7, 0.2993333333333333, 0.0006666666666667` |
| `PHI_PEAK = 3.0e-3` | 波峰孔隙度 | `Peak porosity = 3.0e-3` |
| `WAVE_PEAK_Y = 25000.0` | 峰位（y，= 225 km 深） | `Peak position = 25000` |
| `WAVE_HALF_WINDOW = 10e3` | 只算窗口内 | **删除**，改为 `Truncation threshold = 1e-4`（窗口半径由它导出） |
| `WAVE_AMPLITUDE = 20.0` | 相对振幅（**摆设**：调用方随即 `/1e-3` 约掉） | **删除**，改为导出：`A = φ_peak / f_bg(peak)` |
| `WAVE_DELTA = 53.45` | 压实长度（**硬编码**，与 `calc_compaction_length` 算出的 683 m 不符） | **删除**，改为导出：由材料模型在峰处状态、`φ = φ_0` 求值给出 |
| `background_porosity=1e-3` | 只当缩放因子，**被约掉** | **删除** |
| `shape = clip((φ_r−1)/(A−1),0,1)` | 波形（多余的归一化） | 用 `φ_r ∈ [1,A]` 本身（孤立波习惯） |
| `f_target = max(f_bg + (φ_peak−f_bg(peak))·shape, f_bg)` | 目标（**加法**） | `φ*(x) = f_bg(x)·φ_r(x)`（**乘法**，见 §3.1） |
| `c̄ = (1−f_target)·c_s,bg + f_target·c_l,bg` | tie-line 混合 | `Mixing family = tie line`（不变） |
| `mask = (shape>0) & (f_bg>0)` | 全固区排除 | **自动满足**：`f_bg = 0` 处 `φ* = f_bg·φ_r = 0` |
| `build_initial_condition_PT_profile` 的 P | `2700/3400` 分段静岩 | `Pressure source = adiabatic`（= ASPECT 自己的参考剖面） |
| 同一个函数的 T | `PracticalTemperatureProfile` | **不变**，仍由 `Initial temperature model` 提供 |
| `K_class_list_3_components_recalibrated_June2026` | python 的 K | `.prm` 的 `Melting point/coefficient A/B`（**同一个内核，不可能不一致**） |

---

## 3. 内部逐点算的是什么（用本算例的真实数字）

对每个网格节点 `y`：

```
1. P = adiabatic_conditions->pressure(y)          # = 3.06e9 + 3400·9.81·(150000 − y)
2. T = initial_temperature_manager->initial_temperature(y)   # 读你的温度 IC
3. c̄_bg = (0.7, 0.29933333, 6.66667e-4)                      # 背景
4. (f_bg, c_s,bg, c_l,bg) = Φ(P, T, c̄_bg)                     # 背景平衡
5. 波形参数**全部导出**（只需给 Peak position / Peak porosity / Truncation threshold）：
     φ_0 = f_bg(25000) = 1.054483e-04                 # 峰处的背景孔隙度
     用材料模型在 (P,T,φ=φ_0) 求值 → k, ξ, η, μ_f
     δ = sqrt(k·(ξ + 4η/3)/μ_f) = 34.24 m             # 压实长度
     A = φ_peak/φ_0 = 3.0e-3/1.054483e-4 = 28.450     # 相对振幅
     r = δ·implicit_function(1 + tol) = 2.42 km       # 窗口半径（tol = 1e-4）
6. φ*(y) = f_bg(y)·φ_r(|y−25000|/δ)   （|y−25000| > r 时 φ* = f_bg）   ← **乘法**
7. s = φ*(y)        ← tie-line 族下这**就是**解析解（Φ(s) ≡ s 精确成立）
   c̄(y) = (1−s)·c_s,bg + s·c_l,bg
8. (φ, c_s, c_l) = Φ(P, T, c̄)                                 # 回代，写出去
```

导出量会在 IC 阶段打印：

```
[TEQ solitary wave] phi_0 = f_bg(peak) = 1.054483e-04   (y_peak = 25000 m, 225.0 km)
                    material model at (P=7.2293 GPa, T=1579.83 K, phi=phi_0):
                      k = 1.172616e-20 m^2, xi = 2.142410e+20 Pa s,
                      eta = 2.142410e+20 Pa s, mu_f = 5.000000e-03 Pa s
                    delta = 34.24 m,  A = 28.450,  window = +-2419.4 m  (tol = 1e-4)
                    FWHM = 2853.9 m
```

实测（`prototype_mixing_families.py` / `prototype_equilibrium_ic.py` 用的同一套内核）：

| `y` | 深度 | P | T | f_bg | c_s,bg | c_l,bg |
|---|---|---|---|---|---|---|
| 25000（峰） | 225 km | 7.22925 GPa | 1579.83 K | 1.0545e-4 | (0.7000687, 0.2993618, 5.6951e-4) | (0.0488463, 0.0292559, 0.9218978) |
| 50000 | 200 km | 6.39540 GPa | 1572.10 K | 1.1301e-4 | (0.7000736, 0.2993628, 5.6359e-4) | (0.0487711, 0.0386174, 0.9126116) |
| 100000 | 150 km | 4.72770 GPa | 1556.75 K | 1.4491e-4 | (0.7000943, 0.2993657, 5.3997e-4) | (0.0491069, 0.0760133, 0.8748798) |
| 140000 | 110 km | 3.39354 GPa | 1303.00 K | **0** | (0.7, 0.2993333, 6.6667e-4) | (0.0079837, 0.0101733, 0.1381902) |

波峰处：

```
φ* = 3.000e-03            ⟹  s = 0.003000000000
c̄  = (0.698115, 0.2985515, 0.0033335)     Σ = 1.000000000000000
   → 1000.05 ppm CO2（背景 200.0 ppm）
回代 Φ(c̄) = 3.0000000000000009e-03        |Φ − φ*| = 8.7e-19
```

三条会立刻看到的结果：

1. **窗外逐位等于背景**：`y = 100000` 处 `max|c̄ − c̄_bg| = 1.1e-16`。
   波只在 ±10 km 窗口内存在，与 python 的 mask 行为逐位一致。
2. **窗内相成分恒等于背景**：波峰处回代得到
   `c_s = c_s,bg`、`c_l = c_l,bg`（上表数值逐位相同）——
   这就是 tie-line 族的"均匀相"性质：**波只携带熔融分数异常，不携带任何相成分异常**。
3. **全固区保持背景**：`y = 140000`（110 km，`f_bg = 0`）→ 目标 φ* = 0 →
   `c̄ = c̄_bg`、`φ = 0`。

精度：整个 IC 与"ASPECT 在节点上重算的平衡值"一致到二分容差
（`Solving tolerance = 1e-14`），而不是现在外部流程的 **2.16e-6**。

### 3.1 孤立波是怎么"挂"到随深度变化的背景上的

背景之所以"不水平"，是因为 `P(z)`（静岩）与 `T(z)`（地温）都随深度变，
于是平衡态 `(f_bg, c_s,bg, c_l,bg)` **逐深度都不同**——每个深度有**自己的一条 tie-line**。

拿本算例实测（200 ppm CO2 背景）：

| 深度 | P [GPa] | T [K] | `f_bg` | `c_s,bg` (dun, morb, cmorb) | `c_l,bg` | tie-line 方向 `c_l−c_s` |
|---|---|---|---|---|---|---|
| 140 km | 4.3942 | 1553.70 | 1.575e-04 | (0.7001025, 0.2993663, 5.3117e-04) | (0.049357, 0.089925, 0.860719) | (−0.650746, **−0.209442**, +0.860188) |
| 150 km | 4.7277 | 1556.75 | 1.449e-04 | (0.7000943, 0.2993657, 5.3997e-04) | (0.049107, 0.076013, 0.874880) | (−0.650987, **−0.223352**, +0.874340) |
| 175 km | 5.5616 | 1564.41 | 1.246e-04 | (0.7000811, 0.2993641, 5.5480e-04) | (0.048767, 0.052574, 0.898658) | (−0.651314, **−0.246790**, +0.898104) |
| 200 km | 6.3954 | 1572.10 | 1.130e-04 | (0.7000736, 0.2993628, 5.6359e-04) | (0.048771, 0.038617, 0.912612) | (−0.651303, **−0.260745**, +0.912048) |
| 225 km | 7.2293 | 1579.83 | 1.054e-04 | (0.7000687, 0.2993618, 5.6951e-04) | (0.048846, 0.029256, 0.921898) | (−0.651222, **−0.270106**, +0.921328) |
| 250 km | 8.0631 | 1587.60 | 1.002e-04 | (0.7000652, 0.2993611, 5.7371e-04) | (0.048944, 0.022699, 0.928358) | (−0.651122, **−0.276663**, +0.927784) |

读数：

* **固相端几乎不动**：`c_s,bg` 的 dunite 从 0.7000652 变到 0.7001025（5e-5 相对）；
* **液相端动得很厉害**：`c_l,bg` 的 morb 从 0.0227（250 km）到 0.0899（140 km），**4 倍**；
* 所以 tie-line 这个"线段"是**逐深度在转动、伸缩**的——这就是"背景不水平"的确切含义。

### 3.1.1 无量纲形状的约定

按孤立波理论习惯，无量纲形状就是 `φ_r ∈ [1, A]`——**尾部趋于 1（背景），峰处等于 A**。
不用再归一化成 `[0,1]`；嵌入是**乘法**：`φ*(z) = f_bg(z)·φ_r(z)`。
（旧 python 里那个 `shape = (φ_r−1)/(A−1)` 只是为了配合"加法 + 独立设定峰值"的写法，
现在峰值由 `φ_peak = f_bg(z_peak)·A` 唯一确定，归一化就多余了。）

### 3.1.2 置入方式

我们要的不是"在成分空间里加一个固定的扰动"，而是
"**在每个深度上，沿着那个深度自己的 tie-line 走 `Δφ` 那么远**"。写成公式：

```
c̄(z) = (1 − φ*(z))·c_s,bg(z) + φ*(z)·c_l,bg(z)
     = c_s,bg(z) + φ*(z)·[c_l,bg(z) − c_s,bg(z)]

⟹  Δc̄(z) = c̄(z) − c̄_bg(z) = Δφ(z) · [ c_l,bg(z) − c_s,bg(z) ]
   其中 Δφ(z) = φ*(z) − f_bg(z) = f_bg(z) · (φ_r(z) − 1)
```

注意 `Δφ` 是**乘法**形式（当地背景 × 当地无量纲振幅减 1），
不是"峰值超出量 × 归一化形状"的加法形式。两者在峰处相同，
但乘法形式让 `f_bg → 0` 处 `Δφ → 0` 自动成立（全固区无需 mask），
且波幅随当地背景自然缩放。本算例波很窄（±2.4 km），
两种形式在数值上差别很小（背景在 2 km 内只变 ~6%）。

即：**扰动的"大小"是熔融分数超出量 `Δφ(z)`，扰动的"方向"是当地 tie-line 向量
`c_l,bg(z) − c_s,bg(z)`，两者都随深度变。** 波峰处

```
Δφ = 3.0e-3 − 1.0545e-4 = 2.8946e-3
c_l−c_s = (−0.651222, −0.270106, +0.921328)
Δc̄ = (−1.8850e-3, −7.8185e-4, +2.6669e-3)
```

**为什么这样就"挂"住了**：因为 `Φ(s) ≡ s` 沿 tie-line **精确**成立
（设计稿 §5.4 有代数证明：把 `c_s = K c_l` 代进平衡方程，
左端恒等于 `Σc_s − Σc_l`，而平衡条件就是它 = 0）。所以

1. **窗外**：`shape = 0` ⟹ `Δφ = 0` ⟹ `c̄ = c̄_bg(z)`，**逐位**回到背景
   （实测 `max|Δc̄| = 1.1e-16`）；
2. **窗内**：回代得到的 `φ` 精确等于 `φ*(z)`，
   而 `c_s, c_l` 精确等于**当地**的 `c_s,bg(z), c_l,bg(z)`——
   波**只携带熔融分数异常，不携带任何相成分异常**；
3. 因为每一步都是相对**当地**背景做的，波在整条熔化区里都"贴着"背景，
   不会在某个深度上偏离 tie-line。

**和"直接在成分空间加扰动"的区别**（实测，同一套内核）：

| 做法 | 波峰 | 峰外 1 km（y=24000） | 峰外 2 km | 相成分扰动 |
|---|---|---|---|---|
| **A：逐深度 tie-line**（本案） | err **8.7e-19** | err **8.7e-19** | err **−2.2e-19** | `Δc_s = 0`，`Δc_l = 0`（精确） |
| B：把峰处的 `Δc̄` 固定下来、只按 `shape(z)` 缩放 | err 8.7e-19 | err **−7.6e-07** | err **−4.4e-07** | `Δc_s = 5.1e-07`，`Δc_l = 4.9e-08` |

B 在峰处也对（构造使然），但一离开峰值就既偏了熔融分数、又把相成分搞脏了。
本算例的波很窄（±2 km），背景在这么短的距离里变化不大，所以 B 的误差只有 3e-4 相对量级；
波越宽、跨的深度越大，B 就越不可用——而 A **在任意宽度下都是精确的**。

### 3.1.3 波宽是导出的 ⟹ 对背景成分很敏感

`δ` 由材料模型给出，而 `k ∝ φ_0³`、当前 `.prm` 又设了
`Exponential melt weakening factor = 0`（`ξ = ξ_0·exp(−α_φ φ)` 退化为常数），
所以 **`δ ∝ φ_0^{3/2}`**。实测（峰位 225 km，`ξ = η = 2.1424e20`、`μ_f = 5e-3`、`K_ref = 1e-8`）：

| 背景 CO₂ | `φ_0 = f_bg(peak)` | `δ` [m] | `A = φ_peak/φ_0` | 窗口半径 @tol=1e-4 | FWHM |
|---|---|---|---|---|---|
| **200 ppm（case6）** | 1.0545e-04 | **34.2** | 28.45 | **2.42 km** | 2854 m |
| 300 ppm | 4.672e-04 | 319 | 6.42 | 8.35 km | 6805 m |
| 500 ppm | 1.191e-03 | 1299 | 2.52 | 25.7 km | 14329 m |
| 1000 ppm | 3.000e-03 | 5195 | **1.000** | —（A=1，波不存在） | — |
| 2000 ppm | 6.617e-03 | 17019 | **0.453** | —（A<1，非法） | — |

背景 CO₂ 从 200 提到 500 ppm，波宽涨 **38 倍**——物理上自洽（熔体多 ⟹ 渗透率高 ⟹
压实长度长 ⟹ 波更宽），但意味着**波的形态不再由你直接控制**。
1000 ppm 时 `φ_peak` 恰好等于背景孔隙度，`A = 1`，波直接消失；再多就不合法（`A < 1`）。

所以默认全部导出，但留一个逃生口给"要一条宽度受控的波"的对照实验：

```prm
        set Derive compaction length = false
        set Compaction length        = 53.45   # m
```

（`A` 不需要逃生口——它由你给的**绝对**峰值孔隙度唯一确定，这正是你要的语义。）

另外三条必须在实现里加守卫：

* `φ_0 ≤ 0`（峰位在固相线以下）→ **报错**，提示把峰位放进熔化区或提高背景不相容组分；
* `A ≤ 1` → **报错**（波不存在）；
* 窗口半径超出模型域或超出熔化区 → **只警告**（因为 `f_bg = 0` 处 `φ* = 0`，构造本身安全），
  但要打印 `r` 与熔化区范围供判断。

---

**最后一句必须说清楚**：这个构造只借用孤立波**无量纲的形状** `shape(z) ∈ [0,1]`，
它**不是**分层介质中压实波方程的解。真正的"孤立波"是均匀背景下的解析解
（`SolitaryWaveSolution`），这里只是把它的形状当作一条**目标熔融分数剖面**来用，
再由平衡映射反推出全岩成分。所以这个 IC 是**运动学指定**，不是动力学自洽的初值——
它随后会被自身的压实/迁移演化掉，这正是你想研究的对象。

---

### 3.1.4 附：case6 的背景孔隙度为什么是 ~1e-4 而不是 2e-4

（2026-09-30 核查，因为新配方里 `φ_0 = f_bg(peak)` 直接决定 `δ` 和 `A`。）

**结论：case6 的背景从来就是 1.0–1.6e-4，没有"被改坏"；2.0004e-4 是 region3 的
`mr_*`（mr_06 / mr_new_F）算例。** 实测所有 IC 文件（非波区 `|y−25000| > 20 km`）：

| IC 文件 | 非波区 φ 范围 | 波峰 |
|---|---|---|
| `initial_composition_for_carbon_case06.txt`（旧） | 0 – 1.6148e-04 | 3.000e-03 |
| `..._case06_uniphase.txt`（10 m） | 0 – 1.5848e-04 | 3.000e-03 |
| `..._case06_uniphase_5m / _20m.txt` | 0 – 1.5849e-04 | 3.000e-03 |
| `region3_carbon_test/..._mr_new_F*.txt` | **2.000352e-04 – 2.000352e-04（平的）** | 3.000e-03 |

`postprocess/solitary_wave/report.md:157` 对 case6 的描述也是"约 1.0 – 1.6 × 10⁻⁴"。

**两者不是同一套物理**：

| | `mr_*`（region3，零熔融标定） | case6 |
|---|---|---|
| 熔融曲线 | `A = B = 0`（与压力无关） | `A, B ≠ 0`（与压力相关） |
| 温度 | 常数 1269.89 K（996.74 °C） | 地温 1219–1588 K |
| 压力 | 4.39 GPa | 3.06–8.06 GPa（静岩） |
| 背景 c̄ | (0.7, 0.2993333, 6.6667e-4) | 同 |
| 背景 φ | 2.0004e-4（平） | 1.0–1.6e-4（随深度） |

把 case6 的 P/T 剖面配不同的熔融曲线，实测（200 ppm 背景）：

| 曲线 | 140 km | 225 km | 250 km |
|---|---|---|---|
| 现在（June2026，`A,B≠0`） | 1.575e-04 | **1.054e-04** | 1.002e-04 |
| `.prm` 里的 `cmorb_30`（617/38e-9/−4e-18） | 1.615e-04 | 1.094e-04 | 1.041e-04 |
| 更早的 `cmorb_20`（640/30.1e-9/−1.88e-18, P_thr=6e9） | **0** | **0** | **0**（全列不熔） |
| 假如 `A = B = 0`（压力无关） | 0.311 | **0.324** | 0.328 |

读数：

* **是"压力相关性"把背景压到 1e-4 的**：7.23 GPa 把熔点抬到
  morb 1000 → 1638.7 °C、dunite 1780 → 2003.8 °C、cmorb 617 → 714.8 °C，
  而 225 km 的温度只有 1306.68 °C —— 只有 cmorb 在熔。
  同一温度下若换成压力无关曲线，会是 **32 % 熔融**。
* **化学参数确实调整过，但方向相反**：早期用 `cmorb_20` 时 case6 这条 P/T
  剖面**根本不熔**（熔化区都不存在）；换成 `cmorb_30`/June2026 之后才有 1e-4。
* **温度确实能解释 1e-4 → 2e-4**：要让 225 km 的背景到 2e-4，
  T 需从 1306.68 °C 升到 **1328.55 °C（+21.9 K）**；
  但 case6 的温度剖面一直是现在这条，从没高过 22 K。
* 另一个可能的混淆源：`postprocess/solitary_wave/report.md` 提到 case6 在
  LAB 附近（≈138 km）有个"浅部小包 φ ≈ 1.9e-4"——那个数接近 2e-4，
  但它是**局部结构**（且多半是演化后的场），不是背景；
  当前 IC 在 125–150 km 的极大值只有 1.585e-4 @ 140 km。

---

### 3.1.5 你要找的"背景 ≈1.5e-4 @225 km"的 IC：`..._984ppm` 那一族

**找到了，而且它精确对应你当初的设计意图（相对振幅 20、绝对峰值 3e-3 ⟹ 背景 1.5e-4）。**

两份内容相同的 IC（非波区取样 → 225 km 处内插）：

| IC 文件 | 背景 @235 km | 背景 @215 km | **背景 @225 km** | 峰值 @225 km | **比值** | 背景/峰值 CO₂ |
|---|---|---|---|---|---|---|
| `melt_equilibrium_comparison/bench_aspect/initial_composition_for_carbon_case6_984ppm.txt` | 1.4636e-04 | 1.5186e-04 | **1.4911e-04** | 3.000e-03 | **20.12** | 200 / 983.7 ppm |
| `region5_carbon/case5_practical_0624/initial_composition_for_carbon_case5_984ppm.txt` | 1.4636e-04 | 1.5186e-04 | **1.4911e-04** | 3.000e-03 | **20.12** | 200 / 983.7 ppm |
| 现在的 `..._case06_uniphase.txt`（对照） | 1.0511e-04 | 1.0991e-04 | 1.0751e-04 | 3.000e-03 | **27.9** | 200 / 999.0 ppm |

配对的地温文件是 `initial_temperature_for_carbon_case6_141km.txt`
（`bench_aspect/base.prm:183`）与 `..._case5_141km.txt`。

**差别就是地温，不是化学**：反推两个 T 文件的 `T_LAB`（用 `T_p = T_LAB·exp(−αgh_LAB/C_p)`）：

| T 文件 | `T_potential` | `h_LAB` | 反推 `T_LAB` | T @225 km |
|---|---|---|---|---|
| `..._141km.txt`（984ppm 那族用的） | 1520.462 K | 141.0 km | **1563.11 K** | 1589.09 K |
| `..._case06_uniphase.txt`（现在用的） | 1511.605 K | 139.99 km | **1553.70 K** | 1579.83 K |

`1563.15 K` 正是 `carbon_melting/legacy/reference_equilibrium_calculate/practical_T_profile.py:5`
的默认 `T_LAB`；而现在的 `build_initial_condition_PT_profile()` 显式传了
`T_LAB=1553.70`（`equilibrium/practical_T_profile.py` 的默认则是 1573.15）。
**`T_LAB` 降了 9.4 K ⟹ 225 km 处温度降 9.26 K ⟹ 背景孔隙度降 28 %。**

实测核对（200 ppm，225 km，`P = 7.2293 GPa`）：

```
老地温 T = 1315.94 °C  ->  f_bg = 1.4715e-04,  A = 20.39   （984ppm IC: 1.491e-4, A = 20.1 ✓）
新地温 T = 1306.68 °C  ->  f_bg = 1.0545e-04,  A = 28.45
```

次要因素：那两份 IC 用的 cMORB 系数是 `component_data` 的 `617/38.0e-9/−4.00e-18`，
uniphase 用的是 `calibrating_June2026` 的 `617.2768/3.79917767e-8/−4.00165378e-18`，
贡献约 2–3 %（远小于温度的 40 %）。

#### ⚠️ 一个必须知道的强敏感性

在当前地温下，225 km 的背景孔隙度对背景 CO₂ **极其敏感**（固相线附近的"悬崖"）：

| 背景 CO₂ | 200 | **212.3** | 250 | 300 | 500 | 1000 |
|---|---|---|---|---|---|---|
| `f_bg(225 km)` | 1.054e-04 | **1.500e-04** | 2.863e-04 | 4.672e-04 | 1.191e-03 | 3.000e-03 |
| `A = 3e-3/f_bg` | 28.45 | **20.0** | 10.48 | 6.42 | 2.52 | 1.00 |

也就是说：

* 想用**当前**地温 + 200 ppm 复刻"背景 1.5e-4 / A = 20"，只需把背景 CO₂ 从
  200 ppm 提到 **212.3 ppm**（+6 %），`f_bg` 就翻了 42 %；
* 或者把地温换回老的 `T_LAB = 1563.15 K`，200 ppm 下自然就是 `A = 20.4`；
* 反过来，在新配方（`δ ∝ φ_0^{3/2}`）下这也意味着**波宽对背景 CO₂ 的敏感性被放大**：
  200 → 212 ppm 就让 `δ` 涨 1.75 倍。

---

## 4. 一个必须知道的细节：绝热剖面建好**之前**，插件会被调用一次

`.prm` 里的 `set Composition reference profile = initial composition`
（`compute_profile.cc:158`）让 ASPECT 在建绝热参考剖面时，
对每个参考点调用 `initial_composition_manager->initial_composition(p, c)`
——也就是调用你的新插件。这个调用发生在 `adiabatic_conditions->initialize()`
（`core.cc:372`）内部，**早于 `set_initial_temperature_and_compositional_fields()`（1999）**，
而此刻 `adiabatic_conditions->is_initialized() == false`。

插件的行为定死为：

* `is_initialized() == false` ⟹ 直接返回**背景成分 c̄_bg**（它不需要 P），
  派生场返回 0；
* `is_initialized() == true` ⟹ 走完整的反演流程。

对本算例**没有任何影响**：材料模型的密度是 `3400·(1−αΔT)·(1+βΔP)`，
而 `Thermal expansion coefficient = 0`、`Solid compressibility = 0`，
所以参考密度恒为 3400，与参考成分无关。
（若将来密度依赖成分，用背景成分做参考剖面本来也是合理选择。）

**想彻底避开这次调用**，把那一行换成

```prm
    set Composition reference profile = function
    subsection Function
      set Function expression = 0.7; 0.3-0.0006666666666667; 0.0006666666666667
    end
```

即可——参考剖面与 IC 插件完全解耦。

---

## 5. 跑完怎么确认

插件会在 `signals.post_set_initial_state` 打印一次汇总（`Report initial equilibrium = true`）：

```
[TEQ initial condition] specification = melt fraction, family = tie line, P = adiabatic
  nodes                       : 30001   (interior 29599)
  points with phi* = 0        : 7313    (f_bg = 0, 深度 < ~136.6 km)
  points solved analytically  : 22688   (tie-line 族；含窗口内 4001 个点)
  points solved by bisection  : 0
  points unreachable          : 0
  max |Phi(c_bar) - phi*|     : 8.7e-19
  max |Sigma c_bar - 1|       : 1.1e-16
  max |phi_stored - phi*|     : ~1e-15   (内存 double；读 .vtu 是 float32，只能到 ~1e-10)
  cMORB at wave peak          : 3.3335e-03   (1000.1 ppm CO2)
```

`max |phi_stored − phi*|` 就是"IC 到底落没落在 ASPECT 的平衡态上"的权威判据。
另外 `set Dump generated state to file = true` 会写一份 ASCII，
可以和你现在的 `initial_composition_for_carbon_case06_uniphase.txt` 逐点 diff——
预期最大差 ~1e-6（差异全部来自那条 0.72 MPa 的 `Surface pressure` 舍入，
以及 python 侧用的是另一份 cMORB 系数），而不是"文件读错"。

---

## 6. 想改东西时改哪里

| 我想…… | 改什么 |
|---|---|
| 换背景 CO₂ 浓度（比如 500 ppm） | `Uniform composition` 的第三个分量 = `500e-6/0.30 = 1.6667e-3`，第二个 = `0.3 − 1.6667e-3` |
| 换波峰位置（比如 200 km ⟹ y = 50000） | `Peak position = 50000` |
| 换波峰幅度 / 宽度 | `Peak porosity` / `Compaction length`、`Relative amplitude` |
| 不要窗口截断 | `Half window = 0`（表示不限窗） |
| **直接在固相线下方造熔融** | 换族：`Mixing family = background to component cmorb`（tie-line 族在那里无效，见设计稿 §5.4.2） |
| 目标熔融分数用一条外部剖面 | `Target melt fraction` 的 `Model = depth profile` + 一个只含 `y, phi` 的文本（**y 是 ASPECT 原生坐标，自下而上**） |
| 干脆自己指定 c̄ | `Specification = bulk composition`，`Background bulk composition` 换成 `depth profile`/`function` |
| 指定相成分 + 熔融分数反推 c̄ | `Specification = phase compositions`（会做质量平衡 + 一致性校验，不自洽就报错） |
| 换 K 曲线 / 熔点 | 只改 `.prm` 里那几行——**IC 自动跟着变**，不会再出现"两边不是同一条曲线" |

---

## 7. 外部 ASCII 路径仍然可用

两条路并行（你的决定）。想回到今天的方式，`Initial composition model` 照旧写：

```prm
subsection Initial composition model
  set Model name = ascii data
  subsection Ascii data model
    set Data directory = ./
    set Data file name = initial_composition_for_carbon_case06_uniphase.txt
  end
end
```

要让两条路给出同一个结果，需要（见设计稿 §11）：
`.prm` 用未舍入的 cMORB 系数与 `Surface pressure = 3.06072e9`，
并把 python 侧那两份同名 K 列表统一。
**已存在的算例按你的要求一律不动**，这些只用在新建的 deck 上。

---

## 8. 现状：**P1–P3 已实现并验收**

内置模块已经可用。上面的 `.prm`（§1）就是能跑的版本，唯一的差别是：

* `Peak position` 用 **ASPECT 原生竖直坐标**（2D 盒子里是 y，自下而上），
  所以 225 km 深写 `25000` —— 与 §1 一致；
* `Target melt fraction model` 目前支持 `background | uniform | depth profile | function`；
* 背景成分目前支持 `uniform | depth profile`（`function` 待补）；
* `Pressure source` 只有 `adiabatic`；背景成分的 `function` model、压力源选项、
  `vector` / `component i to j` 族待补（P4.3–P4.5）。
* **`Specification = phase compositions` 已实现**（P4.1）：给定固相/液相成分与
  熔融分数，由相质量平衡算 c̄，再回代校验；不自洽时按
  `If inconsistent = error | keep melt fraction | keep liquid` 处理。
* **`Dump generated state to file` 已实现**（P4.2）：把生成的状态写成
  `# POINTS: nx ny` 的 ASCII，可直接与外部产物比对、也能被 `ascii data` 读回去。

验收算例（**没有改动任何已有算例**）：

```
region5_carbon/case6_practical_0724/verify_internal_ic/
  internal_ic.prm     # 与 verify_teq_pressure/adiabatic.prm 只差输出目录和 composition 段
  run.log             # 插件的参数与自检输出
  analyze.py          # 三项独立验证
  README.md           # 完整说明
```

关键数字（详见该目录 `README.md`）：

| 检查 | 结果 |
|---|---|
| 生成场 vs 独立重算的 `Phi(P,T,c_bar)` | **5.3e-10**（受 vtu float32 限制；插件内存残差 **7.2e-13**） |
| 生成场 vs 解析孤立波目标 `f_bg·φ_r` | 1.2e-08 |
| **窗外**（`|y−25000|>15 km`）与外部 ASCII IC 逐点比较 | **逐位相同** |
| 波峰 | φ = 3.000000e-03，A = 28.45，**δ = 34.2384 m（导出）**，窗口 ±2419.42 m |
| 3 步（45 yr）含 IC 生成 | 16 s（8 MPI） |

插件在 IC 阶段自己打印导出的波形参数：

```
[TEQ initial condition] solitary wave parameters:
    peak position (vertical coordinate) = 25000 m
    background melt fraction at peak   = 0.000105448
    relative amplitude A               = 28.45
    compaction length delta            = 34.2384 m
    window radius (threshold 0.0001)      = 2419.42 m
[TEQ initial condition] generated thermodynamic equilibrium initial state:
    specification = melt fraction, mixing family = tie line, pressure = adiabatic
    evaluated points             : 1.48502e+06
    points solved for a target   : 1.12303e+06
    points kept at the background: 361986
    max |Phi(c_bar) - phi*|      : 7.18297e-13
    max |sum(c_bar) - 1|         : 2.84217e-13
```

注意 `δ = 34.2384 m` 与外部生成器硬编码的 `53.45 m` 不同 —— 这是预期的：
新配方按 `δ = sqrt(k(φ_0)(ξ + 4η/3)/μ_f)` 在峰处背景状态下导出，
所以**窗内的波比原来窄**（FWHM 2.85 km vs 3.18 km），而**窗外逐位不变**。

**两条使用注意**（实现细节，见验证目录 README）：

1. 插件在**绝热参考剖面建好之前**会被调用一次（因为 deck 里
   `Composition reference profile = initial composition`）。此时它直接返回背景成分。
2. 内置插件与 `Boundary composition model = initial composition` 一起用时，
   边界条件会自动取到同一套内部平衡值 —— 这正是我们想要的，
   但插件在 step 1 之后还会被边界条件再调用（用的是缓存下来的温度 IC 指针）。
