# ASPECT 的数值方法：空间离散与时间离散

> 本文基于本仓库的 ASPECT 源码（`VERSION = 2.5.0`，`include/`、`source/`、`CMakeLists.txt`）以及
> `~/dealii-candi/deal.II-v9.5.1` 的 deal.II 源码整理。
> 文中所有结论都标注了文件与行号，便于自行核对。
> 阅读顺序建议：先看 §0 和 §6 的两张总表建立框架，再回头看细节。

---

## 0. 一句话总览

ASPECT 是一个**基于 deal.II 的、以连续/间断有限元做空间离散、以隐式 BDF 做时间离散**的地幔对流软件：

- **空间上**：2D 四边形 / 3D 六面体网格上的混合有限元（mixed FEM）。
  Stokes 用 Taylor–Hood 单元 $Q_k^d \times Q_{k-1}$（默认 $Q_2^d\times Q_1$），
  温度与成分默认用连续 $Q_2$，也可切成间断 $DGQ_2$。
  网格是**自适应的 h 型加密**（没有 p 自适应），误差指标来自物理量而非 Kelly 估计子。
- **时间上**：Stokes 方程**没有时间导数**（准静态，每个时间层解一次椭圆型鞍点问题）；
  只有温度和成分方程做时间推进，用的是**隐式变步长 BDF1/BDF2**，
  而速度与材料参数取**上一/上两个时间层的外推值**，因此整体是"半隐式/IMEX + 算子分裂"。
  时间步长由 **CFL 条件**控制。
- **稳定化**：连续单元用 **entropy viscosity**（熵黏性人工扩散）；间断单元用 **上风 + 对称内罚
  （SIPG）**，此外还有可选的保界限幅器。
- **线性求解**：Stokes 用分块 Schur 补预条件 + AMG（Trilinos ML）+ FGMRES；对流方程用 GMRES + ILU。
- **并行**：`parallel::distributed::Triangulation`（p4est）做分布式自适应网格。

**最重要的一点（很多人会误解）**：deal.II 虽然在 `include/deal.II/base/time_stepping.h` 里有
`TimeStepping` 命名空间（显式/隐式 Runge–Kutta 等），但 ASPECT **完全没有使用它**
（对 `source/`、`include/` 全量 grep `dealii::TimeStepping` 无命中）。
ASPECT 里的 `aspect::TimeStepping` 是它自己的**时间步长控制插件系统**（`source/time_stepping/`），
和 ODE 积分器无关。**ASPECT 的时间离散（BDF1/BDF2、外推、时间步控制）是手写在组装代码里的。**

---

## 1. 要求解的方程与总体策略

方程见 `doc/sphinx/user/methods/basic-equations/index.md`，代码里的实现对应关系：

**(a) 可压缩 Stokes（准静态、无 $\partial \mathbf u/\partial t$）**

$$
-\nabla\cdot\left[2\eta\left(\varepsilon(\mathbf u)-\tfrac13(\nabla\cdot\mathbf u)\mathbf 1\right)\right]+\nabla p=\rho\mathbf g,
\qquad
\nabla\cdot(\rho\mathbf u)=0 .
$$

**(b) 温度方程**

$$
\rho C_p\left(\frac{\partial T}{\partial t}+\mathbf u\cdot\nabla T\right)-\nabla\cdot k\nabla T
= \rho H + 2\eta\,\varepsilon':\varepsilon' + \alpha T\,\mathbf u\cdot\nabla p + \rho T\Delta S\left(\frac{\partial X}{\partial t}+\mathbf u\cdot\nabla X\right)
$$

**(c) 成分场（无扩散的对流方程）**

$$
\frac{\partial c_i}{\partial t}+\mathbf u\cdot\nabla c_i = q_i
$$

关键结构性事实：

1. **(a) 里没有时间导数**。整个模拟中，Stokes 不是"推进"出来的，而是在每个时间层上、
   用当前（或其外推/非线性迭代得到的）$T,c_i$ 重新解一次边界值问题。
   所以严格说 ASPECT 的时间离散只作用于 (b)(c)——见 `assemblers/stokes.cc` 里根本没有
   `time_step` 因子，而 `assemblers/advection.cc` 里到处都是 `time_step`。
2. **(b)(c) 是同一个"对流–扩散/纯对流"骨架**，代码里统一抽象为 `AdvectionField`，
   由 `assemble_advection_system()` 组装（`source/simulator/assembly.cc:1092`）。

---

## 2. 空间离散

### 2.1 网格：只允许四边形/六面体，纯 h 自适应

- 网格由 `geometry_model->create_coarse_mesh(triangulation)` 产生，随后断言
  **所有单元必须是超立方体**（四边形/六面体）：

  ```cpp
  // source/simulator/core.cc:455-457
  Assert (triangulation.all_reference_cells_are_hyper_cube(),
          ExcMessage ("ASPECT only supports meshes that are composed of quadrilateral "
                      "or hexahedral cells."))
  ```

- 初始加密：`initial_global_refinement` 次全局加密；每次先把所有单元标 refine，
  再让 mesh refinement 插件有机会"取消"若干单元（`core.cc:1957-1970`）。
- 初始自适应加密：`initial_adaptive_refinement` 次，每次重新解、重新估计误差再加密。

**几何映射（mapping）**由 `construct_mapping()` 决定（`core.cc:92-121`）：

| 情形 | 映射 |
|---|---|
| 几何模型有曲边/曲单元（如球壳） | `MappingQCache<dim>(4)`，即 4 次等参映射并缓存 |
| 平直 Box、无初始地形 | `MappingCartesian<dim>`（JxW 为常数，最快） |
| 有初始地形导致初始网格已变形 | `MappingQ1<dim>` |
| 开启自由表面/网格变形 | 上述映射被替换为 `MappingQ1Eulerian`（ALE） |

这对空间离散的影响：**高阶映射 + 低阶解的等参/超参组合**。曲面网格上求积和形函数梯度都
经映射 Jacobian 变换，因此 `MappingQ4` 是为了让球面几何的几何误差不污染解的精度。

### 2.2 有限元空间：单元类型与多项式阶数

核心代码在 `source/simulator/introspection.cc`，`construct_default_variables()`
（`introspection.cc:185-232`）：

```cpp
variables.push_back(VariableDeclaration<dim>("velocity",
                     std::make_shared<FE_Q<dim>>(parameters.stokes_velocity_degree), dim, n_velocity_blocks));

if (parameters.use_equal_order_interpolation_for_stokes == false)
  variables.push_back(VariableDeclaration<dim>("pressure",
                     internal::new_FE_Q_or_DGP<dim>(parameters.use_locally_conservative_discretization,
                                                    parameters.stokes_velocity_degree-1), 1, 1));
else
  variables.push_back(VariableDeclaration<dim>("pressure",
                     std::make_shared<FE_Q<dim>>(parameters.stokes_velocity_degree), 1, 1));

variables.push_back(VariableDeclaration<dim>("temperature",
                     internal::new_FE_Q_or_DGQ<dim>(parameters.use_discontinuous_temperature_discretization,
                                                    parameters.temperature_degree), 1, 1));
variables.push_back(VariableDeclaration<dim>("compositions",
                     internal::new_FE_Q_or_DGQ<dim>(parameters.use_discontinuous_composition_discretization,
                                                    parameters.composition_degree),
                     parameters.n_compositional_fields, parameters.n_compositional_fields));
```

对应的参数默认值（`source/simulator/parameters.cc:907-1030`）：

| 变量 | 默认 | 可选 |
|---|---|---|
| 速度 | $Q_2^d$（`FE_Q(2)`），`Stokes velocity polynomial degree = 2` | $Q_k^d$ |
| 压力 | $Q_1$（`FE_Q(1)` = $Q_{k-1}$） | `Use locally conservative discretization=true` → $P_{-(k-1)}$，即 `FE_DGP(k-1)`（单元间间断）；`Use equal order interpolation for Stokes=true` → $Q_k$ |
| 温度 | $Q_2$（`FE_Q(2)`） | `Use discontinuous temperature discretization=true` → `FE_DGQ(2)` |
| 成分 | $Q_2$（每成分一个块） | `FE_DGQ(p)`，$p\ge 0$，$DGQ_0$ 也合法 |

几个要点，都是"空间离散概念"层面的：

1. **Taylor–Hood 满足 LBB（inf–sup）条件**。代码注释明确写了这是为了满足
   Babuška–Brezzi 条件（`parameters.cc:907-923`），不要随意把速度阶数设为 1（会不稳定）。
2. **$P_{-(k-1)}$ 局部守恒**：压力用间断的 $P_{-(k-1)}$（`FE_DGP`）时，可证明每个单元
   $\int_{\partial K}\mathbf u_h\cdot\mathbf n = 0$，即逐单元质量守恒，代价是更多自由度
   （`parameters.cc:944-978`）。
3. **等阶 $Q_k$–$Q_k$**：默认不稳定，需要用 Dohrmann–Bochev 型压力稳定化。
   这一步在组装里看得最清楚（`assemblers/stokes.cc:425-439`）：

   ```cpp
   // If we are using the equal order Q1-Q1 element, then we also need
   // to put the stabilization term into the (P,P) block of the matrix:
   data.local_matrix(i,j) += ( - (one_over_eta * pressure_scaling * pressure_scaling *
                                  (scratch.phi_p[i] - average_pressure_shape_function[i]) *
                                  (scratch.phi_p[j] - average_pressure_shape_function[j]))) * JxW;
   ```

   即向 $(P,P)$ 块加一个 $\frac{1}{\eta}(q-\bar q, q'-\bar q')$ 型的压力 Laplacian 稳定项。
4. **连续 vs 间断（CG vs DG）**：$Q$ 类跨单元 $C^0$ 连续；$DGQ$ 完全间断、单元间只通过
   面上的数值通量耦合。这直接决定后面组装要不要面循环，以及稳定化用哪种机制。

### 2.3 求积（quadrature）

`introspection.cc:107-128`：

```cpp
quadratures.velocities         = reference_cell.get_gauss_type_quadrature<dim>(stokes_velocity_degree+1);
quadratures.pressure           = /* equal order ? k+1 : k */;
quadratures.temperature        = ... (temperature_degree+1);
quadratures.compositional_fields = ... (composition_degree+1);
quadratures.system             = ... (max({vel,temp,comp}) + 1);
```

- 都是 **Gauss–Legendre 张量积求积**（`get_gauss_type_quadrature`），阶数取"多项式阶数 + 1"，
  即对 $Q_p$ 单元用 $(p+1)$ 点一维 Gauss，足以精确积分 $2p$ 次多项式质量矩阵并留有余量。
- Stokes 组装实际用的就是 `quadratures.velocities`（`assembly.cc:704`）。
- **对流方程单独用更高的求积阶数**，因为 $u\cdot\nabla\varphi$ 的被积多项式次数更高：

  ```cpp
  // source/simulator/assembly.cc:1119-1133
  // We have to assemble the term u.grad phi_i * phi_j, which is
  // of total polynomial degree  stokes_deg + 2*temp_deg - 1 ...
  const unsigned int advection_quadrature_degree =
      advection_field.polynomial_degree(introspection) + (parameters.stokes_velocity_degree+1)/2;
  ```

- 熵黏性残差、时间步长估计用的是 `QIterated<dim>(QTrapezoid<1>(), degree)`（梯形法则迭代），
  即**在单元节点上取值**而非 Gauss 点（`helper_functions.cc:718`、`time_stepping/convection_time_step.cc:33`）。

### 2.4 弱形式与组装

组装用 deal.II 的 `WorkStream::run`（按单元并行 + 预条件式流水线），
`local_assemble_*` 把单元贡献写进 `data.local_matrix/local_rhs`，
`copy_local_to_global_*` 再用 `constraints.distribute_local_to_global` 叠加到全局矩阵。
Scratch/CopyData 的定义在 `source/simulator/assemblers/interface.cc`。

**Stokes 的弱形式**（`assemblers/stokes.cc:409-421`）：

$$
\underbrace{(2\eta\,\varepsilon(\mathbf u),\varepsilon(\mathbf v))}_{}
\;-\;\underbrace{(p,\nabla\cdot\mathbf v)}_{}
\;-\;\underbrace{(q,\nabla\cdot\mathbf u)}_{}
\;=\;(\rho\mathbf g,\mathbf v)
$$

```cpp
data.local_matrix(i,j) += ( (eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                            - (pressure_scaling * scratch.div_phi_u[i] * scratch.phi_p[j])
                            - (pressure_scaling * scratch.phi_p[i] * scratch.div_phi_u[j]) ) * JxW;
```

其中 `pressure_scaling` 是把速度/压力不同量纲归一化的对角标度
（`assembly.cc:677` 的 `compute_pressure_scaling_factor()`），
它同时出现在 $(U,P)$、$(P,U)$ 块和等阶稳定项里，是数值上很重要的一步。

RHS 是浮力项 `density * gravity * phi_u`（`stokes.cc:379`），
另可选弹性、指定膨胀、traction 面项等（`assemblers/*`）。

**对流–扩散的弱形式**（`assemblers/advection.cc:86-254`）。令 $\varphi_i$ 为该标量场的形状函数，
$L$ 为潜热 LHS 项。（代码里两边都乘了 $\Delta t$，下面为便于阅读已除以 $\Delta t$；
$\beta_\text{BDF}$ 是 §3.2 的 BDF 质量项系数。）

$$
\underbrace{\int \kappa\,\nabla T^{n+1}\!\cdot\!\nabla\varphi_i}_{\text{扩散（}\kappa\text{取 max(物理, 人工)）}}
+\underbrace{\int (\rho c_p+L)\,\mathbf u^*\!\cdot\!\nabla T^{n+1}\,\varphi_i}_{\text{对流（速度取外推值 }\mathbf u^*\text{）}}
+\underbrace{\int (\rho c_p+L)\,\beta_\text{BDF}\, T^{n+1}\,\varphi_i}_{\text{时间导数（隐式 BDF）}}
=\text{RHS}_i
$$

代码（`advection.cc:216-227`）：

```cpp
data.local_matrix(i,j) += (
     (time_step * diffusion_constant * (scratch.grad_phi_field[i] * scratch.grad_phi_field[j]))
   + ((time_step * (scratch.phi_field[i] * (current_u * scratch.grad_phi_field[j])))
      + (bdf2_factor * scratch.phi_field[i] * scratch.phi_field[j])) * (density_c_P + latent_heat_LHS)
   ) * JxW;
```

注意 `diffusion_constant = std::max(conductivity, artificial_viscosity)`（`advection.cc:178-182`）——
**物理扩散与人工扩散取大者**，这是熵黏性的核心设计（见 §2.6）。

**DG 的面项**（`assemblers/advection.cc` 的 `AdvectionSystemInteriorFace`）：对 $DGQ_p$
在内部面上是**对称内罚 Galerkin（SIPG）+ 上风通量**，同时组装 4 个块
（int-int, int-ext, ext-int, ext-ext，`advection.cc:1086-1210`）：

$$
\big[-\tfrac12(\kappa\nabla T\cdot\mathbf n,\,[\![\varphi]\!])-\tfrac12(\kappa\nabla\varphi\cdot\mathbf n,\,[\![T]\!])+\underbrace{\gamma_D\,(\,[\![T]\!],[\![\varphi]\!])}_{\text{罚项}}-\underbrace{(\mathbf u\cdot\mathbf n\,T^{\text{up}},\,[\![\varphi]\!])}_{\text{上风}}\big]
$$

罚参数（`advection.cc:1010-1019`）：

```cpp
const double penalty = parameters.discontinuous_penalty * temperature_degree^2
                       / approximate_face_measure(face) * conductivity / (density_c_P + latent_heat_LHS);
```

即 $\gamma_D \sim \dfrac{\texttt{Discontinuous penalty}\cdot p^2\,\kappa}{h\,(\rho c_p+L)}$，
典型的 $p^2/h$ 内罚标度，`Discontinuous penalty` 默认 10。
对**成分方程** `penalty` 恒为 0（无扩散项），只有上风通量。
面上还要处理悬挂节点/粗细交界面（subface 循环，`advection.cc:1223-1247`）和周期边界。

### 2.5 自由度、约束与边界条件

`setup_dofs()`（`core.cc:1378-1488`）：

```cpp
dof_handler.distribute_dofs(finite_element);
DoFRenumbering::hierarchical (dof_handler);                    // 保证 restart 可复现
DoFRenumbering::component_wise (dof_handler, ...);             // 按块聚簇 → 分块矩阵
...
geometry_model->make_periodicity_constraints(dof_handler, constraints);   // 必须最先
DoFTools::make_hanging_node_constraints (dof_handler, constraints);       // 悬挂节点
compute_initial_velocity_boundary_constraints(constraints);               // 速度 Dirichlet/自由滑移
constraints.close();
```

- 所有线性约束统一放在 deal.II 的 `AffineConstraints<double>` 里。悬挂节点约束是自适应网格的必需品：
  粗单元边上的解必须由细单元插值决定，否则 $C^0$ 空间不协调。
- 非齐次 Dirichlet（温度、速度）每个时间步在 `compute_current_constraints()` 里用
  `VectorTools::interpolate_boundary_values` 重新生成（`core.cc:647-698`），
  因为边界值可能随时间/温度变化。
- **流入/流出边界**：若不用 Dirichlet 温度，ASPECT 会把 outflow 边界的 boundary id 临时偏移，
  以免被当作 Dirichlet（`core.cc:665-673`）。
- DG 的温度边界条件不走强约束，而是通过面上罚项弱施加（`AdvectionSystemBoundaryFace`，
  `advection.cc:645` 起）。

### 2.6 稳定化（这是"空间离散"里最关键的一环）

对流占优时 $Q_p$ 中心差分型 Galerkin 会产生振荡，ASPECT 提供两种机制：

**(A) 连续单元：entropy viscosity（熵黏性人工扩散）**，默认方法。
核心在 `source/simulator/entropy_viscosity.cc:118-238`：

- 每个单元一个常数 $\nu_h|_K$：

$$
\nu_h|_K=\min\left(\nu_h^{\max}|_K,\ \nu_h^{E}|_K\right)
$$

- 一阶最大耗散（保证单调性）：

```cpp
const double maximum_viscosity = parameters.stabilization_beta[...] *
                                 max_advection_prefactor * max_velocity * cell_diameter;
```

即 $\nu_h^{\max}=\beta\,(\rho c_p)\,\big(\|\mathbf u\|+\gamma h_K\|\varepsilon(\mathbf u)\|\big)h_K$
（$\gamma$ 项在 `entropy_viscosity.cc:193-195` 加入，用应变率增强剪切带的黏性）。

- 熵黏性（利用解的局部正则性）：

```cpp
if (parameters.stabilization_alpha == 2)
  entropy_viscosity = c_R * h^2 * max_residual / global_entropy_variation;
else
  entropy_viscosity = c_R * h * L_Omega * max_velocity * max_residual / (global_u_infty * global_field_variation);
```

  对应文档里的 $v_h^E|_K=\alpha_E\dfrac{h^2\|r_E\|_{\infty,K}}{\|E-E_\text{avg}\|_{\infty,\Omega}}$，
  熵 $E=\tfrac12(T-T_m)^2$，残差 $r_E=\partial_tE+(T-T_m)(\mathbf u\cdot\nabla T-k\triangle T-F)$
  （`compute_residual()`，`advection.cc:259-311`）。
  参数 `alpha`（默认 2）、`cR`（默认 0.11）、`beta`（默认 0.052）、`gamma`（默认 0.0）

- 前两个时间步、或熵变化为零时，直接退回到 $\nu_h^{\max}$（`entropy_viscosity.cc:213-218`）。
- **DG 不需要**熵黏性，函数第一行就 `return 0.`（`:126-128`）——因为上风通量本身提供了稳定化。
- 还有可选的 `Use artificial viscosity smoothing`（对 $\nu_h$ 做邻居平滑，避免指标噪声）。

**(B) DG 单元**：上风 + SIPG（见 §2.4），另可选**保界限幅器**
（`Use limiter for discontinuous ... solution`）。实现见
`helper_functions.cc:1332` 起的 `apply_limiter_to_dg_solutions()`：
在每个单元上用"一维 Gauss × 一维 Gauss–Lobatto"组合出的节点集求极值，
把解限制在 `Global temperature maximum/minimum`（或成分的 min/max）范围内，
是 Zhang–Shu 型 bound-preserving limiter 的思路。代码注释还说明它只对 Cartesian 映射做过测试
（`core.cc:462-468`）。

**(C) 可选 SUPG**（`Stabilization method = SUPG`，实验性）：
在权函数上加 $\tau\,\mathbf u\cdot\nabla\varphi_i$，组装里需要 Hessian（`advection.cc:96-97, 200-252`），
此时 `diffusion_constant = conductivity`（不再用人工黏性），$\tau$ 来自熵黏性量级。

---

## 3. 时间离散

### 3.1 主循环：Stokes 是"每层解一次的稳态问题"

时间主循环在 `source/simulator/core.cc:2009-2108`。一分钟版本：

```
while (true) {
  start_timestep();                 // 更新 BC、物性、信号
  solve_timestep();                 // ← 真正的时间层求解
  time_stepping_manager.update();   // 可能要求重复本步/网格加密
  maybe_refine_mesh(dt_new, ...);
  if (should_repeat_time_step()) { time -= time_step; time_step = dt_new; continue; }
  advance_time(dt_new);             // old_old ← old ← solution；t += dt；n++
  if (should_terminate) break;
}
```

`advance_time()`（`helper_functions.cc:688-709`）只做一件事：**滚动历史解向量**

```cpp
old_time_step = time_step;  time_step = step_size;  time += time_step;  ++timestep_number;
if (timestep_number == 1) { old_old_solution = solution; old_solution = solution; }
else                      { old_old_solution = old_solution;  old_solution = solution; }
```

所以 ASPECT 保存 **两个历史时间层**（`old_solution`、`old_old_solution`）+ 对应的时间步长，
这就是 BDF2 的全部"状态"。

### 3.2 温度和成分：隐式变步长 BDF1 / BDF2

时间导数在 `assemblers/advection.cc:71-75` 和 `:142-153` 里显式写出：

```cpp
const bool   use_bdf2_scheme = (this->get_timestep_number() > 1
                                && this->get_parameters().use_bdf2_for_advection_equations);
const double time_step     = this->get_timestep();
const double old_time_step = this->get_old_timestep();
const double bdf2_factor = (use_bdf2_scheme)
                           ? ((2*time_step + old_time_step) / (time_step + old_time_step)) : 1.0;
...
const double field_term_for_rhs
  = (use_bdf2_scheme
     ? ( scratch.old_field_values[q] * (1 + time_step/old_time_step)
         - scratch.old_old_field_values[q] * (time_step*time_step)
           / (old_time_step * (time_step + old_time_step)) )
     : scratch.old_field_values[q])
    * (density_c_P + latent_heat_LHS);
```

令 $r=\Delta t_n/\Delta t_{n-1}$，把上面两式对照**变步长 BDF2** 的标准形式：

$$
\frac{1}{\Delta t_n}\Big[\underbrace{\frac{1+2r}{1+r}}_{\texttt{bdf2\_factor}}T^{n+1}
-\underbrace{(1+r)}_{1+\Delta t_n/\Delta t_{n-1}}T^{n}
+\underbrace{\frac{r^2}{1+r}}_{\Delta t_n^2/[\Delta t_{n-1}(\Delta t_n+\Delta t_{n-1})]}T^{n-1}\Big]
$$

完全一致。当 $r=1$（等步长）退化为熟悉的 $\frac{3T^{n+1}-4T^n+T^{n-1}}{2\Delta t}$。

- **BDF1（后向 Euler）**：`use_bdf2_for_advection_equations=false`，或 `timestep_number <= 1` 时
  （前两步没有两个历史层）。此时 `bdf2_factor = 1`、RHS 用 $T^n$，即
  $(\rho c_p+L)\frac{T^{n+1}-T^n}{\Delta t}+\cdots$，一阶、无条件稳定、强耗散。
- **BDF2**（默认开启）是二阶、A-稳定（非 L-稳定，但对本问题足够），
  **变步长系数**保证了网格自适应/时间步变化时不会降到一阶。
- 注意时间步因子 `time_step` 乘在扩散与对流项上（因为整个方程乘了 $\Delta t$），
  而 BDF 质量项系数是 `bdf2_factor`。这是"乘 $\Delta t$"形式下最易看错的地方。

**成分方程用同一套系数**（`advection.cc:530-539` 的 Darcy 分支、
`source/simulator/melt.cc:687-850` 的熔融分支也都重复了同样的 BDF2 逻辑）。

### 3.3 速度与物性：外推（显式处理）——整体是 IMEX/半隐式

上一个时间层的解推进时，$\mathbf u$ 用哪个时间层？答案：**外推值**。

`initialize_current_linearization_point()`（`helper_functions.cc:1828-1854`）：

```cpp
current_linearization_point = old_solution;
if (timestep_number > 1 && parameters.use_extrapolated_current_linearization_point)
  {
    distr_solution.sadd((1 + time_step/old_time_step),
                        -time_step/old_time_step, distr_old_solution);
    current_linearization_point = distr_solution;
  }
```

即

$$
\mathbf u^{*}=(1+r)\,\mathbf u^{n}-r\,\mathbf u^{n-1},\qquad r=\frac{\Delta t_n}{\Delta t_{n-1}}
$$

（对等步长就是线性外推 $2\mathbf u^n-\mathbf u^{n-1}$，二阶精度）。

对流组装里读的正是它（`assembly.cc:865-866`）：

```cpp
scratch.finite_element_values[introspection.extractors.velocities]
     .get_function_values(current_linearization_point, scratch.current_velocity_values);
```

**材料模型/热源也在 `current_linearization_point` 上求值**（`assembly.cc:878-895`），
即 $\eta,\rho,c_p,k,H$ 全部显式（外推）处理。

所以整体时间格式是：

| 项 | 时间层 | 处理方式 |
|---|---|---|
| $\partial T/\partial t$ | $T^{n+1},T^n,T^{n-1}$ | **隐式 BDF1/BDF2** |
| $\mathbf u\cdot\nabla T$ | $T^{n+1}$ 隐式，$\mathbf u^*$ 外推 | **半隐式（IMEX）** |
| $\nabla\cdot k\nabla T$ | $T^{n+1}$，$k$ 外推 | 隐式（系数显式） |
| $\rho,\eta,g$（Stokes） | 外推/迭代点 | 显式或非线性迭代内更新 |

这样每个时间步只需解**线性**（或几个线性化）系统，无需处理 $\partial_t$ 带来的时间耦合。

**网格变形（ALE）**：如果开启自由表面，对流速度要减去网格速度
（`advection.cc:156-158`）：`current_u -= scratch.mesh_velocity_values[q]`，
即在对流项里解真是相对于运动网格的速度——这是 ALE 形式的时间离散修正。

### 3.4 时间步长控制：CFL + 一系列钳制

时间步长由一组**插件**给出，取所有插件的最小值（`source/time_stepping/interface.cc:88-156`）：

```cpp
for (const auto &plugin : active_plugins)
  new_time_step = std::min(new_time_step, plugin->execute());
new_time_step = Utilities::MPI::min(new_time_step, mpi_communicator);
new_time_step = std::max(minimum_time_step_size, new_time_step);
new_time_step = std::min(new_time_step, parameters.maximum_time_step);
if (get_timestep() != 0)
  new_time_step = std::min(new_time_step,
                           get_timestep() * (1 + maximum_relative_increase_time_step));
if (timestep_number == 0)
  new_time_step = std::min(new_time_step, maximum_first_time_step);
// 再由 termination 插件缩减，可能触发 repeat_step
```

两个内置插件：

**对流时间步**（`source/time_stepping/convection_time_step.cc:30-91`）：

$$
\Delta t_\text{conv}=\frac{\texttt{CFL}}{\;p_T\cdot\max_K\dfrac{\|\mathbf u\|_\infty}{h_K}\;}
$$

```cpp
min_convection_timestep = this->get_parameters().CFL_number
                          / (this->get_parameters().temperature_degree * max_global_speed_over_meshsize);
```

其中 $h_K$ = `cell->minimum_vertex_distance()`，分子里多一个多项式阶数 $p_T$
（高阶单元的有效网格尺度更小，所以步长要更小）。

**传导时间步**（可选，`Use conduction timestep=true`，`conduction_time_step.cc`）：

$$
\Delta t_\text{cond}=\texttt{CFL}\cdot\frac{h_K^2}{\kappa},\qquad \kappa=\frac{k}{\rho c_p}
$$

```cpp
min_local_conduction_timestep = std::min(min_local_conduction_timestep,
    CFL_number * pow(cell->minimum_vertex_distance(), 2.) / thermal_diffusivity);
```

**注意**：因为时间离散是**隐式**的，CFL 在这里主要不是稳定性限制，而是
**精度/精度换效率**的控制旋钮（`parameters.cc:130-144` 的说明明确写了
"对显式格式需要 $c\le1$，对隐式格式可以选 $c>1$"）。
默认 `CFL number = 1.0`。

此外 `Reaction` 机制允许插件要求"回到上一时刻、减小步长重算"
（`refine_and_repeat_step` / `repeat_step`，`time_stepping/interface.h`），
主循环在 `core.cc:2067-2085` 处理 repeat（并把粒子/网格位移一并回滚）。

### 3.5 时间层内的耦合：算子分裂与非线性迭代

`solve_timestep()`（`core.cc:1805`）先 `initialize_current_linearization_point()`，
按 `Nonlinear solver scheme` 分派到 11 种方案（`solver_schemes.cc`）。
默认是 **`single Advection, single Stokes`**（`parameters.cc:204`），也就是最便宜的
**一阶算子分裂/交错格式**（`solver_schemes.cc:815-837`）：

```
solve_single_advection_single_stokes():
    assemble_and_solve_temperature();
    assemble_and_solve_composition();
    fill_prescribed_fields();
    assemble_and_solve_stokes();
    fill_prescribed_fields();     // 本仓库的 TEQ 改动：Stokes 更新 p_f 后再刷一次
```

温度→成分→Stokes 各解一次，不在同一时间层内迭代耦合。其余方案（`iterated Advection and Stokes`、
`... defect correction`、`... Newton`）把上面的序列套进一个非线性迭代，
直到相对残差 < `Nonlinear solver tolerance`（默认 $10^{-5}$，`parameters.cc:259`）：

- **Picard / 定点迭代**：$x_{k+1}=A(x_k)^{-1}b$，每次重装矩阵。
- **Defect correction（缺陷修正）**：写成 Newton 形式但 Jacobian 用 $J\approx A(x)$，
  数学上等价于 Picard，但线性容差可以放松（1e-2 量级）。
- **Newton**：$J\delta x=-F$，$J=A(x)+c_kA'(x)$，含 SPD 稳定化、residual scaling、
  line search、fail-safe（`newton.cc`、`doc/sphinx/user/methods/nonlinear-stokes-solvers.md`）。

**算子分裂处理"反应"**（`Use operator splitting=true`）：
`compute_reactions()`（`helper_functions.cc:1604`）在求解前对温度/成分做
**显式子步进**（forward Euler），子步数

```cpp
number_of_reaction_steps = max(int(time_step / reaction_time_step),
                               max(reaction_steps_per_advection_step, 1));
reaction_time_step_size  = time_step / number_of_reaction_steps;
```

即"缓变项隐式大步 + 快变反应显式小步"的 IMEX 分裂。

---

## 4. 线性求解与预条件（与时空离散配套）

- **对流/扩散（温度、成分）**：分块 GMRES（重启长度 50）+ 块内 ILU(0) 预条件，
  相对容差默认 $10^{-12}$（`solver.cc:357-408`）。
- **Stokes**：分块鞍点系统
  $\begin{pmatrix}A&B^T\\B&0\end{pmatrix}$，
  用分块下三角 Schur 补预条件
  $\begin{pmatrix}A&0\\B&-S\end{pmatrix}^{-1}$（`solver.cc:171` 的 `BlockSchurPreconditioner`），
  外层 FGMRES。
  - $A$ 块：**代数多重网格**（`LinearAlgebra::PreconditionAMG`，底层 Trilinos ML；
    默认 Chebyshev 平滑 2 次，`AMG smoother type/sweeps`），
    并传入 `constant_modes`（速度零空间/常数模式）帮助 AMG 处理近零空间
    （`assembly.cc:409-434`）。可选 `use_full_A_block_preconditioner`。
  - $S$ 块：默认用**压力质量矩阵的 ILU**（`PreconditionILU`，`assembly.cc:419-422`）；
    只有开启熔融输运时才换成 AMG（因为那时压力块本身是椭圆算子）。
  - 先用较松容差做"cheap"步（`Number of cheap Stokes solver steps = 200`，
    A 块容差 1e-2、S 块 1e-6），失败再用"expensive"步（`solver.cc:750-809`）。
- **直接求解器**：`Use direct solver for Stokes system=true` 时用
  `TrilinosWrappers::SolverDirect`（Amesos），适合小规模。
- **矩阵无关几何多重网格**：`Stokes solver type = block GMG`
  （默认 `block AMG`），只支持速度阶数 2 或 3，用 `stokes_matrix_free.cc`
  的 $\mathcal O(n)$ 无矩阵算子（`core.cc:405-420`）。

---

## 5. ASPECT 与 deal.II 的分工

理解"哪些是 deal.II、哪些是 ASPECT 自己写的"对建立概念很重要：

| 功能 | 由谁提供 | 位置 |
|---|---|---|
| 网格（三角剖分）、自适应加密/粗化、悬挂节点、分布式 p4est 网格 | **deal.II** | `Triangulation`, `parallel::distributed::Triangulation`, `GridRefinement`, `SolutionTransfer` |
| 有限元基函数（$Q_k$、$DGQ_k$、$DGP_k$、`FESystem`）、映射 | **deal.II** | `FE_Q`, `FE_DGQ`, `FE_DGP`, `FESystem`, `MappingQ*` |
| 求积、`FEValues`/`FEFaceValues`（形函数值/梯度/散度/Hessian、JxW） | **deal.II** | `Quadrature`, `FEValues` |
| 自由度枚举、分块、`AffineConstraints`（悬挂节点/Dirichlet/周期） | **deal.II** | `DoFHandler`, `DoFTools`, `VectorTools`, `AffineConstraints` |
| 分布式线性代数（Trilinos）、Krylov 求解器、AMG/ILU 预条件 | **deal.II 封装 Trilinos** | `LinearAlgebra::BlockSparseMatrix`, `TrilinosWrappers::*`, `PreconditionAMG` |
| **方程组的弱形式、单元组装、边界/体源项** | **ASPECT** | `source/simulator/assemblers/*.cc`, `assembly.cc` |
| **时间离散：BDF1/BDF2、外推、历史解滚动** | **ASPECT** | `assemblers/advection.cc`, `helper_functions.cc:688,1828` |
| **时间步长控制（CFL、重复步、插件）** | **ASPECT** | `source/time_stepping/` |
| **熵黏性/SUPG/限幅器** | **ASPECT** | `entropy_viscosity.cc`, `helper_functions.cc:1332`, `assemblers/advection.cc` |
| **网格加密判据（物理量指标、合并策略）** | **ASPECT 插件** | `source/mesh_refinement/` |
| **非线性求解器（Picard/缺陷修正/Newton）** | **ASPECT** | `solver_schemes.cc`, `newton.cc` |

**结论**：deal.II 提供"离散化的零件"（网格、基函数、求积、自由度、线性代数），
**不提供 PDE 的时间积分器**。`deal.II` 的 `TimeStepping::ExplicitRungeKutta` /
`ImplicitRungeKutta` 是给"已知 $y'=f(t,y)$ 的 ODE"用的，ASPECT 没用
（`include/deal.II/base/time_stepping.h` 存在，但 ASPECT 全代码零引用）。
ASPECT 的时间离散完全是在组装循环里手写的——这也是为什么想理解时间离散，
必须去读 `assemblers/advection.cc` 的那几行系数，而不是找某个"时间积分器类"。

### 5.1 deal.II 侧的关键机制（与 ASPECT 的单元选择直接相关）

- **`FE_Q<p>`（连续张量积 Lagrange）**：$Q_p$ 空间是每个坐标方向 1 维 Lagrange
  多项式的张量积。支撑点（support points）分布：
  $p\le 2$ 用**等距点**（$Q_1$ 在 $\{0,1\}$，$Q_2$ 在 $\{0,0.5,1\}$）；
  $p\ge 3$ 默认用 $(p+1)$ 阶 **Gauss–Lobatto 点**。
  文档明确解释了原因：等距点在 $p$ 大时会出现 Runge 现象、质量矩阵条件数爆炸
  （$p=10$ 时约 $2.6\times10^6$，而 Gauss–Lobatto 约 400）。
  自由度编号按"顶点 → 棱 → 面 → 体"重排，所以 `DoFRenumbering::component_wise`
  之后矩阵天然分块。参见 `include/deal.II/fe/fe_q.h:34-95`。
- **`FE_DGQ<p>`（间断张量积）**：同样基于张量积 Lagrange，但**单元之间不共享自由度**，
  边界上的自由度属于该单元独有。它**没有悬挂节点约束**（因为本来就不连续），
  这就是 ASPECT 在 DG 模式下把稳定化交给面通量、并在面上做 subface 循环的原因。
  参见 `include/deal.II/fe/fe_dgq.h:42-97,248-282`。
- **`FE_DGP<p>`**：单元内**完全不连续的完备单项式**空间（$P_p$），
  自由度不在界面上，因此每个单元的自由度都能独立满足 $\int_{\partial K}\mathbf u_h\cdot\mathbf n=0$，
  这正是 ASPECT "locally conservative" 压力单元的理论来源。
- **悬挂节点约束**：`DoFTools::make_hanging_node_constraints` 在 `AffineConstraints`
  里为每个粗单元交界面上的细单元自由度写入线性关系（细边中点值 = 粗边端点值的平均等），
  使 $C^0$ 空间在自适应网格上仍协调。ASPECT 在 `setup_dofs()` 里调用它
  （`core.cc:1458`）。DG 单元不产生这类约束。
- **`MappingQCache`**：把 `MappingQGeneric(4)` 在每个单元上的变换数据缓存起来，
  避免每次 `FEValues::reinit` 重算高阶映射；代价是网格一变就要
  `map->initialize(MappingQGeneric<dim>(4), triangulation)` 重建
  （`core.cc:458-460, 1706-1707, 1968-1969`）。
- **`parallel::distributed::SolutionTransfer`**：自适应网格加密后把旧网格上的
  解向量插值到新网格。ASPECT 在 `refine_mesh()` 里对 `solution`、`old_solution`
  （以及网格变形相关的向量）注册，加密后 `interpolate()` 再 `constraints.distribute()`
  保证协调（`core.cc:1654-1750`）。注意它只做**插值/投影，不保证守恒**——
  粗化时质量守恒在 ASPECT 里靠后续的有限元 L2 投影与物理量重算实现。
- **`FEValues`/`FEFaceValues`**：把参考单元上的形函数值、梯度、散度、Hessian
  经映射 Jacobian 变换到物理单元，同时给出 `JxW`（Jacobian×权重）和法向。
  所有组装循环（`local_assemble_*`）都是这个模式；`update_flags` 控制只算需要的量
  （例如只在 SUPG 时才 `update_hessians`，`assembly.cc:1145-1149`），是性能关键。

---

## 6. 速查表

### 6.1 空间离散

| 项 | ASPECT 的选择 | 代码 |
|---|---|---|
| 网格 | 2D 四边形 / 3D 六面体，分布式 h 自适应 | `core.cc:455`, `core.cc:1592` |
| 映射 | `MappingQCache(4)`（曲边）/ `MappingCartesian`（平直）/ `MappingQ1Eulerian`（变形） | `core.cc:92-121` |
| 速度 | $Q_k^d$，默认 $k=2$ | `introspection.cc:192-196` |
| 压力 | $Q_{k-1}$（默认 $Q_1$）；或 $P_{-(k-1)}$；或 $Q_k$+稳定化 | `introspection.cc:198-212` |
| 温度 | $Q_2$ 或 $DGQ_2$ | `introspection.cc:215-221` |
| 成分 | $Q_2$/成分 或 $DGQ_2$/$DGQ_0$ | `introspection.cc:223-229` |
| 求积 | Gauss–Legendre，阶数 = 多项式阶数+1 | `introspection.cc:107-128` |
| Stokes 弱形式 | $(2\eta\varepsilon(\mathbf u),\varepsilon(\mathbf v))-(p,\nabla\!\cdot\!\mathbf v)-(q,\nabla\!\cdot\!\mathbf u)=(\rho\mathbf g,\mathbf v)$ | `stokes.cc:409-421` |
| 对流弱形式 | $\Delta t(\kappa\nabla T,\nabla\varphi)+\Delta t(\rho c_p\,\mathbf u^*\!\cdot\!\nabla T,\varphi)+\text{BDF}(T,\varphi)$ | `advection.cc:216-227` |
| DG 稳定化 | SIPG 罚 $\gamma_D\sim \texttt{penalty}\,p^2\kappa/(h\rho c_p)$ + 上风 | `advection.cc:1010-1208` |
| CG 稳定化 | entropy viscosity，$\nu=\min(\nu^{\max},\nu^E)$，逐单元常数 | `entropy_viscosity.cc:118-238` |
| 约束 | `AffineConstraints`：周期→悬挂节点→Dirichlet | `core.cc:1447-1463` |

### 6.2 时间离散

| 项 | ASPECT 的选择 | 代码 |
|---|---|---|
| Stokes | 无时间导数，每层解一次椭圆鞍点问题 | `assemblers/stokes.cc` |
| 温度/成分时间导数 | 隐式变步长 BDF2（默认）/ BDF1 | `advection.cc:71-75,142-153` |
| BDF2 系数 | 质量项 $(1+2r)/(1+r)$，RHS $(1+r)T^n-\frac{r^2}{1+r}T^{n-1}$，$r=\Delta t_n/\Delta t_{n-1}$ | 同上 |
| 对流速度 | 外推 $\mathbf u^*=(1+r)\mathbf u^n-r\mathbf u^{n-1}$ | `helper_functions.cc:1828-1854`, `assembly.cc:865` |
| 物性/热源 | 在同一外推点求值（显式） | `assembly.cc:878-895` |
| 时间步 | $\min$ of 插件；CFL 控制 | `time_stepping/interface.cc:88-156` |
| 对流 CFL | $\Delta t=\texttt{CFL}/(p_T\max\|\mathbf u\|/h)$ | `convection_time_step.cc:83` |
| 传导 CFL（可选） | $\Delta t=\texttt{CFL}\,h^2/\kappa$ | `conduction_time_step.cc:99-103` |
| 时间层耦合 | 默认 `single Advection, single Stokes`（一阶交错） | `solver_schemes.cc:815` |
| 非线性 | Picard / 缺陷修正 / Newton，容差 $10^{-5}$ | `solver_schemes.cc`, `newton.cc` |
| 反应 | 算子分裂 + 显式子步 | `helper_functions.cc:1604-1625` |

---

## 7. 常见误解与注意点

1. **"ASPECT 用 deal.II 的时间积分器"——错。** deal.II 只提供 ODE 积分器且 ASPECT 未使用；
   时间离散是手写的 BDF。
2. **"CFL 是稳定性限制"——在本格式里主要不是。** 温度/成分是隐式离散，
   CFL 控制精度与时间步增长率；但外推的速度项和显式反应子步仍会带来精度/Courant 约束。
3. **"网格加密用 Kelly 误差估计子"——默认不是。**
   默认 `Strategy = thermal energy density`（`mesh_refinement/interface.cc:344`），
   是物理量指标（温度、成分、黏度、应变率等），可由多个指标归一化后加权合并；
   加密用 `refine_and_coarsen_fixed_fraction`（默认，按误差占比）或 `fixed_number`（按单元比例）。
   频率由 `Time steps between mesh refinement`（默认 10）控制。
4. **"BDF2 一开始就用"——不是。** 第 0/1 步以及重启后的头两步只能用 BDF1；
   而且 BDF2 可以用 `Use BDF2 for advection equations=false` 关闭。
5. **变步长 BDF2 的系数不是常数 3/2、2、1/2**，而是上面含 $r$ 的表达式；
   网格自适应会让 $\Delta t$ 变化，若用常系数会掉到一阶甚至不稳定。
6. **DG 与 CG 的稳定化机制不同**：DG 靠上风+SIPG（无需熵黏性），CG 靠熵黏性；
   混用时（温度 DG、成分 CG）两套机制会分别生效。
7. **压力标度 `pressure_scaling` 不是物理量**，是数值条件数处理；
   它出现在 Stokes 的 $(U,P)$、$(P,U)$ 块和等阶稳定项中，改等阶选项时要注意。
8. **本仓库的 TEQ 改动**：`solve_single_advection_single_stokes()` 里 Stokes 之后再调用一次
   `fill_prescribed_fields()`（`solver_schemes.cc:822-827`），
   以保证孔隙度用当前步的 $p_f$ 求值——这是纯粹的一致性修正，不改变上述离散格式。

---

## 8. 熔体输运（melt transport）：额外值得补充的方程与离散技巧

> 代码：`include/aspect/melt.h`（552 行）、`source/simulator/melt.cc`（2041 行）、
> `source/simulator/assemblers/advection.cc` 的 `DarcySystem`、
> `source/mesh_refinement/compaction_length.cc`。
> 文档：`doc/sphinx/user/methods/melt-transport.md`。
> 理论出处：Keller, May & Kaus (2013)；Dannberg & Heister (2016)。

前面 §1–§7 讲的是**单相**（固体地幔）问题。开启 `Include melt transport` 后，
问题从"Stokes + 对流扩散"变成**两相流**，数值方法上有 5 类实质性的新东西。

### 8.1 多出来的未知量与方程

`MeltHandler::edit_finite_element_variables()`（`melt.cc:1689-1722`）往有限元系统里插入 3 个变量：

| 新变量 | 单元 | 说明 |
|---|---|---|
| `fluid pressure` $p_f$ | $Q_{k-1}$（默认 $Q_1$） | 熔体（流体）压力 |
| `compaction pressure` $p_c$ | $Q_{k-1}$ **或** `FE_DGP(k-1)`（默认**间断**，`Use discontinuous compaction pressure=true`） | 压实压力，$p_c=(1-\phi)(p_s-p_f)$ |
| `fluid velocity` $\mathbf u_f$ | $Q_k^d$（默认 $Q_2^d$） | 熔体速度，**诊断量**（见 §8.3） |

此外必须有一个名叫 `porosity` 的成分场 $\phi$（否则报错，`parameters.cc:1858-1865`）。

方程被替换/扩充为（`melt-transport.md:18-174`，代码在 `MeltStokesSystem::execute`，`melt.cc:396-566`）：

1. **动量方程**（多了两个压力梯度，浮力用体积密度 $\bar\rho=(1-\phi)\rho_s+\phi\rho_f$）：
   $-\nabla\cdot[2\eta(\varepsilon(\mathbf u_s)-\tfrac13(\nabla\cdot\mathbf u_s)\mathbf 1)]+\nabla p_f+\nabla p_c=\bar\rho\,\mathbf g$
2. **质量守恒被 Darcy 型方程替换**（不再是 $\nabla\cdot\mathbf u=0$）：
   $\nabla\cdot\mathbf u_s-\nabla\cdot(K_D\nabla p_f)-K_D\nabla p_f\cdot\frac{\nabla\rho_f}{\rho_f}
   =-\nabla\cdot(K_D\rho_f\mathbf g)+\Gamma(\frac1{\rho_f}-\frac1{\rho_s})-\cdots$
   （对应 `compute_fluid_pressure_rhs()`，`melt.cc:329-385`）
3. **压实方程**：$\nabla\cdot\mathbf u_s+\dfrac{p_c}{\xi}=0$，$\xi$ 是压实黏度（`MeltOutputs::compaction_viscosities`）
4. **Darcy 定律（诊断）**：$\mathbf u_f=\mathbf u_s-\dfrac{K_D}{\phi}(\nabla p_f-\rho_f\mathbf g)$
5. **孔隙度演化**（置换普通成分方程）：
   $\dfrac{\partial\phi}{\partial t}+\mathbf u_s\cdot\nabla\phi-\phi\,(\nabla\cdot\mathbf u_s+\beta_s\rho_s\mathbf g\cdot\mathbf u_s)=\dfrac{\Gamma}{\rho_s}$

### 8.2 空间离散上的新东西

**(a) 每单元的压力标度 $p_c^{\text{scale}}$。**
$p_c$ 与 $p_f$ 量级可能差很多，所以代码逐单元计算
$p_c^{\text{scale}}=\sqrt{\overline{K_D}/K_D^{\text{ref}}}$（`melt.cc:109-150`），
在组装时把 $p_c$ 乘/除这个标度（`melt.cc:551-554`）。这与 §2.4 的 `pressure_scaling` 是两回事，
是熔体系统特有的**局部无量纲化**。

**(b) "熔体单元 / 非熔体单元"的离散切换。**
如果某单元的 Darcy 系数 $K_D$ 低于
`Melt scaling factor threshold`（默认 $10^{-7}$）×参考值，就认为该单元没有熔体：
**把 $p_c$ 自由度约束为 0、丢掉所有熔体项**，局部退回普通 Stokes
（`melt.cc:1545-1614` 的 `add_current_constraints()` + `PcNonZeroDofsAssembler`；
在 `core.cc:743-744` 每个时间步调用）。
注意这是**依赖当前解 $\phi$ 的约束集**：熔体区边界会随时间移动，
所以约束（以及可能的分块稀疏模式）每个时间步都可能变。

**(c) 混合物单元与块结构。**
熔体分支在 `setup_system_matrix_coupling()` 里显式指定耦合（`core.cc:921-956`）：
速度-速度、$p_f$-$p_f$、$p_c$-$p_c$、速度-$p_f$、速度-$p_c$ 以及它们的转置。
同时熔体组装器会**替换** Stokes 的自由度列表
（`assembly.cc:509-515` 的注释："assemblers below can modify this list of dofs, if they in fact
assemble a different system than the standard Stokes system (e.g. in models with melt transport)"；
实际过滤在 `melt.cc:159-177` 的 `is_velocity_or_pressures`）。

**(d) 硬性限制**（`MeltHandler::initialize()`，`melt.cc:1875-1891`）：
- **DG 温度/成分与熔体不兼容**（"The additional terms in the temperature systems have not been ported to the DG formulation"）——所以 §2.6 的 DG + 上风稳定化在熔体模型里用不了，只能用 CG + 熵黏性；
- **不支持 DG 的 $p_f$**（不能同时用 `Use locally conservative discretization`）；
- 质量守恒只能是 `incompressible` 或 `isentropic_compression`（`melt.cc:1734-1739`）。

**(e) 自由表面稳定化**有熔体专用版本 `apply_free_surface_stabilization_with_melt()`（`melt.cc:1920`）。

### 8.3 时间离散上的新东西

**(a) 孔隙度方程不是标准对流方程。**
`MeltAdvectionSystem::execute()`（`melt.cc:676-938`）在标准对流骨架上加了：
- **LHS 源项** `melt_transport_LHS = ∇·u_s + β_s ρ_s g·u_s`（`melt.cc:823-837`），
  对应 §8.1 第 5 式中从 $(1-\phi)\nabla\cdot\mathbf u_s$ 展开后移过来的 $-\phi(\cdots)$ 项；
- **RHS 源项** `compute_melting_RHS()`：$\Gamma/\rho_s+\nabla\cdot\mathbf u_s+\beta_s\rho_s\mathbf g\cdot\mathbf u_s$（`melt.cc:640-671`）；
- **$1-\phi$ 的展开**：代码注释写得很清楚
  "the full advection term would be $(1-\phi)(\nabla\cdot\mathbf u+\kappa\rho\mathbf u\cdot\mathbf g)$,
  but we expanded the term and move the part $(-\phi)(\cdots)$ over to the LHS"（`melt.cc:662-663`）。
- **BDF1/BDF2 系数与普通成分完全相同**（`melt.cc:792-851`），所以 §3.2 的时间离散结论照旧适用。

**(b) $\nabla\cdot\mathbf u_s$ 取单元平均值。**
`melt.cc:716-719`：

```cpp
double divergence_u = 0.0;
for (unsigned int q=0; q<n_q_points; ++q)
  divergence_u += scratch.current_velocity_divergences[q] * 1./n_q_points;
```

即把散度在单元内取**常数**（单元平均）再用于孔隙度方程。这是一个稳定性处理
（避免 $\nabla\cdot\mathbf u_s$ 的逐点噪声进入孔隙度源项）。注意在熵黏性残差里
（`compute_residual`，`melt.cc:1008-1021`）用的是**逐点**散度，两者有意不同。

**(c) 熔体速度 $\mathbf u_f$ 不在隐式系统里，而是每步投影出来的。**
`compute_melt_variables()`（`melt.cc:1138-1329`，在 `solver.cc:905` 每次 Stokes 求解后调用）：
1. 用 **L2 投影（质量矩阵问题）**把 $\mathbf u_f=\mathbf u_s-\frac{K_D}{\phi}(\nabla p_f-\rho_f\mathbf g)$
   投到 $Q_k^d$ 空间，用 **CG + AMG** 解（`melt.cc:1304-1323`）；
2. 可选 `Average melt velocity=true`（默认）：$K_D/\phi$ 用**单元几何平均**
   （`melt.cc:1239-1248`），避免熔体/无熔体界面上 $K_D,\phi\to0$ 造成的速度尖峰；
3. 非熔体单元里直接令 $\mathbf u_f=\mathbf u_s$（`melt.cc:1288-1292`）。

**(d) 固体压力 $p_s$ 是代数重构的**（`melt.cc:1331-1400`）：
在压力支撑点上逐点算 $p_s=p_f+\dfrac{p_c}{1-\phi}$（代码里含 $p_c^{\text{scale}}$），
不需要解方程。

**(e) 时间步长必须考虑熔体速度。**
`convection_time_step.cc:60-68` 在熔体模型里把 $\|\mathbf u_f\|$ 也纳入 $\max$：

```cpp
if (this->get_parameters().include_melt_transport)
  { ... max_local_velocity = std::max(max_local_velocity, fluid_velocity_values[q].norm()); }
```

由于熔体渗透速度通常远大于固体速度，**熔体模型的时间步长常常由 $\mathbf u_f$ 的 CFL 决定**，
而不是 $\mathbf u_s$。这正是文档里"melt velocity … is only used for postprocessing purposes
and for computing the time step length"的含义。

**(f) 温度/成分的对流速度可以选择。**
- `Heat advection by melt=true`：温度对流用**孔隙度加权的固/液速度组合**
  （`density_c_P_solid=(1-\phi)\rho_s c_p`、`density_c_P_melt=\phi\rho_f c_p`，
  `melt.cc:856-864`）——温度方程 LHS 多出一整套熔体对流项；
- 成分场的 `fem_melt_field` 方法：用 $\mathbf u_f$ 而非 $\mathbf u_s$ 平流该成分（`melt.cc:865-873`）；
- **bulk concentration advection**：把化学组分按"体浓度 + 相内浓度"分解，
  在普通成分方程上额外加 $\nabla\cdot\mathbf u$ 与相分离通量的项
  （`melt.cc:810-846`，受编译宏 `ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS` 控制；
  本仓库的 `melt.h:31-36` 显示该宏在编译时会被打印出来）。

### 8.4 非线性耦合：熔体模型几乎必须迭代

$\eta(\phi)$、$\xi(\phi)$、$K_D(\phi)$、$\Gamma(T,p,\phi)$ 都是强非线性且互相耦合，
所以常规做法是 §3.5 的 `iterated Advection and Stokes` / `defect correction` / `Newton`，
而不是默认的 `single Advection, single Stokes`。相关强制设定：

- 熔体强制打开完整 A 块预条件：`parameters.use_full_A_block_preconditioner = true`（`core.cc:389-393`）；
- $p_f$ 块本身变成椭圆算子，所以 $M_p$ 预条件从 ILU 换成 **AMG**（`assembly.cc:416-422`）；
- `p_c` 的约束依赖解 → 每次非线性迭代/时间步都要重新生成约束。

### 8.5 自适应网格：需要专用的加密判据

§2.2 的默认物理量指标对熔体问题不够，于是有专门的
`compaction length` 插件（`source/mesh_refinement/compaction_length.cc:71-90`）：
用压实长度

$$
\delta_c=\sqrt{\dfrac{(\eta+\tfrac43\xi)\,k}{\eta_f}}
$$

与单元尺寸 $h_K$ 比较，$\delta_c<\texttt{cells\_per\_compaction\_length}\cdot h_K$ 就加密，
$\delta_c<2\cdot(\cdots)$ 就禁止粗化，并且只在 `is_melt_cell` 的单元上操作。
`Mesh cells per compaction length` 默认 1.0。这是典型的"离散尺度必须解析物理边界层"的例子。

### 8.6 一句话总结熔体带来的补充

| 维度 | 单相问题（§1–§7） | 加上熔体后 |
|---|---|---|
| 方程 | Stokes + 对流扩散 | +Darcy 质量守恒、+压实方程、+孔隙度演化 |
| 未知量 | $u,p,T,c_i$ | +$p_f,p_c,u_f,\phi$（其中 $u_f$ 是投影诊断量） |
| 空间离散 | $Q_2^d\times Q_1$ + CG/DG 可选 | +$Q_1$/$DGP_1$ 压力、+每单元 $p_c$ 标度、+熔体单元切换（解相关约束） |
| 稳定化 | 熵黏性 / DG 上风 | 熔体下**禁用 DG**，只能熵黏性；散度取单元平均 |
| 时间离散 | 隐式 BDF1/BDF2，速度外推 | 系数不变；但 $u_f$ 通过质量矩阵投影更新，$\Gamma$ 作源项，Δt 受 $u_f$ 的 CFL 支配 |
| 非线性 | 可选迭代 | 几乎必须迭代；A 块预条件全开、$p_f$ 块用 AMG |
| 网格 | 物理量指标 | +压实长度指标 |

**最值得记住的三点**：
1. 熔体把"椭圆 Stokes + 双曲对流"变成了"**三个耦合的压力/散度方程 + 一个带源项的对流方程**"，$p_f$ 与 $p_c$ 是两个独立的压力未知量；
2. **$\mathbf u_f$ 不是隐式求解的未知量，而是每步用 L2 投影从 Darcy 定律重构的诊断量**——所以它只影响时间步长和对流速度选择；
3. **熔体与 DG 不兼容**（温度/成分），这直接限定了稳定化方案只能用熵黏性。

---

## 附：一页纸的"概念图"

```
空间：      连续体               →  网格 (quad/hex)        →  有限元空间          →  求积/组装
            Ω ⊂ R^d                 K_h, h 自适应            Q2^d×Q1 / Q2 / DGQ2     Gauss, Σ_K ∫_K (...)

时间：      t ∈ [0,T]            →  时间层 t^n, Δt^n        →  历史解 (n, n-1)     →  隐式 BDF + 外推
                                    CFL 控制 Δt             old/old_old solution     u 显式外推, T 隐式

每一步：    ① 外推 u* = (1+r)u^n − r u^{n-1}，物性在 u* 处求值
            ② 解温度：BDF2 隐式 + 熵黏性(或 DG 上风)
            ③ 解成分：同上
            ④ 解 Stokes：2ηε(u):ε(v) − p∇·v − q∇·u = ρg·v   (无时间导数！)
            ⑤ 自适应网格 / 时间步控制 → advance_time()
```

**空间离散 = "怎么把一个 PDE 变成一个大稀疏线性方程组"；
时间离散 = "怎么把 ∂/∂t 变成相邻两个（三个）时间层之间的代数关系"。**
ASPECT 把前者几乎全部交给 deal.II，把后者完全握在自己手里。

```
若开启熔体输运（§8），每步变成：
            ① 外推 u_s*，物性在 u_s* 处求值
            ② 解孔隙度 φ：BDF2 隐式 + 单元平均 ∇·u_s 源项 + 熵黏性
            ③ 解其他成分；温度可用 φ 加权的固+液速度对流
            ④ 解"熔体 Stokes"：动量(∇p_f+∇p_c) + Darcy 质量守恒 + 压实方程
               （非熔体单元把 p_c 约束为 0，退回普通 Stokes）
            ⑤ 投影重构 u_f（质量矩阵 + CG/AMG）与代数重构 p_s
            ⑥ 自适应网格（含压实长度判据）/ 时间步控制（Δt 受 u_f 的 CFL 支配）
```
