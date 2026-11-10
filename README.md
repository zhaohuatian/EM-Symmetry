# EM-Symmetry

用于对称结构电磁过程的半解析光学计算 MATLAB 仓库，面向光散射、光吸收和电偶极子辐射。利用平移对称性或球对称性，结合平面波分解、传输矩阵、散射矩阵及球面波多极展开，逐步积累简洁、可复用的光学计算程序。

## 当前功能

| 模块 | 功能 | 主要接口 |
|---|---|---|
| `01_PlaneWaveMultilayer` | 平面多层结构 TE/TM 反射、透射、吸收与内部场；TMM/SMM | `mlpw.solve`, `mlpw.field` |
| `02_DipoleMultilayer` | 多层介质中电偶极子的上下远场、角分布、逃逸辐射功率与局域总衰减 | `mldipole.farfield`, `mldipole.radiative`, `mldipole.decay` |
| `03_DipoleSphere` | 均匀球外电偶极子的 Mie 散射、相干远场、辐射及吸收 | `mlsphere.mie`, `mlsphere.farfield`, `mlsphere.decay` |

第二模块依赖第一模块，第三模块独立。三个程序均为 CPU MATLAB 实现，不要求 GPU 或额外工具箱。开发环境为 MATLAB R2025a，更早版本兼容性尚未确认。

未来将逐步补充周期结构 RCWA、标量衍射及其他半解析散射方法；这些属于规划，当前尚未实现。

## 使用

在仓库根目录启动 MATLAB：

```matlab
setup;
run_examples;
% 或分别运行 demo_planewave、demo_layered_dipole、demo_sphere
```

薄膜反射、透射和吸收示例：

```matlab
setup;
C=struct('lambda',600e-9,'epsilon',[1 2.25 1], ...
 'thickness',100e-9,'kparallel',0,'polarization','TE');
S=mlpw.solve(C,'smm');
disp([S.reflectance S.transmittance S.absorptance]);
F=mlpw.field(S,linspace(-200e-9,300e-9,1000));
```

各模块 `docs/API.md` 和 `docs/THEORY.md` 提供参数、输出与理论说明，函数内帮助提供调用细节。基本示例直接显示图形，不保存生成图片或数据。

## 物理约定与边界

长度为 m，输入相对介电常数，时间约定为 exp(-i omega t)。电偶极矩为复峰值相量，单位 C m；绝对电场为 V/m，角功率密度为 W/sr，增强因子默认相对真空。

平面偶极子常规路线支持各向同性非磁性层，要求源层透明、正介电常数，源严格位于层内。另提供对齐对角各向异性张量的互易远场路线，不包含各向异性局域总衰减。平面结构的总衰减与逃逸辐射之差可以包含导模和吸收，不能一概解释为吸收。

球模块目前支持透明宿主中的均匀、各向同性、非磁性球及球外电偶极子，不支持同心多层球或球内源。角积分网格、局域反馈积分及球面波截断需针对具体结构检查收敛，尤其是小间距、共振及掠射方向。

## 轻量源码组织

`modules/` 保留三个核心函数包、基本使用示例与必要理论文档；`setup.m` 添加路径，`run_examples.m` 运行示例。不收录验证脚本、验证数据、历史版本、大型图片及开发辅助记录。计算结果建议保存在已忽略的 `results/` 或独立数据目录。

## 作者与来源

Zhaohua Tian。相关原始项目：[TransferMatrix](https://github.com/knifelees3/TransferMatrix)、[DipoleFarfieldinStratifiedMediaMATLAB](https://github.com/knifelees3/DipoleFarfieldinStratifiedMediaMATLAB)、[DipoleEmissionNearNanoSphere](https://github.com/knifelees3/DipoleEmissionNearNanoSphere)。本仓库整理并扩展这些方向的计算接口。

采用 [MIT License](LICENSE)。
