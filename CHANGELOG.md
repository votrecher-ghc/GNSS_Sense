# StarDial 变更记录

本文档记录 StarDial 代码与论文的同步修改。论文侧对应记录位于
`D:\论文\StarDial-src\CHANGELOG.md`。

## 2026-08-12 12:27:00（北京时间，UTC+08:00）

### Fig. 16 紧凑布局

- 将 Attack-defense rates 图由约 `1056×836 pt` 的偏高画布调整为 `1056×575 pt` 的宽扁画布，论文中的显示宽度保持不变，纵向占用减少约 31%。
- 图例改为绘图区上方的单行双项布局，横轴类别标签整体下移，避免标签进入柱状绘图区。
- 柱高、类别、坐标范围和配色保持不变；新增 `gesture_analysis/scientific_graphing/compact_attack_defense_rates.py` 以便复现紧凑矢量 PDF，并同步修改 MATLAB 源绘图布局。

## 2026-08-12 11:11:00（北京时间，UTC+08:00）

### 实验结果图字体与 Fig. 17 布局协调

- 调整 `standardize_environment_figure_set.py` 与 `square_environment_metrics.py`：Fig. 17 四个子图统一轴标题、刻度、图例和色条的视觉层级，两个混淆矩阵的横轴类别标签下移并改为向右下方旋转，避免标签进入热力图单元格。
- 新增 `gesture_analysis/scientific_graphing/harmonize_evaluation_pdf_fonts.py`，根据论文中的实际插入宽度逐图调整实验结果 PDF 的文字比例，使半栏图与全栏图在 Overleaf 中接近 Fig. 11(b) 的视觉字号，而不是机械使用同一个源字号。
- 时空认证通过率图的短轴标题与刻度适当放大；较长的纵轴标题单独采用可完整显示的字号，避免 `Authentication pass rate (%)` 在 PDF 边界处被截断。
- 论文侧同步更新 8 张实验结果 PDF，并将 open-field EER 正文数值由 1.96% 修正为与 ROC 实验图一致的 2.79%。

## 2026-08-11 18:47:00（北京时间，UTC+08:00）

### Fig. 17 四图统一

- 新增 `gesture_analysis/scientific_graphing/standardize_environment_figure_set.py`，将两个混淆矩阵、认证指标柱状图和 ROC 图统一为 `800×800 pt` 的正方形矢量画布。
- 四图统一采用轴标题 50 pt、刻度 27 pt、图例/色条 23 pt 的字号层级，并协调绘图区边距与占比。
- ROC 图从原矢量 PDF 中提取并重绘全部三条曲线，保持曲线点、EER 标注和实验数据不变。
- `export_paper_figures_data_driven.m` 中同步加入方形柱状图布局和统一字号设置，便于后续重新导出时保持一致。
- 论文侧 Fig. 17 的四个子图容器统一为 `0.47\columnwidth`，并全部按 `width=\linewidth` 插入。

## 2026-08-11 18:01:00（北京时间，UTC+08:00）

### Fig. 17 子图尺寸修正

- 修正此前仅统一 LaTeX 高度但可见绘图区仍不一致的问题。
- 将 Fig. 17(c) 的认证指标柱状图及其导出逻辑改为方形画布和方形绘图区，使其与 Fig. 17(d) 的 ROC 图具有一致的可见尺寸。
- 新增 `gesture_analysis/scientific_graphing/square_environment_metrics.py`，用于在无 MATLAB 环境中复现相同数据和配色的方形柱状图。

## 2026-08-11 17:36:18（北京时间，UTC+08:00）

### 论文图表排版

- 将 Fig. 17(c) 的认证指标柱状图和 Fig. 17(d) 的 ROC 曲线图改为相同固定高度，并保持各自纵横比，消除源 PDF 画布比例不同造成的可见尺寸不一致问题。
- 论文侧修改位置：`D:\论文\StarDial-src\Evaluation.tex` 中 Fig. 17 的第二行子图。

## 2026-08-11 16:27:09（北京时间，UTC+08:00）

### 代码与实验

- 新建独立的 StarDial 在线认证流程，仅执行认证所需的信号预处理、衍射特征提取、几何反演、手势模板匹配、物理一致性验证和联合判决。
- 增加运行开销测试，覆盖认证时延、CPU 时间、内存、吞吐量、冷启动和输入规模变化。
- 当前测得平均认证时延为 100.26 ms，P95 为 104.39 ms，P99 为 106.75 ms；平均 CPU 时间为 97.03 ms，峰值内存为 109.52 MiB，吞吐量为 9.97 次/秒。
- Fig. 10 的原始绘图逻辑定位于 `gesture_analysis/scientific_graphing/export_paper_figures_data_driven.m`。
- 将轨迹重建画廊布局由 4×3 调整为 6×2，并采用共享坐标标签、紧凑无边框图例和统一字体规范，以减少论文中的纵向占用。
- 新增论文绘图规范 `gesture_analysis/scientific_graphing/PAPER_FIGURE_STYLE.md`、紧凑布局脚本 `gesture_analysis/scientific_graphing/reflow_trajectory_gallery.py` 和运行开销测试目录 `runtime_benchmark/`。
- 上述代码修改对应本地 Git 提交：`015b4d3 feat: benchmark authentication and compact trajectory gallery`。

### 论文

- 在 `Evaluation.tex` 末尾新增 `Computational Overhead` 小节及实验结果表格。
- 实验设置与论文已有表述统一为 Intel Core i5、16 GB RAM、25 Hz 采样率和 12 类手势。
- 删除 Python/NumPy、精简实现、合成观测等实现限定，以及论文未统一声明的 4 秒窗口、16 颗卫星和 160 个采样点设置。
- Fig. 10 改用 `figures/eval/traj/traj_gallery_data_driven_compact.pdf`，LaTeX 引用位于 `Evaluation.tex:117`。
- 新图尺寸由约 961×802 pt 调整为约 960×410 pt。

### 目录说明

- 代码仓库：`D:\论文\code\SatLock`
- 论文目录：`D:\论文\StarDial-src`（当前未初始化为 Git 仓库）
