# StarDial 变更记录

本文档记录 StarDial 代码与论文的同步修改。论文侧对应记录位于
`D:\论文\StarDial-src\CHANGELOG.md`。

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
