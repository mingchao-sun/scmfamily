# scmfamily — 合成控制（jamovi 模块）

十种合成控制估计量及比较分析，用于比较案例研究和面板数据。 英文完整说明见 [README.md](README.md)。

## 分析（jamovi 菜单：SCM）

| 分析 | 内容 |
|---|---|
| **合成控制** | Abadie、Diamond、Hainmueller (2010)，安慰剂推断 |
| **惩罚合成控制** | Abadie 与 L'Hour (2021) |
| **去均值合成控制** | Ferman 与 Pinto (2021) |
| **增强合成控制** | Ben-Michael 等 (2021)，`augsynth` 的移植 |
| **弹性网合成控制** | Doudchenko 与 Imbens (2016) |
| **贝叶斯合成控制** | Kim、Lee、Gupta (2020) |
| **交互固定效应** | Xu (2017) |
| **矩阵补全** | Athey 等 (2021) |
| **合成双重差分** | Arkhangelsky 等 (2021)，`synthdid` 的移植 |
| **错时增强合成控制** | `multisynth` 的移植 |
| **比较估计量** | 最多九种方法并列 |

每个分析的输出部分都有"方法说明"和"数据与设置指南"。

SDID、增强合成控制和错时增强合成控制与官方包在检验的设定下完全一致（差为 0）；其余方法与 `augsynth`、`gsynth`、`MCPanel`、`pensynth` 等对照（`tests/test_methods.R`）。要求平衡面板。示例数据：加州 Proposition 99、西德统一、最低工资与青少年就业。

## 安装
- **jamovi 模块库**：在 jamovi 右上角打开模块按钮，进入 *jamovi library*，搜索并安装（模块上架后）。
- **侧载**：下载 `.jmo` 文件，在同一菜单选 *Install from file*。`.jmo` 只适用于构建它的操作系统、处理器类型和 jamovi 版本。
- **从源码**：需要 R 和 `jmvtools`：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（路径按本机 jamovi 修改）。

## 引用与许可
引用：Sun, M. *scmfamily* (jamovi module). https://github.com/mingchao-sun/scmfamily 。请同时引用所用方法的原始文献（输出中列出）。

许可：见 `DESCRIPTION`（GPL）；第三方代码与数据见 `NOTICE`。问题与建议：https://github.com/mingchao-sun/scmfamily/issues
