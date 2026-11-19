# scmfamily — 合成コントロール（jamovi モジュール）

比較事例研究とパネルデータのための 10 種類の合成コントロール推定量と比較分析。 詳しい英語の説明は [README.md](README.md)、中国語の要約は [README.zh.md](README.zh.md) を参照してください。

## 分析（jamovi メニュー：SCM）

| 分析 | 内容 |
|---|---|
| **合成コントロール** | Abadie、Diamond、Hainmueller (2010)。プラセボ推測 |
| **ペナルティ付き合成コントロール** | Abadie と L'Hour (2021) |
| **デミーン合成コントロール** | Ferman と Pinto (2021) |
| **拡張合成コントロール** | Ben-Michael ら (2021)。`augsynth` の移植 |
| **弾性ネット合成コントロール** | Doudchenko と Imbens (2016) |
| **ベイズ合成コントロール** | Kim、Lee、Gupta (2020) |
| **交互作用固定効果** | Xu (2017) |
| **行列補完** | Athey ら (2021) |
| **合成差の差** | Arkhangelsky ら (2021)。`synthdid` の移植 |
| **時期のずれる拡張合成コントロール** | `multisynth` の移植 |
| **推定量の比較** | 最大 9 手法を並べて表示 |

各分析の出力には「方法の説明」と「データと設定ガイド」があります。

SDID、拡張合成コントロール、時期のずれる拡張合成コントロールは、検証した設定で公式パッケージと完全に一致（差 0）。他の手法は `augsynth`、`gsynth`、`MCPanel`、`pensynth` などと照合しています（`tests/test_methods.R`）。均衡パネルが必要です。サンプルデータ：カリフォルニア州 Proposition 99、西ドイツ統一、最低賃金と十代雇用。

## インストール
- **jamovi ライブラリ**：jamovi 右上のモジュールボタンから *jamovi library* を開き、検索してインストール（モジュール公開後）。
- **サイドロード**：`.jmo` ファイルをダウンロードし、同じメニューの *Install from file* から導入します。`.jmo` はビルドした OS、プロセッサ、jamovi のバージョンでのみ動作します。
- **ソースから**：R と `jmvtools` が必要です：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（パスはお使いの jamovi に合わせてください）。

## 引用とライセンス
引用：Sun, M. *scmfamily* (jamovi module). https://github.com/mingchao-sun/scmfamily 。使用した手法の原論文（出力に記載）も併せて引用してください。

ライセンス：`DESCRIPTION` を参照（GPL）。第三者のコードとデータは `NOTICE` に記載。不具合や提案：https://github.com/mingchao-sun/scmfamily/issues
