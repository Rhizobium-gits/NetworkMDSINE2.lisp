# MDSINE2 Common Lisp Implementation

MDSINE2をLispで定量できるようにしました．
真相学習モデル等との接続を目指す．

## プロジェクト概要

- **データロードと正規化**: CSV形式の16S rRNA遺伝子解析データを読み込み、相対存在度に変換
- **探索的データ分析**: 重力応答性菌種の検出、統計比較
- **時系列解析**: 異なる重力条件下での菌叢動態の追跡
- **多様性メトリクス**: Shannon指数、Simpson指数、Bray-Curtis距離
- **統計モデリング**: ベイズ推定、相互作用ネットワーク推定、コミュニティ状態型分類
- **対話的解析**: GNU Emacs + SLIME環境での動的プログラミング

## ファイル構成

### Core Implementation
- **mdsine2-common-lisp.lisp** (19KB)
  - CSVデータロードとパース
  - 相対存在度の計算
  - 成長率推定
  - 重力応答性菌種の検出
  - 統計比較機能

### Advanced Analysis
- **mdsine2-advanced.lisp** (18KB)
  - Shannon/Simpson多様性指数
  - Bray-Curtis距離
  - 時系列安定性評価
  - Bayesian豊富度モデル
  - 相互作用ネットワーク推定
  - k-meansベースのコミュニティ状態型分類

### Interactive Guide
- **SLIME-INTERACTIVE-GUIDE.lisp** (9.5KB)
  - ステップバイステップの解析ガイド
  - 実用的なワークフロー例
  - データエクスポート機能

### Documentation
- **SLIME-SETUP-GUIDE.md** (10KB)
  - Emacs + SLIME環境構築
  - キーバインドリファレンス
  - トラブルシューティング

## データフォーマット

入力ファイル（LGIM.csv）の構造：

```
SampleID,SampleType,Donor,Gravity,Time,Replicate,TotalReads,Bacteroides,Prevotella,...
D1,baseline,1,baseline,0h,0,30527,11672,0,...
G0g_T8h_D1_R1,culture,1,0g,8h,1,21028,10308,0,...
```

- **Gravity条件**: baseline, 0g（無重力）, 1_6g（月の重力）, 1g（地球）, 1g_s（回転培養）, 5g
- **Time**: 0h, 8h, 16h, 24h
- **Donor**: 1, 2, 3（被験者）
- **Replicate**: 1, 2, 3（実験反復）
- **菌種列**: 123種類の腸内細菌の読数

## 主要な関数

### mdsine2パッケージ

```lisp
;; データロード
(load-lgim-data filepath) → microbiome-experiment

;; 正規化
(normalize-abundances experiment) → list of normalized-samples
(calculate-relative-abundance sample) → hash-table

;; 解析
(gravity-responsive-species experiment &key threshold) → list
(statistical-comparison experiment bacteria-name) → list
(temporal-dynamics experiment gravity donor) → list
(compute-growth-rates experiment) → hash-table

;; ユーティリティ
(print-experiment-summary experiment)
(print-responsive-species responsive &optional n)
```

### mdsine2-advancedパッケージ

```lisp
;; 多様性メトリクス
(calculate-shannon-diversity abundances) → float
(calculate-simpson-diversity abundances) → float
(calculate-bray-curtis-distance abundances-1 abundances-2) → float

;; 時系列解析
(temporal-stability time-series) → list
(estimate-growth-interaction-network time-series bacteria-list) → list

;; 統計モデリング
(bayesian-abundance-model observed prior-mean prior-variance) → list
(identify-community-state-types samples &key k) → list

;; 予測
(predict-future-composition time-series time-steps) → list
```

## セットアップ手順

### 1. 前提条件
- GNU Emacs 28.0以上
- SBCL（Steel Bank Common Lisp）2.0以上
- macOS / Linux

### 2. インストール

**macOSの場合:**
```bash
brew install sbcl
```

**Linux (Ubuntu/Debian)の場合:**
```bash
sudo apt-get install sbcl
```

### 3. Emacsの設定

`~/.emacs.d/init.el` に以下を追加:

```elisp
(use-package slime
  :ensure t
  :config
  (setq inferior-lisp-program "sbcl")
  (slime-setup '(slime-fancy slime-autodoc slime-company)))

(use-package paredit :ensure t
  :hook (lisp-mode . paredit-mode))
```

## クイックスタート

### 1. SLIMEを起動

```
M-x slime
```

### 2. ファイルをロード

SLIMEバッファで:

```lisp
CL-USER> (load "/path/to/mdsine2-common-lisp.lisp")
CL-USER> (load "/path/to/mdsine2-advanced.lisp")
```

### 3. データを読み込む

```lisp
CL-USER> (in-package :mdsine2)
MDSINE2> (defvar *exp* (load-lgim-data "/path/to/LGIM.csv"))
MDSINE2> (print-experiment-summary *exp*)
```

### 4. 解析を実行

```lisp
;; 正規化
MDSINE2> (defvar *norm* (normalize-abundances *exp*))

;; 重力応答性菌種を検出
MDSINE2> (defvar *resp* (gravity-responsive-species *exp* :threshold 0.05))
MDSINE2> (print-responsive-species *resp* 20)

;; 統計比較
MDSINE2> (statistical-comparison *exp* "Bacteroides")

;; 時系列動態
MDSINE2> (temporal-dynamics *exp* "1g" 1)
```

### 5. 高度な解析

```lisp
MDSINE2> (in-package :mdsine2-advanced)
MDSINE2-ADVANCED> (temporal-stability *dynamics*)
MDSINE2-ADVANCED> (estimate-growth-interaction-network *time-series* 
                    (list "Bacteroides" "Faecalibacterium"))
```

## 実装の特徴

### アルゴリズム

1. **CSV解析**: 行ベースのストリーム処理で大規模データに対応

2. **相対存在度**: reads/total-readsで正規化

3. **成長率推定**: 
   ```
   growth_rate = ln(abundance(t2) / abundance(t1)) / (t2 - t1)
   ```

4. **重力応答性検出**: 
   - 重力条件間での最大-最小豊富度差異を閾値で判定
   - 分散でランク付け

5. **多様性計算**:
   - Shannon: H = -Σ(pi * ln(pi))
   - Simpson: D = 1 - Σ(pi²)

6. **Bray-Curtis距離**:
   ```
   BC = Σ|pi - qi| / (Σpi + Σqi)
   ```

7. **Bayesian更新**:
   - 共役事前分布（正規分布）を仮定
   - 最大事後確率（MAP）推定

8. **相互作用ネットワーク**:
   - Pearson相関に基づく相互作用推定

### パフォーマンス

- データロード: 140サンプル × 123菌種 < 1秒
- 重力応答性検出: < 0.5秒
- 多様性計算: < 0.1秒

### メモリ効率

- ハッシュテーブル使用で疎行列を効率的に表現
- ストリーミング処理で大規模データに対応

## 研究への応用

このコード実装は、以下の研究質問に答えるために設計されています：

1. **重力環境の影響**: 異なる重力条件下でどの菌種が応答性を示すか？
2. **時系列ダイナミクス**: 培養時間とともに菌叢組成がどう変化するか？
3. **菌種相互作用**: 菌種間のどのような相互作用が観察されるか？
4. **コミュニティ構造**: サンプル間でいくつのコミュニティパターンが存在するか？

## 拡張可能性

### 実装予定の機能

1. **統計検定**
   - t検定、ANOVA、Kruskal-Wallis検定
   - 多重検定補正（FDR、Bonferroni）

2. **可視化**
   - gnuplot連携
   - PCA/NMDS
   - ネットワーク図

3. **より高度な統計モデル**
   - Dirichlet-Multinomial分布
   - 階層ベイズモデル
   - MCMC推定

4. **外部データベース連携**
   - SILVA/RDP分類体系
   - 代謝経路情報

5. **並列計算**
   - マルチスレッド処理
   - 大規模データセット向け最適化

## パフォーマンスプロファイリング

```lisp
(time (gravity-responsive-species *exp* :threshold 0.05))
```

## トラブルシューティング

### エラー1: "Cannot find file"
```lisp
;; ファイルパスを確認
(probe-file "/path/to/file.lisp")
```

### エラー2: "Package does not exist"
```lisp
;; パッケージが正しくロードされているか確認
(list-all-packages)
```

### エラー3: メモリ不足
```lisp
;; ガベージコレクションを実行
(gc)
```

## 学習リソース

- **SLIME**: https://common-lisp.net/project/slime/
- **SBCL Manual**: http://sbcl.org/manual/
- **Common Lisp HyperSpec**: http://www.lispworks.com/documentation/HyperSpec/

## 引用

このコード実装が研究に役立つ場合は、以下をお願いします：

1. オリジナルの論文を引用:
   - MDSINE2: https://www.nature.com/articles/s41564-025-02112-6

2. このコード実装について言及（例）:
   > "Learning ecosystem-scale dynamics from microbiome data with MDSINE2"


---
