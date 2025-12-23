# MDSINE2 Common Lisp Implementation - Emacs SLIME セットアップガイド

## 概要

このプロジェクトは、重力条件下での腸内マイクロバイオータ動態解析（MDSINE2）をCommon Lispで実装しています。GNU EmacsとSLIMEを使用して対話的に解析を進行できます。

## 必要なインストール

### 1. Common Lisp処理系のインストール

#### macOS (Homebrew)
```bash
brew install sbcl  # Steel Bank Common Lisp
```

#### Linux (Ubuntu/Debian)
```bash
sudo apt-get install sbcl
```

#### 確認
```bash
sbcl --version
```

### 2. Emacs パッケージのセットアップ

`~/.emacs.d/init.el` (または `~/.config/emacs/init.el`) に以下を追加:

```elisp
;; use-package の設定
(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
(package-initialize)

(unless (package-installed-p 'use-package)
  (package-refresh-contents)
  (package-install 'use-package))

;; SLIME の設定
(use-package slime
  :ensure t
  :config
  (setq inferior-lisp-program "sbcl")
  ;; slime-fancy: フル機能SLIME (推奨)
  ;; slime-autodoc: 関数の引数をミニバッファに表示
  ;; slime-company: 補完機能
  (slime-setup '(slime-fancy slime-autodoc slime-company)))

;; Company モード (補完)
(use-package company
  :ensure t
  :config
  (global-company-mode 1)
  (setq company-backends '((company-dabbrev company-files))))

;; Paredit (括弧操作の便利化)
(use-package paredit
  :ensure t
  :hook (lisp-mode . paredit-mode)
  :config
  (define-key paredit-mode-map (kbd "C-c C-c") nil))
```

### 3. Emacs起動時にロード

Emacsを起動後、以下を実行:

```
M-x package-refresh-contents
M-x package-install RET slime RET
M-x package-install RET paredit RET
M-x package-install RET company RET
```

## SLIMEの起動と基本操作

### 1. SLIMEの起動

Emacs内で:
```
M-x slime
```

すると別ウィンドウでREPLが起動します。

### 2. 基本的なキーバインド

| キーバインド | 動作 |
|------------|------|
| `C-c C-c` | カーソル位置のS式を評価 |
| `C-c C-k` | 現在のバッファをコンパイル |
| `C-c C-l` | ファイルをロード |
| `C-c C-d h` | 関数のドキュメント |
| `M-p` / `M-n` | REPL履歴の前後 |
| `C-c C-b` | 現在の式をハイライト |
| `C-c M-q` | S式を整形 |

### 3. ファイルのロード

SLIMEバッファで:

```lisp
CL-USER> (load "/home/claude/mdsine2-common-lisp.lisp")
```

または、`mdsine2-common-lisp.lisp`のバッファを開いて:
```
C-c C-k
```

## 実装ファイル構成

```
/home/claude/
├── LGIM.csv                          # 実験データ
├── mdsine2-common-lisp.lisp          # メインの解析モジュール
├── mdsine2-advanced.lisp             # 高度な統計解析
├── SLIME-INTERACTIVE-GUIDE.lisp      # インタラクティブガイド
└── SLIME-SETUP-GUIDE.md              # このファイル
```

## クイックスタート

### ステップ1: データロード

Emacsで `SLIME-INTERACTIVE-GUIDE.lisp` を開き、先頭のセクションを順に実行:

```lisp
(in-package :mdsine2)
(defvar *experiment* (load-lgim-data "/home/claude/LGIM.csv"))
(print-experiment-summary *experiment*)
```

各ステップを `C-c C-c` で実行すると、段階的に解析が進みます。

### ステップ2: 探索的データ分析

```lisp
(defvar *normalized* (normalize-abundances *experiment*))
(defvar *responsive* (gravity-responsive-species *experiment* :threshold 0.05))
(print-responsive-species *responsive* 20)
```

### ステップ3: 統計比較

```lisp
(statistical-comparison *experiment* "Bacteroides")
(statistical-comparison *experiment* "Faecalibacterium")
```

### ステップ4: 時系列解析

```lisp
(defvar *dynamics-1g* (temporal-dynamics *experiment* "1g" 1))
(defvar *dynamics-0g* (temporal-dynamics *experiment* "0g" 1))

;; 多様性の計算
(in-package :mdsine2-advanced)
(temporal-stability *dynamics-1g*)
```

## 主要な関数リファレンス

### mdsine2 パッケージ

#### データロードと正規化
- `(load-lgim-data filepath)` → microbiome-experiment
  - CSVファイルを読み込み、実験データ構造を返す
  
- `(normalize-abundances experiment)` → list
  - 全サンプルを相対存在度に正規化

#### 探索的解析
- `(gravity-responsive-species experiment &key threshold)` → list
  - 重力応答性を示す菌種を検出
  
- `(statistical-comparison experiment bacteria-name)` → list
  - 指定菌種について重力条件間での統計比較

- `(temporal-dynamics experiment gravity-condition donor)` → list
  - 時系列動態データを取得

#### ユーティリティ
- `(print-experiment-summary experiment)`
  - 実験概要を表示

- `(print-responsive-species responsive &optional n)`
  - 重力応答性菌種をリスト表示

### mdsine2-advanced パッケージ

#### 多様性メトリクス
- `(calculate-shannon-diversity abundances-hash)` → float
  
- `(calculate-simpson-diversity abundances-hash)` → float
  
- `(calculate-bray-curtis-distance abundances-1 abundances-2)` → float

#### 時系列解析
- `(temporal-stability time-series-abundances)` → list
  - 動態の安定性を評価

#### 相互作用ネットワーク
- `(estimate-growth-interaction-network time-series bacteria-list &key threshold)` → list
  - 菌種間の相互作用を推定

#### 統計モデル
- `(bayesian-abundance-model observed-data prior-mean prior-variance)` → list
  - ベイズ事後分布を計算

- `(identify-community-state-types normalized-samples &key k)` → list
  - k-meansでコミュニティ状態型を分類

## 実用的な解析ワークフロー

### Workflow 1: 重力条件間の比較解析

```lisp
(in-package :mdsine2-advanced)

;; 1. グループ抽出
(defvar *group-0g* 
  (remove-if-not (lambda (s) (string= (getf s :gravity) "0g"))
                 *normalized*))

;; 2. 多様性計算
(defvar *shannon* 
  (mapcar (lambda (s) (calculate-shannon-diversity (getf s :relative-abundances)))
          *group-0g*))

;; 3. 要約統計
(format t "Mean Shannon: ~,3F~%" 
        (/ (apply #'+ *shannon*) (length *shannon*)))

;; 4. Effect Size
(gravity-effect-size *group-0g* 
                     (remove-if-not (lambda (s) (string= (getf s :gravity) "1g")) *normalized*)
                     "Bacteroides")
```

### Workflow 2: 時系列安定性の比較

```lisp
;; 各重力条件での安定性を比較
(defvar *stability-results* '())

(dolist (gravity '("0g" "1_6g" "1g" "1g_s" "5g"))
  (let ((ts (mdsine2:temporal-dynamics *experiment* gravity 1)))
    (push (list :gravity gravity
                :stability (temporal-stability ts))
          *stability-results*)))

;; 結果表示
(dolist (result *stability-results*)
  (format t "~A: ~A~%"
          (getf result :gravity)
          (getf result :stability)))
```

### Workflow 3: データエクスポートとR連携

```lisp
(defun export-abundance-data (filename)
  "Rで使用可能なCSV形式でエクスポート"
  (with-open-file (stream filename :direction :output :if-exists :supersede)
    (format stream "SampleID,Gravity,TimePoint,Donor,Bacteroides,Faecalibacterium~%")
    (loop for sample in *normalized*
          do (format stream "~A,~A,~A,~A,~,6F,~,6F~%"
                    (getf sample :sample-id)
                    (getf sample :gravity)
                    (getf sample :time-point)
                    (getf sample :donor)
                    (gethash "Bacteroides" (getf sample :relative-abundances) 0)
                    (gethash "Faecalibacterium" (getf sample :relative-abundances) 0)))))

(export-abundance-data "/home/claude/abundant-species.csv")
```

## デバッグテクニック

### 1. トレース実行

```lisp
(trace load-lgim-data)
(defvar *exp* (load-lgim-data "/home/claude/LGIM.csv"))
(untrace)
```

### 2. 式の評価を段階的に実行

```lisp
;; 括弧の対応を確認
(paredit-mode 1)

;; 式をステップ実行
C-c C-s  ;; sbcl固有（SLDB）
```

### 3. エラーメッセージの理解

SLIMEは詳細なエラースタックトレースを表示します。`q`で閉じられます。

## パフォーマンス最適化

### 大規模データの処理

```lisp
;; コンパイル時の最適化レベルを上げる
(declaim (optimize (speed 3) (safety 0)))

;; 処理時間の計測
(time (gravity-responsive-species *experiment* :threshold 0.05))
```

### メモリ使用量の監視

```lisp
;; SLIMEバッファで
,show processes  ;; プロセス情報表示
(gc)             ;; ガベージコレクション実行
```

## よくある問題と解決

### 問題1: SBCLが見つからない

```bash
# インストール確認
which sbcl

# パスを設定 (~/.emacs.d/init.el)
(setq inferior-lisp-program "/usr/local/bin/sbcl")
```

### 問題2: ロードエラー

```lisp
;; ファイルパスを確認
(probe-file "/home/claude/mdsine2-common-lisp.lisp")

;; 相対パスではなく絶対パスを使用
```

### 問題3: 日本語文字化け

```elisp
;; init.el で文字コード指定
(set-language-environment "UTF-8")
(set-default-coding-systems 'utf-8)
```

## 推奨される次のステップ

1. **統計検定の追加実装**
   - t検定、ANOVA
   - FDR補正

2. **可視化機能の拡張**
   - gnuplot連携
   - PCA/NMDS

3. **Rとの連携**
   - Rデータフレーム形式でのエクスポート
   - ggplot2での可視化

4. **並列計算**
   - マルチスレッド処理（SBCL）
   - 大規模データセットの高速化

## 参考リソース

- SLIME Documentation: https://common-lisp.net/project/slime/
- Common Lisp HyperSpec: http://www.lispworks.com/documentation/HyperSpec/
- SBCL Manual: http://sbcl.org/manual/

## ご質問・フィードバック

このコード実装について質問がある場合は、以下を参照してください：

- 元のMDSINE2論文: https://doi.org/10.1186/s40168-018-0612-3
- 16S rRNA解析のベストプラクティス

---

**作成日**: 2025年12月23日
**対応環境**: GNU Emacs 28.0+, SBCL 2.0+, macOS/Linux
**ライセンス**: MIT

EOF
