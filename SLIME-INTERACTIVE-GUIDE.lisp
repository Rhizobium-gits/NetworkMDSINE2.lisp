;;;; MDSINE2 SLIME Interactive Guide
;;;; GNU Emacs + SLIME環境での使用ガイド

;;;; ============================================================
;;;; Emacs SLIME Setup Instructions
;;;; ============================================================

;; 1. ~/.emacs.d/init.el に以下を追加:
;;
;; (use-package slime
;;   :ensure t
;;   :config
;;   (setq inferior-lisp-program "sbcl")
;;   (slime-setup '(slime-fancy slime-autodoc slime-company)))
;;
;; (use-package company
;;   :ensure t)

;; 2. ターミナルでこれらのファイルをロード:
;;
;; M-x slime  (SLIMEを起動)
;;
;; その後、SLIMEバッファで:
;;
;; CL-USER> (load "/home/claude/mdsine2-common-lisp.lisp")
;; CL-USER> (load "/home/claude/mdsine2-advanced.lisp")

;;;; ============================================================
;;;; Quick Start: Interactive Commands for SLIME REPL
;;;; ============================================================

;; === Step 1: データロード ===
;; C-c C-c で以下のS式を評価:

(in-package :mdsine2)

(defvar *experiment* (load-lgim-data "/home/claude/LGIM.csv"))

(print-experiment-summary *experiment*)

;; === Step 2: 実験概要の確認 ===
;; 実行結果:
;; === LGIM Experiment Summary ===
;; Total Samples: 140
;; Bacteria Species: 123
;; Gravity Conditions: (0g 1_6g 1g 1g_s 5g baseline)
;; Time Points: (0h 8h 16h 24h)
;; Donors: (1 2 3)

;; === Step 3: 正規化と基本統計 ===

(defvar *normalized* (normalize-abundances *experiment*))

;; 最初のサンプルを確認:
(first *normalized*)

;; === Step 4: 重力応答性菌種の検出 ===

(defvar *responsive* (gravity-responsive-species *experiment* :threshold 0.05))

(print-responsive-species *responsive* 20)

;; === Step 5: 特定菌種の統計比較 ===

;; Bacteroidesについて重力条件間での比較
(statistical-comparison *experiment* "Bacteroides")

;; Faecalibacteriumについて
(statistical-comparison *experiment* "Faecalibacterium")

;; === Step 6: 時系列動態の分析 ===

;; ドナー1, 1g重力下での時系列
(defvar *dynamics-1g* (temporal-dynamics *experiment* "1g" 1))

;; ドナー1, 無重力下での時系列
(defvar *dynamics-0g* (temporal-dynamics *experiment* "0g" 1))

;; === Step 7: Advanced Analysis モジュールの利用 ===

(in-package :mdsine2-advanced)

;; 多様性メトリクスの計算
;; (最初のサンプルの豊富度を取得)
(defvar *first-abundances* 
  (getf (first (mdsine2:normalize-abundances *experiment*)) :relative-abundances))

(calculate-shannon-diversity *first-abundances*)

(calculate-simpson-diversity *first-abundances*)

;; === Step 8: サンプル間距離の計算 ===

(defvar *sample-2nd* 
  (getf (second (mdsine2:normalize-abundances *experiment*)) :relative-abundances))

(calculate-bray-curtis-distance *first-abundances* *sample-2nd*)

;; === Step 9: 成長率と相互作用ネットワーク ===

(defvar *network* 
  (estimate-growth-interaction-network 
   *dynamics-1g* 
   (list "Bacteroides" "Faecalibacterium" "Roseburia" "Blautia")
   :threshold 0.05))

;; === Step 10: 時間安定性の評価 ===

(temporal-stability *dynamics-1g*)

;;;; ============================================================
;;;; Advanced Analysis Workflows
;;;; ============================================================

;; === Workflow A: 重力条件間の詳細比較 ===

(in-package :mdsine2-advanced)

;; 1. 各重力条件のデータ抽出
(defun extract-gravity-group (gravity-string)
  (remove-if-not (lambda (s) (string= (getf s :gravity) gravity-string))
                 *normalized*))

(defvar *group-0g* (extract-gravity-group "0g"))
(defvar *group-1-6g* (extract-gravity-group "1_6g"))
(defvar *group-1g* (extract-gravity-group "1g"))
(defvar *group-1g-s* (extract-gravity-group "1g_s"))
(defvar *group-5g* (extract-gravity-group "5g"))

;; 2. 各グループ内での多様性
(defvar *shannon-0g* 
  (mapcar (lambda (s) (calculate-shannon-diversity (getf s :relative-abundances)))
          *group-0g*))

(defvar *shannon-1g* 
  (mapcar (lambda (s) (calculate-shannon-diversity (getf s :relative-abundances)))
          *group-1g*))

;; 平均多様性
(format t "Shannon Diversity - 0g: ~,3F~%" 
        (/ (apply #'+ *shannon-0g*) (length *shannon-0g*)))
(format t "Shannon Diversity - 1g: ~,3F~%" 
        (/ (apply #'+ *shannon-1g*) (length *shannon-1g*)))

;; 3. Effect sizeの計算
(gravity-effect-size *group-0g* *group-1g* "Bacteroides")
(gravity-effect-size *group-0g* *group-1g* "Faecalibacterium")

;; === Workflow B: 時系列動態の統合解析 ===

;; 1. 各重力条件での時系列を取得
(defvar *ts-0g-d1* (mdsine2:temporal-dynamics *experiment* "0g" 1))
(defvar *ts-1-6g-d1* (mdsine2:temporal-dynamics *experiment* "1_6g" 1))
(defvar *ts-1g-d1* (mdsine2:temporal-dynamics *experiment* "1g" 1))

;; 2. 安定性比較
(format t "~%Temporal Stability Comparison (Donor 1):~%")
(format t "0g: ~A~%" (temporal-stability *ts-0g-d1*))
(format t "1_6g: ~A~%" (temporal-stability *ts-1-6g-d1*))
(format t "1g: ~A~%" (temporal-stability *ts-1g-d1*))

;; === Workflow C: Community State Type 分類 ===

;; コミュニティ状態型を3つのカテゴリに分類
(defvar *cst* (identify-community-state-types *normalized* :k 3))

(format t "~%Community State Type Assignment:~%")
(format t "Assignments: ~A~%" (getf *cst* :assignments))
(format t "Iterations to convergence: ~A~%" (getf *cst* :iterations))

;; === Workflow D: Bayesian Prior Estimation ===

;; Bacteroidesについて、全サンプルの豊富度分布を事前分布として利用
(defvar *bacteroides-abundances*
  (mapcar (lambda (s) (gethash "Bacteroides" (getf s :relative-abundances) 0))
          *normalized*))

(defvar *prior-mean* (/ (apply #'+ *bacteroides-abundances*) 
                        (length *bacteroides-abundances*)))
(defvar *prior-variance* (/ (apply #'+ (mapcar (lambda (x) (expt (- x *prior-mean*) 2))
                                              *bacteroides-abundances*))
                            (length *bacteroides-abundances*)))

;; 新しい観測データに対してBayesian更新
(defvar *new-obs* (subseq *bacteroides-abundances* 0 10))

(bayesian-abundance-model *new-obs* *prior-mean* *prior-variance*)

;;;; ============================================================
;;;; Data Export and Visualization
;;;; ============================================================

;; === R連携のためのデータエクスポート ===

(defun export-for-r (filename)
  "R用にデータをCSVで出力"
  (with-open-file (stream filename :direction :output :if-exists :supersede)
    (format stream "Sample,Gravity,Time,Donor,Shannon_Diversity,Simpson_Diversity~%")
    (loop for sample in *normalized*
          for sha in (mapcar (lambda (s) (calculate-shannon-diversity 
                                         (getf s :relative-abundances)))
                            *normalized*)
          for sim in (mapcar (lambda (s) (calculate-simpson-diversity 
                                         (getf s :relative-abundances)))
                            *normalized*)
          do (format stream "~A,~A,~A,~A,~,4F,~,4F~%"
                    (getf sample :sample-id)
                    (getf sample :gravity)
                    (getf sample :time-point)
                    (getf sample :donor)
                    sha
                    sim)))
  (format t "Data exported to ~A~%" filename))

;; 使用:
(export-for-r "/home/claude/mdsine2-diversity-metrics.csv")

;;;; ============================================================
;;;; SLIME Integration Tips
;;;; ============================================================

;; 1. インクリメンタル開発:
;;    - ファイルの関数定義にカーソルを置いて C-c C-c で評価
;;    - REPL で変数に結果を束縛して即座に実験

;; 2. 対話的なデバッグ:
;;    - (trace 関数名) で関数の実行をトレース
;;    - (untrace) でトレースを解除

;; 3. ドキュメント検索:
;;    - C-c C-d c で関数の定義を検索
;;    - C-c C-d d でドキュメントを表示

;; 4. マクロ展開:
;;    - C-c C-m で S-式のマクロ展開を表示

;; 5. プロファイリング:
;;    - (time (解析処理)) で実行時間を計測

;;;; ============================================================
;;;; Statistical Testing Functions (to be extended)
;;;; ============================================================

(defun t-test-two-sample (sample1 sample2)
  "2標本t検定を実行（簡易版）"
  (let* ((n1 (length sample1))
         (n2 (length sample2))
         (mean1 (/ (apply #'+ sample1) n1))
         (mean2 (/ (apply #'+ sample2) n2))
         (var1 (/ (apply #'+ (mapcar (lambda (x) (expt (- x mean1) 2)) sample1)) n1))
         (var2 (/ (apply #'+ (mapcar (lambda (x) (expt (- x mean2) 2)) sample2)) n2))
         (pooled-var (/ (+ (* (- n1 1) var1) (* (- n2 1) var2)) (+ n1 n2 -2)))
         (t-stat (/ (- mean1 mean2) (sqrt (* pooled-var (+ (/ 1 n1) (/ 1 n2)))))))
    (list :t-statistic t-stat
          :mean1 mean1
          :mean2 mean2
          :n1 n1
          :n2 n2)))

;;;; ============================================================
;;;; Next Steps for Further Development
;;;; ============================================================

;; 1. Full Bayesian MCMC implementation for growth rate estimation
;; 2. Phylogenetic distance matrix integration
;; 3. Functional prediction from 16S data (PICRUSt-like)
;; 4. Statistical significance testing with multiple correction
;; 5. Visualization with gnuplot integration
;; 6. Integration with external databases (SILVA, RDP)

;;;; EOF
