;;;; MDSINE2 Advanced Analysis Module
;;;; 重力条件下での腸内マイクロバイオータ時系列解析
;;;; Gravity-dependent microbiota dynamics with Bayesian inference

(defpackage :mdsine2-advanced
  (:use :cl)
  (:export
   #:estimate-growth-interaction-network
   #:identify-community-state-types
   #:calculate-shannon-diversity
   #:calculate-bray-curtis-distance
   #:temporal-stability
   #:gravity-effect-size
   #:bayesian-abundance-model
   #:predict-future-composition
   #:generate-analysis-report))

(in-package :mdsine2-advanced)

;;; ============================================================
;;; Diversity Metrics
;;; ============================================================

(defun calculate-shannon-diversity (abundances-hash)
  "Shannon多様性指数を計算
   H = -Σ(pi * ln(pi)) where pi is relative abundance"
  (let ((total (reduce #'+ (loop for v being the hash-values of abundances-hash
                                 collect v)))
        (shannon 0))
    (maphash (lambda (bacteria count)
               (let ((pi (/ (float count) (float total))))
                 (when (> pi 0)
                   (setf shannon (- shannon (* pi (log pi)))))))
             abundances-hash)
    shannon))

(defun calculate-simpson-diversity (abundances-hash)
  "Simpson多様性指数を計算
   D = 1 - Σ(pi^2)"
  (let ((total (reduce #'+ (loop for v being the hash-values of abundances-hash
                                 collect v)))
        (sum-pi-sq 0))
    (maphash (lambda (bacteria count)
               (let ((pi (/ (float count) (float total))))
                 (setf sum-pi-sq (+ sum-pi-sq (expt pi 2)))))
             abundances-hash)
    (- 1 sum-pi-sq)))

(defun calculate-bray-curtis-distance (abundances-1 abundances-2)
  "Bray-Curtis距離を計算
   BC = Σ|pi - qi| / (Σpi + Σqi)
   サンプル間の違いを定量化"
  (let ((total-1 (reduce #'+ (loop for v being the hash-values of abundances-1
                                   collect v)))
        (total-2 (reduce #'+ (loop for v being the hash-values of abundances-2
                                   collect v)))
        (sum-abs-diff 0))
    
    ;; 両方のサンプルに含まれるすべての菌について
    (let ((all-bacteria (make-hash-table :test #'equal)))
      (maphash (lambda (b c) (setf (gethash b all-bacteria) t)) abundances-1)
      (maphash (lambda (b c) (setf (gethash b all-bacteria) t)) abundances-2)
      
      (maphash (lambda (bacteria _)
                 (let* ((count-1 (gethash bacteria abundances-1 0))
                        (count-2 (gethash bacteria abundances-2 0))
                        (pi (/ (float count-1) (float total-1)))
                        (qi (/ (float count-2) (float total-2))))
                   (setf sum-abs-diff (+ sum-abs-diff (abs (- pi qi))))))
               all-bacteria))
    
    (if (= (+ total-1 total-2) 0)
        0
        (/ sum-abs-diff (+ (/ (float total-1) (float total-1))
                          (/ (float total-2) (float total-2)))))))

;;; ============================================================
;;; Temporal Stability Analysis
;;; ============================================================

(defun temporal-stability (time-series-abundances)
  "時系列マイクロバイオータの安定性を評価
   各時間点間のBray-Curtis距離の平均と変動から安定性を計算"
  (if (< (length time-series-abundances) 2)
      (list :stability 1.0 :message "Insufficient data")
      (let ((distances '()))
        ;; 連続する時間点間の距離を計算
        (loop for i from 0 below (- (length time-series-abundances) 1)
              do (let ((d (calculate-bray-curtis-distance 
                          (getf (nth i time-series-abundances) :abundances)
                          (getf (nth (1+ i) time-series-abundances) :abundances))))
                   (push d distances)))
        
        (let ((mean-distance (/ (apply #'+ distances) (length distances)))
              (max-distance (apply #'max distances))
              (min-distance (apply #'min distances)))
          (list :mean-distance mean-distance
                :max-distance max-distance
                :min-distance min-distance
                :stability (- 1 mean-distance))))))

;;; ============================================================
;;; Gravity Effect Size Estimation
;;; ============================================================

(defun gravity-effect-size (gravity-group-1 gravity-group-2 bacteria-name)
  "2つの重力条件間での指定菌種のeffect sizeを計算（Cohen's d）"
  (let* ((abundances-1 (mapcar (lambda (sample)
                                 (gethash bacteria-name 
                                         (getf sample :relative-abundances) 0))
                               gravity-group-1))
         (abundances-2 (mapcar (lambda (sample)
                                 (gethash bacteria-name 
                                         (getf sample :relative-abundances) 0))
                               gravity-group-2))
         (mean-1 (/ (apply #'+ abundances-1) (length abundances-1)))
         (mean-2 (/ (apply #'+ abundances-2) (length abundances-2)))
         (var-1 (/ (apply #'+ (mapcar (lambda (x) (expt (- x mean-1) 2)) abundances-1))
                  (length abundances-1)))
         (var-2 (/ (apply #'+ (mapcar (lambda (x) (expt (- x mean-2) 2)) abundances-2))
                  (length abundances-2)))
         (pooled-sd (sqrt (/ (+ var-1 var-2) 2))))
    
    (if (= pooled-sd 0)
        0
        (/ (- mean-1 mean-2) pooled-sd))))

;;; ============================================================
;;; Interaction Network Inference
;;; ============================================================

(defun estimate-growth-interaction-network (time-series-data bacteria-list &key (threshold 0.1))
  "時系列データから菌種間の相互作用ネットワークを推定
   相互作用強度 = 成長率の変化に対する他菌種豊富度の影響"
  (let ((interaction-matrix (make-hash-table :test #'equal))
        (growth-rates (make-hash-table :test #'equal)))
    
    ;; 各菌種の成長率を計算
    (loop for bacteria in bacteria-list
          do (let ((rates '()))
               (loop for i from 0 below (- (length time-series-data) 1)
                     do (let* ((t1 (nth i time-series-data))
                              (t2 (nth (1+ i) time-series-data))
                              (abun1 (gethash bacteria (getf t1 :abundances) 1e-6))
                              (abun2 (gethash bacteria (getf t2 :abundances) 1e-6)))
                         (when (and (> abun1 0) (> abun2 0))
                           (push (log (/ abun2 abun1)) rates))))
               (when rates
                 (setf (gethash bacteria growth-rates)
                       (list :mean (/ (apply #'+ rates) (length rates))
                            :values rates)))))
    
    ;; 相互作用を推定（簡易版：豊富度との相関）
    (loop for bacteria-i in bacteria-list
          do (let ((gi (getf (gethash bacteria-i growth-rates) :mean 0)))
               (loop for bacteria-j in bacteria-list
                     do (when (not (string= bacteria-i bacteria-j))
                          (let ((correlation (estimate-interaction-strength 
                                            time-series-data bacteria-i bacteria-j)))
                            (when (> (abs correlation) threshold)
                              (setf (gethash (format nil "~A->~A" bacteria-i bacteria-j)
                                            interaction-matrix)
                                    correlation)))))))
    
    (list :growth-rates growth-rates
          :interactions interaction-matrix)))

(defun estimate-interaction-strength (time-series bacteria-from bacteria-to)
  "菌種間の相互作用強度を推定（Pearson相関）"
  (let ((values-from '())
        (values-to '()))
    (loop for timepoint in time-series
          do (let ((vf (gethash bacteria-from (getf timepoint :abundances) 0))
                   (vt (gethash bacteria-to (getf timepoint :abundances) 0)))
               (push vf values-from)
               (push vt values-to)))
    
    (if (< (length values-from) 2)
        0
        (pearson-correlation (reverse values-from) (reverse values-to)))))

(defun pearson-correlation (x y)
  "Pearson相関係数を計算"
  (let* ((n (length x))
         (mean-x (/ (apply #'+ x) n))
         (mean-y (/ (apply #'+ y) n))
         (sum-xy (apply #'+ (mapcar #'* 
                                   (mapcar (lambda (v) (- v mean-x)) x)
                                   (mapcar (lambda (v) (- v mean-y)) y))))
         (sum-x2 (apply #'+ (mapcar (lambda (v) (expt (- v mean-x) 2)) x)))
         (sum-y2 (apply #'+ (mapcar (lambda (v) (expt (- v mean-y) 2)) y))))
    
    (if (or (= sum-x2 0) (= sum-y2 0))
        0
        (/ sum-xy (sqrt (* sum-x2 sum-y2))))))

;;; ============================================================
;;; Bayesian Abundance Model
;;; ============================================================

(defun bayesian-abundance-model (observed-data prior-mean prior-variance)
  "単純なベイズ豊富度モデル
   観測データと事前分布から事後分布を推定"
  (let* ((n (length observed-data))
         (data-mean (/ (apply #'+ observed-data) n))
         (data-variance (/ (apply #'+ (mapcar (lambda (x) (expt (- x data-mean) 2)) 
                                             observed-data))
                          n))
         ;; Bayesian update: posterior = (prior_precision * prior_mean + likelihood_precision * data_mean) / (prior_precision + likelihood_precision)
         (prior-precision (/ 1 prior-variance))
         (likelihood-precision (/ n data-variance))
         (posterior-mean (/ (+ (* prior-precision prior-mean)
                              (* likelihood-precision data-mean))
                           (+ prior-precision likelihood-precision)))
         (posterior-variance (/ 1 (+ prior-precision likelihood-precision))))
    
    (list :posterior-mean posterior-mean
          :posterior-variance posterior-variance
          :posterior-sd (sqrt posterior-variance)
          :data-mean data-mean
          :data-variance data-variance)))

;;; ============================================================
;;; Community State Type Classification
;;; ============================================================

(defun identify-community-state-types (normalized-samples &key (k 3))
  "k-means風のアルゴリズムでコミュニティ状態型を特定
   各サンプルの菌叢構成パターンを分類"
  (let* ((sample-vectors (mapcar #'sample-to-abundance-vector normalized-samples))
         (centroids (initialize-centroids sample-vectors k))
         (assignments (make-hash-table :test #'equal))
         (converged nil)
         (iteration 0))
    
    ;; 簡易k-means
    (loop until converged
          do (progn
               (setf iteration (1+ iteration))
               ;; サンプルを最近接centroidに割り当て
               (setf assignments (assign-samples-to-centroids sample-vectors centroids))
               ;; centroidを再計算
               (let ((new-centroids (recalculate-centroids sample-vectors assignments k)))
                 (setf converged (centroids-converged centroids new-centroids))
                 (setf centroids new-centroids)))
          when (> iteration 100)
          do (setf converged t))
    
    (list :centroids centroids
          :assignments assignments
          :iterations iteration)))

(defun sample-to-abundance-vector (normalized-sample)
  "正規化サンプルを豊富度ベクトル（リスト）に変換"
  (let ((abundances (getf normalized-sample :relative-abundances))
        (vector '()))
    (maphash (lambda (bacteria abundance)
               (push abundance vector))
             abundances)
    (reverse vector)))

(defun initialize-centroids (vectors k)
  "ランダムにk個のcentroidを初期化"
  (loop for i from 0 below k
        collect (nth (random (length vectors)) vectors)))

(defun assign-samples-to-centroids (vectors centroids)
  "各ベクトルを最近接centroidに割り当て"
  (let ((assignments (make-hash-table :test #'equal)))
    (loop for i from 0 below (length vectors)
          do (let* ((vector (nth i vectors))
                   (nearest-idx 0)
                   (min-dist (euclidean-distance vector (nth 0 centroids))))
               (loop for j from 1 below (length centroids)
                     do (let ((dist (euclidean-distance vector (nth j centroids))))
                          (when (< dist min-dist)
                            (setf min-dist dist)
                            (setf nearest-idx j))))
               (setf (gethash i assignments) nearest-idx)))
    assignments))

(defun euclidean-distance (v1 v2)
  "2つのベクトル間のユークリッド距離"
  (sqrt (apply #'+ (mapcar (lambda (x y) (expt (- x y) 2)) v1 v2))))

(defun recalculate-centroids (vectors assignments k)
  "各クラスタの新しいcentroidを計算"
  (let ((new-centroids '()))
    (loop for c from 0 below k
          do (let ((members '()))
               ;; クラスタcのメンバーを集める
               (maphash (lambda (idx cluster)
                          (when (= cluster c)
                            (push (nth idx vectors) members)))
                       assignments)
               ;; メンバーが存在する場合のみcentroidを計算
               (if members
                   (let* ((n (length members))
                          (dim (length (first members)))
                          (new-centroid (loop for d from 0 below dim
                                            collect (/ (apply #'+ 
                                                            (mapcar (lambda (m) (nth d m)) 
                                                                   members))
                                                       n))))
                     (push new-centroid new-centroids))
                   (push (nth c (car (last assignments))) new-centroids))))
    (reverse new-centroids)))

(defun centroids-converged (old-centroids new-centroids &key (tolerance 0.001))
  "centroidが収束したかチェック"
  (let ((all-converged t))
    (loop for i from 0 below (length old-centroids)
          do (let ((dist (euclidean-distance (nth i old-centroids) (nth i new-centroids))))
               (when (> dist tolerance)
                 (setf all-converged nil))))
    all-converged))

;;; ============================================================
;;; Predictive Modeling
;;; ============================================================

(defun predict-future-composition (time-series-data time-steps-ahead)
  "単純な外挿法を使って将来の菌叢組成を予測"
  (if (< (length time-series-data) 2)
      (list :error "Insufficient time points")
      (let* ((latest-sample (first (last time-series-data)))
             (previous-sample (first (last time-series-data 2)))
             (latest-abundances (getf latest-sample :abundances))
             (previous-abundances (getf previous-sample :abundances))
             (predicted (make-hash-table :test #'equal)))
        
        ;; 線形外挿：豊富度の変化率を計算して予測
        (maphash (lambda (bacteria latest-count)
                   (let* ((previous-count (gethash bacteria previous-abundances 1))
                          (change-rate (/ (float latest-count) (float previous-count)))
                          (predicted-count (* latest-count (expt change-rate time-steps-ahead))))
                     (setf (gethash bacteria predicted) predicted-count)))
                 latest-abundances)
        
        (list :predicted-composition predicted
              :method "exponential-extrapolation"
              :time-steps time-steps-ahead))))

;;; ============================================================
;;; Report Generation
;;; ============================================================

(defun generate-analysis-report (experiment-data analysis-results)
  "完全な解析レポートを生成"
  (with-output-to-string (output)
    (format output "~%")
    (format output "=====================================~%")
    (format output "MDSINE2 Analysis Report~%")
    (format output "Gravity-dependent Gut Microbiota Dynamics~%")
    (format output "=====================================~%~%")
    
    (format output "Analysis Summary:~%")
    (format output "- Total samples analyzed: ~A~%" 
            (length (getf experiment-data :samples)))
    (format output "- Gravity conditions: ~A~%" 
            (getf experiment-data :gravity-conditions))
    (format output "- Time points: ~A~%" 
            (getf experiment-data :time-points))
    (format output "- Bacterial species detected: ~A~%~%" 
            (length (getf experiment-data :bacteria-names)))
    
    (format output "Key Findings:~%")
    (loop for result in (getf analysis-results :findings)
          do (format output "  • ~A~%" result))
    (format output "~%")
    
    (format output "Recommendations:~%")
    (format output "  • Continue time-series monitoring for 48-72h~%")
    (format output "  • Perform functional metagenomic analysis~%")
    (format output "  • Validate gravity-responsive species in space conditions~%~%")
    
    (format output "Report generated: ~A~%" (get-universal-time))))

;;; EOF
