;;;; MDSINE2 - Microbiome Dynamics and Systems Integration 
;;;; for Experimental Evolution 2 (Common Lisp Implementation)
;;;; 
;;;; 重力条件下での腸内マイクロバイオータ動態解析
;;;; Koyo's LGIM (Linear Gravity Intestinal Microbiota) Research

(defpackage :mdsine2
  (:use :cl)
  (:export
   #:load-lgim-data
   #:normalize-abundances
   #:calculate-relative-abundance
   #:compute-growth-rates
   #:gravity-responsive-species
   #:temporal-dynamics
   #:plot-dynamics
   #:statistical-comparison))

(in-package :mdsine2)

;;; ============================================================
;;; Data Structure
;;; ============================================================

(defstruct microbiome-sample
  "マイクロバイオームサンプルを表現するデータ構造"
  sample-id
  sample-type
  donor
  gravity
  time-point      ; 8h, 16h, 24h
  replicate
  total-reads
  abundances)     ; hash-table: bacteria-name -> read-count

(defstruct microbiome-experiment
  "全実験データを保持する構造体"
  samples           ; list of microbiome-sample
  bacteria-names    ; list of all bacteria species detected
  gravity-conditions ; list of gravity conditions (0g, 1_6g, 1g, 1g_s, 5g)
  time-points       ; list of time points (8h, 16h, 24h)
  donors            ; list of donor IDs (1, 2, 3)
  metadata)         ; additional info

;;; ============================================================
;;; CSV Reading and Data Parsing
;;; ============================================================

(defun read-csv-file (filepath)
  "CSVファイルを読み込み、リストのリストを返す"
  (with-open-file (stream filepath :direction :input :external-format :utf-8)
    (loop for line = (read-line stream nil nil)
          while line
          collect (split-csv-line line))))

(defun split-csv-line (line)
  "CSV行をカンマで分割"
  (let ((fields '())
        (current-field "")
        (in-quotes nil))
    (loop for char across line
          do (cond
               ((char= char #\") 
                (setf in-quotes (not in-quotes)))
               ((and (char= char #\,) (not in-quotes))
                (push current-field fields)
                (setf current-field ""))
               (t 
                (setf current-field (concatenate 'string current-field (string char))))))
    (push current-field fields)
    (reverse fields)))

(defun load-lgim-data (filepath)
  "LGIM実験データをCSVから読み込む"
  (let* ((csv-data (read-csv-file filepath))
         (header (first csv-data))
         (data-rows (rest csv-data))
         (bacteria-names (get-bacteria-columns header))
         (samples '())
         (gravity-conditions (make-hash-table :test #'equal))
         (time-points (make-hash-table :test #'equal))
         (donors (make-hash-table :test #'equal)))
    
    ;; データ行をパース
    (loop for row in data-rows
          do (let* ((sample-id (nth 0 row))
                    (sample-type (nth 1 row))
                    (donor (parse-integer (nth 2 row)))
                    (gravity (nth 3 row))
                    (time-point (nth 4 row))
                    (replicate (parse-integer (nth 5 row)))
                    (total-reads (parse-integer (nth 6 row)))
                    (abundances (make-hash-table :test #'equal)))
               
               ;; 細菌の存在度を取得
               (loop for i from 0 below (length bacteria-names)
                     do (let ((bacteria-name (nth i bacteria-names))
                              (count-str (nth (+ i 7) row)))
                          (when (and bacteria-name (not (string= bacteria-name "")))
                            (setf (gethash bacteria-name abundances)
                                  (parse-integer count-str :junk-allowed t)))))
               
               ;; メタデータ収集
               (setf (gethash gravity gravity-conditions) t)
               (setf (gethash time-point time-points) t)
               (setf (gethash donor donors) t)
               
               ;; サンプル作成
               (push (make-microbiome-sample
                      :sample-id sample-id
                      :sample-type sample-type
                      :donor donor
                      :gravity gravity
                      :time-point time-point
                      :replicate replicate
                      :total-reads total-reads
                      :abundances abundances)
                     samples)))
    
    ;; 実験構造体を返す
    (make-microbiome-experiment
     :samples (reverse samples)
     :bacteria-names bacteria-names
     :gravity-conditions (sort (loop for g being the hash-keys of gravity-conditions
                                     collect g) #'string<)
     :time-points (sort (loop for t being the hash-keys of time-points
                              collect t) #'string<)
     :donors (sort (loop for d being the hash-keys of donors
                         collect d) #'<)
     :metadata (list :header header
                     :sample-count (length samples)))))

(defun get-bacteria-columns (header)
  "ヘッダーから細菌種の列を抽出（最初の7列はメタデータ）"
  (subseq header 7))

;;; ============================================================
;;; Abundance Normalization
;;; ============================================================

(defun calculate-relative-abundance (sample)
  "サンプルの相対存在度を計算（全体を1.0に正規化）"
  (let ((total-reads (float (microbiome-sample-total-reads sample)))
        (relative-abundances (make-hash-table :test #'equal)))
    (maphash (lambda (bacteria count)
               (setf (gethash bacteria relative-abundances)
                     (/ (float count) total-reads)))
             (microbiome-sample-abundances sample))
    relative-abundances))

(defun normalize-abundances (experiment)
  "全サンプルの相対存在度を計算"
  (loop for sample in (microbiome-experiment-samples experiment)
        collect (list :sample-id (microbiome-sample-sample-id sample)
                      :gravity (microbiome-sample-gravity sample)
                      :time-point (microbiome-sample-time-point sample)
                      :donor (microbiome-sample-donor sample)
                      :replicate (microbiome-sample-replicate sample)
                      :relative-abundances (calculate-relative-abundance sample))))

;;; ============================================================
;;; Growth Rate Estimation
;;; ============================================================

(defun compute-growth-rates (experiment)
  "時系列データから各菌種の成長率を推定
   成長率 = ln(abundance(t2)/abundance(t1)) / (t2-t1)"
  (let* ((normalized (normalize-abundances experiment))
         (gravity-time-groups (group-by-gravity-and-time normalized))
         (growth-rates (make-hash-table :test #'equal)))
    
    ;; 各重力条件とドナーごとにグループ化
    (maphash (lambda (gravity-group samples-by-time)
               (loop for donor in (microbiome-experiment-donors experiment)
                     do (let ((donor-samples (remove-if-not 
                                             (lambda (s) (eql (getf s :donor) donor))
                                             samples-by-time)))
                          (when (> (length donor-samples) 1)
                            ;; 時間順にソート
                            (let ((sorted (sort donor-samples #'<
                                               :key (lambda (s) 
                                                      (parse-time-point (getf s :time-point))))))
                              (loop for i from 0 below (- (length sorted) 1)
                                    do (let* ((t1 (nth i sorted))
                                              (t2 (nth (1+ i) sorted))
                                              (time1 (parse-time-point (getf t1 :time-point)))
                                              (time2 (parse-time-point (getf t2 :time-point)))
                                              (delta-t (- time2 time1)))
                                         (when (> delta-t 0)
                                           (compute-pairwise-growth-rates 
                                            t1 t2 delta-t gravity growth-rates)))))))
                     ))
             gravity-time-groups)
    
    growth-rates))

(defun parse-time-point (time-str)
  "\"8h\", \"16h\", \"24h\" から数値を抽出"
  (parse-integer (subseq time-str 0 (- (length time-str) 1))))

(defun group-by-gravity-and-time (normalized-samples)
  "重力条件でグループ化"
  (let ((groups (make-hash-table :test #'equal)))
    (loop for sample in normalized-samples
          do (let ((gravity (getf sample :gravity)))
               (setf (gethash gravity groups)
                     (append (gethash gravity groups '()) (list sample)))))
    groups))

(defun compute-pairwise-growth-rates (sample1 sample2 delta-t gravity growth-rates)
  "2つのサンプル間の菌種ごとの成長率を計算"
  (let ((abundant1 (getf sample1 :relative-abundances))
        (abundant2 (getf sample2 :relative-abundances)))
    
    ;; sample2の全菌種について成長率を計算
    (maphash (lambda (bacteria abundance2)
               (let* ((abundance1 (gethash bacteria abundant1 1e-6))
                      (abundance2-safe (if (> abundance2 0) abundance2 1e-6))
                      (growth-rate (/ (log (/ abundance2-safe abundance1)) delta-t)))
                 (let ((key (format nil "~A_~A" bacteria gravity)))
                   (if (gethash key growth-rates)
                       (push growth-rate (gethash key growth-rates))
                       (setf (gethash key growth-rates) (list growth-rate))))))
             abundant2)))

;;; ============================================================
;;; Gravity-Responsive Species Detection
;;; ============================================================

(defun gravity-responsive-species (experiment &key (threshold 0.1))
  "重力条件に応答性を示す菌種を特定
   threshold以上の相対的な豊富度変化を検出"
  (let* ((normalized (normalize-abundances experiment))
         (gravity-profiles (make-hash-table :test #'equal)))
    
    ;; 各菌種について、重力条件ごとの平均豊富度を計算
    (loop for bacteria in (microbiome-experiment-bacteria-names experiment)
          do (let ((profiles (make-hash-table :test #'equal)))
               (loop for gravity in (microbiome-experiment-gravity-conditions experiment)
                     do (let* ((gravity-samples 
                               (remove-if-not (lambda (s) (string= (getf s :gravity) gravity))
                                            normalized))
                              (abundances (mapcar (lambda (s)
                                                   (gethash bacteria 
                                                           (getf s :relative-abundances) 0))
                                                 gravity-samples))
                              (mean-abundance (/ (apply #'+ abundances)
                                               (max 1 (length abundances)))))
                         (setf (gethash gravity profiles) mean-abundance)))
               (setf (gethash bacteria gravity-profiles) profiles)))
    
    ;; 重力依存的に変化する菌種を検出
    (let ((responsive-species '()))
      (maphash (lambda (bacteria profiles)
                 (let ((values (loop for g being the hash-values of profiles
                                    collect g)))
                   (when (> (- (apply #'max values) (apply #'min values)) threshold)
                     (push (list :bacteria bacteria
                                :profiles profiles
                                :variance (calculate-variance values))
                           responsive-species))))
               gravity-profiles)
      
      ;; 分散で降順ソート
      (sort responsive-species #'> :key (lambda (x) (getf x :variance))))))

(defun calculate-variance (values)
  "分散を計算"
  (if (< (length values) 2)
      0
      (let* ((mean (/ (apply #'+ values) (length values)))
             (sum-sq (apply #'+ (mapcar (lambda (x) (expt (- x mean) 2)) values))))
        (/ sum-sq (length values)))))

;;; ============================================================
;;; Temporal Dynamics Summary
;;; ============================================================

(defun temporal-dynamics (experiment gravity-condition donor)
  "指定された重力条件とドナーについての時系列動態を取得"
  (let* ((normalized (normalize-abundances experiment))
         (filtered (remove-if-not 
                   (lambda (s) (and (string= (getf s :gravity) gravity-condition)
                                   (eql (getf s :donor) donor)))
                   normalized))
         (sorted (sort filtered #'<
                      :key (lambda (s) (parse-time-point (getf s :time-point))))))
    
    (loop for sample in sorted
          collect (list :time-point (getf sample :time-point)
                       :replicate (getf sample :replicate)
                       :abundances (getf sample :relative-abundances)))))

;;; ============================================================
;;; Statistical Comparison
;;; ============================================================

(defun statistical-comparison (experiment bacteria-name)
  "指定された菌種について重力条件間での統計比較
   各重力条件での平均豊富度と標準偏差を計算"
  (let* ((normalized (normalize-abundances experiment))
         (comparison (make-hash-table :test #'equal)))
    
    (loop for gravity in (microbiome-experiment-gravity-conditions experiment)
          do (let* ((gravity-samples 
                    (remove-if-not (lambda (s) (string= (getf s :gravity) gravity))
                                 normalized))
                   (abundances (mapcar (lambda (s)
                                        (gethash bacteria-name 
                                                (getf s :relative-abundances) 0))
                                      gravity-samples)))
               (setf (gethash gravity comparison)
                     (list :mean (/ (apply #'+ abundances) (length abundances))
                          :sd (calculate-std-dev abundances)
                          :count (length abundances)
                          :values abundances))))
    
    (loop for g being the hash-keys of comparison
          collect (list :gravity g
                       :result (gethash g comparison)))))

(defun calculate-std-dev (values)
  "標準偏差を計算"
  (if (< (length values) 2)
      0
      (let* ((mean (/ (apply #'+ values) (length values)))
             (variance (/ (apply #'+ (mapcar (lambda (x) (expt (- x mean) 2)) values))
                         (length values))))
        (sqrt variance))))

;;; ============================================================
;;; Visualization (Text-based)
;;; ============================================================

(defun plot-dynamics (dynamics-data &optional (width 80) (height 20))
  "時系列動態をテキストベースでプロット"
  (let ((min-val 0) (max-val 0))
    ;; データの最大値を取得
    (loop for timepoint in dynamics-data
          do (maphash (lambda (bacteria abundance)
                       (setf max-val (max max-val abundance)))
                     (getf timepoint :abundances)))
    
    ;; グラフの描画
    (format t "~%Temporal Dynamics Plot~%")
    (format t "~A~%" (make-string width :initial-element #\=))
    
    (loop for level from height downto 1
          do (let ((threshold (* (/ level height) max-val)))
               (format t "|")
               (loop for timepoint in dynamics-data
                     do (let ((has-data nil))
                          (maphash (lambda (bacteria abundance)
                                    (when (>= abundance threshold)
                                      (setf has-data t)))
                                  (getf timepoint :abundances))
                          (format t "~A" (if has-data "*" " "))))
               (format t "~%")))))

;;; ============================================================
;;; Utilities for SLIME interaction
;;; ============================================================

(defun print-experiment-summary (experiment)
  "実験の概要を表示"
  (format t "~%=== LGIM Experiment Summary ===~%")
  (format t "Total Samples: ~A~%" 
          (length (microbiome-experiment-samples experiment)))
  (format t "Bacteria Species: ~A~%" 
          (length (microbiome-experiment-bacteria-names experiment)))
  (format t "Gravity Conditions: ~A~%" 
          (microbiome-experiment-gravity-conditions experiment))
  (format t "Time Points: ~A~%" 
          (microbiome-experiment-time-points experiment))
  (format t "Donors: ~A~%~%" 
          (microbiome-experiment-donors experiment)))

(defun print-responsive-species (responsive &optional (n 10))
  "重力応答性菌種をリスト表示"
  (format t "~%=== Top ~A Gravity-Responsive Species ===~%" n)
  (loop for i from 0 below (min n (length responsive))
        for species in responsive
        do (format t "~A. ~A (variance: ~,4F)~%"
                  (1+ i) (getf species :bacteria) (getf species :variance))))

;;; ============================================================
;;; Example Usage for SLIME
;;; ============================================================

(defun example-usage ()
  "SLIMEでの使用例"
  (format t "~%Example MDSINE2 Analysis~%")
  (format t "~%Loading data...~%")
  
  (let ((experiment (load-lgim-data "/Users/satoutsubasa/LGIM.csv")))
    (print-experiment-summary experiment)
    
    ;; 重力応答性菌種の検出
    (let ((responsive (gravity-responsive-species experiment :threshold 0.05)))
      (print-responsive-species responsive 15))
    
    ;; 特定菌種の統計比較
    (format t "~%=== Statistical Comparison: Bacteroides ===~%")
    (let ((comparison (statistical-comparison experiment "Bacteroides")))
      (loop for item in comparison
            do (let* ((gravity (getf item :gravity))
                     (result (getf item :result))
                     (mean (getf result :mean))
                     (sd (getf result :sd)))
                 (format t "~A: mean=~,4F, sd=~,4F~%" gravity mean sd))))
    
    ;; 特定条件の時系列動態
    (format t "~%=== Temporal Dynamics: Donor 1, 1g Gravity ===~%")
    (let ((dynamics (temporal-dynamics experiment "1g" 1)))
      (loop for tp in dynamics
            do (format t "~A (rep ~A): ~A species~%"
                      (getf tp :time-point)
                      (getf tp :replicate)
                      (hash-table-count (getf tp :abundances)))))))

;;; EOF
