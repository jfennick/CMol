(ns infix);copied from clojure google group

(def +precedence+ {'rem 5, '* 4, '/ 3, '+ 2, '- 1})

;; highest level of precedence
(def +highest-precedence+ (apply max (map val +precedence+)))

(defn operator?
  "Check if is valid operator"
  ([sym]
     (not (nil? (get +precedence+ sym)))))

(defn find-lowest-precedence
  "Find the operator with lowest precedence; search from left to right"
  ([seq]
     ;; loop through terms in the sequence
     (loop [idx 0
            seq seq
            lowest-idx nil
            lowest-prec +highest-precedence+]
       ;; nothing left to process
       (if (empty? seq)
         ;; return lowest found
         lowest-idx
         ;; otherwise check if current term is lower
         (let [prec (get +precedence+ (first seq))]
           ;; is of lower or equal precedence
           (if (and prec (<= prec lowest-prec))
             (recur (inc idx) (rest seq)
                    idx prec)
             ;; is of high precedence therefore skip for now
             (recur (inc idx) (rest seq)
                    lowest-idx lowest-prec)))))))

(defn infix-to-prefix
  "Convert from infix notation to prefix notation"
  ([seq]
     (cond
      ;; handle term only
      (not (sequential? seq)) seq
      ;; handle sequence containing one term (i.e. handle parens)
      (= (count seq) 1) (infix-to-prefix (first seq))
      ;; handle all other cases
      true (let [lowest (find-lowest-precedence seq)]
             (if (nil? lowest) ;; nothing to split
               seq
               ;; (a b c) bind a to hd, c to tl, and b to op
               (let [[hd tl] (split-at lowest seq)
                     op (first tl)
                     tl (rest tl)]
                 ;; recurse
                 (list op (infix-to-prefix hd) (infix-to-prefix tl))))))))
                 
(defmacro formula
  "Formula macro translates from infix to prefix"
  ([& equation]
     (infix-to-prefix equation)))