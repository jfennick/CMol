(ns test.test-match-network (:use atom-combos geometry gview match-network xyz-parser
				  clojure.set clojure.test clojure.contrib.seq))

(defn chain-networks [& args])

(defn tih2-n [n hmin hmax tmin tmax]
  "Given H-H and Ti-H bounds, matches a Ti atom with n unique H-H pairs."
  (let [ti (parse-xyz "test/5510TiH90.xyz")
	h2 {[1 2] {:names ["H" "H"] :rmsc [hmin hmax]}}
	tih2n (into {} (for [i (range 2 (+ 2 n))]
			 [[1 i] {:names ["Ti" h2] :rmsc [tmin tmax]}]))]
    (match-networks (make-db (set (first ti))) tih2n {})))
(def tih2-match (tih2-n 3 0.8 1.0 1.6 2.0))

(def C6 (make-chain (repeat 6 "C") (repeat 6 0) (repeat 6 1.5) true))
(def C4N2 (make-chain ["C" "C" "C" "N" "C" "N"]
		      (repeat 6 0) (repeat 6 1.5) true))
(def C3N2 (make-chain ["C" "N" "C" "N" "C"]
		      (repeat 5 0) (repeat 5 1.5) true))
(def NH2 (make-branch ["N" "H" "H"] [0 0] [1.2 1.2]))
(def CH3 (make-branch ["C" "H" "H" "H"] [0 0 0] [1.2 1.2 1.2]))

(defn with-name [num nums name names rmaxes hmaps]
  (zipmap (map #(vector num %) nums)
	  (map #(conj %3 {:names [name %1] :rmsc [0 %2]}) names rmaxes hmaps)))

(def adenine {[1 2] {:names [C4N2 C3N2] :rmsc [0 10] :eq1 [[1] [2]] :eq2 [[1] [5]] :plot false}
	      [1 3] {:names [C4N2 "H"] :rmsc [0 1.2] :a1 [5]}
	      [1 4] {:names [C4N2 NH2] :rmsc [0 1.6] :a1 [3] :a2 [1]}
	      [2 5] {:names [C3N2 "H"] :rmsc [0 2.2] :a1 [2]}
	      [2 6] {:names [C3N2 "H"] :rmsc [0 2.2] :a1 [3]}})

(def guanine {[1 2] {:names [C4N2 C3N2] :rmsc [0 10] :eq1 [[1] [2]] :eq2 [[1] [5]] :plot false}
	      [1 3] {:names [C4N2 "H"] :rmsc [0 1.2] :a1 [4]}
	      [1 4] {:names [C4N2 NH2] :rmsc [0 1.6] :a1 [5] :a2 [1]}
	      [1 7] {:names [C4N2 "O"] :rmsc [0 1.6] :a1 [3]}
	      [2 5] {:names [C3N2 "H"] :rmsc [0 2.2] :a1 [2]}
	      [2 6] {:names [C3N2 "H"] :rmsc [0 2.2] :a1 [3]}})

(def cytosine{[1 2] {:names [C4N2 "H"] :rmsc [0 1.2] :a1 [1]}
	      [1 3] {:names [C4N2 "H"] :rmsc [0 1.2] :a1 [2]}
	      [1 4] {:names [C4N2 "H"] :rmsc [0 2.2] :a1 [6]}
	      [1 5] {:names [C4N2 "O"] :rmsc [0 1.6] :a1 [5]}
	      [1 6] {:names [C4N2 NH2] :rmsc [0 1.6] :a1 [3] :a2 [1]}})

(def thymine {[1 2] {:names [C4N2 "H"] :rmsc [0 1.2] :a1 [1]}
	      [1 3] {:names [C4N2 "H"] :rmsc [0 1.2] :a1 [4]}
	      [1 4] {:names [C4N2 "H"] :rmsc [0 2.2] :a1 [6]}
	      [1 5] {:names [C4N2 "O"] :rmsc [0 1.6] :a1 [3]}
	      [1 6] {:names [C4N2 "O"] :rmsc [0 1.6] :a1 [5]}
	      [1 7] {:names [C4N2 CH3] :rmsc [0 1.6] :a1 [2] :a2 [1]}})

(def base-pairs #{adenine guanine cytosine thymine})
(def any-dna (make-chain (repeat 30 base-pairs)
			 (repeat 30 3.0) (repeat 30 4.5) false))
;base-pairs is a wildcard.  It will match any dna base pair.

(def bases-map {\A adenine \G guanine \C cytosine \T thymine})
(defn make-dna-def [string]
  (make-chain (map bases-map string) (repeat 3.0) (repeat 4.5) false))

(def A (parse-xyz "test/dna.xyz"))
(def adb (make-db (set (first A))))
(def dups {:dups remove-dup-matches :dupset #{base-pairs adenine guanine cytosine thymine}})
(def any-dna-match (time (match-networks adb any-dna dups)))

(def pair-map (zipmap (map any-dna-match [adenine guanine cytosine thymine]) ["A" "G" "C" "T"]))

(defn sequence-dna [match]
  (apply str (map #(first (remove nil? (map (fn [[s n]] (if (s %) n)) pair-map)))
		  (vals match))))

;(def Asorted (sort #(< (magnitude (get-com %1)) (magnitude (get-com %2))) Amatch))
;(def Acoords (ref (map :coordinates (first A))))
;(tree-view Asorted (let [view (viewer-3d (plot-coords Acoords))] (tsg (second view) Acoords)))

;(def Tsorted (sort #(< (magnitude (get-com %1)) (magnitude (get-com %2))) Tmatch))

(def backbone {[1 2] {:names ["P" "O"] :rmsc [0 1.5]}
	       [1 3] {:names ["P" "O"] :rmsc [0 1.5]}
	       [1 4] {:names ["P" "O"] :rmsc [0 1.6]}
	       [4 5] {:names ["O" "C"] :rmsc [0 1.5]}
	       [5 6] {:names ["C" "H"] :rmsc [0 1.1]}
	       [5 7] {:names ["C" "H"] :rmsc [0 1.1]}
	       [5 8] {:names ["C" "C"] :rmsc [0 1.6]}
	       [8 9] {:names ["C" "H"] :rmsc [0 1.1]}
	       [8 10] {:names ["C" "O"] :rmsc [0 1.5]}
	       [10 11] {:names ["O" "C"] :rmsc [0 1.5]}
	       [11 12] {:names ["C" "H"] :rmsc [0 1.1]}
	       [11 13] {:names ["C" "C"] :rmsc [0 1.7]}
	       [13 14] {:names ["C" "H"] :rmsc [0 1.1]}
	       [13 15] {:names ["C" "H"] :rmsc [0 1.1]}
	       [13 16] {:names ["C" "C"] :rmsc [0 1.7]}
	       [16 17] {:names ["C" "H"] :rmsc [0 1.1]}
	       [16 18] {:names ["C" "O"] :rmsc [0 1.5]}
	       [8  16] {:names ["C" "C"] :rmsc [0 1.6]}})

;(def Bmatch (time (match-networks adb backbone {})))
;(def Bsorted (sort #(< (magnitude (get-com %1)) (magnitude (get-com %2))) Bmatch))

(deftest match-network
  (is (= [15 15 15 15 120 60 90 90 4]
	 (map #(count (any-dna-match %))
	      [adenine guanine cytosine thymine C4N2 C3N2 CH3 NH2 any-dna])))
  (is (= ["CTTAGCCGATCTTAGCCGATCTTAGCCGAT"
	  "GAATCGGCTAGAATCGGCTAGAATCGGCTA"
	  "ATCGGCTAAGATCGGCTAAGATCGGCTAAG"
	  "TAGCCGATTCTAGCCGATTCTAGCCGATTC"]
	 (map sequence-dna (any-dna-match any-dna))))
  (is (= 0 (count (tih2-match "TiH2-3")))))