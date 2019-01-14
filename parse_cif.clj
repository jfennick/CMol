(ns parse-cif (:import java.net.URL)
    (:use [infix :only [infix-to-prefix]]
	  [utils :only [read-all-string transpose]]
	  [clojure.contrib.io :only [read-lines]]))

(def codurl "http://fireball-dft.org/cod/cod/")
(def codnames (read-lines (URL. (str codurl "files"))))
(def cifs (map #(URL. (str codurl (.substring #^String % 2))) codnames));skip initial "./"


(def regex #"['\"].*(\\'|\\\")*.*['\"]|\S+")

(defn split [#^String s]
  "Splits string by whitespace, excluding whitespace inside of quotes."
    (map first (re-seq regex s)))

(defn remove-comment [#^String s]
  "If s contains #, removes it and all remaining characters."
  (if (. s contains "#")
    (. s substring 0 (. s indexOf "#"))
    s))

(defn wrap [a b]
  "Like conj, but wraps a in a vector."
  (conj (if (vector? a) a [a]) b))

(defn no-blanks [lines]
  "Removes lines that are solely whitespace"
  (remove #(= "" #^String %) lines))

(defn comment? [#^String s]
  (. s startsWith ";"))

(defn remove-semi [#^String s]
  (.replaceAll s ";" ""))

(defn cif-name? [#^String s]
  (. s startsWith "_"))

(defn loop? [#^String s]
  (= s "loop_"))

(defn no-end-quote? [#^String line];this function working?
  (let [#^LString sp (. line split " ")
	   one (first sp)
	   two (. #^String (apply str (interpose " " (rest sp))) trim)]
    (and (cif-name? one) (. two startsWith "'") (not (. two endsWith "'")))))

(defn reduce-both
  "Like reduce, but f must return both a val and seq of remaining items."
  [f val coll]
  (loop [[acc rem] (f val coll)]
    (if (or (nil? rem) (empty? rem)) acc
	(recur (f acc rem)))))

(defn concat-multi [lines]
  (reduce-both (fn [xs [y & ys]]
		 (if-not (comment? y)
		   [(conj xs y) ys]
		   (if (comment? (first ys));special case for back-to-back comments
		     [(conj (vec (butlast xs)) (remove-semi (str (last xs) " '" y " " (first ys) "'"))) (next ys)]
		     (let [[a b] (split-with (complement comment?) ys)]
		       [(conj (vec (butlast xs)) (remove-semi (apply print-str (concat [(last xs) "'" y] a [(first b) "'"]))))
			(next b)]))))
	       [] lines))
		   
(defn concat-width [data width]
  "data is a seq of seqs. concats consecutive lines until line >= width"
  (->> (reduce-both (fn [xs [y & ys]]
		      (if (nil? ys);Infinite loop when out of data without this nil check
			[(conj xs y) ys]
			(if (>= (count y) width)
			  [(conj xs y) ys]
			  [xs (vec (cons (concat y (first ys)) (next ys)))])))
		    [] (map split data))
       (map #(apply print-str %))))

(defn fix-split [lines]
  "Concatenates data split across multiple lines, both standalone and in loop_ sections."
  (reduce-both (fn [xs [y & ys]]
		 (if-not (loop? y)
		   (let [[y2 & more] ys]
		     (if (and (cif-name? y)
			      (= 1 (count (split y)));line contains only a cif-name, and no value
			      (not (or (cif-name? y2) (loop? y2) (comment? y2))))
		       [(conj xs (str y " " y2)) more];concatenate subsequent lines
		       [(conj xs y) ys]));do nothing
		   (let [[name val] (split (first ys))]
		     (if val [xs ys];unnecessary loop_ section
			 (let [[a b] (split-with cif-name? ys)
			       [data rem] (split-with #(and (not (loop? %)) (not (cif-name? %))) b);loop_ cif-name [no-data] loop_ causes infinite loop
			       done (vec (concat (conj xs y) a (concat-width data (count a))))]
			   [done rem])))))
	       [] lines))

(defn fix-line [#^String line]
  (let [journal (re-matches #"_journal volume\s+(\S+)" line)]
    (cond (= line "_chemical_formula_structural  C31 H32 Br2 Ni O2 P2 x CH2Cl2")
	  "_chemical_formula_structural  'C31 H32 Br2 Ni O2 P2 x CH2Cl2'"
	  
	  (re-matches #"data_[0-9]+" line) ""
	  (re-matches #"loop_(\s+\S+)+" line) (seq (.split line "\\s+"))
	  journal (second journal)
	  (no-end-quote? line) (str line "'")
	  true line)))

(defn get-lines [url]
  "Removes comments, blank lines, concats multiline and split values, and fixes specific errors."
  (->> (read-lines url)
       (map #(. #^String % trim))
       (map remove-comment)
       (map fix-line)
       (reduce #(if (sequential? %2) (vec (concat %1 %2)) (conj %1 %2)) [])
       no-blanks
       concat-multi
       fix-split))

(defn add-to-hash [hash-ref new-hash]
  "Takes a ref to a hashmap and merges with new-hash"
  (dosync (ref-set hash-ref (merge-with wrap @hash-ref new-hash))))

(defn parse-cif [url]
  "Parses the given cif url into a hashmap of name item pairs."
  (let [hash (ref {})
	names (ref []) ;a vector of "recent" strings starting with _
	build (ref false)] ;true if the last line started with a _
    (dorun (for [line (get-lines url)]
	     (cond (nil? line) nil ;shouldn't be any nils, but make sure
		   
		   (= line "loop_") (dosync (ref-set build true)
					(ref-set names []))
		   
		   (cif-name? line)
		   (if @build
			 (dosync (alter names conj line));build up the list of cif-names
		         ;add single value to hashmap
		     (let [tmp1 (split line)
			   tmp2 (if (= (count tmp1) 1);if key, but not value specified, default value to nil
				  (concat tmp1 '(nil))
				  tmp1)]
		       (add-to-hash hash (apply hash-map tmp2))))

		   ;add multiple values to hashmap
		   (not (empty? @names))
		   (dosync (dorun (map #(ref-set hash (merge-with wrap @hash {%1 %2}))
				       @names (split line)))
			   (ref-set build false))
		   true (println (apply str url " Default:" line))
		   )))
    @hash))

(def num-keys ["_atom_site_fract_"
	       "_atom_site_aniso_U_"
	       "_atom_site_U_iso_or_equiv"
	       "_atom_site_symmetry_multiplicity"
	       "_atom_site_occupancy"
	       "_atom_site_attached_hydrogens"
	       "_cell_angle_"
	       "_cell_length_"
	       ])
(comment "_atom_site_Wyckoff_symbol" "_atom_site_disorder_assembly"
	 "_atom_site_refinement_flags" "_atom_site_disorder_group")

(defn get-num-keys [cif-map]
  "Returns all the numeric keys in cif-map"
  (filter (fn [#^String key] (some (fn [#^String x] (. key contains x))
				   num-keys))
	  (keys cif-map)))

(defn get-nums [cif-map]
  (let [nums (get-num-keys cif-map)]
    (zipmap nums (map cif-map nums))))

(defn parse-numerics [cif-map]
  "Processes numeric keys (coordinates, etc), replacing strings with parsed numbers."
  (let [parser #(-> (condp = %, "." "0", "?" "Infinity", %)
		    (.replaceAll "\\(" "")
		    (.replaceAll "\\)" "")
		    Float.)
;remove parens around uncertain trailing digit(s)
;replace unknowns with Infinity (equally uninformative, but still parses)
;use Float. instead of read-string because sometimes numbers do not have leading zero's.
	nks (get-num-keys cif-map)]
    (->> nks
	 (map #(let [val (get cif-map %)]
		 (if (sequential? val)
		   (map parser val);use array to save memory?
		   (parser val))))
	 (interleave nks)
	 (apply hash-map)
	 (merge cif-map))))

(defn get-fractional [cif-map]
  "Returns a list of the fractional coordinates"
  (let [s "_atom_site_fract_"]
    (transpose (map #(get cif-map (str s %))
		    ["x" "y" "z"]))))

(defmacro capture-vars [vars expr]
  "This intentionally captures vars, so we can include vars as part of expr."
  `(fn [~@vars] ~expr))

(defn get-op-fn [#^String op-str]
  (->>  (if (re-matches #"[+-][xyz]" op-str) (str "0" op-str) op-str)
	(interpose " ")
	(apply str)
	read-all-string
	infix-to-prefix
	((fn [q] `(capture-vars ~'[x y z] ~q)))
	eval))

(defn get-symmetry-ops [cif-map]
  "Returns a sequence of fn 3-tuples representing the symmetry operations."
  (let [f #(-> (.replaceAll % "'" "") (.split ","))]
    (->> (get cif-map "_symmetry_equiv_pos_as_xyz")
	 (map f)
	 (map #(map get-op-fn %)))))

(defn get-frac-coords [cif-map]
  (for [op (get-symmetry-ops cif-map)
	coord (get-fractional cif-map)]
    (map #(apply % coord) op)))

(defn parse-cifs []
  (remove empty? (map #(try (-> % parse-cif parse-numerics get-frac-coords)
			  (catch Exception e
			    (println (str "Error: " %))))
		     cifs)))