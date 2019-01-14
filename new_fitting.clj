(ns new-fitting (:use atom-combos geometry jgap-wrapper utils match-network potentials fortran-parser
		      clojure.set (clojure.contrib seq combinatorics)))

;"keys are atom/set targets, can be any regular expression"
(def potentials
     '{[".*" ".*"] {:name harmonic :params [1.5 -0.01] :rmsc [[0.5] [3.0]]}
       ["C" "C"] {:name harmonic :params [1.4 -0.002] :rmsc [[0.5] [2.0]]}
       ["C" "H"] {:name harmonic :params [1.2 -0.003] :rmsc [[0.5] [1.5]]}
       ["H" "C"] {:name harmonic :params [1.2 -0.003] :rmsc [[0.5] [1.5]]}
       ["H" "H"] {:name harmonic :params [0.8 -0.00465] :rmsc [[0.5] [1.2]]}
       :offset 0
       })

(defn fill-wildcards [potential structure-database]
  "Returns a map of all matches of potential's target in structure-database, with wildcards replaced with real targets."
  (let [[pkey pval] potential
	potkeys (if (= :offset pkey) [pkey]
		    (filter #(reduce (fn [x y] (and x y));if all match individually
				     (map (fn [a b] (re-matches (re-pattern a) b)) pkey %))
			    (selections (keys structure-database) (count pkey))))]
    (zipmap potkeys (repeat pval))))

(defn fill-potentials [potentials structure-database]
  "Fills in all wildcards in all potentials, returning a seq of 'raw' potentials"
  (reduce merge (map #(fill-wildcards % structure-database) potentials)))

;(def raw-potentials (fill-potentials potentials @stdb))

(defn sum-potential [potential stdb]
  "Retrieves all sets of objs (satisfying radius and angle bounds) from structure-database, and returns the sum of the potential's value on each object."
  (let [[pkey pval] potential]
    (if (= :offset pkey) pval
    ;todo replace cartesian-product with bsptree code
	(let [allatms (apply cartesian-product (map stdb pkey))
	      cutoffs (:rmsc pval)
	      atms (if-not cutoffs allatms (filter #(objs-pred % cutoffs) allatms))
	      ]
	  (reduce + (map #((resolve (:name pval)) (:params pval) %)
			 atms))))))

(defn all-potentials [raw-potentials structure-database]
  "Calculates the total energy of the system."
  ;(let [raw-potentials (fill-potentials potentials structure-database)]
  (reduce + (map #(sum-potential % structure-database) raw-potentials)))

(defn get-params [raw-potentials]
;todo: add flag to return/not return params/cutoffs
  "Returns a seq of all the parameters where each item is [[params] [cutoffs]]"
  (concat (map #(vector (:params %) (get % :rmsc []))
	       (vals (dissoc raw-potentials :offset)))
	  [(get raw-potentials :offset)]))

(defn update-params [params raw-potentials]
  "updates the values of the parameters in raw-potentials, returning new raw-potentials."
  (into {:offset (last params)}
	(map #(let [[param cutoff] %1
		    [pkey pval] %2]
		[pkey (merge pval {:params param} {:rmsc cutoff})])
	     (unflatten (get-params raw-potentials) params)
	     (dissoc raw-potentials :offset))))

(defn relative-error [actual approx]
  (/ (Math/abs (- approx actual)) (Math/abs actual)))

;(nlcg #(Math/pow (- (first energies) (all-potentials (update-params % raw-potentials) structure-database)) 2) (flatten (get-params raw-potentials)) 0.01 10 0.01)

(defn bonds-fit-func [funcs widths energy]
  (proxy [org.jgap.FitnessFunction] []
    (evaluate [chromo]
	      (let [al (get-alleles chromo)
		    args (partition-width al widths)
		    result (reduce #'+ (map #(apply %1 %2) funcs args))]
		(/ (relative-error energy result))))))

;(def numbonds (count all-bonds))
;(def energies (:cohesive (parse-fortran "/home/jacob/5510TiH90/f.log" fireball-parsers)))
;(def harmonic-params
;     (evolve (bonds-fit-func (harmonics all-bonds) (repeat numbonds 2) (first energies))
;		(flatten (repeat numbonds [0.7 -1000]))
;		(flatten (repeat numbonds [3 1000])) 10 10))
;(get-alleles harmonic-params);returns a seq of Ro and K values
