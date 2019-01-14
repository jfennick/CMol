(ns fitting (:use testing jgap-wrapper utils))

(defn count-params [fn]
  (count (first (:arglists (meta fn)))))

(defn total-params [fns]
  (reduce #'+ (map #'count-params fns)))

(defn harmonic [set Ro K]
  (reduce #'+ (map #(let [c1 (get-com (. %1 getElement 0))
			  c2 (get-com (. %1 getElement 1))
			  dr (get-dr c1 c2)
			  rmag (. dr magnitude)]
		      (* 0.5 K (Math/pow (- rmag Ro) 2)))
		   (get-private-field set "set"))))

(defn harmonics []
  (map #(partial #'harmonic %1) all-bonds))

(defn relative-error [actual approx]
  (/ (Math/abs (- approx actual)) (Math/abs actual)))

(defn bonds-fit-func [funcs widths energy]
  (proxy [org.jgap.FitnessFunction] []
    (evaluate [chromo]
	      (let [al [1 10 1 10 1 10 1 10];(get-alleles chromo)
		    args (partition-width al widths)
		    result (reduce #'+ (map #(apply %1 %2) funcs args))]
		(/ (relative-error energy result))))))

(def numbonds (count all-bonds))
(def energies (map #(nth (. %1 y) 0) (. wph parseCohesiveEnergy)))
(def harmonic-params
     (evolve (bonds-fit-func (harmonics) (replicate numbonds 2) (first energies))
		(flatten (replicate numbonds [0.7 -1000]))
		(flatten (replicate numbonds [3 1000])) 10 10))
