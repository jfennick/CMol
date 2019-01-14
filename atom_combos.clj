;this file calculates and stores the one two and three center combinations of connected atoms.
(ns atom-combos (:use [geometry :only [magnitude get-com]]
		      [clojure.set :only [select]]))

(defn get-symbols [timestep]
  (reduce #(conj %1 %2) #{}
	  (map :species timestep)))

(defn get-species [timestep species]
  "Returns the set of all atoms of the given species in the current timestep."
  (select #(= species (:species %)) timestep))

(defn get-all-species [timestep symbols]
  "Returns a seq of sets of each species of atoms in the current timestep."
  (map #(get-species timestep %1) symbols))


(defn distance-pred [rmin rmax]
  "Returns a predicate of two args that returns true if the two coordinates are between rmin and rmax."
  #(let [dr (map - %1 %2) mag (magnitude dr)]
     (and (> mag rmin) (< mag rmax))))

(defn com-pred [obj1 obj2 rmin rmax]
  "Returns true if get-com of obj1 and obj2 are between rmin and rmax."
  ((distance-pred rmin rmax) (get-com obj1) (get-com obj2)))

(defn objs-pred [objs cutoffs]
  "Given n objs and n-1 cutoffs, returns true if every pair of objs satisfies com-pred.
   cutoffs are of the form [[mins] [maxes]]"
  (every? true? (apply map (fn [[a b] c d] (com-pred a b c d)) (partition 2 1 objs) cutoffs)))

;todo allow rmax regression
;(def all-bonds (map #(let [best (evolve (apply bond-fit-func %) [0.5 0] [3.0 20] 2 1000)
;			   al (get-alleles best)]
;		       (bondlength a b (first al) (+ (first al) (second al)) true))
;		    bond-combos))
