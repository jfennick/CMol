(ns ga (:use [clojure.contrib.seq :only [reductions]]))

(defn find-index
  "Returns the (0-based) index of sorted coll where val < coll[index]"
  [val coll]
  (dec (count (take-while #(<= % val) coll))));todo: binary search may be faster

(defn tournament-select
  [pop fit-fn k p]
  (let [indexed (zipmap pop (map fit-fn pop))]
    ))

(defn roulette-select
  "Select num individuals from pop, with an individual's selection likelihood proportional to fitness."
  [pop fit-fn num]
  (println "roulette")
  (let [pop-fits (map fit-fn pop)
	fit-min (apply min pop-fits), fit-max (apply max pop-fits), fit-range (- fit-max fit-min)
	rescaled-fits (map #(- (- % fit-max)) pop-fits)
	inc-fits (reductions + 0 rescaled-fits)]
    (take num (repeatedly #(let [val (rand fit-range)
				 index (find-index val inc-fits)
				 max-index (dec (count pop))]
			     (if (or (> index max-index) (neg? index))
			       (throw (RuntimeException. (str "index out of bounds:" index ":" max-index ":" val ":" (vec inc-fits)))))
			     (nth pop (min index max-index)))))))

(defn remove-atmost [n f coll]
  (loop [num-removed 0, kept [], left coll]
    (if (or (> num-removed n) (empty? left))
      (concat kept left)
      (if (f (first left))
	(recur (inc num-removed) kept (rest left))
	(recur num-removed (conj kept (first left)) (rest left))))))

(defn map-meta
  "Maps f across coll, but keeps the original elements of coll as metadata."
  [f & colls]
  (apply map (fn [& args] (vary-meta (apply f args) assoc :args args)) colls))

(defn my-contains? [coll element equal-fn]
  (every? true? (map #(equal-fn element %) coll)))

(defn lazy-finite-approx-set
  "Use this function to return a seq of size n of approximately
   equal elements from an initially infinite stream of elements."
  [n equal-fn infinite-seq]
  (loop [coll [], [element & remaining] infinite-seq]
    (if (or (= n (count coll)) (empty? remaining))
      coll
      (if (my-contains? coll element equal-fn)
	(recur coll remaining)
	(recur (conj coll element) remaining)))))

(defn lazy-finite-set [n infinite-seq]
  (->> infinite-seq distinct (take n) set))

(defn ga
  "Returns and infinite seq of populations.
init-fn: takes no arguments and returns an infinite seq new population members
fit-fn: takes a population member and outputs 'annotated' member
fit-keyfn: take an 'annotated' member and outputs fitness
mut-fn: takes a population member and returns a mutated member.
sel-fn: takes a population, a fitness function, and a number to select. Returns selected members.
cross-fn: combines two or more parents into one child
equal-fn is a key-fn for population members.  It is used to ensure population diversity.
pop-size, elite-size, mut-prob, and cross-prob should be self-explanatory."
  [hmap]
  (let [{:keys [init-fn fit-fn fit-keyfn mut-fn lat-fn sel-fn cross-fn equal-fn
		mut-prob lat-prob cross-prob cross-num pop-size elite-size init-pop]
	 :or {pop-size 20, init-pop (take pop-size (init-fn)), mut-prob 0.01, cross-prob 0.9, lat-prob 0.01
	      cross-num 2, elite-size 1, equal-fn =, sel-fn roulette-select, lat-fn identity}} hmap]
    (iterate
     (fn [pop](println "iterating population")
       (let [temp-pop (take pop-size (concat (remove nil? (map fit-fn pop)) (init-fn)))
	     elites (lazy-finite-approx-set elite-size equal-fn (sort-by fit-keyfn < temp-pop))
	     ;TODO: add opposite of elite: unconditionally reject ~40% of worst structures
	     new-fn #(let [lol (println "new-fn")
			   temp1 (if (> cross-prob (rand))
				  (let [selected (sel-fn temp-pop fit-keyfn cross-num)]
				    (apply cross-fn selected))
				  (first (sel-fn temp-pop fit-keyfn 1)))
			   temp2 (if (> mut-prob (rand))
				   (mut-fn temp1) temp1)
			   temp3 (if (> lat-prob (rand))
				   (lat-fn temp2) temp2)
			   ]temp3)
	     chosen (->> (repeatedly new-fn)
;if population is not very diverse, equal-fn may infinitely remove crossover'd/mutated members.
;Currently, detect and give up. Generate brand-new members with some probability?
			 (remove-atmost pop-size #(my-contains? elites % equal-fn));pop-size is arbitrary in this case
			 (take (- pop-size (count elites))))]
	 (concat elites chosen))) init-pop)))
;Each element of iterate is itself lazy, so you need (dorun n (map dorun (ga ...))) to actually force the results.