(ns match-network (:use [atom-combos :only [com-pred distance-pred get-all-species get-symbols]]
			[bsptree :only [add-to-tree get-bounds make-bspnode]]
			[geometry :only [normalize-mag get-moi]]
			[jgap-wrapper :only [evolve get-alleles]]
			[utils :only [unflatten map-to-seq mat-mat-mult mat-vect-mult transpose]]
			[tree-utils :only [tree-apply]]
			[xyz-parser :only [parse-xyz]]
			[clojure.set :only [intersection]]
			[clojuratica-config :only [eigenvectors invert-matrix]]
			[clojure.contrib.generic.functor :only [fmap]]))

(defn make-db [timestep]
  "Creates a new database, and initializes with species"
  (let [symbols (get-symbols timestep)
	species (get-all-species timestep symbols)]
    (zipmap symbols species)))
;todo: allow for multiple timesteps

(defn make-bspdb [timestep]
  "Creates a new bsp database, and initializes with species"
  (let [[a b] (get-bounds timestep get-com)
	db (make-bspnode nil a b 10)
	lol (add-to-tree db timestep get-com)]
    db))

(defn get-names [network]
  "Returns a hashmap of [number #{set of names}]"
  (->> (mapcat #(let [[k v] %, [i1 i2] k, [n1 n2] (:names v)]
		  [{i1 [n1]} {i2 [n2]}])
	       network)
       (apply merge-with concat)
       (#(zipmap (keys %) (map set (vals %))))))

(defn check-names [network]
  "Returns an empty coll if each name corresponds to a unique number."
  (filter #(not= 1 (count (second %)))
	  (get-names network)))

;todo: check-lengths

(defn do-meta [match]
  ;(println (str "match1 " match))
  (vary-meta
    (disj (set (map (fn [[k v]]
		      (if-not (or (symbol? v) (coll? v))
			[k v] ;only add metadata to supported types
			(vary-meta (if (and (map? v) (not (:coordinates v)))
				     (do-meta v) v)
			  assoc :key k))) match))
	  nil);remove nil
    assoc :coll (empty match)))

(defn undo-meta [match]
  ;(println (str "match2 " match))
  (into (:coll (meta match))
	(map #(let [k (:key (meta %))]
		(vector (if k k (first %))
			(if k
			  (vary-meta (if (:coll (meta %))
				       (undo-meta %) %)
				     dissoc :key)
			  (second %)))) match)))

(defn remove-dup-matches [matches]
  "This is useful for networks where order doesn't matter (i.e. carbon rings)."
  (map undo-meta (set (map do-meta matches))))

(defn atom-unique [match bond]
  "Add empty? to form predicate to test if the new trial match is atomistically unique"
  (intersection (into #{} (flatten (map-to-seq match)))
		(if-not (:coordinates bond);if atom
		  (into #{} (flatten (map-to-seq bond)))
		  #{bond})))

(defn equal-nths [e1 hmap]
  "like atom-unique, but also tests if match and bond have equal nths."
  (fn [match bond]
    (let [eq1 (:eq1 hmap) eq2 (:eq2 hmap)
	  unique (atom-unique match bond)
	  inter (map #(set [(get-in e1 %1) (get-in bond %2)]) eq1 eq2)
	  one (every? #(= 1 (count %)) inter)
	  two (= unique (apply union inter))]
      ;(println (str "one " one " two " two " uni " (count unique) " int " (count inter)))
      (and one two))))

(defn make-chain [names mins maxes ring]
  "names is a seq of strings"
  (let [chain (into {} (map #(vector (vec %1) {:names %2 :rmsc [%3 %4]})
			    (partition 2 1 (iterate inc 1)) (partition 2 1 names) mins maxes))
	end (if ring
	      {[1 (count names)] {:names [(first names) (last names)]
				  :rmsc [(last mins) (last maxes)]}})]
    (conj chain end)))

(defn make-branch [[name & names] mins maxes]
  "names is a seq of strings"
  (into {} (map #(vector [1 %1] {:names [name %2] :rmsc [%3 %4]})
		(iterate inc 2) names mins maxes)))

(defn get-at [a item]
  (get-com (if a (get-in item a) item)))

(defn get-nums [network]
  (set (apply concat (keys network))))

(defn check-nums [network]
  (let [nums (get-nums network)
	r (set (range 1 (inc (apply max nums))))]
    (= nums r)))

(defn get-from-db [stdb name]
  (cond (set? name) (mapcat stdb (filter name (keys stdb)))
	(instance? java.util.regex.Pattern name)
	(mapcat stdb (filter #(re-matches name %)
			     (filter string? (keys stdb))))
	true (get stdb name)))

(defn get-match-schema [match]
  (zipmap (keys match)
	  (map #(if (:coordinates %) (:species %)
		    (get-match-schema %))
	       (vals match))))
       
(defn get-query-schema [query]
  (let [match (get-names query)]
    (zipmap (keys match)
	    (map #(let [o (first %)]
		    (if (string? o) o
			(get-query-schema o)))
		 (vals match)))))

(defn get-wildcard [db match]
  "Returns the key in db which is the definition of match."
  (let [ms (get-match-schema match)]
    (first (filter #(= ms (get-query-schema %))
		   (remove string? (keys db))))))

(defn update-mol [mol f]
  "Maps f across the :coordinates of mol recursively, returning an 'updated' structure."
  (fmap #(if-not (:coordinates %)
	   (update-mol (if-let [o (:obj %)] o %) f)
	   (update-in % [:coordinates] f)) mol))

(defn move-mol [mol coord]
  "Recursively updates mol's :coordinates such that its new center of mass is equal to coord."
  (let [shift (map - (get-com mol) coord)]
    (update-mol mol #(map - % shift))))

(defn shift-mol [mol coord]
  "Recursively translates mol's :coordinates by coord."
  (update-mol mol #(map + % coord)))

(defn get-pai [mol]
  "Returns the principal axes of inertia of mol."
  (if (:coordinates mol) [[1 0 0] [0 1 0] [0 0 1]]
      (->> (get-moi mol) (eigenvectors user/math) transpose)))

(defn get-dpai [a b]
  "Returns the 'delta' between the principal axes of inertia,
   i.e. the rotation matrix necessary to rotate a to align with b"
  (let [apai (get-pai a) bpai (get-pai b)
	inv  (invert-matrix user/math apai)]
    (mat-mat-mult bpai inv)))

(defn align-mol [mol moi]
  "Recursively updates mol's :coordinates such that the principal axes of inertia of mol and moi are equal."
  (let [inv (->> (get-pai mol) (invert-matrix user/math))
	mat (->> (eigenvectors user/math moi) transpose)
	m (mat-mat-mult mat inv)]
    (update-mol mol #(mat-vect-mult m %))))

(defn rotate-mol [mol moi]
  "Recursively multiplies mol's :coordinates by the rotation matrix moi."
  (update-mol mol #(mat-vect-mult moi %)))

(defn rotate-mol-about [mol moi com]
  (-> mol
      (shift-mol (map - com))
      (rotate-mol moi)
      (shift-mol com)))

(defn spin-mol [mol moi]
  "Like rotate-mol, but moves mol to/from its center of mass, thus mol will rotate 'in-place'."
  (-> mol
      (move-mol [0 0 0])
      (rotate-mol moi)
      (move-mol (get-com mol))))

(declare match-networks)
(defn add-to-db [stdb name options]
;todo: figure out if this will work with lazyiness/multithreading.  i.e., need to use a lock/etc. to avoid recomputation?
  (if (contains? stdb name) stdb
      (if (set? name)
	(reduce #(add-to-db %1 %2 options) stdb name)
	(match-networks stdb name options))))

(defn do-dups [n2 options matches]
  ((if ((:dupset options identity) n2)
     (:dups options identity) identity) matches))

(defn rot-mat-x [t]
  (let [s (Math/sin t) c (Math/cos t)]
    [[1 0 0] [0 c s] [0 (- s) c]]))

(defn rot-mat-y [t]
  (let [s (Math/sin t) c (Math/cos t)]
    [[c 0 (- s)] [0 1 0] [s 0 c]]))

(defn rot-mat-z [t]
  (let [s (Math/sin t) c (Math/cos t)]
    [[c s 0] [(- s) c 0] [0 0 1]]))

(defn rotate-nthrest [mol mat n]
  "Multiplies the nth component and all remaining components of mol, returning a new mol."
  (let [vm (vals mol), obj (nth vm n), com (get-com obj)]
    (zipmap (keys mol)
	    (concat (take n vm)
		    (map #(rotate-mol-about % mat com) (drop n vm))))))

(defn save-chain [match]
  (let [mov-rot (fn [dcom dpai] #(-> % (shift-mol (map - dcom)) (rotate-mol dpai)))
	one (first (vals match))
	com (get-com one) pai (get-pai one)]
    (loop [left (map (mov-rot com (invert-matrix user/math pai)) (vals match))
	   trans [] result []]
      (if (= 1 (count left))
	[(zipmap (keys match) (concat result left))
	 (cons [com pai] trans)]
	(let [[a b & more] left
	      dcom (map - (get-com b) (get-com a))
	      dmoi (get-pai b)]
	  (recur (map (mov-rot dcom (invert-matrix user/math dmoi)) (next left))
		 (conj trans [dcom dmoi])
		 (conj result a)))))))

(defn build-chain [[match trans]]
  (let [rot-mov (fn [dcom dpai] #(-> % (rotate-mol dpai) (shift-mol dcom)))
	[com pai] (first trans)
	builder (fn builder [mol tr]
		  (let [f (apply rot-mov (first tr))]
		    (if (= 1 (count mol))
		      (map f mol)
		      (map f (cons (first mol) (builder (next mol) (next tr)))))))]
    (->> (builder (next (vals match)) (next trans))
	 (cons (first (vals match)))
	 (map (apply rot-mov (first trans)))
	 (zipmap (keys match)))))

(defn repeat-com1 [[match t]]
  (let [[coms mois] (transpose t)
	[c & cs] coms
	[m & ms] mois]
    [match (transpose [(cons c (repeat 29 (first cs))) mois])]))

(defn save-network [match]
  (zipmap (keys match)
	  (map #(if-let [c (:coordinates %)] {:obj %}
			(let [temp (move-mol % [0 0 0])
			      mat (get-pai temp)
			      inv (->> mat (invert-matrix user/math) (map normalize-mag))
			      rotated (update-mol temp (fn [coord] (mat-vect-mult inv coord)))]
			  {:dcom (get-com %) :dmoi mat :obj (save-network rotated)}))
	       (vals match))))

(defn build-network [match]
  (zipmap (keys match)
	  (map #(if-let [c (:coordinates (:obj %))] (:obj %)
			(let [net (build-network (:obj %))
			      mat (map normalize-mag (:dmoi %))
			      rotated (update-mol net (fn [coord] (mat-vect-mult mat coord)))]
			  (move-mol rotated (:dcom %))))
	       (vals match))))

(defn match-networks [stdb network options]
  "Returns a set of hmaps containing [num obj] that match the given network."
  (let [names (check-names network)]
    (if-not (empty? names)
      (println (str "One or more names have the same number: " (print-str names)))
      (let [sorted (apply sorted-map (apply concat network))
	    initval (second (first sorted))
	    nums (get-nums sorted)
	    cs (apply create-struct (sort nums))]
	(loop [db (add-to-db stdb (first (:names initval)) options)
	       matches (map #(hash-map (apply min nums) %);struct-map ~4% slower.
			    (get-from-db db (first (:names initval))))
	       hmap (seq sorted)]
	  (println (str (count matches) " " (count (keys db))));(first hmap)))
	  ;(println (str "keys  " (keys (first matches))))
	  ;(println (str "keys= " (apply = (map keys matches))))
	  (if (or (nil? hmap) (empty? matches))
	    (conj db {network (set matches)});return new database with network added
	    (let [[k v] (first hmap)
		  [i1 i2] k, [n1 n2] (:names v)
		  [rmin rmax] (:rmsc v)
		  e2 (get (first matches) i2)
		  dbtemp (add-to-db db n2 options)
		  result (if e2;if both atoms already exist in match, ensure distance constraints are met
			   (filter #(com-pred (% i1) (% i2) rmin rmax) matches)
			   (->> (for [match matches]
				  (let [e1 (get match i1)
					p1 (get-at (:a1 v) e1)
					bonds (if-let [dist (:dist-fn options)]
						(dist rmin rmax p1 n2);dist-fn must return all atoms of type n2
						;else brute-force neighbor search
						(filter #((distance-pred rmin rmax) p1 (get-at (:a2 v) %))
							(get-from-db dbtemp n2)))
					options2 (if (and (:eq1 v) (:eq2 v)) (conj options {:constraint-fn (equal-nths e1 v)}) options)]
				    (for [bond bonds :when ((:constraint-fn options2 (comp empty? atom-unique)) match bond)]
				      ((:combiner-fn options #(conj %1 {i2 %2})) match bond))))
				(apply concat) (do-dups network options)))]
	      (recur dbtemp result (next hmap)))))))))

(defn maxes-fit-func [stdb network options]
  (proxy [org.jgap.FitnessFunction] []
    (evaluate [chromo]
	      (let [al (get-alleles chromo)
		    [fmins fmaxes] (apply map list (partition 2 al))
		    new-network nil;todo: (unflatten mins fmins) (unflatten maxes fmaxes)
		    matches (match-networks stdb new-network options)]
		(count (get matches new-network))))))

(defn fit-maxes [stdb options names mins maxes size maxiter]
  "Finds maxes that maximize the number of matches in chain-networks"
  (let [mymaxes (if (sequential? maxes) maxes
		    (unflatten mins (repeat maxes)))
	mymins (if (not (sequential? maxes)) mins
		   (tree-apply #(* 0.8 %) false false maxes))
	fit-func (maxes-fit-func stdb options names mins maxes)
	result (evolve fit-func (flatten mymins) (flatten mymaxes) size maxiter)
	al (get-alleles result)
	[fmins fmaxes] (apply map list (partition 2 al))]
    (unflatten maxes fmaxes)))

(defn match-angles [bonds anglesmin anglesmax]
  "bonds is the output of match-bonds, and angles are seqs of two less elements than each tuple of bonds."
  );todo
