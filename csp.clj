(ns csp (:use
         [atomic-structure-output :only [write-bas]]
         [geometry :only [align-lattice cross-product dot-product get-com invert-matrix3 magnitude mat-mat-mult
			       mat-vect-mult normalize-mag parameters->vectors vectors->parameters point-plane-distance]]
	      [supercell :only [supercell-atoms]]
	      [fireball-wrapper :only [cgopt kpts normalize-weights pressure required-files run-fireball]]
	      [gulp-wrapper :only [run-gulp]]
	      [symmetry-wrapper :only [get-symmetry]]
	      [fortran-parser :only [parse-fortran fireball-parsers get-lattice-params gulp-parsers]]
	      [utils :only [inter-cat-tree slash]]
	      [xyz-parser :only [parse-xyz]]
              [xyzSTATS.xyzcore :only [atoms]]
	      [ga :only [ga]]
	      [clojure.contrib.probabilities.random-numbers :only [rand-stream]]
	      [clojure.contrib.probabilities.monte-carlo :only [random-stream normal]]
	      [clojure.contrib.seq :only [frequencies rand-elt]]))

(defn rand-point []
  "Returns a random point [0..1] in fractional coordinates."
  (take 3 (repeatedly rand)))

(defn fpart [num]
  "Returns the fractional part of num."
  (- num (int num)))

(defn reimage-point [point]
  "Returns all components to the range [0..1]"
  (map #(fpart (first (drop-while neg? (iterate inc %)))) point))

(defn rand-vec []
  "Returns a random vector [-1..1]"
  (take 3 (repeatedly #(dec (rand 2)))))

(defn rand-lat-vectors []
  "Returns a 3x3 matrix of random lattice vectors."
  (take 3 (repeatedly rand-vec)))

(defn rand-lat-params []
  "Returns a seq of random values of [a b c alpha beta gamma]."
  (let [pi Math/PI, pi2 (/ pi 2)
	alpha (rand pi2), beta (rand pi2), abmax (+ alpha beta), abmin (Math/abs (- alpha beta))
	gamma (+ abmin (rand (- abmax abmin)))];respect the triangle inequality!
    [(rand) (rand) (rand) alpha beta gamma]))

(defn lat-volume [matrix]
  "Returns the volume of the unit cell defined by the seq of lattice vectors in matrix."
  (let [[a b c] matrix]
    (Math/abs (dot-product a (cross-product b c)))))

(defn scale-lat-vectors [volume lat-matrix]
  (let [v (lat-volume lat-matrix)
	m (Math/cbrt (/ volume v)) f #(* m %)]
    (map #(map f %) lat-matrix)))

(defn scale-lat-params [volume params]
  (let [v (lat-volume (apply parameters->vectors params))
	m (Math/cbrt (/ volume v)) f #(* m %)]
    (concat (map f (take 3 params)) (drop 3 params))))

(defn frac-cart [points matrix]
  "Converts points from fractional to cartesian coordinates."
  (map #(mat-vect-mult matrix %) points))

(defn afrac-cart [atms matrix]
  (let [points (frac-cart (map :coordinates atms) matrix)]
    (map #(assoc %1 :coordinates %2) atms points)))

(defn cart-frac [points matrix]
  "Converts points from cartesian to fractional coordinates."
  (let [mat (invert-matrix3 matrix)]
    (map #(reimage-point (mat-vect-mult mat %)) points)))

(defn acart-frac [atms matrix]
  (let [points (cart-frac (map :coordinates atms) matrix)]
    (map #(assoc %1 :coordinates %2) atms points)))

(defn average-lattice-vectors [lat1 lat2];perform weighted average with random weight
  (let [w (rand), y (- 1 w), f #(+ (* w %1) (* y %2))
	vol (lat-volume lat1)]
    (scale-lat-vectors vol (map #(map f %1 %2) lat1 lat2))))

(defn average-lattice-params [lat1 lat2];perform weighted average with random weight
  (let [w (rand), y (- 1 w), f #(+ (* w %1) (* y %2))
	vol (lat-volume (apply parameters->vectors lat1))]
    (scale-lat-params vol (map f lat1 lat2))))

(defn sample-distribution [rt]
  (clojure.contrib.generic.collection/seq (random-stream rt rand-stream)))

(defn strain-matrix []
  "Returns a random strain matrix per Equation 3 of J. Chem. Phys. 124, 244704 (2006)."
  (let [[e1 e2 e3 e4 e5 e6] (take 6 (filter #(min 1 (max -1  %)) (sample-distribution (normal 0 1))))
	one (int 1), two (int 2), e4half (/ e4 two), e5half (/ e5 two), e6half (/ e6 two)]
    [[(+ one e1) e6half e5half]
     [e6half (+ one e2) e4half]
     [e5half e4half (+ one e3)]]))

(defn mutate-lattice [structure]
  (let [strain (strain-matrix);strain matrix is column vectors.  need to transpose %?
	mutate #(map (partial mat-vect-mult strain) %)]
    (if-let [params (:lat-params structure)]
      (let [vectors (apply parameters->vectors params)
	    new-params (-> vectors mutate align-lattice vectors->parameters)]
	(println "old params:" (vec params))
	(println "new params:" (vec new-params));too much change may or may not indicate a problem
	(assoc structure :lat-params (scale-lat-params (lat-volume vectors) new-params)))
      (let [vectors (:lat-vectors structure)]
	(assoc structure :lat-vectors (->> vectors mutate (scale-lat-vectors (lat-volume vectors))))))))

(defn crossover-lattices [system1 system2 atoms]
  (if (:lat-params system1)
;    (assoc system1 :atoms atoms
;	   :lat-params (vectors->parameters (align-lattice (average-lattice-vectors (apply parameters->vectors (:lat-params system1))
;										    (apply parameters->vectors (:lat-params system2))))))
    (assoc system1 :lat-params  (average-lattice-params  (:lat-params system1)  (:lat-params system2))  :atoms atoms)
    (assoc system1 :lat-vectors (average-lattice-vectors (:lat-vectors system1) (:lat-vectors system2)) :atoms atoms)))

(defn plane-crossover [system1 system2]
  "Splits the unit cell with a randomly generated plane.
   Returns all atoms above the plane from atoms1 and all below for atoms2."
  (let [atoms (map (fn [atoms1 atoms2]
		     (let [plane [(rand-point) (normalize-mag (rand-point))]
			   dist #(point-plane-distance % plane)
			   pred (comp pos? dist :coordinates)]
		       (concat (filter pred atoms1)
			       (remove pred atoms2))))
		   (:atoms system1) (:atoms system2))]
    (crossover-lattices system1 system2 atoms)))

(defn subgrid-crossover [system1 system2]
  "Divides fractional space into a grid.  For each element of the grid,
   randomly selects all atoms from either system1 or system2."
  (let [atoms (map (fn [atoms1 atoms2]
		     (let [divisor (inc (rand 4)), step (/ 1.0 divisor)
			   ranges (partition 2 1 (concat (range 0 1 step) [1]))
			   between (fn [val [min max]] (and (>= val min) (< val max)))]
		       (apply concat (for [x ranges, y ranges, z ranges]
				       (let [atoms (if (> 0.5 (rand)) atoms1 atoms2)
					     bounded #(map between % [x y z])]
					 (filter #(every? bounded %) atoms))))))
		   (:atoms system1) (:atoms system2))]
    (crossover-lattices system1 system2 atoms)))

(defn check-distances [atoms lattice dist]
  "Returns a seq of all bond lengths less than dist."
  (let [temp (afrac-cart (supercell-atoms 1 atoms lattice) lattice)]
    ;todo: use bsptree
    (sort (set (for [a temp, b temp
		     :let [d (magnitude (map - (get-com a) (get-com b)))]
		     :when (and (not= a b) (< d dist))] d)))))

(defn add-atom [atms lattice dist species max-iter]
  "Adds an atom at a random coordinate. If intersecting, will retry upto max-iter times,
   and the least intersecting structure will be returned."
  (let [aseq (take max-iter (repeatedly #(conj atms (struct atoms species (rand-point)))))
	indexed (zipmap (map #(let [d (check-distances % lattice dist)]
				(if (empty? d) 0 (apply min d))) aseq) aseq)
	good (first (filter #(zero? (first %)) indexed))]
    (if good (second good);if a non-intersecting structure is found, return it
	(get indexed (apply min (keys indexed))))))

(defn remove-atom [atoms species]
  (let [temp (filter #(= species (:species %)) atoms)]
    (remove #(= (rand-elt temp) %) atoms)))

(defn fix-atoms [atoms lattice formula]
  "Fixes chemical composition of atoms. formula is a hmap of [species count]"
  (let [species (set (keys formula))
	freqs (frequencies (map :species atoms))
	diffs (merge-with - formula freqs)
	sep (into {} (map (fn [s] (vector s (filter #(= s (:species %)) atoms))) species))
	];separate atoms by species
    (loop [s (keys diffs), nums (vals diffs), atms atoms]
      (if (nil? s) atms
	  (let [num (first nums)
		new-atms (cond (zero? num) atms
			       (pos? num) (nth (iterate #(add-atom % lattice 0.75 (first s) 10) atms) num)
;Atoms are added in order of species.  This may or may not be disadvantageous for latter species.  Try round robin?
			       (neg? num) (nth (iterate #(remove-atom % (first s)) atms) (Math/abs num)))]
	    (recur (next s) (next nums) new-atms))))))

(defn remove-elt [elt coll]
  (remove #(= elt %) coll))

(defn permute-atoms [atoms formula]
  "Randomaly swaps the identities of two atoms. formula is a hmap of [species count]"
  (let [species (set (keys formula))
	freqs (frequencies (map :species atoms))
	diffs (merge-with - formula freqs)
	sep (into {} (map (fn [s] (vector s (filter #(= s (:species %)) atoms))) species))
	];separate atoms by species
    (let [species1 (rand-elt (seq species))
	  species2 (rand-elt (seq (disj species species1)))
	  swap (fn [s1 s2 hmap]
		 (let [elt (rand-elt (get hmap s2))]
		   (assoc hmap s2 (replace {elt (assoc elt :species s1)} (get hmap s2)))))]
      (->> sep
	   (swap species1 species2)
	   (swap species2 species1)
	   vals (apply concat)))))

(defn fingerprint [atoms lattice dr rmax]
  "Returns the pair correlation function - 1, discretized into bins dr wide."
  (let [temp (afrac-cart (supercell-atoms 1 atoms lattice) lattice)
	num (count temp), rho (/ num (lat-volume lattice))]
    (->> (for [a temp, b (afrac-cart atoms lattice)
	       :let [d (magnitude (map - (get-com a) (get-com b)))]
	       :when (and (not= a b) (< d rmax))]
	   {(int (/ d dr)) 1})
	 (apply merge-with +)
	 ((fn [x] (dissoc 0)));prevent divide by zero
	 (map (fn [[k v]] [k (dec (/ v num rho (* 4 Math/PI k k dr)))])))))

(defn rms-fingerprint [f1 f2]
  "Returns the rms distance between two fingerprints."
  (magnitude (map - (map second f1) (map second f2))))

(defn cos-fingerprint [f1 f2]
  "Returns the cosine distance between two fingerprints."
  (let [one (map second f1), two (map second f2)]
    (- 1 (/ (dot-product one two) (* (magnitude one) (magnitude two))))))

(defn quasi-entropy [prints dist-fn]
  "Returns the sum of the pairwise fingerprint distances."
  (reduce + (for [a prints, b prints] (dist-fn a b))))

(defn fb-fitness-fn [hmap]
  "Generates a structure from the given hmap, runs fireball, then parses and returns the logfile and answer.xyz"
  (let [directory (str "/srv/cache/jfennick" slash (System/nanoTime))
	files (-> required-files (conj [cgopt "cgopt.optional"]) (conj [kpts "gamma.kpts"]) (conj [pressure "pressure.optional"]))
	fpath "/home/users/jfennick/workspace/progs-cp/src/fireball.x"
	params (:lat-params hmap)
	lattice (if params (apply parameters->vectors params) (:lat-vectors hmap))
	cart (afrac-cart (last (:atoms hmap)) lattice)
	temp (try (run-fireball fpath directory files (write-bas cart) (inter-cat-tree ["\n" "\t"] lattice) "f.log" hmap)
		  (catch RuntimeException e))
	log (parse-fortran (str directory slash "f.log") fireball-parsers)
	ks [:runtime :etot]
	bad-keys (remove log ks)]
    (if-not (empty? bad-keys);Ensure fireball ran and parsed correctly, so we don't get weird errors later.
      (println "Unparseable keys:" (vec bad-keys))
      (let [xyz (parse-xyz (str directory slash "answer.xyz"))
	    new-lat (:lat-vectors log);Fireball seems to use row vectors, not column vectors.  need to transpose?
	    frac (map acart-frac xyz new-lat)
	    alat (align-lattice (last new-lat))
	    #_(fingerprint (last frac) a-lat 0.01 10)]
	(if params
	  (assoc hmap :log log :atoms frac :lat-params (vectors->parameters alat))
	  (assoc hmap :log log :atoms frac :lat-vectors alat))))))

(defn gulp-fitness-fn [hmap]
  (let [dir (str (:directory hmap) slash (System/nanoTime))
	temp (try (run-gulp dir hmap) (catch RuntimeException e))
	log (parse-fortran (str dir slash "output.gout") gulp-parsers)
	ks [:total-time :final-lat-params :final-enthalpy]
	bad-keys (remove log ks)]
    (if-not (empty? bad-keys);Ensure gulp ran and parsed correctly, so we don't get weird errors later.
      (println "Unparseable keys:" (vec bad-keys))
      (let [xyz (parse-xyz (str dir slash "outputmovie.xyz"))
	    temp1 (first (:final-lat-params log)), f #(* Math/PI (/ % 180))
	    temp2 (concat (take 3 temp1) (map f (drop 3 temp1)))
	    lat-params (repeat (count (first xyz)) temp2)
	    ;when doing optimization, only final lattice parameters are written out.
	    lat-vectors (map #(apply parameters->vectors %) lat-params)
	    frac (map acart-frac xyz lat-vectors)]
	(get-symmetry (assoc hmap :log log :atoms frac :lat-params (last lat-params) ;:lat-vectors (last lat-vectors)
			     #_(fingerprint (last frac) (last lat-vectors) 0.01 10)))))))

(defn init-population [options]
  (->> #(let [v (:volume options), init (:lat-init options), type (:lat-type options)
	      lat (if (= init 'params)
		    (let [x (scale-lat-params v (rand-lat-params))]
		      (if (= type 'params) x (apply parameters->vectors x)))
		    (let [x (scale-lat-vectors v (rand-lat-vectors))]
		      (if (= type 'params) (vectors->parameters x) x)))
	      temp (if (= 3 (count lat)) lat (apply parameters->vectors lat))
	      atoms (fix-atoms [] temp (:formula options))]
	  (-> (if (= type 'params)
		(assoc options :atoms [atoms] :lat-params lat)
		(assoc options :atoms [atoms] :lat-vectors lat))
		get-symmetry))
       repeatedly (map (:fit-fn options)) (remove nil?)))

(def formula {"C" 7 "Si" 7})
(def gulp-ops {:keywords '#{conp optimise numdiag spatial}
	       :directory "/srv/cache/jfennick"
	       :gulpx "/home/users/jfennick/gulp3.4/gulp"
	       :fit-fn gulp-fitness-fn :fit-keyfn #(last (get-in % [:log :final-enthalpy]))
	       :icluster 0 :volume 300 :formula formula :lat-type 'params :lat-init 'params
	       :potentials {'brenner []}
	       :opt-vals '{ensemble [npt 0.05 0.05]
			   temperature 300
			   ftol [opt 6], gtol [opt 6]
			   timestep 0.0001
			   production 1.0
			   pressure 0.000101325
			   output [movie xyz outputmovie.xyz] }})
(def fb-ops {:fdata "/home/users/jfennick/Fdata3_Hs3.8_Bs4.4p4.9"
	     :iquench 1, :max-scf 40, :icluster 0, :itheory 1, :kpts "gamma.kpts", :timesteps "1,1000"
	     :deltat 0.1, :iimage 1, :cell-mass 1000.0, :idynmat 1;ipressure has the same index as idynmat
	     :formula formula, :volume 155.1, :big4 "1 1 1 1", :dynamics 2
	     :fit-fn fb-fitness-fn :fit-keyfn #(last (get-in % [:log :etot]))
	     ;:kpoints (normalize-weights (let [mesh (range -0.5 0.51 1)]
		});			   (for [a mesh, b mesh, c mesh] [a b c 1]))) })
(def init-options {:deltat 0.1 :timesteps "1,100" :cell-mass 1000.0});use ultra-conservative initial settings

(def ga-ops {:mut-prob 0.5, :mut-fn #(assoc % :atoms [(permute-atoms (last (:atoms %)) formula)])
	     :lat-prob 0.5, :lat-fn mutate-lattice :cross-prob 0.95 :pop-size 4
	     :cross-fn (comp #(let [lat (if-let [params (:lat-params %)] (apply parameters->vectors params) (:lat-vectors %))
				    atoms (fix-atoms (last (:atoms %)) lat formula)]
				(assoc % :atoms [atoms])) plane-crossover)
	     :init-fn #(init-population (merge gulp-ops init-options)) })
(def options (merge ga-ops gulp-ops))

(defn evolve [options init-options];deprecated. use ga instead.
  "Performs genetic algorithm using the given options."
  (lazy-seq
    (let [pop (:population options (take (:pop-size options) (init-population (merge options init-options))))
	  temp1 (take (:pop-size options)
		      (repeatedly #(if (> (rand) (:cross-prob options))
				     (rand-elt pop)
				     ((:cross-fn options) (rand-elt pop) (rand-elt pop)))))
	  temp2 (take (:pop-size options)
		      (repeatedly #(if (> (rand) (:mut-prob options))
				     (rand-elt temp1)
				     ((:mut-fn options) (rand-elt temp1)))))
	  temp3 (take (:pop-size options)
		      (repeatedly #(if (> (rand) (:lat-p options))
				     (rand-elt temp2)
				     ((:lat options) (rand-elt pop)))))]
	(cons pop (evolve (assoc options :population (map (:fit-fn options) temp3)) init-options)))))