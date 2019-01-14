(ns gulp-wrapper (:use [clojure.contrib.except :only [throw-if]]
		       [clojure.contrib.shell :only [sh]]
		       [clojure.string :only [join]]
                       [empirical-data :only [atomic-symbols]]
		       [fortran-parser :only [parse-fortran gulp-parsers]]
                       [geometry :only [tolerance?]]
                       [not-clima :only [min-max-coordinates update-mol-coords shift]]
		       [utils :only [endline inter-cat-tree slash squared-diff]]))

(def all-keywords
     '#{neb nodneb oldunits nomodcoord spatial nomolecularinternalke nononanal zero_potential
	montecarlo nomcediff marvinSE noelectrostatics zsisa dcharge qeq minimum_image
	noquicksearch prt_eam prt_two newda qok dipole oldelastic conserved hexagonal
	c6 nofirst_point fix_molecule norepulsive_cutoff libff libdump save restore nod2sym
	md nolist_md sm bulk_noopt noanisotropic_2b defect regi_before optimise gradient
	transition_state rfo distance eem qbond gasteiger pot nodpsym conv conp noflags shell
	breathe nobreathe cellonly isotropic orthorhombic free_energy static_first nozeropt
	noenergy nodensity_out unfix debug efg single bond angle linmin nosymmetry nodsymmetry
	fit fbfgs genetic umorse preserve_Q simultaneous hill voigt property molecule outcon
	compare unit lbfgs numdiag dfp positive conjugate hessian nosderv molq molmec
	noautobond relax phonon dynamical_matrix eigenvectors lower_symmetry meanke nokpoints
	kfull nofrequency full noksymmetry average norecip noreal noexclude storevectors
	operators predict anneal global Cost })

(defn write-keywords [options]
  (let [key-set (:keywords options)
	unknown (remove all-keywords key-set)]
    (if-not (empty? unknown)
      (println (apply str "Unknown keywords: " (interpose " " unknown))))
    (apply str (interpose " " (filter all-keywords key-set)))))

(defn write-vectors [options]
  (cond
    (= (count (:lat-vectors options)) 3)
  (apply str "vectors\n" (inter-cat-tree ["\n" "\t"] (:lat-vectors options)))

    (= (count (:lat-vectors options)) 2)
  (apply str "svectors\n" (inter-cat-tree ["\n" "\t"] (:lat-vectors options)))

    (= (count (:lat-vectors options)) 1)
  (apply str "pvector\n" (inter-cat-tree ["\n" "\t"] (:lat-vectors options)))))

(defn write-params
  "gulp uses degrees for all angles, but we assume radians (to make math easier)."
  [options]
  (cond
    (= (count (:lat-params options)) 6)
    (let [params (:lat-params options)
          angles (map #(/ (* 180 %) Math/PI) (drop 3 params))]
      (apply str "cell " endline (interpose " " (map float (concat (take 3 params) angles)))))

    (= (count (:lat-params options)) 3)
    (let [params (:lat-params options)
          angles (map #(/ (* 180 %) Math/PI) (drop 2 params))]
      (apply str "scell " endline (interpose " " (map float (concat (take 2 params) angles)))))

    (= (count (:lat-params options)) 1)
    (apply str "pcell" endline (:lat-params options))))

(defn write-potentials [options]
  (let [default-fn (fn [& args] (apply str (interpose " " args)))
         pots {'brenner (constantly "brenner \n")
                  'brenner1 (constantly "brenner 1\n")}]
    ;todo: check for unknown potentials
    (apply str (map #(apply (get pots % default-fn) (get (:potentials options) %))
              (keys (:potentials options))))))

(defn gulp-atom-string
  "This will change the atoms struct map to a string for one particular atom in
a mol.  This code still needs work (ion-radius and occupancy).  This is designed
to be used by gulp-mol-string.  In order to tell GULP that a particular atom should
be forced to not move in a (x, y, z, or all) direction(s) you should place the
following in the :name value of the atom (:xstop, :ystop, :zstop, :allstop). If
an atom is a shell atom place the key :shell in the atom's :name set.  A clear
way of adding these keys to :name is using not-clima/update-mol-name."
   [atom]
  (let [species_temp (:species atom)
        species (if (string? species_temp) species_temp (atomic-symbols species_temp))
        coords (join " " (map float (:coordinates atom)))
        charge (if (nil? (:charge atom)) 0 (:chage atom))
        occupancy 1
        ion-radius 0.0
        name (:name atom)
        coreshell (if (:shell name) "shell" "core")
        xstop (if (or (:allstop name)(:xstop name)) 0 1)
        ystop (if (or (:allstop name)(:ystop name)) 0 1)
        zstop (if (or (:allstop name)(:zstop name)) 0 1)]
    (join " " (vector species coreshell coords charge occupancy ion-radius
                                xstop ystop zstop))))

(defn gulp-mol-string
  "This will change the atoms struct of a mol into a string.  This is designed
to be used by write-coordinates.
Usage: (gulp-mol-string water) => (str O 0.0 0.0 0.0 0 1 0.0 0 0 0 \nH 0.0 0.86 0.0 0 1 0.0 0 0 0 \nH 0.0 0.0 0.86 0 1 0.0 0 0 0)"
   [mol]
  (join endline (map #(gulp-atom-string %) mol)))

(defn write-coordinates
  "While this states that this writes the coordinates to the GULP input file, it
actually writes the coordinates and the lattice parameters/vectors.  When the
system of interest is a unit cell (or supercell), besides specifying the
parameters/vectors using :lat-params/:lat-vectors and the coordinates using :atoms,
the user should specify the keyval :coord-type with string vals of fractional or
cartesian .  If the system of interest is a molecule, :lat-params and
 :lat-vectors should not be specified, and :coord-type need not be.


Now Deprecated:  :icluster"
   [options]
  (throw-if (and (not (:lat-params options)) (not (:lat-vectors options)) (= "fractional" (:coord-type options)))
    ":lat-params and :lat-vectors are not set while :coord-type = fractional.")
    (throw-if (and (or (:lat-params options) (:lat-vectors options))
                (not (or (= "fractional" (:coord-type options)) (= "cartesian" (:coord-type options)))))
    ":lat-params or :lat-params are set while :coord-type is not set properly.")

  (let [coords (gulp-mol-string (:atoms options))]; assume fractional/cartesian conversion has already been done
    (cond (:lat-params options)
      (str (write-params options) endline (:coord-type options) endline coords)

      (:lat-vectors options)
      (str (write-vectors options) endline (:coord-type options) endline coords)

      (and (not (:lat-vectors options)) (not (:lat-params options)))
      (str "cartesian" endline coords))))

(defn write-options-vals
  "Converts each option/value pair into a string.  Returns a concatenated string separated by newlines."
  [options]
  (->> (:opt-vals options)
    (map (fn [[opt val]] (if (sequential? val)
                           (apply str opt " " (interpose " " val))
                           (str opt " " val))))
    (interpose endline)
    (apply str)))

(defn build-input [options]
  (let [builders [write-keywords write-coordinates write-potentials write-options-vals]]
    (apply str (map #(str (% options) endline) builders))))

(defn run-gulp
  ([dir options]
    (. (java.io.File. dir) mkdirs)
    (spit (str dir slash "runGULP")
      (str "ulimit -s unlimited" endline
        (:gulpx options) " < input.gin > output.gout"))
    ;make script executable
    (doto (java.io.File. (str dir slash "runGULP"))
      (.setExecutable true))
    (spit (str dir slash "input.gin")
      (build-input options))
    (sh "./runGULP" :dir dir))

  ([dir options folder-name]
    (. (java.io.File. dir) mkdirs)
    (spit (str dir slash folder-name)
      (str (:gulpx options) " < " folder-name ".gin > " folder-name ".gout"))
    ;make script executable
    (doto (java.io.File. (str dir slash folder-name))
      (.setExecutable true))
    (spit (str dir slash folder-name ".gin")
      (build-input options))
    (sh (str "./" folder-name) :dir dir)))
