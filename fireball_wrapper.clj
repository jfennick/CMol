(ns fireball-wrapper (:use [atomic-structure-output :only [write-bas write-xyz]]
                           [geometry :only [normalize-sum]]
			   [utils :only [endline inter-cat-tree read-table slash temp-dir]]
			   [tree-utils :only [tree-apply]]
			   [fortran-parser :only [parse-fortran fireball-parsers]]
			   [empirical-data :only [atomic-numbers]]
			   [xyz-parser :only [parse-xyz]]
			   [not-clima :only [update-mol-coords]]
			   [viewer-2d :only [find-window plot-points viewer-2d]]
			   [clojure.contrib.shell :only [sh]]
                           [clojure.string :only [join]])
    (:import java.io.File))

(defn pad-line [val desc]
  "Returns val padded to 35 spaces, then description."
  (let [s (str val)]
    (apply str (concat s (repeat (- 35 (count s)) " ") desc))))

(defmacro def-input
  "Creates defns for fireball input files."
  [name keys vals descs]
  `(defn ~name
     ([] (~name {}))
     (~'[hmap]
	(let [~'keys ~keys ~'vals ~vals ~'descs ~descs
	      ~'newmap (merge (zipmap ~'keys ~'vals) ~'hmap)
	      ~'lines  (map (fn ~'[key desc]
			      (pad-line (get ~'newmap ~'key) ~'desc)) ~'keys ~'descs)]
	  (join ~endline ~'lines)))))

(def-input barrier
  [:strength :rms-diff :steps :push :slow :fraction :bfile]
  [0.005 0.01 0 0 32000 1 "final.bas"]
  ["strength for pushing atoms"
   "rms diff between current/final configurations"
   "quench velocities every X steps"
   "push atoms heading in right direction"
   "slow down atoms that are going too fast in right direction"
   "save this fraction of the bad forces"
   "Final basis file"])

(def-input cgopt
     [:max-disp :scale-disp :etot-tol :ftot-tol :max-steps :max-iter :refine]
     [0.10 0.5 0.00001 0.01 100 3 0]
     ["max. displacement"
      "scale displacement by this factor if next energy > current energy"
      "etot tolerance"
      "ftot tolerance"
      "max. number of CG steps"
      "max. number of internal iterations"
      "refine conjugate gradient results with quenching (0/1 N/Y)"])

(def-input diagnostics
  [:2overlap :2vna_ontopl :2vna_ontopr :2vna_atom-atom :2non-local
   :2xc_ontop :2xc_atom-atom :2xc_correction :2z-dipole :2y-dipole
   :2x-dipole :2coulomb :2kinetic :2Hubbard :2den_ontopl
   :2den_ontopr :2den_atom :2denS_ontopl :2denS_ontopr :2denS_atom
   :2overlapS :2coulomb-Hubbard :2srEwald :2lrEwald :3eutral-atom
   :3xc :3OLSXC :3OLSXC-spheric :test]
  [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "1, 0.0"]
  ["2 center overlap"   "2 center vna_ontopl"   "2 center vna_ontopr"
   "2 center vna_atom-atom"   "2 center non-local"   "2 center xc_ontop"
   "2 center xc_atom-atom"   "2 center xc_correction"   "2 center z-dipole"
   "2 center y-dipole"   "2 center x-dipole"   "2 center coulomb"
   "2 center kinetic"   "2 center extended Hubbard"   "2 center den_ontopl"
   "2 center den_ontopr"   "2 center den_atom"   "2 center denS_ontopl"
   "2 center denS_ontopr"   "2 center denS_atom"   "2 center overlapS"
   "2 center coulomb (extended Hubbard)"   "2 center short range Ewald"
   "2 center long range Ewald"   "3 center neutral-atom"
   "3 center exchange-correlation"   "3 center average density OLSXC"
   "3 center average density OLSXC (spheric)"   "itestrange, rangetest"])

(def-input dos
  [:scale :natom :nsteps :step :write :min-max :green]
  [1 "1        5" 150 "-10      0.1" 0 "-0.51   -0.51" 0.05]
  ["scale factor of coord / ratio of lattice parameter"
   "natom_beg, natom_end for dos calculation"
   "number of energy steps"
   "first energy and step for dos calculation"
   "1/0 yes/no write the tip_e_str.inp"
   "minimum and maximum energies for the writing the tip"
   "(eta) imaginary part for green function in calculation"])

(def-input dyn
  [:disp :vec :output :flag]
  [0.03 "1 1 1" "dyn.dat" 0]
  ["elementary displacement"
   "dimension vector"
   "output file"
   "flag on list of free atoms"])

(defn nh [freqs]
  "Returns the output of nh.optional Use for Nose-Hoover thermostat"
  (let [num (count freqs)
	ints (range 1 (inc num))
	lines (map #(pad-line %1 (str "characteristic frequency of link " %2)) freqs ints)]
    (join endline (cons (pad-line num "number of links on NH chain") lines))))

(def-input options
  [:itheory :ispin :itheory-xc :max-scf :iqout :qstate :iquench :iensemble :ibarrier
   :big4 :ifixc :ifixn :icluster :iordern :iumbrella :ivdw :igauss :iimage :ipath :idynmat
   :iexfield :ithermoint :ireduce :iendtemp]
  [0 0 0 20 1 0.0 0 0 0 "0 1 1 1" 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
  ["itheory (0=Harris, 1=DOGS, 2=Ext-hubbard)"   
   "ispin"
   "itheory_xc (0=HXC, 1=SNXC, 2=OLSXC)"
   "max_scf_iterationis (irrelevant for Harris)"
   "iqout (Lowdin = 1, Mulliken = 2)"
   "qstate"   
   "iquench"   
   "iensemble"   
   "ibarrier"
   "The big 4 constraints"   
   "ifixcharges"   
   "ifixneigh"
   "icluster (0=Periodic, 1=Cluster)"   
   "iordern"
   "iumbrella"   
   "ivdw (Van der Waals)"
   "igauss (Gaussian Approach)"   
   "iimage (Reimage atoms)"
   "ipathintegral"   
   "idynmat (Dynamic Matrix)"   
   "iexfield"
   "ithermoint"   
   "ireducekpts"   
   "iendTemp"])

(def-input output
  [:cdcoefs :charges :density :eigen :fermi :fpieces :hampieces
   :components :neigh :neigh-com :xyz :dos :hop :atom :pop :hs :vel]
  [0 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0]
  ["iwrtcdcoefs (coefficients for wavefunctions)"
   "iwrtcharges (writes out Lowdin or Mulliken)"
   "iwrtdensity (write out density matrix)"
   "iwrteigen (writes out the energy eigenvalues)"
   "iwrtefermi (write out fermi occupations)"
   "iwrtfpieces (writes out force components)"
   "iwrthampiece (write pieces of H matrix)"
   "iwrtcomponents (write components of ebs energy)"
   "iwrtneigh (write out neighbor map)"
   "iwrtneigh_com (write out common neighbor map)"
   "iwrtxyz (write out xyz movie file)"
   "iwrtdos (write out dos files)"
   "iwrthop (write out hopping values for STM)"
   "iwrtatom (write out the Atomo_i files)"
   "iwrtpop (write out population file)"
   "iwrtHS (write out H & S file)"
   "iwrtvel (write out VELOCITY.dat file)"])

(def-input pressure
  [:gpa :cell-mass :dynamics :velocity :lquench :betaq :alphaq :volfile :s-kinfile :l-kinfile :latvecfile]
  [0.000101325 1000.0 3 1 1 0.6 0.4 "nameofvolume_file" "nameofskinetic_file" "nameoflkinetic_file" "nameoflatvec_file"]
  ["Pressure in GPa"   "Mass of cell"
   "Type of cell dynamics (1, 2 or 3) (Parrinello-Rahman, Projected force, or Wentzcovitch)"
   "Scaled velocity quench"   "Lattice quench"   "betaQ"   "alphaQ"
   "nameofvolume_file"   "nameofskinetic_file"
   "nameoflkinetic_file"   "nameoflatvec_file"])

(def-input quench
  [:etol :ftol :twant :taurelax]
  [0.001 0.05 0.0 50.0]
  ["energy tolerance"
   "force tolerance"
   "Temperature wanted"
   "taurelax for simulated annealing"])

(defn sa [temp disp flag atoms seed]
  "Returns SA.optional. atoms is a seq of 2-tuples of (atom# 1 or 0)"
  (let [lines (map #(pad-line (str (first %) ",    " (second %))
			      (str (if (= 0 (second %)) "do not ")
				   "anneal atom " (first %))) atoms)
	headers (list (pad-line temp "simulated annealing temperature")
		      (pad-line disp "maximum displacement per step")
		      (pad-line flag "1 to anneal all atoms, 0 to specify which atoms to anneal"))
	]
    (join endline (concat headers lines (list (pad-line seed "integer seed value for the random number generator"))))))

(def-input scf
  [:sigma :bmixfactor :fermitemp]
  [0.00001 0.08 0.0]
  ["sigmasmall (Default=0.0001)"
   "BMIX mixing factor (Default=0.08, or 8%)"
   "Fermi temperature (Smearing)"])

(def-input script
  [:basfile :latfile :kpts :timesteps :deltat :xvfile :acfile :tinit :tend]
  ["input.bas" "far.lvs" "automatic" "1,1" 0.25 "fire.xv" "fire.ac" 0.0 0.0]
  ["basis file"   "lattice vector filename"
   "kpt preference"   "initstep,finalstep"
   "time step (fs)"   "xv file name"
   "ac file name"   "initial temperature"
   "end temperature"])

(defn umbrella [filename time atoms]
  "atoms is a seq of 4-tuples of (atom1# atom2# bond-distance K)"
  (let [num (count atoms)
	ints (range 1 (inc num))
	lines (map #(list (pad-line (str (first %1) " " (second %1))
				    (str "atomic pair " %2))
			  (pad-line (nth %1 2) "bond distance")
			  (pad-line (nth %1 3) "K")) atoms ints)]
    (join endline (concat (list (pad-line filename "umbrella output file")
				     (pad-line time "equilibration time"))
			       (apply concat lines)))))

(defn kpts [hmap]
  "Points is a seq of 4-tuples: x,y,z coordinates in k-space, then the weight.  The weights must sum to 1."
  (let [points (get hmap :kpoints [[0 0 0 1]])]
    (str (count points) endline (inter-cat-tree [endline "\t"] points))))

(defn normalize-weights [points]
  (let [weights (normalize-sum (map #(nth % 3) points))]
    (map #(concat (take 3 %1) (list (double %2)))
	 points weights)))

(defn normalize-kpts [kptfile]
  "Read in a kptfile and overwrite it with normalized weights."
  (spit kptfile (kpts {:kpoints (normalize-weights (next (read-table kptfile)))})))
	
;todo vdw, fragments

(defn wrap [x]
  (if (sequential? x) x [x]))

(defn fragments [hmap]
  ;switch is 0 or 1 - fix atoms or fix center-of-mass
  (let [switch (get hmap :fswitch 0)
	indices (get hmap :fragments)
	ifn #(vec (take 4 (concat (wrap %) (repeat 1))))
	contents (vec (map #(vec (cons [(count %)] (map ifn %))) indices))]
    (str switch endline (count indices) endline (inter-cat-tree ["\t" "\n" "\t"] contents))))

(defn fdata [hmap]
  (get hmap :fdata))

(def required-files
  [[options "options.input"]
   [output "output.input"]
   [quench "quench.optional"]
   [script "script.input"]
   [fdata "Fdata.optional"]])

(defn run-fireball [fireballx directory files bas lvs logfile hmap]
  "Writes files with settings from hmap into the given directory and runs fireball.x
   files is a seq of 2-tuples of (fn filename). bas and lvs are the contents of the bas and lvs files."
  (. (File. directory) mkdirs)
  ;todo cleanup directory
  (doseq [file files]
    (spit (str directory slash (second file))
      ((first file) hmap)))
  ;write bas, lvs, and script files
  (spit (str directory slash (get hmap :basfile "input.bas")) bas)
  (spit (str directory slash (get hmap :lvsfile "far.lvs")) lvs)
  (spit (str directory slash "runFireball")
        (str "ulimit -s unlimited" endline fireballx))
  ;make script executable
  (doto (File. (str directory slash "runFireball"))
    (.setExecutable true))
  ;execute fireball
  (spit (str directory slash logfile)
    ;for large logfiles, may need to read and write streams directly
    (sh "./runFireball" :dir directory)))

(defn scale-nums [scale nums]
  "Multiplies each leaf of nums by scale. nums may be any tree"
  (tree-apply #(* scale %) false false nums))

(defn scale-atoms [scale atoms]
  "scales the :coordinates of a seq of atoms struct-maps"
  (map #(merge % {:coordinates (scale-nums scale (get % :coordinates))}) atoms))


(defn eq-lat-opt [smin smax step fireballx files atoms lattice logfile hmap]
  "Performs equilibrium lattice optimization on atoms.
   Scales the unit cell from smin to smax in step, returning
   a vector of [etot scale atoms lattice] of the structure with minimum etot."
  (let [scale-factors (range smin (+ smax step) step)
	result (map #(try
		      (let [sdir (temp-dir)
			    fdir (File. sdir)
			    bas  (scale-atoms % atoms)
			    sbas (write-bas bas)
			    lvs  (scale-nums % lattice)
			    slvs (inter-cat-tree [endline "\t"] lvs)
			    temp (run-fireball fireballx sdir files sbas slvs logfile hmap)
			    parse (parse-fortran (str sdir slash logfile) fireball-parsers)
			    xyz (parse-xyz (str sdir slash "answer.xyz"))]
			(println (str "scale: " % "etot/atom: " (:etot/atom parse)))
		        ;(delete-dir fdir)
		        ;for optimized structure, min and last element should be equal
			[(apply min (:etot/atom parse)) % (last xyz) lvs])
			;todo: when optimizing, return parsed coordinates, not simply scaled coordinates
		      ;if fireball crashes, return worst possible value
		      (catch RuntimeException e [Double/POSITIVE_INFINITY % 0 0]))
		    scale-factors)
	sorted (sort #(< (first %1) (first %2)) result)
	coords (ref (partition 2 (interleave scale-factors (map first result))))]
    (viewer-2d (plot-points coords) (find-window coords))
    (first sorted)))