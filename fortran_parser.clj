(ns fortran-parser (:use [merge-utils :only [wrap-merge]]
			 [utils :only [empty-str last-widths read-table transpose]]
			 [clojure.contrib.io :only [read-lines]]
                         [clojure.contrib.seq :only [find-first]])
                   (:require clojure.contrib.string))

(defn parse-fortran [logfile parsers]
  "parsers is a seq of (predicate parser key-symbol).
   This function returns a hashmap."
  (when (.exists (java.io.File. #^String logfile))
    (loop [lines (read-lines logfile)
	   hmap {}]
      (if (empty? lines)
	hmap;end of logfile
	(let [tuples (filter #((first %) lines) parsers)
	      keys (map #(nth % 2) tuples)
	      vals (map #((second %) lines) tuples)]
	  (if (< 1 (count tuples)) (println "count" (count tuples)))
	  (recur (next lines) (wrap-merge hmap (apply hash-map (interleave keys vals)))))))))

(defn string-starts [#^String str lines]
  "Returns true if (first lines) starts with the given string literal."
  (let [#^String one (first lines)]
    (.. one trim (startsWith str))))

(def pat #"\s*((\s*\S+)+)\s*=\s*(.*)")

(defn string-matches [lines]
  (println (first lines))
  (re-matches pat (first lines)))

(defn string-contains [#^String str lines]
  "Returns true if (first lines) contains the given string literal."
  (let [#^String one (first lines)]
    (.. one trim (contains str))))

(defn star-str [#^String s]
  (every? #(= \* %) (.trim s)))

(defn dash-str [#^String s]
  (every? #(= \- %) (.trim s)))

(defn str-skip [#^String s #^String str]
  "Skips past the first instance of str. Returns the remainder of s, assuming str found"
  (let [#^String st (. s substring (+ (count str) (. s indexOf str)))] st))

#_(defn read-next [#^String str]
  "splits string by whitespace, then reads and returns the next element."
  (read-string (last (. str split "\\s+"))))

(defn read-next
  "splits string by whitespace, then reads and returns the first element that is a number."
  [str]
    (read-string (first (clojure.contrib.string/grep #"\b[1234567890]+" (. str split "\\s+")))))

(defn read-rest [#^String s]
  "splits s by whitespace and returns all remaining elements."
    (map read-string (.. s trim (replace "D" "E") (split "\\s+"))))

(def fb-key-vals1 '[("Cohesive Energy per atom  =" :cohesive/atom) ("CohesiveE =" :cohesive) ("sigma =" :sigma) ("fxcot_ca  =" :fxcot_ca) ("Atomic Energy =" :atomic) ("Etot/atom =" :etot/atom) ("ETOT =" :etot) ("uxcdcc =" :uxdcc) ("etotxc_1c =" :etotxc_1c) ("uii - uee =" :uii-uee) ("ebs =" :ebs) ("Fermi Level =" :fermi) ("deltaE/atom  (meV) =" :deltae/atom) ("grand total energy per atom =" :grandtotal/atom) ("total energy per atom =" :total/atom) ("Grand Total = nuclear kinetic + potential =" :grandtotal) ("nuclear kinetic energy =" :nuclear) ("average temperature       =" :taverage) ("T_instantaneous " :tinstantaneous) ("Density of States" :dos) ("Deviation (rms) of input/output charges =" :rmsdeviation) ("Deviation (max) of input/output charges =" :maxdeviation) ("Fi_max=" :fi_max) ("Etot  RES =" :etotres) ("Fmax  RES =" :fmaxres) ("deltae, etotnew, etotold =" :deltae) ("gg0=" :gg0) ("gg1=" :gg1) ("orth =" :orth) ("beta =" :beta) ("sqrh =" :sqrh) ("frs  =" :frs) ("alpha=" :alpha) ("dx=" :dx) ("drmax=" :drmax) ("etot0  =" :etot0) ("etot1  =" :etot1) ("Fproj  =" :fproj) ("alpha =" :alpha2) ("A    =" :a) ("B    =" :b) ("C    =" :c) ("gamma_min=" :gamma_min) ("comp_KE     =" :comp_ke) ("comp_VNA    =" :comp_vna) ("comp_VXC    =" :comp_vxc) ("comp_VNL    =" :comp_vnl) ("comp_VXC_1C =" :comp_vxc1c) ("comp_VCA    =" :comp_vca) ("comp_VXC_CA =" :comp_vxc_ca) ("comp_EWDLR  =" :comp_ewdlr) ("comp_EWDSR  =" :comp_ewdsr) ("comp_EBS =" :comp_ebs) ("EBS =" :ebs2) ("error =" :error) ("vdw =" :vdw) ("dt =" :deltat)])

(def fb-key-vals3 '[("center of mass position =" :com-p) ("center of mass velocity =" :com-v) ("center of mass angular momentum =" :com-m)])

(def gulp-key-vals1 [["Brenner potentials         =" :brenner-potE]
                     ["Calculation of Brenner" :brenner-time]
                     ["Calculation of real space" :real-time]
                     ["Cell angle : alpha (o) =" :lat-param-al]
                     ["Cell angle : beta  (o) =" :lat-param-be]
                     ["Cell angle : gamma (o) =" :lat-param-ga]
                     ["Cell parameter : a (A) =" :lat-param-a]
                     ["Cell parameter : b (A) =" :lat-param-b]
                     ["Cell parameter : c (A) =" :lat-param-c]
                     ["Cell volume :   (A**3) =" :cell-volume]
                     ["Dimensionality =" :icluster]
                     ["Entropy (eV/K) =" :entropy]
                     ["Equilibration time        =" :equilibration-time]
                     ["Final enthalpy =" :final-enthalpy (comp first read-rest)]
                     ["Final energy =" :final-energy]
                     ["Force start time          =" :start-time]
                     ["Friction for temperature bath =" :friction]
                     ["Global summation overhead" :summation-time]
                     ["Heat capacity -const vol. (eV/K) =" :heat-capacity-const-vol]
                     ["Helmholtz free-energy        = " :Helmholtz-free-energy]
                     ["Interatomic potentials     =" :interatomic-potE]
                     ["Kinetic energy    (eV) =" :kineticE]
                     ["Monopole - monopole (real) =" :monopole-potE]
                     ["No. of degrees of freedom =" :degrees-of-freedom]
                     ["No. of mobile ions        =" :num-mobile-ions]
                     ["Number of CPUs =" :num-cpus]
                     ["Potential energy  (eV) =" :potentialE]
                     ["Production time           =" :production-time]
                     ["Sampling frequency        =" :sampling-freq]
                     ["Scaling frequency         =" :scaling-freq]
                     ["Scaling time              =" :scaling-time]
                     ["Temperature       (K)  =" :temperature]
                     ["Temperature of configuration =" :init-temp]
                     ["Time limit =" :time-limit]
                     ["Time step                 =" :deltat]
                     ["Total CPU time" :total-time]
                     ["Total energy      (eV) =" :totalE]
                     ["Total lattice energy       =" :etot-lattice]
                     ["Total number atoms/shells =" :num-atoms]
                     ["Total number of configurations input =" :num-configs]
                     ["Write frequency           =" :write-freq]])

(defn bulk-parser [lines skip end-pred item-parser]
  "Used to parse aggregate quantities.  both functions take (first lines) an as argument, and item-parser must return a seq"
  (loop [result []
	 temp (skip lines)]
    (if (or (nil? temp) (end-pred (first temp)))
      result
      (recur (concat result (item-parser (first temp)))
	     (next temp)))))

(defmacro str-fn [& body]
  "Returns a type-hinted fn of one String arg. Intentionally captures str."
  `(fn ~'[#^String str] ~@body))

(defn contains= [#^String s]
  (. s contains "="))

(def force-comp
     '[;band structure force components
       ("The kinetic force:" ft)
       ("The neutral atom force:" fna)
       ("The neutral atom force, 3-center:" f3na)
       ("The non-local force:" fnl)
       ("The exchange-correlation force:" fxc)
       ("The exchange-correlation force, 3c-center:" f3xc)
       ("The exchange-correlation force atm:" fxcatm)
       ("The exchange-correlation force ontop:" fxcot)
       ;total force components
       ("The band-structure force:" fbs)
       ("The short-range force:" dusr)
       ("The exchange-correlation double-counting force:" dxcv)
       ("The overlap-repulsive force:" fro)
       ("The grand total force (eV/A):" ftot)])

(defn force-comp-parser [symbol lines]
  (let [s1 (str symbol)
	s2 (repeat (- 10 (count s1)) " ")
	s3 (apply str (concat s1 (apply str s2) "="))]
    (bulk-parser lines next empty-str
		 (str-fn (list (read-rest (str-skip str s3)))))))

(defn lattice-parser [lines]
  (bulk-parser lines next star-str (str-fn (list (read-rest (str-skip str "Corrected LATTICE VECTORS"))))))

(def eigen-strings ["The overlap eigenvalues:" "The energy eigenvalues:"])
(def eigen-keys [:overlap_eigen :energy_eigen])

(defn eigen-parser [lines]
  (bulk-parser lines nnext empty-str read-rest))

(def ham-comp-key-vals '[("s: overlap" :s) ("t: kinetic" :t) ("vna: neutral atom" :vna) ("vxc: exchange/correlation" :vxc) ("vxc_1c: one-center exchange/correlation" :vxc_1c) ("Total Hamiltonian: h" :total_h)])

(defn ham-comp-parser [lines]
  (bulk-parser lines nnext empty-str (str-fn (list (read-rest str)))))

(defn occupy-parser [#^String string lines]
  (bulk-parser lines next empty-str (str-fn (read-rest (str-skip str string)))))

(defn shell-charges-parser [lines]
  (bulk-parser lines nnext contains= #(list (drop 3 (read-rest %)))))

(defn sum-charges-parser [lines]
  (bulk-parser lines nnext contains= (str-fn (list (read-next str)))))

(defn make-parsers [#^String str key item-parser]
  (list #(string-starts str %)
	#(let [#^String one (first %)]
	   (try
	    (item-parser (str-skip (. one trim) str))
	    (catch RuntimeException e (first %))))
	key))

(def gulp-parsers `(~@(map #(make-parsers (first %) (second %) (nth % 2 read-next)) gulp-key-vals1)
		    ~[#(string-starts "Final cell parameters and derivatives :" %)
		      (fn [lines] (bulk-parser lines #(nthnext % 3) dash-str (comp list second read-rest))) :final-lat-params]))

(def fireball-parsers `(
	       ~(make-parsers "FIREBALL RUNTIME :" :runtime read-string)
	       ~(make-parsers "Time step =" :scfsteps #(nth (read-rest %) 4))
	       ~(make-parsers "T_want =" :twant #(first (read-rest %)))
	       ~@(map #(make-parsers (first %) (second %) read-next) fb-key-vals1)
	       ~@(map #(make-parsers (first %) (second %) read-rest) fb-key-vals3)
	       ~(list #(string-starts "Corrected LATTICE VECTORS" %) lattice-parser :lat-vectors)
;	       ~(list #(string-contains "fermi ioccupy_k" %) #(occupy-parser "ioccupy =" %) :ioccupy)
;	       ~(list #(string-contains "fermi foccupy" %) #(occupy-parser "foccupy =" %) :foccupy)
;	       ~@(map (fn [str key] (list #(string-starts str %) eigen-parser key)) eigen-strings eigen-keys)
;	       ~@(map (fn [comp] (list #(string-starts (first comp) %)
;				       #(force-comp-parser (second comp) %)
;				       (keyword (str (second comp))))) force-comp)
;	       ~@(map (fn [comp]
;			(list #(string-starts (first comp) %) ham-comp-parser (second comp))) ham-comp-key-vals)
;disable by default to save memory
;	       ~(list #(let [#^String one (first %)]
;			 (. one contains "charges for each atom:")) sum-charges-parser :sum_charges)
;	       ~(list #(string-starts "Atom #    Type    Shells   Charges" %) shell-charges-parser :shell_charges)
     ))

(defn remove-scfsteps [hmap]
  "Given a parsed logfile, finds all quantities that are
   written at each timestep, and remove all but the last timestep."
  (let [scf (:scfsteps hmap)
	num (reduce + scf)
	scfkeys (filter #(= num (count (get hmap %))) (keys hmap))
	scfvals (map #(last-widths (get hmap %) scf) scfkeys)]
    (merge hmap (zipmap scfkeys scfvals))))

(defn get-lattice-params
  "Note: gulp prints out lattice parameters differently when doing md and optimization.  Use this for md."
  [log]
  (let [ks [:lat-param-a :lat-param-b :lat-param-c :lat-param-al :lat-param-be :lat-param-ga]
	f #(/ (* % Math/PI) 180), p (map #(get log %) ks)]
    (transpose (concat (take 3 p) (map #(map f %) (drop 3 p))))))

;usage
;(def logfile (parse-fortran "path-to-f.log" fireball-parsers))
;(keys logfile)
;(get logfile :any-key)

(defn parse-zeta [filename]
  "Parses Zeta.dat"
  ;currently assumes all atoms and 1 1 1 dyn vector
  (let [zeta (read-table filename)
	zetap (partition (second (first zeta))
			 (next zeta))]
    (transpose zetap)))