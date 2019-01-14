(ns clojuratica-config (:use utils clojuratica.clojuratica)
    (:import [com.wolfram.jlink MathLinkFactory]))

(System/setProperty "java.library.path" "/opt/Wolfram/MathematicaPlayer/7.0/SystemFiles/Libraries/Linux")

(def kernelpath "'/opt/Wolfram/MathematicaPlayer/7.0/Executables/MathKernel'")

(defn get-math [kernelpath]
  (let [kernel-link (MathLinkFactory/createKernelLink
		     (str "-linkmode launch -linkname " kernelpath))]
    (.discardAnswer kernel-link)
    (comp (get-parser kernel-link) (get-evaluator kernel-link))))

(defn sparse-eigensystem [math indices vals m n]
  "Given 3 vectors and the matrix dimensions, returns the eigenvalues and eigenvectors using mathematica."
  (math ["vals" vals "indices" indices "m" m "n" n]
    "Eigensystem[SparseArray[indices->N[vals],{m,n}]]"))

(defn invert-matrix [math mat]
  (math ["mat" mat] "Inverse[mat]"))

(defn eigenvalues [math mat]
  (math ["mat" mat] "Eigenvalues[mat]"))

(defn eigenvectors [math mat]
  (math ["mat" mat] "Eigenvectors[mat]"))

(defn eigensystem [math mat]
  (math ["mat" mat] "Eigensystem[mat]"))

(defn diagonalize [math mat]
  (math ["mat" mat] "vecs=Eigenvectors[mat];Inverse[Transpose[vecs]].mat.Transpose[vecs]"))

;(def m (read-table "/home/users/jfennick/matrix"))
;(dorun m)

;(def rows (map first m))
;(def cols (map second m))
;(def values (doall (map #(nth % 2) m)))
;(def indices (doall (map list rows cols)))