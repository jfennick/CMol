(ns geometry (:use [utils :only [map-to-seq transpose]]
		   [empirical-data :only [atomic-mass]]
                   [clojure.contrib.generic.math-functions :only [rint]]))

;(def ^:const PI 3.141592653589793) ;clojure 1.3 version
(def PI 3.141592653589793) ;clojure 1.2 version

(defn magnitude [coords]
  "Calculates the sqrt of the sum of the squares of coords"
  (Math/sqrt (reduce + (map #(* % %) coords))))

(defn normalize-mag [lst]
  "Divides each element of lst by the magnitude of lst, such that the new magnitude is equal to 1."
  (let [mag (magnitude lst)]
    (if (= 0.0 mag)
      (recur (take (count lst) (repeat 1)))
      (map #(/ % mag) lst))))

(defn normalize-sum [lst]
  "Divides each element of lst by the sum of lst, such that the new sum is equal to 1."
  (let [sum (reduce + lst)]
    (if (= 0.0 sum)
      (recur (take (count lst) (repeat 1)))
      (map #(/ % sum) lst))))

(defn euclidean [x y]
  "Returns the euclidean distance between two coordinates, represented as seqs of numbers."
  (magnitude (map - x y)))

(defn get-com [obj]
  "Returns the center of mass coordinate of an atom, or a nested hmap of atoms"
  (if-let [c (:coordinates obj)] c
    (let [coords (flatten (map-to-seq obj))
	  ;uniques (set coords);assume double-counting is okay
	  ;perform a weighted sum of the coordinates
	  result (reduce #(list (+ (first %1) (atomic-mass %2))
				(map (fn [x1 x2]
				       (+ x1 (* x2 (atomic-mass %2))))
				     (second %1) (:coordinates %2)))
			 [0 [0 0 0]] coords)]
      ;divide result by total weight and return
      (map #(/ % (first result)) (second result)))))

(defmacro icomp [expr coll]
  "helper macro for moment of inertia tensor components. captures the symbols x y z"
  `(reduce + (map (fn [[~'a ~'b ~'c]] (let [~'x (float ~'a) ~'y (float ~'b) ~'z (float ~'c)] ~expr)) ~coll)))

(defn get-moi [obj]
  "Calculates the moment of inertia tensor of an atom, or a nested hmap of atoms."
  (let [temp (map :coordinates (flatten (map-to-seq obj)))
	com (get-com obj)
	coords (map #(map - com %) temp)
	xx (icomp (+ (* y y) (* z z)) coords)
	yy (icomp (+ (* x x) (* z z)) coords)
	zz (icomp (+ (* x x) (* y y)) coords)
	xy (- (icomp (* x y) coords))
	xz (- (icomp (* x z) coords))
	yz (- (icomp (* y z) coords))
	yx xy, zx xz, zy yz];tensor is symmetric; only 6 independent components
    [[xx xy xz] [yx yy yz] [zx zy zz]]))

(defn com-dr [one two]
  "returns the magnitude of the distance between the center of masses of objects one and two."
  (let [c1 (get-com one) c2 (get-com two)]
       (magnitude (map - c1 c2))))

(defn com-angle [one two three]
  "returns the angle between the center of masses of objects one and two three."
  (let [c1 (get-com one) c2 (get-com two) c3 (get-com three)
	d1 (map - c1 c2) d2 (map - c3 c2)
	dot (reduce + (map * d1 d2))]
    (Math/acos (/ dot (* (magnitude d1) (magnitude d2))))))

(defn vec-times-scalar [scalar vect]
  (map * (repeat scalar) vect))

(defn dot-product [x y]
  (reduce + (map * x y)))

(defn cross-product [pnt1 pnt2]
  "pnt1 and pnt2 are 3-tuples."
  (let [[x1 y1 z1] pnt1
	[x2 y2 z2] pnt2]
    [(- (* y1 z2) (* z1 y2))
     (- (* z1 x2) (* x1 z2))
     (- (* x1 y2) (* y1 x2))]))

(defn point-plane-distance [point plane]
  "Returns the (signed) distance.  Positive distance means point and normal
vector are on the same side."
  (let [[location normal] plane
	d (- (dot-product location normal))]
    (+ d (dot-product normal point))))

(defn mat-vect-mult [mat vect]
  "Performs matrix vector multplication.
   mat must be a seq of row vectors, not column vectors"
  (map #(reduce + (map * % vect)) mat))

(defn mat-mat-mult [mat1 mat2]
  "Performs matrix matrix multiplication.
   mat1 and mat2 are seqs of row vectors, not column vectors."
  (let [t (transpose mat2)]
    (map #(mat-vect-mult t %) mat1)))

(defn adjugate3 [matrix]
  (let [[[m11 m12 m13]
	 [m21 m22 m23]
	 [m31 m32 m33]] matrix]
    [[(- (* m22 m33) (* m32 m23)) (- (* m31 m23) (* m21 m33)) (- (* m21 m32) (* m31 m22))]
     [(- (* m32 m13) (* m12 m33)) (- (* m11 m33) (* m31 m13)) (- (* m31 m12) (* m11 m32))]
     [(- (* m12 m23) (* m22 m13)) (- (* m21 m13) (* m11 m23)) (- (* m11 m22) (* m21 m12))]]))

(defn determinant3 [matrix]
  (dot-product (first matrix) (first (adjugate3 matrix))))

(defn invert-matrix3 [matrix]
  (let [det (double (determinant3 matrix)), f #(/ % det)]
    (map #(map f %) (adjugate3 (transpose matrix)))))

(defn get-rotation-matrix [vectr angle]
  "Returns the Rodrigues' rotation matrix associated with the vector and angle.
           | 00  01  02 |
     rmat: | 10  11  12 |
           | 20  21  22 |"
  (let [[a b c] (normalize-mag vectr)
	x (float a) y (float b) z (float c) one (int 1);primitive locals = 10X faster!
	sina (Math/sin angle) cosa (Math/cos angle) dcosa (- one cosa)
	r00 (+ (* x x) (* (- one (* x x)) cosa))
	r01 (- (* x y dcosa) (* z sina))
	r02 (+ (* x z dcosa) (* y sina))
	r10 (+ (* x y dcosa) (* z sina))
	r11 (+ (* y y) (* (- one (* y y)) cosa))
	r12 (- (* y z dcosa) (* x sina))
	r20 (- (* x z dcosa) (* y sina))
	r21 (+ (* y z dcosa) (* x sina))
	r22 (+ (* z z) (* (- one (* z z)) cosa))]
    [[r00 r01 r02] [r10 r11 r12] [r20 r21 r22]]))

(defn mean [lst]
  (/ (reduce + lst) (count lst)))

(defn parameters->vectors
  "Converts a unit cell specified in lattice parameters into lattice vectors.
   Assumes the a vector is aligned with the cartesian x axis, and the b vector is in the cartesian xy plane.
   Formula found at http://en.wikipedia.org/wiki/Fractional_coordinates
   Returns a seq of column vectors, [a b c].
   where alpha is the angle between axis b and c
       beta is the angle between axis a and c
       gamma is the angle between axis a and b"
  [a b c alpha beta gamma]
  (let [cosa (Math/cos alpha);TODO: add primitive locals for performance
	cosb (Math/cos beta), sinb (Math/sin beta)
	cosg (Math/cos gamma), sing (Math/sin gamma)
	b1 (* b cosg), b2 (* b sing), c1 (* c cosb)
	cosbg (* cosb cosg), c2 (* c (/ (- cosa cosbg) sing))
	bca (* b c (Math/sin alpha))
	bc (- (* b1 c2) (* b2 c1)), c3 (Math/sqrt (Math/abs (/ (- (* bca bca) (* bc bc)) (* b b))))]
    [[a 0 0]
     [b1 b2 0]
     [c1 c2 c3]]))

(defn find-angle [a b]
  "Returns the angle between the vectors a and b."
  (Math/acos (/ (dot-product a b) (* (magnitude a) (magnitude b)))))

(defn align-lattice*
  "Rotates an arbitrary lattice so that
   the a vector is aligned with the cartesian x axis, and the b vector is in the cartesian xy plane."
  [lattice]
  (let [[a1 b1 c1] lattice
	angle1 (find-angle a1 [1 0 0])
	vector1 (cross-product a1 [1 0 0])
	matrix1 (get-rotation-matrix vector1 angle1)
	;dot-product is unsigned, so check if we rotated in the correct direction
	lattice2 (if (> 1E-6 (find-angle (mat-vect-mult matrix1 a1) [1 0 0]))
		   (map #(mat-vect-mult matrix1 %) lattice)
		   (let [m (get-rotation-matrix vector1 (- angle1))]
		     (map #(mat-vect-mult m %) lattice)))
	[a2 b2 c2] lattice2;Note: a2 b2 c2 are each still vectors
	xydist (point-plane-distance b2 [[0 0 0] [0 0 1]])
	b2xy (map - b2 [0 0 xydist])];the vector b2 projected onto the xy plane
    (if (zero? (magnitude b2xy));if b2 is already parallel with the cartesian z axis
      (let [pi2 (/ PI 2), angle (if (neg? xydist) pi2 (- pi2));same check as above
	    m (get-rotation-matrix [1 0 0] angle)]
	(map #(mat-vect-mult m %) lattice2))
      (let [angle2 (find-angle b2 b2xy), vector2 [1 0 0]
	    matrix2 (get-rotation-matrix vector2 angle2)]
	(if (> 1E-6 (find-angle (mat-vect-mult matrix2 b2) b2xy));same check as above
	  (map #(mat-vect-mult matrix2 %) lattice2)
	  (let [m (get-rotation-matrix vector2 (- angle2))]
	    (map #(mat-vect-mult m %) lattice2)))))))

(defn align-lattice [lattice]
  "align-lattice* works, but may give vectors that are slightly off numerically.
   aligning twice (or more) should give better results."
  ;Note:nth is indexed from zero, but the first iteration is x, not (f x).
  (nth (iterate align-lattice* lattice) 2))

(defn vectors->parameters
  "Converts a unit cell specified in lattice vectors to lattice parameters.
   Assumes the a vector is aligned with the cartesian x axis, and the b vector is in the cartesian xy plane."
  [lattice]
  (let [[a b c] lattice]
    (concat (map magnitude lattice)
	    (map find-angle [b c a] [c a b]))))

(defmacro tolerance? [value tol]
  "This function is a two-tuple predicate.
   It is TRUE if abs(val) < abs(tol), FALSE otherwise.

Usage: (tolerance? 2 1.0E-8) => false
       (tolerance? 1.0E-10 1.0E-8) => true"
   `(< (Math/abs ~value) (Math/abs ~tol)))
