(ns potentials (:use geometry))
;todo: write a macro to cast params to doubles

(defn harmonic [[Ro K] objs]
  "1/2K(dr - Ro)^2"
  (let [dr (apply com-dr objs)]
    (* 0.5 K (Math/pow (- dr Ro) 2))))

(defn general [[A rho C m n] objs]
  "Ae^(-dr/rho)/dr^m - c/dr^n"
  (let [dr (apply com-dr objs)]
    (- (/ (* A (Math/exp (/ (- dr) rho)))
	  (Math/pow dr m))
       (/ C (Math/pow dr n)))))

(defn buckingham [[A rho C] objs]
  "Ae^(-dr/rho) - c/dr^6"
  (general [A rho C 6 0] objs))

(defn lennard-jones [[A B m n] objs]
  "A/dr^m - B/dr^n"
  (let [dr (apply com-dr objs)]
    (- (/ A (Math/pow dr (if m m 12)))
       (/ B (Math/pow dr (if n n 6))))))

(defn inverted-gauss [[A B r0] objs]
  "-Ae^(-B(dr-r0)^2)"
  (let [dr (apply com-dr objs)]
    (- (* A (Math/exp (- (* B (Math/pow (- dr r0) 2))))))))

(defn covalent-exponential [[A D r0] objs]
  "-De^(-A(dr-r0)^2/(2dr))"
  (let [dr (apply com-dr objs)]
    (- (* D (Math/exp (- (/ (* A (Math/pow (- dr r0) 2)) (* 2 dr))))))))

(defn fermi-dirac [[A B r0] objs]
  "A/(1+exp(Bdr))"
  (let [dr (apply com-dr objs)]
    (/ A (+ 1 (Math/exp (* B (- dr r0)))))))

(defn mei-davenport [[phi0 delta gamma r0] objs]
  "-phi0*(1 + delta*(dr/r0 - 1))*e^(-gamma*(dr/r0 - 1))"
  (let [dr (apply com-dr objs)
	dr1 (- (/ dr r0) 1)]
    (- (* phi0 (+ 1 (* delta dr1))
	  (Math/exp (- (* gamma dr1)))))))

(defn ljbuffered [[A B r0 m n] objs]
  "A/(dr + r0)^m - B/(dr + r0)^n"
  (let [dr (apply com-dr objs)]
    (- (/ A (Math/pow (+ dr r0) (if m m 12)))
       (/ B (Math/pow (+ dr r0) (if n n 6))))))

(defn polynomial [coefs objs]
  "sum of cn*dr^n"
  (let [dr (apply com-dr objs)]
    (reduce + (map #(* %1 (Math/pow dr %2)) coefs (iterate inc 0)))))

(defn rydberg [[A B r0] objs]
  "-A*[1+B*((dr/r0)-1)]*e^(-B*((dr/r0)-1))"
  (let [dr (apply com-dr objs)
	dr1 (- (/ dr r0) 1)]
    (- (* A (+ 1 (* B dr1))
	  (Math/exp (- (* B dr1)))))))
  