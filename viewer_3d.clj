(ns viewer-3d (:use bsptree geometry match-network [clojure.contrib combinatorics seq]
		    [penumbra window opengl])
    (:import (javax.swing JFrame) (java.awt Canvas Color Graphics)
	     (java.awt.event MouseAdapter MouseEvent MouseMotionListener MouseWheelEvent MouseWheelListener)))

(defn mouse-drag [[[dx dy] [x y]] button state]
  (assoc state
    :rot-x (+ (:rot-x state) dy)
    :rot-y (+ (:rot-y state) dx)))

(defn key-scale [key state]
  (cond (= "w" key)
	(update-in state [:scale] * 0.9)
	(= "s" key)
	(update-in state [:scale] / 0.9)
	true state))

(defn display [[delta time] state]
  (apply scale (repeat 3 (:scale state 1.0)))
  (rotate (:rot-x state) 1 0 0)
  (rotate (:rot-y state) 0 1 0)
  ((:plot-fn state identity) state))

(defn viewer-3d [callbacks init-state]
  (start (merge {:display display :mouse-drag mouse-drag :key-type key-scale} callbacks)
	 (merge {:rot-x 0 :rot-y 0} init-state)))

(defn get-box-lines [l u]
  "Returns the coordinates of the lines representing a wireframe box."
  (let [ooi (assoc l 2 (u 2)) oio (assoc l 1 (u 1))
	oii (assoc u 0 (l 0)) ioo (assoc l 0 (u 0))
	ioi (assoc u 1 (l 1)) iio (assoc u 2 (l 2))]
    [[l ioo] [l oio] [l ooi] [ooi ioi] [ooi oii] [oio iio] [oio oii] [oii u] [ioo ioi] [ioo iio] [ioi u] [iio u]]))

(defn get-bsp-lines [nodes]
  (let [lbounds (doall (map #(vec @(:lbounds %)) nodes))
	ubounds (doall (map #(vec @ (:ubounds %)) nodes))
	temp (map #(get-box-lines %1 %2) lbounds ubounds)]
    (into #{} (apply concat temp))))

(defn get-network-lines [db network matches]
  "network and matches are per match-network.  draws lines between elements of each match"
  (let [get-lines (fn get-lines [match net]
		    (map (fn [[[i1 i2] v]]
			   (let [[n1 n2] (:names v)]
			     [(if (map? n1) (get-lines (match i1) n1) [])
			      (if (set? n1) (get-lines (match i1) (get-wildcard db (match i1))) [])
			      (if (map? n2) (get-lines (match i2) n2) [])
			      (if (set? n2) (get-lines (match i2) (get-wildcard db (match i2))) [])
			      (if (:plot v true)
				{:1 (get-at (:a1 v) (match i1))
				 :2 (get-at (:a2 v) (match i2))}
				[])]))
			 net))]
	(->> (map #(set (map vals (flatten (get-lines % network)))) matches)
	     (apply concat))))

(defn plot-lines [state]
  (doseq [[a b] (:lines state)]
	(draw-lines (apply vertex a) (apply vertex b))))

(defn plot-coords [state]
  (doseq [c (:coords state)]
	(draw-points (apply vertex c))))

(defn plot-cell [coords matrix]
  (fn [state]
  	(plot-coords (assoc state :coords coords))
  	(plot-lines (assoc state :lines (map #(vector [0 0 0] %) matrix)))))

(defn plot-coms [objs]
  (plot-coords (map get-com @objs)))

(defn plot-pais [objs]
  (fn [state]
    (let [coms (map get-com @objs)
	  mois (map #(map normalize-mag (get-pai %)) @objs)]
      (doseq [com coms [a b c] mois]
		(plot-lines (assoc state :lines [[com (map + com a)] [com (map + com b)] [com (map + com c)]]))))))

(defn plot-objs [objs]
  (fn [state]
    ((juxt (plot-coms objs) (plot-pais objs)) state)))

(defn eval-fn [f argseq]
  "Evaluates the function f for every args in argseq. Returns a seq of (args (f args))."
  (map #(list % (apply f %)) argseq))

(defn get-func-coords [f xgrid ygrid]
  "xgrid and ygrid are [min max step].  Returns a plot-fn suitable for plotting the given function."
    (map (fn [[[x y] z]] [x y z])
	 (eval-fn f (cartesian-product (apply range xgrid) (apply range ygrid)))))