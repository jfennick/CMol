(ns viewer-2d (:import (javax.swing JFrame) (java.awt Canvas Color Graphics)
		       (java.awt.event MouseAdapter MouseEvent MouseMotionListener MouseWheelListener)))

(defn eval-fn [f x-vals]
  "Evaluates the function f at every value of x. Returns a seq of (x (f x))."
  (map #(list % (f %)) x-vals))

(defn find-window [coords]
  "finds the minimum and maximum coordinates."
  (let [xs (map first @coords)
	ys (map second @coords)]
    (ref [(reduce min xs) (reduce max xs) (reduce min ys) (reduce max ys)])))

(defn x-pixel [x xmin xrange width]
  (int (* width (/ (- x xmin) xrange))))

(defn y-pixel [y ymin yrange height]
  (- height (int (* height (/ (- y ymin) yrange)))))

(defn clear-g [#^Graphics g width height]
  (. g setColor Color/black)
  (. g fillRect 0 0 width height)
  (. g setColor Color/white))

(defn make-coords [vals]
  "Returns a ref containing a seq of [index val]"
  (ref (map list (iterate inc 1) vals)))

(defn plot-points [coords]
  "Returns a plot-fn suitable for plotting coordinates.
   coords is a ref containing a seq of x-y pairs."
  (fn [g window width height]
    (let [[xmin xmax ymin ymax] @window
	  xrange (- xmax xmin) yrange (- ymax ymin)]
      (doseq [coord @coords]
	(. g fillOval
	   (x-pixel (first coord) xmin xrange width)
	   (y-pixel (second coord) ymin yrange height)
	   2 2)))))

(defn plot-lines [coords]
 "Returns a plot-fn suitable for plotting lines.
  coords is a ref containing a seq of x-y pairs."
 (fn [g window width height]
   (let [[xmin xmax ymin ymax] @window
	 xrange (- xmax xmin) yrange (- ymax ymin)]
     (doseq [coord (partition 2 1 @coords)]
       (let [a (first coord) b (second coord)]
	 (. g drawLine
	    (x-pixel (first a) xmin xrange width)
	    (y-pixel (second a) ymin yrange height)
	    (x-pixel (first b) xmin xrange width)
	    (y-pixel (second b) ymin yrange height)))))))

(defn plot-func [f]
  "Returns a plot-fn suitable for plotting the given function."
  (fn [g window width height]
    (let [[xmin xmax ymin ymax] @window
	  xrange (- xmax xmin) yrange (- ymax ymin)]
      (doseq [coord (partition 2 1 (eval-fn f (range xmin xmax (/ xrange width))))]
	(let [[a b] coord]
	  (. g drawLine
	     (x-pixel (first a) xmin xrange width)
	     (y-pixel (second a) ymin yrange height)
	     (x-pixel (first b) xmin xrange width)
	     (y-pixel (second b) ymin yrange height)))))))

(defn plot-logfile [log keywords coords-fn]
  "Given a fireball logfile, a seq of keywords, and a plot-fn that takes coords, returns a composition of the coords-fn's."
  (let [fns (map #(coords-fn (make-coords (% log))) keywords)]
    (fn [& args]
      (doseq [f fns]
	(apply f args)))))

(defn window-logfile [log keywords]
  "Given a fireball logfile and a seq of keywords, finds the window dimensions that encompass all coordinates."
  (let [windows (map #(find-window (make-coords (% log))) keywords)]
    (reduce #(let[[x1min x1max y1min y1max] @%1
		  [x2min x2max y2min y2max] @%2]
	       (ref [(min x1min x2min) (max x1max x2max) (min y1min y2min) (max y1max y2max)]))
	    windows)))

(defn viewer-2d [plot-fn window]
  "Opens a 2d viewer in a new JFrame.  plot-fn is fn of args (graphics window width height).
   window is a ref of (xmin xmax ymin ymax)"
  (let [xI (ref 0) yI (ref 0)
	down (ref false)
	jframe (JFrame.)
	canvas (proxy [Canvas] []
		 (paint [graphics] (let [width (. this getWidth) height (. this getHeight)]
				     (clear-g graphics width height)
				     (plot-fn graphics window width height))))
	mouse (proxy [MouseAdapter] []
		(mousePressed [e] (dosync (ref-set xI (. e getX))
					  (ref-set yI (. e getY))
					  (ref-set down true)))
		(mouseReleased [e] (dosync  (ref-set down false))))
	motion (proxy [MouseMotionListener] []
		 (mouseDragged [e] (if @down
				     (let [[xmin xmax ymin ymax] @window
					   xrange (- xmax xmin) yrange (- ymax ymin)
					   width (. this getWidth) height (. this getHeight)
					   dx (/ (* (- @xI (. e getX)) xrange) width)
					   dy (/ (* (- (. e getY) @yI) yrange) height)
					   new-window (map + @window [dx dx dy dy])]
				       (dosync (ref-set window new-window)
					       (ref-set xI (. e getX))
					       (ref-set yI (. e getY)))
				       (.repaint canvas))))
		 (mouseMoved [e]))
	wheel (proxy [MouseWheelListener] []
		(mouseWheelMoved [e] (let [[xmin xmax ymin ymax] @window
					   xp (* 0.10 (- xmax xmin)) yp (* 0.10 (- ymax ymin))
					   diff [xp (- xp) yp (- yp)]]
				       (dosync (if (< 0 (. e getWheelRotation))
						 (ref-set window (map + @window diff))
						 (ref-set window (map - @window diff))))
				       (.repaint canvas))))
	]
    (doto canvas
      (.addMouseListener mouse)
      (.addMouseMotionListener motion)
      (.addMouseWheelListener wheel))
    (doto jframe
      (.add canvas)
      (.setSize 640 480)
      (.setDefaultCloseOperation JFrame/DISPOSE_ON_CLOSE)
      (.pack)
      (.setVisible true))
    (list jframe canvas)))
