(ns bsptree (:use geometry clojure.contrib.seq))

(defstruct bspnode :parent :lbounds :ubounds :nodes :objects :recursive-size :maxobjects)

(defn make-bspnode
  ([parent lbounds ubounds maxobjects]
     (apply #'struct bspnode (concat (list (atom parent) (atom lbounds) (atom ubounds))
				     (map (fn [x] (atom nil)) (range 3)) (list maxobjects))))
  ([maxobjects]
     (apply #'struct bspnode (concat (map (fn [x] (atom nil)) (range 6)) (list maxobjects)))))

(defn inside-node? [coordinate node]
  "Returns node if all the components of the coordinate are within the bounds of the bspnode, else nil."
  (if (every? #'true? (map #(and (>= %1 %2) (< %1 %3))
			   coordinate (deref (:lbounds node)) (deref (:ubounds node))))
    node))

(defn inside-subnode? [coordinate node]
  "Return the subnode containing coordinate, else nil."
  (first (remove #'nil? (map #(inside-node? coordinate %1)
			     (vals (deref (:nodes node)))))))

(defn make-subnode [node object coordinate-func]
  "Adds a subnode to node."
  (let [midbounds (map #(/ (+ %1 %2) 2) (deref (:lbounds node)) (deref (:ubounds node)))
	subkey (map #(< %1 %2) (coordinate-func object) midbounds)
	sublbounds (map #(if %1 %2 %3) subkey (deref (:lbounds node)) midbounds)
	sububounds (map #(if %1 %2 %3) subkey midbounds (deref (:ubounds node)))
	subnode (make-bspnode node sublbounds sububounds 10)]
    (swap! (:nodes node) #'merge {subkey subnode})
    subnode))

(declare split-node)

(defn add-to-node
  "Adds object to node.  If node becomes too large, splits into subnodes."
  ([object node coordinate-func] (add-to-node object node coordinate-func false))
  ([object node coordinate-func update]
     (swap! (:objects node) #'conj object)
     (swap! (:recursive-size node) (fn [x y] (if (nil? x) 1 (inc x))) nil)
     (if update
       (loop [parent (deref (:parent node))]
	 (if (not (nil? parent))
	   (do (swap! (:recursive-size parent) (fn [x y] (if (nil? x) 1 (inc x))) nil)
	       (recur (deref (:parent parent)))))))
     (if (> (count (deref (:objects node))) (:maxobjects node))
       (trampoline split-node node coordinate-func))))

(defn split-node [node coordinate-func]
  "Subdivides a node and adds its objects to the new subnode(s)."
    (dorun (map (fn [object]
		  (if-let [subnode (inside-subnode? (coordinate-func object) node)]
	            ;add to existing subnode (none on first 'iteration')
		    (add-to-node object subnode coordinate-func false)
		    (add-to-node object (make-subnode node object coordinate-func) coordinate-func false)))
		(deref (:objects node))))
    (swap! (:objects node) (fn [old new] new) nil))

(defn not-on-edge? [node1 node2 root]
  "Returns true if node1 is a subnode of node2, and if all the bounds of node1 are 'nontrivially' inside node2, else nil"
  (let [lbounds1 (deref (:lbounds node1))
	ubounds1 (deref (:ubounds node1))
	lbounds2 (deref (:lbounds node2))
	ubounds2 (deref (:ubounds node2))
	lboundsr (deref (:lbounds root))
	uboundsr (deref (:ubounds root))
	l-ignore (map #(= %1 %2) lbounds1 lboundsr)
	u-ignore (map #(= %1 %2) ubounds1 uboundsr)
	;ignore dimension if edge coincides with root edge
	widths1 (map #(- %1 %2) ubounds1 lbounds1)]
    (every? #'true? (map #(and (or %6 (>= %1 (+ %3 %5)))
			       (or %7 (<= (+ %2 %5) %4)));adding the widths of the subnode is somewhat arbitrary, but it should work
				lbounds1 ubounds1 lbounds2 ubounds2 widths1 l-ignore u-ignore))))

;TODO (defn remove-from-node

(defn find-node [object node coordinate-func]
  "Returns the leaf node that contains coordinate. Call with root node."
  (if-let [subnode (inside-subnode? (coordinate-func object) node)]
    (find-node object subnode coordinate-func)
    (if-let [subnodes (deref (:nodes node))]
      ;if not a leaf, add a new subnode
      (trampoline make-subnode node object coordinate-func)
      ;else return node
      node)))

(defn min-node-node-distance [node1 node2]
  "Returns the minimum distance between two nodes."
  (let [lbounds1 (deref (:lbounds node1))
	ubounds1 (deref (:ubounds node1))
	lbounds2 (deref (:lbounds node2))
	ubounds2 (deref (:ubounds node2))]
    (Math/sqrt (reduce #'+ (map #(let [x1 (- %3 %2);lbound2 - ubound1
				       x2 (- %1 %4)];lbound1 - ubound2
				   (cond (> x1 0) (* x1 x1)
					 (> x2 0) (* x2 x2)
					 true 0))
				lbounds1 ubounds1 lbounds2 ubounds2)))))

(defn min-node-coord-distance [node coord]
  "Returns the minimum distance between a node and a coordinate."
    (let [l (deref (:lbounds node))
	  u (deref (:ubounds node))]
      (Math/sqrt (reduce #'+ (map #(let [x1 (- %3 %2)
					 x2 (- %1 %3)]
				     (cond (> x1 0) (* x1 x1)
					   (> x2 0) (* x2 x2)
					   true 0))
				  l u coord)))))

(defn find-nodes-rcutoff [coord root rcutoff]
  "Finds all leaf nodes that are within rcutoff of coord."
  (let [nodes (vals (deref (:nodes root)))
	filtered (filter #(<= (min-node-coord-distance %1 coord) rcutoff) nodes)]
    (concat (if (not nodes) (list root))
	    (map #(find-nodes-rcutoff coord %1 rcutoff) filtered))))

(defn get-count [zeros-count nonzeros k]
  (loop [i 0 num zeros-count non nonzeros]
    (if (or (nil? non) (>= num k))
      i
      (recur (inc i) (+ zeros-count (deref (:recursive-size (second (first non))))) (rest non)))))

(defn find-nodes-rk
  [node root k]
  (loop [nodes (vals (deref (:nodes root)))]
    (let [list (map #(list (min-node-node-distance node %1) %1) nodes)
	  zeros (filter #(= (first %1) 0) list)
	  nonzeros (filter #(not (= (first %1) 0)) list)
	  rofl (println nonzeros)
	  zeros-count (reduce #'+ (map #(deref (:recursive-size %1)) zeros))
	  newnodes (if (or (nil? nonzeros) (>= zeros-count k)) zeros
		       (merge zeros (take (get-count zeros-count nonzeros k) nonzeros)))
	  lol (println "HI1")
	  subnodes (map #(vals (deref (:nodes %1))) newnodes)
	  lol2 (println "HI2")]
      (if (every? #'nil? subnodes)
	newnodes;if all leaf nodes
	(recur (flatten (map #(if %1 %1 %2) subnodes newnodes)))))))

(defn get-leaf-nodes [node]
  (if-let [nodes (vals (deref (:nodes node)))]
    (do (if (not (nil? (deref (:objects node))))
	  ;there should only be objects in leaf nodes
	  (println "OH NODES!" (count (deref (:objects node)))))
	(map #(get-leaf-nodes %1) nodes))
    node))

(defn find-nodes-k
  "Starting at a leaf node, finds the node (a) that contains at least k objects, then finds the node (b) such that (a) is not on an edge of (b), then returns the leaves of (b)."
  [node k root]
  (loop [a node]
    (if (and (not (nil? (deref (:parent a))))
	     (< (deref (:recursive-size a)) k))
      (recur (deref (:parent a)))
      (loop [b a]
	(if (and (not (nil? (deref (:parent b))))
		 (not (not-on-edge? a b root)))
	  ;further optimization possible: if one or more edges of (a) coincide with the edges of the root,
	  ;only the other edges need to be inside
	  (recur (deref (:parent b)))
	  (trampoline get-leaf-nodes b))))))

(defn add-to-tree [root objs coord-func]
  (dorun (map #(let [node (find-node %1 root coord-func)]
		 (add-to-node %1 node coord-func false)) objs)))

(defn get-objects [nodes]
  (flatten (map #(deref (:objects %1)) nodes)))

(defn avg-num-objects [node]
  "Returns the average number of objects stored in the leaf nodes of node."
  (let [nums (map #(count (deref (:objects %1)))
		  (flatten (get-leaf-nodes node)))]
    (float (/ (reduce + nums) (count nums)))))

(defn filter-neighbors-r [objects obj rcutoff coord-func]
  "Returns a list of neighbors within rcutoff from the list of possible neighbors."
  (let [coord (coord-func obj)]
    (filter #(< (euclidean coord (coord-func %1)) rcutoff) objects)))

(defn filter-neighbors-k [obj objects k coordinate-func]
  "Returns a list of the k nearest neighbors from the list of possible neighbors."
  (let [coord (coordinate-func obj)]
    ;sorted-map is emulating a priority queue
    (loop [map (sorted-map)
	   objs objects]
      (if (nil? objs)
	(vals map)
	(let [d (euclidean coord (coordinate-func (first objs)))
	      last (first (last map))]
	  (if (< (count map) k)
	    (recur (merge map (sorted-map d (first objs))) (rest objs))
	    (if (and (> d 0) (< d (ffirst map)))
	      (recur (merge (dissoc map last) (sorted-map d (first objs))) (rest objs))
	      (recur map (rest objs)))))))))

(defn print-bsp [hmap]
  "bspnode contains circular references, so use this to print"
  (print-str "{" :l @(:lbounds hmap) :u @(:ubounds hmap)
	 :objs @(:objects hmap) :size @(:recursive-size hmap) :max (:maxobjects hmap)
	 :nodes (doall (map print-bsp (vals @(:nodes hmap)))) "}" ))

(defn get-objects-rcutoff [obj root coord-func rcutoff]
  (-> (coord-func obj)
      (find-nodes-rcutoff root rcutoff);get all leaf nodes within rcutoff of previous node
      flatten
      (#(remove nil? %))
      get-objects
      (filter-neighbors-r obj rcutoff coord-func);nodes may be within rcutoff, but objects in node may not
      ))

(defn get-all-objects [root]
  (-> (get-leaf-nodes root) flatten get-objects))

(defn filter-names [objs name]
  "Returns all objects that satisfy the given name"
  (remove #(empty? (re-seq  (re-pattern name) (:species %))) objs))

(defn get-bounds [objects coordinate-func]
  "Returns a 2-element list of the min and max coordinates"
  (let [[start & coords] (map coordinate-func objects)]
    (reduce (fn [[a b] c]
	      [(doall (map min a c)) (doall (map max b c))])
	    ;doall is necessary to avoid blowing the stack here
	    [start start] coords)))