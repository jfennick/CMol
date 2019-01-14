(ns tree-utils)

(defn map-same [f coll]
  (into (empty coll) (map f coll))) 

(defn tree-apply [func only-coll recursive tree]
  "This function will call the given function on every element of the tree."
  (if only-coll
      (if (coll? tree)
	(let [result (map-same #(tree-apply func only-coll recursive %) tree)]
	  (if (or recursive (not-any? coll? tree))
	    (func result)
	    result))
	tree)
      (if (coll? tree)
          (map-same #(tree-apply func only-coll recursive %) tree)
          (func tree))))

(defn tree-map [func & trees]
  "This function will map the given function across corresponding elements of the trees.  Trees must have the same structure."
  (cond (not-any? coll? trees)
	(map-same func trees)
	
	(every? coll? trees)
	(map-same #(apply tree-map func %)
	     (partition (count trees)
			(apply interleave trees)))
	
	true
	trees;this should never happen
	))

(defn tree-reduce [func & trees]
  "Reduces func across the corresponding elements of the trees.  Trees must have the same structure."
  (tree-apply #(reduce func %) true false
	      (apply tree-map identity trees)))
