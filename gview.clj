(ns gview (:import (javax.swing JFrame JScrollPane JTree))
    (:use clojure.contrib.seq geometry reflection-utils))

(defn make-node
  ([parent obj] (make-node parent obj pr-str))
  ([parent obj pr-fn]
     (cond (nil? obj)
	   (make-node parent [] pr-fn)
	   
	   (instance? java.awt.Container obj)
	   (make-node parent obj pr-fn #'identity #(. % getComponents))
	   
	   (instance? java.util.ArrayList obj)
	   (make-node parent (seq obj) pr-fn)
	   
	   (or (coll? obj) (seq? obj) (map? obj))
	   (make-node parent obj pr-fn #'identity #'seq)
	   
	   true
	   (make-node parent obj pr-fn #(empty? (get-private-field-objects %))
		      #(get-private-field-objects %))
	   ))
  ([parent obj pr-fn ischildfn getchildfn]
     (proxy [javax.swing.tree.DefaultMutableTreeNode] [obj]
       (toString []          (pr-fn obj))
       (getAllowsChildren [] (ischildfn obj))
       (getChildAt [i]       (make-node this (nth (getchildfn obj) i)))
       (getChildCount []     (count (getchildfn obj)))
       (getIndex [n]         -1)
       (getParent []         parent)
       (isLeaf []            (not (ischildfn obj))))))

(defn tree-view
  "Creates a hierarchical viewer in a new JFrame. Calls func with jtree on selection changed."
  ([obj func] (tree-view obj func pr-str))
  ([obj func pr-fn]
     (let [jtree (JTree. (make-node nil obj pr-fn))]
       (. jtree (addTreeSelectionListener
		 (proxy [javax.swing.event.TreeSelectionListener] []
		   (valueChanged [e] (func jtree)))))
       (doto (JFrame.)
	 (.add (JScrollPane. jtree))
	 (.setTitle (str "tree-view: " (.getName (class obj))))
	 (.setDefaultCloseOperation JFrame/DISPOSE_ON_CLOSE)
	 (.pack)
	 (.setVisible true)))))

(defn get-selected-objects [jtree]
  "Returns a seq of the user objects at the selected leaves of a jtree."
  (map #(.. % getLastPathComponent getUserObject)
       (. jtree getSelectionPaths)))

(defn tsg [canvas coords objs]
  "Tree-subset-grapher For use with tree-view and viewer-3d. This function will highlight the currently selected set of atoms."
  (fn [jtree]
    (let [objects (get-selected-objects jtree)
	  atoms (into #{} (flatten (set-to-seq objects)))]
      (dosync (ref-set coords (map get-com atoms)))
      (dosync (ref-set objs objects))
      (. canvas repaint))))

'(tree-view obj (let [[f c] (viewer-3d (plot-coords coords))]
	      (tsg c coords)))

#_(defn slide-view [n func]
  "Creates a linear viewer in a new JFrame with n marks.  Calls func with getValue when a new value is selected."
     (let [jframe (new javax.swing.JFrame)
	   slider (new javax.swing.JSlider)
	   change (proxy [javax.swing.event.ChangeListener] []
		    (stateChanged [e]
				  (func (. slider getValue))))]
       (. slider (addChangeListener change))
       (. slider (setMinimum 1))
       (. slider (setMaximum (count objs)))
       (. slider (setMinorTickSpacing 1))
       (. slider (setBounds 0 0 600 30))
       (. jframe (setLayout nil))
       (. jframe (add slider))
       (. jframe (setVisible true))))

#_(defn plot-zeta [canvas coords zeta]
  "For use with slide-view and viewer-3d. This function will plot the currently selected mode."
  (fn [n]
    (let [couplings (map #(last (nth % (dec n))) zeta)
	  ;partition couplings with widths
	  ;associate with atomic coordinates & normalize to get colors
	  ]
      (dosync (ref-set coords new-coords))
      (. canvas repaint))))