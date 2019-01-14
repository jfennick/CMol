(ns utils (:import java.io.File)
 (:use [clojure.contrib.io :only [file read-lines append-spit]]
   [clojure.contrib.combinatorics :only [combinations]]
   [clojure.string :only [join split trim trim-newline]]))

(def endline (System/getProperty "line.separator"))
(def slash (System/getProperty "file.separator"))
(def home (System/getProperty "user.home"))
(def tmp (System/getProperty "java.io.tmpdir"))
(defn nanotime [] (System/nanoTime))

(defn delete-dir [#^File dir]
  "Recursively deletes the contents of dir"
  (if (. dir isDirectory)
    (doseq [file (. dir listFiles)]
      (delete-dir file)))
  (. dir delete))

(defn temp-dir []
  "Returns a unique, temporary directory."
  (str tmp slash (nanotime)))

(defn empty-str [#^String s]
  (= "" (. s trim)))

(defn inter-cat [item coll]
  "Interposes item in coll, then concatenates into a string"
  (println  "Use clojure.string/join instead of inter-cat."))

(defn inter-cat-tree [items tree]
  "applies join recursively to tree"
  (if (nil? (next items))
    (join (first items) tree)
    (join (first items)
      (map #(inter-cat-tree (next items) %)
        tree))))

(defn read-all-string [#^String s]
  "reads all items in the string, not just the first one."
  (try
    (read-string (str "(" s ")"))
    (catch RuntimeException e s)))

(defn read-table [filename]
  "Reads a file with whitespace-delimited, tabular data."
  (map read-all-string (read-lines filename)))

(defn read-fortran-table [filename]
  "Reads a file with whitespace-delimited, tabular data."
  (->> (read-lines filename)
    (map #(.trim %))
    (remove #(= "" %))
    (map #(-> % (.replace "D" "E")))
    (map #(-> % (.replace ";" "")))
    (map #(if (and (.startsWith % "(") (not (.endsWith % ")")))
            (str % ")") %))
    (map read-all-string)))

(defn pad-string [#^String s length]
  "If str is less than length, appends spaces to str until it is length."
  (let [len (count s)]
    (if (>= len length) s
      (apply str (concat s (repeat len " "))))))

(defn partition-width
  "Partitions coll into variable width chunks."
  [coll widths]
  (loop [one coll
	 two widths
	 result nil]
    (if (or (nil? one) (nil? two))
      result
      (let [w (first two)
	    p (take w one)]
	(recur (nthnext one w) (next two) (concat result (list p)))))))

(defn last-widths [coll widths]
  "partitions coll into widths, returning the last element of each partition."
  (map last (partition-width coll widths)))

(defn flatten-n [n coll]
  "Like flatten, but only goes n levels deep."
  (if (= n 0)
    coll
    (recur (dec n) (apply concat (map #(if (sequential? %) % (list %)) coll)))))

(defn- unflatten* [tree coll]
  (loop [val-tree []
	 new-tree tree
	 new-coll coll]
    (if (nil? (first new-tree))
      [val-tree new-coll]
      (if (sequential? (first new-tree))
	(let [[a b] (unflatten* (first new-tree) new-coll)]
   (recur (conj val-tree a) (next new-tree) b))
	(recur (conj val-tree (first new-coll)) (next new-tree) (next new-coll))))))

(defn unflatten [tree coll]
  "Returns a new tree with the same shape as tree containing the elements of coll.
coll must be of length #leaves of tree"
  (first (unflatten* tree coll)))

(defn transpose [x]
  (if x (apply map list x) (throw (RuntimeException. "Can't transpose a nil matrix."))))
;If x is nil, the error will be "Wrong number of args passed to: core$map" and the backtrace will contain no user-level code. Oh no!

(defn set-to-seq [obj]
  "Recursively turns all sets into seqs"
  (if (set? obj)
    (map #(set-to-seq %) obj) obj))

(defn map-to-seq [obj]
  "Recursively turns all hmaps into seqs"
  (if-not (:coordinates obj)
    (map #(map-to-seq %) (vals obj)) obj))

(defn filter-user-code [trace]
  (let [java #{"LazySeq.java" "Cons.java" "Compiler.java" "RT.java" "Var.java" "AFn.java" "RestFn.java" "Thread.java"}
	clj #{"core.clj" "basic.clj" "NO_SOURCE_FILE"}
	code (clojure.set/union java clj)]
    (remove #(code (.getFileName %)) trace)))

(defn get-name-number [e]
  (let [s (filter-user-code (.getStackTrace e))]
    (map #(list (.getFileName %) (.getLineNumber %)) s)))

(defn find-cause [exception]
  (let [c (.getCause exception)]
    (if c (distinct (concat (get-name-number exception) (find-cause c)))
      (get-name-number exception))))
