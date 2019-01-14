(ns reflection-utils)

(defn get-name [x]
  (. x getName))

(defn get-names [coll]
  (map #'get-name coll))

(defn get-private-fields [instance]
  (.. instance getClass getDeclaredFields))

(defn get-private-field-objects [instance]
  (remove nil? (map #(let [x (. %1 (setAccessible true))]
		       (if x (get instance)))
		    (get-private-fields instance))))

(defn get-private-methods [instance]
  (.. instance getClass getDeclaredMethods))

(defn get-fields [instance]
  (.. instance getClass getFields))

(defn get-methods [instance]
  (.. instance getClass getMethods))

(defn get-private-field [instance field-name]
  (. (doto (first (filter
		   (fn [x] (.. x getName (equals field-name)))
		   (.. instance getClass getDeclaredFields)))
       (.setAccessible true))
     (get instance)))

(defn get-private-field-names [instance field-names]
  (if-let [x (first field-names)]
	  (get-private-field-names (get-private-field instance x) (rest field-names))
	  instance))

(defn invoke-private-method [obj fn-name-string & args]
  (let [m (first (filter (fn [x] (.. x getName (equals fn-name-string)))
			 (.. obj getClass getDeclaredMethods)))]
    (. m (setAccessible true))
    (println "OBJ" obj)
    (println "ARGS" args)
    (. m (invoke obj args))))
