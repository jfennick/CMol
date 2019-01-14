(ns merge-utils)

(defn always-merge-with [f amap bmap]
  "Like merge-with, but always applies f"
  (loop [bkeys (keys bmap)
         bvals (vals bmap)
         result amap]
    (if (nil? bkeys)
      result
      (let [bkey (first bkeys)
            bval (first bvals)]
        (recur (next bkeys) (next bvals)
               (merge result {bkey (f (get amap bkey) bval)}))))))

(defn wrap-merge [amap bmap]
  "If key doesn't exist in amap, wraps val from bmap in [].
   Useful when incrementally merging a seq of vals into hashmaps when
the vals can be seqs themselves."
  (always-merge-with #(if (nil? %1) [%2] (conj %1 %2)) amap bmap))