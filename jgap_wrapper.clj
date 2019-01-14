(ns jgap-wrapper (:import (org.jgap Configuration Chromosome FitnessFunction Genotype)
			(org.jgap.impl DefaultConfiguration DoubleGene)))

(defn get-conf
  ([] (get-conf true true))
  ([elite constant]
  (Configuration/reset)
  ;apparently jgap only allows one configuration, so we can't evolve in parallel
  (let [conf (new DefaultConfiguration)]
    (. conf (setPreservFittestIndividual elite))
    (. conf (setKeepPopulationSizeConstant constant))
    conf)))

(defn get-pop [#^Configuration conf fit-func genes size]
  (. conf (setFitnessFunction fit-func))
  (. conf (setSampleChromosome (new Chromosome conf genes)))
  (. conf (setPopulationSize size))
  (Genotype/randomInitialGenotype conf))

(defn make-genes [conf mins maxes]
  (into-array (map #(new DoubleGene conf %1 %2) mins maxes)))

(defn get-alleles [best]
  (map #(. %1 getAllele) (. best getGenes)))

(defn evolve [#^FitnessFunction fit-func mins maxes size maxiter]
  (let [conf (get-conf)
	genes (make-genes conf mins maxes)
	#^Genotype pop (get-pop conf fit-func genes size)]
    (dotimes [i maxiter]
      (. pop evolve))
    (. pop getFittestChromosome)))

(defn two-param-fit-func [func]
  (proxy [org.jgap.FitnessFunction] []
    (evaluate [chromo]
	      (let [[one two] (get-alleles chromo)
		    result (func one (+ one two))]
		(/ (. result size) (* two two))))))