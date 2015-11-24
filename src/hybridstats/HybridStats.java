package hybridstats;

import java.io.PrintWriter;
import java.util.Arrays;

import pal.tree.Tree;

// TODO: The two letter abbreviations have become important and should be better integrated 
// with the code (e.g. defined all in one place in an array or static final definitions) 

// 2015-06-02 renamed "TS" (total splits) to "US" (unique splits) for consistency with "UC" (unique cherries)
// (which used to be "TC" until TreeCertainty got added.)
public class HybridStats {
	private int nTrees;
	private int nTaxa;
	int nSplits; // total number of splits. I.e. nTrees*(nTaxa-3) if all trees are fully resolved
	private double topoEntropy; // S1, "TE"
	private int[] cumTopoCounts; // can generate S2, S3 
	private int splitIncompat; // S4, "SI"
	private int[] reducedSplitIncompat; // S7, "SI-#"
	private int consensusDist; // S5, "DC"
	private int nCherries; // S9, "UC"
	private int splitsObs; // S10, "US"
	private double quartetEnt; // S11, "QE"
	private int[] cumulativeSplitCountByFreq; // can generate S12, "RS#"
	private double treeCertainty; // "TC"
	private double treeCertaintyAll; // "TCA"
	private SplitCounts splitCounts;
	
	public HybridStats(Tree[] trees) {
		this(new Forest(trees));
	}
	
	public HybridStats(Forest forest) {
		nTrees = forest.size();
		nTaxa = forest.get(0).getIdCount();
		splitCounts = new SplitCounts(forest);		
		nSplits=splitCounts.totalNumberSplits();
		TreeTopologyCounts topoCounts = new TreeTopologyCounts(forest); 
		topoEntropy=entropy(topoCounts.getCounts());
		cumTopoCounts = topoCounts.cumulativeCounts();
		splitIncompat = splitCounts.weightedPairwiseSplitIncompatibility();
		int nThreshold = nTrees/2;
		reducedSplitIncompat = new int[nThreshold];
		for (int i=0; i<nThreshold; i++) {
			reducedSplitIncompat[i] = splitCounts.weightedPairwiseSplitIncompatibility(i);
		}
		consensusDist = splitCounts.sumRFtoMajRuleTree();
		nCherries = splitCounts.numUniqueCherries();
		splitsObs = splitCounts.numUniqueSplits();
		treeCertainty = splitCounts.treeCertainty();
		treeCertaintyAll = splitCounts.treeCertaintyAll(0); // Possible TODO: use a suitable threshold instead of 0.
		quartetEnt = QuartetEntropy.entropy(forest);
		int[] splitCountByFreq = splitCounts.countByFrequency();
		cumulativeSplitCountByFreq = new int[splitCountByFreq.length];
		cumulativeSplitCountByFreq[0]=0;
		for (int i=1; i<splitCountByFreq.length; i++) {
			cumulativeSplitCountByFreq[i] = cumulativeSplitCountByFreq[i-1] + splitCountByFreq[i-1];
		}
	}
	
	public double getStatByName(String statName) {
		switch (statName) {
			case "1"  : return 1; // allows constant (intercept) term
			case "TE" : return topoEntropy;
			case "SI" : return splitIncompat;
			case "DC" : return consensusDist;
			case "UC" : return nCherries;
			case "US" : return splitsObs;
			case "QE" : return quartetEnt;
			case "TC" : return treeCertainty;
			case "TCA": return treeCertaintyAll;
			default : break;
		}
		// Only legitimate statNames left are "RS<int>" and "SI-<int>"
		if (statName.startsWith("RS")) {
			int index = Integer.valueOf(statName.substring(2));
			return cumulativeSplitCountByFreq[index];
		} else if (statName.startsWith("SI-")) {
			int index = Integer.valueOf(statName.substring(3));
			return reducedSplitIncompat[index];
		} else throw new IllegalArgumentException("Unrecognized stat name '"+statName+"'");
	}
	
	public void printHumanFriendly(PrintWriter out) {
		double maxTopoEntropy = nTrees * Math.log(nTrees);
		out.printf("(S1) Topology entropy = %f (max possible=%f)\n", topoEntropy, maxTopoEntropy);
		out.printf("(S2, S3) Cumulative counts of topologies: [%d", cumTopoCounts[0]);
		int i;
		for (i=1; i<cumTopoCounts.length && cumTopoCounts[i]-cumTopoCounts[i-1]>1; i++) {
			out.printf(" %d",  cumTopoCounts[i]);
		}
		if (i!=cumTopoCounts.length) {
			out.printf("... counts increase by one up to %d", cumTopoCounts.length);
		}
		out.print("]\n");
		//out.printf("(S2, S3) Cumulative counts of topologies: %s\n", Arrays.toString(cumTopoCounts));
		out.printf("(S4) Total pairwise split incompatibility = %d\n", splitIncompat);		
		out.printf("(S5) Sum diff Robinson Foulds distance to majority rule tree = %d (max possible = %d)\n", 
				consensusDist, (nTaxa-3)*nTrees);
		out.printf("(S9) Number of unique cherries = %d (max possible = %d)\n", nCherries, nTaxa*(nTaxa-1)/2);
		out.printf("(S10) Number of unique non-trivial splits observed = %d (c.f. %d for a single fully resolved tree, max %.0f)\n", 
				splitsObs, nTaxa-3, Math.min((nTaxa-3)*nTrees,Math.pow(2, nTaxa-1)-nTaxa-1));
		out.printf("(S11) Quartet entropy = %f\n", quartetEnt);
		int nPairs = nSplits*(nSplits-1)/2; // Number of pairwise split comparisons
		out.printf("(S12) Cumulative number of splits with a given frequency = %s\n", Arrays.toString(cumulativeSplitCountByFreq));
		out.printf("Tree certainty = %f\n", treeCertainty);
		out.printf("Tree certainty all = %f\n", treeCertaintyAll);
		out.print("(S7) Split incompatibilities beyond threshold:\nThresh.   Pairwise incompat.\n");
		int reduced=1;
		for (int threshold=0; threshold<nTrees/2 && reduced>0; threshold++) {
			reduced = reducedSplitIncompat[threshold];
			out.printf("%d/%d (%2.0f%%)      %d/%d\n",
					threshold, nTrees,
					100*((float)threshold)/nTrees,
					reduced,
					nPairs);
		}
	}
	
	
	public static void printRFriendlyHeadings(PrintWriter out, SummaryStatParameters stats) {
		out.print("TE\tSI\tDC\tUC\tUS\tQE\tTC\tTCA"); 
		stats.printHeadings(out);
		out.println();
	}
	
	/**
	 * Print stats in format suitable for input file for R (i.e. tab delimited table.)
	 * @param out 
	 * @param headers: if true, print column headings instead of data
	 * @param sitThresholds: Split Incompatibility with Threshold: list of thresholds.
	 */
	/*
	 * TODO: consider structure here. HybridStats and SummaryStatParameters are intertwined
	 * in an unseemly way. Perhaps move this method to SummaryStatParameters, taking
	 * a HybridStats as argument?
	 */
	public void printRFriendly(PrintWriter out, SummaryStatParameters stats, boolean headers) {
		if (headers) printRFriendlyHeadings(out,stats);
		out.printf("%7.3f\t%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%5.3f", 
				topoEntropy, splitIncompat, consensusDist, nCherries, splitsObs, quartetEnt, treeCertainty, treeCertaintyAll);
		stats.printValues(out, this);
	}
	
	public int getCumulativeSplitCount(int n) { return cumulativeSplitCountByFreq[n]; }
	public int getReducedSplitIncompatibility(int n) { return reducedSplitIncompat[n]; }
	public SplitCounts getSplitCounts() { return splitCounts; }
		
	/**
	 * Return the entropy of an observed multinomial distribution
	 * @param data
	 * @return
	 */
	/*
	 * Entropy = -sum_i (log(p_i^n_i))
	 * where n_i is number in state i, p_i = probability of state i
	 * and we take p_i = n_i/n_tot, which means
	 * Entropy = n_tot log(n_tot) - sum_i (n_i log(n_i))
	 */
	private static double entropy(int[] data) {
		int sum=0;
		double entropy = 0;
		for (int x : data) {
			sum += x;
			entropy -= x * Math.log(x);
		}
		entropy += sum * Math.log(sum);
		return entropy;
	}
}
