package hybridstats;

import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Vector;

import biojavaExtensions.GenericBlock;

/**
 * Info on which summary stats to calculate and report on.
 * Can include compound stats, which are new stats which are 
 * polynomial in the base stats. 
 * @author woodhams
 *
 */
public class SummaryStatParameters {
	public int[] siThresholds=null; // Split incompatibility thresholds. Don't use '0', that is just SI stat.
	public int[] rsThresholds=null; // rare splits thresholds.
	private Vector<CompoundStat> compoundStats;
	
	// Label strings:
	private static final String SI_THRESH = "split incompatibility thresholds";
	private static final String RS_THRESH = "rare splits thresholds";
	public final static Set<String> HYBRID_STATS_VALID_KEYS = new HashSet<String>(); 
	static {
		HYBRID_STATS_VALID_KEYS.add(SI_THRESH);
		HYBRID_STATS_VALID_KEYS.add(RS_THRESH);
	}

	// Defaults:
	private static final int[] SI_THRESH_DEF = new int[]{1,2};
	private static final int[] RS_THRESH_DEF = new int[]{1};
	public static final SummaryStatParameters DEFAULT = new SummaryStatParameters();
	
	/**
	 * NOTE! This method modifies 'block', removing the entries which it has processed.
	 * @param block
	 */
	public SummaryStatParameters(GenericBlock block) {
		// start with the defaults:
		this();
		// then overwrite them if required:
		if (block.hasKey(SI_THRESH)) { 
			siThresholds = parseIntSpecification(block.getValueTrimmed(SI_THRESH));
			block.removeField(SI_THRESH);
		}
		if (block.hasKey(RS_THRESH)) {
			rsThresholds = parseIntSpecification(block.getValueTrimmed(RS_THRESH));
			block.removeField(RS_THRESH);
		}
	}
	
	/**
	 * Generates the default SummaryStatsParameters: RS1, SI-1, SI-2.
	 */
	public SummaryStatParameters() {
		siThresholds = SI_THRESH_DEF;
		rsThresholds = RS_THRESH_DEF;
		compoundStats = new Vector<CompoundStat>();
	}
	
	public Iterator<CompoundStat> compoundStatIterator() {
		return compoundStats.iterator();
	}
	
	public void addCompoundStat(CompoundStat cStat) {
		compoundStats.add(cStat);
	}
	
	// Print headings for the non-standard stats (RS#, SI-# and compounds)
	// Tab-deliminted (friendly for R)
	public void printHeadings(PrintWriter out) {
		if (rsThresholds!=null) {
			for (int t : rsThresholds) {
				out.printf("\tRS%d",t);
			}
		}
		if (siThresholds!=null) {
			for (int t : siThresholds) {
				out.printf("\tSI-%d",t);
			}
		}
		for (CompoundStat compound: compoundStats) {
			out.printf("\t%s", compound.getName());
		}
	}
	
	// tab delimited
	public void printValues(PrintWriter out, HybridStats hStats) {
		if (rsThresholds != null) {
			for (int t : rsThresholds) {
				out.printf("\t%d",hStats.getCumulativeSplitCount(t));
			}
		}
		if (siThresholds != null) {
			for (int t : siThresholds) {
				out.printf("\t%d",hStats.getReducedSplitIncompatibility(t));
			}
		}
		for (CompoundStat compound: compoundStats) {
			out.printf("\t%7f", compound.evaluate(hStats));
		}
		out.println();
	}
	
	public SortedMap<String,Double> calculateCompoundStats(HybridStats hStats) {
		TreeMap<String,Double> map = new TreeMap<String,Double>(); 
		for (CompoundStat compound: compoundStats) {
			map.put(compound.getName(), compound.evaluate(hStats));
		}
		return map;
	}

	/**
	 * Throw an exception if thresholds are inconsistent with number of trees or taxa
	 * @param nTaxa
	 * @param nTrees
	 */
	public void rangeCheck(int nTrees) {
		int max=-1;
		for (int threshold : siThresholds) max = Math.max(max,threshold);
		// Maximum threshold for Split Incompatibility is just under half nTrees 
		if (2*max>=nTrees) throw new IllegalArgumentException("Split incompatibility threshold of "+max+" is too large for "+nTrees+" trees.");
		max=-1;
		for (int threshold : rsThresholds) max = Math.max(max,threshold);
		// Rare splits threshold maximum is nTrees (at which point we're counting all splits, not 'rare' ones.)
		if (max>=nTrees) throw new IllegalArgumentException("Rare split threshold of "+max+" is too large for "+nTrees+" trees.");
	}
	
	// Has similarity with parsing code from GridDimension
	// TODO: abstract this code to somewhere and share it with GridDimension
	private int[] parseIntSpecification(String str) {
		int[] intArray;
		int open  = str.indexOf('{');
		int close = str.indexOf('}');
		if (open>=0 && close >=0 && close > open) {
			// string contains /.*{.*}/
			String temp = str.substring(open+1, close);
			String[] tokens = temp.split("\\|");
			int n = tokens.length;
			intArray = new int[n];
			for (int i=0; i<n; i++) intArray[i] = Integer.valueOf(tokens[i]);
		} else {
			open = str.indexOf('(');
			close  = str.indexOf(')');
			if (open>=0 && close >=0 && close > open) {
				String temp = str.substring(open+1, close);
				String[] tokens = temp.split(":");
				int start = Integer.valueOf(tokens[0]);
				int stop  = Integer.valueOf(tokens[1]);
				int incr = (tokens.length==2) ? 1 : Integer.valueOf(tokens[2]);
				int n=(stop-start)/incr+1;
				intArray = new int[n];
				// remember incr is potentially negative, which causes problems for a naive 'for' loop implementation
				int x=start;
				for (int i=0; i<n; i++) {
					intArray[i]=x;
					x+=incr;
				}
			} else throw new RuntimeException("Could not parse threshold values from '"+str+"'");
		} 
		return intArray;
	}
}
