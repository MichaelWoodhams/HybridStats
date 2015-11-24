package hybridstats;
import java.io.PrintWriter;
import java.util.*;
import java.util.Map.Entry;

import mdwUtils.DoubleList;

import pal.tree.*;
import pal.misc.*;
import palExtensions.ExtRandom;
import palExtensions.IdGroupUtils;
import palExtensions.NeoSplitSystem;
import palExtensions.NeoSplitUtils;
import palExtensions.Split;

/**
 * 
 * 
 * @author woodhams
 *
 */

/*
 * TODO: Consider making this class Iterable
 */
public class SplitCounts implements IdGroup {
	private static final long serialVersionUID = -6924983916412736269L;
	private static final ExtRandom DEFAULT_RNG = new ExtRandom(4); // for shuffling for resolving ties in greedy consensus tree
	/*
	 * I spent considerable time trying to make HashMap<Split,Integer> work, but splits with the same
	 * hashCode() still got entered separately in the HashMap, so use two maps keyed by hex string instead
	 * TODO: fix this. This workaround is ugly.
	 */
	private HashMap<String,Integer> counts;
	private HashMap<String,Split> splits;
	private Vector<String> sortedSplits; // Splits sorted by frequency, ties randomly resolved. Cached result: is null until first needed.
	private Vector<Split> greedySplits; // Derived from sortedSplits, the splits in the greedy consensus tree. Is null until first needed.
	private Vector<Integer> greedySplitIndex; // where in sortedSplits the greedySplits occur.
	private IdGroup idg; // all splits must have the same IdGroup, to ensure consistent ordering of taxa.
	private int nTrees; // when adding splits a tree at a time, how many trees were added?
	private boolean splitsAddedOnlyViaTrees;
	private int nSplits; // total number of splits. Equal to the sum of values in 'counts'.
	private ExtRandom shuffler = DEFAULT_RNG;
	
	/*
	 * I see danger here: NeoSplitSystem is an IdGroup, for which we also have a constructor
	 */
	public SplitCounts(NeoSplitSystem system) {
		this();
		addSplitSystem__(system);
	}
	
	public SplitCounts(Forest forest) {
		this();
		this.addForest(forest);
	}
	
	public SplitCounts(IdGroup idGroup) {
		this();
		setIdGroup(idGroup);
	}
	
	/*
	 * If this constructor is used, 'id' will be set by first call to addSplitSystem.
	 */
	public SplitCounts() {
		counts = new HashMap<String,Integer>();
		splits = new HashMap<String,Split>();
		idg = null;
		nTrees = 0;
		nSplits =0;
		splitsAddedOnlyViaTrees = true;
		sortedSplits = null;
		greedySplits = null;
		greedySplitIndex = null;
	}
	
	public int totalNumberSplits() {
		return nSplits;
	}
	
	/*
	 * Sets idg if not already set.
	 * Throws error if try to set a non-equivalent idg.
	 */
	private void setIdGroup(IdGroup idGroup) {
		if (idg == null) idg = idGroup;
		if (!IdGroupUtils.equals(idg,idGroup)) throw new IllegalArgumentException("Tried to add split on different taxon set");
	}
	
	public void addForest(Forest forest) {
		for (Tree tree : forest) {
			this.addTree(tree);
		}
	}
	
	/**
	 * Using this method invalidates majorityRuleConsensusTree() method
	 * @param splitSys
	 */
	public void addSplitSystem(NeoSplitSystem splitSys) {
		splitsAddedOnlyViaTrees=false;
		addSplitSystem__(splitSys);
	}
	
	private void addSplitSystem__(NeoSplitSystem splitSys) {
		setIdGroup(splitSys.getIdGroup());
		for (Split split : splitSys) {
			nSplits++;
			String hex = split.toHexString();
			if (counts.containsKey(hex)) {
				counts.put(hex, counts.get(hex)+1);
			} else {
				counts.put(hex, 1);
				splits.put(hex,split);
			}
		}
	}
	
	public void addTree(Tree tree) {
		nTrees++;
		NeoSplitSystem splitSys = NeoSplitUtils.getSplits(tree);
		addSplitSystem__(splitSys);
	}
	
	/*
	 * If you want your own random number generator to control the shuffling (e.g. so you can set the seed).
	 */
	public void setRNG(ExtRandom rng) {
		this.shuffler = rng;
	}
	
	public Tree majorityRuleConsensusTree() {
		// Somebody used 'addSplitSystem()' method to supply splits to this count.
		if (!splitsAddedOnlyViaTrees) throw new RuntimeException("Can't determine consensus tree unless splits added only via trees");
		NeoSplitSystem consensus = new NeoSplitSystem(idg);
		int majority = nTrees/2+1;
		for (String hex : splits.keySet()) {
			if (counts.get(hex)>=majority) {
				consensus.add(splits.get(hex));
			}
		}
		return NeoSplitUtils.treeFromSplits(consensus);
	}
	
	/**
	 * Returns the sum of the Robinson Foulds distances between the 
	 * majority rules consensus tree and each tree in the collection.
	 * This does not actually require calculating the majority rule tree.
	 * @return
	 */
	public int sumRFtoMajRuleTree() {
		// Somebody used 'addSplitSystem()' method to supply splits to this count.
		if (!splitsAddedOnlyViaTrees) throw new RuntimeException("Can't determine consensus tree unless splits added only via trees");
		int sumDist = 0;
		for (Integer count : counts.values()) {
			sumDist += Math.min(count, nTrees-count);
		}
		return sumDist;
	}
	
	/**
	 * For each pair of splits: 
	 * If one or both occur <threshold> or fewer times, count 0.
	 * If they are compatible, count 0.
	 * Else reduce the count of each by <threshold>, multiply these counts, add to the sum.
	 * @param threshold
	 * @return
	 */
	public int weightedPairwiseSplitIncompatibility (int threshold) {
		int n = splits.size();
		int sum = 0;
		String[] keySet = new String[n];
		counts.keySet().toArray(keySet);
		for (int i=0; i<n-1; i++) {
			int count1 = counts.get(keySet[i]);
			if (count1<=threshold) continue;
			count1 -= threshold;
			for (int j=i+1; j<n; j++) {
				int count2 = counts.get(keySet[j]);
				if (count2>threshold && !splits.get(keySet[i]).compatible(splits.get(keySet[j]))) {
					sum += count1*(count2-threshold);
				}
			}
		}
		return sum;
	}
	
	/**
	 * Return the total number of incompatible split pairs in the split collection
	 * @return
	 */
	public int weightedPairwiseSplitIncompatibility() {
		return weightedPairwiseSplitIncompatibility(0);
	}

	public void printInternodeCertainties(PrintWriter out) {
		DoubleList<Split,Double> ic = getICs();
		for (int i=0; i<ic.size(); i++) {
			out.printf("IC = %f for split %s\n", ic.getB(i), ic.getA(i).toString());
		}
	}
	
	/**
	 * Returns the sum of Internode Certainties over a greedy consensus tree.
	 * @return
	 */
	public double treeCertainty() {
		double tc=0;
		DoubleList<Split,Double> ic = getICs();
		for (int i=0; i<ic.size(); i++) {
			tc += ic.getB(i);
		}
		return tc;
	}
	
	/**
	 * Returns the sum of Internode Certainty All over a greedy consensus tree.
	 * @return
	 */
	public double treeCertaintyAll(int threshold) {
		double tca=0;
		DoubleList<Split,Double> ic = getICAs(threshold);
		for (int i=0; i<ic.size(); i++) {
			tca += ic.getB(i);
		}
		return tca;
	}

	
	/* The comparator to use whenever we sort splits in this SplitCounts object */
	private final Comparator<String> frequencyComparator = new Comparator<String>() {
		public int compare(String s1, String s2) {
			return counts.get(s2)-counts.get(s1);
		}
	};
	/*
	 * Will redo the sort-by-frequency even if it has already been done.
	 * Clears greedySplits.
	 */
	private void resortSplits() {
		sortedSplits = new Vector<String>(counts.size()); // to be sorted by frequency
		for (String splitStr : counts.keySet()) {
			sortedSplits.add(splitStr);
		}
		shuffler.shuffle(sortedSplits);
		Collections.sort(sortedSplits,frequencyComparator);
		greedySplits = null;
		greedySplitIndex = null;
	}
	/*
	 * Ensure sortedSplits is populated. Don't recalculate if it is.
	 */
	private void sortSplits() {
		if (sortedSplits==null) resortSplits();
	}
	
	/*
	 * Does not recalculate if greedySplits are already cached.
	 */
	private void findGreedySplits() {
		if (greedySplits != null) return;
		sortSplits();
		int nGreedySplits = idg.getIdCount()-3;
		greedySplits = new Vector<Split>(nGreedySplits);
		greedySplitIndex = new Vector<Integer>(nGreedySplits);
		for (int i=0; i<sortedSplits.size() && greedySplits.size() < nGreedySplits; i++) {
			String splitStr = sortedSplits.elementAt(i);
			Split split = splits.get(splitStr);
			if (split.compatible(greedySplits)) {
				greedySplits.add(split);
				greedySplitIndex.add(i);
			}
		}
	}
	
	/*
	 * If recalculate = true, will recalculate greedy consensus tree with
	 * different random tie breaks, even if a greedy consensus tree is already cached.
	 * If recalculate = false, will recalculate only if no tree is already cached.
	 */
	public Tree greedyConsensusTree(boolean recalculate) {
		if (recalculate) resortSplits();
		findGreedySplits();
		return NeoSplitUtils.treeFromSplits(greedySplits);
	}
	
	/*
	 * Return a twin list of the greedy consensus tree splits and their internode certainties
	 */
	public DoubleList<Split,Double> getICs() {
		DoubleList<Split,Double> ic = new DoubleList<Split,Double>();
		List<DoubleList<Split,Integer>> greedyTreeConflicts = findConflictingSplitCounts(2, true);
		for (DoubleList<Split,Integer> splitList : greedyTreeConflicts) {
			ic.add(splitList.getA(0), internodeCertainty(splitList));
		}
		return ic;
	}

	/*
	 * Return a twin list of the greedy consensus tree splits and their ICA (internode certainty all) 
	 * scores, with a threshold for which incompatible splits are included in the ICA calculation
	 */
	public DoubleList<Split,Double> getICAs(int threshold) {
		DoubleList<Split,Double> ica = new DoubleList<Split,Double>();
		List<DoubleList<Split,Integer>> greedyTreeConflicts = findConflictingSplitCounts(threshold, false);
		for (DoubleList<Split,Integer> splitList : greedyTreeConflicts) {
			ica.add(splitList.getA(0), internodeCertainty(splitList));
		}
		return ica;
	}

	
	
	
	/*
	 * If thresholdIsLength is true, threshold = max number of splits to consider. (Use 2 to 
	 * get the classic Internode Certainty: main split plus strongest competitor.)
	 * If thresholdIsLength is false, threshold = min split weight for split to make the list.
	 */

	public List<DoubleList<Split,Integer>> findConflictingSplitCounts(int threshold, boolean thresholdIsLength) {
		findGreedySplits();
		List<DoubleList<Split,Integer>> results = new Vector<DoubleList<Split,Integer>>(greedySplits.size());
		// could be more efficient by storing where the greedy splits appear in the sorted list.
		for (int i : greedySplitIndex) {
			String splitStr = sortedSplits.elementAt(i);
			Split split = splits.get(splitStr);
			if (!thresholdIsLength && counts.get(splitStr)<threshold) break; // ignore splits with frequency below threshold
			DoubleList<Split,Integer> splitList = new DoubleList<Split,Integer>(); 
			splitList.add(split, counts.get(splitStr));
			// Only check splits after this one in sorted list: ones before this one are guaranteed to be
			// compatible, else this split would not be in the greedy list.
			for (int j=i+1; j<sortedSplits.size(); j++) {
				String otherSplitString = sortedSplits.elementAt(j);
				if (!thresholdIsLength && counts.get(otherSplitString)<threshold) break; // ignore conflicting splits with frequency below threshold
				Split otherSplit = splits.get(otherSplitString);
				if (!split.compatible(otherSplit)) {
					splitList.add(otherSplit,  counts.get(otherSplitString));
				} // if !compatible
				if (thresholdIsLength && splitList.size()==threshold) break; // have enough secondary splits now
			} // for otherSplit (j)
			results.add(splitList);		
		} // for i over sortedSplits
		return results;
	}
	
	/*
	 * Salichos Stamatakis and Rokas, MBE v31 p1261 (2014)
	 * If only two splits are listed, returns IC, the Internode Certainty.
	 * If more splits are listed, returns ICA, (IC All).
	 */
	private double internodeCertainty(DoubleList<Split,Integer> splits) {
		int sum = 0;
		int n = splits.size(); // number of splits under consideration
		for (int i=0; i<n; i++) sum += splits.getB(i);
		double ic = 1;
		if (n>1) {
			ic = Math.log(n);
			for (int i=0; i<n; i++) {
				double p = ((double)splits.getB(i))/sum;
				ic += p*Math.log(p);
			}
			ic /= Math.log(n); // convert from natural log to log base n
		}
		return ic;
	}

	
	/**
	 * 
	 * @return f[] where f[i]==c indicates there were c splits which were present in exactly i trees.
	 */
	public int[] countByFrequency() {
		int max=0;
		if (splitsAddedOnlyViaTrees) max = nTrees;
		else for (int n : counts.values()) max = Math.max(max, n);
		int[] freq = new int[max];
		for (int i : counts.values()) freq[i-1]++;
		return freq;
	}
	
	
	public int numUniqueCherries() {
		int count = 0;
		for (Split split : splits.values()) {
			if (split.sizeOfSmaller()==2) count++;
		}
		return count;
	}
	
	/**
	 * Crude display of the object contents:
	 * Splits as hex numbers, then count.
	 */
	public void hexDump(PrintWriter out) {
		for (Entry<String,Integer> entry : counts.entrySet()) {
			out.printf("%s: %d\n", entry.getKey(), entry.getValue());
		}
	}
	
	/**
	 * Replacement for hexDump
	 */
	public void tempDump(PrintWriter out) {
		for (String key : counts.keySet()) {
			out.printf("%s: %d\n", splits.get(key).toString(), counts.get(key));
		}
	}
	
	public int getCount(String key) {
		return (counts.containsKey(key)) ? counts.get(key) : 0;
	}
	
	public int getCount(Split split) {
		return getCount(split.toHexString());
	}
	
	public Split getSplit(String key) {
		return splits.get(key);
	}
	
	public Iterator<String> getHexIterator() {
		return splits.keySet().iterator();
	}
	
	public String[] getHexArray() {
		String[] array = new String[splits.size()];
		splits.keySet().toArray(array);
		return array;
	}
	
	public int numUniqueSplits() {
		return splits.size();
	}
	
	/*
	 * Methods to implement IdGroup, which just pass through to 'idg' member
	 */
	@Override
	public int getIdCount()                  { return idg.getIdCount(); }
	@Override
	public Identifier getIdentifier(int i)    {return idg.getIdentifier(i); }
	@Override
	public void setIdentifier(int i, Identifier id) { idg.setIdentifier(i, id); };
	@Override
	public int whichIdNumber(String name)     {return idg.whichIdNumber(name); }
}
