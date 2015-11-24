package hybridstats;

import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

import pal.misc.IdGroup;
import pal.misc.Identifier;
import pal.tree.Tree;
import palExtensions.ExTreeUtils;
import palExtensions.IdGroupUtils;

public class TreeTopologyCounts implements IdGroup {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private HashMap<String,Integer> counts;
	private IdGroup idg;
	
	public TreeTopologyCounts() {
		counts = new HashMap<String,Integer>();
	}
	
	
	public TreeTopologyCounts(Forest forest) {
		this();
		this.addForest(forest);
	}

	public void addForest(Forest forest) {
		for (Tree tree : forest) {
			this.addTree(tree);
		}
	}
	
	public void addTree(Tree tree) {
		checkIdGroup(tree);
		// Puts tree in canonical order and omits branch lengths
		String str = ExTreeUtils.toTopologyString(tree);
		if (counts.containsKey(str)) {
			counts.put(str, counts.get(str)+1);
		} else {
			counts.put(str, 1);
		}
	}
	
	private void checkIdGroup(IdGroup newGroup) {
		if (idg==null) {
			idg = newGroup;
		} else {
			if (!IdGroupUtils.sameLabels(idg, newGroup)) {
				throw new IllegalArgumentException("Tried to add tree with incompatible IdGroup");
			}
		}
	}
	
	public int[] getCounts() {
		int n=counts.size();
		Integer[] data = new Integer[n];
		counts.values().toArray(data);
		Arrays.sort(data);
		int[] sorted = new int[n];
		for (int i=0; i<n; i++) {
			sorted[i] = data[i];
		}
		return sorted;
	}
	
	public int[] cumulativeCounts() {
		int[] sorted = getCounts();
		int n = sorted.length;
		int[] cumulative = new int[n];
		cumulative[0] = sorted[n-1];
		for (int i=1; i<n; i++) {
			cumulative[i] = cumulative[i-1]+sorted[n-i-1];
		}
		return cumulative;
	}
	
	public int getNumberUniqueTopologies() {
		return counts.size();
	}
	
	public void printSummary(Writer writer) throws IOException {
		for (Entry<String,Integer> entry : counts.entrySet()) {
			writer.write(String.format("%s: %d\n", entry.getKey(), entry.getValue()));
		}
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
