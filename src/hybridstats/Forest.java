package hybridstats;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.RandomAccess;

import org.biojavax.bio.phylo.io.nexus.TreesBlock;
import org.biojavax.bio.phylo.io.nexus.TreesBlock.NewickTreeString;

import biojavaExtensions.NexusUtils;

import pal.tree.Tree;
import pal.tree.TreeParseException;
import palExtensions.ExTreeUtils;

/**
 * A collection of trees.
 * 
 * This is implemented as a wrapper around an ArrayList,
 * with the addition of some extra methods.
 * 
 * @author woodhams
 *
 */

/* Interfaces same as ArrayList */
public class Forest implements List<Tree>, RandomAccess, Cloneable, Serializable {
	private static final long serialVersionUID = 1L;
	private ArrayList<Tree> array;
	
	public Forest() {
		array = new ArrayList<Tree>();
	}
	
	public Forest(Tree[] treeArray) {
		this();
		for (Tree tree : treeArray) array.add(tree);
	}
	
	// TODO: (possibly) allow constraint that all trees must have same taxon set.
	// If so, throw error if input violates this, and ensure all trees have same Identifiers (and IdGroup?)
	public Forest(TreesBlock treesBlock) throws IOException, TreeParseException {
		this();
		for (Object obj : treesBlock.getTrees().values()) {
			String treeString = ((NewickTreeString)obj).getTreeString();
			int equalsIndex = treeString.indexOf('=');
			// transform 'T0=(a,(b,c))' into '(a,(b,c))'
			if (equalsIndex>=0) {
				treeString = treeString.substring(equalsIndex+1);
			}
			if (!treeString.endsWith(";")) treeString=treeString+";";
			try {
				array.add(ExTreeUtils.robustStringToTree(treeString));
			} catch (TreeParseException e) {
				// add useful information, then rethrow.
				System.err.printf("Error parsing tree from: %s\n", treeString);
				throw e;
			}
		}
	}
	
	public Forest(String filename) throws IOException, TreeParseException {
		this();
        FileReader fileReader = new FileReader(filename);
        BufferedReader bufferedReader = new BufferedReader(fileReader);
        String line = null;
        while ((line = bufferedReader.readLine()) != null) {
            array.add(ExTreeUtils.robustStringToTree(line));
        }
        bufferedReader.close();
	}
	
	public TreesBlock toTreesBlock(String blockComment) {
		int nTrees = array.size();
		String[] treeStrings = new String[nTrees];
		for (int i=0; i<nTrees; i++) {
			treeStrings[i] = array.get(i).toString();
		}			
		return NexusUtils.makeTreesBlock(treeStrings, null, blockComment);
	}

	/*
	 * From here on, methods required by the interfaces, which just pass straight through to 'array'.
	 */
	@Override
	public boolean add(Tree tree) { return array.add(tree); }
	@Override
	public void add(int n, Tree tree) { array.add(n,tree); }
	@Override
	public boolean addAll(Collection<? extends Tree> treeCollection) { return array.addAll(treeCollection); }
	@Override
	public boolean addAll(int n, Collection<? extends Tree> treeCollection) { 
		return array.addAll(n,treeCollection);
	}
	@Override
	public void clear() { array.clear(); }
	@Override
	public boolean contains(Object tree) { return array.contains(tree); } 
	@Override
	public boolean containsAll(Collection<?> treeCollection) { return array.containsAll(treeCollection); }
	@Override
	public Tree get(int n) { return array.get(n); }
	@Override
	public int indexOf(Object tree) { return array.indexOf(tree); }
	@Override
	public boolean isEmpty() { return array.isEmpty(); }
	@Override
	public Iterator<Tree> iterator() { return array.iterator(); }
	@Override
	public int lastIndexOf(Object tree) { return array.lastIndexOf(tree); }
	@Override
	public ListIterator<Tree> listIterator() { return array.listIterator(); }
	@Override
	public ListIterator<Tree> listIterator(int n) { return array.listIterator(n); }
	@Override
	public boolean remove(Object tree) { return array.remove(tree); }
	@Override
	public Tree remove(int n) { return array.remove(n); } 
	@Override
	public boolean removeAll(Collection<?> treeCollection) { return array.removeAll(treeCollection); }
	@Override
	public boolean retainAll(Collection<?> treeCollection) { return array.retainAll(treeCollection); }
	@Override
	public Tree set(int n, Tree tree) { return array.set(n, tree); }
	@Override
	public int size() { return array.size(); }
	@Override
	public List<Tree> subList(int from, int to) { return array.subList(from, to); }
	@Override
	public Object[] toArray() { return array.toArray(); }
	@Override
	public <T> T[] toArray(T[] destination) { return array.toArray(destination); }
}
