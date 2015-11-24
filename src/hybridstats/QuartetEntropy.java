package hybridstats;

import java.util.EnumMap;
import java.util.HashMap;
import java.util.Iterator;

import pal.distance.DistanceMatrix;
import pal.misc.IdGroup;
import pal.tree.Node;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.Tree;
import palExtensions.ExTreeUtils;
import palExtensions.Quadruple;
import palExtensions.Quartet.Topology;

/*
 * A static class.
 * 
 * Routines to calculate quartet entropy of a set of trees
 */

public class QuartetEntropy {
	public static double entropy(Forest forest) {
		int nTrees = forest.size();
		Quadruple[] quads = Quadruple.allQuadruples(forest.get(0));
		int nQuads = quads.length;
		HashMap<Quadruple,EnumMap<Topology,Integer>> topologyCounts = new HashMap<>(nQuads);
		for (Quadruple quad : quads) {
			EnumMap<Topology,Integer> map = new EnumMap<>(Topology.class);
			for (Topology topo : Topology.values()) map.put(topo, 0);
			topologyCounts.put(quad,map);
		}
		Iterator<Tree> iter = forest.iterator();
		IdGroup order = quads[0].getIdGroup();
		while (iter.hasNext()) {
			Tree tree = iter.next();
			// I should simply be able to do:
			// DistanceMatrix dist = new TreeDistanceMatrix(tree, order, true, 0);
			// but that is buggy.
			DistanceMatrix dist = ExTreeUtils.treeToDistanceMatrix(tree, order, true, 0);
			for (Quadruple quad : quads) {
				Topology topo = quad.whichQuartet(dist).getTopology();
				EnumMap<Topology,Integer> count = topologyCounts.get(quad);
				count.put(topo, count.get(topo)+1);
			}
		}
		// cache log-of-integer values:
		double[] log = new double[nTrees+1];
		for (int i=1; i<=nTrees; i++) log[i] = Math.log(i);
		double sum = 0;
		for (Quadruple quad : quads) {
			EnumMap<Topology,Integer> count = topologyCounts.get(quad);
			for (Topology topo : Topology.values()) {
				int n = count.get(topo); // number of times this topology was observed
				if (n>0) {
					if (topo == Topology.UNRESOLVED) throw new RuntimeException("haven't figured out how to deal with this yet");
					sum += n*log[n];
				}
			}
		}
		// same result as sum(-p_i log(p_i)) over all quartets where p_i = proportion of trees the quartet appears in.
		double entropy = nQuads*log[nTrees]-sum/nTrees; 
		/*
		 *  Normalize to max entropy, so return value is in range [0,1]
		 */
		return entropy / (nQuads*Math.log(3));
	}
	
	public static void test() {
		// Two simple tests: all trees the same (should give entropy 0)
		// and a forest of 3 different 4 taxon trees.
		Node leaf1, leaf2, leaf3, leaf4, int1, int2, root;
		Forest diverseForest = new Forest();
		leaf1 = new SimpleNode("A",1);
		leaf2 = new SimpleNode("B",1);
		leaf3 = new SimpleNode("C",1);
		leaf4 = new SimpleNode("D",1);
		int1  = new SimpleNode(); int1.addChild(leaf1); int1.addChild(leaf2); int1.setBranchLength(1); 
		int2  = new SimpleNode(); int2.addChild(leaf3); int2.addChild(leaf4); int2.setBranchLength(1);
		root  = new SimpleNode(); root.addChild(int1);  root.addChild(int2);
		Tree treeAB = new SimpleTree(root);
		diverseForest.add(treeAB);
		leaf1 = new SimpleNode(leaf1); // clone
		leaf2 = new SimpleNode(leaf2); 
		leaf3 = new SimpleNode(leaf3); 
		leaf4 = new SimpleNode(leaf4);
		int1  = new SimpleNode(); int1.addChild(leaf1); int1.addChild(leaf3); int1.setBranchLength(1); 
		int2  = new SimpleNode(); int2.addChild(leaf2); int2.addChild(leaf4); int2.setBranchLength(1);
		root  = new SimpleNode(); root.addChild(int1);  root.addChild(int2);
		Tree treeAC = new SimpleTree(root);
		diverseForest.add(treeAC);
		leaf1 = new SimpleNode(leaf1); // clone
		leaf2 = new SimpleNode(leaf2); 
		leaf3 = new SimpleNode(leaf3); 
		leaf4 = new SimpleNode(leaf4);
		int1  = new SimpleNode(); int1.addChild(leaf1); int1.addChild(leaf4); int1.setBranchLength(1); 
		int2  = new SimpleNode(); int2.addChild(leaf2); int2.addChild(leaf3); int2.setBranchLength(1);
		root  = new SimpleNode(); root.addChild(int1);  root.addChild(int2);
		Tree treeAD = new SimpleTree(root);
		diverseForest.add(treeAD);
		System.out.printf("Entropy of tree conflicting quartets = %f\n", entropy(diverseForest));
	}
}
