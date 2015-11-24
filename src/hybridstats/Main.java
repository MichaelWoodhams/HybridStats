package hybridstats;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.seq.io.ParseException;
import org.biojavax.bio.phylo.io.nexus.CharactersBlock;
import org.biojavax.bio.phylo.io.nexus.DataBlock;
import org.biojavax.bio.phylo.io.nexus.DistancesBlock;
import org.biojavax.bio.phylo.io.nexus.NexusBlock;
import org.biojavax.bio.phylo.io.nexus.NexusBlockParser;
import org.biojavax.bio.phylo.io.nexus.NexusComment;
import org.biojavax.bio.phylo.io.nexus.NexusFile;
import org.biojavax.bio.phylo.io.nexus.NexusFileBuilder;
import org.biojavax.bio.phylo.io.nexus.NexusFileFormat;
import org.biojavax.bio.phylo.io.nexus.TaxaBlock;
import org.biojavax.bio.phylo.io.nexus.TreesBlock;

import pal.tree.Tree;
import pal.tree.TreeParseException;
import palExtensions.ExTreeUtils;
import biojavaExtensions.NexusUtils;
import biojavaExtensions.UseableUnknownBlockParser;

public class Main {
	private static final String DEFAULT_IN_FILE = "output.nex";

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		PrintWriter out = new PrintWriter(System.out);
		boolean useLineageTrees = false;
		// Very crude command line parsing: '-l' = use lineage trees, optional input file name
		String filename = DEFAULT_IN_FILE; 
		if (args.length > 0) {
			if (args[0].equals("-l")) {
				useLineageTrees = true;
				if (args.length == 2) filename = args[1];
				if (args.length >2) throw new RuntimeException("Too many command line arguments");
			} else {
				filename = args[0];
				if (args.length>1) throw new RuntimeException("Too many command line arguments");
			}
		}
		Forest forest=null;
		try {
			forest = readTreesFromFile(filename,useLineageTrees);
		} catch (TreeParseException | IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			System.exit(1);
		}

		int nTaxa = forest.get(0).getIdCount();
		int nTrees = forest.size();
		HybridStats stats = new HybridStats(forest);
		out.printf("%d trees on %d taxa read from file %s\n\n",nTrees,nTaxa,filename);
		stats.printHumanFriendly(out);
		
		// and some misc. output that doesn't matter so much:
		if (!useLineageTrees) {
			Forest lineageForest=null;
			try {
				lineageForest = readTreesFromFile(filename,true);
			} catch (TreeParseException | IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
			SplitCounts lineageSplitCounts = new SplitCounts(lineageForest);
			out.printf("(L4) Lineage total pairwise split incompatibility = %d\n", lineageSplitCounts.weightedPairwiseSplitIncompatibility());
			out.printf("(L5) Lineage sum diff Robinson Foulds distance to majority rule tree = %d\n", lineageSplitCounts.sumRFtoMajRuleTree());
			out.printf("(L11) Lineage quartet entropy = %f\n", QuartetEntropy.entropy(lineageForest));
		}
		// These two are redundant recalculations already done in HybridStats object.
		SplitCounts splitCounts = stats.getSplitCounts();
		TreeTopologyCounts topoCounts = new TreeTopologyCounts(forest); 
		if (splitCounts.numUniqueSplits()<100) {
			out.println("\nCounts of splits:");
			splitCounts.tempDump(out);
		}
		if (topoCounts.getNumberUniqueTopologies()<100) {
			out.println("\nCounts of tree topologies:");
			try {
				topoCounts.printSummary(out);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		Tree majRule = splitCounts.majorityRuleConsensusTree();
		out.println("\nMajority rule consensus tree: "+ExTreeUtils.toTopologyString(majRule));
		out.println("\nGreedy consensus tree: "+ExTreeUtils.toTopologyString(splitCounts.greedyConsensusTree(false)));
		// Eclipse refuses to accept the following as valid, for reasons I cannot fathom
		// stats.PrintRFriendly(out, SummaryStatParameters.DEFAULT, true);
		stats.printRFriendly(out, new SummaryStatParameters(), true);
		
		splitCounts.printInternodeCertainties(out);

		out.close();
	}
	
	public static Forest readTreesFromFile(String filename, boolean useLineageTrees) throws TreeParseException, IOException {
		File file = new File(filename);
		NexusFileBuilder builder=new NexusFileBuilder();
		builder.setBlockParser(TaxaBlock.TAXA_BLOCK, new UseableUnknownBlockParser());
		builder.setBlockParser(CharactersBlock.CHARACTERS_BLOCK, new UseableUnknownBlockParser());
		builder.setBlockParser(DataBlock.DATA_BLOCK, new UseableUnknownBlockParser());
		builder.setBlockParser(DistancesBlock.DISTANCES_BLOCK, new UseableUnknownBlockParser());
		builder.setBlockParser(NexusBlockParser.UNKNOWN_BLOCK,	new UseableUnknownBlockParser());
		
		try {
			NexusFileFormat.parseFile(builder, file);
		} catch (ParseException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		NexusFile nexusFile = builder.getNexusFile();
		@SuppressWarnings("unchecked")
		Iterator<NexusBlock> blockIter = nexusFile.blockIterator();
		TreesBlock treesBlock = null;
		TreesBlock firstTreesBlock = null;
		int nTreesBlocks = 0;
		String regex = useLineageTrees ? "\\[Randomly selected lineage.*" : "\\[Randomly selected coalescent.*";
		while (blockIter.hasNext()) {
			NexusBlock block = blockIter.next();
			if (block.getBlockName().equalsIgnoreCase("trees")) {
				TreesBlock tempTreesBlock = (TreesBlock)block;
				nTreesBlocks++;
				if (nTreesBlocks == 1) firstTreesBlock = tempTreesBlock;
				@SuppressWarnings("rawtypes") // forced on me by Biojava.
				List comments = tempTreesBlock.getComments();
				if (comments.size()>0) {
					String comment = NexusUtils.toString((NexusComment)comments.get(0));
					if (comment.matches(regex)) {
						if (treesBlock!=null) throw new RuntimeException("Too many applicable trees blocks - can't happen");
						treesBlock = tempTreesBlock;
					} // if matches regex
				} // if comment.size>0
			} // if trees block
		} // block iteration
		if (treesBlock==null) {
			if (nTreesBlocks==1) {
				System.err.println("Only one trees block found, so using that one");
				treesBlock = firstTreesBlock;
			} else {
				throw new RuntimeException("Required trees block not found");
			}
		}
		return new Forest(treesBlock);
	}
}
