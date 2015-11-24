package hybridstats;

import java.io.PrintWriter;

import palExtensions.Split;
/*
 * A pair of splits (on same taxon set) each with a weight.
 * Intended purpose is for Internode Certainty calculation (Solachos and Rokas, 2014)
 */

public class PairWeightedSplits {
	public Split split1;
	public Split split2;
	public double wt1;
	public double wt2;
	
	public PairWeightedSplits (Split s1, Split s2, double w1, double w2) {
		split1 = s1;
		split2 = s2;
		wt1 = w1;
		wt2 = w2;
	}

	public double internodeCertainty() {
		if (wt2==0) {
			return 1;
		} else {
			double p1 = wt1/(wt1+wt2);
			double p2 = wt2/(wt1+wt2);
			return 1+(p1*Math.log(p1)+p2*Math.log(p2))/Math.log(2);
		}
	}
	
	public void printInternodeCertainty(PrintWriter out) {
		if (split2 != null) {
			out.printf("Split %s (%.0f) vs %s (%.0f) = %f\n", split1.toString(), wt1, split2.toString(), wt2, internodeCertainty());
		} else {
			out.printf("Split %s (%.0f) vs nothing = 1.0\n", split1.toString(), wt1);
		}
	}
}
