package hybridstats;

public class CompoundCoefficient {
	private final double coefficient;
	private final String[] variables;
	
	public CompoundCoefficient(String[] vars, double coefficient) {
		variables = vars;
		this.coefficient = coefficient;
	}
	
	/**
	 * Parse compound coefficient from string, e.g. 
	 * 1.491e-03*QE*QE*TE
	 */
	public CompoundCoefficient(String input) {
		String[] bits = input.split("\\*");
		variables = new String[bits.length-1];
		try {
			coefficient = Double.valueOf(bits[0]);
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException("NumberFormatException trying to parse compound coefficient '"+input+"'");
		}
		for (int i=0; i<variables.length; i++) {
			variables[i]=bits[i+1].trim();
		}
	}
	
	public double evaluate(HybridStats stats) {
		double result = coefficient;
		for (String var : variables) {
			result *= stats.getStatByName(var);
		}
		return result;
	}

	public String toString() {
		StringBuffer buf = new StringBuffer(Double.toString(coefficient));
		for (String var : variables) {
			buf.append("*").append(var);
		}
		return buf.toString();
	}
}
