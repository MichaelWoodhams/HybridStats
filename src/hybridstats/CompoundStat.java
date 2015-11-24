package hybridstats;


import java.util.Vector;

import org.biojava.bio.seq.io.ParseException;

/*
 * A real-valued function of the summary statistics stored in a
 *  HybridStats object. (Functionality required for ABC project.)
 *  
 *  For now, it is a linear combination of HybridStats entries. In future,
 *  more complex functions may be allowed.
 */
public class CompoundStat {
	private String name;
	private Vector<CompoundCoefficient> coefficients;
	private double error; // std deviation of this stat
	
	/**
	 * Creates a CompoundStat with name and std deviation, but no coefficients.
	 * Use addCoefficient to supply them.
	 * @param name
	 * @param stdev
	 */
	public CompoundStat(String name, double stdev) {
		this.name = name;
		coefficients = new Vector<CompoundCoefficient>();
	}
	
	/**
	 * Parse a CompoundStat from inputLine
	 * @param inputLine
	 * 
	 * E.g. format of a compound stat:
	 * fitCoal(1.25) 8.138e-01 : TE*1.491e-03 : SI*2.012e-07 : DC*-4.405e-04 : TS*-7.387e-03
	 * (":" acts like "+", but would be harder to parse if we'd used "+".)
	 */
	//private static final Pattern NAME_STDDEV_PATTERN = Pattern.compile("(.*)\\((.+)\\)\\s*(.*)");
	public CompoundStat(String inputLine) throws ParseException {
		String[] firstPass = inputLine.split("[()]");
		//matcher matcher = NAME_STDDEV_PATTERN.matcher(inputLine);
		if (firstPass.length!=3) throw new ParseException("Could not parse compound stat specifier '"+inputLine+"'");
		name = firstPass[0];;
		error = Double.valueOf(firstPass[1]);
		String[] coefficientStrings = firstPass[2].split("\\s*:\\s*");
		coefficients = new Vector<CompoundCoefficient>();
		for (String coefString : coefficientStrings) {
			coefficients.add(new CompoundCoefficient(coefString));
		}
	}
	
	public String toString() {
		StringBuffer buf = new StringBuffer(name);
		buf.append("(").append(error).append(")");
		String separator = " ";
		for (CompoundCoefficient coef : coefficients) {
			buf.append(separator).append(coef.toString());
			separator = " : ";
		}
		return buf.toString();
	}
	
	public String getName() {
		return name;
	}
	
	public void addCoefficient(String coefString) {
		coefficients.add(new CompoundCoefficient(coefString));
	}
	
	public double evaluate(HybridStats stats) {
		double sum = 0;
		for (CompoundCoefficient coefficient : coefficients) {
			sum += coefficient.evaluate(stats);
		}
		return sum;
	}
	
	public void evaluateRange(HybridStats stats, Double[] range, double scale) {
		double eval = evaluate(stats);
		range[0] = eval-scale*error;
		range[1] = eval+scale*error;
	}
}
