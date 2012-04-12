/**
 * Class to abstract away some of the probability computations
 * (so it can be done in log space or normal space)
 * @author arthiramachandran
 *
 */
public class Probability
{
	// all in log space
	// stores logProbability
	double value;
	
	public Probability(double value)
	{
		this.value = Math.log(value);
		if (Double.isNaN(this.value) || Double.isInfinite(this.value))
			this.value = Double.NEGATIVE_INFINITY;
	}
	
	public Probability multiply(Probability p)
	{
		Probability result = new Probability(this.value + p.value);
		return result;
	}
	
	public Probability add(Probability p)
	{
		double max = Math.max(this.value, p.value);
		double min = Math.min(this.value, p.value);
		if ( Double.isInfinite(min) || (max - min) >= 15.7)
			return new Probability(max);
		else
			return new Probability((max + Math.log(1 + Math.exp(-(max - min)))));
	}
	
	public static double add(double a, double d)
	{
		return (a+d);
	}

	public static Double multiply(double a, double b)
	{
		return (a*b);
	}

	// computes log (e^a + e^b)
	// a and b should be log values...
	public static double logSum(double a, double b)
	{
		double max = Math.max(a, b);
		double min = Math.min(a, b);
		if ( Double.isInfinite(min) || (max - min) >= 15.7)
			return max;
		else
			return (max + Math.log(1 + Math.exp(-(max - min))));
	}
}
