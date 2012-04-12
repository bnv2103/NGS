
public class TestLogSum
{
	public static void main(String arg[])
	{
		double a = -15;
		double b = -5;
		
		double max = Math.max(a, b);
		double min = Math.min(a, b);
		if ( Double.isInfinite(min) || (max - min) >= 15.7)
			System.out.println(max);
		else
			System.out.println(max + Math.log(1 + Math.exp(-(max - min) )));
		System.out.println(Math.log(Math.exp(a) + Math.exp(b)));
	}
}
