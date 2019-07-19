import java.util.*;

public class Boltzmann {
	//public static void main(String args[]) {
	//	System.out.println(boProb(-100.0,10.0));
	//}
	private static Random rnd = new Random();

	/**computes the Boltzmann probablility factor
 	@param deltaE double, the temperature difference
	@param temperature double
	@return double, the probability factor
	Assumption: deltaE >0 */
	private static double boProb(double deltaE, double temperature) {
		return Math.exp(- deltaE / temperature);
	}

	/**accepts or rejects the new situation according to the Boltzmann probablity
	@param newE double, the new energy
	@param oldE double, the old energy
	@param temperature double
	@return boolean, accept the new situation or not*/
	public static boolean accept(double newE, double oldE, double temperature) {
		//calculate the Boltzmann probability
		double deltaE = newE - oldE;
		if (deltaE <= 0) return true;
		double prob = boProb(deltaE, temperature);
		//get a random number between 0 and 1
		double x = rnd.nextDouble();
		//and accept or reject
		if (x <= prob) return true;
		else return false;
	}
}
