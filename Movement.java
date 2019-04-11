public class Movement {
	private int ind;
	private double direction;
	private double moveDist;

	public Movement(int i, double d, double dist) {
		ind = i;
		direction = d;
		moveDist = dist;
	}	

	public int getInd() {
		return ind;
	}

	public double getDirection() {
		return direction;
	}

	public double getMoveDist() {
		return moveDist;
	}
	
	/** get the reverse movement where the direction is the opposite 
 * 	@return the reverse movement */
	public Movement reverse() {
		int resultI = ind;
		double resultDist = moveDist;
		double resultD = (Math.PI + direction)%(2*Math.PI);
		Movement result = new Movement(resultI, resultD, resultDist);
		return result;
	}
}
