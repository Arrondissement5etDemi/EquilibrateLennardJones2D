public class Movement {
	private int ind;
	private double dx;
	private double dy;

	public Movement(int i, double x, double y) {
		ind = i;
		dx = x;
		dy = y;
	}	

	public int getInd() {
		return ind;
	}

	public double getDx() {
		return dx;
	}

	public double getDy() {
		return dy;
	}
	
	/** get the reverse movement where the direction is the opposite 
 * 	@return the reverse movement */
	public Movement reverse() {
		Movement result = new Movement(ind, -dx, -dy);
		return result;
	}
}
