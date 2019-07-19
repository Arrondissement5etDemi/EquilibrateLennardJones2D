import java.util.*;

public class Particle {
        //data
	private double x,y;

	//constructor
        public Particle(double xx, double yy) {
		x = xx;
                y = yy;
        }

        //gets the coordinates
        public double getx() {
		return x;
	}
    
	public double gety() {
		return y;
	}

	//modifiers
	public void setX(double xx) {
		x = xx;
	}

	public void setY(double yy) {
		y = yy;
	}

	//distance to another city
	public double distanceto(Particle another) {
		double anotherX = another.getx();
                double anotherY = another.gety();
                return Math.sqrt(Math.pow((x - anotherX),2) + Math.pow((y - anotherY),2));
	}

	/**the minimum distance between two particles in PBC
 *      @param another Particle, another particle 
 *      @param sideLength double, the side length of the unit cell
 *      @return the minimum distance to that particle among all its replicates */
	public double distanceto(Particle another, double sideLength) {
                double thatX = another.getx();
                double thatY = another.gety();
                double miniDx = minIn3(thatX,thatX+sideLength,thatX-sideLength,x);
                double miniDy = minIn3(thatY,thatY+sideLength,thatY-sideLength,y);
                return Math.sqrt(miniDx*miniDx + miniDy*miniDy);
        }

	public static double minIn3(double a, double b, double c, double center) {
		return Math.min(Math.abs(a-center),Math.min(Math.abs(b-center),Math.abs(c-center)));
	}


	public boolean equals(Particle another) {
		return (x == another.getx() && y == another.gety());
	}

	public String toString() {
		return Double.toString(x) + " " + Double.toString(y);
	}
}
