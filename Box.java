import java.util.*;


//this is a d*d square box with n Particles in the periodic boundary condition.
public class Box {
	private int n; //the number of particles
	private double d; //the dimension of the box
	private Particle[] partiArr;
	private double cutoff; //the cutoff distance for calculating the energy at a particle
	
	/** constructs a box with lattice-like initial configuration 
 * 	@param nGiven int, the given numer of particles
 * 	@param dGiven double, the given size of the box
 * 	@param co     double, the given cutoff distance
 * 	@return Box, the constructed box */
	public Box(int nGiven, double dGiven, double co) {
		n = nGiven;
		d = dGiven;
		cutoff = co;
		partiArr = new Particle[n];
		double delta = d/45.0;
		for (int i = 0; i < n; i++) {
			/**generate a candidate particle p*/
			int coordX = i%45;
			int coordY = i/45;
			double x = (double)coordX * delta;
			double y = (double)coordY * delta;
			Particle p = new Particle(x,y);
			partiArr[i] = p;	
		}
	}

	/** a box whose particle positions are predetermined
 *      @param nGiven int, the given numer of particles
 *      @param dGiven double, the given size of the box
 *      @param co     double, the given cutoff distance
 *      @return Box, the constructed box */
	public Box(int nGiven, double dGiven, double co, double[][] coordinates) {
		if (coordinates.length < nGiven) {
			throw IllegalArgumentException("the coordinate array size is smaller than nGiven.");	
		}
		n = nGiven;
                d = dGiven;
                cutoff = co;
		partiArr = new Particle[n];
		for (int i = 0; i < n; i++) {
			double x = coordinates[i][0];
			double y = coordinates[i][1];
			Particle p = new Particle(x,y);
			partiArr[i] = p;
		}
	}

	/** gets d, the dimension of the box
 * 	@return double */
	public double getD() {
		return d;
	}

	/** gets n, the number of particles
 *      @return int */
        public int getN() {
                return n;
        }

	/**computes the sum energy of all particles in the box
 * 	@return double, the ensemble energy */
	public double getEnergy() {
		double sum = 0;
		for (int i=0; i < n; i++) {
			sum = sum + partiE(i);
		}
		return sum/2.0; /**don't doublecount pair interactions!*/
	}

	/**gets the energy felt by a single particle in the box using minimum image
 * 	@param ind integer, the index of the particle in partiArr 
 * 	@return double, the energy felt by the particle
 * 	@assume Box side length > 2 * cutcoff */
	public double partiE(int ind) {
		if (ind >= n) {
			throw new IllegalArgumentException("input index must be < the number of particles");
		}
		if (2*cutoff >= d) {
			throw new RuntimeException("box dimension must > 2*cutoff for method partiE to work");	
		}
		Particle thisParti = partiArr[ind];
		/**For every particle in the box, compute their countribution to the particle*/
		double result = 0;
		for (int i = 0; i < n ; i++) {
			Particle thatParti = partiArr[i];
			double distToDupli = minDist(thisParti, thatParti); /**minimum image*/
			/**get the pair energy*/
			if (distToDupli <= cutoff && ind != i) {
				double pairEnergy = LjPotential.lj(distToDupli);
				result = result + pairEnergy;
			}
		
		}
		return result;
	}

	/**the minimum distance between two particles in PBC
 *     @param thisP Particle
 *     @param another Particle, another particle 
 *     @return the minimum distance to that particle among all its replicates in 2D */
	public double minDist(Particle thisP, Particle another) {
                double thatX = another.getx();
                double thatY = another.gety();
                double miniDx = minIn3(thatX,thatX+d,thatX-d,thisP.getx());
                double miniDy = minIn3(thatY,thatY+d,thatY-d,thisP.gety());
                return Math.sqrt(miniDx*miniDx + miniDy*miniDy);
        }

	public static double minIn3(double a, double b, double c, double center) {
		return Math.min(Math.abs(a-center),Math.min(Math.abs(b-center),Math.abs(c-center)));
	}                                                                                      }

	/**moves a random particle in a random direction by a random distance between 0 and maxDist
 * 	@param ind int, the index of the particle moved
 *	@param maxDist double, the maximum distance to move for a movement
 * 	@return Movement, what movement has occured. */
	public Movement move(int ind, double maxDist) {
		//randomly choose a direction
		double direction = getRandomNumberInRange(0.0,2*Math.PI);
		//randomly choose a distance to move
		double dist = getRandomNumberInRange(0,maxDist);
		Movement result = new Movement(ind, direction, dist);
		move(result);
		return result;
	}
	
	/**moves a particlar particle in a particular direction by a distance
 * 	@param m, Movement, including the index of the particle and the direction to move */
	public void move(Movement m) {
		int ind = m.getInd();
		double direction = m.getDirection();
		double moveDist = m. getMoveDist();
		if (direction >= 2*Math.PI || direction < 0) {
			throw new IllegalArgumentException("input direction must be between 0 and 2*PI");
		}
		if (ind >= n) {
                        throw new IllegalArgumentException("input index must be < the number of particles");
                }
		/**obtain the coordinates of the particle at ind*/
		double x = partiArr[ind].getx();
		double y = partiArr[ind].gety();
		/**obtain the amount of x and y to move*/
		double deltaX = Math.cos(direction)*moveDist;
		double deltaY = Math.sin(direction)*moveDist;
		/**get new coordinates of the partilce at ind*/
		partiArr[ind].setX(positiveModulo(x + deltaX,d));
		partiArr[ind].setY(positiveModulo(y + deltaY,d));	
	}

	/** returns the array of particles in the Box
 * 	@return a newly constructed array of particles that is a deep copy of partiArr */
	public Particle[] toArray() {
		Particle[] partiArrCopy = new Particle[n];
		for (int i = 0; i < n; i++) {
			double x = partiArr[i].getx();
			double y = partiArr[i].gety();
			Particle pCopy = new Particle(x,y);
			partiArrCopy[i] = pCopy;
		}
		return partiArrCopy;
	}

	/** returns the content of partiArray as a string
 * 	@return the content of partiArr */
	public String toString() {
		String result = "";
		for (int i = 0; i < n; i++) {
			result = result + partiArr[i].toString() + "\n";
		}
		return result;
	}

	/**auxillary functions*/
	public static double getRandomNumberInRange(double min, double max) {
		if (min >= max) {
			throw new IllegalArgumentException("max must be greater than min");
	        }		
		Random r = new Random();	
       		return r.nextDouble()*(max - min) + min;
	}
	
		
	private static double positiveModulo(double x, double d) {
		double result = x%d;
		if (result < 0) {
			result = result + d;
		}
		return result;
	}

	private static double expoRandomNumber(double lambda) {
		/**get a uniform random variable */
		Random r = new Random();
		double u = r.nextDouble();
		/**create the exponentially distributed random number*/
		double result = Math.log(1 - u)/(-lambda);
		return result;
	}
}
