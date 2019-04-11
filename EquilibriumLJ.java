import java.util.*;
import java.io.*;
/**A class that generates equilibrium LJ configs in 2D*/

public class EquilibriumLJ {
	public static void main(String[] args) throws IOException{
		double rho = Double.parseDouble(args[0]);
                double t = Double.parseDouble(args[1]);
		Box equil = equilibrate(2000,rho,t); /**generates the initial equilibrated config*/
		for (int i= 0; i<1; i++) {
			equil = equilibrate(equil,t); /**generates more configs*/
			/**System.out.println(equil.getD());*/
			/**WRITE THE CONGIGS TO FILES, ONE FILE FOR EACH CONFIG!*/
			String rho100 = Integer.toString((int)Math.round(rho*100));
			String t100 = Integer.toString((int)Math.round(t*100));
			String path = "rho"+rho100+"t"+t100+"/file"+Integer.toString(i);
			/**String path = "rho"+rho100+"/HSfile"+Integer.toString(i); FOR HS TESTING*/
			FileWriter f = new FileWriter(path);
			f.write(equil.toString());
			f.close();
		}
	}

/**BENCHMARK TESTING OF 2d HARD SPHERE SYSTEMS
	generates an equilibrium HS in 2D with given density
 *      @param n int, the number of particles in the fundamental simu box
 *      @param rho double, the number density
 *      @return Particle[] the list of the particles 
        public static Box equilibrateHS(int n, double rho) {
                double d = Math.sqrt(n*(1.0/rho));
                Box pandora = new Box(n,d,2.5);
                double displace = 0.01*d;
                int sweeps = 1200;
                for (int i=1; i<= sweeps*n; i++) {
                        int ind = (int) Math.floor(Box.getRandomNumberInRange(0.0,n));
                        Movement proposedMove = pandora.move(ind,displace);
			boolean accept = true;
			Particle[] pa = pandora.toArray();
			Particle thisParti = pa[ind];
			for (int j = 0; j < n; j++) {
				Particle thatParti = pa[j];
				if (j!=ind && thisParti.distanceto(thatParti,d)<1) {
					accept = false;
					break;
				}
			}
                        if (!accept) {
                                pandora.move(proposedMove.reverse());
                        }
                }
                return pandora;
        }
	
	generates an equilibrium HS in 2D from an existing config
 *  	@param pandora Box, the given already-equilibrated configuration
 *      @return Particle[] the list of the particles 
        public static Box equilibrateHS(Box pandora) {
                int n = pandora.getN();
                double d = pandora.getD();
                double displace = 0.01*d;
                int sweeps = 400;
                for (int i=1; i<= sweeps*n; i++) {
			int ind = (int) Math.floor(Box.getRandomNumberInRange(0.0,n));
                        Movement proposedMove = pandora.move(ind,displace);
                        boolean accept = true;
                        Particle[] pa = pandora.toArray();
                        Particle thisParti = pa[ind];
                        for (int j = 0; j < n; j++) {
                                 Particle thatParti = pa[j];
                                 if (j!=ind && thisParti.distanceto(thatParti,d)<1) {
                                        accept = false;
                                        break;
                                }
                        }
                        if (!accept) {
                                pandora.move(proposedMove.reverse());
                        }
                }
                return pandora;
        }
*/

	/**generates an equilibrium LJ in 2D with given density and temperature 
 * 	@param n int, the number of particles in the fundamental simu box
 * 	@param rho double, the number density
 * 	@param temperature double, the scaled temperature 
 * 	@return Particle[] the list of the particles */
	public static Box equilibrate(int n, double rho, double temperature) {
		double d = Math.sqrt(n*(1.0/rho));	
		Box pandora = new Box(n,d,LjPotential.rc);/**create a box with n particles, dimension d, cutoff radius rc*/
		double displace = 0.6;
		int sweeps = 3000;
		for (double t = 1.0; t >= temperature; t = t*0.98) {
			for (int i=1; i<= sweeps*n; i++) {
                        	int ind = (int) Math.floor(Box.getRandomNumberInRange(0.0,n));
                        	double oldE = pandora.partiE(ind);
                        	Movement proposedMove = pandora.move(ind,displace);
                        	double newE = pandora.partiE(ind);
				if (!Boltzmann.accept(newE,oldE,t)) {
                                	pandora.move(proposedMove.reverse());
				}
			}
			System.out.println(t + " " + pandora.getEnergy()/n);
		}
		return pandora;		
	}

	/**generates an equilibrium LJ in 2D from an existing config at a given temperature
 * 	@param pandora Box, the given already-equilibrated configuration
 *      @param temperature double, the scaled temperature
 *      @return Particle[] the list of the particles */
        public static Box equilibrate(Box pandora, double temperature) {
		int n = pandora.getN();
                double d = pandora.getD();
                double displace = 0.6;
                int sweeps = 600;
                for (double t = temperature; t >= temperature; t = t*0.98) {
                        for (int i=1; i<= sweeps*n; i++) {
                                int ind = (int) Math.floor(Box.getRandomNumberInRange(0.0,n));
                                double oldE = pandora.partiE(ind);
                                Movement proposedMove = pandora.move(ind,displace);
                                double newE = pandora.partiE(ind);
                                if (!Boltzmann.accept(newE,oldE,t)) {
                                        pandora.move(proposedMove.reverse());
                                }
                        }
                        System.out.println(t + " " + pandora.getEnergy()/n);
                }
                return pandora;
        }

}
