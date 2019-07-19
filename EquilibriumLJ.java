import java.util.*;
import java.nio.file.*;
import java.io.*;
/**A class that generates equilibrium LJ configs in 2D*/

public class EquilibriumLJ {
	public static void main(String[] args) throws IOException{
		double rho = Double.parseDouble(args[0]);
                double t = Double.parseDouble(args[1]);
		int fileStart = Integer.parseInt(args[2]);
		int fileEnd = Integer.parseInt(args[3]);
		double rc = Double.parseDouble(args[4]);
		int rcInt = (int) rc;
		String rhoStr = args[0];
                String tStr = args[1];
		String dirName = "cutoff"+rcInt+"n800/rho"+rhoStr+"t"+tStr;
		if (!Files.exists(Paths.get(dirName))) { 
			Files.createDirectory(Paths.get(dirName));
		}
		Box equil = equilibrate(800,rho,t,rc);
		for (int i= fileStart; i<fileEnd; i++) {
			equil = equilibrate(equil,t);
			String p = dirName+"/file"+Integer.toString(i);
			/**String path = "rho"+rho100+"/HSfile"+Integer.toString(i);*/

			FileWriter f = new FileWriter(p);
			f.write(equil.toString());
			f.close();
		}
	}

	/**generates an equilibrium HS in 2D with given density
 *      @param n int, the number of particles in the fundamental simu box
 *      @param rho double, the number density
 *      @return Particle[] the list of the particles */
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

	/**generates an equilibrium HS in 2D from an existing config
 *  	@param pandora Box, the given already-equilibrated configuration
 *      @return Particle[] the list of the particles */
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


	/**generates an equilibrium LJ in 2D with given density and temperature 
 * 	@param n int, the number of particles in the fundamental simu box
 * 	@param rho double, the number density
 * 	@param temperature double, the scaled temperature
 * 	@param rc double, the cutoff of the potential 
 * 	@return Particle[] the list of the particles */
	public static Box equilibrate(int n, double rho, double temperature, double rc) {
		double d = Math.sqrt(n*(1.0/rho));	
		Box pandora = new Box(n,d,rc);
		double displace = 0.2;
		int sweeps = 10000;
		double tempEnergy = 0;
		for (double t = temperature; t >= temperature; t = t*0.98) {
			for (int i=1; i<= sweeps; i++) {
				for (int j=0; j < n; j++) {
                        		int ind = j;
                        		double oldE = pandora.partiE(ind);
                        		Movement proposedMove = pandora.move(ind,displace);
                        		double newE = pandora.partiE(ind);
					if (!Boltzmann.accept(newE,oldE,t)) {
                                		pandora.move(proposedMove.reverse());
					}
				}
				/**checking convergence
				if (i%500 == 0) {
					System.out.println(i+" "+tempEnergy/500);
					tempEnergy = 0;
				}
				else {
					tempEnergy = tempEnergy + pandora.getEnergy()/n;
				}
				*/
				/**checking speed*/
				if (i%500 == 0) {
                                        System.out.println(i);
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
                double displace = 0.2;
                int sweeps = 1000;
                for (double t = temperature; t >= temperature; t = t*0.98) {
                        for (int i=1; i<= sweeps; i++) {
				for (int j = 0; j < n; j++) {
                                	int ind = j;
                                	double oldE = pandora.partiE(ind);
                                	Movement proposedMove = pandora.move(ind,displace);
                                	double newE = pandora.partiE(ind);
                                	if (!Boltzmann.accept(newE,oldE,t)) {
                                        	pandora.move(proposedMove.reverse());
                                	}
				}
                        }
                }
                return pandora;
        }

}
