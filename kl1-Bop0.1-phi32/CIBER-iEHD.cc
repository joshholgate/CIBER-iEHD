/*  CIBER-iEHD2.1.cc - a botched-together simulation of plasma-liquid interactions (set up for flow past drop)
    Trying to use memory more efficiently	 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

// grid and global definitions
#define IMAX 94			/* number of cells in r-direction */
#define JMAX 1050 		/* number of cells in z-direction */
#define INV_H 30.0		/* cell width */
#define TOL 1e-6		/* successive-over-relaxation (SOR) error tolerance */
#define RELAX 1.9		/* relaxation parameter for SOR */
#define TIME_STEP 0.1		/* time between writing data to file */
#define FINAL_TIME 3.0		/* termination time of simulation */
#define PI 3.141592653589793

#define SHEATH			/* turn this on for a sheath simulation (default is flow past objects) */
//#define CARTESIAN			/* turn this on for cartesian x-y coordiates (default is axisymmetric r-z) */
					/* FIXME: liquid region is currently for (r,z) only!!! */

// liquid region definitions
#define INV_RE 0.1			/* inverse Reynolds number */
#define INV_WE 10.0			/* inverse Weber number */
#define INV_WEP 1.0			/* inverse plasma Weber number */
#define INV_FRSQR 0.0			/* inverse square Froude number */
#define INIT_L 0.0			/* initial (normalized) angular momentum of drop */
#define DLDT 0.0			/* rate of increase/decrease of angular momentum */

// plasma region definitions
#define FLOW_VELOCITY 1.0		/* (magnitude of) velocity of cold ion flow (set as 1 for Bohm sheath) */
#define GRAIN_POTENTIAL -32.0		/* initial potential of sphere / surface */
#define GRAIN_Z 25.0			/* z-position of grain (bottom of simulation domain is z=0) OR height of sheath */
#define GRAIN_RAD 1.0			/* radius of spherical dust grain */

/*  Need a grid  */
typedef struct Grid {
	// liquid variables
	double u[IMAX+1][JMAX+2];	/* radial liquid velocity */
	double v[IMAX+2][JMAX+1];	/* axial liquid velocity */
	double p[IMAX+2][JMAX+2];	/* liquid pressure */
	double F[IMAX+2][JMAX+2];	/* intermediate-step radial liquid velocity */
	double G[IMAX+2][JMAX+2];	/* intermediate-step axial liquid velocity */
	double RHS[IMAX+2][JMAX+2];	/* intermediate values of RHS of pressure-Poisson equation */
	double vel_pot[IMAX+2][JMAX+2];	/* fictional velocity potential */
	double LS[IMAX+2][JMAX+2];	/* level-set function */
	double ANG_MOM;			/* (dimensionless) angular momentum of droplet */
	double MOMENT_OF_INERTIA;	/* (dimensionless) moment of inertia around z-axis */
	double LONGEST_AXIS;		/* elongation of droplet */
	double ANG_VEL;			/* (dimensionless) angular velocity of droplet */
	double INV_ROSQR;		/* gives information on current angular velocity */

	// plasma variables
	double dt;
	double E_pot[IMAX+2][JMAX+2];	/* electric potential */
	double ION_U[IMAX+1][JMAX+2];
	double ION_V[IMAX+2][JMAX+1];
	double ION_N[IMAX+2][JMAX+2];
	double GRAIN_PHI;
} Grid;


/*  Define a function to initialize the grid values  */
void init_liquid( Grid* grid ) {
	for (int i=0; i<=IMAX; i++) {
		for (int j=0; j<=JMAX+1; j++) {
			grid->u[i][j] = 0.0;
			grid->F[i][j] = 0.0;
		}
	}

	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX; j++) {
			grid->v[i][j] = 0.0;
			grid->G[i][j] = 0.0;
		}
	}

	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX+1; j++) {
			grid->p[i][j] = 0.0;
			grid->RHS[i][j] = 0.0;
			grid->E_pot[i][j] = 0.0;
		}
	}

#ifdef SHEATH
	// Set up LS for 'surface' problem
	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX+1; j++) {
			grid->LS[i][j] = -( (double(j)-0.5)/INV_H - 10.0 )
						 + 0.2*cos(1.0*(double(i)-0.5)/INV_H);
		}
	}
#else
	double r,z;
	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX+1; j++) {
			// Sphere
			r = ((double(i)-0.5)/INV_H);
			z = ((double(j)-0.5)/INV_H - GRAIN_Z);
			grid->LS[i][j] = GRAIN_RAD - sqrt(r*r+z*z);
		}
	}
#endif

	grid->ANG_MOM = INIT_L;
	grid->MOMENT_OF_INERTIA = 0.891094;		// Dimensionless moment of inertia of a sphere
	grid->LONGEST_AXIS = 1.0;
	grid->ANG_VEL = 0.0;
	grid->INV_ROSQR = 0.0;
}


/*  Define a function to initialize the grid values  */
void init_plasma( Grid* grid ) {
	grid->GRAIN_PHI = GRAIN_POTENTIAL;
	double R, Z, UN, HALF;
	UN = sqrt( FLOW_VELOCITY*FLOW_VELOCITY - 2.0*grid->GRAIN_PHI);
	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX+1; j++) {
#ifdef SHEATH
			Z = GRAIN_Z;
			grid->E_pot[i][j] = -6.0/(Z*Z);
			grid->ION_U[i][j] = 0.0;
			grid->ION_V[i][j] = -FLOW_VELOCITY - 6.0/(Z*Z);
#else
			grid->E_pot[i][j] = 0.0;
			/* Initialize velocity so that upwinding points away from grain - use solution for potential flow with 
			   sink (normal velocity -1) at radius 1 */
			R = (double(i)-0.5)/INV_H;
			Z = (double(j)-GRAIN_Z*INV_H-0.5)/INV_H;
			HALF = 0.5/INV_H;
			if ( i!=IMAX+1) {
				grid->ION_U[i][j] = -(R+HALF)*UN*GRAIN_RAD*GRAIN_RAD/((R+HALF)*(R+HALF)+Z*Z)
											/sqrt((R+HALF)*(R+HALF)+Z*Z)
						+ 1.5*GRAIN_RAD*GRAIN_RAD*GRAIN_RAD*Z*(R+HALF)*FLOW_VELOCITY
						/((R+HALF)*(R+HALF)+Z*Z)/((R+HALF)*(R+HALF)+Z*Z)/sqrt((R+HALF)*(R+HALF)+Z*Z);
			}
			if ( j!=JMAX+1) {
				grid->ION_V[i][j] = -(Z+HALF)*UN*GRAIN_RAD*GRAIN_RAD/(R*R+(Z+HALF)*(Z+HALF))
											/sqrt(R*R+(Z+HALF)*(Z+HALF))
						- (1.0 + 0.5*(R*R-2.0*(Z+HALF)*(Z+HALF))*GRAIN_RAD*GRAIN_RAD*GRAIN_RAD
							/(R*R+(Z+HALF)*(Z+HALF))/(R*R+(Z+HALF)*(Z+HALF))
										/sqrt(R*R+(Z+HALF)*(Z+HALF)))*FLOW_VELOCITY;
			}
#endif
			grid->F[i][j] = grid->ION_U[i][j];
			grid->G[i][j] = grid->ION_V[i][j];
		}
	}

	// Ion velocity/density zero inside object (so no density flux leaves it)
	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX+1; j++) {
			grid->ION_N[i][j] = 0.0;
#ifdef SHEATH
			Z = GRAIN_Z;
			if (grid->LS[i][j] < 0.0) grid->ION_N[i][j] = 1.0 - 6.0/(Z*Z);
#else
			if (grid->LS[i][j] < 0.0) grid->ION_N[i][j] = 1.0;	
#endif
		}
	}
}


/*  Define a function to write grid-centred values to file labelled with time t */
void write_grid( Grid* grid, double t) {
	double Ex, Ey;
	ostringstream fileNameStream("");
	fileNameStream << "output/t" << t << "_gridvalues.dat";
	string fileName = fileNameStream.str();
	ofstream myfile;
	myfile.open ( fileName.c_str() );

	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX+1; j++) {
			/* Reading out: r / z / u / v / p / level-set fn / electric potential / Er / Ez */
			/* Make sure no infs or nans are read out */
			Ex = (grid->E_pot[i-1][j]-grid->E_pot[i+1][j])*0.5*INV_H;
			Ex = (grid->E_pot[i][j-1]-grid->E_pot[i][j+1])*0.5*INV_H;
			myfile << (0.5+float(i-1))/INV_H << "\t" << (0.5+float(j-1))/INV_H << "\t" 
				<< (isinf(grid->u[i][j]) || isinf(grid->u[i][j]) ? 0.0 : grid->u[i][j]) << "\t"
				<< (isinf(grid->v[i][j]) || isinf(grid->v[i][j]) ? 0.0 : grid->v[i][j]) << "\t"
				<< (isinf(grid->p[i][j]) || isinf(grid->p[i][j]) ? 0.0 : grid->p[i][j]) << "\t"
				<< (isinf(grid->LS[i][j]) || isinf(grid->LS[i][j]) ? 0.0 : grid->LS[i][j]) << "\t"
				<< (isinf(grid->E_pot[i][j]) || isinf(grid->E_pot[i][j]) ? 0.0 : grid->E_pot[i][j]) << "\t"
				<< (isinf(Ex) || isinf(Ex) ? 0.0 : Ex) << "\t"
				<< (isinf(Ey) || isinf(Ey) ? 0.0 : Ey) << "\t" 
				<< (isinf(grid->ION_U[i][j]) || isinf(grid->ION_U[i][j]) ? 0.0 : grid->ION_U[i][j]) << "\t"
				<< (isinf(grid->ION_V[i][j]) || isinf(grid->ION_V[i][j]) ? 0.0 : grid->ION_V[i][j]) << "\t"
				<< (isinf(grid->ION_N[i][j]) || isinf(grid->ION_N[i][j]) ? 0.0 : grid->ION_N[i][j]) << endl;
		}
		myfile << endl;
	}
	myfile.close();
}


void write_timeseries( Grid* grid, double t) {
	/* also write fluid volume to file */
	double volume = 0.0;
	double r_j = 0.0;
	for (int j=1; j<=JMAX-1; j++) {
		if ( grid->LS[0][j] < 0.0 ) {
			r_j = 0.0;
		} else {
			for (int i=1; i<=IMAX; i++) {
				if (grid->LS[i+1][j] < 0.0) {
					r_j = (1.0/INV_H) * (double(i) - 0.5 +
						abs(grid->LS[i][j])/(abs(grid->LS[i][j])+abs(grid->LS[i+1][j])));
					break;
				}
				if (i == IMAX) r_j = 0.5;
			}
		}
 		volume += PI*r_j*r_j/INV_H;
	}

	ostringstream fileNameStream("");
	fileNameStream << "output/macroscopic.dat";
	string fileName = fileNameStream.str();
	ofstream myfile2;
	myfile2.open ( fileName.c_str(), ios_base::app);
#ifdef SHEATH
	/* amplitude of surface perturbation if doing surface instability */
	double ymin = 0.0, ymax = 0.0;
	for (int j=0; j<JMAX; j++) {
		if ( grid->LS[1][j]*grid->LS[1][j+1] < 0.0 ) ymax = double(j) - 0.5
						+ abs(grid->LS[1][j])/(abs(grid->LS[1][j])+abs(grid->LS[1][j+1]));
		if ( grid->LS[IMAX][j]*grid->LS[IMAX][j+1] < 0.0 ) ymin = double(j) - 0.5
				+ abs(grid->LS[IMAX][j])/(abs(grid->LS[IMAX][j])+abs(grid->LS[IMAX][j+1]));
	}
	double amplitude = ( ymax - ymin )/INV_H;

	/* and energy in the plasma region */
	double energy=0.0, maxenergy=0.0, xm=0.0, ym = 0.0, avenergy=0.0, npoints=0.0;
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			if (grid->LS[i][j] <  0.0) {
				energy = 0.5*(0.25*(grid->ION_V[i][j]+grid->ION_V[i][j-1])*(grid->ION_V[i][j]+grid->ION_V[i][j-1])
					+ 0.25*(grid->ION_U[i][j]+grid->ION_U[i-1][j])*(grid->ION_U[i][j]+grid->ION_U[i-1][j]) )
								+ grid->E_pot[i][j];
				energy = abs(energy-0.5);
				if (energy > maxenergy) {
					maxenergy = energy;
					xm = (double(i)-0.5)/INV_H;
					ym = (double(j)-0.5)/INV_H;
				}
				avenergy += energy;
				npoints += 1.0;
			}
		}
	}
	avenergy /= npoints;

	myfile2 << t << "\t" << amplitude << "\t" << volume << "\t"
				<< avenergy << "\t" << maxenergy << "\t" << xm << "\t" << ym << endl;
#else
	myfile2 << t << "\t" << grid->ANG_MOM << "\t" << grid->MOMENT_OF_INERTIA << "\t"
				<< grid->ANG_VEL << "\t" << grid->LONGEST_AXIS << "\t" << volume << endl;
#endif
	myfile2.close();
}


/*  This function sets DT to give numerically stable calculations
    Also have DT << Fr*sqrt(H) but this is almost always satisfied by CFL for realistic Fr and small H */
double set_dt( Grid* grid) {
	double max_velocity = 0.0;
	for (int i=1; i<=IMAX; i++) {			// loop such that all interior u and v covered
		for (int j=1; j<=JMAX; j++) {
			max_velocity = max( abs(grid->u[i][j]), max_velocity);	// make sure magnitude is used!
			max_velocity = max( abs(grid->v[i][j]), max_velocity);
		}
	}
	double dt = 0.5 * min ( 1.0/(INV_H*max_velocity), 0.25/(INV_RE*INV_H*INV_H) );
	dt = min (dt, 0.5/sqrt(INV_WE*INV_H*INV_H*INV_H*2.0*PI) );	// Need to resolve capillary waves (Brackbill 1992)
	dt = min (dt, 0.5/sqrt(INV_WEP*INV_H*INV_H*INV_H*2.0*PI) );	// And electrocapillary waves

	return dt;
}


/*  This function sets DT in sheath region to give numerically stable calculations according to CFL condition */
void set_plasma_dt( Grid* grid) {
	double max_velocity = 0.0;
	for (int i=1; i<=IMAX; i++) {			// loop such that all interior u and v covered
		for (int j=1; j<=JMAX; j++) {
			if (grid->LS[i][j] <  0.0) {
				max_velocity = max( abs(grid->ION_U[i][j]), max_velocity);
				max_velocity = max( abs(grid->ION_V[i][j]), max_velocity);
			}
		}
	}
	grid->dt = 0.25 /(INV_H*max_velocity);
}


/*  Function to solve Poisson-Boltzmann equation in plasma region */
void solve_PB_eqn( Grid* grid) {
	double temp, err, maxerr, R;			/* Use min norm for residual and tolerance */
	double S, E, W;					/* Switches to zero for points on the boundary */

#ifndef SHEATH
	for (int i=0; i<=IMAX+1; i++) {
		grid->E_pot[i][JMAX+1] = 0.0;
	}
#endif
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			if (grid->LS[i][j] >= 0.0) {
				grid->E_pot[i][j] = grid->GRAIN_PHI;
			}
		}
	}
	do {
		// Successive-over-relaxation.
		maxerr = 0.0;
		for (int i=1; i<=IMAX; i++) {
			for (int j=1; j<=JMAX; j++) {
				if (grid->LS[i][j] <  0.0) {	// only solve in vacuum
					E = 1.0;
					W = 1.0;
					S = 1.0;
					if (i == 1) W = 0.0;
					if (i == IMAX) E = 0.0;
					if (j == 1) S = 0.0;

#ifdef CARTESIAN
					temp = (1.0-RELAX)*grid->E_pot[i][j] + RELAX *
						( E*grid->E_pot[i+1][j] + W*grid->E_pot[i-1][j]
							+ grid->E_pot[i][j+1] + S*grid->E_pot[i][j-1] 
							+ (grid->ION_N[i][j] - exp(grid->E_pot[i][j]))/INV_H/INV_H )
							 	/(1.0+S+E+W) ;
#else
					R = double(i)-0.5;
					temp = (1.0-RELAX)*grid->E_pot[i][j] + RELAX *
						( (R+0.5)*E*grid->E_pot[i+1][j] + (R-0.5)*W*grid->E_pot[i-1][j]
							+ R*grid->E_pot[i][j+1] + R*S*grid->E_pot[i][j-1] 
							+ R*(grid->ION_N[i][j] - exp(grid->E_pot[i][j]))/INV_H/INV_H )
							 	/((1.0+S)*R+(R+0.5)*E+(R-0.5)*W) ;
#endif
					err = abs(grid->E_pot[i][j] - temp);
					if (err > maxerr) maxerr = err;
					grid->E_pot[i][j] = temp;
				}
			}
		}
	} while (maxerr > TOL);

	for (int j=1; j<=JMAX-1; j++) {
		grid->E_pot[0][j] = grid->E_pot[1][j];
		grid->E_pot[IMAX+1][j] = grid->E_pot[IMAX][j];
	}
#ifndef SHEATH
	for (int i=0; i<=IMAX+1; i++) {
		grid->E_pot[i][0] = grid->E_pot[i][1];
	}
#endif
}



/*  Function to calculate ion density with second order upwind differences */
bool calc_ion_n( Grid* grid) {

	double R, DNURDR, DNVDZ, RFLUXPLUS, RFLUXMINUS, ZFLUXPLUS, ZFLUXMINUS;
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			if (grid->LS[i][j] <  0.0) {	// only solve in vacuum
				RFLUXPLUS = max(grid->ION_U[i][j],0.0)*grid->ION_N[i][j]
								+ min(grid->ION_U[i][j],0.0)*grid->ION_N[i+1][j];
				RFLUXMINUS = max(grid->ION_U[i-1][j],0.0)*grid->ION_N[i-1][j]
								+ min(grid->ION_U[i-1][j],0.0)*grid->ION_N[i][j];

				ZFLUXPLUS = max(grid->ION_V[i][j],0.0)*grid->ION_N[i][j]
								+ min(grid->ION_V[i][j],0.0)*grid->ION_N[i][j+1];
				ZFLUXMINUS = max(grid->ION_V[i][j-1],0.0)*grid->ION_N[i][j-1]
								+ min(grid->ION_V[i][j-1],0.0)*grid->ION_N[i][j];

				DNVDZ = ZFLUXPLUS - ZFLUXMINUS;

#ifdef CARTESIAN
				DNURDR = RFLUXPLUS - RFLUXMINUS;
				grid->F[i][j] = grid->ION_N[i][j] - grid->dt*INV_H*(DNURDR + DNVDZ);
#else
				R = double(i)-0.5;
				DNURDR = RFLUXPLUS*(R+0.5) - RFLUXMINUS*(R-0.5);
				grid->F[i][j] = grid->ION_N[i][j] - grid->dt*INV_H*(DNURDR/R + DNVDZ);
#endif
			}
		}
	}

	/* FIXME: Is this section really necessary??? */
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			grid->ION_N[i][j] = min(grid->F[i][j], 1000.0);
			if (grid->LS[i][j] >= 0.0) grid->ION_N[i][j] = 0.0;
			if ( grid->ION_N[i][j] < 0.0 ) grid->ION_N[i][j] = 0.25*(grid->ION_N[i+1][j]+grid->ION_N[i-1][j]
											+grid->ION_N[i][j+1]+grid->ION_N[i][j-1]);
			if ( isnan(grid->ION_N[i][j]) == true ) {
				cout << "N is nan at (" << (double(i)-0.5)/INV_H << ", " << (double(j)-0.5)/INV_H << ")" << endl; 
				return true;
			}
		}
	}

	/* tidy up boundaries (fixed inflow, quasineutral Boltzmann approximation at outer/bottom edges)*/
#ifndef SHEATH
	for (int i=1; i<=IMAX; i++) {
		grid->ION_N[i][0] = exp(grid->E_pot[i][0]);
		grid->ION_N[i][JMAX+1] = 1.0;
	}
#endif
	for (int j=0; j<=JMAX+1; j++) {
		grid->ION_N[0][j] = grid->ION_N[1][j];
#if defined(CARTESIAN) || defined(SHEATH)
		grid->ION_N[IMAX+1][j] = grid->ION_N[IMAX][j];
#else
		grid->ION_N[IMAX+1][j] = exp(grid->E_pot[IMAX+1][j]);
#endif
	}

	return false;
}


/*  Function to calculate ion velocity (convective derivative calculation taken from iEHD) */
bool calc_ion_uv( Grid* grid) {
	double UDUDR, VDUDZ, UDVDR, VDVDZ, DPHIDR, DPHIDZ, DUMINUS, DUPLUS, DVMINUS, DVPLUS;
	double max_vel = sqrt(FLOW_VELOCITY*FLOW_VELOCITY-2.0*grid->GRAIN_PHI);
	for (int i=1; i<=IMAX-1; i++) {				/* Looping over U points only */
		for (int j=1; j<=JMAX; j++) {
			if ( grid->LS[i][j] <  0.0 || grid->LS[i+1][j] < 0.0 ) {	// only solve in vacuum

				DUMINUS = (3.0*grid->ION_U[i][j] - 4.0*grid->ION_U[i-1][j] + grid->ION_U[i-2][j])*0.5*INV_H;
				if ( i==1 ) DUMINUS = (grid->ION_U[1][j]-grid->ION_U[0][j])*INV_H;
				DUPLUS = (-grid->ION_U[i+2][j] + 4.0*grid->ION_U[i+1][j] - 3.0*grid->ION_U[i][j])*0.5*INV_H;
				if ( i==IMAX-1 ) DUPLUS = (grid->ION_U[IMAX][j]-grid->ION_U[IMAX-1][j])*INV_H;
				UDUDR = max(grid->ION_U[i][j],0.0)*DUMINUS + min(grid->ION_U[i][j],0.0)*DUPLUS;

				DUMINUS = (3.0*grid->ION_U[i][j] - 4.0*grid->ION_U[i][j-1] + grid->ION_U[i][j-2])*0.5*INV_H;
				if ( j==1 ) DUMINUS = (grid->ION_U[i][1]-grid->ION_U[i][0])*INV_H;
				DUPLUS = (-grid->ION_U[i][j+2] + 4.0*grid->ION_U[i][j+1] - 3.0*grid->ION_U[i][j])*0.5*INV_H;
				if ( j==JMAX ) DUPLUS = (grid->ION_U[i][JMAX+1]-grid->ION_U[i][JMAX])*INV_H;
				VDUDZ = max((grid->ION_V[i][j]+grid->ION_V[i+1][j]+grid->ION_V[i][j-1]
								+grid->ION_V[i+1][j-1])*0.25,0.0) * DUMINUS
					+ min((grid->ION_V[i][j]+grid->ION_V[i+1][j]+grid->ION_V[i][j-1]
								+grid->ION_V[i+1][j-1])*0.25,0.0) * DUPLUS;

				DPHIDR = ( grid->E_pot[i+1][j] - grid->E_pot[i][j] ) * INV_H;

				grid->F[i][j] = grid->ION_U[i][j] - grid->dt*( UDUDR + VDUDZ + DPHIDR );
			}
		}
	}

	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX-1; j++) {
			if ( grid->LS[i][j] <  0.0 || grid->LS[i][j+1] < 0.0 ) {	// only solve in vacuum

				DVMINUS = (3.0*grid->ION_V[i][j] - 4.0*grid->ION_V[i][j-1] + grid->ION_V[i][j-2])*0.5*INV_H;
				if ( j==1 ) DVMINUS = (grid->ION_V[i][1]-grid->ION_V[i][0])*INV_H;
				DVPLUS = (-grid->ION_V[i][j+2] + 4.0*grid->ION_V[i][j+1] - 3.0*grid->ION_V[i][j])*0.5*INV_H;
				if ( j==JMAX-1 ) DVPLUS = (grid->ION_V[i][JMAX]-grid->ION_V[i][JMAX-1])*INV_H;
				VDVDZ = max(grid->ION_V[i][j],0.0)*DVMINUS + min(grid->ION_V[i][j],0.0)*DVPLUS;

				DVMINUS = (3.0*grid->ION_V[i][j] - 4.0*grid->ION_V[i-1][j] + grid->ION_V[i-2][j])*0.5*INV_H;
				if ( i==1 ) DVMINUS = (grid->ION_V[1][j]-grid->ION_V[0][j])*INV_H;
				DVPLUS = (-grid->ION_V[i+2][j] + 4.0*grid->ION_V[i+1][j] - 3.0*grid->ION_V[i][j])*0.5*INV_H;
				if ( i==IMAX ) DVPLUS = (grid->ION_V[IMAX+1][j]-grid->ION_V[IMAX][j])*INV_H;
				UDVDR = max((grid->ION_U[i][j]+grid->ION_U[i][j+1]+grid->ION_U[i-1][j+1]
								+grid->ION_U[i-1][j])*0.25,0.0) * DVMINUS
					+ min((grid->ION_U[i][j]+grid->ION_U[i][j+1]+grid->ION_U[i-1][j+1]
								+grid->ION_U[i-1][j])*0.25,0.0) * DVPLUS;
	
				DPHIDZ = ( grid->E_pot[i][j+1] - grid->E_pot[i][j] ) * INV_H;

				grid->G[i][j] = grid->ION_V[i][j] - grid->dt*( VDVDZ + UDVDR + DPHIDZ );
			}
		}
	}

	// FIXME: is this bit really necessary?
	for (int i=1; i<=IMAX-1; i++) {				/* Looping over U points only */
		for (int j=1; j<=JMAX; j++) {
			if ( grid->LS[i][j] <  0.0 || grid->LS[i+1][j] < 0.0 ) {
				grid->ION_U[i][j] = grid->F[i][j];
				if ( isnan(grid->ION_U[i][j]) == true ) {
					cout << "U is nan at (" << (double(i)-0.5)/INV_H << ", "
									<< (double(j)-0.5)/INV_H << ")" << endl; 
					return true;
				}
			}
		}
	}

	// FIXME: or this?
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX-1; j++) {
			if ( grid->LS[i][j] <  0.0 || grid->LS[i][j+1] < 0.0 ) {
				grid->ION_V[i][j] = grid->G[i][j];
				if ( isnan(grid->ION_V[i][j]) == true ) {
					cout << "V is nan at (" << (double(i)-0.5)/INV_H << ", "
									<< (double(j)-0.5)/INV_H << ")" << endl; 
					return true;
				}
			}
		}
	}

	/* top/bottom boundary conditions */
#ifndef SHEATH
	for (int i=0; i<=IMAX; i++) {
		grid->ION_U[i][JMAX+1] = 0.0;
		grid->ION_U[i][0] = grid->ION_U[i][1];
		grid->ION_V[i][0] = -FLOW_VELOCITY + grid->E_pot[i][1]/FLOW_VELOCITY;
		grid->ION_V[i][JMAX] = -FLOW_VELOCITY;
	}
	grid->ION_V[IMAX+1][0] = -FLOW_VELOCITY + grid->E_pot[IMAX+1][1]/FLOW_VELOCITY;
#endif
#if defined(CARTESIAN) || defined(SHEATH)
	for (int j=1; j<=JMAX; j++) {
		grid->ION_U[0][j] = 0.0;
		grid->ION_U[IMAX][j] = 0.0;
		grid->ION_V[0][j] = grid->ION_V[1][j];
		grid->ION_V[IMAX+1][j] = grid->ION_V[IMAX][j];
	}
#else
	double V;
	for (int j=1; j<=JMAX; j++) {
		grid->ION_U[0][j] = 0.0;
		V = 0.5*( grid->ION_V[IMAX][j-1] + grid->ION_V[IMAX][j] );
		grid->ION_U[IMAX][j] = (double(IMAX)-1.0)*grid->ION_U[IMAX-1][j]/double(IMAX)
					+ (double(IMAX)-0.5)*(1.0/V-V )
							*0.5*(grid->E_pot[IMAX][j+1] - grid->E_pot[IMAX][j-1])/double(IMAX);
		grid->ION_V[0][j] = grid->ION_V[1][j];
		grid->ION_V[IMAX+1][j] = -FLOW_VELOCITY + 0.5*(grid->E_pot[IMAX+1][j] + grid->E_pot[IMAX+1][j+1])/FLOW_VELOCITY;
	}
#endif
	
	return false;
}



/*  Need a function to progress the the grid by one timestep:
	1. Update F and G and calculate RHS of pressure Poisson equation
	2. Solve pressure Poisson equation with SOR
	3. Update u and v						*/

bool step_grid( Grid* grid, double DT) {
	grid->ANG_MOM += DLDT*DT;

	/* Find longest axis by assuming it is along z axis. Also find normalized moment of inertia. */
	grid->MOMENT_OF_INERTIA = 0.0;
	double maj_rad_upper = 0.0, maj_rad_lower = 0.0, r_j = 0.0;
	for (int j=1; j<=JMAX-1; j++) {

		if ( (grid->LS[1][j]*grid->LS[1][j+1] < 0.0) && (j > JMAX/2) ) {
			 maj_rad_upper = double(j) + abs(grid->LS[1][j])/(abs(grid->LS[1][j])+abs(grid->LS[1][j+1]));
		}
		if ( (grid->LS[1][j]*grid->LS[1][j+1] < 0.0) && (j < JMAX/2) ) {
			 maj_rad_lower = double(j) + abs(grid->LS[1][j])/(abs(grid->LS[1][j])+abs(grid->LS[1][j+1]));
		}

		if ( grid->LS[0][j] < 0.0 ) {
			r_j = 0.0;
		} else {
			for (int i=1; i<=IMAX; i++) {
				if (grid->LS[i+1][j] < 0.0) {
					r_j = (double(i) - 0.5 +
						abs(grid->LS[i][j])/(abs(grid->LS[i][j])+abs(grid->LS[i+1][j])));
					break;
				}
			}
		}
 		grid->MOMENT_OF_INERTIA += 0.5*(double(j)-double(JMAX)*0.5)*(double(j)-double(JMAX)*0.5)*r_j*r_j 
							+ 0.125*r_j*r_j*r_j*r_j;
	}

	/* Here follows a horrendous conversion between different dimensionless conventions - my deepest apologies for this */ 
	grid->LONGEST_AXIS = 0.5*(maj_rad_upper-maj_rad_lower)/INV_H;
	grid->MOMENT_OF_INERTIA *= 3.341602735 / (INV_H*INV_H*INV_H*INV_H*INV_H);
	grid->ANG_VEL = grid->ANG_MOM / grid->MOMENT_OF_INERTIA;
	grid->INV_ROSQR = 8.0 * grid->ANG_VEL * grid->ANG_VEL * INV_WE;

	double DUPLUS, DUMINUS, UDUDR, VDUDZ, R, LAPLU;		/* Variables for derivatives and Laplacian */
	for (int i=1; i<=IMAX-1; i++) {				/* Looping over interior points only */
		for (int j=1; j<=JMAX; j++) {

			DUMINUS = (3.0*grid->u[i][j] - 4.0*grid->u[i-1][j] + grid->u[i-2][j])*0.5*INV_H;
			if ( i==1 ) DUMINUS = (grid->u[1][j]-grid->u[0][j])*INV_H;
			DUPLUS = (-grid->u[i+2][j] + 4.0*grid->u[i+1][j] - 3.0*grid->u[i][j])*0.5*INV_H;
			if ( i==IMAX-1 ) DUPLUS = (grid->u[IMAX][j]-grid->u[IMAX-1][j])*INV_H;
			UDUDR = max(grid->u[i][j],0.0)*DUMINUS + min(grid->u[i][j],0.0)*DUPLUS;

			DUMINUS = (3.0*grid->u[i][j] - 4.0*grid->u[i][j-1] + grid->u[i][j-2])*0.5*INV_H;
			if ( j==1 ) DUMINUS = (grid->u[i][1]-grid->u[i][0])*INV_H;
			DUPLUS = (-grid->u[i][j+2] + 4.0*grid->u[i][j+1] - 3.0*grid->u[i][j])*0.5*INV_H;
			if ( j==JMAX ) DUPLUS = (grid->u[i][JMAX+1]-grid->u[i][JMAX])*INV_H;
			VDUDZ = max((grid->v[i][j]+grid->v[i+1][j]+grid->v[i][j-1]+grid->v[i+1][j-1])*0.25,0.0)
											*DUMINUS
				+ min((grid->v[i][j]+grid->v[i+1][j]+grid->v[i][j-1]+grid->v[i+1][j-1])*0.25,0.0)
											*DUPLUS;

			R = double(i);
			LAPLU = ( (R+0.5)*grid->u[i+1][j]/R + (R-0.5)*grid->u[i-1][j]/R + grid->u[i][j+1]
					+ grid->u[i][j-1] - 4.0*grid->u[i][j] )*INV_H*INV_H - grid->u[i][j]/(R*R);

			grid->F[i][j] = grid->u[i][j] + DT*(INV_RE*LAPLU - UDUDR - VDUDZ);
		}
	}

	/* and for the boundary values (only need DU/DR) */
	for (int j=0; j<=JMAX; j++) {
		grid->F[0][j] = grid->u[0][j];
		grid->F[IMAX][j] = grid->u[IMAX][j];
	}

	double DVPLUS, DVMINUS, UDVDR, VDVDZ, LAPLV;	/* Variables for derivatives and Laplacian */
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX-1; j++) {

			DVMINUS = (3.0*grid->v[i][j] - 4.0*grid->v[i][j-1] + grid->v[i][j-2])*0.5*INV_H;
			if ( j==1 ) DVMINUS = (grid->v[i][1]-grid->v[i][0])*INV_H;
			DVPLUS = (-grid->v[i][j+2] + 4.0*grid->v[i][j+1] - 3.0*grid->v[i][j])*0.5*INV_H;
			if ( j==JMAX-1 ) DVPLUS = (grid->v[i][JMAX]-grid->v[i][JMAX-1])*INV_H;
			VDVDZ = max(grid->v[i][j],0.0)*DVMINUS + min(grid->v[i][j],0.0)*DVPLUS;


			DVMINUS = (3.0*grid->v[i][j] - 4.0*grid->v[i-1][j] + grid->v[i-2][j])*0.5*INV_H;
			if ( i==1 ) DVMINUS = (grid->v[1][j]-grid->v[0][j])*INV_H;
			DVPLUS = (-grid->v[i+2][j] + 4.0*grid->v[i+1][j] - 3.0*grid->v[i][j])*0.5*INV_H;
			if ( i==IMAX ) DVPLUS = (grid->v[IMAX+1][j]-grid->v[IMAX][j])*INV_H;
			UDVDR = max((grid->u[i][j]+grid->u[i][j+1]+grid->u[i-1][j+1]+grid->u[i-1][j])*0.25,0.0)
											*DVMINUS
				+ min((grid->u[i][j]+grid->u[i][j+1]+grid->u[i-1][j+1]+grid->u[i-1][j])*0.25,0.0)
											*DVPLUS;

			R = double(i)-0.5;
			LAPLV = ( (R+0.5)*grid->v[i+1][j]/R + (R-0.5)*grid->v[i-1][j]/R + grid->v[i][j+1]
					+ grid->v[i][j-1] - 4.0*grid->v[i][j] )*INV_H*INV_H;

			grid->G[i][j] = grid->v[i][j] + DT*(INV_RE*LAPLV - VDVDZ-UDVDR - INV_FRSQR
						+ (double(j)-1.0-double(JMAX)*0.5)*grid->INV_ROSQR/INV_H);
		}
	}

	/* and for the boundary values (only need DV/DZ) */
	for (int i=0; i<=IMAX; i++) {
		grid->G[i][0] = grid->v[i][0];
		grid->G[i][JMAX] = grid->v[i][JMAX];
	}

	/* Calculate RHS of pressure Poisson equation from updated F and G */
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			R = double(i)-0.5;
			grid->RHS[i][j] = ( (R+0.5)*grid->F[i][j]/R - (R-0.5)*grid->F[i-1][j]/R + grid->G[i][j]
							- grid->G[i][j-1] ) * INV_H/DT;
		}
	}

	/* And now initialize exterior pressure values as those given by Young-Laplace and Maxwell stress */
	double DXM, DXP, DLSDR, DLSDZ, DLSDRR, DLSDZZ, DLSDRZ, NR, NZ, CURV, ER, EZ, UR, UZ, RAM;
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			if (grid->LS[i][j] < 0.0) {
				if ( grid->LS[i+1][j] >= 0.0 || grid->LS[i-1][j] >= 0.0
						|| grid->LS[i][j+1] >= 0.0 || grid->LS[i][j-1] >= 0.0) {

					DXM = (3.0*grid->LS[i][j] - 4.0*grid->LS[i-1][j]
									 + grid->LS[i-2][j])*0.5*INV_H;
					if ( i==1 ) DXM = (grid->LS[1][j]-grid->LS[0][j])*INV_H;
					DXP = (-grid->LS[i+2][j] + 4.0*grid->LS[i+1][j]
									 - 3.0*grid->LS[i][j])*0.5*INV_H;
					if ( i==IMAX ) DXP = (grid->LS[IMAX+1][j]-grid->LS[IMAX][j])*INV_H;
					DLSDR = min( DXM, 0.0) + max( DXP, 0.0);

					DXM = (3.0*grid->LS[i][j] - 4.0*grid->LS[i][j-1]
									 + grid->LS[i][j-2])*0.5*INV_H;
					if ( j==1 ) DXM = (grid->LS[i][1]-grid->LS[i][0])*INV_H;
					DXP = (-grid->LS[i][j+2] + 4.0*grid->LS[i][j+1]
									 - 3.0*grid->LS[i][j])*0.5*INV_H;
					if ( j==JMAX ) DXP = (grid->LS[i][JMAX+1]-grid->LS[i][JMAX])*INV_H;
					DLSDZ = min( DXM, 0.0) + max( DXP, 0.0);

					/* I've had a few issues with DLSDR = DLSDZ = 0... The following code
										is a quick and dirty fix */
					if (DLSDR == 0.0) {
						if ( i==1 ) {
							DLSDR = (-2.0*grid->LS[i-1][j]-3.0*grid->LS[i][j]
								+6.0*grid->LS[i+1][j]-grid->LS[i+2][j])*INV_H/6.0;
						} else if ( i==IMAX ) {
							DLSDR = (grid->LS[i-2][j]-6.0*grid->LS[i-1][j]
							      +3.0*grid->LS[i][j]+2.0*grid->LS[i+1][j])*INV_H/6.0;
						} else {
							DLSDR = (grid->LS[i-2][j]-8.0*grid->LS[i-1][j]
							       +8.0*grid->LS[i+1][j]-grid->LS[i+2][j])*INV_H/12.0;
						}
						if ( j==1 ) {
							DLSDZ = (-2.0*grid->LS[i][j-1]-3.0*grid->LS[i][j]
								+6.0*grid->LS[i][j+1]-grid->LS[i][j+2])*INV_H/6.0;
						} else if ( j==JMAX ) {
							DLSDZ = (grid->LS[i][j-2]-6.0*grid->LS[i][j-1]
							      +3.0*grid->LS[i][j]+2.0*grid->LS[i][j+1])*INV_H/6.0;
						} else {
							DLSDZ = (grid->LS[i][j-2]-8.0*grid->LS[i][j-1]
							       +8.0*grid->LS[i][j+1]-grid->LS[i][j+2])*INV_H/12.0;
						}
					}

					DLSDRR = (grid->LS[i+1][j]-2.0*grid->LS[i][j]
									+grid->LS[i-1][j])*INV_H*INV_H;
					DLSDZZ = (grid->LS[i][j+1]-2.0*grid->LS[i][j]
									+grid->LS[i][j-1])*INV_H*INV_H;
					DLSDRZ = (grid->LS[i+1][j+1]-grid->LS[i+1][j-1]-grid->LS[i-1][j+1]
									+grid->LS[i-1][j-1])*0.25*INV_H*INV_H;
					R = double(i)-0.5;
					NR =  - DLSDR / sqrt(DLSDR*DLSDR+DLSDZ*DLSDZ);
					NZ =  - DLSDZ / sqrt(DLSDR*DLSDR+DLSDZ*DLSDZ);
					CURV = - NR*INV_H/R    // R should be in length units, not multiples of i
						+ (DLSDZ*DLSDZ*DLSDRR+DLSDR*DLSDR*DLSDZZ-2.0*DLSDR*DLSDZ*DLSDRZ)
						/( sqrt(DLSDR*DLSDR+DLSDZ*DLSDZ)*(DLSDR*DLSDR+DLSDZ*DLSDZ) );

					if ( DLSDR < 0.0 ) {
						ER = (-grid->E_pot[i+2][j] + 4.0*grid->E_pot[i+1][j]
									 - 3.0*grid->E_pot[i][j])*0.5*INV_H;
						if ( i==IMAX ) ER = (grid->E_pot[IMAX+1][j]
									 - grid->E_pot[IMAX][j])*INV_H;
					} else {
						ER = (3.0*grid->E_pot[i][j] - 4.0*grid->E_pot[i-1][j]
									 + grid->E_pot[i-2][j])*0.5*INV_H;
						if ( i==1 ) ER = (grid->E_pot[1][j]-grid->E_pot[0][j])*INV_H;
					}

					if ( DLSDZ < 0.0 ) {
						EZ = (-grid->E_pot[i][j+2] + 4.0*grid->E_pot[i][j+1]
									 - 3.0*grid->E_pot[i][j])*0.5*INV_H;
						if ( j==JMAX ) EZ = (grid->E_pot[i][JMAX+1]
									 - grid->E_pot[i][JMAX])*INV_H;
					} else {
						EZ = (3.0*grid->E_pot[i][j] - 4.0*grid->E_pot[i][j-1]
									 + grid->E_pot[i][j-2])*0.5*INV_H;
						if ( j==1 ) EZ = (grid->E_pot[i][1]-grid->E_pot[i][0])*INV_H;
					}

					UR = 0.5*(grid->ION_U[i][j] + grid->ION_U[i-1][j]);
					UZ = 0.5*(grid->ION_U[i][j] + grid->ION_U[i][j-1]);
					RAM = grid->ION_N[i][j]*(NR*UR + NZ*UZ)*(NR*UR + NZ*UZ);

					grid->p[i][j] = - INV_WE*CURV
							- ( (ER*ER+EZ*EZ)*0.5 - RAM - exp(grid->E_pot[i][j]) )*INV_WEP;

					if ( isnan(grid->p[i][j]) == true ) {
						cout << "P_out at boundary cell (i,j) = (" << i << "," << j <<
							") is nan!" << endl;
						cout << "Divided by DLSDR*DLSDR+DLSDZ*DLSDZ = " <<
							DLSDR*DLSDR+DLSDZ*DLSDZ << " at this point..." << endl;
						cout << "DLSDR = " << DLSDR << ", DLSDZ = " << DLSDZ << endl; 
					}
				} else grid->p[i][j] = 0.0;
			}
		}
	}

	// Introduce variables for pressure Poisson equation
	double temp, err, maxerr;			/* Use min norm for residual and tolerance */
	double N, S, E, W;				/* Switches to zero for points on the boundary */
	do {
		/* Griebel recommends "copying the pressure values along the boundary to their neighboring cells 			   in the boundary strip prior to each iteration step." So that's what I do here.		*/
		for (int i=1; i<=IMAX; i++) {
			grid->p[i][0] = grid->p[i][1];
			grid->p[i][JMAX+1] = grid->p[i][JMAX];
		}
		for (int j=1; j<=JMAX+1; j++) {
			grid->p[0][j] = grid->p[1][j];
			grid->p[IMAX+1][j] = grid->p[IMAX][j];
		}

		/* Then SOR */
		maxerr = 0.0;
		for (int i=1; i<=IMAX; i++) {
			for (int j=1+(i&1); j<=JMAX; j+=2) {
				if (grid->LS[i][j] >=  0.0) {	/* only solve in fluid! */
					N = 1.0;
					S = 1.0;
					E = 1.0;
					W = 1.0;
					if (i == 1) W = 0.0;
					if (i == IMAX) E = 0.0;
					if (j == 1) S = 0.0;
					if (j == JMAX) N = 0.0;
					R = double(i)-0.5;
					temp = (1.0-RELAX)*grid->p[i][j] + RELAX*( (R+0.5)*E*grid->p[i+1][j]
							+ (R-0.5)*W*grid->p[i-1][j] + R*N*grid->p[i][j+1]
							+ R*S*grid->p[i][j-1] - R*grid->RHS[i][j]/(INV_H*INV_H))
							 	/(R*N+R*S+(R+0.5)*E+(R-0.5)*W);
					err = fabs(grid->p[i][j] - temp);
					if (err > maxerr) maxerr = err;
					grid->p[i][j] = temp;
					if ( isnan(temp) == true ) {
						cout << "Error in pressure Poisson solver at (i,j) = (" <<
							i << "," << j << ")" << endl;
						return false;
					}
				}
			}
		}
		maxerr = 0.0;
		for (int i=1; i<=IMAX; i++) {
			for (int j=2-(i&1); j<=JMAX; j+=2) {
				if (grid->LS[i][j] >=  0.0) {	/* only solve in fluid! */
					N = 1.0;
					S = 1.0;
					E = 1.0;
					W = 1.0;
					if (i == 1) W = 0.0;
					if (i == IMAX) E = 0.0;
					if (j == 1) S = 0.0;
					if (j == JMAX) N = 0.0;
					R = double(i)-0.5;
					temp = (1.0-RELAX)*grid->p[i][j] + RELAX*( (R+0.5)*E*grid->p[i+1][j]
							+ (R-0.5)*W*grid->p[i-1][j] + R*N*grid->p[i][j+1]
							+ R*S*grid->p[i][j-1] - R*grid->RHS[i][j]/(INV_H*INV_H))
							 	/(R*N+R*S+(R+0.5)*E+(R-0.5)*W);
					err = fabs(grid->p[i][j] - temp);
					if (err > maxerr) maxerr = err;
					grid->p[i][j] = temp;
					if ( isnan(temp) == true ) {
						cout << "Error in pressure Poisson solver at (i,j) = (" <<
							i << "," << j << ")" << endl;
						return false;
					}
				}
			}
		}
	} while (maxerr > TOL);

	/* Update velocities, and we're done! */
	for (int i=1; i<=IMAX-1; i++) {
		for (int j=1; j<=JMAX; j++) {
			if ( grid->LS[i][j] >= 0.0 || grid->LS[i+1][j] >= 0.0 ) {	/* only solve in fluid! */
				grid->u[i][j] = grid->F[i][j]-(grid->p[i+1][j]-grid->p[i][j])*DT*INV_H;
			}
		}
	}
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX-1; j++) {
			if (grid->LS[i][j] >= 0.0 || grid->LS[i][j+1] >= 0.0 ) {	/* only solve in fluid! */
				grid->v[i][j] = grid->G[i][j]-(grid->p[i][j+1]-grid->p[i][j])*DT*INV_H;
			}
		}
	}

	return true;
}


void velocity_extension_and_BCs( Grid* grid) {		/* solving div(u) everywhere in vacuum to preserve volume
							 in advection of LS and give sensible free-surface BCs */
#ifndef SHEATH
	// Set average drop velocity to zero - this is often nescessary to keep the drop rotating around its axis
	double average_v = 0.0;
	double N_values = 0.0;
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX-1; j++) {
			if ( 0.5*(grid->LS[i][j]+grid->LS[i][j+1]) >= 0.0 ) {
				average_v += grid->v[i][j];
				N_values += 1.0;
			}
		}
	}
	average_v /= N_values;
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX-1; j++) {
			if ( 0.5*(grid->LS[i][j]+grid->LS[i][j+1]) >= 0.0 ) {
				grid->v[i][j] -= average_v;
			}
		}
	}
#endif

	// set free slip along r=0
	for (int j=0; j<=JMAX+1; j++) {
		grid->u[0][j] = 0.0;				// u = 0
		grid->u[IMAX][j] = 0.0;
	}

	double temp, err, maxerr;
	double N, S, E, W, R;
	do {
		maxerr = 0.0;
		for (int i=1; i<=IMAX; i++) {
			for (int j=1+(i&1); j<=JMAX; j+=2) {
				N = 1.0;
				S = 1.0;
				E = 1.0;
				W = 1.0;
				if ( i==1 ) W = 0.0;
				if ( i==IMAX ) E = 0.0;
				if ( grid->LS[i][j+1] >= 0.0 ) N = 0.0;
				if ( grid->LS[i][j-1] >= 0.0 ) S = 0.0;
				if ( grid->LS[i+1][j] >= 0.0 ) E = 0.0;
				if ( grid->LS[i-1][j] >= 0.0 ) W = 0.0;

				if ( (N+S+E+W) < 0.5 ) {
					// Cell surrounded by fluid: do nothing
				} else {		// Use RHS to store temporary values
					R = double(i)-0.5;
					temp = (1.0-RELAX)*grid->vel_pot[i][j] + RELAX *
					       ( (R+0.5)*E*grid->vel_pot[i+1][j] + (R-0.5)*W*grid->vel_pot[i-1][j]
							+ R*N*grid->vel_pot[i][j+1] + R*S*grid->vel_pot[i][j-1]
							+ (1.0-E)*(R+0.5)*grid->u[i][j]/INV_H
							- (1.0-W)*(R-0.5)*grid->u[i-1][j]/INV_H
							+ (1.0-N)*R*grid->v[i][j]/INV_H
							- (1.0-S)*R*grid->v[i][j-1]/INV_H )
							 	/ ( R*N + R*S + (R+0.5)*E + (R-0.5)*W );

					err = fabs(grid->vel_pot[i][j] - temp);
					if (err > maxerr) maxerr = err;
					grid->vel_pot[i][j] = temp;
				}
			}
		}
		maxerr = 0.0;
		for (int i=1; i<=IMAX; i++) {
			for (int j=2-(i&1); j<=JMAX; j+=2) {
				N = 1.0;
				S = 1.0;
				E = 1.0;
				W = 1.0;
				if ( i==1 ) W = 0.0;
				if ( i==IMAX ) E = 0.0;
				if ( grid->LS[i][j+1] >= 0.0 ) N = 0.0;
				if ( grid->LS[i][j-1] >= 0.0 ) S = 0.0;
				if ( grid->LS[i+1][j] >= 0.0 ) E = 0.0;
				if ( grid->LS[i-1][j] >= 0.0 ) W = 0.0;

				if ( (N+S+E+W) < 0.5 ) {
					// Cell surrounded by fluid: do nothing
				} else {		// Use RHS to store temporary values
					R = double(i)-0.5;
					temp = (1.0-RELAX)*grid->vel_pot[i][j] + RELAX *
					       ( (R+0.5)*E*grid->vel_pot[i+1][j] + (R-0.5)*W*grid->vel_pot[i-1][j]
							+ R*N*grid->vel_pot[i][j+1] + R*S*grid->vel_pot[i][j-1]
							+ (1.0-E)*(R+0.5)*grid->u[i][j]/INV_H
							- (1.0-W)*(R-0.5)*grid->u[i-1][j]/INV_H
							+ (1.0-N)*R*grid->v[i][j]/INV_H
							- (1.0-S)*R*grid->v[i][j-1]/INV_H )
							 	/ ( R*N + R*S + (R+0.5)*E + (R-0.5)*W );

					err = fabs(grid->vel_pot[i][j] - temp);
					if (err > maxerr) maxerr = err;
					grid->vel_pot[i][j] = temp;
				}
			}
		}
		if ( isnan(maxerr) == true ) cout << "Error in velocity potential solver" << endl;
	} while (maxerr > TOL);

	for (int i=1; i<=IMAX-1; i++) {
		for (int j=1; j<=JMAX; j++) {
			if ( grid->LS[i][j] < 0.0 || grid->LS[i+1][j] < 0.0 ) {
				grid->u[i][j] = (grid->vel_pot[i+1][j]-grid->vel_pot[i][j])*INV_H;
			}
		}
	}
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX-1; j++) {
			if ( grid->LS[i][j] < 0.0 || grid->LS[i][j+1] < 0.0 ) {
				grid->v[i][j] = (grid->vel_pot[i][j+1]-grid->vel_pot[i][j])*INV_H;
			}
		}
	}

	/* Free-slip BCs at r=0 and r=R; no-slip BCs at z = +-L */
	for (int j=0; j<=JMAX+1;j++) {				// free-slip
		grid->u[0][j] = 0.0;				// u = 0
		grid->v[0][j] = grid->v[1][j];			// dv/dr=0
		grid->u[IMAX][j] = 0.0;
		grid->v[IMAX+1][j] = grid->v[IMAX][j];
	}

	for (int i=0; i<=IMAX+1; i++) {  			// no-slip
		grid->v[i][JMAX] = 0.0;				// v = 0
		grid->u[i][JMAX+1] = (-1.0)*grid->u[i][JMAX];	// u = 0 (by averaging)
		grid->v[i][0] = 0.0;
		grid->u[i][0] = (-1.0)*grid->u[i][1];
	}

}


void advect_LS( Grid* grid, double DT) {
	// Use 2nd order upwinding and use RHS to store temporary values
	double DXP, DXM, DYP, DYM, U, V;		// [d(LS)/dx]^+ etc
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {

			DYM = (3.0*grid->LS[i][j] - 4.0*grid->LS[i][j-1] + grid->LS[i][j-2])*0.5*INV_H;
			if ( j==1 ) DYM = (grid->LS[i][1]-grid->LS[i][0])*INV_H;
			DYP = (-grid->LS[i][j+2] + 4.0*grid->LS[i][j+1] - 3.0*grid->LS[i][j])*0.5*INV_H;
			if ( j==JMAX ) DYP = (grid->LS[i][JMAX+1]-grid->LS[i][JMAX])*INV_H;

			DXM = (3.0*grid->LS[i][j] - 4.0*grid->LS[i-1][j] + grid->LS[i-2][j])*0.5*INV_H;
			if ( i==1 ) DXM = (grid->LS[1][j]-grid->LS[0][j])*INV_H;
			DXP = (-grid->LS[i+2][j] + 4.0*grid->LS[i+1][j] - 3.0*grid->LS[i][j])*0.5*INV_H;
			if ( i==IMAX ) DXP = (grid->LS[IMAX+1][j]-grid->LS[IMAX][j])*INV_H;

			U = (grid->u[i][j] + grid->u[i-1][j])/2.0;
			V = (grid->v[i][j] + grid->v[i][j-1])/2.0;
			/* use RHS to store temporary values of DLSDT */
			grid->RHS[i][j] = max(U,0.0)*DXM + min(U,0.0)*DXP + max(V,0.0)*DYM + min(V,0.0)*DYP;
		}
	}
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			grid->LS[i][j] -= DT*grid->RHS[i][j];
		}
	}

	/* impose 90 degree BCs */
	for (int i=1; i<=IMAX; i++) {
		grid->LS[i][0] = grid->LS[i][1];
		grid->LS[i][JMAX+1] = grid->LS[i][JMAX];
	}
	for (int j=0; j<=JMAX+1; j++) {
		grid->LS[0][j] = grid->LS[1][j];
		grid->LS[IMAX+1][j] = grid->LS[IMAX][j];
	}
}


int main() {

	// set up grid at t=0
	double t = 0.0;
	double plasma_t = 0.0;
	double integer_t = 0.0;
	double fluid_params_t = 0.0;
	Grid grid;
	init_liquid(&grid);
	init_plasma(&grid);


	// Set up sheath region
	bool nan_detected = false;
	do {
		set_plasma_dt(&grid);
		solve_PB_eqn(&grid);
		nan_detected = calc_ion_n(&grid);
		nan_detected = max( calc_ion_uv(&grid), nan_detected);
		plasma_t += grid.dt;
		cout << plasma_t << endl;
	} while (nan_detected == false && plasma_t < 10.0);

	write_grid(&grid, 0);
	double DT = set_dt(&grid);

	bool did_step_work = true;

	// main loop
	do {
		// write data to file if necessary
		if ( t >= integer_t* TIME_STEP ) {
			write_grid(&grid, integer_t);	// NB real time is t, integer_t will have rounding error
			integer_t += 1.0;		// but integer_t is better for filenames
		}
		if ( t >= fluid_params_t ) {
			write_timeseries(&grid, t);
			fluid_params_t += 0.01;
		}

		// mini sheath loop
		plasma_t = 0.0;
		do {
			set_plasma_dt(&grid);
			solve_PB_eqn(&grid);
			nan_detected = calc_ion_n(&grid);
			nan_detected = max( calc_ion_uv(&grid), nan_detected);
			plasma_t += grid.dt;
		} while (nan_detected == false && plasma_t < 0.1);

		// update liquid region
		did_step_work = step_grid(&grid, DT);
		if ( did_step_work == false ) {
			cout << "Terminating simulation: Poisson solver returned nan" << endl;
			write_grid(&grid, t);
			return 0;
		}
		velocity_extension_and_BCs(&grid);
		advect_LS(&grid, DT);

		// re-set dt at end of loop ready for next iteration
		DT = set_dt(&grid);
		t += DT;

		// tell me where the simulation is up to
		cout << "T = " << t << endl;
	} while (nan_detected == false && t <= FINAL_TIME+DT);

	return 0;
}
