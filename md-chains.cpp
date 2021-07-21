
/*----------------------------------------------------------------

                    Simple MD-Simulation.

        Used in parallel to test new installations into
    the University of Reading's Polymer Physics group's overall
              code to test different components.

  ----------------------------------------------------------------*/

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "md-chains.h"
#include "correl.h"
#include "rand.h"

#if defined UMBRELLA_SAMPLING
#define USE_MPI
#define SIMPLE_SPRNG
#include "sprng.h"
#include "mpi.h"
#endif

/*----------------------------------------------------------------*/
/* Centre of mass calculation and save it into a structure */
void centre_of_mass() {

  for (int i = 0; i < poly.N; i++) 
  {
    /* (Re)initialise the centre of mass structure. */
    for (unsigned int j = 0; j < 3; j++)
    {
      mole[i].cm[j] = 0.0;
    }

    /* Calculate the centre of mass for molecule i */
    for (int j = 0; j < poly.MPC; j++)
    {
      int ip = j + i * poly.MPC;

      for (unsigned int k = 0; k < 3; k++)
      {
        mole[i].cm[k] += part[ip].r[k];
      }
    }
    
    /* Normalise the centre of mass of molecule i */
    for (unsigned int j = 0; j < 3; j++)
    {
      mole[i].cm[j] /= (double)poly.MPC;
    }
  }

}

/*----------------------------------------------------------------*/
/* GLOBAL PARAMETERS. */

long idum;
long *idumPtr = &idum;

int istep, MAX_blist, numtasks, rankConst, rc, swapFlag;
double box_l[3], volume, Temp, dt, Gamma;
bool WARMEDUP;

int acceptedTemperingCounter[50];
int allTemperingCounter[50];
bool BIN_FLAG;

Sigma sT;

#define PI 3.14159265358979323846264338328

/*----------------------------------------------------------------*/

/* Initialise chain. */
void initChain();

/* Rescale the coordinates and velocity of the particles */
void rescaleParticle();

/* Calculate all the forces - no Lennard-Jones */
void forces_calcul();

/* The Velocity-Verlet integrator */
void integrator();

/* Euler method to integration */
void euler();

/* Calculate the average bond length throughout the simulation which will
   be printed out at the end rather than saved continuously... */
double calcBondLength();

double calcE2Esq();

void paraPrint(FILE *sys_p_Ptr, int ns_warm, int nstep, int nrelax);

void write_out_n();

void harmonic_CM_forces_U(void);

int temperingSweep(double E, double *th, int N, double *TA, int istep);

void r_neg1();

void umbrellaShift(double shift);

/*----------------------------------------------------------------*/
/* Main function - initialisation + MD-loop. */
int main(int argc, char **argv) {

int anaCount;
double ttime, bondDist, e2eDist, *TA;
char *file_name;

  #if defined USE_MPI
  rc = MPI_Init(&argc, &argv);
  rc = MPI_Comm_rank(MPI_COMM_WORLD, &rankConst);
  rc = MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  TA = (double *) malloc(sizeof(double) * numtasks);

  if (rankConst == 0) {
    std::cout << "There is " << numtasks << " task(s) being computed right now.\n";
  }
  #endif


  if (argc < 2) {
    std::cout << "Please put in the number of polymers.\n";
  } else {
    poly.N = atoi(argv[1]);

    if (argc < 3) {
      std::cout << "Please put in the number of monomers within each molecule.\n";
    } else {
      poly.MPC = atoi(argv[2]);

      if (argc < 4) {
        std::cout << "Please put in the file-name of the system.\n";
      } else {
        file_name = argv[3];
      }
    }
  }

  std::cout << "The number of particles in the system is equal to " 
            << poly.N * poly.MPC << std::endl;

          
  /* Set for each processor so we can distinguish each one + debug. */
  #if defined UMBRELLA_SAMPLING
  sprintf(file_name, "%s_%d", file_name, rankConst);

  init_sprng(DEFAULT_RNG_TYPE, idum, SPRNG_DEFAULT);
  print_sprng();
  idum = -(long)time(NULL);

  /* Effectively parallelises idum for each case. */
  double V = sprng();
  idum *= V;
  std::cout << "idum.parallelized................................ = U(0,1) = " << idum << std::endl;

  BIN_FLAG = false;
  um.k = 2.0e+04;
  um.size = 5000;
  um.spring_len = 1.0;
  um.range = {0.0, 10.0};

  um.d_eq = um.range[0]
          + ((double)rankConst + 0.5) * fabs(um.range[1] - um.range[0]) / (double)numtasks;

  std::cout << "The equilibrium distance for the harmonic centre of mass force is "
            << um.d_eq << ". For simulation number " << rankConst << "." << std::endl;

  um.xVal = new double [um.size] ();
  um.binn = new int [um.size] ();

  for (unsigned int i = 0; i < um.size; i++)
  {
    um.xVal[i] = um.range[0] + (i + 0.5) * (um.range[1] - um.range[0]) / (double)um.size;
  }

  #else
  idum = -(long)time(NULL);
  #endif

  /* Global parameters. */
  MAX_blist = 1;
  dt    = 1e-03;
  Gamma = 1.0;
  Temp  = 1.0;

  /* Minimising global parameters whenever possible. */
  const int nstep    = 6e+07;
  const int warmstep = 1e+07;
  const int nrelax   = 1e+07;
  const double Rho   = 1.e-03;

  volume = (double)(poly.N * poly.MPC) / Rho;
  box_l[0] = box_l[1] = box_l[2] = pow(volume,1./3.);
  std::cout << "Box length is " << box_l[0] << std::endl;

  /* Dynamically allocate the particle data array. */
  part = (Particle *) malloc(sizeof(Particle) * poly.N * poly.MPC);
  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    part[i].nb = 0;
    part[i].bl = new int [MAX_blist];

    for (int j = 0 ; j < MAX_blist; j++)
    {
      part[i].bl[j] = -1;
    }
  }

  mole = (polymer_whole_struct *) malloc(sizeof(polymer_whole_struct) * poly.N);

  char SysPara_FileName[100];
  FILE *SysPara_Ptr;
  sprintf(SysPara_FileName, "sys_para_%s.dat", file_name); 
  if ((SysPara_Ptr = fopen(SysPara_FileName, "w")) == NULL) {
    printf("File %s could not be opened!\n", SysPara_FileName);
    exit(1);
  
  } else {
    paraPrint(SysPara_Ptr,warmstep,nstep,nrelax);
    fclose(SysPara_Ptr);
  }

  /* Initialisation of the particles. */
  initChain();

  /* Rescale. */
  #if defined UMBRELLA_SAMPLING
  umbrellaShift(2.5);
  #endif
  rescaleParticle();

  /* Initial force calulation + Warm-up loop -
     normally used for Lennard Jones interactions in the main code. */
  WARMEDUP = false;
  forces_calcul();

  std::cout << "Begin warming up.\n";
  istep = 0, anaCount = 0;
  do {

    istep++;
    if ( !(istep % (int)1e+05) ) {
      std::cout << "The warm-up step is " << istep << std::endl;
    }

    /* Integrate the degrees of freedom of the system. */
    #if defined VERLET
      integrator();
    #else
      euler();
    #endif
    
  } while (istep < warmstep);
  WARMEDUP = true;
  std::cout << "The warm-up has been completed!!\n\n";


  /*
  s1 = (int) ceil(log2((nstep-nrelax) / p) + 1);
  s1 = (int)( log((nstep - nrelax) / p) / log(2.0) + 1.0 );
  */

  Corrl stress;
  stress.constitute(1, 1.0 / (3.0 * (double)(poly.N * (poly.MPC-1))), "stress");
  //stress.constitute(1, 1./(3. * box_l[0] * box_l[1] * box_l[2]), "stress");

  Corrl *e2e = new Corrl [poly.N];
  for (int i = 0; i < poly.N; i++) e2e[i].constitute(1,1.0,"end-to-end");

  Corrl *vel = new Corrl [poly.N];
  for (int i = 0; i < poly.N; i++) vel[i].constitute(2,1.0,"velocity");

  Corrl *rot = new Corrl [poly.N];
  for (int i = 0; i < poly.N; i++) rot[i].constitute(2,1.0,"rotation");

  Corrl c2c;
  c2c.constitute(2,1.0,"com-to-core");


  std::cout << "Begin the MD-loop.\n";

  time_t time1, time2; 
  time1 = time(NULL);

  istep = 0, ttime = 0.0, bondDist = 0.0, e2eDist = 0.0;
  do {
    istep++;
    ttime += dt;

    if ( !(istep % (int)1e+05) ) {
      std::cout << "The MD-step is " << istep << std::endl;
    }

    /* Integrate the degrees of freedom of the system. */
    #if defined VERLET
    integrator();
    #else
    euler();
    #endif

    /* Any analysis is done here underneath here. */
    if (istep > nrelax) {
      // Stress correl calc at every time-step due to high the fluctuations.
      if (istep % stress.grab_wait() == 0) mainStressCorrl(stress);

      if (poly.N * poly.MPC > 0) {
        for (int i = 0; i < poly.N; i++)
        {
          if (istep % e2e[i].grab_wait() == 0) {
            mainEndToEndCorrl(e2e[i], i * poly.MPC, (i+1) * poly.MPC - 1);
          }

          if (istep % vel[i].grab_wait() == 0) {
            mainVelocityCorrl(vel[i],i);
          }
        }
      }
      
      anaCount++;
      bondDist += calcBondLength();
      e2eDist += calcE2Esq();
    }

    #if defined UMBRELLA_SAMPLING
      if (istep == nrelax) BIN_FLAG = true;

      if (BIN_FLAG) {
        if (istep % (int)2.5e+06 == 0) {
          temperingSweep(0.0, (double *)part, poly.N * poly.MPC * sizeof(Particle) / sizeof(double), TA, istep);
        }

        if (istep % (int)5e+06 == 0) {
          write_out_n();

          #if defined DEBUG
          double dist, distAv = 0.0;

          for (int i = 0; i < poly.N-1; i++)
          {
            for (int j = i+1; j < poly.N; j++)
            {
              dist = 0.0;
              for (unsigned int k = 0; k < 3; k++)
              {
                dist += SQR(mole[i].cm[k] - mole[j].cm[k]);
              }

              distAv += sqrt(dist);
            }
          }

          if (poly.N > 1) {
            double div = 0.0;
            for (double i = 0; i < poly.N-1; i++) div *= i+1;
            distAv /= div;
          }

          std::cout << "THE AVERAGE DISTANCE IS " << distAv << std::endl;
          #endif
        }
      }
    #endif

  } while (istep < nstep);

  /* Calculate the simulation time accordingly... */
  time2 = time(NULL);
  fprintf(stderr, "The total run-time of the md-loop is %le\n", difftime(time2,time1));

  std::cout << "The MD-loop has been completed.\nPrinting out the relevant data...\n\n";


  /* OPERATOR TAU - print out to a file for later plotting using MATLAB file. */
  FILE *acf_FilePtr;
  char acfFileName[100];
  sprintf(acfFileName, "acf_%s.dat", file_name);
  if ((acf_FilePtr = fopen(acfFileName, "w")) == NULL) {
    printf ("File %s has not been generated! \n", acfFileName);
    exit(1);

  } else {
    std::cout << "Opened the auto-correlation function save file\n";

    stress.calcCorrl(acf_FilePtr);

    if (poly.N * poly.MPC > 0) {
      collateCorrlData(acf_FilePtr, e2e, poly.N);
      collateCorrlData(acf_FilePtr, vel, poly.N);
    }

  }
  fclose(acf_FilePtr);

  std::cout << "The average bond-length is " << bondDist / (double)anaCount << std::endl;
  std::cout << "The average end-to-end distance is " << e2eDist / (double)anaCount << std::endl;


  /* Deallocate memory-step for clean-up. */
  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    delete [] part[i].bl;
  }

  #if defined USE_MPI
  rc = MPI_Finalize();
  delete[] um.xVal;
  delete[] um.binn;
  #endif

  free(part);
  free(mole);

  printf("Fin simulation run.\n");
  return 0;
}

/*----------------------------------------------------------------*/
/* Initialise chain. */
void initChain() {

int ip;
double radius = 1.0, d_constr, theta, phi, dist, dr[3];
bool checkFlag;

  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    do {
      checkFlag = true;

      if ( !(i % poly.MPC) ) {

        for (unsigned int j = 0; j < 3; j++)
        {
          part[i].r[j] = (ran2(idumPtr) - 0.5) * box_l[j];
        }

      } else {
        d_constr = radius * (ran2(idumPtr) + 0.5);
        theta    = acos(2.0 * ran2(idumPtr) - 1.0);
        phi      = ran2(idumPtr) * 2.0 * PI;

        part[i].r[0] = part[i-1].r[0] + d_constr * sin(theta) * cos(phi);
        part[i].r[1] = part[i-1].r[1] + d_constr * sin(theta) * sin(phi);
        part[i].r[2] = part[i-1].r[2] + d_constr * cos(theta);

        /* Make sure it lies in the box. */
        for (unsigned int j = 0; j < 3; j++)
        {
          if (fabs(part[i].r[j]) >= 0.5 * box_l[j]) {
            checkFlag = false;
            break;
          }
        }
      }

      /* Make sure that there exists no overlap */
      if (checkFlag && i != 0) {
        for (ip = i-1; ip >= 0; ip--)
        {
          dist = 0.0;
          for (unsigned int j = 0; j < 3; j++)
          {
            dr[j] = part[i].r[j] - part[ip].r[j];
            dist += dr[j] * dr[j];
          }
          dist = sqrt(dist);

          if (dist < 0.4) {
            checkFlag = false;
            break;
          }
        }
      }
    } while (!checkFlag);

    for (unsigned int j = 0; j < 3; j++)
    {
      part[i].v[j] = sqrt(Temp) * ran2(idumPtr);
    }
  }

  /* Set the bond list. */
  for (int i = 0; i < poly.N; i++)
  {
    for (int j = 0; j < poly.MPC-1; j++)
    {
      ip = i * poly.MPC + j;

      part[ip].nb = 1;
      part[ip].bl[0] = ip+1;
    }
  }

}

/*------------------------------------------------------*/
/* Rescale the coordinates and velocity of the particles */
void rescaleParticle() {

double v_sum[3] = { 0.0, 0.0, 0.0 };
double kine_eng = 0.0 ;

  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      v_sum[j] += part[i].v[j];
    }
  }

  for (unsigned int j = 0; j < 3; j++)
  {
    v_sum[j] /= (double)(poly.N * poly.MPC);
  }

  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      part[i].v[j] -= v_sum[j];
      kine_eng += SQR(part[i].v[j]);
      part[i].v[j] *= dt;
    }
  }

  std::cout << "The kinetic energy is equal to " << 0.5 * kine_eng << std::endl;

}

/*---------------------------------------------------------------------------*/
/* Harmonic potential - strictly for Gaussian chains! */
void harmonic_force() {

int    j, l, t1, t2;
double d[3], dist2, spring_const;

  double bondl = 1.0;
  spring_const = 3.0 * Temp / bondl / bondl;
  
  for (t1 = 0; t1 < poly.N * poly.MPC; t1++)
  {
    if (part[t1].nb != 0) {
      for (j = 0; j < part[t1].nb; j++)
      {
        t2 = part[t1].bl[j];
        
        if (t1 < t2 && t2 != -1) {
          dist2 = 0.0;
          for (l = 0; l < 3; l++)
          {
            d[l]   = part[t1].r[l] - part[t2].r[l];
            dist2 += SQR(d[l]);

            part[t1].f[l] -= spring_const * d[l];
            part[t2].f[l] += spring_const * d[l];
          }

          if (WARMEDUP) addToStress(d,spring_const);
        }

      }
    }
  }

}

/*----------------------------------------------------------------------*/
/* FENE potential [U_fene = -1/2 * k * R_0^2 * ln(1 - r^2/R_0^2)] */
void fene_force(void) {
/* Expects in the FOLDED particle coordinates?
   The routine increments the forces with the fene
   contributions. If the pointer energy is given, it is assumed that the 
   FENE contribution to the energy is to be computed as well. */

int    i, j, ii, jj, t1, t2;
double d[3], r2, fac, pre1, pre2;
 
  pre1 = 0.5 * fene.k * fene.r2;

  energy.FENE = 0.0;

  for (ii = 0; ii < poly.N; ii++)
  {
    for (jj = 0; jj < poly.MPC; jj++)
    {
      t1 = jj + ii * poly.MPC;

      if (part[t1].nb != 0) {
        for (i = 0; i < part[t1].nb; i++)
        {
          t2 = part[t1].bl[i];
          
          if (t1 < t2) {
    
            r2 = 0.0;
            for (j = 0; j < 3; j++)
            {
              d[j]  = part[t1].r[j] - part[t2].r[j];
              d[j] -= dround(d[j]/box_l[j]) * box_l[j];
              r2   += SQR(d[j]);
            }
      
            fac = fene.k / (1.0 - r2 / fene.r2);
            energy.FENE -= pre1 * log(1.0 - r2 / fene.r2);
            
            for (j = 0; j < 3; j++)
            {
              part[t1].f[j] -= fac * d[j]; 
              part[t2].f[j] += fac * d[j]; 
            }

            if (WARMEDUP) addToStress(d,fac);
          }
        }
      }
    }
  }

}

/*------------------------------------------------------------*/
/* Calculate all the forces - no Lennard-Jones */
void forces_calcul() {

  centre_of_mass();

  if (WARMEDUP) sT.reinit();
  
  /* Initialise the forces s.t. */
  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      part[i].f[j] = 0.0;
    }
  }

  /* Bonded forces consideration. */
  harmonic_force();
  //fene_force();

  #if defined UMBRELLA_SAMPLING
  harmonic_CM_forces_U();
  #endif

  /* Unbonded interactions are calulcated below ... */
  r_neg1();

  /* Thermal + frictional forces for each particle... */
  static double Constant = sqrt(24.0 * Gamma * Temp / dt);

  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      part[i].f[j] -= Gamma * part[i].v[j] / dt;
      part[i].f[j] -= Constant * (ran2(idumPtr) - 0.5);
    }
  }

  /* Transfer the forces to REDUCED UNITs by multiplying dt^2/2 */
  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      part[i].f[j] *= 0.5 * dt * dt;
    }
  }

}

/*---------------------------------------------------*/
/* The Velocity-Verlet integrator */
void integrator() {

  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      part[i].v[j] += part[i].f[j];
      part[i].r[j] += part[i].v[j];
    }
  }

  /* Calculate the forces at t + dt ... */
  forces_calcul();

  /* Propagate the velocity to t + dt */
  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      part[i].v[j] += part[i].f[j];
    }
  }

}

/*---------------------------------------------------*/
/* Euler's method of integration */
void euler() {

  static double C = sqrt(12.0 * 2.0 * Gamma * Temp / dt);

  if (WARMEDUP) sT.reinit();

  /* Initialise the forces. */
  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      part[i].f[j] = 0.0;
    }
  }

  forces_calcul();

  for (int i = 0; i < poly.N * poly.MPC; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      part[i].v[j] = (part[i].f[j] + C * (ran2(idumPtr) - 0.5)) / Gamma;
      part[i].r[j] += dt * part[i].v[j];
    }
  }

  centre_of_mass();

}

/*----------------------------------------------------------------*/
/* Calculate the average bond length throughout the simulation which
   will be printed out at the end rather than saved continuously. */
double calcBondLength() {

int countBonds = 0;
double dist = 0.0;

  for (int t1 = 0 ; t1 < poly.N * poly.MPC; t1++)
  {
    if (part[t1].nb > 0) {
      for (int i = 0; i < part[t1].nb; i++)
      {
        unsigned int t2 = part[t1].bl[i];

        double dist2Tmp = 0.0;
        for (unsigned int j = 0; j < 3; j++)
        {
          dist2Tmp += SQR(part[t2].r[j] - part[t1].r[j]);
        }

        dist += sqrt(dist2Tmp);
        countBonds++;
      }
    }
  }

  return dist / (double)countBonds;
}

/*----------------------------------------------------------------*/
/* End-to-end distance of the chain - printed out at the end. */
double calcE2Esq() {

  double dist2 = 0.0;
  for (int i = 0; i < poly.N; i++)
  {
    int t1 = i * poly.MPC;
    int t2 = (i + 1) * poly.MPC - 1;  

    if (t1 != t2) {
      double dist2Tmp = 0.0;
      for (unsigned int j = 0; j < 3; j++)
      {
        dist2Tmp += SQR(part[t2].r[j] - part[t1].r[j]);
      }

      dist2 += dist2Tmp;
    }
  }

  return dist2 / (double)poly.N;
}


/*----------------------------------------------------------------*/
/* Print paramters in the same vein as the main md-code. */
void paraPrint(FILE *sys_p_Ptr, int ns_warm, int nstep, int nrelax) {

  fprintf(sys_p_Ptr, "Number.of.data.points.for.analysis.....= 0\n");
  fprintf(sys_p_Ptr, "----N/A----\n");
  fprintf(sys_p_Ptr, "----N/A----\n");
  fprintf(sys_p_Ptr, "----N/A----\n");
  fprintf(sys_p_Ptr, "----N/A----\n");
  fprintf(sys_p_Ptr, "----N/A----\n");
  fprintf(sys_p_Ptr, "Time.step..............................= %le\n",dt);
  fprintf(sys_p_Ptr, "Number.of.time.steps.for.warming.up....= %d\n", ns_warm);
  fprintf(sys_p_Ptr, "Total.number.of.time.steps.............= %d\n", nstep);
  fprintf(sys_p_Ptr, "Total.number.of.steps.for.relaxation...= %d\n", nrelax);
  fprintf(sys_p_Ptr, "----N/A----\n");
  fprintf(sys_p_Ptr, "----N/A----\n");
  fprintf(sys_p_Ptr, "Time.steps.to.tape.the.configuration...= 0\n");

}

#if defined UMBRELLA_SAMPLING
/*--------------------------------------------------------------------*/
/* Tempering sweep function. */
int temperingSweep(double E, double *th, int N, double *TA, int istep) {

static int direction = -1;
static int numtasks;
const static int dM = 1;
static double T;
static double *thbuf;

  std::cout << "Begin tempering sweep for "
            << rankConst << " at " << istep << std::endl;

  if (rankConst == -1) {
    thbuf = new double[N+2];
    T = TA[rankConst];
  }

  static MPI_Status status;

  direction *= -1;
  int partner =
    ((rankConst / dM) % 2) ? (rankConst + dM * direction) : (rankConst - dM * direction);

  MPI_Barrier(MPI_COMM_WORLD);


  if (partner >= 0 && partner < numtasks) {
   
    const double r1 = sprng();
    th[N] = E;
    th[N+1] = r1;
    const int l = N+2;
    const int ID1 = partner;
    const int ID2 = rankConst;

    MPI_Sendrecv(th, l, MPI_DOUBLE, partner, ID1, thbuf, l, MPI_DOUBLE,
                 partner, ID2, MPI_COMM_WORLD, &status);

    const double Epartner = thbuf[N];
    const double r2 = thbuf[N+1];

    const double bE = -(1.0 / T - 1.0 / TA[partner]) * (E - Epartner);
    const double P  = (bE > 0) ? exp (-bE) : 1.0;
    const double r  = (rankConst > partner) ? r1 : r2;

    if (r < P) {
      int ic;
      for (ic = 0; ic < N; ic++) th[ic] = thbuf[ic];
      acceptedTemperingCounter[dM]++;
      swapFlag = 1;
    }

    allTemperingCounter[dM]++;


    std::cout << " istep: " << istep << " accept: " << acceptedTemperingCounter[dM] 
              << " out of " << allTemperingCounter[dM]<< " ....... rank "
              << rankConst << std::endl;
    
    return 1;
  }

  return 0;
}

/*--------------------------------------------------------------------*/
/* Write out the bin. This unit will write out the collective
   bin variable to file. This is initiated given a certain LONG
   time-step (can see the exact in main.c within the md-loop). */
void write_out_n() {

  double EXPECTED   = um.xiBar / (double)um.count;
  double EXPECTED2  = um.xi2Bar / (double)um.count;
  double SDEVIATION = sqrt(EXPECTED2 - SQR(EXPECTED));

/*
  std::ofstream binFile;
  std::string binFileName = "n_" + IntToStr(rankConst) + ".dat";
  binFile.open(binFileName.c_str());

  binFile << boost::format("%1% %2% %3% %4% %5% %6% %7%\n")
            % Temp % um.spring_len % um.k % um.range[0] 
            % um.d_eq % EXPECTED % sqrt(VARIANCE);
  
  for (unsigned int i = 0; i < um.size; i++)
  {
    binFile << um.xVal[i] << " " << um.binn[i] << std::endl;
  }
  
  binFile.close();
*/

  FILE *bin_FilePtr;
  char bin_FileName[100];
	sprintf(bin_FileName, "n_%d.dat", rankConst);
  if ((bin_FilePtr = fopen(bin_FileName, "w")) == NULL) {
    printf("File %s has not been generated!\n", bin_FileName);
    exit(1);
  
  } else {

    std::cout << "Begin to write out the umbrella bins.\n";

    fprintf(bin_FilePtr, "%f %f %f %f %f %f %f %f\n",
            Temp, um.k, um.K, um.range[0], um.range[1],
            um.d_eq, EXPECTED, SDEVIATION);
    
    for (unsigned int i = 0; i < um.size; i++)
    {
      fprintf(bin_FilePtr, "%f %d\n", um.xVal[i], um.binn[i]);
    }
    
    fclose(bin_FilePtr);

  }

}

#endif

/*------------------------------------------------------------------------*/
/* Add stress to the tensor, (the class) Sigma sT */
void addToStress(double d[3], double spring_const)
{
  sT.addIn1 ( -d[0] * spring_const * d[1] );
  sT.addIn2 ( -d[1] * spring_const * d[2] );
  sT.addIn3 ( -d[2] * spring_const * d[0] );
}

/*------------------------------------------------------------------------*/
/* Harmonic force wrt. CENTRE OF MASS used for Umbrella Sampling */
void harmonic_CM_forces_U() {

unsigned int t1, t2, bin;
double d[3], dist;
double tmp[3];

  static double increment = fabs(um.range[1] - um.range[0]) / (double)um.size;

  if (poly.N > 1) {
    for (int i = 0; i < poly.N - 1; i++)
    {
      for (int j = i+1; j < poly.N; j++)
      {

        dist = 0.0;
        for (unsigned int k = 0; k < 3; k++)
        {
          d[k]  = mole[i].cm[k] - mole[j].cm[k];
          dist += SQR(d[k]);
        }
        dist = sqrt(dist);

        if (BIN_FLAG) {
          bin = (unsigned int)floor((dist - um.range[0]) / increment);

          if (bin < um.size) um.binn[bin]++;

          if (um.range[0] <= dist && dist <= um.range[1]) {
            um.xiBar  += dist;
            um.xi2Bar += SQR(dist);
            um.count ++;
          }
        }
        
        /* Add to the CoM force to each particle: */
        for (int l = 0; l < poly.MPC; l++)
        {
          t1 = i * poly.MPC + l;
          t2 = j * poly.MPC + l;

          part[t1].f[0] -= um.k * (d[0] - um.d_eq) / (double)poly.MPC;
          part[t2].f[0] += um.k * (d[0] - um.d_eq) / (double)poly.MPC;
        
          for (unsigned int k = 1; k < 3; k++)
          {
            part[t1].f[k] -= um.K * d[k] / (double)poly.MPC;
            part[t2].f[k] += um.K * d[k] / (double)poly.MPC;
          }

          if (WARMEDUP) {
            tmp[0] = um.k * (d[0] - um.d_eq);
            tmp[1] = um.K * d[1];
            tmp[2] = um.K * d[2];

            addToStress(tmp,1.0);
          }
        }
        
      }
    }
  }

}

/*------------------------------------------------------------------------*/
/* r-1 interactions - unbonded interactions between monomers. */
void r_neg1() {
double fr1, dist, frac, d[3], force[3];

  if (poly.N * poly.MPC > 1) {
    for (int t1 = 0; t1 < poly.N * poly.MPC - 1; t1++)
    {
      for (int t2 = t1 + 1; t2 < poly.N * poly.MPC; t2++)
      {
        
        dist = 0.0;
        for (unsigned int k = 0; k < 3; k++)
        {
          d[k]  = part[t1].r[k] - part[t2].r[k];
          dist += SQR(d[k]);
        }
        dist = sqrt(dist);

        frac = 1.0 / dist;
        fr1  = frac / SQR(dist);

        for (unsigned int k = 0; k < 3; k++)
        {
          part[t1].f[k] += fr1 * d[k];
          part[t2].f[k] -= fr1 * d[k];
        }

        if (WARMEDUP) {
          addToStress(d,-fr1);
        }

      }
    }
  }

}

/*--------------------------------------------------------------------*/
/* Umbrella shift - initiate at the origin of the box then push away -
   add the radius of gyration to it to make sure its not penetrating */
void umbrellaShift(double shift) {
int ip, t1, t2;

  if (poly.N != 2) {
    printf("Umbrella shift not optimised for poly.N != 2");
    exit(1);
  }

  centre_of_mass();

  for (int i = 0; i < poly.N; i++)
  {
    for (int j = 0; j < poly.MPC; j++)
    {
      ip = i * poly.MPC + j;

      for (unsigned int k = 0; k < 3; k++)
      {
        part[ip].r[k] -= mole[i].cm[k];
      }
    }
  }

  double r_G = 0.0;
  for (int i = 0; i < poly.N; i++)
  {
    for (int j = 0; j < poly.MPC; j++)
    {
      ip = i * poly.MPC + j;

      for (unsigned int k = 0; k < 3; k++)
      {
        r_G += SQR(part[i].r[j] - mole[i].cm[k]);
      }
    }
  }

  for (t1 = 0; t1 < poly.MPC; t1++)
  {
    part[t1].r[0] -= shift + r_G;

    t2 = t1 + poly.MPC;
    part[t2].r[0] += shift + r_G;

    for (unsigned int k = 1; k < 3; k++)
    {
      part[t1].r[k] = part[t2].r[k];
    }
  }

  printf("Shift the molecules for the umbrella analysis.\n");

}