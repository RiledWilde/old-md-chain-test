
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "md-chains.h"
#include "correl.h"

/* Main correlation function for the stress tensor */
void mainStressCorrl(Corrl &cor);

/* Main correlation function for the rotational vector ON THE FLY */
void mainRotationalCorrl(Corrl &cor);

/* Main correlation function for the velocity (3-components) */
void mainVelocityCorrl(Corrl &cor);

/* Main correlation function for the end-to-end vector - unit to
   calculate the end to end autocorrelation relaxation function */
void mainEndToEndCorrl(Corrl &cor);

/* Main correlation function comparing centre of mass to the centre of
   the dendrimer - see analysis file for function comparing both */
void mainComToCoreCorrl(Corrl &cor);

/*
template<typename Base, typename T>
inline bool instanceof(const T*) { return std::is_base_of<Base, T>::value; }
*/

/* Collate the data within one averaged- correlator to save to file */
void collateCorrlData(FILE *FilePtr, Corrl cor[], int iter);

/* Calculate the correlation function */
void printCorrl(Corrl cor, FILE *acf_FilePtr, int normArr[], double corrlArr[]);

/*--------------------------------------------------------------------*/
/* Main correlation function for the stress tensor ON THE FLY. */
void mainStressCorrl(Corrl &cor) {

  cor.transfer(s_T);

  //std::cout << s_T << std::endl;
  //std::cout << cor << std::endl;

  cor.OperatorTau(false,true);
  
}

/*--------------------------------------------------------------------*/
/* Main correlation function for the rotational vector ON THE FLY */
void mainRotationalCorrl(Corrl &cor, bool normalise, int im1) {

int ip, count = 0;
double Q[3], abs = 0.0;
bool flagCont = false;

  cor.reinit();

  /* Create the list of the particle types we want to consider */
  static int list[] = {0,1};
  static size_t sizelist = SO(list);
  
  for (int i = 0; i < poly.MPC; i++)
  {
    ip = im1 * poly.MPC + i;

    flagCont = false;
    for (unsigned int k = 0; k < (unsigned int)sizelist; k++)
    {
      if (list[k] == part[ip].part_type) {
        flagCont = true;
        break;
      }
    }


    if (flagCont) {
      
      if (normalise) {
        abs = 0.0;
        for (unsigned int k = 0; k < 3; k++)
        {
          Q[k] = part[ip].r[k] - mole[im1].centre[k];
          abs += SQR(Q[k]);
        }

        if (abs > 0.0) {
          for (unsigned int k = 0; k < 3; k++)
          {
            cor.t_T[k] += Q[k] / sqrt(abs);
          }

          count++;
        }

      } else {
        for (unsigned int k = 0; k < 3; k++)
        {
          cor.t_T[k] += part[ip].r[k] - mole[im1].centre[k];
        }
        
        count++;
      }

    }
  }

  for (unsigned int k = 0; k < 3; k++)
  {
    if (count != 0) cor.t_T[k] /= (double)count;
  }

  cor.OperatorTau(true,true);

}

/*--------------------------------------------------------------------*/
/* Main correlation function for the velocity (3-components), for every
   monomer in the molecules ... */
void mainVelocityCorrl(Corrl &cor, int im1) {

  cor.reinit();

  for (int i = 0; i < poly.MPC; i++)
  {
    int ip = im1 * poly.MPC + i;

    cor.t_T[0] += part[ip].v[1] * part[ip].v[0];
    cor.t_T[1] += part[ip].v[2] * part[ip].v[0];
    cor.t_T[2] += part[ip].v[2] * part[ip].v[1];
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    cor.t_T[i] /= (double)poly.MPC;
  }

  cor.OperatorTau(true,true);

}

/*--------------------------------------------------------------------*/
/* Main correlation function for the end-to-end vector - unit to
   calculate the end to end autocorrelation relaxation function */
void mainEndToEndCorrl(Corrl &cor, int t1, int t2) {

  cor.reinit();

  /* Absolute end to end calculation located here */
  for (unsigned int j = 0; j < 3; j++)
  {
    cor.t_T[j] += part[t1].r[j] - part[t2].r[j];
  }

  cor.OperatorTau(true,true);

}

/*--------------------------------------------------------------------*/
/* Main correlation function comparing centre of mass to the centre of
   the dendrimer - see analysis file for function comparing both */
void mainComToCoreCorrl(class Corrl &cor) {

  cor.reinit();

  CentreOfMassToCentre(cor.t_T);

  cor.OperatorTau(true,true);

}

/*--------------------------------------------------------------------*/
/* Push correlation right */
void Corrl::pushCorrl(unsigned int is, bool smooth) {
unsigned int icor;

  /* Push the is correlator to the right and update the is+1 accumulator */
  for (unsigned int ip = p-1; ip > 0; ip--)
  {
    icor = p * is + ip;

    for (unsigned int i = 0; i < 3; i++)
    {
      tcor_T[i][icor] = tcor_T[i][icor-1];
    }
  }

  /* Which level are we on? */
  icor = p * is;
  if (smooth == false || is == 0) {
    for (unsigned int i = 0; i < 3; i++)
    {
      tcor_T[i][icor] = t_T[i];
    }
  
  } else {
    for (unsigned int i = 0; i < 3; i++)
    {
      tcor_T[i][icor] = tacc_T[i][is] / (double)m;
    }
  }

  /* Accumulator the next second correlator - as long as
     it is not the last level of the correlator. */
  if (is < s-1) {

    tcou[is+1]++;

    /* Smoothed out average (normal version employed by most people) */
    if (smooth == true) {

      for (unsigned int i = 0; i < 3; i++)
      {
        tacc_T[i][is+1] += tcor_T[i][icor];
      }

    /* Will mean the systematic error is void, and contains an acceptable
       statistical error. A.L. said this can be used to get correlation
       functions using higher moments of displacement efficiently. */
    } else if (tcou[is+1] == m) {
      
      for (unsigned int i = 0; i < 3; i++)
      {
        tacc_T[i][is+1] = tcor_T[i][icor];
      }
    }
  }


  /* RESET the element values of the arrays s.t. ... */
  tcou[is] = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    tacc_T[i][is] = 0.0;
  }

}

/*--------------------------------------------------------------------*/
/* Calculate the correlation function */
void Corrl::calc_ACF_Corrl(unsigned int is, unsigned int ip0, bool normalised) {
int icor;

  int icor0 = p * is;
  double normVal = 0.0;

  /* Normlisation scalar tcor_T**2 ... */
  if (normalised == true)
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      normVal += SQR(tcor_T[i][icor0]);
    }

    if (normVal == 0.0) return;
  }
  
  /* The T-3x1-array is used as follows ... */
  for (unsigned int ip = ip0; ip < p; ip++)
  {
    if ((int)(istep / iwait) > (int)(ip * pow((double)m,(double)is))) {

      icor = p * is + ip;

      if (tcor_T[0][icor] > -1.0E10 && tcor_T[1][icor] > -1.0E10 && tcor_T[2][icor] > -1.0E10) {
        
        NormT[icor]++;
      
        if (normalised == true) {
          if (normVal != 0.0) {
            for (unsigned int i = 0; i < 3; i++)
            {
              acf_T[i][icor] += tcor_T[i][icor] * tcor_T[i][icor0] / normVal;
            }
          } else { break; }


        } else {
          for (unsigned int i = 0; i < 3; i++)
          {
            acf_T[i][icor] += tcor_T[i][icor] * tcor_T[i][icor0];
          }
        }

      }
    } else { break; }

  }

}

/*--------------------------------------------------------------------*/
/* Calculate the correlation function + save... */
void Corrl::calcCorrl(FILE *acf_FilePtr) {

int icor, ip0;
double cor_time;

  /* Calculate the correlation function s.t. ... */
  for (unsigned int i = 0; i < s * p; i++)
  {
    if (NormT[i] != 0) {

      for (unsigned int j = 0; j < 3; j++)
      {
        acf_TA[i] += acf_T[j][i];
      }

      acf_TA[i] *= fact1 / (double)NormT[i];

    }
  }

  /* Write out variables for matlab to read */
  fprintf(acf_FilePtr, "%s %d %d %d\n", correlator.c_str(), s, p, m);

  /* Write out correlation function... */
  for (unsigned int is = 0; is < s; is++)
  {

    if (is == 0) ip0 = 0;
    else ip0 = (int)(p / m);

    for (unsigned int i = ip0; i < p; i++)
    {
      cor_time = (double)(i * iwait) * dt * pow((double)m,(double)is);

      icor = p * is + i;
      fprintf(acf_FilePtr, "%10f %10d %10f\n", cor_time, NormT[icor], acf_TA[icor]);
    }

  }

}

/*--------------------------------------------------------------------*/
/* Collate the data within one averaged- correlator to save to file */
void collateCorrlData(FILE *FilePtr, Corrl cor[], int iter) {

  static unsigned int tmpSize = correl_p * correl_s;
  int *tmpNorm1 = new int [tmpSize] ();
  double *tmpCorrl1 = new double [tmpSize] ();

  for (int i = 0; i < iter; i++)
  {
    cor[i].calcCorrl2(tmpNorm1, tmpCorrl1);
  }

  printCorrl(cor[0], FilePtr, poly.N, tmpNorm1, tmpCorrl1);

  delete[] tmpNorm1;
  delete[] tmpCorrl1;
}

/*--------------------------------------------------------------------*/
/* Calculate the correlation function */
void Corrl::calcCorrl2(int normArr[], double corrlArr[]) {

  /* Calculate the correlation function s.t. ... */
  for (unsigned int i = 0; i < s * p; i++)
  {
    if (NormT[i] != 0) {

      for (unsigned int j = 0; j < 3; j++)
      {
        acf_TA[i] += acf_T[j][i];
      }

      acf_TA[i] *= fact1 / (double)NormT[i];

      /* Transfer into the overall arrays. */
      normArr[i]  += NormT[i];
      corrlArr[i] += acf_TA[i];

    }
  }

}

/*--------------------------------------------------------------------*/
/* Calculate the correlation function - input the first cor each time */
void printCorrl(Corrl cor, FILE *acf_FilePtr, int Const, int normArr[],
                double corrlArr[]) {

int ip0, icor;
double cor_time;

  /* Normalise if available ... */
  if (Const != 0) {
    for (unsigned int i = 0; i < cor.s * cor.p; i++)
    {
      normArr[i]  /= (int)ceil((double)Const);
      corrlArr[i] /= (int)ceil((double)Const);
    }
  }

  /* Write out variables for matlab to read */
  fprintf(acf_FilePtr, "%s %d %d %d\n", cor.correlator.c_str(), cor.s, cor.p, cor.m);

  /* Write out correlation function... */
  for (unsigned int is = 0; is < cor.s; is++)
  {

    if (is == 0) ip0 = 0;
    else ip0 = (int)(cor.p / cor.m);
    
    for (unsigned int i = ip0; i < cor.p; i++)
    {
      cor_time = (double)(i * cor.iwait) * dt * pow((double)cor.m,(double)is);

      icor = cor.p * is + i;
      fprintf(acf_FilePtr, "%10f %10d %10f\n", cor_time, normArr[icor], corrlArr[icor]);
    }

  }

}