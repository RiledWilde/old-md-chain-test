#pragma once

#ifndef CORREL_H
#define CORREL_H

/*----------------------------------------------------------------*/
/* Correlation class for oeprator-tau method */
class Corrl {
   
   private:

      unsigned int s = 16, p = 36, m = 2;

      // Used to distinguish between each correlation function in Matlab.
      std::string correlator;

      int iwait;
      unsigned int *tcou;
      int *NormT;
      
      /* The normal tensor accumulation and calculation for operator tau: */
      double t_T[3], *acf_TA, *(tacc_T[3]), *(tcor_T[3]), *(acf_T[3]), fact1;


   public:

    /* Reintialise the tensor below ... */
    void reinit()
    {
      for (unsigned int i = 0; i < 3; i++)
      {
        t_T[i] = 0.0;
      }
    }

    Corrl() {

      reinit();

      fact1 = 1.0;
      iwait = 1;

      NormT  = new int [s*p] ();
      tcou   = new unsigned int [s*p] ();
      acf_TA = new double [s*p] ();

      for (unsigned int i = 0; i < 3; i++)
      {
        tacc_T[i] = new double [s] ();
        acf_T[i]  = new double [s * p] ();
        tcor_T[i] = new double [s * p];

        for (unsigned int j = 0; j < s * p; j++)
        {
          tcor_T[i][j] = -2.0e+10;
        }
      }

    }
    
    ~Corrl() {
      for (unsigned int i = 0; i < 3; i++)
      {
        delete[] acf_T[i];
        delete[] tcor_T[i];
        delete[] tacc_T[i];
      }

      delete[] acf_TA;
      delete[] NormT;
      delete[] tcou;
    }

    Corrl& operator=(const Corrl&) = default;
    Corrl(const Corrl&) = default;

    Corrl& operator=(Corrl&&) = default;
    Corrl(Corrl&&) = default;

    void constitute(int value, double prefactor, const char name[])
    {
      iwait = value;
      fact1 = prefactor;
      correlator = name;
    }

    int grab_wait() { return iwait; }

    void operator/=(double divisionIn)
    {
      for (unsigned int i = 0; i < 3; i++)
      {
        this->t_T[i] = this->t_T[i] / divisionIn;
      }
    }

    void operator*=(double multiplicationIn)
    {
      for (unsigned int i = 0; i < 3; i++)
      {
        this->t_T[i] = multiplicationIn * this->t_T[i];
      }
    }

    template <class TypeTmp> void transfer(TypeTmp cTmp)
    {
      t_T[0] = cTmp.out1();
      t_T[1] = cTmp.out2();
      t_T[2] = cTmp.out3();
    }

    friend std::ostream& operator<<(std::ostream& os, Corrl& cor)
    {
      os << "( " << cor.t_T[0] << ", " << cor.t_T[1] << ", " << cor.t_T[2] << " )";
      return os;
    }

    /* Push correlation right */
    void pushCorrl(unsigned int is, bool smooth);

    /* Calculate the correlation function - note that the normalisation is
        for looking at the relaxation analysis of the ACF. */
    void calc_ACF_Corrl(unsigned int is, unsigned int ip0, bool normalised);

    /* Calculate the correlation function */
    void calcCorrl(FILE *acf_FilePtr);
    void calcCorrl2(int normArr[], double corrlArr[]);
    friend void printCorrl(Corrl cor, FILE *acf_FilePtr, int Const, int normArr[], double corrlArr[]);

    /* Add to stress external function. */
    template <class TypeTmp1> 
    friend void addToStress(TypeTmp1 cor, double d[3], double spring_const);

    /* Friend functions to grab variables from it defined outside of this scope ... */
    friend void mainStressCorrl(Corrl &cor);
    friend void mainVelocityCorrl(Corrl &cor, int im1);
    friend void mainEndToEndCorrl(Corrl &cor, int t1, int t2);

    /* Operator Tau method calling is below. */
    void OperatorTau(bool normalised, bool smooth) {
    unsigned int ip0;

      for (unsigned int is = 0; is < s; is++)
      {
        if (is == 0) {
          pushCorrl(is, smooth);
          ip0 = 0;
          calc_ACF_Corrl(is, ip0, normalised);

        } else if (tcou[is] == m) {
          pushCorrl(is, smooth);
          ip0 = p / m;
          calc_ACF_Corrl(is, ip0, normalised);
        }
      }
    }

};

template<class T, class U = T> T exchange(T& obj, U&& new_value)
{
  T old_value = std::move(obj);
  obj = std::forward<U>(new_value);
  return old_value;
}

/* Stress elements - used alongside the correl method above in the wider code. */
class Sigma { 

   private:

     double t1, t2, t3;

   public:
      
    Sigma(double xval=0.0, double yval=0.0, double zval=0.0)
        : t1(xval), t2(yval), t3(zval) { }
    ~Sigma() = default;

    Sigma& operator=(const Sigma& rhs)
    {
      if (this == &rhs) return *this;
      t1 = rhs.t1; t2 = rhs.t2; t3 = rhs.t3;
      return *this;
    }

    Sigma(const Sigma& src) : t1(src.t1), t2(src.t2), t3(src.t3) { }

    Sigma& operator=(Sigma&& rhs)
    {
      if (this == &rhs) return *this;
      t1 = rhs.t1; t2 = rhs.t2; t3 = rhs.t3;
      rhs.t1 = 0.0; rhs.t2 = 0.0; rhs.t3 = 0.0;
      return *this;
    }

    Sigma(Sigma&& src)
    {
      t1 = exchange(src.t1,0.0);
      t2 = exchange(src.t2,0.0);
      t3 = exchange(src.t3,0.0);
    }

    void reinit() { t1 = 0.0; t2 = 0.0; t3 = 0.0; }

    void operator /= (double divisionIn)
    {
      this->t1 = this->t1 / divisionIn;
      this->t2 = this->t2 / divisionIn;
      this->t3 = this->t3 / divisionIn;
    }

    void operator *= (double multiplicationIn)
    {
      this->t1 = multiplicationIn * this->t1;
      this->t2 = multiplicationIn * this->t2;
      this->t3 = multiplicationIn * this->t3;
    }

    /* Overload the << operator for quicker print out */
    friend std::ostream& operator << (std::ostream& os, Sigma& c)
    {
      os << "( " << c.t1 << ", " << c.t2 << ", " << c.t3 << " )";
      return os;
    }

    void addIn1(double valIn) { this->t1 += valIn; }
    void addIn2(double valIn) { this->t2 += valIn; }
    void addIn3(double valIn) { this->t3 += valIn; }

    double out1() { return this->t1; }
    double out2() { return this->t2; }
    double out3() { return this->t3; }

};
extern Sigma sT

#endif