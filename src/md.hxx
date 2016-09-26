#ifndef HXX_MD
#define HXX_MD

#include <iostream>
#include <fstream>
#include <string>

#include <cassert>
#include <cmath>

#include <sys/time.h>

using namespace std;

typedef float DorF;
//typedef double DorF;

#define PBC
//#define NEWTON3

class MD{
protected:
  DorF *r,*rx,*ry,*rz;
  DorF *v,*vx,*vy,*vz;
  DorF *f,*fx,*fy,*fz;

  int  nmol;
  DorF dens;
  DorF temp;

  DorF dt,dth;
  DorF rcut;
  DorF cs,csh;
  DorF pot,kin,tot;

public:
  MD();
  ~MD();

  class Timer{
  private:
    unsigned int   count;
    double elapsed_time,max,min;
    struct timeval timer_start,timer_end;
  public:
    Timer(){}
    ~Timer(){}
    void flush(){elapsed_time=max=0.0;min=1e32;count=0;}
    void start(){
      gettimeofday(&timer_start,NULL);
    }
    void end(){
      gettimeofday(&timer_end,NULL);
      double e = (timer_end.tv_sec - timer_start.tv_sec) + (timer_end.tv_usec - timer_start.tv_usec) * 0.000001;
      elapsed_time += e;
      if(max<e) max = e;
      if(min>e) min = e;
      count++;
    }
    void print(ostream &os){
      os << "# elapsed time[sec]:\t" << elapsed_time << endl;
      os << "# average time[sec]:\t" << elapsed_time / (double)count << endl;
      os << "# maximum time[sec]:\t" << max << endl;
      os << "# minimum time[sec]:\t" << min << endl;
    }
  };
  Timer timer;

  virtual void Initialize(const int,const DorF,const DorF);
  void SetCoorFCC();
  void SetVelRandom();

  void IntegrateCoor();
  void IntegrateVel();
  virtual void CalcForce();

  void KillMomentum();
  void VelocityScaling();

  void CalcKineticEnergy();
  void CalcPotentialEnergy();
  void CalcHamiltonian(){tot = pot + kin;};

  void DisplayEnergies(ostream&);
  void DisplayConditions(ostream&);

  void OutputCDV(const string);
};

#endif
