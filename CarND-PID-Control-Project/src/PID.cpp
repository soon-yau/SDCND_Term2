#include "PID.h"
#include <cmath>
#include <iostream>
#include <iomanip>   
/*
* TODO: Complete the PID class.
*/
PID::PID() {}

PID::~PID() {}

void PID::Reset()
{
  bestCte=100000;
  intCte=0;
  prevCte=0;  

  runTwiddle=true;
  twiddle_i=0;
  iteration=0;

  delta[PID_P]=0.05;
  delta[PID_I]=0.001;
  delta[PID_D]=0.05;
}

void PID::Init(double Kp, double Ki, double Kd) {
  Reset();

  // add delta to P coef as part of twaddle algorithm initialization
  coef[PID_P]=Kp+delta[PID_P];
  coef[PID_I]=Ki;
  coef[PID_D]=Kd;

}

double PID::TotalError() {
  return (fabs(delta[0])+fabs(delta[1])+fabs(delta[2]));
}

double PID::CalcSteerValue(double cte)
{
    double diffCte=cte-prevCte;
    prevCte=cte;
    intCte+=cte;

    double steer_value=-coef[PID_P]*cte-coef[PID_D]*diffCte-coef[PID_I]*intCte;
    std::cout<<"coef:"<<coef[PID_P]<<"," <<coef[PID_I]<<","<<coef[PID_D] <<std::endl
             <<"cte:"<<cte<<" diffCte"<<diffCte<<" intCte:"<<intCte
              <<" steer:"<<steer_value<<std::endl;

    return steer_value;
}
/*
Twiddle is implemented in a simple state machine where twiddle_i
is the index of PID coefficient, i.e. twiddle_i=0 refers to P, 
twiddle_i=1 is I and twiddle_i=2 is D. The state changes based on twiddle_i, 
its previous value and cte.
*/
double PID::Twiddle(double cte)
{
  static int  prev_twiddle_i=-1;
  // iteration is used to prevent calibration stopping too early
  // e.g. when cte is small in straight track
  if (TotalError()<0.0001 && iteration>1000)
  {
    std::cout<<"Total Error:"<<TotalError()<<std::endl;
    std::cout<<delta[0]<<','<<delta[1]<<','<<delta[2]<<','<<std::endl;
    runTwiddle=false;
    return 0;
  }
  if (init==true)
  {
    init=false;
    return 0;
  }
  int twiddle_state;
  
  if (cte<bestCte)
  {
    twiddle_state=0;    
    bestCte=1.1*cte;  // prevent the learning to stop too soon
    delta[twiddle_i]*=1.1;

    // change state
    prev_twiddle_i=twiddle_i;
    twiddle_i=(twiddle_i+1)%3;
    coef[twiddle_i]+=delta[twiddle_i];
  }
  else if (prev_twiddle_i!=twiddle_i)
  {
    twiddle_state=1;    
    coef[twiddle_i]-=2*delta[twiddle_i];
    prev_twiddle_i=twiddle_i;
  }
  else
  {
    twiddle_state=2;    
    coef[twiddle_i]+=delta[twiddle_i];
    delta[twiddle_i]*=0.9;

  // change state
    prev_twiddle_i=twiddle_i;
    twiddle_i=(twiddle_i+1)%3;
    coef[twiddle_i]+=delta[twiddle_i];
  }
  iteration++;
  std::cout<<"iter:"<<iteration<<" i:"<<twiddle_i<<" state:"<<twiddle_state<<" delta:"<<delta[0]<<","<<delta[1]<<","<<delta[2]<<std::endl;
}