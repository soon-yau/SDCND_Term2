#ifndef PID_H
#define PID_H

enum PidEnum {PID_P,PID_I,PID_D};

class PID {
public:
  // 
  int twiddle_i=0;
  int iteration=0;
  bool runTwiddle=true;
  bool init=true;
  double intCte=0;
  double prevCte=0;
  double bestCte;  
  /*
  * Errors
  */
  //double p_error;
  //double i_error;
  //double d_error;
  double delta[3];

  /*
  * Coefficients
  */ 
  //double Kp;
  //double Ki;
  //double Kd;
  double coef[3];

  /*
  * Constructor
  */
  PID();

  /*
  * Destructor.
  */
  virtual ~PID();


  // Initialize PID.
  void Init(double Kp, double Ki, double Kd);

  // reset PID states and delta but keeping the PID coefficients
  void Reset();

  // calculate the sum of absolute value of each coefficient delta
  double TotalError();

  // run twiddle algorithm iteration to update the delta and coefficient
  // but does not calculate the steering value
  double Twiddle(double cte);

  // calculate the steering value from PID coefficients and cte values
  double CalcSteerValue(double cte);

};

#endif /* PID_H */
