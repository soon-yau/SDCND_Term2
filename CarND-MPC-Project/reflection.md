# Reflection on my MPC C++ implementation

## The Model
I use 6 states as described in lessons:
* px - car's position in x-axis
* py - car's position in y-axis
* phi - car's angle
* v - car's velocity in mph
* cte - cross track error
* epsi - error of phi to reference
The actuators are
* steering, the unit is in radian from the optimizer
* throttle

After receiving the path trajectory (ptsx and ptsy), they are first converted from map coordinates to car coordinate. I derived the conversion equations myself. 
They are fitted to 3 degree of polynomial functions and the coefficents are obtained.
Since the path is in car coordinate, therefore px, py and phi becomes zero. We find the CTE by finding the error of the car position (px,py) to the supposed trajectory. This is done by using x=0 to the polynimal equations to find y, and CTE is y-py.
Similarly, epsi is the error between the car's angle (which is zero) and the gradient of path at x=0. 

The update equations were taken from the lessons and MPC quiz solution, nothing fancy there. The multiplication factors to the cost e.g. COST_CTE 1000.0  were obtained by trial and error.

Finally, the steering obtained from the optimizer is normalized to plus minus 1.0 that corresponds to plus minus 25 degrees and multiplied with -1 sign. This is because in car coordinate system, x-axis is toward the front of the car, and y-axis is to the left. Steering of positive angle would mean steering to the left which is a negative value assumed by the simulator.

## Timestep Length and Elapsed Duration
This is probably the toughest path, I spent a lot of time understanding the effect of using different N and dt. There is a latency of 100ms (0.1s) for the simulator to react after receiving the steering and throttle commands. 
Therefore, if dt is smaller than the latency, that means the first actuator values from optimizer were out of date when the car responded. One potential solution is to factor in this latency and take the later actuators value from optimizer instead of the first one.
However, this won't work in practise. If dt is small, then N will need to be large for the model to accurately predict the path ahead and this will result in higher computational time. I tried dt=0.05 and N=100 (to plan for 5s ahead), with optimizer won't be able to solve it within specified time of 0.5s. 
Beside, it is tedious to have to work out which actuators predictions to take. 

It is therefore better to choose dt that is larger than the latency, note that actual latency is the 0.1s plus time taken for optimizer to solve the optimization problem.
I tried N=25, dt=0.5 (to predict 12.5 ahead) but found the optimizer will need to take more than 1.0s to handle one of the road curve. The latency is too high and the car failed to respond to the fast changing landscape and crashed.
After trying different combinations, I found N=10 and dt=0.2 to has low optimizer computation cost where dt is close to total latency and 5.0s of planning ahead is just about right.

## Polynomial Fitting
I tried different degree of polynomial by plotting and comparing the waypoints and the fitted curve. We need at least 2 degree (quadratic) to approximate the road curve and I ended up using 3 degrees. The equations for CTE calculation FG_eval() has been updated to reflect that:

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + 
                      coeffs[2] * pow(x0,2) +
                      coeffs[3] * pow(x0,3);  // fits 3 degree of polynomial
                      
      AD<double> psides0 = CppAD::atan(coeffs[1]+
                                      2*x0*coeffs[2]+
                                      3*x0*x0*coeffs[3]); // find gradient of 3 degree polynomial function at x0

## Notes on Other Implementation
* main.cpp:SolveWrapper() to take in states and waypoitns from simulator that is to passed to MPC solver
* I took waypoints and other simulator states of challenging curve in main.cpp:test() to try different MPC hyperparameters, look at optimizer solver's status and plot the results for visualizaton. There are similar plottings in MPC:Solve() which has been commented out.





