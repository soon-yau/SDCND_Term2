#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <ctime>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"
#include "matplotlibcpp.h"
// for convenience
using json = nlohmann::json;
using namespace std;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }
double rad2steer(double x){ return (x*180/(M_PI*25));}
// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

vector<double> map2carCoord(double ptsx, double ptsy, double px, double py, double phi)
{
  double theta = atan2(ptsy-py,ptsx-px)-phi;
  double l = sqrt(pow(ptsx-px,2)+pow(ptsy-py,2));
  double x = l*cos(theta);
  double y = l*sin(theta);
  return{x,y};
}


tuple<vector<vector<double>>, Eigen::VectorXd>
 SolveWrapper(MPC &mpc, vector<double> &ptsx, vector<double> &ptsy,
                   double px, double py, double psi, double v)
{
  vector<double> next_x_vals;
  vector<double> next_y_vals;
  
  Eigen::VectorXd eigen_ptsx(ptsx.size());
  Eigen::VectorXd eigen_ptsy(ptsy.size());
  for (int i=0; i<ptsx.size(); i++)
  {
    auto coords=map2carCoord(ptsx[i], ptsy[i], px, py, psi);            
    next_x_vals.push_back(coords[0]);
    next_y_vals.push_back(coords[1]);

    eigen_ptsx(i)=coords[0];
    eigen_ptsy(i)=coords[1];
  }

  // vehicle coord
  auto coeffs = polyfit(eigen_ptsx, eigen_ptsy, 2);
  double cte = polyeval(coeffs, 0);
  double epsi = -atan(coeffs[1]);
  //cout<<"cte:"<<cte<<" Coef:"<<coeffs[0]<<','<<coeffs[1]<<','<<coeffs[2]<<','<<coeffs[2]<<endl;
  Eigen::VectorXd state(6);                    
  state << 0, 0, 0, v, cte, epsi;
  
  
  auto vars = mpc.Solve(state, coeffs);

  double steer_value=(vars[0][0]);
  double throttle_value=vars[0][1];

  return(make_tuple(vector<vector<double>>{vars[0], next_x_vals, next_y_vals, vars[1], vars[2]}, coeffs));
}


void test(MPC &mpc)
{
/*  
  //42["telemetry",{"ptsx":[91.19354,107.3083,114.2635,122.5583,126.6183,129.1083],"ptsy":[-139.979,-134.209,-130.6506,-123.779,-118.159,-108.669],"psi_unity":1.290459,"psi":0.2803371,"x":91.38155,"y":-140.1407,"steering_angle":-0.02231064,"throttle":0.250027,"speed":29.16532}]
  //cte:0.20909 Coef:0.20909,0.327313,-0.024286,-0.024286

42["telemetry",{"ptsx":[107.3083,114.2635,122.5583,126.6183,129.1083,129.1283],"ptsy":[-134.209,-130.6506,-123.779,-118.159,-108.669,-100.349],"psi_unity":1.068266,"psi":0.50253,"x":111.7632,"y":-131.0737,"steering_angle":0.02676496,"throttle":0.2500176,"speed":29.11715}]
delta t:0.311836
cte:-0.554005 Coef:-0.554005,-0.00295137,-0.00126532,-0.00126532  
*/
  double ptsx_samples[]={91.19354,107.3083,114.2635,122.5583,126.6183,129.1083};
  double ptsy_samples[]={-139.979,-134.209,-130.6506,-123.779,-118.159,-108.669};
  vector<double>ptsx(ptsx_samples,ptsx_samples+sizeof(ptsx_samples)/sizeof(double));
  vector<double>ptsy(ptsy_samples,ptsy_samples+sizeof(ptsy_samples)/sizeof(double));

  double px=91.38155;
  double py=-140.1407;
  double psi=0.2803371;
  double v=29.16532;
  

  auto t=SolveWrapper(mpc, ptsx, ptsy, px, py, psi, v);
  vector<vector<double>> ret=get<0>(t);
  double steer_value = ret[0][0];
  double throttle_value = ret[0][1];
  vector<double>next_x_vals=ret[1];
  vector<double>next_y_vals=ret[2];          
  vector<double>mpc_x_vals=ret[3];
  vector<double>mpc_y_vals=ret[4];
  auto coeffs=get<1>(t);

  vector<double>fitted_ptsy;
  for (auto x:next_x_vals)
  {
    fitted_ptsy.push_back(polyeval(coeffs, x));
  }

  cout<<"Test complete"<<endl;
  cout<<fitted_ptsy.size()<<" "<<next_x_vals.size();


  matplotlibcpp::subplot(1, 1, 1);
  matplotlibcpp::title("Trajectory");
  matplotlibcpp::plot(next_x_vals,next_y_vals, "g-");
  matplotlibcpp::plot(next_x_vals, fitted_ptsy, "b"); 
  //matplotlibcpp::subplot(2, 1, 2);
  
  //matplotlibcpp::title("Delta (Radians)");
  matplotlibcpp::plot(mpc_x_vals, mpc_y_vals, "r*"); 
  matplotlibcpp::show();

}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;
  //test(mpc);  

  float solver_latency=0;
  h.onMessage([&mpc, &solver_latency](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double a = j[1]["throttle"];
          double delta = j[1]["steering_angle"];
          
          
          /*
          * TODO: Calculate steeering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
       
          // predict new states 100ms into future plus solver's latency
          // solver latency varies with machine, you can adjust it accordingly, it is in the cout
          double dt=0.1+0.03;
          clock_t prevT=clock();
          v = v * 0.44704;  //convert miles per hour to meter per second
          double px_ = px + v * cos(psi)*dt;
          double py_ = py + v * sin(psi)*dt;
          double psi_ = psi - v*delta*dt/Lf;
          double v_ = v + a*dt;
          //cout<<"delta x:"<<px_-px<<" delta y:"<<py_-py<<" delta psi:"<< psi_-psi<<" delta v:"<<a*dt<<endl;
          
          auto t=SolveWrapper(mpc, ptsx, ptsy, px_, py_, psi_, v_);
        
          vector<vector<double>> ret=get<0>(t);
          vector<double>next_x_vals=ret[1];
          vector<double>next_y_vals=ret[2];          
          vector<double>mpc_x_vals=ret[3];
          vector<double>mpc_y_vals=ret[4];
          bool success=mpc_x_vals.size()!=0;
          double steer_value = success?ret[0][0]:delta;
          double throttle_value = success?ret[0][1]:a;

          json msgJson;
          msgJson["steering_angle"] = -rad2steer(steer_value);
          msgJson["throttle"] = throttle_value;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
//          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          solver_latency=float(clock()-prevT)/CLOCKS_PER_SEC;
          cout<<"solver_latency:"<<solver_latency<<endl;          
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

