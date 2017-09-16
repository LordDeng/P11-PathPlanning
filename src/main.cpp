#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "utils/spline.h"
#include "utils/coordinate_transforms.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// returns the id (in sensor_fusion) of the occuping vehicle
// -1 if no vehicle
int isForwardLaneOccupied(const vector<vector<double>>& sensor_fusion,
                           double car_s,
                           int lane) {
  double lane_center = 2.0 + (double)lane * 4.0;
  for(size_t i = 0; i < sensor_fusion.size(); i++) {
    float d = sensor_fusion[i][6];
    if(d < (lane_center + 2.0) && d > (lane_center - 2)) {
      double vx = sensor_fusion[i][3];
      double vy = sensor_fusion[i][4];
      double v_mag = sqrt(vx * vx + vy * vy);
      double check_car_s = sensor_fusion[i][5];

      check_car_s += 0.1 * v_mag;
      if(check_car_s > car_s && (check_car_s - car_s) < 18) {
        return (int)sensor_fusion[i][0];
      }
    }
  }
  return -1;
}

// returns the id (in sensor_fusion) of the occuping vehicle
// -1 if no vehicle
int isBackwardLaneOccupied(const vector<vector<double>>& sensor_fusion,
                           double car_s,
                           int lane) {
  double lane_center = 2.0 + (double)lane * 4.0;
  for(size_t i = 0; i < sensor_fusion.size(); i++) {
    float d = sensor_fusion[i][6];
    if(d < (lane_center + 2.0) && d > (lane_center - 2)) {
      double vx = sensor_fusion[i][3];
      double vy = sensor_fusion[i][4];
      double v_mag = sqrt(vx * vx + vy * vy);
      double check_car_s = sensor_fusion[i][5];

      check_car_s += 0.1 * v_mag;
      if(check_car_s < car_s && (car_s - check_car_s) < 18) {
        return (int)sensor_fusion[i][0];
      }
    }
  }
  return -1;
}

double getSpeedById(const vector<vector<double>>& sensor_fusion,
                    int car_id) {
  for(size_t i = 0; i < sensor_fusion.size(); i++) {
    if(sensor_fusion[i][0] == car_id) {
      double vx = sensor_fusion[i][3];
      double vy = sensor_fusion[i][4];
      double v_mag = sqrt(vx * vx + vy * vy);
      return v_mag;
    }
  }
  return 1000.0;
}

enum class State {
  Ready,
  KeepLane,
  PrepareChangeLeft,
  PrepareChangeRight,
  ChangeLaneLeft,
  ChangeLaneRight,
};

std::ostream& operator<<(std::ostream& os, const State& state)
{
  switch(state) {
    case State::Ready:
      os << "Ready             ";
      break;
    case State::KeepLane:
      os << "KeepLane          ";
      break;
    case State::PrepareChangeLeft:
      os << "PrepareChangeLeft ";
      break;
    case State::PrepareChangeRight:
      os << "PrepareChangeRight";
      break;
    case State::ChangeLaneLeft:
      os << "ChangeLaneLeft    ";
      break;
    case State::ChangeLaneRight:
      os << "ChangeLaneRight   ";
      break;
  }
  return os;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  int lane = 1;
  double ref_vel = 0.0; // km/h
  State state = State::Ready;
  h.onMessage([&map_waypoints_x,&map_waypoints_y,
               &map_waypoints_s,&map_waypoints_dx,
               &map_waypoints_dy,
               &lane,
               &ref_vel,
               &state](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
            // format: [car_id, x, y, vx, vy, s, d]
            vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];
            int prev_size = previous_path_x.size();
            double ego_lane_center = 2.0 + (double)lane * 4.0;
            if(prev_size > 0) {
              car_s = end_path_s;
            }

            int lane0_next_car_id = isForwardLaneOccupied(sensor_fusion, car_s, 0);
            int lane1_next_car_id = isForwardLaneOccupied(sensor_fusion, car_s, 1);
            int lane2_next_car_id = isForwardLaneOccupied(sensor_fusion, car_s, 2);

            bool lane0_occupied = lane0_next_car_id != -1;
            bool lane1_occupied = lane1_next_car_id != -1;
            bool lane2_occupied = lane2_next_car_id != -1;

            bool decelerate;
            if(lane == 0) {
              decelerate = lane0_occupied;
            } else if(lane == 1) {
              decelerate = lane1_occupied;
            } else {
              decelerate = lane2_occupied;
            }

            State previous_state = state;
            if(state == State::Ready) {
              state = State::KeepLane;
            } // State::Ready
            else if(state == State::KeepLane) {
              if(decelerate) {
                if(lane == 0) {
                  if(lane1_occupied) {
                    state = State::KeepLane;
                  } else {
                    state = State::PrepareChangeRight;
                  }
                }
                else if(lane == 1) {
                  if(lane0_occupied && lane2_occupied) {
                    state = State::KeepLane;
                  }
                  else if(lane0_occupied && !lane2_occupied) {
                    state = State::PrepareChangeRight;
                  }
                  else if(!lane0_occupied && lane2_occupied) {
                    state = State::PrepareChangeLeft;
                  }
                  else if(!lane0_occupied && !lane2_occupied) {
                    // resolve by going to left lane
                    state = State::PrepareChangeLeft;
                  }
                }
                else if(lane == 2) {
                  if(lane1_occupied) {
                    state = State::KeepLane;
                  } else {
                    state = State::PrepareChangeLeft;
                  }
                }
              } // decelerate
              else {
              } // !decelerate
            } // State::KeepLane
            else if(state == State::PrepareChangeLeft) {
              if(lane == 1) {
                if(isBackwardLaneOccupied(sensor_fusion, car_s, 0) == -1) {
                  state = State::ChangeLaneLeft;
                } else {
                  state = State::PrepareChangeLeft;
                }
              }
              else if(lane == 2) {
                if(isBackwardLaneOccupied(sensor_fusion, car_s, 1) == -1) {
                  state = State::ChangeLaneLeft;
                } else {
                  state = State::PrepareChangeLeft;
                }
              }
            } // State::PrepareChangeLeft
            else if(state == State::PrepareChangeRight) {
              if(lane == 0) {
                if(isBackwardLaneOccupied(sensor_fusion, car_s, 1) == -1) {
                  state = State::ChangeLaneRight;
                } else {
                  state = State::PrepareChangeRight;
                }
              }
              else if(lane == 1) {
                if(isBackwardLaneOccupied(sensor_fusion, car_s, 2) == -1) {
                  state = State::ChangeLaneRight;
                } else {
                  state = State::PrepareChangeRight;
                }
              }
            } // State::PrepareChangeRight
            else if(state == State::ChangeLaneLeft) {
              lane -= 1;
              state = State::KeepLane;
            } // State::ChangeLaneLeft
            else if(state == State::ChangeLaneRight) {
              lane += 1;
              state = State::KeepLane;
            } // State::ChangeLaneRight

            if(state != previous_state) {
              cout << previous_state << " --> " << state << endl;
            }
            if(decelerate) {
              ref_vel -= 0.6;
            }
            else {
              double max_vel = 79.5;
              double next_vehicle_speed;
              if(lane == 0) {
                next_vehicle_speed = 1.61 * getSpeedById(sensor_fusion, lane0_next_car_id);
              }
              else if(lane == 1) {
                next_vehicle_speed = 1.61 * getSpeedById(sensor_fusion, lane0_next_car_id);
              }
              if(lane == 2) {
                next_vehicle_speed = 1.61 * getSpeedById(sensor_fusion, lane0_next_car_id);
              }
              max_vel = min(max_vel, next_vehicle_speed);
              if(ref_vel < max_vel){
                ref_vel += 0.6;
              }
            }

            vector<double> anchor_x;
            vector<double> anchor_y;

            double ref_x = car_x;
            double ref_y = car_y;
            // double ref_yaw = deg2rad(car_yaw);
            double ref_yaw = car_yaw;

            if(prev_size < 2) {
              // create points tangent to the car
              double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              anchor_x.push_back(prev_car_x);
              anchor_x.push_back(car_x);

              anchor_y.push_back(prev_car_y);
              anchor_y.push_back(car_y);

            }
            else {
              // use previous path points as reference
              // using the current car_x, car_y, car_yaw would not be as smooth
              ref_x = previous_path_x[prev_size - 1];
              ref_y = previous_path_y[prev_size - 1];

              double prev_car_x = previous_path_x[prev_size - 2];
              double prev_car_y = previous_path_y[prev_size - 2];
              ref_yaw = atan2(ref_y - prev_car_y, ref_x - prev_car_x);

             anchor_x.push_back(prev_car_x);
             anchor_x.push_back(ref_x);

             anchor_y.push_back(prev_car_y);
             anchor_y.push_back(ref_y);
            }

            // create evenly spaced points in frenet coords - every 30m
            for(int i = 1; i < 4; i++) {
              vector<double> anchor_xy = getXY(car_s + 30.0 * i, ego_lane_center, map_waypoints_s, map_waypoints_x, map_waypoints_y);
              anchor_x.push_back(anchor_xy[0]);
              anchor_y.push_back(anchor_xy[1]);
            }
            globalToVehicle(&anchor_x, &anchor_y, ref_x, ref_y, ref_yaw);

            tk::spline s;
            s.set_points(anchor_x, anchor_y);

            // create points from spline in vehicle coordinates
            vector<double> spline_x;
            vector<double> spline_y;
            double target_x = 30.0;
            double target_y = s(target_x);
            double target_dist = sqrt(pow(target_x, 2) + pow(target_y, 2));
            double x_point = 0.0;
            for(int i = 0; i < 50 - prev_size; i++) {
              double N = target_dist / (0.02 * ref_vel / 3.6);  // kmh -> m/s
              x_point += target_x / N;
              double y_point = s(x_point);

              spline_x.push_back(x_point);
              spline_y.push_back(y_point);
            }

            vehicleToGlobal(&spline_x, &spline_y, ref_x, ref_y, ref_yaw);
            vector<double> next_x_vals;
            vector<double> next_y_vals;

            // add previous points
            for(int i = 0; i < prev_size; i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            // add points from spline
            for(int i = 0; i < 50 - prev_size; i++) {
              next_x_vals.push_back(spline_x[i]);
              next_y_vals.push_back(spline_y[i]);
            }

            json msgJson;
            msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
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
















































































