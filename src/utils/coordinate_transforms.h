#ifndef COORDINATE_TRANSFORM_H
#define COORDINATE_TRANSFORM_H
#include <vector>
#include <math.h>

using namespace std;
// CLOCKWISE rotation:
// https://en.wikipedia.org/wiki/Rotation_matrix
void globalToVehicle(vector<double>* global_coords_x, vector<double> *global_coords_y, double ref_x, double ref_y, double ref_yaw);
// inverse to globalToVehicle
void vehicleToGlobal(vector<double>* vehicle_coords_x, vector<double> *vehicle_coords_y, double ref_x, double ref_y, double ref_yaw);

#endif /* COORDINATE_TRANSFORM_H */
