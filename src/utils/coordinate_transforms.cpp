#include "coordinate_transforms.h"

void globalToVehicle(vector<double>* global_coords_x,
                     vector<double>* global_coords_y,
                     double ref_x,
                     double ref_y,
                     double ref_yaw) {
  for(size_t i = 0; i < global_coords_x->size(); i++) {
    double shift_x = (*global_coords_x)[i] - ref_x;
    double shift_y = (*global_coords_y)[i] - ref_y;

    (*global_coords_x)[i] = shift_x * cos(ref_yaw) + shift_y * sin(ref_yaw);
    (*global_coords_y)[i] = -shift_x * sin(ref_yaw) + shift_y * cos(ref_yaw);
  }
}

void vehicleToGlobal(vector<double>* vehicle_coords_x,
                     vector<double> *vehicle_coords_y,
                     double ref_x,
                     double ref_y,
                     double ref_yaw) {
  for(size_t i = 0; i < vehicle_coords_x->size(); i++) {
    double local_x = (*vehicle_coords_x)[i];
    double local_y = (*vehicle_coords_y)[i];
    (*vehicle_coords_x)[i] = local_x * cos(ref_yaw) - local_y * sin(ref_yaw) + ref_x;
    (*vehicle_coords_y)[i] = local_x * sin(ref_yaw) + local_y * cos(ref_yaw) + ref_y;
  }
}

