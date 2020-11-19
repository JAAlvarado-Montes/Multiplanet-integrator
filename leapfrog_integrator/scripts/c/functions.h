double initial_velocity(double M, double m, double r0, double a){
  double v0;
  v0 = pow((M + m) * (2.0 / r0 - 1.0 / a), 0.5);

  return v0;
}

double mag_vec_soft(double a, double b, double soft){
  double mag;
  mag = pow(a*a + b*b + soft*soft, 0.5);

  return mag;
}


double mag_vec(double a, double b){
  double mag;
  mag = pow(a*a + b*b, 0.5);

  return mag;
}

double *grav_force(double m, double M, double sep_vec[1][2], double pos){
  static double g_force[1][2];

  g_force[0][0] = - m * M * sep_vec[0][0] / pow(pos, 3);
  g_force[0][1] = - m * M * sep_vec[0][1] / pow(pos, 3);

  return *g_force;
}

double *grav_relat_force(double m, double M, double ls, double sep_vec[1][2], double pos, double vel_vec[1][2]){
  static double grav_rel_f[1][2];
  double vel_squared = vel_vec[0][0] * vel_vec[0][0] + vel_vec[0][1] * vel_vec[0][1];
  double dot_prod = sep_vec[0][0] * vel_vec[0][0] + sep_vec[0][1] * vel_vec[0][1];

  grav_rel_f[0][0] = m * M / (pow(ls, 2) * pow(pos, 3)) * ((4. * M / pos - vel_squared) * sep_vec[0][0] + 4. * dot_prod * vel_vec[0][0]);
  grav_rel_f[0][1] = m * M / (pow(ls, 2) * pow(pos, 3)) * ((4. * M / pos - vel_squared) * sep_vec[0][1] + 4. * dot_prod * vel_vec[0][1]);

  return *grav_rel_f;
}

double timestep(double forces[1][2], double M, double R, double dt_std){
  double dt;

  dt = 0.3 * dt_std * pow(R / (mag_vec(forces[0][0], forces[0][1]) / M), 0.5);

  return dt;
}

double find_min(double array[], int array_size){
  int c;
 	double minimum;
 
  minimum = array[0];
 
  for(c=1; c < array_size; c++){
    if(array[c] < minimum){
      minimum = array[c];
      }
    } 
 
  return minimum;
}
