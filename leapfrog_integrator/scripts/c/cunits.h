double *get_canon_units(double a, double b, char var[]){

  char var1[] = "uM";
  char var2[] = "uL";
  char var3[] = "uT";
  double G = 6.6740831e-11;

  double uM, uL, uT;

  static double array[3];

  if(strcmp(var, var1) == 0){
    // I will find the unit of mass
    uL = a;
    uT = b;
    uM = uL * uL *uL/(uT * uT * G);
    }	

  if(strcmp(var, var2) == 0){
    // I will find the unit of longitude
    uM = a;
    uT = b;
    uL = pow(G * uM * uT * uT, 1.0 / 3.0);
    }

  if(strcmp(var, var3) == 0){
    // I will find the unit of time
    uM = a;
    uL = b;
    uT = pow(uL * uL * uL / (G * uM), 0.5);
    }

  array[0] = uM;
  array[1] = uL;
  array[2] = uT;
    
  return array;
}