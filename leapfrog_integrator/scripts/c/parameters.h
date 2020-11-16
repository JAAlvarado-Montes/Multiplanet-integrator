typedef struct Inpar_args{

  // inner planet properties
  double m1;
  double R1;
  double gyr_rad_inner;
  double a_inner;
  double e_inner;
  double P_rot_inner;
  double alpha_inner;

  // outer planet properties
  double m2;
  double R2;
  double gyr_rad_outer;
  double a_outer;
  double e_outer;
  double P_rot_outer;
  double alpha_outer;

  // Inner B properties
  double m0;
  double R0;

  // General properties
  const char *sim_name;
  double total_t;
  double dt_store;
  double dt_std;
  double softening;

  // Initial conditions
  
  // Canonical units info.
  double uM;
  double uL;
  double uT;
} Inpar;

Inpar args(){
  
  Inpar rest; //defining the REturnSTucture

  config_t cfg;
  config_setting_t *setting;


  config_init(&cfg);
  
  // Read the file. If there is an error, report it and exit. 
  if(! config_read_file(&cfg, "config.cfg"))
    {
      fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
	      config_error_line(&cfg), config_error_text(&cfg));
      config_destroy(&cfg);
    }

  ////////////////////////////////////////////////////////////
  //
  // Getting and storing General Properties and initial conditions
  //
  ////////////////////////////////////////////////////////////

  const char *sim_name;
  double total_t;
  double dt_store;
  double dt_std;
  double softening;
  
  
  config_lookup_string(&cfg, "sim_name", &sim_name);
  config_lookup_float(&cfg,  "total_t", &total_t);
  config_lookup_float(&cfg,  "dt_store", &dt_store);
  config_lookup_float(&cfg,  "dt_std", &dt_std);
  config_lookup_float(&cfg,  "softening", &softening);

  ////////////////////////////////////////////////////////////
  //
  // Getting information about inner planet
  //
  ////////////////////////////////////////////////////////////
  
  double m1;
  double R1;
  double gyr_rad_inner;
  double a_inner;
  double e_inner;
  double P_rot_inner;
  double alpha_inner;
  
  
  setting = config_lookup(&cfg, "bodies.inner_planet");
  if(setting != NULL){
    int count = config_setting_length(setting);
    int i;
    
    for(i = 0; i < count; ++i){
      config_setting_t *object = config_setting_get_elem(setting, i);
      
      config_setting_lookup_float(object,  "mass",    &m1);
      config_setting_lookup_float(object,  "radius",  &R1);
      config_setting_lookup_float(object,  "gyr_rad", &gyr_rad_inner);
      config_setting_lookup_float(object,  "a",       &a_inner);
      config_setting_lookup_float(object,  "e",       &e_inner);
      config_setting_lookup_float(object,  "P_rot",   &P_rot_inner);
      config_setting_lookup_float(object,  "alpha",   &alpha_inner);
    }
  }
  

  ////////////////////////////////////////////////////////////
  //
  // Getting information about outer planet
  //
  ////////////////////////////////////////////////////////////
  
  double m2;
  double R2;
  double gyr_rad_outer;
  double a_outer;
  double e_outer;
  double P_rot_outer;
  double alpha_outer;
  
  
  setting = config_lookup(&cfg, "bodies.outer_planet");
  if(setting != NULL){
    int count = config_setting_length(setting);
    int i;
    
    for(i = 0; i < count; ++i){
      config_setting_t *object = config_setting_get_elem(setting, i);
      
      config_setting_lookup_float(object,  "mass",    &m2);
      config_setting_lookup_float(object,  "radius",  &R2);
      config_setting_lookup_float(object,  "gyr_rad", &gyr_rad_outer);
      config_setting_lookup_float(object,  "a",       &a_outer);
      config_setting_lookup_float(object,  "e",       &e_outer);
      config_setting_lookup_float(object,  "P_rot",   &P_rot_outer);
      config_setting_lookup_float(object,  "alpha",   &alpha_outer);
    }
  }
 
  ////////////////////////////////////////////////////////////
  //
  // Getting information about star
  //
  ////////////////////////////////////////////////////////////
  
  double m0;
  double R0;
  setting = config_lookup(&cfg, "bodies.star");
  if(setting != NULL){
    int count = config_setting_length(setting);
    int i;
    
    for(i = 0; i < count; ++i){
      config_setting_t *object = config_setting_get_elem(setting, i);
      
      config_setting_lookup_float(object,  "mass",    &m0);
      config_setting_lookup_float(object,  "radius",  &R0);
    }
  }
  
  ////////////////////////////////////////////////////////////
  //
  // Finding derived time canonical unit
  //
  ////////////////////////////////////////////////////////////
  double uM, uL, uT;
  uM = MSUN;
  uL = RSUN;
  uT = get_canon_units(uM, uL, "uT")[2];

  rest.m1       = m1 * MJUP / uM;
  rest.R1       = R1 * RJUP / uL;
  rest.gyr_rad_inner   = gyr_rad_inner;
  rest.a_inner         = a_inner * AU / uL;
  rest.e_inner         = e_inner;
  rest.P_rot_inner     = P_rot_inner * DAY / uT;
  rest.alpha_inner     = alpha_inner;
  rest.m2              = m2 * MJUP / uM;
  rest.R2              = R2 * RJUP / uL;
  rest.gyr_rad_outer   = gyr_rad_outer;
  rest.a_outer         = a_outer * AU / uL;
  rest.e_outer         = e_outer;
  rest.P_rot_outer     = P_rot_outer * DAY / uT;
  rest.alpha_outer     = alpha_outer;
  rest.m0       = m0 * MSUN / uM;
  rest.sim_name  = sim_name;
  rest.total_t     = total_t * YEAR / uT;
  rest.dt_store  = dt_store * DAY / uT;
  rest.dt_std  = dt_std * MIN / uT;
  rest.uM        = uM;
  rest.uL        = uL;
  rest.uT        = uT;
  
  return rest;
  config_destroy(&cfg);

}