#include "globvars.h"

int main(void)
{	
  double uM = MSUN;
  double uL = RSUN;
  double uT = get_canon_units(uM, uL, "uT")[2];

  printf("uM is: %f kg\n", uM);
  printf("uL is: %f m\n", uL);
  printf("uT is: %f s\n", uT);

  FILE *ftxv;

  char src[100];
  char dest[100];
  char file_name[100];
  strcpy(src, "../../sims/");
  strcpy(dest,"1y.dat");
  strcpy(file_name, strcat(src,dest));
  ftxv = fopen(file_name, "w");

  // ######################################################
  // # #### INTEGRATION TIME, TIMESTEP, AND SOFTENING #####
  // ######################################################
  int NBody = 3;
  double total_t = 1 * YEAR / uT;
  double dt_std = 0.2 * MIN / uT;
  double dt_store = 1. * DAY / uT;
  double softening = 0.4;

  printf("Integration time = %f y \n", total_t * uT / YEAR);
  printf("timestep = %f d\n", dt_std * uT / DAY);

  // #################################################
  // # #### SET INITIAL CONDITIONS OF THE SYSTEM #####
  // #################################################
  double m0 = MSUN / uM;
  double m1 = MSAT / uM;
  double m2 = MJUP / uM;

  double R0 = RSUN / uL;
  double R1 = RJUP / uL;
  double R2 = RJUP / uL;	

  double c = c_light * uT / uL;

  double e_inner = 0.0;
  double a_inner = 0.017 * AU / uL;
  double ini_pos_inner_x = a_inner * (1 - e_inner * e_inner) / (1 + e_inner);
  double ini_pos_inner_y = 0.0;
  double ini_vel_inner_x = 0.0;
  double ini_vel_inner_y = initial_velocity(m0, m1, ini_pos_inner_x, a_inner);

  double e_outer = 0.2;
  double a_outer = 0.08 * AU / uL;
  double ini_pos_outer_x = a_outer * (1 - e_outer * e_outer) / (1 + e_outer);
  double ini_pos_outer_y = 0.0;
  double ini_vel_outer_x = 0.0;
  double ini_vel_outer_y = initial_velocity(m0, m2, ini_pos_outer_x, a_outer);


  double ini_pos_star_x = 0.0;
  double ini_pos_star_y = 0.0;


  // #############################################################
  // # #### COMPUTE INITIAL COORDINATES OF THE CENTER OF MASS ####
  // #############################################################
  double COM_x = (m0 * ini_pos_star_x + m1 * ini_pos_inner_x + m2 * ini_pos_outer_x) / (m0 + m1 + m2);
  double COM_y = (m0 * ini_pos_star_y + m1 * ini_pos_inner_y + m2 * ini_pos_outer_y) / (m0 + m1 + m2);
  double COM[1][2];
  COM[0][0] = COM_x;
  COM[0][1] = COM_y;

  double star_positions[1][2];
  star_positions[0][0] = ini_pos_star_x;
  star_positions[0][1] = ini_pos_star_y;

  double inner_positions[1][2];
  inner_positions[0][0] = ini_pos_inner_x - COM_x;
  inner_positions[0][1] = ini_pos_inner_y - COM_y;

  double outer_positions[1][2];
  outer_positions[0][0] = ini_pos_outer_x - COM_x;
  outer_positions[0][1] = ini_pos_outer_y - COM_y;

  double star_velocities[1][2];
  star_velocities[0][0] = 0.0;
  star_velocities[0][1] = 0.0;

  double inner_velocities[1][2];
  inner_velocities[0][0] = ini_vel_inner_x;
  inner_velocities[0][1] = ini_vel_inner_y;

  double outer_velocities[1][2];
  outer_velocities[0][0] = ini_vel_outer_x;
  outer_velocities[0][1] = ini_vel_outer_y;

  // #########################################
  // # #### SEPARATION VECTORS FOR INNER ####
  // #########################################
  double SI_separation_vectors[1][2];
  SI_separation_vectors[0][0] = inner_positions[0][0] - star_positions[0][0];
  SI_separation_vectors[0][1] = inner_positions[0][1] - star_positions[0][1];

  double IS_separation_vectors[1][2];
  IS_separation_vectors[0][0] = star_positions[0][0] - inner_positions[0][0];
  IS_separation_vectors[0][1] = star_positions[0][1] - inner_positions[0][1];

  double SI_separations[1];
  SI_separations[0] = mag_vec_soft(SI_separation_vectors[0][0], SI_separation_vectors[0][1], softening);

  // ########################################
  // # #### SEPARATION VECTORS FOR OUTER ####
  // ########################################
  double SO_separation_vectors[1][2];
  SO_separation_vectors[0][0] = outer_positions[0][0] - star_positions[0][0];
  SO_separation_vectors[0][1] = outer_positions[0][1] - star_positions[0][1];

  double OS_separation_vectors[1][2];
  OS_separation_vectors[0][0] = star_positions[0][0] - outer_positions[0][0];
  OS_separation_vectors[0][1] = star_positions[0][1] - outer_positions[0][1];

  double SO_separations[1];
  SO_separations[0] = mag_vec_soft(SO_separation_vectors[0][0], SO_separation_vectors[0][1], softening);

  // ###############################################
  // # #### SEPARATIONS VECTORS FOR INNER-OUTER ####
  // ###############################################
  double IO_separation_vectors[1][2];
  IO_separation_vectors[0][0] = outer_positions[0][0] - inner_positions[0][0];
  IO_separation_vectors[0][1] = outer_positions[0][1] - inner_positions[0][1];

  double OI_separation_vectors[1][2];
  OI_separation_vectors[0][0] = inner_positions[0][0] - outer_positions[0][0];
  OI_separation_vectors[0][1] = inner_positions[0][1] - outer_positions[0][1];

  double IO_separations[1];
  IO_separations[0] = mag_vec_soft(IO_separation_vectors[0][0], IO_separation_vectors[0][1], softening);

  // ################################
  // # #### FORCES OF THE SYSTEM ####
  // ################################
  double star_forces[1][2];
  star_forces[0][0] =	grav_force(m0, m1, IS_separation_vectors, SI_separations[0])[0] + grav_force(m0, m2, OS_separation_vectors, SO_separations[0])[0];
  star_forces[0][1] = grav_force(m0, m1, IS_separation_vectors, SI_separations[0])[1] + grav_force(m0, m2, OS_separation_vectors, SO_separations[0])[1];

  double inner_forces[1][2];
  inner_forces[0][0] = grav_force(m1, m0, SI_separation_vectors, SI_separations[0])[0] + grav_force(m1, m2, OI_separation_vectors, IO_separations[0])[0] + grav_relat_force(m1, m0, c, SI_separation_vectors, SI_separations[0], inner_velocities)[0] * (m0 + m1) / m0;
  inner_forces[0][1] = grav_force(m1, m0, SI_separation_vectors, SI_separations[0])[1] + grav_force(m1, m2, OI_separation_vectors, IO_separations[0])[1] + grav_relat_force(m1, m0, c, SI_separation_vectors, SI_separations[0], inner_velocities)[1] * (m0 + m1) / m0;

  double outer_forces[1][2];
  outer_forces[0][0] = grav_force(m2, m0, SO_separation_vectors, SO_separations[0])[0] + grav_force(m2, m1, IO_separation_vectors, IO_separations[0])[0] + grav_relat_force(m1, m0, c, SI_separation_vectors, SI_separations[0], inner_velocities)[0] * m2 / m0;
  outer_forces[0][1] = grav_force(m2, m0, SO_separation_vectors, SO_separations[0])[1] + grav_force(m2, m1, IO_separation_vectors, IO_separations[0])[1] + grav_relat_force(m1, m0, c, SI_separation_vectors, SI_separations[0], inner_velocities)[1] * m2 / m0;

  double new_star_forces[1][2];
  double new_inner_forces[1][2];
  double new_outer_forces[1][2];

  double t, t_save;
  double progress;
  int i, j;

  double dt, dt0, dt1, dt2;

  double timestep_array[NBody + 1];

  t = 0.0;
  t_save = 0.0;
  fprintf(ftxv, "%10s\t%34s\t%34s\t%34s\t%34s\t%34s\t%34s\t%34s\t%34s\t%34s\t%34s\t%34s\t%34s\t%34s\n", 
                 "#Time", "dt", "Star_x", " Star_y","Inner_x", " Inner_y", " Outer_x", "Outer_y", "Star_vx", 
                 "Star_vy", "Inner_vx", "Inner_vy", "Outer_vx", "Outer_vy");

  while (t <= total_t){

    // Compute timestep for next calculations		
    timestep_array[0] = dt_std;
    timestep_array[1] = timestep(star_forces, m0, R0, dt_std);
    timestep_array[2] = timestep(inner_forces, m1, R1, dt_std);
    timestep_array[3] = timestep(outer_forces, m2, R2, dt_std);
    // Choose the minimum of the above timesteps
    dt = find_min(timestep_array, NBody + 1);

    if (t > t_save){
      fprintf(ftxv, "%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\n", 
              t, dt, star_positions[0][0], star_positions[0][1], inner_positions[0][0], inner_positions[0][1],
              outer_positions[0][0], outer_positions[0][1], star_velocities[0][0], star_velocities[0][1],
              inner_velocities[0][0], inner_velocities[0][1], outer_velocities[0][0], outer_velocities[0][1]);


      t_save += dt_store;
    }


    progress = (t / total_t) * 100.0;
    printf("Progress: %.6f per cent \n", progress);
    // printf("%.5f %.5f\n", dt_std, dt);
    for (i=0; i<2; i++){		
      // Update positions
      star_positions[0][i] = star_positions[0][i] + star_velocities[0][i] * dt + 0.5 * star_forces[0][i] / m0 * dt*dt;
      inner_positions[0][i] = inner_positions[0][i] + inner_velocities[0][i] * dt + 0.5 * inner_forces[0][i] / m1 * dt*dt;
      outer_positions[0][i] = outer_positions[0][i] + outer_velocities[0][i] * dt + 0.5 * outer_forces[0][i] / m2 * dt*dt;

      COM[0][i] = (m0 * star_positions[0][i] + m1 * inner_positions[0][i] + m2 * outer_positions[0][i]) / (m0 + m1 + m2);

      // Reset positions to common frame of reference
      star_positions[0][i] = star_positions[0][i] - COM[0][i];
      inner_positions[0][i] = inner_positions[0][i] - COM[0][i];
      outer_positions[0][i] = outer_positions[0][i] - COM[0][i];

      // Calculate separation vectors
      SI_separation_vectors[0][i] = inner_positions[0][i] - star_positions[0][i];
      IS_separation_vectors[0][i] = -1.0 * SI_separation_vectors[0][i];

      SO_separation_vectors[0][i]= outer_positions[0][i] - star_positions[0][i];
      OS_separation_vectors[0][i] = -1.0 * SO_separation_vectors[0][i];

      IO_separation_vectors[0][i] = outer_positions[0][i] - inner_positions[0][i];
      OI_separation_vectors[0][i] = -1.0 * IO_separation_vectors[0][i];
    }

    // Compute magnitude of the separation vectors
    SI_separations[0] = mag_vec_soft(SI_separation_vectors[0][0], SI_separation_vectors[0][1], softening);
    SO_separations[0] = mag_vec_soft(SO_separation_vectors[0][0], SO_separation_vectors[0][1], softening);
    IO_separations[0] = mag_vec_soft(IO_separation_vectors[0][0], IO_separation_vectors[0][1], softening);

    for (j=0; j<2; j++){
      // Compute new forces
      new_star_forces[0][j] = grav_force(m0, m1, IS_separation_vectors, SI_separations[0])[j] + grav_force(m0, m2, OS_separation_vectors, SO_separations[0])[j];
      new_inner_forces[0][j] = grav_force(m1, m0, SI_separation_vectors, SI_separations[0])[j] + grav_force(m1, m2, OI_separation_vectors, IO_separations[0])[j] + grav_relat_force(m1, m0, c, SI_separation_vectors, SI_separations[0], inner_velocities)[j] * (m0 + m1) / m0;
      new_outer_forces[0][j] = grav_force(m2, m0, SO_separation_vectors, SO_separations[0])[j] + grav_force(m2, m1, IO_separation_vectors, IO_separations[0])[j] + grav_relat_force(m1, m0, c, SI_separation_vectors, SI_separations[0], inner_velocities)[j] * m2 / m0;

      // Update velocities
      star_velocities[0][j] = star_velocities[0][j] + 0.5 * (star_forces[0][j] + new_star_forces[0][j]) / m0 * dt;
      inner_velocities[0][j] = inner_velocities[0][j] + 0.5 * (inner_forces[0][j] + new_inner_forces[0][j]) / m1 * dt;
      outer_velocities[0][j] = outer_velocities[0][j] + 0.5 * (outer_forces[0][j] + new_outer_forces[0][j]) / m2 * dt;

      star_forces[0][j] = new_star_forces[0][j];
      inner_forces[0][j] = new_inner_forces[0][j];
      outer_forces[0][j] = new_outer_forces[0][j];
    }

    t += dt;
  }

  printf("prueba %.30f\n", inner_forces[0][0]);

  fclose(ftxv);
  return 0;
}
