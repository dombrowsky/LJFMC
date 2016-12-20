#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <iostream>
#include <cmath>
#include <vector>
#include "include/pcg_random.hpp"

// evaluation function
double LJ_calc (double x, double xs, double ep){
  double x_rel = xs / x;
  double LJP = ep *( pow( x_rel , 12 ) - 2 * pow( x_rel , 6 ) );
  return LJP;
}
long double ran_num_gen (double size){
  pcg_extras::seed_seq_from<std::random_device> seed_source;
  // initialize random number machine
  pcg32 rng(seed_source);
  // tune uniform distributer
  std::uniform_real_distribution<long double> uniform_dist(0, size);
  // get number
  long double out = uniform_dist(rng);
  return out;
}

double calculateSD(std::vector<double> data){
  double sum, mean, standardDeviation;
  int i;
  int length;
  length = data.size();
  for(i = 0; i < length; ++i){
    sum += data[i];
  }
  mean = sum/length;
  for(i = 0; i < length; ++i){
    standardDeviation += pow(data[i] - mean, 2);
  }
  return sqrt(standardDeviation / length);
}
int main(void){
  // Containers
  int max_steps= 2000;
  float Temperatur;
  int n_particles, n_steps, i, j, k, num, num_two, klm, nanana, descriptor, counter, step;
  double Box_Size, r, epsilon, r_t, dx, dy, dz, dx_1, dx_2, dy_1, dy_2, dz_1, dz_2, LJ_Sum_init, LJ_Sum, criterion, deltaE;
  std::vector<double> Energies(max_steps);
  std::vector<double> E_ref(20);
  double SD, mean;
  long double boltzman_factor_1, boltzman_factor;
  long double kB = 1.38064852 * pow(10,-23);
  FILE *fenergyout;
  FILE *foutput;
  // ask for input variables
  printf("****************** Lennard-Jones fluid Monte Carlo simulations ******************\n");
  printf("************************* Author: Maximilian Dombrowsky *************************\n");
  printf("Please enter simulation box size in nm:\n");
  scanf("%lf",&Box_Size);
  printf("Please enter number of particles in your desired system:\n");
  scanf("%d",&n_particles);
  printf("Please enter simulation temperature:\n");
  scanf("%f",&Temperatur);
  printf("Please enter equilibrium distance between your desired particles:\n");
  scanf("%lf",&r);
  printf("Please enter equilibrium potential energy of your desired particles:\n");
  scanf("%lf",&epsilon);
  printf("Please enter number of simulation steps:\n");
  scanf("%d",&n_steps);
  printf("Begin Simulation \n");
  // begin system creation
  // initialize particle container in x, y and z with size n_particles
  long double x_1[n_particles], y_1[n_particles],z_1[n_particles];
  long double x_2[n_particles], y_2[n_particles],z_2[n_particles];
  
  for(step=0;step<n_steps;step++){
    // generate random particle coordinates
    for(klm=0; klm<n_particles;klm++){
      x_1[klm] = ran_num_gen(Box_Size);
      y_1[klm] = ran_num_gen( Box_Size);
      z_1[klm] = ran_num_gen(Box_Size);
    }
    //calculate Lennard Jones Energy between two particles
    for(num=0; num<n_particles-1;num++){
      for(num_two=num+1; num_two<n_particles;num_two++){
        // Periodic boundary Conditions correction
        dx_1 = x_1[num] - x_1[num_two];
        dx_2 = dx_1 - Box_Size;
        dx = std::min(fabs(dx_1), fabs(dx_2));
        dy_1 = y_1[num] - y_1[num_two];
        dy_2 = dy_1 - Box_Size;
        dy = std::min(fabs(dy_1), fabs(dy_2));
        dz_1 = y_1[num] - z_1[num_two];
        dz_2 = dz_1 - Box_Size;
        dz = std::min(fabs(dz_1), fabs(dz_2));
        // minimum distance between two particles
        r_t = sqrt(( dx * dx ) + ( dy * dy ) + ( dz * dz ));
        // Lennard-Jones Potential between two particles
        LJ_Sum_init += LJ_calc( r_t, r, epsilon );
      }
    }
    nanana = 0;
    // descriptor = 1;
    counter = 0;
    do{
      printf("%i\n",nanana);
        for(klm=0; klm<n_particles;klm++){
        x_2[klm] = ran_num_gen(Box_Size);
        y_2[klm] = ran_num_gen( Box_Size);
        z_2[klm] = ran_num_gen(Box_Size);
      }
      for(num=0; num<n_particles-1;num++){
        for(num_two=num+1; num_two<n_particles;num_two++){
          // Periodic boundary Conditions correction
          dx_1 = x_2[num] - x_2[num_two];
          dx_2 = dx_1 - Box_Size;
          dx = std::min(fabs(dx_1), fabs(dx_2));
          dy_1 = y_2[num] - y_2[num_two];
          dy_2 = dy_1 - Box_Size;
          dy = std::min(fabs(dy_1), fabs(dy_2));
          dz_1 = y_2[num] - z_2[num_two];
          dz_2 = dz_1 - Box_Size;
          dz = std::min(fabs(dz_1), fabs(dz_2));
          // minimum distance between two particles
          r_t = sqrt(( dx * dx ) + ( dy * dy ) + ( dz * dz ));
          // Lennard-Jones Potential between two particles
          LJ_Sum += LJ_calc( r_t, r, epsilon );
        }
      }
      foutput = fopen("coords.dat","w");
      fprintf(foutput,"# xyz-style\n");
      for(i=0; i<n_particles; i++){
        fprintf(foutput,"%Lf,%Lf,%Lf", x_1[i], y_1[i], z_1[i]);
      }
      fclose(foutput);
      // calculate boltzman factor
      criterion = ran_num_gen(1);
      deltaE = LJ_Sum - LJ_Sum_init;
      boltzman_factor_1 = ( -deltaE / (kB * Temperatur) );
      boltzman_factor = exp(boltzman_factor_1);
      // Ensemble 2 is accepted if LJ_Sum is smaller than LJ_Sum_init
      if(LJ_Sum < LJ_Sum_init){
        for(i=0; i<n_particles;i++){
          x_1[i] = x_2[i];
          y_1[i] = y_2[i];
          z_1[i] = z_2[i];
        }
        Energies[nanana] = LJ_Sum;
        LJ_Sum_init = LJ_Sum;
        if(counter > 10){
          counter -= 10;
        }
        else{
          counter = 0;
        }
        nanana++;
      }
      // if LJ_Sum is bigger than LJ_Sum_init then the boltzman factor will be between 0 and 1
      // Ensemble 2 will be accepted if bolzman factor is bigger than random variable between 0 and 1 
      // => can oversample local minima
      else if(boltzman_factor > criterion){
        for(i=0; i<n_particles;i++){
          x_1[i] = x_2[i];
          y_1[i] = y_2[i];
          z_1[i] = z_2[i];
        }
        Energies[nanana] = LJ_Sum;
        LJ_Sum_init = LJ_Sum;
        if(counter > 10){
          counter -= 10;
        }
        else{
          counter = 0;
        }
        nanana++;
      }
      else{
        counter++;
      }
    } while (nanana < max_steps || counter < 50);
    
  // write coordinates to file
    foutput = fopen("coords.dat","a");
    fprintf(foutput,"# step %i\n", (step + 1));
    for(i=0; i<n_particles; i++){
        fprintf(foutput,"%Lf,%Lf,%Lf", x_1[i], y_1[i], z_1[i]);
    }
    fclose(foutput);
    // write energies to file
    fenergyout = fopen("energy.dat","a");
    fprintf(fenergyout,"# energy \n");
    for(i=0; i<nanana; i++){
      fprintf(fenergyout,"%lf,",Energies[i]);
    }
    fclose(fenergyout);
  }
  return 0;
}

// Graveyard
// stop loop if energy didn't change significantly in the last 20 accepted steps 
//if(nanana > 100){
//  for(i=0; i<20; i++){
//    E_ref[i] = Energies[ ( nanana - i ) ];
//    mean += Energies[ ( nanana - i ) ];
//  }
//  SD = calculateSD(E_ref);
//  mean = mean / 20;
//  if(SD < mean / 10){
//    descriptor = 0;
//  }
//}
// random number generator
//long double ran_num_gen (double size){
//  long double out;
//  std::random_device rd;
//  std::mt19937 gen(rd());
//  std::uniform_real_distribution<> dis(0,size);
//  out = dis(gen);
//  return out;
//}

