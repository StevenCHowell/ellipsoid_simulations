#ifndef H_NPTMC
#define H_NPTMC

#include "ellipsoid.h"

#include <sasmol.h>
#include <sasio.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406 
#endif

class RunParameters;

class RunParameters{

  public:
      RunParameters() = default ;

      const bool flag_restart = true;
      const int Nstart = 1;

      float box_length = (RBEADTOBEAD+2.*RADIUS)*N1d*1.011*2.0;        // initial box length

      int number_of_steps = 100;

      const int dcd_save_frequency = 10 ;     // how often to save frames

      float delta_move_max[3] = {(RBEADTOBEAD+2.*RADIUS)*1.011, 2*M_PI, 2*M_PI/3.0}; // max step size for translation, rotation, and bending
      float delta_move[3] = {RADIUS*1.0, 2*M_PI*0.1, 2*M_PI/3.0*0.1}; // step size for translation, rotation, and bending
      float delta_boxlength = box_length/10. ;              // volume move step size (percent)

      int interval_move = 1 ;		// update interval for translation, rotation, and bending
      int interval_boxlength = 100 ;		// update interval for box length


      float temperature = 100./119.8 ;                // temperature (119.8 Kelvin in reduced units)

      float pressure = 0.5 ;                    // pressure (bar)
      float goal_pressure = 100000.;                    // pressure (41.9 MPa (419 bar) in reduced units)

      float boxlinv = 0.0 ;                     // inverse box length
      float density = 0.0 ;                     // density of box

      ~RunParameters() = default ;
} ;

extern void stop_here() ;

#endif
