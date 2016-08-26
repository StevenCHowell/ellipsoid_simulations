#include <sasmol.h>
#include <sasio.h>

#include <iostream>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406 
#endif

#define N1d 12
#define RADIUS 1.0
#define RSQCUT RADIUS*RADIUS*4
#define RBEADTOBEAD RADIUS*3.1934
#define RCENTERTOBEAD RADIUS*1.8437

#define RHH RSQCUT // cut off radius
#define RHL1 RSQCUT // cut off radius
#define RHL2 RSQCUT // cut off radius
#define RL1L1 RSQCUT // cut off radius
#define RL1L2 RSQCUT // cut off radius
#define RL2L2 RSQCUT // cut off radius

// globals aren't good but I will move it later
Eigen::Array<float,1,3> v_H(RCENTERTOBEAD, 0.0, 0.0);
Eigen::Array<float,1,3> v_L1(-RCENTERTOBEAD/2.0, RCENTERTOBEAD*1.732/2.0, 0.0);
Eigen::Array<float,1,3> v_L2(-RCENTERTOBEAD/2.0, -RCENTERTOBEAD*1.732/2.0, 0.0);


class RunParameters;

class RunParameters{

  public:
      RunParameters() = default ;
      std::string runname = "run_1" ;
      const std::string dcd_output_filename = "run_1.dcd" ;

      const std::string input_filename = "run_0.pdb" ;
      const std::string output_filename = "final.pdb" ;

      float box_length = (RBEADTOBEAD+2.*RADIUS)*N1d*1.011*2.0;        // initial box length

      int number_of_steps = 10000 ;

      const int dcd_save_frequency = 10 ;     // how often to save frames

      float delta_translation = 0.20*RADIUS ;          // translational move step size
      float delta_rotation = 0.0 ;              // molecular rotation move step size
      float delta_volume = 0.01 ;              // volume move step size (percent)

	float translation_ratio = 0.0 ;	// parameter to re-assess frequency of translation move step size
	float volume_ratio = 0.0 ;		//


      float temperature = 300./119.8 ;                // temperature (119.8 Kelvin in reduced units)

      float pressure = 0.5 ;                    // pressure (bar)
      float goal_pressure = 100./419.;                    // pressure (41.9 MPa (419 bar) in reduced units)

      float volume = 0.0 ;                      // volume of box
      float boxlinv = 0.0 ;                     // inverse box length
      float density = 0.0 ;                     // density of box
      float eta = 0.0 ;                     	// packing fraction of box

      float rcut = 0.0 ;
      float rmin = 0.7 ;

      ~RunParameters() = default ;
} ;

extern void stop_here() ;
