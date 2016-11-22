#include <sasmol.h>
#include <sasio.h>

#include <iostream>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406 
#endif

class RunParameters;

class RunParameters{

  public:
      RunParameters() = default ;
      std::string runname = "run_1" ;
      const std::string dcd_output_filename = "run_1.dcd" ;

      const std::string input_filename = "run_0.pdb" ;
      const std::string output_filename = "final.pdb" ;

      const int start_type = 2 ;                // 1 == xyz, 2 == restart

      float box_length = 3000.0;        // initial box length

      int number_of_steps = 10000 ;

      const int dcd_save_frequency = 10 ;     // how often to save frames

      float delta_translation = 40.0 ;          // translational move step size
      float delta_rotation = 0.0 ;              // molecular rotation move step size
      float delta_volume = 0.020 ;              // volume move step size (percent)

	float translation_ratio = 0.0 ;	// parameter to re-assess frequency of translation move step size
	float volume_ratio = 0.0 ;		//



      float temperature = 0.71 ;                // temperature (85 Kelvin in reduced units)

      float pressure = 0.5 ;                    // pressure (bar)
      float goal_pressure = 0.2 ;                    // pressure (bar)

      float volume = 0.0 ;                      // volume of box
      float boxlinv = 0.0 ;                     // inverse box length
      float density = 0.0 ;                     // density of box

      ~RunParameters() = default ;
} ;

extern void stop_here() ;

extern float get_random_float(float &a, float &b) ;
