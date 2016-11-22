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

    float sigma_Ar = 3.405 ;                 // used to put distances in reduced/unitless representation
    float diameter_of_Ar = 3.76 / sigma_Ar ; // van der Waals diameter (# in A, division converts to unitless)

    /*** set these parameters ***/
    public:
        RunParameters() = default ;
        std::string runname = "run_2" ;
        const std::string dcd_output_filename = "run_2.dcd" ;

    const int dcd_save_frequency = 10 ;               // how often to save frames
    const std::string input_filename = "run_1.pdb" ;  // input pdb
    const std::string output_filename = "run_2.pdb" ; // output pdb (last frame)
    float box_length = 52.7545533335145 / sigma_Ar ;  // initial box length (# in A, division converts to unitless) = 15.4932608909 unitless
    int number_of_steps = 40010 ;                     // number of steps
    float delta_translation = 0.20 ;                  // translational move step size
    float temperature = 0.71 ;                        // temperature (85 Kelvin in units of 119.8 K)
    float goal_pressure = 0.2 ;                       // pressure (unitless) <--- ??? do not believe these units ???
    float r_lj = diameter_of_Ar / 2.0 ;               // Lennard-Jones radius

    /*** do NOT set here, only initialize, these are overwritten ***/
    float translation_ratio = 0.0 ;          // set in nptmc.cpp to number_of_steps/100.0
    float volume_ratio = 0.0 ;               // set in nptmc.cpp to number_of_steps/100.0
    float pressure = 0.0 ;                   // pressure (bar)
    float volume = 0.0 ;                     // volume of box
    float inv_box_length = 0.0 ;             // inverse box length
    float density = 0.0 ;                    // density of box
    float atf = 0.0 ;                        // packing fraction of box
    float r_cutoff = 0.0 ;                   // max distance for potential and force calculations

    /*** these are no longer used ***/
    // const int start_type = 2 ;               // 1 == xyz, 2 == restart
    // float delta_rotation = 0.0 ;             // molecular rotation move step size
    // float delta_volume = 0.020 ;             // volume move step size (percent)
    // float rmin = 0.7 ;

    ~RunParameters() = default ;
} ;

extern void stop_here() ;

extern void sumup(sasmol::SasMol &mol,RunParameters &par,float &v12,float &v6,float &w12,float &w6) ;

extern void lj_energy(sasmol::SasMol &mol, float &x1, float &y1, float &z1, int &thisi, float &box_length,float &v12, float &v6, float &w12, float &w6) ;

extern float get_random_float(float &a, float &b) ;
