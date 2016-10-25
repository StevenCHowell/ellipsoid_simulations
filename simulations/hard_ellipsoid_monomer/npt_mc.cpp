#include "npt_mc.h"
#include <dcdio.h>
#include <util.h>
#include <fstream>
#include <Eigen/Geometry>


bool Q_profiling_overall=false;
bool Q_profiling_detail=false;

/********* methods        ******************/
void stop_here(){ exit(0) ;} ;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
float dist_2atoms_pbc(const Eigen::Array<float,1,3> & R1, const Eigen::Array<float,1,3> & R2, const float box_length)
{
  float rxij = R2(0,0)-R1(0,0) ;
  float ryij = R2(0,1)-R1(0,1) ;
  float rzij = R2(0,2)-R1(0,2) ;
  rxij=rxij-box_length*(round(rxij/box_length)) ;
  ryij=ryij-box_length*(round(ryij/box_length)) ;
  rzij=rzij-box_length*(round(rzij/box_length)) ;

  return rxij*rxij+ryij*ryij+rzij*rzij;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool check_conflict_1atom(const int natoms, const Eigen::Array<float,1,3> & R1_H, const Eigen::Array<float,1,3> & R1_L1, const Eigen::Array<float,1,3> & R1_L2, const int thisi, const float box_length, const Eigen::ArrayXXf & coor_H, const Eigen::ArrayXXf & coor_L1, const Eigen::ArrayXXf & coor_L2)
{
  Eigen::Array<float,1,3> R2_H, R2_L1, R2_L2;

  /*
    std::cout<<"R1_H: "<<R1_H<<std::endl;
    std::cout<<"R1_L1: "<<R1_L1<<std::endl;
    std::cout<<"R1_L2: "<<R1_L2<<std::endl;
  */
  struct timeval t1, t2;
  gettimeofday(&t1,NULL) ;

  // check this molecule first
  if (dist_2atoms_pbc(R1_H, R1_L1, box_length)<RHL1) return true;
  if (dist_2atoms_pbc(R1_H, R1_L2, box_length)<RHL2) return true;
  if (dist_2atoms_pbc(R1_L1, R1_L2, box_length)<RL1L2) return true;

  // then check against other molecules
  for(int j = 0 ; j < natoms; ++j)
    {
      if( j != thisi)
        {
          R2_H << coor_H(j,0),coor_H(j,1),coor_H(j,2) ;
          R2_L1 << coor_L1(j,0),coor_L1(j,1),coor_L1(j,2) ;
          R2_L2 << coor_L2(j,0),coor_L2(j,1),coor_L2(j,2) ;
          //util::profiling(t1,t2,"grab coor") ;
          /*
            std::cout<<"R2_H: "<<R2_H<<std::endl;
            std::cout<<"R2_L1: "<<R2_L1<<std::endl;
            std::cout<<"R2_L2: "<<R2_L2<<std::endl;
            std::cout<<dist_2atoms_pbc(R1_H, R2_H, box_length)<<std::endl;
            std::cout<<dist_2atoms_pbc(R1_H, R2_L1, box_length)<<std::endl;
            std::cout<<dist_2atoms_pbc(R1_H, R2_L2, box_length)<<std::endl;
            std::cout<<dist_2atoms_pbc(R1_L1, R2_H, box_length)<<std::endl;
            std::cout<<dist_2atoms_pbc(R1_L1, R2_L1, box_length)<<std::endl;
            std::cout<<dist_2atoms_pbc(R1_L1, R2_L2, box_length)<<std::endl;
            std::cout<<dist_2atoms_pbc(R1_L2, R2_H, box_length)<<std::endl;
            std::cout<<dist_2atoms_pbc(R1_L2, R2_L1, box_length)<<std::endl;
            std::cout<<dist_2atoms_pbc(R1_L2, R2_L2, box_length)<<std::endl;
          */

          if (dist_2atoms_pbc(R1_H, R2_H, box_length)<RHH) return true;
          if (dist_2atoms_pbc(R1_H, R2_L1, box_length)<RHL1) return true;
          if (dist_2atoms_pbc(R1_H, R2_L2, box_length)<RHL2) return true;
          if (dist_2atoms_pbc(R1_L1, R2_H, box_length)<RHL1) return true;
          if (dist_2atoms_pbc(R1_L1, R2_L1, box_length)<RL1L1) return true;
          if (dist_2atoms_pbc(R1_L1, R2_L2, box_length)<RL1L2) return true;
          if (dist_2atoms_pbc(R1_L2, R2_H, box_length)<RHL2) return true;
          if (dist_2atoms_pbc(R1_L2, R2_L1, box_length)<RL1L2) return true;
          if (dist_2atoms_pbc(R1_L2, R2_L2, box_length)<RL2L2) return true;
          //util::profiling(t1,t2,"calc distance") ;
        }

    } // end of loop j over mol._natoms()
  //util::profiling(t1,t2,"check_conflict_1atom") ;

  return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool check_conflict_scale(const int natoms, const Eigen::ArrayXXf & coor_center, const Eigen::ArrayXXf & coor_H, const Eigen::ArrayXXf & coor_L1, const Eigen::ArrayXXf & coor_L2, const float scale, const float boxlnew)
{
  Eigen::Array<float,1,3> R1_center, R1_H, R1_L1, R1_L2, R2_center, R2_H, R2_L1, R2_L2;

  for(int i = 0 ; i < natoms-1 ; ++i)
    {
      R1_center << coor_center(i,0), coor_center(i,1), coor_center(i,2) ;
      R1_H << coor_H(i,0), coor_H(i,1), coor_H(i,2) ;
      R1_L1 << coor_L1(i,0), coor_L1(i,1), coor_L1(i,2) ;
      R1_L2 << coor_L2(i,0), coor_L2(i,1), coor_L2(i,2) ;
      R1_H = R1_center*scale + (R1_H-R1_center) ;
      R1_L1 = R1_center*scale + (R1_L1-R1_center) ;
      R1_L2 = R1_center*scale + (R1_L2-R1_center) ;
      for(int j = i+1 ; j < natoms; ++j)
        {
          R2_center << coor_center(j,0), coor_center(j,1), coor_center(j,2) ;
          R2_H << coor_H(j,0), coor_H(j,1), coor_H(j,2) ;
          R2_L1 << coor_L1(j,0), coor_L1(j,1), coor_L1(j,2) ;
          R2_L2 << coor_L2(j,0), coor_L2(j,1), coor_L2(j,2) ;
          R2_H = R2_center*scale + (R2_H-R2_center) ;
          R2_L1 = R2_center*scale + (R2_L1-R2_center) ;
          R2_L2 = R2_center*scale + (R2_L2-R2_center) ;

          if (dist_2atoms_pbc(R1_H, R2_H, boxlnew)<RHH) return true;
          if (dist_2atoms_pbc(R1_H, R2_L1, boxlnew)<RHL1) return true;
          if (dist_2atoms_pbc(R1_H, R2_L2, boxlnew)<RHL2) return true;
          if (dist_2atoms_pbc(R1_L1, R2_H, boxlnew)<RHL1) return true;
          if (dist_2atoms_pbc(R1_L1, R2_L1, boxlnew)<RL1L1) return true;
          if (dist_2atoms_pbc(R1_L1, R2_L2, boxlnew)<RL1L2) return true;
          if (dist_2atoms_pbc(R1_L2, R2_H, boxlnew)<RHL2) return true;
          if (dist_2atoms_pbc(R1_L2, R2_L1, boxlnew)<RL1L2) return true;
          if (dist_2atoms_pbc(R1_L2, R2_L2, boxlnew)<RL2L2) return true;
        } // end of loop j over mol._natoms()

    } // end of loop i over mol._natoms()

  return false;
}


/********* main           ******************/

int main(){

  std::cout << "\n\n\n" ;

  util::compile_time(__FILE__, __func__, __DATE__, __TIME__ ) ;

  util::pp("welcome to a new beginning") ;

  RunParameters par ;

  int frame = 0 ;

  // set up sasmol object
  sasmol::SasMol mol_single;
  util::pp(">>> reading pdb") ;
  mol_single.read_pdb("ellipsoid.pdb") ;

  util::pp(">>> creating copies of sasmol object") ;
  sasmol::SasMol mol;
  //sh note: N1d is a global variable defined in npt_mc.h
  mol_single.duplicate_molecule(mol, frame, N1d, par.box_length/N1d) ;

  util::pp(">>> centeringing composite system of copies") ;
  mol.center(frame) ;

  util::pp(">>> writing pdb of composite system to local dir") ;
  mol.write_pdb("ellipsoids.pdb",frame) ;
  //sh: need to write ANISOU line!!!
  mol.calc_mass() ;

  //sh: how do I want to represent ellipsoid orientation?!?
  //sh: need to randomize the orientations of starting ellipsoids!!!

  std::cout << "total mass = " << mol._total_mass() << std::endl ;
  std::cout << "number of atoms = " << mol._natoms() << std::endl ;
  std::cout << "starting box length = " << par.box_length << std::endl ;

  FILE *dcdoutfile = sasio::open_write_dcd_file(par.dcd_output_filename, mol_singlets._natoms(), par.number_of_steps) ;
  //sh: cannot use dcd to store orientation of ellipsoids!?!

  par.volume = pow(par.box_length,3.0) ;
  par.inv_box_length = 1.0/par.box_length ;
  par.density = mol_singlets._natoms()/par.volume ;

  par.translation_ratio = par.number_of_steps/100.0 ;
  par.volume_ratio = par.number_of_steps/100.0 ;

  float dboxmx = par.box_length/40.0 ;
  float drmax = par.delta_translation ;

  float beta = 1.0/par.temperature ;         // i.e., kb=1.0

  float acm = 0.0, acatma = 0.0, acboxa = 0.0, acp = 0.0, acd = 0.0 ;
  float acpsq = 0.0, acdsq = 0.0, flv = 0.0, flp = 0.0, fld = 0.0 ;

  float ps = par.density * par.temperature;

  std::cout << "density = " << par.density << std::endl ;
  std::cout << "initial pressure = " << ps << std::endl ;

  int m = 0, v = 0, tm = 0, tv = 0 ;

  int fr = 1 ;

  std::cout << "\nSTEP\t\tVN\t%ACC_MOV\t%ACC_VOL\t ETA\t\t<ETA>\t\tDENSITY\t\t<DENSITY>\tPRESSURE\t<PRESSURE>\tBOXL\n" << std::endl ;
  //print '%i\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%10.6f\t%8.6f\n' % (i+1,acatma*100.0/m,acboxa*100.0/v,density,acd/(m+v),pressure,acp/(m+v),boxl)

  ///// MAIN LOOP

  int frames_since_last_dcd_save = 0 ;

  float rxiold, ryiold, rziold ;
  float rxinew, ryinew, rzinew ;

  Eigen::Array<float,1,3> Rinew_H, Rinew_L1, Rinew_L2, Rinew_center;

  float neg_1 = -1.0, pos_1 = 1.0, zero = 0.0 ;

  float ratbox,dpv,dvol,delthb,rrbox,boxlnew ;

  int accepted_moves = 0 ;
  int count = 0 ;
  int this_step ;
  int overall_frame = 0 ;
  float ratio, bratio ;

  float d_theta;
  int rot_axis;
  Eigen::Vector3f v_axis;

  std::ofstream fout("box_length.txt",std::ofstream::out) ;
  std::cout << "STARTING MC SIMULATION : " << par.number_of_steps << " steps " << std::endl ;
  fout << "\nSTEP\t\tbox_length\t\tdensity\t\tpressure"<<std::endl;

  struct timeval t0_overall, t1_overall, t2_overall;
  struct timeval t1_detail, t2_detail;

  for(this_step = 0 ; this_step < par.number_of_steps ; ++this_step)
    {
      if (Q_profiling_overall) gettimeofday(&t0_overall,NULL) ;
      std::cout << this_step << " " << std::endl;

      if(this_step == 0)
        {
          frames_since_last_dcd_save = par.dcd_save_frequency ;
        }
      else
        {
          frames_since_last_dcd_save += 1 ;
        }

      if (Q_profiling_overall) gettimeofday(&t1_overall,NULL) ;
      for(int i = 0 ; i < mol._natoms() ; ++i)
        {
          count++ ;
          m=m+1 ; tm=tm+1 ;

          if (Q_profiling_detail) gettimeofday(&t1_detail,NULL) ;
          int dice = int(util::get_random_float(0.0,3.0)) ;
          // translate the singlet
          if (dice == 0)
            {
              rxiold = mol._x()(i,frame) ;
              ryiold = mol._y()(i,frame) ;
              rziold = mol._z()(i,frame) ;

              rxinew = rxiold + (util::get_random_float(neg_1,pos_1)) * drmax ;
              ryinew = ryiold + (util::get_random_float(neg_1,pos_1)) * drmax ;
              rzinew = rziold + (util::get_random_float(neg_1,pos_1)) * drmax ;

              rxinew = rxinew-par.box_length*(round(rxinew*par.inv_box_length)) ;
              ryinew = ryinew-par.box_length*(round(ryinew*par.inv_box_length)) ;
              rzinew = rzinew-par.box_length*(round(rzinew*par.inv_box_length)) ;

              Rinew_center = coor_center.row(i) ;
              Rinew_H = coor_H.row(i)-Rinew_center;
              Rinew_L1 = coor_L1.row(i)-Rinew_center;
              Rinew_L2 = coor_L2.row(i)-Rinew_center;
              Rinew_center << rxinew, ryinew, rzinew;
              Rinew_H = Rinew_center + Rinew_H;
              Rinew_L1 = Rinew_center + Rinew_L1;
              Rinew_L2 = Rinew_center + Rinew_L2;
              if (Q_profiling_detail) util::profiling(t1_detail,t2_detail,"shift singlet") ;
            }
          // rotate the singlet
          else if (dice==1)
            {
              rxinew = mol._x()(i,frame) ;
              ryinew = mol._y()(i,frame) ;
              rzinew = mol._z()(i,frame) ;

              d_theta = (util::get_random_float(neg_1,pos_1)) * PI;
              //d_theta = -PI/2.;
              rot_axis = int(util::get_random_float(0.0,3.0)) ;
              if (rot_axis==0) v_axis<<1.0,0.0,0.0;
              else if (rot_axis==1) v_axis<<0.0,1.0,0.0;
              else v_axis<<0.0,0.0,1.0;
              Eigen::AngleAxisf rot(d_theta, v_axis) ;

              Rinew_center = coor_center.row(i) ;
              Rinew_H = coor_H.row(i)-Rinew_center;
              Rinew_L1 = coor_L1.row(i)-Rinew_center;
              Rinew_L2 = coor_L2.row(i)-Rinew_center;
              Rinew_H = Rinew_center + (rot*(Rinew_H.matrix().transpose())).transpose().array() ;
              Rinew_L1 = Rinew_center + (rot*(Rinew_L1.matrix().transpose())).transpose().array() ;
              Rinew_L2 = Rinew_center + (rot*(Rinew_L2.matrix().transpose())).transpose().array() ;
              if (Q_profiling_detail) util::profiling(t1_detail,t2_detail,"rotate singlet") ;
            }
          // deform the singlet
          else
            {
              rxinew = mol._x()(i,frame) ;
              ryinew = mol._y()(i,frame) ;
              rzinew = mol._z()(i,frame) ;

              Rinew_center = coor_center.row(i) ;
              Rinew_H = coor_H.row(i)-Rinew_center;
              Rinew_L1 = coor_L1.row(i)-Rinew_center;
              Rinew_L2 = coor_L2.row(i)-Rinew_center;

              v_axis = (Rinew_H.matrix().transpose()).cross((Rinew_L1.matrix().transpose())) ;
              v_axis = v_axis/v_axis.norm() ;
              d_theta = (util::get_random_float(neg_1,pos_1)) * PI/6.0; // ZHL hardwired to rotate a small amount
              Eigen::AngleAxisf rot(d_theta, v_axis) ;

              rot_axis = int(util::get_random_float(0.0,3.0)) ;
              if (rot_axis==0)
                {
                  Rinew_H = Rinew_center + (rot*(Rinew_H.matrix().transpose())).transpose().array() ;
                  Rinew_L1 = Rinew_center + Rinew_L1;
                  Rinew_L2 = Rinew_center + Rinew_L2;
                }
              else if (rot_axis==1)
                {
                  Rinew_H = Rinew_center + Rinew_H;
                  Rinew_L1 = Rinew_center + (rot*(Rinew_L1.matrix().transpose())).transpose().array() ;
                  Rinew_L2 = Rinew_center + Rinew_L2;
                }
              else
                {
                  Rinew_H = Rinew_center + Rinew_H;
                  Rinew_L1 = Rinew_center + Rinew_L1;
                  Rinew_L2 = Rinew_center + (rot*(Rinew_L2.matrix().transpose())).transpose().array() ;
                }

              if (Q_profiling_detail) util::profiling(t1_detail,t2_detail,"deform singlet") ;
            }
          /*
            std::cout<<"v_axis: "<<v_axis<<std::endl;
            std::cout<<"d_theta: "<<d_theta<<std::endl;
            std::cout<<"Rinew_H: "<<Rinew_H<<std::endl;
            std::cout<<"Rinew_L1: "<<Rinew_L1<<std::endl;
            std::cout<<"Rinew_L2: "<<Rinew_L2<<std::endl;
          */

          if ( ! check_conflict_1atom(mol._natoms(), Rinew_H,Rinew_L1,Rinew_L2,i,par.box_length, coor_H, coor_L1, coor_L2) )
            {
              mol._x()(i,frame) = rxinew ;
              mol._y()(i,frame) = ryinew ;
              mol._z()(i,frame) = rzinew ;
              //util::profiling(t1,t2,"assign new sasmol") ;
              coor_center.row(i) <<rxinew, ryinew, rzinew;
              coor_H.row(i) = Rinew_H;
              coor_L1.row(i) = Rinew_L1;
              coor_L2.row(i) = Rinew_L2;
              //util::profiling(t1,t2,"assign new coor") ;
              acatma=acatma+1.0 ;
              accepted_moves ++ ;
              //                std::cout << std::endl << "deltv acc" << std::endl ;
            }

          par.pressure = par.density * par.temperature;

          acm = acm + 1.0 ;
          acp = acp + par.pressure ;
          acd = acd + par.density ;

          acpsq = acpsq + pow(par.pressure,2.0) ;
          acdsq = acdsq + pow(par.density,2.0) ;

        } // end of loop i over mol._natoms()
      if(Q_profiling_overall) util::profiling(t1_overall,t2_overall,"update per atoms") ;

      v=v+1 ; tv=tv+1 ;

      boxlnew = par.box_length + (util::get_random_float(neg_1,pos_1)) * dboxmx ;

      ratbox = par.box_length / boxlnew ;
      rrbox = 1.0 / ratbox ;

      dpv = par.goal_pressure * (pow(boxlnew,3.0) - par.volume) ;
      dvol = 3.0 * par.temperature * mol_singlets._natoms() * log(ratbox) ; //log == ln
      delthb = beta * ( dpv + dvol ) ;

      if ( ! check_conflict_scale(mol._natoms(), coor_center, coor_H, coor_L1, coor_L2, rrbox, boxlnew) )
        {
          if(delthb < 75.0)
            {
              if(delthb <= 0.0)
                {
                  mol._x() *= rrbox ;
                  mol._y() *= rrbox ;
                  mol._z() *= rrbox ;

                  coor_H = coor_center*rrbox + (coor_H-coor_center) ;
                  coor_L1 = coor_center*rrbox + (coor_L1-coor_center) ;
                  coor_L2 = coor_center*rrbox + (coor_L2-coor_center) ;
                  coor_center *= rrbox;

                  par.box_length = boxlnew ;
                  acboxa = acboxa + 1.0 ;
                }
              else if(exp(-delthb) > util::get_random_float(zero,pos_1))
                {
                  mol._x() *= rrbox ;
                  mol._y() *= rrbox ;
                  mol._z() *= rrbox ;

                  coor_H = coor_center*rrbox + (coor_H-coor_center) ;
                  coor_L1 = coor_center*rrbox + (coor_L1-coor_center) ;
                  coor_L2 = coor_center*rrbox + (coor_L2-coor_center) ;
                  coor_center *= rrbox;

                  par.box_length = boxlnew ;
                  acboxa = acboxa + 1.0 ;
                }
            }
        }
      if (Q_profiling_overall) util::profiling(t1_overall, t2_overall, "update scale") ;

      par.inv_box_length = 1.0/par.box_length ;
      par.volume = pow(par.box_length,3.0) ;
      par.density=mol_singlets._natoms()/par.volume ;

      par.pressure = par.density * par.temperature;

      acm = acm + 1.0 ;
      acp = acp + par.pressure ;
      acd = acd + par.density ;

      acpsq = acpsq + (pow(par.pressure,2.0)) ;
      acdsq = acdsq + (pow(par.density,2.0)) ;

      if(frames_since_last_dcd_save == par.dcd_save_frequency)
        {
          overall_frame += 1 ;
          for (int i=0; i<mol_singlets._natoms()/3; ++i)
            {
              mol_singlets._x()(i*3,frame) = coor_H(i,0) ;
              mol_singlets._y()(i*3,frame) = coor_H(i,1) ;
              mol_singlets._z()(i*3,frame) = coor_H(i,2) ;
              mol_singlets._x()(i*3+1,frame) = coor_L1(i,0) ;
              mol_singlets._y()(i*3+1,frame) = coor_L1(i,1) ;
              mol_singlets._z()(i*3+1,frame) = coor_L1(i,2) ;
              mol_singlets._x()(i*3+2,frame) = coor_L2(i,0) ;
              mol_singlets._y()(i*3+2,frame) = coor_L2(i,1) ;
              mol_singlets._z()(i*3+2,frame) = coor_L2(i,2) ;
            }
          mol_singlets.write_dcd_step(dcdoutfile,frame,overall_frame) ;
          frames_since_last_dcd_save = 0 ;
          fout << this_step << " "<<par.box_length<<" "<<par.density <<" "<<par.pressure<<std::endl;
        }

      //        boxfile.write("%i\t%f\n" % (step,boxl))
      //        pressfile.write("%i\t%f\t%f\n" % (step,pressure,acp/(tm+tv)))
      //        densityfile.write("%i\t%f\t%f\n" % (step,density,acd/(tm+tv)))
      //        boxfile.flush() ; pressfile.flush() ; densityfile.flush()

      if(fmodf(float(this_step),(par.number_of_steps/float(100.0)))==0)
        {
          if(fmodf(float(fr),20.0)==0)
            {
              std::cout << "\nSTEP\t\t%ACC_MOV\t%ACC_VOL\t\tDENSITY\t\t<DENSITY>\tPRESSURE\t<PRESSURE>\tBOXL\n" << std::endl ;
              std::cout << this_step << " "<< acatma*100.0/m <<" "<< acboxa*100.0/v <<" "<< par.density <<" "<< acd/(tm+tv) <<" "<< par.pressure <<" "<< acp/(tm+tv) <<" "<< par.box_length << std::endl ;
            }

          fr += 1 ;
        }

      if(fmodf(float(this_step),par.translation_ratio)==0)
        {
          ratio = acatma/(mol_singlets._natoms()*par.translation_ratio) ;
          if (ratio > 0.5)
            {
              drmax=drmax*1.05 ;
            }
          else
            {
              drmax=drmax*0.95 ;
            }
          acatma = 0.0 ;
          m=0 ;
        }
      if(fmodf(this_step,par.volume_ratio)==0)
        {
          bratio = acboxa/par.volume_ratio ;
          if (bratio > 0.5)
            {
              dboxmx=dboxmx*1.05 ;
            }
          else
            {
              dboxmx=dboxmx*0.95 ;
            }
          acboxa = 0.0 ;
          v=0 ;
        }

      std::cout << "\n\npercent mcmoves accepted " << acatma*100.0/tm << std::endl ;
      std::cout << "\n\npercent mcvols accepted " << acboxa*100.0/tv << std::endl ;

      std::cout << "goal pressure = " << par.goal_pressure << std::endl ;




      //    std::cout << std::endl << par.density << std::endl ;
      util::profiling(t0_overall, t2_overall, "MC step") ;

    } // end of loop step over number_of_steps

  std::cout << "\nnumber of moves = " << count << std::endl ;
  std::cout << "accepted moves = " << accepted_moves << std::endl ;

  std::cout << "ratio = " << accepted_moves / float(count) << std::endl ;

  mol.write_pdb(par.output_filename,frame) ;
  fout.close() ;





  close_dcd_write(dcdoutfile) ;

  //    stop_here() ;

  std::cout << std::endl << "\n\n\n" << std::endl ;
  util::print_run_details() ;
  std::cout << " >>> DONE <<< " << std::endl << std::endl ;
  return 0 ;


}
