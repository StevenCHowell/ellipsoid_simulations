#include "nptmc.h"
#include <dcdio.h>
#include <util.h>
#include <fstream>

using namespace Eigen ;

/********* methods        ******************/

float get_random_float(float &a, float &b)
{
  return ((b-a)*((float)rand()/RAND_MAX))+a;
}

bool check_conflict_scale(sasmol::SasMol &mol, const float scale, const float boxlnew)
{

  float xi,yi,zi,xj,yj,zj;
  float rxij,ryij,rzij,rijsq;
  float radius_i, radius_j;

  int frame = 0 ;

  sasmol::SasMol mol_seg_i, mol_seg_j;
  const std::vector<boost::dynamic_bitset<>> & masks = mol._mask_segnames();
  const int number_of_segnames = masks.size();
  std::vector<sasmol::SasMol> v_mol_segs;
  for (int i_seg=0; i_seg<number_of_segnames; ++i_seg)
    {
      mol.copy_molecule_using_mask(mol_seg_i, masks[i_seg], frame);
      auto com_i = mol_seg_i.calc_com(frame);
      for (auto & tmp : com_i) tmp *= scale;
      mol_seg_i.move_to(frame, com_i);
      v_mol_segs.push_back(mol_seg_i);
    }

  sasmol::SasMol * p_mol_seg_i, * p_mol_seg_j;
  for (int i_seg=0; i_seg<number_of_segnames-1; ++i_seg)
    {
      p_mol_seg_i = &v_mol_segs[i_seg];
      for (int i=0; i<p_mol_seg_i->_natoms(); ++i)
        {
          xi = p_mol_seg_i->_x()(i, frame);
          yi = p_mol_seg_i->_y()(i, frame);
          zi = p_mol_seg_i->_z()(i, frame);
          radius_i = std::stof(p_mol_seg_i->_atom_beta()[i]);
          for (int j_seg=i_seg+1; j_seg<number_of_segnames; ++j_seg)
            {
              p_mol_seg_j = &v_mol_segs[j_seg];
              for (int j=0; j<p_mol_seg_j->_natoms(); ++j)
                {
                  xj = p_mol_seg_j->_x()(j, frame);
                  yj = p_mol_seg_j->_y()(j, frame);
                  zj = p_mol_seg_j->_z()(j, frame);
                  radius_j = std::stof(p_mol_seg_j->_atom_beta()[j]);

                  rxij = xj-xi ;
                  ryij = yj-yi ;
                  rzij = zj-zi ;

                  rxij=rxij-boxlnew*(round(rxij/boxlnew)) ;
                  ryij=ryij-boxlnew*(round(ryij/boxlnew)) ;
                  rzij=rzij-boxlnew*(round(rzij/boxlnew)) ;

                  rijsq=rxij*rxij+ryij*ryij+rzij*rzij ;

                  if (rijsq < pow((radius_i+radius_j),2.0) ) return true;
                }
            }
        }
    }

  return false;
}

////////////////////////////////////////////////////////
/// check for overlap based on the candidate coordinates from the mask sasmol object
////////////////////////////////////////////////////////
bool checkOverlapUsingMask(sasmol::SasMol & mol, sasmol::SasMol & mol_mask_candidate, const boost::dynamic_bitset<> & mask, const int frame, const float box_length)
{
  double x1, y1, z1, x2, y2, z2, rxij, ryij, rzij, rijsq;
  int i_mask_candidate=0;

  for (int i = 0; i<mask.size(); ++i)
    {
      if(mask[i])
        {
          x1 = mol_mask_candidate._x()(i_mask_candidate,frame);
          y1 = mol_mask_candidate._y()(i_mask_candidate,frame);
          z1 = mol_mask_candidate._z()(i_mask_candidate,frame);
          for(int j = 0 ; j < mol._natoms() ; ++j)
            {
              if( ! mask[j] )
                {
                  x2 = mol._x()(j,frame) ;
                  y2 = mol._y()(j,frame) ;
                  z2 = mol._z()(j,frame) ;

                  rxij = x2-x1 ;
                  ryij = y2-y1 ;
                  rzij = z2-z1 ;

                  rxij=rxij-box_length*(round(rxij/box_length)) ;
                  ryij=ryij-box_length*(round(ryij/box_length)) ;
                  rzij=rzij-box_length*(round(rzij/box_length)) ;

                  rijsq=rxij*rxij+ryij*ryij+rzij*rzij ;

                  if (rijsq < pow(std::stof(mol._atom_beta()[i])+std::stof(mol._atom_beta()[j]),2.0)) return true;
                }

            } // end of loop j over mol._natoms()
          ++i_mask_candidate;
        }
    }
  return false;
}

/********* main           ******************/

int main(){

  std::cout << "\n\n\n" ;

  util::compile_time(__FILE__, __func__, __DATE__, __TIME__ );

  util::pp("welcome to a new beginning") ;

  RunParameters par ;

  // set up sasmol object
  sasmol::SasMol mol_single;
  util::pp(">>> reading pdb") ;
  mol_single.read_pdb("triplet.pdb");
  int frame = 0;
  sasmol::SasMol mol;
  mol_single.duplicate_molecule(mol, frame, 10, 150.0);
  mol.center(frame);
  mol.write_pdb("triplets.pdb",frame);
  mol.initialize_children();

  util::pp(">>> writing pdb") ;

  //mol.write_pdb(par.output_filename,frame) ;

  mol.calc_mass() ;

  std::cout << "total mass = " << mol._total_mass() << std::endl ;
  std::cout << "number of atoms = " << mol._natoms() << std::endl ;
  std::cout << "starting box length = " << par.box_length << std::endl ;

  FILE *dcdoutfile = sasio::open_write_dcd_file(par.dcd_output_filename, mol._natoms(), par.number_of_steps);

  //    boxfile=open(runname+'_box.dat','w')
  //    pressfile=open(runname+'_press.dat','w')
  //    densityfile=open(runname+'_density.dat','w')

  par.volume = pow(par.box_length,3.0) ;
  par.boxlinv = 1.0/par.box_length ;
  par.density = mol._natoms()/par.volume ;

  par.translation_ratio = par.number_of_steps/100.0 ;
  par.volume_ratio = par.number_of_steps/100.0 ;

  float dboxmx = par.box_length/40.0 ;
  float drmax = par.delta_translation ;

  float beta = 1.0/par.temperature ;         // i.e. kb=1.0

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

  //    Array<float,Dynamic,Dynamic> rxiold ;
  //     rxiold.setZero(mol._natoms(), 1);

  float rxiold, ryiold, rziold ;
  float rxinew, ryinew, rzinew ;

  float neg_1 = -1.0, pos_1 = 1.0, zero = 0.0 ;

  float ratbox,dpv,dvol,delthb,rrbox,boxlnew ;

  int accepted_moves = 0 ;
  int count = 0 ;
  int this_step ;
  int overall_frame = 0 ;
  float ratio, bratio ;

  sasmol::SasMol mol_seg;
  std::vector<float> d_move(3);
  float d_theta;
  int rot_axis;

  std::ofstream fout("box_length.txt",std::ofstream::out);
  std::cout << "STARTING MC SIMULATION : " << par.number_of_steps << " steps " << std::endl ;
  fout << "\nSTEP\t\tbox_length\t\tdensity\t\tpressure"<<std::endl;

  for(this_step = 0 ; this_step < par.number_of_steps ; ++this_step)
    {
      std::cout << this_step << " " << std::flush ;

      if(this_step == 0)
        {
          frames_since_last_dcd_save = par.dcd_save_frequency ;
        }
      else
        {
          frames_since_last_dcd_save += 1 ;
        }
      for (auto & mask : mol._mask_segnames())
        {
          count++ ;
          m=m+1 ; tm=tm+1 ;

          mol.copy_molecule_using_mask(mol_seg, mask, frame);
          d_move = {(get_random_float(neg_1,pos_1)) * drmax, (get_random_float(neg_1,pos_1)) * drmax, (get_random_float(neg_1,pos_1)) * drmax};
          mol_seg.translate(frame, d_move);
          d_theta = (util::get_random_float(neg_1,pos_1)) * PI;
          rot_axis = int(util::get_random_float(0.0,3.0));
          if (rot_axis==0) mol_seg.rotate(frame, "x", d_theta);
          else if (rot_axis==1) mol_seg.rotate(frame, "y", d_theta);
          else if (rot_axis==2) mol_seg.rotate(frame, "z", d_theta);
          else util::Error("failed to generate rotation axis");
          if ( ! checkOverlapUsingMask(mol, mol_seg, mask, frame, par.box_length) )
            {
              mol.set_sasmolcoor_using_mask(mol_seg, frame, mask);
              acatma=acatma+1.0 ;
              accepted_moves ++ ;
            }

          par.pressure = par.density * par.temperature;

          acm = acm + 1.0 ;
          acp = acp + par.pressure ;
          acd = acd + par.density ;

          acpsq = acpsq + pow(par.pressure,2.0) ;
          acdsq = acdsq + pow(par.density,2.0) ;

        } // end of loop i over mol._natoms()

      v=v+1 ; tv=tv+1 ;

      boxlnew = par.box_length + (get_random_float(neg_1,pos_1)) * dboxmx ;

      ratbox = par.box_length / boxlnew ;
      rrbox = 1.0 / ratbox ;

      dpv = par.goal_pressure * (pow(boxlnew,3.0) - par.volume) ;
      dvol = 3.0 * par.temperature * mol._natoms() * log(ratbox) ; //log == ln
      delthb = beta * ( dpv + dvol ) ;

      if ( ! check_conflict_scale(mol, rrbox, boxlnew) )
        {
          if(delthb < 75.0)
            {
              if(delthb <= 0.0)
                {
                  mol._x() *= rrbox ;
                  mol._y() *= rrbox ;
                  mol._z() *= rrbox ;

                  par.box_length = boxlnew ;
                  acboxa = acboxa + 1.0 ;
                }
              else if(exp(-delthb) > get_random_float(zero,pos_1))
                {
                  mol._x() *= rrbox ;
                  mol._y() *= rrbox ;
                  mol._z() *= rrbox ;

                  par.box_length = boxlnew ;
                  acboxa = acboxa + 1.0 ;
                }
            }
        }

      par.boxlinv = 1.0/par.box_length ;
      par.volume = pow(par.box_length,3.0) ;
      par.density=mol._natoms()/par.volume ;

      par.pressure = par.density * par.temperature;

      acm = acm + 1.0 ;
      acp = acp + par.pressure ;
      acd = acd + par.density ;

      acpsq = acpsq + (pow(par.pressure,2.0)) ;
      acdsq = acdsq + (pow(par.density,2.0)) ;

      if(frames_since_last_dcd_save == par.dcd_save_frequency)
        {
          overall_frame += 1 ;
          mol.write_dcd_step(dcdoutfile,frame,overall_frame) ;
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
          ratio = acatma/(mol._natoms()*par.translation_ratio) ;
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

    } // end of loop step over number_of_steps

  std::cout << "\nnumber of moves = " << count << std::endl ;
  std::cout << "accepted moves = " << accepted_moves << std::endl ;

  std::cout << "ratio = " << accepted_moves / float(count) << std::endl ;

  mol.write_pdb(par.output_filename,frame) ;
  fout.close();

  close_dcd_write(dcdoutfile) ;

  std::cout << std::endl << "\n\n\n" << std::endl ;
  util::print_run_details() ;
  std::cout << " >>> DONE <<< " << std::endl << std::endl ;
  return 0 ;


}
