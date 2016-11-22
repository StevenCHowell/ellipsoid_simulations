#include "npt_mc.h" // located here
#include <dcdio.h>  // located /usr/local/sasmol_cpp/include
#include <util.h>   // located /usr/local/sasmol_cpp/include
#include <fstream>

using namespace Eigen ;

/********* methods    ******************/
void stop_here()
{
  exit(0) ;
}

float get_random_float(float &a, float &b)
{
  return ((b - a)*((float)rand()/RAND_MAX)) + a;
}

void lj_energy(sasmol::SasMol &mol, float &x1, float &y1, float &z1, int &thisi, float &box_length, float &v12, float &v6, float &w12, float &w6)
{
  float x2, y2, z2 ;
  float sr2, sr6, vij12, vij6, r2_cutoff ;
  float rxij, ryij, rzij, r2ij ;

  v12 = 0.0, v6 = 0.0, w12 = 0.0, w6 = 0 ;

  r2_cutoff = pow((box_length/2.0), 2.0) ;

  int frame = 0 ;

  for(int j = 0 ; j < mol._natoms() ; ++j)
    {
      if( j != thisi)
        {
          x2 = mol._x()(j, frame) ;
          y2 = mol._y()(j, frame) ;
          z2 = mol._z()(j, frame) ;

          rxij = x2 - x1 ;
          ryij = y2 - y1 ;
          rzij = z2 - z1 ;

          rxij = rxij - box_length*(round(rxij/box_length)) ;
          ryij = ryij - box_length*(round(ryij/box_length)) ;
          rzij = rzij - box_length*(round(rzij/box_length)) ;

          r2ij = rxij*rxij + ryij*ryij + rzij*rzij ;

          if (r2ij < r2_cutoff)
            {
              sr2 = 1.0/r2ij ;
              sr6 = sr2*sr2*sr2 ;
              vij12 = sr6*sr6 ;
              vij6 = -sr6 ;
              v12 = v12 + vij12 ;
              v6 = v6 + vij6 ;
              w12 = w12 + vij12 ;
              w6 = w6 + vij6*0.5 ;
            }
        }

    } // end of loop j over mol._natoms()

  v12 = 4.0*v12 ;
  v6 = 4.0*v6 ;
  w12 = 48.0*w12/3.0 ;
  w6 = 48.0*w6/3.0 ;

  return ;
}

void sum_up(sasmol::SasMol &mol, RunParameters &par, float &v12, float &v6, float &w12, float &w6)
{
  int frame = 0 ;
  float r2_cutoff = par.r_cutoff * par.r_cutoff ;

  float r2ij, sr2, sr6, vij12, vij6 ;
  float rxi, ryi, rzi, rxj, ryj, rzj ;
  float rxij, ryij, rzij ;

  v12 = 0.0, v6 = 0.0, w12 = 0.0, w6 = 0.0 ;

  for(int i = 0 ; i < mol._natoms() - 1 ; ++i)
    {
      rxi = mol._x()(i, frame); ryi = mol._y()(i, frame) ; rzi = mol._z()(i, frame) ;
      for(int j = i + 1 ; j < mol._natoms() ; ++j)
        {
          rxj = mol._x()(j, frame); ryj = mol._y()(j, frame) ; rzj = mol._z()(j, frame) ;

          rxij = rxi - rxj ;
          ryij = ryi - ryj ;
          rzij = rzi - rzj ;

          rxij = rxij - par.box_length*(round(rxij*par.inv_box_length)) ;
          ryij = ryij - par.box_length*(round(ryij*par.inv_box_length)) ;
          rzij = rzij - par.box_length*(round(rzij*par.inv_box_length)) ;

          r2ij = rxij*rxij + ryij*ryij + rzij*rzij ;

          if(r2ij < r2_cutoff)
            {
              sr2 = 1.0/r2ij ;
              sr6 = sr2*sr2*sr2 ;
              vij12 = sr6*sr6 ;
              vij6 = -sr6 ;
              v12 += vij12 ;
              v6 += vij6 ;
              w12 += vij12 ;
              w6 += vij6*0.5 ;
            }
        }
    }

  v12 = 4.0*v12 ;
  v6 = 4.0*v6 ;
  w12 = 48.0*w12/3.0 ;
  w6 = 48.0*w6/3.0 ;

  return ;
}

/************************* main *****************************/

int main()
{
  /*** Intro print Statements ***/
  std::cout << "\n\n\n" ;

  util::compile_time(__FILE__, __func__, __DATE__, __TIME__ );

  util::pp("welcome to a new beginning") ;

  /*** instantiate a RunParameters object labeled par (defined in the header file) ***/
  RunParameters par ;

  /*** instantiate a sasmol object labeled mol,
       then populate it with coordinates from par.input_filename ***/
  sasmol::SasMol mol ;

  util::pp(">>> reading pdb") ;

  mol.read_pdb(par.input_filename) ;

  //util::pp(">>> writing pdb") ;
  int frame = 0 ; // used when writing dcd
  //mol.write_pdb(par.output_filename, frame) ;

  // gather and print info about the run
  // calc_mass(), _total_mass(), and _natoms() are not defined in python sasmol???
  mol.calc_mass() ;

  std::cout << "total mass = " << mol._total_mass() << std::endl ;
  std::cout << "number of atoms = " << mol._natoms() << std::endl ;
  std::cout << "starting box length = " << par.box_length << std::endl ;

  // slightly different from similar method in python
  FILE *dcdoutfile = sasio::open_write_dcd_file(par.dcd_output_filename, mol._natoms(), par.number_of_steps);

  par.volume = pow(par.box_length, 3.0) ;
  par.inv_box_length = 1.0/par.box_length ;
  par.density = mol._natoms()/par.volume ;

  par.translation_ratio = par.number_of_steps/100.0 ; // frequency to re-assess translation move step size
  par.volume_ratio = par.number_of_steps/100.0 ;      // frequency to re-assess volume move step size ?

  par.r_cutoff = 0.5*par.box_length ; // why l/2? is this arbitrary?
  std::cout << "cutoff (r_cutoff) = " << par.r_cutoff << std::endl ;

  float dboxmx = par.box_length/40.0 ;  // box move step size ?
  float drmax = par.delta_translation ; // translational move step size
  float beta = 1.0/par.temperature ;         // reduced units: kb = 1.0

  float sr3 = pow((1.0/par.r_cutoff), 3.0) ;
  float sr9 = pow(sr3, 3.0) ;

  float vlrc12 = 8.0*M_PI*par.density*mol._natoms()*sr9/9.0 ;
  float vlrc6 = -8.0*M_PI*par.density*mol._natoms()*sr3/3.0 ;
  float vlrc = vlrc12 + vlrc6 ;

  float wlrc12 = 4.0*vlrc12 ;
  float wlrc6  = 2.0*vlrc6 ;
  float wlrc = wlrc12 + wlrc6 ;

  float acm = 0.0, acatma = 0.0, acboxa = 0.0, acv = 0.0, acp = 0.0, acd = 0.0 ;
  float acvsq = 0.0, acpsq = 0.0, acdsq = 0.0, flv = 0.0, flp = 0.0, fld = 0.0 ;

  float v12 = 0.0, v6 = 0.0, w12 = 0.0, w6 = 0.0 ;

  sum_up(mol, par, v12, v6, w12, w6) ;

  float vs = ( v12 + v6 + vlrc ) / mol._natoms() ;
  float ws = ( w12 + w6 + wlrc ) / mol._natoms() ;

  float ps = par.density * par.temperature + ( w12 + w6 + wlrc ) / par.volume ;

  v12 = v12 + vlrc12 ;
  v6 = v6 + vlrc6 ;
  w12 = w12 + wlrc12 ;
  w6 = w6 + wlrc6 ;

  util::pp(">>> initial energies") ;

  std::cout << "sr3 = " << sr3 << std::endl ;
  std::cout << "sr9 = " << sr9 << std::endl ;
  std::cout << "vlrc12 = " << vlrc12 << std::endl ;
  std::cout << "vlrc6 = " << vlrc6 << std::endl ;
  std::cout << "wlrc12 = " << wlrc12 << std::endl ;
  std::cout << "wlrc6 = " << wlrc6 << std::endl ;
  std::cout << "vlrc = " << vlrc << std::endl ;
  std::cout << "wlrc = " << wlrc << std::endl ;

  std::cout << "v12 = " << v12 << std::endl ;
  std::cout << "v6 = " << v6 << std::endl ;
  std::cout << "w12 = " << w12 << std::endl ;
  std::cout << "w6 = " << w6 << std::endl ;

  std::cout << "density = " << par.density << std::endl ;
  std::cout << "initial energy/atom = " << vs << std::endl ;
  std::cout << "initial w/atom = " << ws << std::endl ;
  std::cout << "initial pressure = " << ps << std::endl ;

  float molvolume = ((mol._natoms()*4.0*M_PI*(pow(par.r_lj, 3.0)))/3.0) ; // volume of all the atoms: N 4/3 pi r^3

  int m = 0, v = 0, tm = 0, tv = 0 ;

  float sumeta = 0.0 ;

  int fr = 1 ;

  //std::cout << "\nSTEP\t\tVN\t%ACC_MOV\t%ACC_VOL\t ETA\t\t<ETA>\t\tDENSITY\t\t<DENSITY>\tPRESSURE\t<PRESSURE>\tBOXL\n" << std::endl ;
  //print '%i\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%10.6f\t%8.6f\n' % (i + 1, acatma*100.0/m, acboxa*100.0/v, density, acd/(m + v), pressure, acp/(m + v), boxl)

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN LOOP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  int frames_since_last_dcd_save = 0 ;

  //    Array<float, Dynamic, Dynamic> rxiold ;
  //     rxiold.setZero(mol._natoms(), 1);

  float rxiold, ryiold, rziold ;
  float rxinew, ryinew, rzinew ;

  float neg_1 = -1.0, pos_1 = 1.0, zero = 0.0 ;

  float v12old, v6old, w12old, w6old ;
  float v12new, v6new, w12new, w6new ;

  float delv12, delv6, delw12, delw6, deltv, deltvb ;
  float vn ;

  float ratbox, dpv, dvol, delthb, rat6, rat12, rcutn, rrbox, boxlnew ;

  int accepted_moves = 0 ;
  int count = 0 ;
  int this_step ;
  int overall_frame = 0 ;
  float ratio, bratio ;

  std::ofstream fout("output/box_length.txt", std::ofstream::out);
  std::cout << "STARTING MC SIMULATION : " << par.number_of_steps << " steps " << std::endl ;
  fout << "#step  box_length  density  pressure"<<std::endl;

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
      for(int i = 0 ; i < mol._natoms() ; ++i)
        {
          count++ ;
          m = m + 1 ; tm = tm + 1 ;

          rxiold = mol._x()(i, frame) ;
          ryiold = mol._y()(i, frame) ;
          rziold = mol._z()(i, frame) ;

          lj_energy(mol, rxiold, ryiold, rziold, i, par.box_length, v12old, v6old, w12old, w6old) ;

          rxinew = rxiold + (get_random_float(neg_1, pos_1)) * drmax ;
          ryinew = ryiold + (get_random_float(neg_1, pos_1)) * drmax ;
          rzinew = rziold + (get_random_float(neg_1, pos_1)) * drmax ;

          rxinew = rxinew - par.box_length*(round(rxinew*par.inv_box_length)) ;
          ryinew = ryinew - par.box_length*(round(ryinew*par.inv_box_length)) ;
          rzinew = rzinew - par.box_length*(round(rzinew*par.inv_box_length)) ;

          lj_energy(mol, rxinew, ryinew, rzinew, i, par.box_length, v12new, v6new, w12new, w6new) ;

          delv12 = v12new - v12old ;
          delv6 = v6new - v6old ;
          delw12 = w12new - w12old ;
          delw6 = w6new - w6old ;
          deltv = delv12 + delv6 ;
          deltvb = beta * deltv ;

          if(deltvb < 75.0)
            {
              if(deltv <= 0.0)
                {
                  v12 = v12 + delv12 ;
                  v6 = v6 + delv6 ;
                  w12 = w12 + delw12 ;
                  w6 = w6 + delw6 ;

                  mol._x()(i, frame) = rxinew ;
                  mol._y()(i, frame) = ryinew ;
                  mol._z()(i, frame) = rzinew ;
                  acatma = acatma + 1.0 ;
                  accepted_moves++ ;
                  //                std::cout << std::endl << "deltv acc" << std::endl ;
                }
              else if(exp(-deltvb) > get_random_float(zero, pos_1))
                {
                  v12 = v12 + delv12 ;
                  v6 = v6 + delv6 ;
                  w12 = w12 + delw12 ;
                  w6 = w6 + delw6 ;
                  mol._x()(i, frame) = rxinew ;
                  mol._y()(i, frame) = ryinew ;
                  mol._z()(i, frame) = rzinew ;
                  acatma = acatma + 1.0 ;
                  accepted_moves++ ;
                  //                std::cout << std::endl << "boltz acc" << std::endl ;
                }
            }

          vn = ( v12 + v6 ) / mol._natoms() ;
          //vn = ( v12 + vlrc12 + v6 + vlrc6 ) / mol._natoms() ;
          par.pressure = par.density * par.temperature + ( w12 + w6 ) / par.volume ;
          //par.pressure = par.density * par.temperature + ( w12 + wlrc12 + w6 + wlrc6 ) / par.volume ;

          acm = acm + 1.0 ;
          acv = acv + vn ;
          acp = acp + par.pressure ;
          acd = acd + par.density ;

          acvsq = acvsq + pow(vn, 2.0) ;
          acpsq = acpsq + pow(par.pressure, 2.0) ;
          acdsq = acdsq + pow(par.density, 2.0) ;

        } // end of loop i over mol._natoms()


      v = v + 1 ; tv = tv + 1 ;

      boxlnew = par.box_length + (get_random_float(neg_1, pos_1)) * dboxmx ;

      ratbox = par.box_length / boxlnew ;
      rrbox = 1.0 / ratbox ;
      rcutn = par.r_cutoff * rrbox ;

      rat6 = pow(ratbox, 6.0) ;
      rat12 = rat6 * rat6 ;

      v12new = v12 * rat12 ;
      v6new = v6 * rat6 ;
      w12new = w12 * rat12 ;
      w6new = w6 * rat6 ;

      deltv = v12new + v6new - v12 - v6 ;
      dpv = par.goal_pressure * (pow(boxlnew, 3.0) - par.volume) ;
      dvol = 3.0 * par.temperature * mol._natoms() * log(ratbox) ; //log == ln
      delthb = beta * ( deltv + dpv + dvol ) ;

      if(delthb < 75.0)
        {
          if(delthb <= 0.0)
            {
              v12 = v12new ;
              v6  = v6new ;
              w12 = w12new ;
              w6  = w6new ;

              mol._x() *= rrbox ;
              mol._y() *= rrbox ;
              mol._z() *= rrbox ;

              par.box_length = boxlnew ;
              par.r_cutoff = rcutn ;
              acboxa = acboxa + 1.0 ;
            }
          else if(exp(-delthb) > get_random_float(zero, pos_1))
            {
              v12 = v12new ;
              v6  = v6new ;
              w12 = w12new ;
              w6  = w6new;

              mol._x() *= rrbox ;
              mol._y() *= rrbox ;
              mol._z() *= rrbox ;

              par.box_length = boxlnew ;
              par.r_cutoff = rcutn ;
              acboxa = acboxa + 1.0 ;
            }
        }

      par.inv_box_length = 1.0/par.box_length ;
      par.volume = pow(par.box_length, 3.0) ;
      par.density = mol._natoms()/par.volume ;

      //#vlrc12 = 8.0*pi*density*natoms*sr9/9.0
      //#vlrc6 = -8.0*pi*density*natoms*sr3/3.0
      //#wlrc12 = 4.0*vlrc12
      //#wlrc6 =    2.0*vlrc6
      vn = ( v12 + v6 ) /mol._natoms() ;
      //#vn = ( v12 + vlrc12 + v6 + vlrc6) /mol._natoms() ;
      par.pressure = par.density * par.temperature + ( w12 + w6 ) / par.volume ;
      //#pressure = density * temperature + ( w12 + wlrc12 + w6 + wlrc6) / volume

      par.atf = molvolume/pow(par.box_length, 3.0) ;

      acm = acm + 1.0 ;
      acv = acv + vn ;
      acp = acp + par.pressure ;
      acd = acd + par.density ;

      sumeta = sumeta + par.atf ;

      acvsq = acvsq + (pow(vn, 2.0)) ;
      acpsq = acpsq + (pow(par.pressure, 2.0)) ;
      acdsq = acdsq + (pow(par.density, 2.0)) ;

      if(frames_since_last_dcd_save == par.dcd_save_frequency)
        {
          overall_frame += 1 ;
          mol.write_dcd_step(dcdoutfile, frame, overall_frame) ;
          fout <<this_step<<" "<<par.box_length<<" "<<par.density <<" "<<par.pressure<<std::endl;
          frames_since_last_dcd_save = 0 ;
        }

      //        boxfile.write("%i\t%f\n" % (step, boxl))
      //        pressfile.write("%i\t%f\t%f\n" % (step, pressure, acp/(tm + tv)))
      //        densityfile.write("%i\t%f\t%f\n" % (step, density, acd/(tm + tv)))
      //        boxfile.flush() ; pressfile.flush() ; densityfile.flush()

      if(fmodf(float(this_step), (par.number_of_steps/float(100.0))) == 0)
        {
          if(fmodf(float(fr), 20.0) == 0)
            {
              std::cout << "\nSTEP\t\tVN\t%ACC_MOV\t%ACC_VOL\t ETA\t\t<ETA>\t\tDENSITY\t\t<DENSITY>\tPRESSURE\t<PRESSURE>\tBOXL\n" << std::endl ;
              std::cout << this_step << " "<<vn <<" "<< acatma*100.0/m <<" "<< acboxa*100.0/v <<" "<< par.atf <<" "<< sumeta/tv <<" "<< par.density <<" "<< acd/(tm + tv) <<" "<< par.pressure <<" "<< acp/(tm + tv) <<" "<< par.box_length << std::endl ;
            }

          fr += 1 ;
        }

      if(fmodf(float(this_step), par.translation_ratio) == 0)
        {
          ratio = acatma/(mol._natoms()*par.translation_ratio) ;
          if (ratio > 0.5)
            {
              drmax = drmax*1.05 ;
            }
          else
            {
              drmax = drmax*0.95 ;
            }
          acatma = 0.0 ;
          m = 0 ;
        }
      if(fmodf(this_step, par.volume_ratio) == 0)
        {
          bratio = acboxa/par.volume_ratio ;
          if (bratio > 0.5)
            {
              dboxmx = dboxmx*1.05 ;
            }
          else
            {
              dboxmx = dboxmx*0.95 ;
            }
          acboxa = 0.0 ;
          v = 0 ;
        }

      std::cout << "\npercent mcmoves accepted " << acatma*100.0/tm << std::endl ;
      std::cout << "percent mcvols accepted " << acboxa*100.0/tv << std::endl ;

      std::cout << par.pressure << "\t = current pressure" << std::endl ;
      std::cout << par.goal_pressure << "\t\t = goal pressure" << std::endl ;




    } // end of loop step over number_of_steps

  mol.write_dcd_step(dcdoutfile, frame, overall_frame) ;
  fout <<this_step<<" "<<par.box_length<<" "<<par.density<<" "<<par.pressure<<std::endl;
  fout.close();

  std::cout << "\nnumber of moves = " << count << std::endl ;
  std::cout << "accepted moves = " << accepted_moves << std::endl ;

  std::cout << "ratio = " << accepted_moves / float(count) << std::endl ;

  mol.write_pdb(par.output_filename, frame) ;

  close_dcd_write(dcdoutfile) ;

  //    stop_here() ;

  std::cout << std::endl << "\n\n\n" << std::endl ;
  util::print_run_details() ;
  std::cout << " >>> DONE <<< " << std::endl << std::endl ;
  return 0 ;


}
