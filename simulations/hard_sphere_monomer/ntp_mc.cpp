#include "nptmc.h"
#include <dcdio.h>
#include <util.h>
#include <fstream>

#define r_cutoff 0.5 // cut off radius

using namespace Eigen ;

/********* methods        ******************/
void stop_here(){ exit(0) ;} ;

float get_random_float(float &a, float &b)
{
  return ((b-a)*((float)rand()/RAND_MAX))+a;
}

bool check_conflict_scale(sasmol::SasMol &mol, const float scale, const float boxlnew)
{

	float x1,y1,z1,x2,y2,z2 ;
	float rxij,ryij,rzij,r2ij ;

	int frame = 0 ;
    const float r2_cutoff = (2*r_cutoff)*(2*r_cutoff);

	for(int i = 0 ; i < mol._natoms()-1 ; ++i)
    {
        x1 = mol._x()(i,frame)*scale ;
        y1 = mol._y()(i,frame)*scale ;
        z1 = mol._z()(i,frame)*scale ;
	    for(int j = i+1 ; j < mol._natoms() ; ++j)
	    {
			x2 = mol._x()(j,frame)*scale ;
			y2 = mol._y()(j,frame)*scale ;
			z2 = mol._z()(j,frame)*scale ;

			rxij = x2-x1 ;
			ryij = y2-y1 ;
			rzij = z2-z1 ;

			rxij=rxij-boxlnew*(round(rxij/boxlnew)) ;
			ryij=ryij-boxlnew*(round(ryij/boxlnew)) ;
			rzij=rzij-boxlnew*(round(rzij/boxlnew)) ;

			r2ij=rxij*rxij+ryij*ryij+rzij*rzij ;

			if (r2ij < r2_cutoff) return true;
	    } // end of loop j over mol._natoms()

	} // end of loop i over mol._natoms()

	return false;
}

bool check_conflict_1atom(sasmol::SasMol &mol, float &x1, float &y1, float &z1, int &thisi, float &box_length)
{

	float x2,y2,z2 ;
	float rxij,ryij,rzij,r2ij ;

	int frame = 0 ;
    const float r2_cutoff = (2*r_cutoff)*(2*r_cutoff);

	for(int j = 0 ; j < mol._natoms() ; ++j)
	{
        if( j != thisi)
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

			r2ij=rxij*rxij+ryij*ryij+rzij*rzij ;

			if (r2ij < r2_cutoff) return true;
        }

	} // end of loop j over mol._natoms()

	return false;
}

/********* main           ******************/

int main(){

    std::cout << "\n\n\n" ;

    util::compile_time(__FILE__, __func__, __DATE__, __TIME__ );

    util::pp("welcome to a new beginning") ;

    RunParameters par ;

    sasmol::SasMol mol ;

    util::pp(">>> reading pdb") ;

    mol.read_pdb(par.input_filename) ;

    util::pp(">>> writing pdb") ;

    int frame = 0 ;
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
        for(int i = 0 ; i < mol._natoms() ; ++i)
        {
            count++ ;
            m=m+1 ; tm=tm+1 ;

            rxiold = mol._x()(i,frame) ;
            ryiold = mol._y()(i,frame) ;
            rziold = mol._z()(i,frame) ;

            rxinew = rxiold + (get_random_float(neg_1,pos_1)) * drmax ;
            ryinew = ryiold + (get_random_float(neg_1,pos_1)) * drmax ;
            rzinew = rziold + (get_random_float(neg_1,pos_1)) * drmax ;

            rxinew = rxinew-par.box_length*(round(rxinew*par.boxlinv)) ;
            ryinew = ryinew-par.box_length*(round(ryinew*par.boxlinv)) ;
            rzinew = rzinew-par.box_length*(round(rzinew*par.boxlinv)) ;

			if ( ! check_conflict_1atom(mol,rxinew,ryinew,rzinew,i,par.box_length) )
            {
                    mol._x()(i,frame) = rxinew ;
                    mol._y()(i,frame) = ryinew ;
                    mol._z()(i,frame) = rzinew ;
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

//    stop_here() ;

    std::cout << std::endl << "\n\n\n" << std::endl ;
    util::print_run_details() ;
    std::cout << " >>> DONE <<< " << std::endl << std::endl ;
    return 0 ;


}
