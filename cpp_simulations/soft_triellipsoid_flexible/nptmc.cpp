#include "nptmc.h"

#include <dcdio.h>
#include <util.h>
#include <fstream>
#include <Eigen/Geometry>


bool Q_profiling_overall=true;
bool Q_profiling_detail=false;

/********* methods        ******************/
void stop_here(){ exit(0) ;} ;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initHLCoor(sasmol::SasMol & mol, EigenArrayXXfRow & coor_center, EigenArrayXXfRow & coor_H, EigenArrayXXfRow & coor_L1, EigenArrayXXfRow & coor_L2,
                EigenArrayXXfRow & u0_H, EigenArrayXXfRow & u1_H, EigenArrayXXfRow & u2_H,
                EigenArrayXXfRow & u0_L1, EigenArrayXXfRow & u1_L1, EigenArrayXXfRow & u2_L1,
                EigenArrayXXfRow & u0_L2, EigenArrayXXfRow & u1_L2, EigenArrayXXfRow & u2_L2)
{
    coor_center = mol._coor();

    float d_theta;
    int rot_axis;
    Eigen::Vector3f v_axis;
    for (int i=0; i<mol._natoms(); ++i)
    {
        d_theta = (util::get_random_float(-1.0, 1.0)) * PI;
        rot_axis = int(util::get_random_float(0.0,3.0));
        if (rot_axis==0) v_axis<<1.0,0.0,0.0;
        else if (rot_axis==1) v_axis<<0.0,1.0,0.0;
        else v_axis<<0.0,0.0,1.0;
        Eigen::AngleAxisf rot(d_theta, v_axis);
        coor_H.row(i) = coor_center.row(i) + (rot*(VECTOR_CENTER_TO_H.matrix().transpose())).transpose().array();
        coor_L1.row(i) = coor_center.row(i) + (rot*(VECTOR_CENTER_TO_L1.matrix().transpose())).transpose().array();
        coor_L2.row(i) = coor_center.row(i) + (rot*(VECTOR_CENTER_TO_L2.matrix().transpose())).transpose().array();
        u0_H.row(i) = (rot*(U0_H.matrix().transpose())).transpose().array();
        u1_H.row(i) = (rot*(U1_H.matrix().transpose())).transpose().array();
        u2_H.row(i) = (rot*(U2_H.matrix().transpose())).transpose().array();
        u0_L1.row(i) = (rot*(U0_L1.matrix().transpose())).transpose().array();
        u1_L1.row(i) = (rot*(U1_L1.matrix().transpose())).transpose().array();
        u2_L1.row(i) = (rot*(U2_L1.matrix().transpose())).transpose().array();
        u0_L2.row(i) = (rot*(U0_L2.matrix().transpose())).transpose().array();
        u1_L2.row(i) = (rot*(U1_L2.matrix().transpose())).transpose().array();
        u2_L2.row(i) = (rot*(U2_L2.matrix().transpose())).transpose().array();
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
float dist_2atoms_pbc(const EigenArray13fRow & R1, const EigenArray13fRow & R2, const float box_length)
{
    float rxij = R2(0,0)-R1(0,0);
    float ryij = R2(0,1)-R1(0,1);
    float rzij = R2(0,2)-R1(0,2);
    rxij=rxij-box_length*(round(rxij/box_length)) ;
    ryij=ryij-box_length*(round(ryij/box_length)) ;
    rzij=rzij-box_length*(round(rzij/box_length)) ;
                         
	return rxij*rxij+ryij*ryij+rzij*rzij;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// For sphere
////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool check_conflict_1atom(const int natoms, const EigenArray13fRow & R1_H, const EigenArray13fRow & R1_L1, const EigenArray13fRow & R1_L2, const int thisi, const float box_length, const EigenArrayXXfRow & coor_H, const EigenArrayXXfRow & coor_L1, const EigenArrayXXfRow & coor_L2)
{
    EigenArray13fRow R2_H, R2_L1, R2_L2;

    /*
	std::cout<<"R1_H: "<<R1_H<<std::endl;
	std::cout<<"R1_L1: "<<R1_L1<<std::endl;
	std::cout<<"R1_L2: "<<R1_L2<<std::endl;
    */
	struct timeval t1, t2;
	gettimeofday(&t1,NULL);

    // check this molecule first
    if (dist_2atoms_pbc(R1_H, R1_L1, box_length)<RHL1) return true;
    if (dist_2atoms_pbc(R1_H, R1_L2, box_length)<RHL2) return true;
    if (dist_2atoms_pbc(R1_L1, R1_L2, box_length)<RL1L2) return true;

    // then check against other molecules
    for(int j = 0 ; j < natoms; ++j)
    {
        if( j != thisi)
        { 
			R2_H << coor_H(j,0),coor_H(j,1),coor_H(j,2);
			R2_L1 << coor_L1(j,0),coor_L1(j,1),coor_L1(j,2);
			R2_L2 << coor_L2(j,0),coor_L2(j,1),coor_L2(j,2);
	        //util::profiling(t1,t2,"grab coor");
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
	        //util::profiling(t1,t2,"calc distance");
        }

    } // end of loop j over mol._natoms() 
	//util::profiling(t1,t2,"check_conflict_1atom");

    return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// For ellipsoid
////////////////////////////////////////////////////////////////////////////////////////////////////////////
float lj_energy_ellipsoid(const int natoms, const EigenArray13fRow & R1_H, const EigenArray13fRow & R1_L1, const EigenArray13fRow & R1_L2,
                          const EigenArray13fRow & u0inew_H, const EigenArray13fRow & u1inew_H, const EigenArray13fRow & u2inew_H,
                          const EigenArray13fRow & u0inew_L1, const EigenArray13fRow & u1inew_L1, const EigenArray13fRow & u2inew_L1, 
                          const EigenArray13fRow & u0inew_L2, const EigenArray13fRow & u1inew_L2, const EigenArray13fRow & u2inew_L2, 
                          const int thisi, const float box_length, const EigenArrayXXfRow & coor_H, const EigenArrayXXfRow & coor_L1, const EigenArrayXXfRow & coor_L2,
                          const EigenArrayXXfRow & u0_H, const EigenArrayXXfRow & u1_H, const EigenArrayXXfRow & u2_H,
                          const EigenArrayXXfRow & u0_L1, const EigenArrayXXfRow & u1_L1, const EigenArrayXXfRow & u2_L1, 
                          const EigenArrayXXfRow & u0_L2, const EigenArrayXXfRow & u1_L2, const EigenArrayXXfRow & u2_L2)
{
    const Eigen::Vector3f v_R1_H = R1_H.transpose().matrix();
    const Eigen::Vector3f v_R1_L1 = R1_L1.transpose().matrix();
    const Eigen::Vector3f v_R1_L2 = R1_L2.transpose().matrix();
    Eigen::Vector3f v_R2_H, v_R2_L1, v_R2_L2;

    // get the internal energy for this molecule first
    float energy = 0.0;
    energy += pot_energy(v_R1_H, axis_lengths_H, u0inew_H.transpose().matrix(), u1inew_H.transpose().matrix(), u2inew_H.transpose().matrix(),
                         v_R1_L1, axis_lengths_L1, u0inew_L1.transpose().matrix(), u1inew_L1.transpose().matrix(), u2inew_L1.transpose().matrix(), sigma, epsilon, box_length);

    // then check against other molecules
    for(int j = 0 ; j < natoms; ++j)
    {
        if( j != thisi)
        { 
			v_R2_H << coor_H(j,0),coor_H(j,1),coor_H(j,2);
			v_R2_L1 << coor_L1(j,0),coor_L1(j,1),coor_L1(j,2);
			v_R2_L2 << coor_L2(j,0),coor_L2(j,1),coor_L2(j,2);
            energy += pot_energy(v_R1_H, axis_lengths_H, u0inew_H.transpose().matrix(), u1inew_H.transpose().matrix(), u2inew_H.transpose().matrix(),
                                 v_R2_H, axis_lengths_H, u0_H.row(j).transpose().matrix(), u1_H.row(j).transpose().matrix(), u2_H.row(j).transpose().matrix(), sigma, epsilon, box_length);
            energy += pot_energy(v_R1_H, axis_lengths_H, u0inew_H.transpose().matrix(), u1inew_H.transpose().matrix(), u2inew_H.transpose().matrix(),
                                 v_R2_L1, axis_lengths_L1, u0_L1.row(j).transpose().matrix(), u1_L1.row(j).transpose().matrix(), u2_L1.row(j).transpose().matrix(), sigma, epsilon, box_length);
            energy += pot_energy(v_R1_H, axis_lengths_H, u0inew_H.transpose().matrix(), u1inew_H.transpose().matrix(), u2inew_H.transpose().matrix(),
                                 v_R2_L2, axis_lengths_L2, u0_L2.row(j).transpose().matrix(), u1_L2.row(j).transpose().matrix(), u2_L2.row(j).transpose().matrix(), sigma, epsilon, box_length);
            energy += pot_energy(v_R1_L1, axis_lengths_L1, u0inew_L1.transpose().matrix(), u1inew_L1.transpose().matrix(), u2inew_L1.transpose().matrix(),
                                 v_R2_H, axis_lengths_H, u0_H.row(j).transpose().matrix(), u1_H.row(j).transpose().matrix(), u2_H.row(j).transpose().matrix(), sigma, epsilon, box_length);
            energy += pot_energy(v_R1_L1, axis_lengths_L1, u0inew_L1.transpose().matrix(), u1inew_L1.transpose().matrix(), u2inew_L1.transpose().matrix(),
                                 v_R2_L1, axis_lengths_L1, u0_L1.row(j).transpose().matrix(), u1_L1.row(j).transpose().matrix(), u2_L1.row(j).transpose().matrix(), sigma, epsilon, box_length);
            energy += pot_energy(v_R1_L1, axis_lengths_L1, u0inew_L1.transpose().matrix(), u1inew_L1.transpose().matrix(), u2inew_L1.transpose().matrix(),
                                 v_R2_L2, axis_lengths_L2, u0_L2.row(j).transpose().matrix(), u1_L2.row(j).transpose().matrix(), u2_L2.row(j).transpose().matrix(), sigma, epsilon, box_length);
            energy += pot_energy(v_R1_L2, axis_lengths_L2, u0inew_L2.transpose().matrix(), u1inew_L2.transpose().matrix(), u2inew_L2.transpose().matrix(),
                                 v_R2_H, axis_lengths_H, u0_H.row(j).transpose().matrix(), u1_H.row(j).transpose().matrix(), u2_H.row(j).transpose().matrix(), sigma, epsilon, box_length);
            energy += pot_energy(v_R1_L2, axis_lengths_L2, u0inew_L2.transpose().matrix(), u1inew_L2.transpose().matrix(), u2inew_L2.transpose().matrix(),
                                 v_R2_L1, axis_lengths_L1, u0_L1.row(j).transpose().matrix(), u1_L1.row(j).transpose().matrix(), u2_L1.row(j).transpose().matrix(), sigma, epsilon, box_length);
            energy += pot_energy(v_R1_L2, axis_lengths_L2, u0inew_L2.transpose().matrix(), u1inew_L2.transpose().matrix(), u2inew_L2.transpose().matrix(),
                                 v_R2_L2, axis_lengths_L2, u0_L2.row(j).transpose().matrix(), u1_L2.row(j).transpose().matrix(), u2_L2.row(j).transpose().matrix(), sigma, epsilon, box_length);
	        //util::profiling(t1,t2,"calc distance");
        }

    } // end of loop j over mol._natoms() 

    return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool check_conflict_scale(const int natoms, const EigenArrayXXfRow & coor_center, const EigenArrayXXfRow & coor_H, const EigenArrayXXfRow & coor_L1, const EigenArrayXXfRow & coor_L2, const float scale, const float boxlnew)
{
    EigenArray13fRow R1_center, R1_H, R1_L1, R1_L2, R2_center, R2_H, R2_L1, R2_L2;

    for(int i = 0 ; i < natoms-1 ; ++i)
    {
        R1_center << coor_center(i,0), coor_center(i,1), coor_center(i,2);
        R1_H << coor_H(i,0), coor_H(i,1), coor_H(i,2);
        R1_L1 << coor_L1(i,0), coor_L1(i,1), coor_L1(i,2);
        R1_L2 << coor_L2(i,0), coor_L2(i,1), coor_L2(i,2);
        R1_H = R1_center*scale + (R1_H-R1_center);
        R1_L1 = R1_center*scale + (R1_L1-R1_center);
        R1_L2 = R1_center*scale + (R1_L2-R1_center);
        for(int j = i+1 ; j < natoms; ++j)
        {
            R2_center << coor_center(j,0), coor_center(j,1), coor_center(j,2);
            R2_H << coor_H(j,0), coor_H(j,1), coor_H(j,2);
            R2_L1 << coor_L1(j,0), coor_L1(j,1), coor_L1(j,2);
            R2_L2 << coor_L2(j,0), coor_L2(j,1), coor_L2(j,2);
            R2_H = R2_center*scale + (R2_H-R2_center);
            R2_L1 = R2_center*scale + (R2_L1-R2_center);
            R2_L2 = R2_center*scale + (R2_L2-R2_center);
         
			if (dist_2atoms_pbc(R1_H, R2_H, boxlnew)<RHH) {std::cout<<"Conflict found between seg# "<<i<<" and seg# "<<j<<" during volume change"<<std::endl; return true;}
			if (dist_2atoms_pbc(R1_H, R2_L1, boxlnew)<RHL1) {std::cout<<"Conflict found between seg# "<<i<<" and seg# "<<j<<" during volume change"<<std::endl; return true;}
			if (dist_2atoms_pbc(R1_H, R2_L2, boxlnew)<RHL2) {std::cout<<"Conflict found between seg# "<<i<<" and seg# "<<j<<" during volume change"<<std::endl; return true;}
			if (dist_2atoms_pbc(R1_L1, R2_H, boxlnew)<RHL1) {std::cout<<"Conflict found between seg# "<<i<<" and seg# "<<j<<" during volume change"<<std::endl; return true;}
			if (dist_2atoms_pbc(R1_L1, R2_L1, boxlnew)<RL1L1) {std::cout<<"Conflict found between seg# "<<i<<" and seg# "<<j<<" during volume change"<<std::endl; return true;}
			if (dist_2atoms_pbc(R1_L1, R2_L2, boxlnew)<RL1L2) {std::cout<<"Conflict found between seg# "<<i<<" and seg# "<<j<<" during volume change"<<std::endl; return true;}
			if (dist_2atoms_pbc(R1_L2, R2_H, boxlnew)<RHL2) {std::cout<<"Conflict found between seg# "<<i<<" and seg# "<<j<<" during volume change"<<std::endl; return true;}
			if (dist_2atoms_pbc(R1_L2, R2_L1, boxlnew)<RL1L2) {std::cout<<"Conflict found between seg# "<<i<<" and seg# "<<j<<" during volume change"<<std::endl; return true;}
			if (dist_2atoms_pbc(R1_L2, R2_L2, boxlnew)<RL2L2) {std::cout<<"Conflict found between seg# "<<i<<" and seg# "<<j<<" during volume change"<<std::endl; return true;}
        } // end of loop j over mol._natoms() 

    } // end of loop i over mol._natoms() 

    return false;
}


/********* main           ******************/

int main(){

    //srand (time(NULL));
    
    RunParameters par ;

    int frame = 0 ;

    // restart parameters
    std::string file_box_length_start= "box_length_"+std::to_string(par.Nstart)+".txt";
    std::string file_com_start = "triplets_com_"+std::to_string(par.Nstart)+".pdb";
    std::string file_triplets_start = "triplets_"+std::to_string(par.Nstart)+".pdb";
    std::string file_triplets_final = "triplets_"+std::to_string(par.Nstart+1)+".pdb";
    std::string file_com_final = "triplets_com_"+std::to_string(par.Nstart+1)+".pdb";
    std::string file_box_length_run= "box_length_"+std::to_string(par.Nstart+1)+".txt";
    std::string file_triplets_dcd_run = "triplets_"+std::to_string(par.Nstart+1)+".dcd";

    // set up sasmol reference object (for writing out coordinate)
    sasmol::SasMol mol;
    sasmol::SasMol mol_triplets;
    const int natoms = N1d*N1d*N1d;
    EigenArrayXXfRow coor_center(natoms, 3);
    EigenArrayXXfRow coor_H(natoms, 3);
    EigenArrayXXfRow coor_L1(natoms, 3);
    EigenArrayXXfRow coor_L2(natoms, 3);
    EigenArrayXXfRow u0_H(natoms, 3);
    EigenArrayXXfRow u1_H(natoms, 3);
    EigenArrayXXfRow u2_H(natoms, 3);
    EigenArrayXXfRow u0_L1(natoms, 3);
    EigenArrayXXfRow u1_L1(natoms, 3);
    EigenArrayXXfRow u2_L1(natoms, 3);
    EigenArrayXXfRow u0_L2(natoms, 3);
    EigenArrayXXfRow u1_L2(natoms, 3);
    EigenArrayXXfRow u2_L2(natoms, 3);

    // set up the mol_triplets coordinates
    if (! par.flag_restart) // start from scratch
    {
        sasmol::SasMol mol_single;
        util::pp(">>> reading pdb") ;
        mol_single.read_pdb("triplet_H.pdb");
        mol_single.duplicate_molecule(mol, frame, N1d, par.box_length/N1d);
        mol.center(frame);

        sasmol::SasMol mol_triplet;
        mol_triplet.read_pdb("triplet.pdb");
        mol_triplet.duplicate_molecule(mol_triplets, frame, N1d, par.box_length/N1d);

        // initialize the heavy and light chain center coordinates
        initHLCoor(mol, coor_center, coor_H, coor_L1, coor_L2, u0_H, u1_H, u2_H, u0_L1, u1_L1, u2_L1, u0_L2, u1_L2, u2_L2);

        // update the coordinates in mol_triplets
        for (int i=0; i<mol_triplets._natoms()/3; ++i)
        {
            mol_triplets._x()(i*3,frame) = coor_H(i,0);
            mol_triplets._y()(i*3,frame) = coor_H(i,1);
            mol_triplets._z()(i*3,frame) = coor_H(i,2);
            mol_triplets._x()(i*3+1,frame) = coor_L1(i,0);
            mol_triplets._y()(i*3+1,frame) = coor_L1(i,1);
            mol_triplets._z()(i*3+1,frame) = coor_L1(i,2);
            mol_triplets._x()(i*3+2,frame) = coor_L2(i,0);
            mol_triplets._y()(i*3+2,frame) = coor_L2(i,1);
            mol_triplets._z()(i*3+2,frame) = coor_L2(i,2);
        }
        // save the triplets pdb file for the starting structure
        mol_triplets.write_pdb("triplets.pdb",frame);
    }
    else
    {
        mol.read_pdb(file_com_start);
        mol_triplets.read_pdb(file_triplets_start);
        // set the coordinates in coor_H/L12
        for (int i=0; i<mol_triplets._natoms()/3; ++i)
        {
            coor_center(i,0) = mol._x()(i,frame);
            coor_center(i,1) = mol._y()(i,frame);
            coor_center(i,2) = mol._z()(i,frame);
            coor_H(i,0) = mol_triplets._x()(i*3,frame);
            coor_H(i,1) = mol_triplets._y()(i*3,frame);
            coor_H(i,2) = mol_triplets._z()(i*3,frame);
            coor_L1(i,0) = mol_triplets._x()(i*3+1,frame);
            coor_L1(i,1) = mol_triplets._y()(i*3+1,frame);
            coor_L1(i,2) = mol_triplets._z()(i*3+1,frame);
            coor_L2(i,0) = mol_triplets._x()(i*3+2,frame);
            coor_L2(i,1) = mol_triplets._y()(i*3+2,frame);
            coor_L2(i,2) = mol_triplets._z()(i*3+2,frame);
            /// ZHL to do: read the ANISO
        }
        std::ifstream fin(file_box_length_start);
        std::string dump_string;
        std::stringstream ss;
        while (getline(fin,dump_string)) ss.str(dump_string);
        int dump_int;
        ss>>dump_int>>par.box_length>>par.density>>par.pressure>>par.delta_move[0]>>par.delta_move[1]>>par.delta_move[2]>>par.delta_boxlength;
        fin.close();
    }
    par.boxlinv = 1.0/par.box_length ;
    par.density = mol_triplets._natoms()/pow(par.box_length,3.0) ;
    float beta = 1.0/par.temperature ;         // i.e. kb=1.0

    std::cout << "number of atoms = " << mol._natoms() << std::endl ;
    std::cout << "starting box length = " << par.box_length << std::endl ;
    std::cout << "starting density = " << par.density << std::endl ;

    FILE *dcdoutfile = sasio::open_write_dcd_file(file_triplets_dcd_run, mol_triplets._natoms(), par.number_of_steps);

    int attempted_move[3] = {0, 0, 0}; // attempted moves for translation, rotation, and bending
    int accepted_move[3] = {0, 0, 0}; //acceptance rate for translation, rotation, and bending
    int attempted_boxlength = 0;
    int accepted_boxlength = 0;
    float acp = 0.0, acd = 0.0 ;

    std::cout << "density = " << par.density << std::endl ;
    std::cout << "initial pressure = " << par.density * par.temperature << std::endl ;

    int fr = 1 ;

///// MAIN LOOP

    int frames_since_last_dcd_save = 0 ;

    float rxiold, ryiold, rziold ;
    float rxinew, ryinew, rzinew ;
    EigenArray13fRow Riold_H, Riold_L1, Riold_L2, Riold_center; 
    EigenArray13fRow Rinew_H, Rinew_L1, Rinew_L2, Rinew_center; 
    EigenArray13fRow u0iold_H, u1iold_H, u2iold_H;
    EigenArray13fRow u0iold_L1, u1iold_L1, u2iold_L1;
    EigenArray13fRow u0iold_L2, u1iold_L2, u2iold_L2;
    EigenArray13fRow u0inew_H, u1inew_H, u2inew_H;
    EigenArray13fRow u0inew_L1, u1inew_L1, u2inew_L1;
    EigenArray13fRow u0inew_L2, u1inew_L2, u2inew_L2;

    float deltv, deltvb;
    
    float ratbox,dpv,dvol,delthb,rrbox,boxlnew ;

    int count = 0 ;
    int this_step ;
    int overall_frame = 0 ;
    float ratio, bratio ;

    float d_theta;
    int rot_axis;
    Eigen::Vector3f v_axis;

    std::ofstream fout(file_box_length_run,std::ofstream::out);
    std::cout << "STARTING MC SIMULATION : " << par.number_of_steps << " steps " << std::endl ;
    fout << "\nSTEP\t\tbox_length\t\tdensity\t\tpressure\t\tdelta_translation\t\tdelta_rotation\t\tdelta_bending\t\tdelta_boxlength"<<std::endl;

    struct timeval t0_overall, t1_overall, t2_overall;
    struct timeval t1_detail, t2_detail;


    for(this_step = 0 ; this_step < par.number_of_steps ; ++this_step)
    {
        std::cout<<"================================================================================"<<std::endl<<"Starting step: "<<this_step<<std::endl<<std::endl;
        if (Q_profiling_overall) gettimeofday(&t0_overall,NULL);

        if(this_step == 0)
        { 
            frames_since_last_dcd_save = par.dcd_save_frequency ;
        }
        else
        {
            frames_since_last_dcd_save += 1 ;    
        }    

        if (Q_profiling_overall) gettimeofday(&t1_overall,NULL);
        for(int i = 0 ; i < mol._natoms() ; ++i)
        {
            count++ ;
    
            if (Q_profiling_detail) gettimeofday(&t1_detail,NULL);
            int dice = int(util::get_random_float(0.0,3.0));
            if (dice<0||dice>2) dice = 0; //ZHL hardwired for a bug
            attempted_move[dice] ++;

            // get the temporary coordinates and u0/1/2
            Riold_center = coor_center.row(i);
            Riold_H = coor_H.row(i);
            Riold_L1 = coor_L1.row(i);
            Riold_L2 = coor_L2.row(i);
            Rinew_center = coor_center.row(i);
            Rinew_H = coor_H.row(i)-Rinew_center;
            Rinew_L1 = coor_L1.row(i)-Rinew_center;
            Rinew_L2 = coor_L2.row(i)-Rinew_center;
            u0iold_H = u0inew_H = u0_H.row(i);
            u1iold_H = u1inew_H = u1_H.row(i);
            u2iold_H = u2inew_H = u2_H.row(i);
            u0iold_L1 = u0inew_L1 = u0_L1.row(i);
            u1iold_L1 = u1inew_L1 = u1_L1.row(i);
            u2iold_L1 = u2inew_L1 = u2_L1.row(i);
            u0iold_L2 = u0inew_L2 = u0_L2.row(i);
            u1iold_L2 = u1inew_L2 = u1_L2.row(i);
            u2iold_L2 = u2inew_L2 = u2_L2.row(i);

            // translate the triplet
            if (dice == 0)
            {
                rxiold = mol._x()(i,frame) ;        
                ryiold = mol._y()(i,frame) ;        
                rziold = mol._z()(i,frame) ;        

                rxinew = rxiold + (util::get_random_float(-1.0,1.0)) * par.delta_move[dice] ;
                ryinew = ryiold + (util::get_random_float(-1.0,1.0)) * par.delta_move[dice] ;
                rzinew = rziold + (util::get_random_float(-1.0,1.0)) * par.delta_move[dice] ;

                rxinew = rxinew-par.box_length*(round(rxinew*par.boxlinv)) ;
                ryinew = ryinew-par.box_length*(round(ryinew*par.boxlinv)) ;
                rzinew = rzinew-par.box_length*(round(rzinew*par.boxlinv)) ;
                
                Rinew_center << rxinew, ryinew, rzinew;
                Rinew_H = Rinew_center + Rinew_H;
                Rinew_L1 = Rinew_center + Rinew_L1;
                Rinew_L2 = Rinew_center + Rinew_L2;
                if (Q_profiling_detail) util::profiling(t1_detail,t2_detail,"shift triplet");
            }
            // rotate the triplet
            else if (dice==1)
            {
                rxinew = mol._x()(i,frame) ;        
                ryinew = mol._y()(i,frame) ;        
                rzinew = mol._z()(i,frame) ;        
                
                d_theta = (util::get_random_float(-1.0,1.0)) * par.delta_move[dice];
                rot_axis = int(util::get_random_float(0.0,3.0));
                if (rot_axis==0) v_axis<<1.0,0.0,0.0;
                else if (rot_axis==1) v_axis<<0.0,1.0,0.0;
                else v_axis<<0.0,0.0,1.0;
                Eigen::AngleAxisf rot(d_theta, v_axis);

                Rinew_H = Rinew_center + (rot*(Rinew_H.matrix().transpose())).transpose().array();
                Rinew_L1 = Rinew_center + (rot*(Rinew_L1.matrix().transpose())).transpose().array();
                Rinew_L2 = Rinew_center + (rot*(Rinew_L2.matrix().transpose())).transpose().array();
                u0inew_H = (rot*(u0inew_H.matrix().transpose())).transpose().array();
                u1inew_H = (rot*(u1inew_H.matrix().transpose())).transpose().array();
                u2inew_H = (rot*(u2inew_H.matrix().transpose())).transpose().array();
                u0inew_L1 = (rot*(u0inew_L1.matrix().transpose())).transpose().array();
                u1inew_L1 = (rot*(u1inew_L1.matrix().transpose())).transpose().array();
                u2inew_L1 = (rot*(u2inew_L1.matrix().transpose())).transpose().array();
                u0inew_L2 = (rot*(u0inew_L2.matrix().transpose())).transpose().array();
                u1inew_L2 = (rot*(u1inew_L2.matrix().transpose())).transpose().array();
                u2inew_L2 = (rot*(u2inew_L2.matrix().transpose())).transpose().array();
                if (Q_profiling_detail) util::profiling(t1_detail,t2_detail,"rotate triplet");
            }
            // deform the triplet
            else
            {
                rxinew = mol._x()(i,frame) ;        
                ryinew = mol._y()(i,frame) ;        
                rzinew = mol._z()(i,frame) ;        
                
                v_axis = (Rinew_H.matrix().transpose()).cross((Rinew_L1.matrix().transpose()));
                v_axis = v_axis/v_axis.norm();
                d_theta = (util::get_random_float(-1.0,1.0)) * par.delta_move[dice];
                Eigen::AngleAxisf rot(d_theta, v_axis);

                rot_axis = int(util::get_random_float(0.0,3.0));
                if (rot_axis==0)
                {
                    Rinew_H = Rinew_center + (rot*(Rinew_H.matrix().transpose())).transpose().array();
                    Rinew_L1 = Rinew_center + Rinew_L1;
                    Rinew_L2 = Rinew_center + Rinew_L2;
                    u0inew_H = (rot*(u0inew_H.matrix().transpose())).transpose().array();
                    u1inew_H = (rot*(u1inew_H.matrix().transpose())).transpose().array();
                    u2inew_H = (rot*(u2inew_H.matrix().transpose())).transpose().array();
                }
                else if (rot_axis==1)
                {
                    Rinew_H = Rinew_center + Rinew_H;
                    Rinew_L1 = Rinew_center + (rot*(Rinew_L1.matrix().transpose())).transpose().array();
                    Rinew_L2 = Rinew_center + Rinew_L2;
                    u0inew_L1 = (rot*(u0inew_L1.matrix().transpose())).transpose().array();
                    u1inew_L1 = (rot*(u1inew_L1.matrix().transpose())).transpose().array();
                    u2inew_L1 = (rot*(u2inew_L1.matrix().transpose())).transpose().array();
                }
                else
                {
                    Rinew_H = Rinew_center + Rinew_H;
                    Rinew_L1 = Rinew_center + Rinew_L1;
                    Rinew_L2 = Rinew_center + (rot*(Rinew_L2.matrix().transpose())).transpose().array();
                    u0inew_L2 = (rot*(u0inew_L2.matrix().transpose())).transpose().array();
                    u1inew_L2 = (rot*(u1inew_L2.matrix().transpose())).transpose().array();
                    u2inew_L2 = (rot*(u2inew_L2.matrix().transpose())).transpose().array();
                }

                if (Q_profiling_detail) util::profiling(t1_detail,t2_detail,"deform triplet");
            }

            deltv = lj_energy_ellipsoid(mol._natoms(), Rinew_H,Rinew_L1,Rinew_L2, u0inew_H, u1inew_H, u2inew_H, u0inew_L1, u1inew_L1, u2inew_L1, u0inew_L2, u1inew_L2, u2inew_L2, i,par.box_length, coor_H, coor_L1, coor_L2, u0_H, u1_H, u2_H, u0_L1, u1_L1, u2_L1, u0_L2, u1_L2, u2_L2)
                 - lj_energy_ellipsoid(mol._natoms(), Riold_H,Riold_L1,Riold_L2, u0iold_H, u1iold_H, u2iold_H, u0iold_L1, u1iold_L1, u2iold_L1, u0iold_L2, u1iold_L2, u2iold_L2, i,par.box_length, coor_H, coor_L1, coor_L2, u0_H, u1_H, u2_H, u0_L1, u1_L1, u2_L1, u0_L2, u1_L2, u2_L2);

			deltvb = beta * deltv ;

			if(deltvb < 75.0)
			{
				if ( deltv <= 0.0 || (exp(-deltvb) > util::get_random_float(0.0, 1.0)) )
                {
                    mol._x()(i,frame) = rxinew ;
                    mol._y()(i,frame) = ryinew ;
                    mol._z()(i,frame) = rzinew ;
                    //util::profiling(t1,t2,"assign new sasmol");
                    coor_center.row(i) <<rxinew, ryinew, rzinew;
                    coor_H.row(i) = Rinew_H;
                    coor_L1.row(i) = Rinew_L1;
                    coor_L2.row(i) = Rinew_L2;
                    u0_H.row(i) = u0inew_H;
                    u1_H.row(i) = u1inew_H;
                    u2_H.row(i) = u2inew_H;
                    u0_L1.row(i) = u0inew_L1;
                    u1_L1.row(i) = u1inew_L1;
                    u2_L1.row(i) = u2inew_L1;
                    u0_L2.row(i) = u0inew_L2;
                    u1_L2.row(i) = u1inew_L2;
                    u2_L2.row(i) = u2inew_L2;
                    //util::profiling(t1,t2,"assign new coor");
                    accepted_move[dice] += 1.0 ;
                }
            }

            par.pressure = par.density * par.temperature;

            acp = acp + par.pressure ;
            acd = acd + par.density ;
            
        } // end of loop i over mol._natoms()

        if(Q_profiling_overall) util::profiling(t1_overall,t2_overall,"update atom moves");
            
        // change the volume

        const float rand_volchange = util::get_random_float(-1.0,1.0);
        if (rand_volchange<0.0) attempted_boxlength ++;
        std::cout<<"rand_volchange: "<<rand_volchange<<std::endl;
        boxlnew = par.box_length + rand_volchange * par.delta_boxlength ;
        
        ratbox = par.box_length / boxlnew ;
        rrbox = 1.0 / ratbox ;

        dpv = par.goal_pressure * (pow(boxlnew,3.0) - pow(par.box_length,3.0)) ;
        dvol = 3.0 * par.temperature * mol_triplets._natoms() * log(ratbox) ; //log == ln
        delthb = beta * ( dpv + dvol ) ;
        std::cout<<"ZHL delthb: "<<delthb<<" dpv: "<<dpv<<" dvol: "<<dvol<<" old boxlength: "<<par.box_length<<" new boxlength: "<<boxlnew<<" rand: "<<rand_volchange<<std::endl;

        bool flag_conflict;
        if (rand_volchange>=0.0) flag_conflict=false;
        else flag_conflict = check_conflict_scale(mol._natoms(), coor_center, coor_H, coor_L1, coor_L2, rrbox, boxlnew);
        if ( ! flag_conflict )
        {
            if(delthb < 75.0)
            {
                if(delthb <= 0.0)
                {
                    mol._x() *= rrbox ;
                    mol._y() *= rrbox ;
                    mol._z() *= rrbox ;

                    coor_H = coor_center*rrbox + (coor_H-coor_center);
                    coor_L1 = coor_center*rrbox + (coor_L1-coor_center);
                    coor_L2 = coor_center*rrbox + (coor_L2-coor_center);
                    coor_center *= rrbox;

                    par.box_length = boxlnew ;
                    if (rand_volchange<0.0) accepted_boxlength++;
                }
                else if(exp(-delthb) > util::get_random_float(0.0,1.0))
                {
                    mol._x() *= rrbox ;
                    mol._y() *= rrbox ;
                    mol._z() *= rrbox ;

                    coor_H = coor_center*rrbox + (coor_H-coor_center);
                    coor_L1 = coor_center*rrbox + (coor_L1-coor_center);
                    coor_L2 = coor_center*rrbox + (coor_L2-coor_center);
                    coor_center *= rrbox;

                    par.box_length = boxlnew ;
                    if (rand_volchange<0.0) accepted_boxlength++;
                }    
            }
        }
        if (Q_profiling_overall) util::profiling(t1_overall, t2_overall, "update volume");

        par.boxlinv = 1.0/par.box_length ;
        par.density = mol_triplets._natoms()/pow(par.box_length,3.0) ;

        par.pressure = par.density * par.temperature;

        acp = acp + par.pressure ;
        acd = acd + par.density ;

        if(frames_since_last_dcd_save == par.dcd_save_frequency || this_step+1 == par.number_of_steps)
        {
            overall_frame += 1 ;
            for (int i=0; i<mol_triplets._natoms()/3; ++i)
            {
                mol_triplets._x()(i*3,frame) = coor_H(i,0);
                mol_triplets._y()(i*3,frame) = coor_H(i,1);
                mol_triplets._z()(i*3,frame) = coor_H(i,2);
                mol_triplets._x()(i*3+1,frame) = coor_L1(i,0);
                mol_triplets._y()(i*3+1,frame) = coor_L1(i,1);
                mol_triplets._z()(i*3+1,frame) = coor_L1(i,2);
                mol_triplets._x()(i*3+2,frame) = coor_L2(i,0);
                mol_triplets._y()(i*3+2,frame) = coor_L2(i,1);
                mol_triplets._z()(i*3+2,frame) = coor_L2(i,2);
            }
            mol_triplets.write_dcd_step(dcdoutfile,frame,overall_frame) ;
            frames_since_last_dcd_save = 0 ;
            fout << this_step << " "<<par.box_length<<" "<<par.density <<" "<<par.pressure<<" "<<par.delta_move[0]<<" "<<par.delta_move[1]<<" "<<par.delta_move[2]<<" "<<par.delta_boxlength<<std::endl;
        }

        //if(fmodf(float(this_step),par.interval_move)==0)
        {
            std::cout<<std::endl;
            std::cout<<"**************************"<<std::endl;
            std::cout<<"STEP: "<<this_step<<std::endl;
            std::cout<<std::endl;
            std::cout<<"ACC_TRANSLATION: "<<accepted_move[0]*100.0/(attempted_move[0]*par.interval_move)<<" % ("<<accepted_move[0]<<" out of "<<attempted_move[0]*par.interval_move<<")"<<std::endl;
            std::cout<<"ACC_ROTATION:    "<<accepted_move[1]*100.0/(attempted_move[1]*par.interval_move)<<" % ("<<accepted_move[1]<<" out of "<<attempted_move[1]*par.interval_move<<")"<<std::endl;
            std::cout<<"ACC_BENDING:     "<<accepted_move[2]*100.0/(attempted_move[2]*par.interval_move)<<" % ("<<accepted_move[2]<<" out of "<<attempted_move[2]*par.interval_move<<")"<<std::endl;
            std::cout<<"ACC_BOXL:      "<<accepted_boxlength*100.0/attempted_boxlength<<" % ("<<accepted_boxlength<<" out of "<<attempted_boxlength<<")"<<std::endl;
            std::cout<<std::endl;
            std::cout<<"STEP_TRANSLATION: "<<par.delta_move[0]<<std::endl;
            std::cout<<"STEP_ROTATION:    "<<par.delta_move[1]*180./M_PI<<std::endl;
            std::cout<<"STEP_BENDING:     "<<par.delta_move[2]*180./M_PI<<std::endl;
            std::cout<<"STEP_BOXL:      "<<par.delta_boxlength<<std::endl;
            std::cout<<std::endl;
            std::cout<<"DENSITY:    "<<par.density<<std::endl;
            std::cout<<"PRESSURE:   "<<par.pressure<<std::endl;
            std::cout<<"BOX_LENGTH: "<<par.box_length<<std::endl;
            std::cout<<std::endl;
        }

    
        // update translation, rotation and bending move step size
        for (int i=0; i<3; ++i)
        {
            if (float(accepted_move[i])/(float(attempted_move[i])) > 0.5) par.delta_move[i] += (par.delta_move_max[i]-par.delta_move[i])*0.05 ;
            else par.delta_move[i] *= 0.95 ;
            if(fmodf(float(this_step+1),par.interval_move)==0)
            {
                    attempted_move[i] = 0;
                    accepted_move[i] = 0;
            }
        }

        // update volume step size
        if (rand_volchange<0.0)
        {
            if (float(accepted_boxlength)/float(attempted_boxlength) > 0.5) par.delta_boxlength *= 1.05;
            else par.delta_boxlength *= 0.95;
        }
        if(fmodf(float(this_step+1),par.interval_boxlength)==0)
        {
            accepted_boxlength = 0 ;
            attempted_boxlength = 0 ;
        }

        if (Q_profiling_overall) util::profiling(t0_overall, t2_overall, "MC step");
        std::cout<<"================================================================================"<<std::endl<<std::endl;

    } // end of loop step over number_of_steps

    std::cout << "ending box length = " << par.box_length << std::endl ;
    std::cout << "ending density = " << par.density << std::endl ;

    mol.write_pdb(file_com_final,frame) ;
    mol_triplets.write_pdb(file_triplets_final,frame) ;
    fout.close();
    close_dcd_write(dcdoutfile) ;

//    stop_here() ;

    util::print_run_details() ;
    std::cout << " >>> DONE <<< " << std::endl << std::endl ;
    return 0 ;


}


