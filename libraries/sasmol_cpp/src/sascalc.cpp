#include <sascalc.h>
#include <sasmol.h>
#include <map>
#include <Eigen/Dense>

void
initialize_amu(std::map<std::string,float> &amu)
{
        amu["H"]  = 1.007947 ;
        amu["He"] = 4.0026022 ;
        amu["Li"] = 6.9412 ; 
        amu["Be"] = 9.0121823 ; 
        amu["B"]  = 10.8117 ;
        amu["C"]  = 12.01078 ;
        amu["N"]  = 14.00672 ;
        amu["O"]  = 15.99943 ;
        amu["F"]  = 18.99840325 ;
        amu["Ne"] = 20.17976 ;
        amu["Na"] = 22.9897702 ; 
        amu["Mg"] = 24.30506 ;
        amu["Al"] = 26.9815382 ; 
        amu["Si"] = 28.08553 ;
        amu["P"]  = 30.9737612; 
        amu["S"]  = 32.0655 ; 
        amu["Cl"] = 35.4532 ; 
        amu["Ar"] = 39.9481 ;
        amu["K"]  = 39.09831 ; 
        amu["Ca"] = 40.0784 ; 
        amu["Sc"] = 44.9559108 ;
        amu["Ti"] = 47.8671 ;
        amu["V"]  = 50.94151 ; 
        amu["Cr"] = 51.99616 ; 
        amu["Mn"] = 54.9380499 ; 
        amu["Fe"] = 55.8452 ;
        amu["Co"] = 58.9332009 ; 
        amu["Ni"] = 58.69342 ; 
        amu["Cu"] = 63.5463 ; 
        amu["Zn"] = 65.4094 ;
        amu["Ga"] = 69.7231 ; 
        amu["Ge"] = 72.641 ; 
        amu["As"] = 74.921602 ; 
        amu["Se"] = 78.963 ;
        amu["Br"] = 79.9041 ; 
        amu["Kr"] = 83.7982 ; 
        amu["Rb"] = 85.46783 ;
        amu["Sr"] = 87.621 ;
        amu["Y"]  = 88.905852 ; 
        amu["Zr"] = 91.2242 ; 
        amu["Nb"] = 92.906382 ;
        amu["Mo"] = 95.942 ; 
        amu["Tc"] = 98.0 ; 
        amu["Ru"] = 101.072 ; 
        amu["Rh"] = 102.905502 ; 
        amu["Pd"] = 106.421 ; 
        amu["Ag"] = 107.86822 ; 
        amu["Cd"] = 112.4118 ; 
        amu["In"] = 114.8183 ; 
        amu["Sn"] = 118.7107 ;
        amu["Sb"] = 121.7601 ; 
        amu["Te"] = 127.603 ; 
        amu["I"]  = 126.904473 ; 
        amu["Xe"] = 131.2936 ;
        amu["Cs"] = 132.905452 ; 
        amu["Ba"] = 137.3277 ; 
        amu["La"] = 138.90552 ; 
        amu["Ce"] = 140.1161 ; 
        amu["Pr"] = 140.907652 ; 
        amu["Nd"] = 144.243 ; 
        amu["Pm"] = 145.0 ; 
        amu["Sm"] = 150.363 ; 
        amu["Eu"] = 151.9641 ; 
        amu["Gd"] = 157.253 ; 
        amu["Tb"] = 158.925342 ; 
        amu["Dy"] = 162.5001 ;
        amu["Ho"] = 164.930322 ; 
        amu["Er"] = 167.2593 ; 
        amu["Tm"] = 168.93421 ; 
        amu["Yb"] = 173.043 ; 
        amu["Lu"] = 174.9671 ; 
        amu["Hf"] = 174.9671 ; 
        amu["Ta"] = 180.94791 ; 
        amu["W"]  = 183.841 ; 
        amu["Re"] = 186.2071 ; 
        amu["Os"] = 190.233 ; 
        amu["Ir"] = 192.2173 ; 
        amu["Pt"] = 195.0782 ; 
        amu["Au"] = 196.966552 ; 
        amu["Hg"] = 200.592 ; 
        amu["Tl"] = 204.38332 ; 
        amu["Pb"] = 207.21 ; 
        amu["Bi"] = 208.980382 ; 
        amu["Po"] = 209.0 ; 
        amu["At"] = 210.0 ; 
        amu["Rn"] = 222.0 ; 
        amu["Fr"] = 223.0 ; 
        amu["Ra"] = 226.0 ; 
        amu["Ac"] = 227.0 ; 
        amu["Th"] = 232.03811 ; 
        amu["Pa"] = 231.035882 ; 
        amu["U"]  = 238.028913 ; 
        amu["D"]  = 2.01410177804 ; 
        amu["1H"] = 1.00782503214 ;
} ;

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::vector<float>
sascalc::Prop::
calc_com(int frame)
{

    if(_total_mass() <= 0.0) {
        calc_mass() ;
        std::cout << "com tot mass = " << _total_mass() << std::endl ;
    }

    float comx = 0.0, comy = 0.0, comz = 0.0 ;

    try
    {
        comx = (_atom_mass()*_x().col(frame)).sum() / _total_mass() ;    
        comy = (_atom_mass()*_y().col(frame)).sum() / _total_mass() ;    
        comz = (_atom_mass()*_z().col(frame)).sum() / _total_mass() ;    
    }
    catch(std::overflow_error e) {
        std::cout << "failed in calc_com" << std::endl ;
        std::cout << "total_mass = " << _total_mass() << std::endl ;
        std::cout << e.what() << " -> ";
		float tmp[]={comx, comy, comz};
        std::vector<float> com(tmp, tmp+3);
        return com ;
    }

	float tmp[]={comx, comy, comz};
    std::vector<float> com(tmp, tmp+3);

    return com ;
}


///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
float
sascalc::Prop::
calc_rg(int &frame){

    std::vector<float> com = calc_com(frame) ;

    float rgx = ((_x().col(frame) - com[0]).pow(2.0)).sum() ;
    float rgy = ((_y().col(frame) - com[1]).pow(2.0)).sum() ;
    float rgz = ((_z().col(frame) - com[2]).pow(2.0)).sum() ;
    
    float rg = sqrt((rgx + rgy + rgz) / _natoms()) ;    

    return rg ;
}


///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sascalc::Prop::
calc_mass()
{
    /*
        http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some

        standard atomic weight is based on the natural istopic composition

        NOTE: deuterium is 'D' 2H1 and '1H' is 1H1, all other elements 
        have their natural abundance weight. These elements are located
        at the end of the dictionary.

    */    
        
    std::map<std::string,float> standard_atomic_weight ;
    initialize_amu(standard_atomic_weight) ;

    float total_mass = 0.0 ;
    std::string this_element ;
    float this_mass ;

    _atom_mass().setZero(_natoms(), 1);

    for(int i = 0 ; i < _natoms() ; ++i)
    {
        this_element = _atom_selement()[i] ;
        this_mass = standard_atomic_weight.find(this_element)->second ;
        total_mass += this_mass ;
        _atom_mass()[i] = this_mass ;
    }

    _total_mass() = total_mass ;
    
    return ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
float
sascalc::Prop::
calc_volume()
{
    /*
     * 0.73 cm**3/g  x  10**24 A**3/cm**3  x  molecular weight g/mole
     *  --------------------------------------------------------------
     *               6.02 x 10**23 molecules/mole
     */

    std::string this_resname;
    float this_mass;
    const float scaling_protein = 0.73e24/6.02e23;
    const float scaling_na = 0.56e24/6.02e23;
    float tot_volume = 0.0;

    std::vector<std::string> protein_resnames={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HSD","HSE","HSP","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
   
    std::vector<std::string> dna_resnames={"NUSA","NUSG","NUSC","NUSU","DA","DG","DC","DT","ADE","GUA","CYT","THY"};

    std::vector<std::string> rna_resnames={"RNUS","RNUA","RUUG","RNUC","A", "C", "G", "U","ADE","CYT","GUA","URA"};

    std::vector<std::string> nucleic_resnames = {"GUA","ADE","CYT","THY","URA","G", "A", "C", "T", "U","DA","DG","DC","DT"};

    for(int i = 0 ; i < _natoms() ; ++i)
    {
        this_resname = _atom_resname()[i] ;
        this_mass = _atom_mass()[i] ;
        if (std::find(protein_resnames.begin(), protein_resnames.end(), this_resname) != protein_resnames.end())
        {
            tot_volume += scaling_protein*this_mass;
        }
        else if (std::find(dna_resnames.begin(), dna_resnames.end(), this_resname) != dna_resnames.end())
        {
            tot_volume += scaling_na*this_mass;
        }
        else if (std::find(rna_resnames.begin(), rna_resnames.end(), this_resname) != rna_resnames.end())
        {
            tot_volume += scaling_na*this_mass;
        }
        else if (std::find(nucleic_resnames.begin(), nucleic_resnames.end(), this_resname) != nucleic_resnames.end())
        {
            tot_volume += scaling_na*this_mass;
        }
        else
        {
            tot_volume += scaling_protein*this_mass;///< @note ZHL hack
        }
    }
    
    return tot_volume;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sascalc::Prop::
calc_pmi(int &frame)
{
    std::vector<float> com = calc_com(frame) ;

    _x().col(frame) -= com[0] ;
    _y().col(frame) -= com[1] ;
    _z().col(frame) -= com[2] ;

    float Ixx = (_atom_mass()*(_y().col(frame)*_y().col(frame) + _z().col(frame)*_z().col(frame))).sum() ;
    float Iyy = (_atom_mass()*(_x().col(frame)*_x().col(frame) + _z().col(frame)*_z().col(frame))).sum() ;
    float Izz = (_atom_mass()*(_x().col(frame)*_x().col(frame) + _y().col(frame)*_y().col(frame))).sum() ;

    float Ixy = (-_atom_mass()*(_x().col(frame)*_y().col(frame))).sum() ;
    float Ixz = (-_atom_mass()*(_x().col(frame)*_z().col(frame))).sum() ;
    float Iyz = (-_atom_mass()*(_y().col(frame)*_z().col(frame))).sum() ;
    float Iyx = (-_atom_mass()*(_y().col(frame)*_x().col(frame))).sum() ;
    float Izx = (-_atom_mass()*(_z().col(frame)*_x().col(frame))).sum() ;
    float Izy = (-_atom_mass()*(_z().col(frame)*_y().col(frame))).sum() ;
    
    Eigen::Matrix3f I ;

    I << Ixx, Ixy, Ixz,
         Iyx, Iyy, Iyz,
         Izx, Izy, Izz;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(I);

    if (eigensolver.info() != Eigen::Success) 
    {
        std::cout << "PMI calculation failed" << std::endl ;
        return ;
    }

    _uk() = eigensolver.eigenvalues() ; // .col(0) ;
    _ak() = eigensolver.eigenvectors() ;

//    std::cout << "The eigenvalues of I are:\n" << eigensolver.eigenvalues() << std::endl;
//    std::cout << "Here's a matrix whose columns are eigenvectors of I \n"
//    << "corresponding to these eigenvalues:\n"
//    << eigensolver.eigenvectors() << std::endl;

    _x().col(frame) += com[0] ;
    _y().col(frame) += com[1] ;
    _z().col(frame) += com[2] ;

    return ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
int
sascalc::Prop::
self_overlap(float &cutoff, int &frame){

    int check = 0 ;     // overlapped == 1
                  // default == no overlap will return == 0 

    float dist ; 
    float this_x , that_x ;
    float this_y , that_y ;
    float this_z , that_z ;
    float dx2, dy2, dz2 ;
    int count = 1 ;

    for(int i = 0 ; i != _natoms()-1 ; ++i)
    {
        this_x = _x()(i,frame) ;
        this_y = _y()(i,frame) ;
        this_z = _z()(i,frame) ;

        for(int j = i+1 ; j != _natoms() ; ++j)
        {
            that_x = _x()(j,frame) ;
            that_y = _y()(j,frame) ;
            that_z = _z()(j,frame) ;

            dx2 = (this_x - that_x)*(this_x - that_x) ;
            dy2 = (this_y - that_y)*(this_y - that_y) ;
            dz2 = (this_z - that_z)*(this_z - that_z) ;

            dist = sqrt(dx2+dy2+dz2) ;

            if(dist < cutoff) 
            {
                std::cout << "count = " << count << std::endl ;
                std::cout << "this_atom = " << _atom_name()[count-1] << std::endl ;
                std::cout << "that_atom = " << _atom_name()[count] << std::endl ;
                std::cout << "dist = " << dist << std::endl ;
                std::cout << "this_x = " << this_x << std::endl ;
                std::cout << "this_y = " << this_y << std::endl ;
                std::cout << "this_z = " << this_z << std::endl ;
                std::cout << "that_x = " << that_x << std::endl ;
                std::cout << "that_y = " << that_y << std::endl ;
                std::cout << "that_z = " << that_z << std::endl ;
                check = 1 ; 
                return check ;
            }
        }    
//        std::cout << count << " " ;
        count++ ;

    }
    std::cout << std::endl ;
    
    return check ;

}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
float
sascalc::Prop::
calc_rmsd(sasmol::SasMol &mol, int &frame){

    float dx = (_x().col(frame)-mol._x().col(frame)).pow(2.0).sum() ;
    float dy = (_y().col(frame)-mol._y().col(frame)).pow(2.0).sum() ;
    float dz = (_z().col(frame)-mol._z().col(frame)).pow(2.0).sum() ;

    float rmsd = sqrt((dx + dy + dz)/_natoms()) ;

    return rmsd ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::vector<std::vector<float > >
sascalc::Prop::
calc_minmax(int &frame){

    float min_x = _x().col(frame).minCoeff() ; float max_x = _x().col(frame).maxCoeff() ;
    float min_y = _y().col(frame).minCoeff() ; float max_y = _y().col(frame).maxCoeff() ;
    float min_z = _z().col(frame).minCoeff() ; float max_z = _z().col(frame).maxCoeff() ;

    std::vector<float> min_vector ;
    std::vector<float> max_vector ;

    min_vector.push_back(min_x) ; min_vector.push_back(min_y) ; min_vector.push_back(min_z) ;
    max_vector.push_back(max_x) ; max_vector.push_back(max_y) ; max_vector.push_back(max_z) ;

    std::vector<std::vector<float > > minmax ;
    minmax.push_back(min_vector) ; minmax.push_back(max_vector) ;

    return minmax ;
}

/*
to do:
    def calc_minmax_all_steps() --> apbs.py & sld.py
    def calculate_residue_charge() --> mmc only
*/


///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sascalc::Prop::
find_surf(int * const idx)
{
	const double r = 1.8; // parameter to derive almn and blmn
	const double d = 2.5; // thickness of the surface layer
	const double rou = -15; // molecular 1 interior parameter
	const double eta = 1.0; // grid step size, 1.0-1.2
    const int N = 90; // Number of grid points
	const bool qsurf = 1; // flag for surf atoms to be passed to construct_lmn
	const double percsurf = 0.2; // percentage of coverage to decide a surf atom in construct_lmn
    double *lmn = (double*)malloc(N*N*N*sizeof(int));
    int i, ix, iy, iz, ix_tmp, iy_tmp, iz_tmp;
    int iX_low, iY_low, iZ_low, iX_high, iY_high, iZ_high;
    double x,y,z,tx,ty,tz;
    double dis;
    for (i=0; i<_natoms(); ++i)
    {
        x=_x()(i);
        y=_y()(i);
        z=_z()(i);
        iX_low = (int)(floor((x-r-d)/eta)) + N/2;
        iY_low = (int)(floor((y-r-d)/eta)) + N/2;
        iZ_low = (int)(floor((z-r-d)/eta)) + N/2;
        iX_high = (int)(ceil((x+r+d)/eta)) + N/2;
        iY_high = (int)(ceil((y+r+d)/eta)) + N/2;
        iZ_high = (int)(ceil((z+r+d)/eta)) + N/2;
        for (ix=iX_low; ix<=iX_high; ix++)
        {
            ix_tmp = ix;
            while(ix_tmp<0) {ix_tmp+=N;}
            while(ix_tmp>=N) {ix_tmp-=N;}
            tx = (ix-N/2)*eta;
            for (iy=iY_low; iy<=iY_high; iy++)
            {
                iy_tmp = iy;
                ty = (iy-N/2)*eta;
                while(iy_tmp<0) {iy_tmp+=N;}
                while(iy_tmp>=N) {iy_tmp-=N;}
                for (iz=iZ_low; iz<=iZ_high; iz++)
                {
                    iz_tmp = iz;
                    tz = (iz-N/2)*eta;
                    while(iz_tmp<0) {iz_tmp+=N;}
                    while(iz_tmp>=N) {iz_tmp-=N;}
                    dis = sqrt(pow((tx-x),2.0)+pow((ty-y),2.0)+pow((tz-z),2.0));
                    if (dis<=(r+d))
                    {
                        if (dis>r)
                        {
                            if (lmn[ix_tmp*N*N + iy_tmp*N + iz_tmp]==0)
                            {
                                lmn[ix_tmp*N*N + iy_tmp*N + iz_tmp]=1;
                            }
                        }
                        else
                        {
                            lmn[ix_tmp*N*N + iy_tmp*N + iz_tmp] = rou;
                        }
                    }

                }
            }
        }
    }

    // identify the surface atoms
    double val;
    double perc,rp=12;
    int count,count_s;
    double epsilon=1.0e-8;
    for (i=0; i<_natoms(); ++i)
    {
        x=_x()(i);
        y=_y()(i);
        z=_z()(i);
        iX_low = (int)(floor((x-r-d)/eta)) + N/2;
        iY_low = (int)(floor((y-r-d)/eta)) + N/2;
        iZ_low = (int)(floor((z-r-d)/eta)) + N/2;
        iX_high = (int)(ceil((x+r+d)/eta)) + N/2;
        iY_high = (int)(ceil((y+r+d)/eta)) + N/2;
        iZ_high = (int)(ceil((z+r+d)/eta)) + N/2;
        count=0;
        count_s=0;
        for (ix=iX_low; ix<=iX_high; ix++)
        {
            ix_tmp = ix;
            while(ix_tmp<0) {ix_tmp+=N;}
            while(ix_tmp>=N) {ix_tmp-=N;}
            tx = (ix-N/2)*eta;
            for (iy=iY_low; iy<=iY_high; iy++)
            {
                iy_tmp = iy;
                ty = (iy-N/2)*eta;
                while(iy_tmp<0) {iy_tmp+=N;}
                while(iy_tmp>=N) {iy_tmp-=N;}
                for (iz=iZ_low; iz<=iZ_high; iz++)
                {
                    iz_tmp = iz;
                    tz = (iz-N/2)*eta;
                    while(iz_tmp<0) {iz_tmp+=N;}
                    while(iz_tmp>=N) {iz_tmp-=N;}
                    dis = sqrt(pow((tx-x),2.0)+pow((ty-y),2.0)+pow((tz-z),2.0));
                    if (dis>r && dis<=(r+d))
                    {
                        val = lmn[ix_tmp*N*N + iy_tmp*N + iz_tmp];
                        if (abs(val-1.)<epsilon)
                        {
                            count_s++;
                        }
                        //if (abs(val-0)<epsilon) cout<<"WRONG"<<endl;
                        if (abs(val-0)<epsilon) printf("aWRONG\n");
                        count++;
                    }
                }
            }
        }
        perc=(double)(count_s)/(double)(count);
        //cout<<"ATOM "<<i<<" portion covered: "<<perc<<endl;
        if (perc>percsurf)
        {
            //cout<<i<<endl;
            idx[i]=1;
        }
        else
        {
            idx[i]=0;
        }
    }

    free(lmn);
}


///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sascalc::Prop::
find_non_aliphatic_Hs(int * const idx)
{
    char previous_nonH_atom_short_name;
    for (int i=0; i<_natoms(); ++i)
    {
        const std::string & atom_full_name = _atom_name()[i];
        if (atom_full_name[0] == 'H') // H's
        {
            if (previous_nonH_atom_short_name!='C') idx[i] = 1;
            else idx[i] = 0;
        }
        else // non-H's
        {
            idx[i] = 0;
            previous_nonH_atom_short_name = atom_full_name[0];
        }
    }
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sascalc::Prop::
find_backbone_Hs(int * const idx)
{
    char previous_nonH_atom_short_name;
    for (int i=0; i<_natoms(); ++i)
    {
        const std::string & atom_full_name = _atom_name()[i];
        if (atom_full_name[0] == 'H' && atom_full_name[1] == 'N') idx[i] = 1;
        else idx[i] = 0;
    }
}
