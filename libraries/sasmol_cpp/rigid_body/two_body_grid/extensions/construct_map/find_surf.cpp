#include <math.h>
#include <stdio.h>


// find the surface atoms for the solute
void find_surf(int * const idx)
{
	const double r = 1.8; // parameter to derive almn and blmn
	const double d = 2.5; // thickness of the surface layer
	const double rou = -15; // molecular 1 interior parameter
	const double eta = 1.0; // grid step size, 1.0-1.2
    const int N = 90; // Number of grid points
	const bool qsurf = 1; // flag for surf atoms to be passed to construct_lmn
	const double percsurf = 0.3; // percentage of coverage to decide a surf atom in construct_lmn
    double *lmn = (double*)malloc(N*N*N*sizeof(int));
    int i, ix, iy, iz, ix_tmp, iy_tmp, iz_tmp;
    int iX_low, iY_low, iZ_low, iX_high, iY_high, iZ_high;
    double x,y,z,tx,ty,tz;
    double dis;
    for (i=0; i<_natoms(); ++i)
    {
        x=_x()[i];
        y=_y()[i];
        z=_z()[i];
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
        x=_x()[i];
        y=_y()[i];
        z=_z()[i];
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
                        if (abs(val-0)<epsilon) printf("WRONG\n");
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
