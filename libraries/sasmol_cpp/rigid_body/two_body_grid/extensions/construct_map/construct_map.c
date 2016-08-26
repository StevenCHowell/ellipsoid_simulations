#include <math.h>
#include <Python.h>
#include <stdio.h>
#include "numpy/arrayobject.h"


PyObject *func_construct_map(PyObject *self, PyObject *args)
{
    PyObject * input;
    PyArrayObject * py_coor;
    double rou;
    int N;
    double r;
    double d;
    double eta;
    int qsurf;
    double percsurf;
    if (!PyArg_ParseTuple(args, "OOdiddd|id", &pdb, &input, &rou, &N, &r, &d, &eta, &qsurf, &percsurf))
		return NULL;
    py_coor = (PyArrayObject*)PyArray_ContiguousFromObject(input,PyArray_DOUBLE, 2,2);
    int dim[3]={N,N,N};
    PyArrayObject * py_lmn = (PyArrayObject *)PyArray_FromDims(3,dim,PyArray_DOUBLE);
    double *lmn = (double*)(py_lmn->data);
    double *coor = (double*)(py_coor->data);
    int i, ix, iy, iz, ix_tmp, iy_tmp, iz_tmp;
    int iX_low, iY_low, iZ_low, iX_high, iY_high, iZ_high;
    double x,y,z,tx,ty,tz;
    double dis;
    for (i=0; i<py_coor->dimensions[0]; i++)
    {
        x=coor[i*3];
        y=coor[i*3+1];
        z=coor[i*3+2];
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

    if (!qsurf)
    {
        return PyArray_Return(py_lmn);
    }
    else
    // identify the surface atoms
    {
        double val;
        double perc,rp=12;
        int count,count_s;
        double epsilon=1.0e-8;
        int Natoms = ((int*)(py_coor->dimensions))[0];
        int dim1[1]={Natoms};
        PyArrayObject * py_idx = (PyArrayObject *)PyArray_FromDims(1,dim1,PyArray_INT);
        int *idx = (int*)(py_idx->data);
        for (i=0; i<py_coor->dimensions[0]; i++)
        {
            //cout<<endl<<"ATOM "<<i<<endl;
            x=coor[i*3];
            y=coor[i*3+1];
            z=coor[i*3+2];
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
        return PyArray_Return(py_idx);
    }
}


static PyMethodDef exampleMethods[] = 
{
	{ "construct_lmn", func_construct_map, METH_VARARGS },
	{ NULL, NULL }
} ;

PyMODINIT_FUNC initconstruct()
{
	PyObject *m ;
    import_array();
	m = Py_InitModule("construct", exampleMethods);
}
