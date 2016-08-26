#include <Python.h>
#include "numpy/arrayobject.h"


PyObject *func_calc_threadholds(PyObject *self, PyObject *args)
{
    PyObject * PO_clmn, * PO_threadholds;
    PyArrayObject * py_clmn, *py_threadholds;
    if (!PyArg_ParseTuple(args, "OO", &PO_clmn, &PO_threadholds))
		return NULL;
    py_clmn = (PyArrayObject*)PyArray_ContiguousFromObject(PO_clmn,PyArray_DOUBLE,1,1);
    py_threadholds = (PyArrayObject*)PyArray_ContiguousFromObject(PO_threadholds,PyArray_DOUBLE,1,1);
    double *clmn = (double*)(py_clmn->data);
    double *threadholds = (double*)(py_threadholds->data);

    int N_threadholds = py_threadholds->dimensions[0];
    int dim[1]={N_threadholds};
    PyArrayObject * py_tot_nf_threadholds = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_DOUBLE);
    double * tot_nf_threadholds = (double*)(py_tot_nf_threadholds->data);

    int i,j;
    int flag_threadholds[N_threadholds];
    for (i=0; i<N_threadholds; i++)
    {
        flag_threadholds[i]=0;
    }
    double tot_cl = 0.0;
    for (i=0; i<py_clmn->dimensions[0]; i++)
    {
        tot_cl += clmn[i];
    }
    //cout<<"tot cl "<<tot_cl<<endl;
    double ratio;
    double sum_cl = 0.0;
    for (i=0; i<py_clmn->dimensions[0]; i++)
    {
        sum_cl += clmn[i];
        //cout<<"sum cl "<<sum_cl<<endl;
        ratio = sum_cl/tot_cl;
        //cout<<"ratio "<<ratio<<endl;
        for (j=0; j<N_threadholds; j++)
        {
            if ( (flag_threadholds[j]==0) && (ratio>threadholds[j]) )
            {
                //cout<<"Before "<<tot_nf_threadholds[j]<<endl;
                tot_nf_threadholds[j] += i;
                //cout<<"After "<<tot_nf_threadholds[j]<<endl;
                flag_threadholds[j]=1;
            }
        }
    }
    Py_DECREF(py_clmn);
    Py_DECREF(py_threadholds);
    return PyArray_Return(py_tot_nf_threadholds);
}


static PyMethodDef exampleMethods[] = 
{
	{ "calc_threadholds", func_calc_threadholds, METH_VARARGS },
	{ NULL, NULL }
} ;

PyMODINIT_FUNC initcalc_threadholds()
{
	Py_InitModule("calc_threadholds", exampleMethods);
    import_array();
}


