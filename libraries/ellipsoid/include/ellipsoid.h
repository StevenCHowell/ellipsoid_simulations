#ifndef H_ELLIPSOID
#define H_ELLIPSOID

#include <vector>
#include <math.h>
#include <Eigen/Dense>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_matrix_float.h>
#include <gsl/gsl_sort_vector_float.h>
#include <gsl/gsl_permute_vector_float.h>

/// workhorse of gsl calcultor for FAB
namespace FAB
{
    struct FAB_params
    {
        const Eigen::Vector3f * p_Avec;
        const Eigen::Vector3f * p_Bvec;
        const Eigen::Matrix3f * p_Amat;
        const Eigen::Matrix3f * p_Bmat;
        Eigen::Vector3f x_of_lam;
        float Ascr;
    };

    static float calc_FAB(const Eigen::Vector3f & Avec, const Eigen::Vector3f & Bvec, const Eigen::Matrix3f & Amat, const Eigen::Matrix3f & Bmat, Eigen::Vector3f & x_of_lam_new, float & lam_new);
    static double func_FAB(double lam, void * params);
}

namespace ellipsoid
{
    /// @brief class of ellipsoid
    /// @par ellispsoid
    class Ellipsoid
    {
        /// data members
        protected:
            float _machine_precision_float = 1.2e-38; ///< @note ZHL hack
            float _machine_inf_float = 3.3e38; ///< @note ZHL hack

            Eigen::Vector3f _vx = Eigen::Vector3f(1.,0.,0.); ///< unit vector x
            Eigen::Vector3f _vy = Eigen::Vector3f(0.,1.,0.); ///< unit vector y
            Eigen::Vector3f _vz = Eigen::Vector3f(0.,0.,1.); ///< unit vector z

            float _ep;
            float _sig;

            Eigen::Vector3f _Avec;

            Eigen::Matrix3f _Amat1;
            Eigen::Matrix3f _Amat2;

        /// interface
        public:
            inline Eigen::Vector3f Avec() const {return _Avec;}
            inline Eigen::Matrix3f Amat1() const {return _Amat1;}
            inline Eigen::Matrix3f Amat2() const {return _Amat2;}
            inline float machine_inf_float() const {return _machine_inf_float;}

            inline void set_Avec(const Eigen::Vector3f & Avec) {_Avec=Avec;}
            inline void set_Amat1(const Eigen::Matrix3f & Amat1) {_Amat1=Amat1;}
            inline void set_Amat2(const Eigen::Matrix3f & Amat2) {_Amat2=Amat2;}

        /// methods
        protected:
            inline Eigen::Matrix3f _ellipsoid_matrix(const Eigen::Vector3f & axis_lengths, const Eigen::Vector3f & u0, const Eigen::Vector3f & u1, const Eigen::Vector3f & u2) const;
        public:
            inline void rotate_Amat(const Eigen::Matrix3f & rot_matrix);
            bool energy_pbc(const ellipsoid::Ellipsoid & e, const float box_length, float & energy) const;

        /// constructor/destructor
        public:
            inline Ellipsoid();///< default constructor
            inline Ellipsoid(const float a1, const float a2, const float a3,
                             const float a11, const float a12, const float a13, const float a21, const float a22, const float a23,
                             const float ep, const float sig); ///< constructor for identical ellipsoids from individual length values
            inline Ellipsoid(const Eigen::Vector3f & coor_A,
                             const Eigen::Vector3f & axis_lengths_A1, const Eigen::Vector3f & axis_lengths_A2,
                             const float ep, const float sig); ///< constructor for identical ellipsoids from length vectors
            virtual ~Ellipsoid();

/*
        /// copy/assign constructor
        public:
            Ellipsoid & operator=(const Ellipsoid &);
            Ellipsoid(const Ellipsoid &);
            */
        
    };
}

/// get the inline definitions
#define ellipsoid_icc
#include "ellipsoid.icc"
#undef ellipsoid_icc

#endif
