#pragma once

#include <iostream>
#include <cmath>

#include "sasatomBase.h"

// forward declaration
namespace sasmol
{
    class SasMol;
}

namespace sascalc 
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class Prop: public sasmol::SasAtomBase
    {
        public:
            int self_overlap(float &cutoff, int &frame); ///< (Brief description)

            void calc_mass() ; ///< (Brief description)

            float calc_volume() ; ///< (Brief description)

            std::vector<float> calc_com(int frame) ; ///< (Brief description)

            float calc_rg(int &frame) ; ///< (Brief description)

            void calc_pmi(int &frame) ; ///< (Brief description)

            float calc_rmsd(sasmol::SasMol &mol, int &frame) ; ///< (Brief description)

            std::vector<std::vector<float> > calc_minmax(int &frame) ; ///< (Brief description)

            void find_surf(int * const idx); ///< (Brief description)
            void find_non_aliphatic_Hs(int * const idx); ///< (Brief description)
            void find_backbone_Hs(int * const idx); ///< (Brief description)
    };
}
