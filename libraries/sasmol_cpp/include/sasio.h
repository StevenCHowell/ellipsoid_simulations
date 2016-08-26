#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "sasatomBase.h"

/********* declarations   ******************/



namespace sasio
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class Files: public sasmol::SasAtomBase
    {
        public:
            void read_xyz(const std::string &filename) ; ///< (Brief description)

            void write_dcd(const std::string &filename); ///< (Brief description)

            void read_pdb(const std::string &filename); ///< (Brief description)
            void read_pdb_anisou_hack(const std::string &filename, std::vector<Eigen::Matrix3f> & Us, std::vector<Eigen::Matrix3f> & As); ///< (Brief description)

            void write_pdb(const std::string &filename, int frame); ///< (Brief description)
            void write_pdb_anisou_hack(const std::string &filename, int frame, const std::vector<Eigen::Matrix3f> & Us, const std::vector<Eigen::Matrix3f> & As); ///< (Brief description)

            void write_dcd_step(FILE *outfile, int frame, int step); ///< (Brief description)
            void write_dcd_step_pbc(FILE *outfile, int frame, int step, const float x0, const float y0, const float z0, const float lx, const float ly, const float lz); ///< @note ZHL use only

            void read_dcd_step(FILE *dcd_infile, int frame, int natoms, int reverseEndian) ; ///< (Brief description)

            void read_dcd(std::string dcd_input_filename); ///< (Brief description)

            void resize_array() ; ///< (Brief description)

            void _set_unique_attributes() ; ///< (Brief description)
    };

    /// @note Where is my home?
    FILE *open_write_dcd_file(const std::string &filename, int natoms, int nset); ///< (Brief description)
    int write_dcd_header(FILE *outfile, const std::string &filename, int natoms, int nset); ///< (Brief description)
    FILE *open_dcd_read(std::string &filename, int &natoms, int &nset, int &reverseEndian) ; ///< (Brief description)
}

