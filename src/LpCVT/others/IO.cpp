//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-2.
//

#include "LpCVT/others/IO.h"
#include "LpCVT/others/macro.h"
#include "LpCVT/common/line_stream.h"

#include <fstream>

namespace Geex{
    bool load_pts(
        const std::string& filename, std::vector<vec3>& pts
        ) {
        std::vector<vec3>().swap(pts);
        std::ifstream in_stream(filename);
        if (!in_stream) {
            WARNING("Could not open " << filename);
            return false;
        }
        LineInputStream in(in_stream);
        while (!in.eof()) {
            in.get_line();
            std::string kw;
            in >> kw;
            if (kw == "v") {
                vec3 v;
                in >> v;
                pts.push_back(v);
            }
        }

        return true;
    }

    bool write_pts(
        const std::string& filename, const std::vector<vec3>& pts
        ){
        std::ofstream out(filename);
        if (!out.good()) {
            WARNING("Could not write " << filename);
            return false;
        }

        VERBOSE("write pts to " << filename);
        for (const auto& pt : pts)
            out << "v" << " " << pt << std::endl;

        return true;
    }

    bool write_vtk(
        const std::string& filename, const std::vector<vec3>& pts
        ) {
        std::ofstream out(filename);
        if (!out.good()) {
            WARNING("Could not write " << filename);
            return false;
        }

        out << "# vtk DataFile Version 2.0" << std::endl;
        out << filename << std::endl;
        out << "ASCII" << std::endl;
        out << "DATASET UNSTRUCTURED_GRID" << std::endl;

        /** Points **/
        out << "POINTS " << pts.size() << " FLOAT" << std::endl;
        for (const auto& p : pts)
            out << p << std::endl;

        /** Cells **/
        out << "CELLS " << pts.size() << " " << 2 * pts.size() << std::endl;
        for (std::size_t i = 0, i_end = pts.size(); i < i_end; ++i)
            out << "1" << " " << i << std::endl;

        /** Cell type **/
        out << "CELL_TYPES " << pts.size() << std::endl;
        for (std::size_t i = 0, i_end = pts.size(); i < i_end; ++i)
            out << "1" << std::endl;

        return true;
    }
}