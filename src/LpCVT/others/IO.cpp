//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-2.
//

#include "LpCVT/others/IO.h"
#include "LpCVT/others/macro.h"

#include <fstream>

namespace Geex{
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
}