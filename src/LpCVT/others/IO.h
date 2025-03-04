//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-2.
//

#ifndef IO_H
#define IO_H

#include <string>
#include <vector>
#include <LpCVT/common/types.h>

namespace Geex{
    /** Loads points from file.
    * [Canjia Huang] file format: "v" << x << y << z **/
    bool load_pts(const std::string& filename, std::vector<vec3>& pts);

    /** write pts to file **/
    bool write_pts(const std::string& filename, const std::vector<vec3>& pts);

    bool write_vtk(const std::string& filename, const std::vector<vec3>& pts);
}
#endif //IO_H
