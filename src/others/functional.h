//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-2.
//

#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H

#include <vector>
#include <LpCVT/common/types.h>
#include <LpCVT/combinatorics/mesh.h>

namespace Geex {
    /** Randomly sample pts on the mesh **/
    void sample_mesh_pts(
        const Mesh& M, std::vector<vec3>& pts,
        int pts_nb,
        bool volume = false);
}

#endif //FUNCTIONAL_H
