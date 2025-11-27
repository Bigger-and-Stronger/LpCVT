//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-2.
//

#include "LpCVT/others/optimizer.h"
#include "LpCVT/others/IO.h"
#include "LpCVT/others/Geex_test_algebra.h"

namespace Geex {
    double compute_f_grad(
        Mesh* M, const std::vector<vec3>& pts,
        const bool volume, const unsigned int p,
        std::vector<double>& grad
        ) {
        std::vector<int> I;
        std::vector<vec3> C;
        std::vector<int> F;
        get_combinatorics(M, pts, I, C, F, volume);
        const unsigned int nb_integration_simplices = static_cast<unsigned int>(I.size()) / 10;

        std::vector<mat3> Mat(nb_integration_simplices);
        for (auto & i : Mat)
            i.load_identity();
        std::vector<plane3> Q(M->nb_facets());
        for (unsigned int i = 0, i_end = M->nb_facets(); i < i_end; ++i)
            Q[i] = M->facet_plane(i);

        std::vector<double>(pts.size() * 3).swap(grad);
        return compute_F_Lp(volume, p, M, I, C, pts, Q, Mat, grad);
    }
}