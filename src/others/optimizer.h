//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-2.
//

#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <LpCVT/common/types.h>
#include <LpCVT/combinatorics/mesh.h>

namespace Geex {
    /** For a given reference mesh and the current sampling points,
     *  calculate the F_{L_p} energy and the gradient for each point.
     * @param[out] grad the gradient for each point
     * @return total F_{L_p}
     */
    double compute_f_grad(
        Mesh* M, const std::vector<vec3>& pts,
        bool volume, unsigned int p,
        std::vector<double>& grad);

    /** optimize the sample points
     * class "solver" needs to have the constructor (M, pts, volume, p),
     * and a "solve()" function to optimize the pts **/
    template<class T> void optimize(const T& solver) {
        solver.solve();
    }
}

#endif //OPTIMIZE_H
