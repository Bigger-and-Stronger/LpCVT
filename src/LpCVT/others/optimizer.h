//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-2.
//

#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <LpCVT/common/types.h>
#include <LpCVT/combinatorics/mesh.h>

namespace Geex {
    class Optimizer {
    public:
        Optimizer() {};

        /** Read reference mesh **/
        void read_mesh(const std::string& filename);

        /** Read initial pts **/
        bool read_pts(const std::string& filename);

        void write_pts(const std::string& filename) const ;

        void optimize();
    private:
        Mesh M_;
        std::vector<vec3> pts_;

        bool volume_ = false;
        unsigned int p_ = 4;
    };
}

#endif //OPTIMIZE_H
