//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-2.
//

#include "LpCVT/others/optimizer.h"
#include "LpCVT/others/IO.h"
#include "LpCVT//others/Geex_test_combinatorics.h"
#include "LpCVT/others/Geex_test_algebra.h"
// Eigen
#include <Eigen/Core>
// LBFGSpp
#include "LBFGSpp/LBFGS.h"

namespace Geex {
    void Optimizer::read_mesh(const std::string& filename
        ) {
        M_.load(filename);
    }

    bool Optimizer::read_pts(const std::string& filename
        ) {
        return load_pts(filename, pts_);
    }

    void Optimizer::write_pts(const std::string& filename
        ) const {
        Geex::write_pts(filename, pts_);
    }

    /** Used for optimizing **/
    class Compute_f_grad {
    public:
        Compute_f_grad(Mesh* M, std::vector<vec3>& pts, bool volume, unsigned int p
            ) : M_(M), pts_(pts), volume_(volume), p_(p) {
        }

        double operator()(
            const Eigen::VectorXd& x, Eigen::VectorXd& grad
            ) {
            for (int i = 0, i_end = pts_.size(); i != i_end; ++i) {
                pts_[i][0] = x[3 * i];
                pts_[i][1] = x[3 * i + 1];
                pts_[i][2] = x[3 * i + 2];
            }

            std::vector<int> I;
            std::vector<vec3> C;
            std::vector<int> F;
            get_combinatorics(M_, pts_, I, C, F, volume_);
            const unsigned int nb_integration_simplices = static_cast<unsigned int>(I.size()) / 10;

            std::vector<mat3> M(nb_integration_simplices);
            for (auto & i : M)
                i.load_identity();
            std::vector<plane3> Q(M_->nb_facets());
            for (unsigned int i = 0, i_end = M_->nb_facets(); i < i_end; ++i)
                Q[i] = M_->facet_plane(i);

            std::vector<double> g(pts_.size() * 3);
            const double f = compute_F_Lp(volume_, p_, M_, I, C, pts_, Q, M, g);

            {
                double gnorm = 0.0;
                for (const auto& i : g)
                    gnorm += i * i;
                VERBOSE("LBFGS iteration f: " << f << ", gnorm: " << gnorm);
            }

            // return
            for (size_t i = 0, i_end = g.size(); i < i_end; ++i)
                grad[i] = g[i];
            return f;
        }

        std::vector<vec3>& pts() { return pts_; }
    private:
        Mesh* M_;
        std::vector<vec3> pts_;
        bool volume_;
        unsigned int p_;
    };

    void Optimizer::optimize(
        ) {
        // Set up parameters
        LBFGSpp::LBFGSParam<double> param;
        param.epsilon = 1e-6;
        param.max_iterations = 100;

        // Create solver and function object
        LBFGSpp::LBFGSSolver<double> solver(param);
        Compute_f_grad fun(&M_, pts_, volume_, p_);

        // Initial guess
        Eigen::VectorXd x = Eigen::VectorXd::Zero(pts_.size() * 3);
        for (std::size_t i = 0, i_end = pts_.size(); i < i_end; ++i) {
            x[3 * i] = pts_[i][0];
            x[3 * i + 1] = pts_[i][1];
            x[3 * i + 2] = pts_[i][2];
        }

        // x will be overwritten to be the best point found
        double fx;
        const int niter = solver.minimize(fun, x, fx);

        VERBOSE(niter << " iterations");
        VERBOSE("f(x) = " << fx);

        // rewrite
        pts_ = fun.pts();
    }
}