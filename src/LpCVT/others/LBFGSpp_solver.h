//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-4.
//

#include "LpCVT/others/optimizer.h"
#include <LpCVT/common/types.h>
#include <LpCVT/combinatorics/mesh.h>
// Eigen
#include <Eigen/Core>
// LBFGSpp
#include "LpCVT/external/LBFGSpp/LBFGS.h"

namespace LBFGSppSolver{
    using namespace Geex;

    /** Used for optimizing verbose**/
    class verbose_f_g {
    public:
        verbose_f_g(const unsigned int iter, const double f, std::vector<double>& g
            ) : iter_(iter), f_(f), g_(g) {}

        void operator()() const {
            double gnorm;
            for (const auto& i : g_)
                gnorm += i * i;
            gnorm = ::sqrt(gnorm);

            std::cout.precision(16);
            VERBOSE(
                "iter=" << iter_ << " "
                "f=" << std::scientific << f_ << "  g=" << gnorm
                );
        }
    private:
        unsigned int iter_;
        double f_;
        std::vector<double>& g_;
    };

    /** Used for optimizing **/
    class LpCVT_iterate {
    public:
        LpCVT_iterate(Mesh* M, std::vector<vec3>& pts, const bool volume, const unsigned int p
            ) : M_(M), pts_(pts), volume_(volume), p_(p) {}

        std::vector<vec3>& pts() const { return pts_; }

        template<class T> void for_each_iter(const T& do_it) {
            do_it();
        }

        double operator()( // TODO: This function should be more efficient by avoiding vector copying.
            const Eigen::VectorXd& x, Eigen::VectorXd& grad
            ) {
            for (std::size_t i = 0, i_end = pts_.size(); i != i_end; ++i) {
                pts_[i][0] = x[3 * i];
                pts_[i][1] = x[3 * i + 1];
                pts_[i][2] = x[3 * i + 2];
            }

            std::vector<double> g;
            const double f = compute_f_grad(M_, pts_, volume_, p_, g);

            /** Verbose fx and grad **/
            if (verbose_ && ++iter_ % verbose_gap_ == 0)
                for_each_iter(verbose_f_g(iter_, f, g));

            for (std::size_t i = 0, i_end = g.size(); i < i_end; ++i)
                grad[i] = g[i];
            return f;
        }
    private:
        Mesh* M_;
        std::vector<vec3>& pts_;
        bool volume_;
        unsigned int p_;

        unsigned int iter_ = 0;
        bool verbose_ = true;
        unsigned int verbose_gap_ = 1;
    };

    class LBFGSpp_solver {
    public:
        LBFGSpp_solver(Mesh* M, std::vector<vec3>& pts, const bool volume, const unsigned int p
            ) : M_(M), pts_(pts), volume_(volume), p_(p) {}

        double solve() const {
            LBFGSpp::LBFGSParam<double> param;
            param.epsilon = param_epsilon_;
            param.max_iterations = param_max_iterations_;

            // Create solver and function object
            LBFGSpp::LBFGSSolver<double> solver(param);
            LpCVT_iterate fun(M_, pts_, volume_, p_);

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

            pts_.swap(fun.pts());
            return fx;
        }
    private:
        Mesh* M_;
        std::vector<vec3>& pts_;
        bool volume_;
        unsigned int p_;

        double param_epsilon_ = 1e-6;
        unsigned int param_max_iterations_ = 100;
    };
}