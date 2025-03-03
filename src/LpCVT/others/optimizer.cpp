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

    /** Used for optimizing verbose**/
    class Iteration_verbose {
    public:
        Iteration_verbose(unsigned int iter, double f, std::vector<double>& g
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

    /** Used for optimizing point animation **/
    class Iteration_point_animation {
    public:
        Iteration_point_animation(const unsigned int iter, std::vector<vec3>& pts
            ) : iter_(iter), pts_(pts) {}

        void operator()() const {
            std::ofstream out(filename_ + "." + std::to_string(iter_));
            out << "x coord, y coord, z coord, scalar" << std::endl;
            for (const auto& p : pts_) {
                out << p[0] << ", " << p[1] << ", " << p[2] << ", " << "1" << std::endl;
            }
            out.close();
        }
    private:
        std::string filename_ = "../data/point_animation.csv";
        unsigned int iter_;
        std::vector<vec3>& pts_;
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

            /** Verbose func and grad **/
            if (verbose_ && ++iter_ % verbose_gap_ == 0) {
                for_each_iter(Iteration_verbose(iter_, f, g));
            }
            for_each_iter(Iteration_point_animation(iter_, pts_));



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
        unsigned int verbose_gap_ = 10;
    };

    class LBFGSpp_solver {
    public:
        LBFGSpp_solver(Mesh* M, std::vector<vec3>& pts, const bool volume, const unsigned int p
            ) : M_(M), pts_(pts), volume_(volume), p_(p) {}

        double operator()() const {
            LBFGSpp::LBFGSParam<double> param;
            param.epsilon = 1e-6;
            param.max_iterations = 10;

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
    };

    void Optimizer::optimize(
    ) {
        solve(LBFGSpp_solver(&M_, pts_, volume_, p_));
    }
}