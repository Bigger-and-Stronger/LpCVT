//
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-2.
//

#include "LpCVT/others/functional.h"
#include "LpCVT/others/macro.h"

#include <random>

namespace Geex {
    /** Used by sample_mesh_pts. **/
    class ComputeAreaIndividually {
    public:
        explicit ComputeAreaIndividually(
            double& result, std::vector<double>& facet_areas
            ) : result_(result), facet_areas_(facet_areas) {
            result = 0.0 ;
        }

        void operator()(const vec3& p1, const vec3& p2, const vec3& p3) const {
            const double area = length(cross(p2 - p1, p3 - p1)) / 2.0;
            facet_areas_.push_back(area);
            result_ += area;
        }
    private:
        double& result_ ;
        std::vector<double>& facet_areas_;
    } ;

    class RandomSampleFacets {
    public:
        RandomSampleFacets(
            std::default_random_engine& randGen, int& pts_nb, std::vector<vec3>& pts
            ) : randGen_(randGen), pts_nb_(pts_nb), pts_(pts) {
            pts_.reserve(pts_.size() + pts_nb_);
        }

        void operator()(
            const vec3& p1, const vec3& p2, const vec3& p3
            ) const {
            for (size_t i = 0; i < pts_nb_; i++) {
                double r1 = randGen_() / static_cast<double>(RAND_MAX);
                double r2 = randGen_() / static_cast<double>(RAND_MAX);
                if (r1 + r2 > 1) {
                    const double tmp = r1; r1 = 1 - r2; r2 = 1 - tmp;
                }
                const double r3 = 1 - r1 - r2;

                pts_.push_back(r1*p1 + r2*p2 + r3*p3);
            }
        }
    private:
        std::default_random_engine& randGen_;
        int& pts_nb_;
        std::vector<vec3>& pts_;
    };

    void sample_mesh_pts(
        const Mesh& M, std::vector<vec3>& pts,
        const int pts_nb,
        const bool volume
        ) {

        if (volume) {
            // TODO: randomly sample the volume bounded by mesh M
        }
        else {
            /** Count the area of the facets of the mesh **/
            double sf_area;
            std::vector<double> facet_areas;
            facet_areas.reserve(M.nb_facets());

            M.for_each_triangle(ComputeAreaIndividually(sf_area, facet_areas));

            /** Randomly sample the mesh **/
            std::vector<int> pts_per_facet(M.nb_facets(), 0);

            std::default_random_engine randGen;
            std::discrete_distribution<int> thisDistribution(facet_areas.begin(), facet_areas.end());
            for (int i = 0; i < pts_nb; ++i) {
                const int sample_tri = thisDistribution(randGen); // current sample facet
                ++pts_per_facet[sample_tri];
            }

            for (unsigned int f = 0, f_end = M.nb_facets(); f < f_end; ++f)
                M.for_each_triangle(f, RandomSampleFacets(randGen, pts_per_facet[f], pts));
        }
    }
}