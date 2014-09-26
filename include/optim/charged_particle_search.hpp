/*
 * file charged_particle_search.hpp
 * Copyright (C) 2011       Simon Knapp
 *
 * This progrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OPT_CHARGED_PARTICLE_HEADER_INCLUDED_PEQOTUPWGOIPW5YUP50
#define OPT_CHARGED_PARTICLE_HEADER_INCLUDED_PEQOTUPWGOIPW5YUP50

// based on the document: A novel heuristic optimization method - charged system search - 2010 kaveh.pdf
// found at: https://svn-d1.mpi-inf.mpg.de/AG1/MultiCoreLab/papers/KavehTalatahari10%20-%20Charged%20System%20Search.pdf
// - equation numbers refer to equations in that paper.

#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <cassert>

#include "optim_utilities.hpp"
#include "particle_searcher.hpp"

#define CPO_OPT_EPS 1e-7
#define CPO_OPT_JIT 0.1

namespace optim {
    namespace CPS {
        template<typename XT, typename YT>
        struct ParticleImpl : Result<XT, YT> {
            typedef std::vector<XT> velocity_type;
            ParticleImpl(void) {}
            ParticleImpl(size_t size) : Result<XT, YT>(size), v(size) { std::fill_n(v.begin(), size, 0.0); }
            template<typename Position>
            void adjustForPositionChange(Position const& b) {
                std::transform(this->x.begin(), this->x.end(), b.begin(), v.begin(), std::minus<XT>());
            }
            velocity_type v;
        };





        template<
            typename XT, 
            typename YT, 
            typename BoundsAdjustor=DefaultBoundsAdjustor<XT>, 
            typename Stopper=DefaultStoppingRule, 
            typename Bounds=Bounds<XT>, 
            typename Generator=UniformGenerator<XT> 
        >
        class ChargedParticleSearcher : public ParticleSearcher<
            ChargedParticleSearcher<XT, YT, BoundsAdjustor, Stopper, Bounds, Generator>,
            XT,
            YT,
            ParticleImpl<XT, YT>,
            Result<XT, YT>,
            BoundsAdjustor,
            Stopper,
            Generator,
            Bounds
        > {
            typedef ParticleSearcher<
                ChargedParticleSearcher<XT, YT, BoundsAdjustor, Stopper, Bounds, Generator>,
                XT,
                YT,
                ParticleImpl<XT, YT>,
                Result<XT, YT>,
                BoundsAdjustor,
                Stopper,
                Generator,
                Bounds
            > Base;

            typedef typename Base::bounds_type bounds_type;
            typedef typename Base::particle_type particle_type;
            typedef typename Base::particle_store_type particle_store_type;

            friend class ParticleSearcher<
                ChargedParticleSearcher<XT, YT, BoundsAdjustor, Stopper, Bounds, Generator>,
                XT,
                YT,
                ParticleImpl<XT, YT>,
                Result<XT, YT>,
                BoundsAdjustor,
                Stopper,
                Generator,
                Bounds
            >;

        public:
            ChargedParticleSearcher(
                bounds_type const& bounds,
                BoundsAdjustor const& boundsAdjustor,
                Stopper const& stopper,
                size_t nParticles,
                XT a
            ) : Base(bounds, boundsAdjustor, stopper, nParticles), a_(a)
            {
                typedef typename bounds_type::const_iterator Iter;
                // from equation 21.
                if(a <= 0.0) {
                    for(Iter i(this->getBounds().begin()); i!=this->getBounds().end(); ++i) a_ = std::max(a_, i->second - i->first);
                    a_ *= 0.1;           
                }
            }

        private:
            // rules 4, 5 (equations 25, 26)
            // - implemented together as a proceedure for the sake of efficiency
            void updateParticle(
                particle_type& j,
                size_t iter, 
                size_t iterMax,
                Generator& gen
            ) {
                typedef typename particle_type::arg_type Position;
                typedef typename particle_type::velocity_type Velocity;
                typedef typename Position::iterator PIter;
                typedef typename Velocity::iterator VIter;
                typedef typename particle_store_type::const_iterator PvIter;

                YT iterFrac(1.0 + (YT)iter/(YT)iterMax), rij, qi, randj1(0.5*gen()), randj2(0.5*gen()), tmp;
                particle_store_type const& particles(this->getParticleStore());
                particle_type const& worst(this->getWorstParticle()), best(this->getBestParticle());
                PIter ix(j.x.begin());
                VIter iv(j.v.begin());
                
                for(size_t currentDim(0); ix!=j.x.end(); ++ix, ++iv, ++currentDim) { // iterate over the x values.
                    tmp = randj2*iterFrac**iv + *ix;
                    for(PvIter i(particles.begin()); i!=particles.end(); ++i) { // iterate over the other particles.
                        if(pijFunc(*i, j, best) /* use rule 3*/) {
                            // use rule 1
                            rij = rijFunc(i->x, j.x, best.x, CPO_OPT_EPS);
                            qi  = qiFunc(*i, worst, best);
                            tmp += 
                                randj1*iterFrac*qi*
                                ((rij<a_) ? (rij/(a_*a_*a_)) : (1.0/(rij*rij)))*
                                (i->x[currentDim]-*ix);
                        }
                    }
                    *iv = tmp - *ix;
                    *ix = tmp;
                }
            }

            //rule1 (qi)
            template<typename Particle>
            typename Particle::result_type qiFunc(Particle const& i, Particle const& worst, Particle const& best) {
                return (i.y - worst.y) / (best.y - worst.y);
            }

            // rule1 (rij)
            template<typename Position>
            typename Position::value_type rijFunc(Position const& i, Position const& j, Position const& best, double eps) {
                typedef typename Position::const_iterator Iter;
                typedef typename Position::value_type distance;
                distance num(0.0), den(0.0), tmp;
                Iter ib(i.begin()), jb(j.begin()), ie(i.end()), bestb(best.begin());
                for(; ib!=ie; ++ib, ++jb, ++bestb) {
                    tmp = *ib - *jb; num += tmp * tmp;
                    tmp = (*ib + *jb) / 2.0 - *bestb; den += tmp * tmp;
                }
                return sqrt(num) / (sqrt(den) + eps);
            }

            // rule 3
            template<typename Particle>
            bool pijFunc(Particle const& i, Particle const& j, Particle const& best) { // j will attract i
                return (j.y >= i.y) || ((i.y - best.y) / (j.y - i.y));
            }

        private:
            XT a_;
        };
    } // end namespace CPS
} // end namespace optim





template<typename XT, typename YT, typename FUNC, typename Bounds>
typename optim::CPS::ChargedParticleSearcher<XT, YT, optim::DefaultBoundsAdjustor<XT>, optim::DefaultStoppingRule, Bounds>::result_type chargedParticleOptimise(
    FUNC const& func,
    Bounds const& bounds,
    size_t nIters,
    size_t nParticles,
    XT cmcr=0.5, 
    XT par=0.1, 
    XT a=-1.0, 
    XT jitter=CPO_OPT_JIT
) {
    return optim::CPS::ChargedParticleSearcher<XT, YT, optim::DefaultBoundsAdjustor<XT>, optim::DefaultStoppingRule, Bounds>(
        bounds,
        optim::DefaultBoundsAdjustor<XT>(cmcr, par, jitter),
        optim::DefaultStoppingRule(),
        nParticles,
        a
    ).run(
        func, 
        nIters
    );
}

#endif //OPT_CHARGED_PARTICLE_HEADER_INCLUDED_PEQOTUPWGOIPW5YUP50
