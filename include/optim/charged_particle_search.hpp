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
#include <iterator>
#include <cassert>
#include <stdexcept>

#include "optim_utilities.hpp"
#include "particle_searcher.hpp"

#define CPO_OPT_EPS 1e-7
#define CPO_OPT_JIT 0.1

namespace optim {
    namespace CPS {
        template<typename XT, typename YT, typename VT>
        struct ParticleImpl : Result<XT, YT> {
            typedef Result<XT, YT> base_type;
            typedef std::vector<VT> velocity_type;
            typedef typename velocity_type::iterator velocity_iterator;
            typedef typename velocity_type::const_iterator const_velocity_iterator;

            ParticleImpl(size_t nArgs)
              : base_type(nArgs), particleNumber_(0), velocity_(nArgs) {
                std::fill(velocity_.begin(), velocity_.end(), 0.0);
            }

            void swap(ParticleImpl & other) {
                base_type::swap(*this);
                this->velocity_.swap(other.velocity_);
            }

            void id(int id) {
                particleNumber_ = id;
            }

            int id(void) const {
                return particleNumber_;
            }

            bool operator<(ParticleImpl const& other) const {
                return this->y > other.y;
            }

            const_velocity_iterator beginVelocity(void) const {
                return velocity_.begin();
            }

            velocity_iterator beginVelocity(void) {
                return velocity_.begin();
            }

            template<typename PositionIter>
            void updateForPositionChange(PositionIter bNewPos, PositionIter eNewPos) {
                // update the velocity
                std::transform(bNewPos, eNewPos, this->beginPosition(), this->beginVelocity(), std::minus<XT>());

                // update the position
                std::copy(bNewPos, eNewPos, this->position_.begin());
            }

        private:
            int
                particleNumber_;

            velocity_type
                velocity_;
        };





        template<
            typename XT, // argument type
            typename YT, // fitness type
            typename VT, // velocity type
            typename BoundsAdjustor=DefaultBoundsAdjustor<XT>, 
            typename Stopper=DefaultStoppingRule, 
            typename Bounds=Bounds<XT>, 
            typename Generator=UniformGenerator<XT> 
        >
        class ChargedParticleSearcher : public ParticleSearcher<
            ChargedParticleSearcher<XT, YT, VT, BoundsAdjustor, Stopper, Bounds, Generator>,
            XT,
            YT,
            ParticleImpl<XT, YT, VT>,
            Result<XT, YT>,
            BoundsAdjustor,
            Stopper,
            Generator,
            Bounds
        > {
            typedef ChargedParticleSearcher<XT, YT, VT, BoundsAdjustor, Stopper, Bounds, Generator> my_type;
            typedef ParticleSearcher<
                my_type,
                XT,
                YT,
                ParticleImpl<XT, YT, VT>,
                Result<XT, YT>,
                BoundsAdjustor,
                Stopper,
                Generator,
                Bounds
            > Base;

            typedef typename Base::bounds_type bounds_type;
            typedef typename Base::particle_type particle_type;
            typedef typename Base::particle_store_type particle_store_type;
            typedef typename particle_type::position_type position_type;
            typedef typename particle_type::velocity_type velocity_type;

            typedef typename particle_store_type::const_iterator particle_store_const_iterator;

            friend class ParticleSearcher<
                my_type,
                XT,
                YT,
                ParticleImpl<XT, YT, VT>,
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
                int seed,
                XT a
            ) : Base(bounds, boundsAdjustor, stopper, nParticles, seed), a_(a) {
                typedef typename bounds_type::const_iterator Iter;
                typedef typename particle_store_type::iterator PIter;

                assert(nParticles >= BEST_STORE_SIZE);
                // from equation 21.
                if(a_ <= 0.0) {
                    for(Iter i(this->getBounds().begin()); i!=this->getBounds().end(); ++i) {
                        a_ = std::max(a_, i->second - i->first);
                    }
                    a_ *= 0.1;
                }

                int
                    particleNumber(0);

                for(PIter b(this->getParticleStore().begin()), e(this->getParticleStore().end()); b!=e; ++b, ++particleNumber) {
                    b->id(particleNumber);
                }
            }

        private:
            // rules 4, 5 (equations 25, 26)
            // - implemented together as a proceedure for the sake of efficiency
            void updateParticle(
                particle_type const& j,
                typename position_type::iterator newPosIter, // assumed to be zeroed already.
                size_t iter, 
                size_t iterMax,
                Generator& gen
            ) const {
                typedef typename position_type::iterator posIter;
                typedef typename position_type::const_iterator constPosIter;
                typedef typename velocity_type::const_iterator constVelIter;
                typedef typename particle_store_type::const_iterator parStoreIter;

                XT
                    kv(0.5 * gen() * (1.0 - static_cast<XT>(iter)/static_cast<YT>(iterMax))),
                    ka(0.5 * gen() * (1.0 + static_cast<XT>(iter)/static_cast<YT>(iterMax))),
                    Fj, rij, qi;

                particle_store_type const&
                    particles(this->getParticleStore());

                particle_type const&
                    worst(this->getWorstParticle()),
                    best(this->getBestParticle());

                posIter
                    njx;

                constPosIter
                    ex(j.endPosition()),
                    ix, jx;

                constVelIter
                    jv;

                int
                    jId(j.id());
                
                // iterate over the that we are not going to drop particles
                for(parStoreIter i(particles.begin()), ie(particles.end()); i!=ie; ++i) {
                    if(jId!=i->id() && pijFunc(*i, j, best, gen) /* use rule 3*/) {
                        // use rule 1
                        rij = rijFunc(*i, j, best, CPO_OPT_EPS);
                        qi  = qiFunc(*i, worst, best);
                        Fj = (rij<a_) ? (qi*rij/(a_*a_*a_)) : (qi/(rij*rij)); // note no qi outside, as this cancels with mj (see text just after eqn 23)

                        if(rij != rij) throw std::runtime_error("rij is nan");
                        if(qi != qi) throw std::runtime_error("qi is nan");
                        if(Fj != Fj) throw std::runtime_error("Fj is nan");

                        // rule 5 (eqns 25, 26) - iterating over the particle's x
                        for(njx=newPosIter, ix=i->beginPosition(), jx=j.beginPosition(); jx!=ex; ++ix, ++jx, ++njx) {
                            *njx += ka*Fj*(*ix-*jx);
                        }
                    }
                }

                for(njx=newPosIter, jx=j.beginPosition(), jv=j.beginVelocity(); jx!=ex; ++njx, ++jx, ++jv) {
                    *njx += kv**jv + *jx;
                }
            }



            particle_type const& getReplacementParticle(void) const {
                typedef typename particle_store_type::const_iterator parStoreIter;

                parStoreIter
                    b(this->beginBestParticleStore()),
                    e(this->endBestParticleStore());

                size_t
                    nInStore(std::distance(b, e)),
                    index(static_cast<size_t>(this->rand() * nInStore));

                assert(nInStore > 0);

                if(index == nInStore) {
                    // can only happen if we got 1.0 from gen()
                    index = nInStore - 1;
                }

                return *(b + index);
            }


            particle_type const& getWorstParticle(void) const {
                return this->getParticleStore().back();
            }


            particle_type const& getBestParticle(void) const {
                return this->getParticleStore().front();
            }


            particle_store_const_iterator beginBestParticleStore(void) const {
                return this->getParticleStore().begin();
            }


            particle_store_const_iterator endBestParticleStore(void) const {
                return this->getParticleStore().begin() + BEST_STORE_SIZE;
            }



            //rule1 (qi)
            static typename particle_type::result_type qiFunc(
                particle_type const& i,
                particle_type const& worst,
                particle_type const& best
            ) {
                return (i.y - worst.y) / (best.y - worst.y);
            }



            // rule1 (rij)
            static typename position_type::value_type rijFunc(
                particle_type const& particle1,
                particle_type const& particle2,
                particle_type const& best,
                double eps
            ) {
                typedef typename position_type::const_iterator Iter;
                typedef typename position_type::value_type distance;

                distance
                    num(0.0),
                    den(0.0),
                    tmp;

                Iter
                    p1b(particle1.beginPosition()),
                    p2b(particle2.beginPosition()),
                    p1e(particle1.endPosition()),
                    bestb(best.beginPosition());

                for(; p1b!=p1e; ++p1b, ++p2b, ++bestb) {
                    tmp = *p1b - *p2b;
                    num += tmp * tmp;
                    tmp = (*p1b + *p2b) / 2.0 - *bestb;
                    den += tmp * tmp;
                }

                return sqrt(num) / (sqrt(den) + eps);
            }



            // rule 3
            static bool pijFunc(
                particle_type const& i,
                particle_type const& j,
                particle_type const& best,
                Generator& gen
            ) { // j will attract i
                return (j.y >= i.y) || (((i.y - best.y) / (j.y - i.y)) > gen());
            }


        private:
            XT
                a_;
        };
    } // end namespace CPS
} // end namespace optim





template<typename XT, typename YT, typename VT, typename FUNC, typename Bounds>
typename optim::CPS::ChargedParticleSearcher<XT, YT, VT, optim::DefaultBoundsAdjustor<XT>, optim::DefaultStoppingRule, Bounds>::result_type
chargedParticleOptimise(
    FUNC const& func,
    Bounds const& bounds,
    size_t nIters,
    size_t nParticles,
    int seed,
    XT cmcr=0.5, 
    XT par=0.1, 
    XT a=-1.0, 
    XT jitter=CPO_OPT_JIT
) {
    typedef optim::CPS::ChargedParticleSearcher<XT, YT, VT, optim::DefaultBoundsAdjustor<XT>, optim::DefaultStoppingRule, Bounds> Searcher;
    return Searcher(
        bounds,
        optim::DefaultBoundsAdjustor<XT>(cmcr, par, jitter),
        optim::DefaultStoppingRule(),
        nParticles,
        seed,
        a
    ).run(
        func, 
        nIters
    );
}

#endif //OPT_CHARGED_PARTICLE_HEADER_INCLUDED_PEQOTUPWGOIPW5YUP50
