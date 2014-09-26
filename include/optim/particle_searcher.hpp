/*
 * file particle_searcher.hpp
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

#ifndef OPTIM_PARTICLE_SEARCHER_HEADER_INCLUDED_ODIJWPEUP5OIJRSDLFGJOH
#define OPTIM_PARTICLE_SEARCHER_HEADER_INCLUDED_ODIJWPEUP5OIJRSDLFGJOH

#include <vector>

#include "optim_utilities.hpp"

#define CPO_CM_FRAC 4

namespace optim {
    template<
        typename Super,
        typename XT, 
        typename YT,
        typename Particle,
        typename Result=Particle,
        typename BoundsAdjustor=DefaultBoundsAdjustor<XT>,
        typename Stopper=DefaultStoppingRule,
        typename Generator=UniformGenerator<XT>,
        typename Bounds=Bounds<XT>,
        typename BestStore=BestStoreImpl<Particle>,
        typename ParticleStore=std::vector<Particle>
    >
    class ParticleSearcher {
    public:
        typedef Particle particle_type;
        typedef Result result_type;
        typedef BestStore best_particle_store_type;
        typedef Bounds bounds_type;
        typedef ParticleStore particle_store_type;

    private:
        typedef typename particle_store_type::iterator ParticleVectorIterator;

    public:
        ParticleSearcher(
            bounds_type const& bounds,
            BoundsAdjustor const& boundsAdjustor,
            Stopper const& stopper,
            size_t nParticles
        ) : bounds_(bounds),
            stopper_(stopper),
            boundsAdjustor_(boundsAdjustor),
            particles_(nParticles),
            bestStore_(nParticles/CPO_CM_FRAC)
        {
            if(nParticles < CPO_CM_FRAC) {
                // should emit a warning (somehow) here.
                particles_ = particle_store_type(CPO_CM_FRAC);
                bestStore_ = best_particle_store_type(1);
            }
        }


        template<typename FUNC>
        particle_type const& run(FUNC func, size_t iters) {
            typedef typename particle_type::arg_type Position;
            typedef typename particle_type::arg_type::iterator PIter;
            typedef typename particle_type::velocity_type::iterator VIter;

            size_t size(bounds_.size()), currParticleIndex, nextWorstIndex;
            YT worstValue(-std::numeric_limits<YT>::infinity()), nextWorstValue(worstValue);
            ParticleVectorIterator pvi;
            Position oldPos(bounds_.size());
            // use rule 2
            generateRandomParticles(particles_, bounds_, gen_);
            // evaluate func for all the particles.
            for(pvi=particles_.begin(), currParticleIndex=0; pvi!=particles_.end(); ++pvi, ++currParticleIndex) {
                pvi->y = func(&pvi->x[0], size);
                // use rule 6
                bestStore_.addToStore(*pvi);
                if(pvi->y > nextWorstValue) {
                    nextWorstValue = pvi->y;
                    nextWorstIndex = currParticleIndex;
                }
                if(pvi->y > worstValue || currParticleIndex == worstParticleIndex_) {
                    worstValue = nextWorstValue;
                    worstParticleIndex_ = nextWorstIndex;
                    nextWorstValue = -std::numeric_limits<YT>::infinity();
                }
            }
            
            // use rule 8 (stopping proceedure) - max number of iterations.
            for(currentIteration_=0; currentIteration_<iters; ++currentIteration_) {
                for(pvi=particles_.begin(), currParticleIndex=0; pvi!=particles_.end(); ++pvi, ++currParticleIndex) {
                    // use rules 4 and 5
                    std::copy(pvi->x.begin(), pvi->x.end(), oldPos.begin());
                    static_cast<Super*>(this)->updateParticle(*pvi, currentIteration_, iters, gen_);
                    if(!boundsAdjustor_.checkBounds(*pvi, bounds_)) {
                        // use rule 7
                        boundsAdjustor_.adjustBounds(*pvi, *this, gen_);
                        pvi->adjustForPositionChange(oldPos);
                    }
                    pvi->y = func(&pvi->x[0], size);
                    // use rule 6
                    if(!bestStore_.addToStore(*pvi) && pvi->y > nextWorstValue) {
                        nextWorstValue = pvi->y;
                        nextWorstIndex = currParticleIndex;
                    }
                    if(pvi->y > worstValue || currParticleIndex == worstParticleIndex_) {
                        worstValue = nextWorstValue;
                        worstParticleIndex_ = nextWorstIndex;
                        nextWorstValue = -std::numeric_limits<YT>::infinity();
                    }
                }
                if(stopper_(*this)) break;
            }
            return bestStore_.getBestParticle();
        }

        particle_store_type const& getParticleStore(void) const {
            return particles_;
        }

        best_particle_store_type const& getBestParticleStore(void) const {
            return bestStore_;
        }

        particle_type const& getWorstParticle(void) const {
            return particles_[worstParticleIndex_];
        }

        particle_type const& getBestParticle(void) const {
            return bestStore_.getBestParticle();
        }

        bounds_type const& getBounds(void) const {
            return bounds_;
        }

    private:
        Generator gen_;
        bounds_type bounds_;
        BoundsAdjustor boundsAdjustor_;
        Stopper stopper_;
        particle_store_type particles_;
        best_particle_store_type bestStore_;
        size_t worstParticleIndex_;
        size_t currentIteration_;

    };
} // end namespace optim

#endif // OPTIM_PARTICLE_SEARCHER_HEADER_INCLUDED_ODIJWPEUP5OIJRSDLFGJOH