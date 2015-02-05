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
#include <algorithm>

#include "optim_utilities.hpp"

#define BEST_STORE_SIZE 5

namespace optim {
    template<
        typename Super,
        typename XT, // argument type
        typename YT, // fitness type
        typename Particle,
        typename Result,
        typename BoundsAdjustor=DefaultBoundsAdjustor<XT>,
        typename Stopper=DefaultStoppingRule,
        typename Generator=UniformGenerator<XT>,
        typename Bounds=Bounds<XT>,
        typename ParticleStore=std::vector<Particle>
    >
    class ParticleSearcher {
        typedef Super super_type;

    public:
        typedef Particle particle_type;
        typedef Result result_type;
        typedef Bounds bounds_type;
        typedef ParticleStore particle_store_type;

    public:
        ParticleSearcher(
            bounds_type const& bounds,
            BoundsAdjustor const& boundsAdjustor,
            Stopper const& stopper,
            size_t nParticles,
            int seed
        ) : gen_(seed),
            bounds_(bounds),
            boundsAdjustor_(boundsAdjustor),
            stopper_(stopper),
            particles_(nParticles, particle_type(bounds.size()))
        {
            size_t
                nArgs(bounds.size());

            for(typename bounds_type::const_iterator b(bounds.begin()), e(bounds.end()); b!=e; ++b) {
                if(b->first >= b->second) throw std::runtime_error("bounds.first >= bound.second");
            }
        }


        template<typename FUNC>
        particle_type const& run(FUNC fit, size_t iters) {
            typedef typename particle_type::position_type position_type;
            typedef typename position_type::iterator PosIter;
            typedef typename particle_type::velocity_type::iterator VIter;
            typedef typename particle_store_type::iterator storeIter;

            size_t
                nArgs(bounds_.size());

            position_type
                newPos(nArgs * particles_.size()); // buffer for new positions

            PosIter
                newPosIter;

            // use rule 2 - fill the particles with random locatiions
            generateRandomParticles(particles_, bounds_, gen_);

            storeIter
                pvi, pve;

            // evaluate the fitness of every particle
            for(pvi=particles_.begin(), pve=particles_.end(); pvi!=pve; ++pvi) {
                pvi->y = fit(pvi->positionData(), nArgs);
            }

            // sort the particles according to fitness
            std::sort(particles_.begin(), particles_.end());

            // run the iterations.
            for(currentIteration_=0; currentIteration_<iters; ++currentIteration_) {

                // initialise the new particle positions to zero
                std::fill(newPos.begin(), newPos.end(), 0.0);

                // get the new positions of the particles
                for(pvi=particles_.begin(), pve=particles_.end(), newPosIter=newPos.begin(); pvi!=pve; ++pvi, newPosIter+=nArgs) {

                    // use rules 4 and 5 - get the particles new position
                    static_cast<Super*>(this)->updateParticle(*pvi, newPosIter, currentIteration_, iters, gen_);
                    if(!boundsAdjustor_.checkBoundsOK(newPosIter, bounds_)) {
                        // use rule 7 - ajust the new position if reqired
                        boundsAdjustor_.adjustBounds(*pvi, newPosIter, *this, gen_);
                    }
                }

                // update the positions of the particles
                for(pvi=particles_.begin(), pve=particles_.end(), newPosIter=newPos.begin(); pvi!=pve; ++pvi, newPosIter+=nArgs) {
                    pvi->updateForPositionChange(newPosIter, newPosIter+nArgs);
                    pvi->y = fit(pvi->positionData(), nArgs);
                }

                std::sort(particles_.begin(), particles_.end());

                // use rule 8 (stopping proceedure)
                if(stopper_(*this)) break;
            }

            return static_cast<super_type const*>(this)->getBestParticle();
        }

        bounds_type const& getBounds(void) const {
            return bounds_;
        }

        particle_type const& getReplacementParticle(void) const {
            return static_cast<super_type const*>(this)->getReplacementParticle();
        }

    protected:
        particle_store_type const& getParticleStore(void) const {
            return particles_;
        }        

        particle_store_type& getParticleStore(void) {
            return particles_;
        }

        double rand(void) const {
            return gen_();
        }

    private:
        size_t
            currentIteration_;

        mutable Generator
            gen_;

        bounds_type
            bounds_;

        BoundsAdjustor
            boundsAdjustor_;

        Stopper
            stopper_;

        particle_store_type
            particles_;
    };
} // end namespace optim

#endif // OPTIM_PARTICLE_SEARCHER_HEADER_INCLUDED_ODIJWPEUP5OIJRSDLFGJOH