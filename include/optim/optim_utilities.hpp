/*
 * file optim_utilities.hpp
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

#ifndef OPTIM_UTILITIES_HEADER_INCLUDED_SDFPGJWPY8U60Y9WU850W49
#define OPTIM_UTILITIES_HEADER_INCLUDED_SDFPGJWPY8U60Y9WU850W49

#include <vector>
#include <utility>
#include <set>
#include <cassert>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

namespace optim {





    template<typename NT, typename G=boost::mt19937>
    struct UniformGenerator {
        UniformGenerator(void) {}
        template<typename UIntType> explicit UniformGenerator(UIntType seed) : generator_(seed) {}
        NT operator()(void) { return dist_(generator_); }

    private:
        G generator_;
        boost::uniform_01<NT> dist_;
    };





    template<typename NT>
    struct Bounds : std::vector< std::pair<NT, NT> > {
        Bounds(size_t nDim) : std::vector< std::pair<NT, NT> >(nDim) { }
    };





    struct DefaultStoppingRule {
        template<typename T>
        bool operator()(T const&) const { return false; }
    };





    template<typename XT, typename YT>
    struct Result {
        typedef std::vector<XT> arg_type;
        typedef YT result_type;
        Result(void) {}
        Result(size_t xDim) : x(xDim) { }
        arg_type x;
        result_type y;
    };





    template<typename Particle, typename Bounds, typename Generator>
    void generateRandomPosition(Particle& particle, Bounds const& bounds, Generator& gen) {
        typedef typename Bounds::const_iterator BvIter;
        typedef typename Particle::arg_type::iterator PIter;
        PIter x(particle.x.begin());
        BvIter b(bounds.begin());
        for(; b!=bounds.end(); ++b, ++x) {
            *x = b->first + (b->second - b->first) * gen();
        }
    }





    template<typename ParticleStore, typename Bounds, typename Generator>
    void generateRandomParticles(ParticleStore& particles, Bounds const& bounds, Generator& gen) {
        typedef typename ParticleStore::value_type Particle;
        typedef typename ParticleStore::iterator PvIter;
        for(PvIter p(particles.begin()); p!=particles.end(); ++p) {
            *p = Particle(bounds.size());
            generateRandomPosition(*p, bounds, gen);
        }
    }





    template<typename NT>
    struct DefaultBoundsAdjustor {
        DefaultBoundsAdjustor(NT CMCR, NT PAR, NT jitter) 
          : CMCR_(CMCR), PAR_(PAR), jitter_(jitter) 
        {
            assert(0.0 <= cmcr && cmcr <= 1.0);
            assert(0.0 <=  par && par  <= 1.0);
            assert(jitter < 1.0);
        }

        template<typename Particle, typename Alg, typename Generator>
        void adjustBounds(Particle& particle, Alg const& alg, Generator& gen) const {
            if(gen() < CMCR_) {
                if(gen() < PAR_) {
                    // not sure if this is exactly what was meant... but it seems to make sense.
                    moveParticleInbounds(particle, alg.getBounds(), jitter_, gen);
                } else {
                    // use a previous good particle.
                    Particle const& p(alg.getBestParticleStore().getReplacementParticle(gen));
                    std::copy(p.x.begin(), p.x.end(), particle.x.begin());
                    moveParticleSlightly(particle, alg.getBounds(), jitter_, gen); // I don't think they did this bit.
                }
            } else {
                // select a new random position.
                generateRandomPosition(particle, alg.getBounds(), gen);
            }
        }

        template<typename Particle, typename Bounds>
        bool checkBounds(Particle const& particle, Bounds const& bounds) {
            typedef typename Bounds::const_iterator BvIter;
            typedef typename Particle::arg_type::const_iterator PIter;
            BvIter b(bounds.begin());
            PIter  p(particle.x.begin());
            for(; b!=bounds.end(); ++b, ++p) if(*p < b->first || *p > b->second) return false;
            return true;
        }

    private:
        // rule 7
        template<typename Particle, typename Bounds, typename Generator>
        void moveParticleSlightly(Particle& particle, Bounds const& bounds, NT jitter, Generator& gen) const {
            typedef typename Bounds::const_iterator BvIter;
            typedef typename Particle::arg_type::iterator PIter;
            BvIter bb(bounds.begin());
            PIter  pb(particle.x.begin());
            for(; bb!=bounds.end(); ++bb, ++pb) {
                *pb += (gen() - 0.5) * jitter;
                if(*pb < bb->first) *pb = bb->first + gen() * (bb->second - bb->first) * 0.5000001 * jitter;
                else if(*pb > bb->second) *pb = bb->second - gen() * (bb->second - bb->first) * 0.5000001 * jitter;
            }
        }

        // rule 7
        template<typename Particle, typename Bounds, typename Generator>
        void moveParticleInbounds(Particle& particle, Bounds const& bounds, NT jitter, Generator& gen) const {
            typedef typename Bounds::const_iterator BvIter;
            typedef typename Particle::arg_type::iterator PIter;
            BvIter bb(bounds.begin());
            PIter  pb(particle.x.begin());
            for(; bb!=bounds.end(); ++bb, ++pb) {
                if(*pb < bb->first) *pb = bb->first + gen() * (bb->second - bb->first) * jitter;
                else if(*pb > bb->second) *pb = bb->second - gen() * (bb->second - bb->first) * jitter;
            }
        }
        
    private:
        NT
            CMCR_,
            PAR_,
            jitter_;
    };





    template<typename Particle>
    struct BestStoreImpl {
        typedef Particle particle_type;

    private:
        struct BPC {
            bool operator()(particle_type const& a, particle_type const& b) const { 
                return a.y < b.y;
            }
        };
        typedef std::set<particle_type, BPC> BestVector;

    public:
        BestStoreImpl(size_t size) : maxSize_(size) {}

        bool addToStore(particle_type const& particle) {
            typedef BestVector::iterator Iter;
            bool inserted(false);
            if(bestVector_.size() < maxSize_) {
                bestVector_.insert(particle);
                inserted = true;
            }
            Iter insertAt(bestVector_.upper_bound(particle));
            if(!inserted && insertAt!=bestVector_.end()) {
                *insertAt = particle;
            }
            return insertAt==bestVector_.begin();
        }

        template<typename Generator>
        particle_type const& getReplacementParticle(Generator& gen) const {
            typedef BestVector::const_iterator Iter;
            size_t index(std::min((size_t)(gen() * bestVector_.size()), bestVector_.size()));
            Iter iter(bestVector_.begin());
            for(; index; --index) ++iter; 
            return *iter;
        }

        particle_type const& getBestParticle(void) const {
            return *bestVector_.begin();
        }

    private:
        size_t maxSize_;
        BestVector bestVector_;
    };





} // end namespace optim

#endif //OPTIM_UTILITIES_HEADER_INCLUDED_SDFPGJWPY8U60Y9WU850W49