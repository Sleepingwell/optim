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
#include <algorithm>
#include <iterator>
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
        typedef XT position_element;
        typedef std::vector<position_element> position_type;
        typedef typename position_type::iterator position_iterator;
        typedef typename position_type::const_iterator const_position_iterator;
        typedef YT result_type;

        Result(size_t nArgs)
            : position_(nArgs) {
        }

        virtual ~Result(void) {
        }

        void swap(Result& other) {
            this->position_.swap(other.position_);
            std::swap(this->y, other.y);
        }

        position_iterator beginPosition(void) {
            return position_.begin();
        }

        position_iterator endPosition(void) {
            return position_.end();
        }

        const_position_iterator beginPosition(void) const {
            return position_.begin();
        }

        const_position_iterator endPosition(void) const {
            return position_.end();
        }

        XT getPositionElement(size_t index) const {
            return position_[index];
        }

        position_element* positionData(void) {
            return &position_[0];
        }

        position_element const* positionData(void) const {
            return &position_[0];
        }

        result_type
            y;

    protected:
        position_type
            position_;
    };





    template<typename BoundsIter, typename PositionIter, typename Generator>
    void generateRandomPosition(
        BoundsIter bBounds,
        BoundsIter eBounds,
        PositionIter bPosition,
        Generator& gen
    ) {
        for(; bBounds!=eBounds; ++bBounds, ++bPosition) {
            *bPosition = bBounds->first + (bBounds->second - bBounds->first) * gen();
        }
    }





    template<typename ParticleStore, typename Bounds, typename Generator>
    void generateRandomParticles(
        ParticleStore& particles,
        Bounds const& bounds,
        Generator& gen
    ) {
        typedef typename ParticleStore::value_type Particle;
        typedef typename ParticleStore::iterator PvIter;

        for(PvIter p(particles.begin()); p!=particles.end(); ++p) {
            generateRandomPosition(bounds.begin(), bounds.end(), p->beginPosition(), gen);
        }
    }





    template<typename NT>
    struct DefaultBoundsAdjustor {
        DefaultBoundsAdjustor(NT CMCR, NT PAR, NT jitter) 
          : CMCR_(CMCR), PAR_(PAR), jitter_(jitter) 
        {
            assert(0.0 <= CMCR && CMCR <= 1.0);
            assert(0.0 <=  PAR && PAR  <= 1.0);
            assert(jitter < 1.0);
        }



        template<typename PositionIter, typename Bounds>
        static bool checkBoundsOK(
            PositionIter posIter,
            Bounds const& bounds
        ) {
            typedef typename Bounds::const_iterator BvIter;

            BvIter
                b(bounds.begin()),
                e(bounds.end());

            for(; b!=e; ++b, ++posIter) if(*posIter < b->first || *posIter > b->second) return false;

            return true;
        }



        template<typename Particle, typename PositionIter, typename Alg, typename Generator>
        void adjustBounds(
            Particle const& particle,
            PositionIter posIter,
            Alg const& alg,
            Generator& gen
        ) const {
            if(gen() < CMCR_) {
                if(gen() < PAR_) {
                    // not sure if this is exactly what was meant... but it seems to make sense.
                    moveParticleInbounds(posIter, alg.getBounds(), jitter_, gen);
                } else {
                    // use a previous good particle.
                    Particle const&
                        p(alg.getReplacementParticle());

                    std::copy(p.beginPosition(), p.endPosition(), posIter);
                    moveParticleSlightly(posIter, alg.getBounds(), jitter_, gen); // I don't think they did this bit.
                }
            } else {
                // select a new random position.
                generateRandomPosition(alg.getBounds().begin(), alg.getBounds().end(), posIter, gen);
            }
        }



    private:
        // rule 7
        template<typename PositionIter, typename Bounds, typename Generator>
        static void moveParticleSlightly(
            PositionIter posIter,
            Bounds const& bounds,
            NT jitter,
            Generator& gen
        ) {
            typedef typename Bounds::const_iterator BvIter;

            BvIter
                bb(bounds.begin());

            for(; bb!=bounds.end(); ++bb, ++posIter) {
                *posIter += (gen() - 0.5) * jitter;
                if(*posIter < bb->first) *posIter = bb->first + gen() * (bb->second - bb->first) * 0.5000001 * jitter;
                else if(*posIter > bb->second) *posIter = bb->second - gen() * (bb->second - bb->first) * 0.5000001 * jitter;
            }
        }

        // rule 7
        template<typename PositionIter, typename Bounds, typename Generator>
        static void moveParticleInbounds(
            PositionIter posIter,
            Bounds const& bounds,
            NT jitter,
            Generator& gen
        ) {
            typedef typename Bounds::const_iterator BvIter;

            BvIter
                bb(bounds.begin());

            for(; bb!=bounds.end(); ++bb, ++posIter) {
                if(*posIter < bb->first) *posIter = bb->first + gen() * (bb->second - bb->first) * jitter;
                else if(*posIter > bb->second) *posIter = bb->second - gen() * (bb->second - bb->first) * jitter;
            }
        }
        
    private:
        NT
            CMCR_,
            PAR_,
            jitter_;
    };





    //template<typename Particle>
    //struct BestStoreImpl {
    //    typedef Particle particle_type;

    //private:
    //    struct BPC {
    //        bool operator()(particle_type const& a, particle_type const& b) const { 
    //            return a.y > b.y;
    //        }
    //    };
    //    typedef std::set<particle_type, BPC> BestVector;



    //public:
    //    BestStoreImpl(size_t size) : maxSize_(size) {}



    //    bool addToStore(particle_type const& particle) {
    //        typedef typename BestVector::iterator Iter;

    //        bool
    //            inserted(false);

    //        if(bestVector_.size() < maxSize_) {
    //            bestVector_.insert(particle);
    //            inserted = true;
    //        }

    //        Iter
    //            insertAt(bestVector_.lower_bound(particle));

    //        bool
    //            isBest(insertAt == bestVector_.begin());

    //        if(!inserted && insertAt!=bestVector_.end()) {
    //            // remember, set iterators are const!
    //            // *insertAt = particle;
    //            bestVector_.erase(*insertAt);
    //            bestVector_.insert(particle);
    //        }

    //        return isBest;
    //    }



    //template<typename StoreIter, typename Generator>
    //particle_type const& getReplacementParticle(
    //    StoreIter b,
    //    StoreIter e,
    //    Generator& gen
    //) const {
    //    typedef typename BestVector::const_iterator Iter;

    //    size_t
    //        nInStore(std::distance(b, e)),
    //        index(static_cast<size_t>(gen() * nInStore));

    //    assert(nInStore > 0);

    //    if(index == nInStore) {
    //        // can only happen if we got 1.0 from gen()
    //        index = nInStore - 1;
    //    }

    //    return *(b + index);
    //}



        //particle_type const& getBestParticle(void) const {
        //    return *bestVector_.begin();
        //}



    //private:
    //    size_t
    //        maxSize_;

    //    BestVector
    //        bestVector_;
    //};





} // end namespace optim

#endif //OPTIM_UTILITIES_HEADER_INCLUDED_SDFPGJWPY8U60Y9WU850W49