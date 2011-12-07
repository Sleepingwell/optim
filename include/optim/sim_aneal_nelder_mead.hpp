/*
 * file sim_aneal_nelder_mead.hpp
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

#ifndef SIM_ANEAL_NELDER_MEAD_HEADER_INCLUDED_PWEOUPWSBOIJPIIJPHR
#define SIM_ANEAL_NELDER_MEAD_HEADER_INCLUDED_PWEOUPWSBOIJPIIJPHR

#include <cmath>
#include <cassert>
#include <numeric>

#include "optim_utilities.hpp"

#define GET_PSUM \
    for (n=0; n<ndim; ++n) {\
        for (sum=0.0,m=0; m<mpts; ++m) sum += p[m][n];\
        psum[n]=sum;\
    }

namespace optim {
    namespace detail {
        template<typename NT, typename Generator, typename Bounds>
        NT amotsa(
            Generator& gen,
            NT **p, 
            Bounds const& bounds, 
            NT y[], 
            NT psum[], 
            int ndim, 
            NT pb[], 
            NT *yb, 
            NT (*funk)(NT [], Bounds const&), 
            int ihi, 
            NT *yhi, 
            NT fac, 
            NT tt, 
            NT ptry[], 
            NT jitter
        ) {
            //Extrapolates by a factor fac through the face of the simplex across from the high point, tries
            //it, and replaces the high point if the new point is better.
            int j;
            NT fac1, fac2, yflu, ytry, tmpj;
            fac1=(1.0f-fac)/ndim;
            fac2=fac1-fac;
            for(j=0; j<ndim; ++j) {
                tmpj = psum[j]*fac1-p[ihi][j]*fac2;
                //if(tmpj < bounds[j].first) tmpj = bounds[j].first + gen()*jitter;
                //else if(tmpj > bounds[j].first) tmpj = bounds[j].second - gen()*jitter;
                ptry[j] = tmpj;
            }
            ytry=(*funk)(ptry, bounds);
            if(ytry <= *yb) { //Save the best-ever.
                for(j=0; j<ndim; ++j) pb[j]=ptry[j];
                *yb=ytry;
            }
            yflu=ytry-tt*log(gen()); //We added a thermal fluctuation to all the current vertices, 
                                      //but we subtract it here, so as to give the simplex a thermal 
                                      //Brownian motion: It likes to accept any suggested change.
            if(yflu < *yhi) {
                y[ihi]=ytry;
                *yhi=yflu;
                for(j=0; j<ndim; ++j) {
                    psum[j] += ptry[j]-p[ihi][j];
                    p[ihi][j]=ptry[j];
                    //std::cout << ptry[j] << ",";
                }
                //std::cout << ytry;
                //std::cout << '\n';
            }
            
            return yflu;
        }
    } // end namespace detail


    /*
    Multidimensional minimization of the function funk(x) where x[1..ndim] is a vector in
    ndim dimensions, by simulated annealing combined with the downhill simplex method of Nelder
    and Mead.

    - The input matrix p[1..ndim+1][1..ndim] has ndim+1 rows, each an ndimdimensional
    vector which is a vertex of the starting simplex.

    - the vector y[1..ndim+1], whose components must be pre-initialized to the values of funk evaluated
    at the ndim+1 vertices (rows) of p.

    - ftol, the fractional convergence tolerance to be achieved in the function value for an early return;

    - iter, and temptr. The routine makes iter function evaluations at an annealing temperature
    temptr, then returns. You should then crease temptr according to your annealing schedule, 
    reset iter, and call the routine again (leaving other arguments unaltered between calls). If
    iter is returned with a positive value, then early convergence and return occurred. If you 
    initialize yb to a very large value on the firt call then yb and pb[1..ndim] will 
    subsequently return the best function value and point ever encountered (even if it is 
    no longer a point in the simplex).
     */
    template<typename NT, typename Bounds>
    void amebsa(
        NT **p, 
        Bounds const& bounds, 
        NT y[], 
        int ndim, 
        NT pb[], 
        NT *yb, 
        NT ftol, 
        NT (*funk)(NT [], Bounds const&), 
        int *iter, 
        NT temptr, 
        NT jitter=1e-7f
    ) {
        UniformGenerator<NT> gen;
        int i,ihi,ilo,j,m,n,mpts(ndim+1);
        NT tt,rtol,sum,swap,yhi,ylo,ynhi,ysave,yt,ytry,*psum,*ptry,tmpj;
        bool recalc;
        ptry = new NT[ndim];
        psum = new NT[ndim];
        tt = -temptr;
        //Ensure all p are within bounds
        for(i=0; i<mpts; ++i) {
            recalc = false;
            for(j=0; j<ndim; ++j) {
                if(p[i][j] < bounds[j].first) {
                    recalc = true;
                    p[i][j] = bounds[j].first + gen()*jitter;
                } else if(p[i][j] > bounds[j].second) {
                    recalc = true;
                    p[i][j] = bounds[j].second - gen()*jitter;
                } 
                if(recalc) {
                    y[i] = (*funk)(p[i], bounds);
                }
            }
        }
        GET_PSUM
        for(;;) {
            ilo=0; //Determine which point is the highest (worst),
            ihi=1; //next-highest, and lowest (best)
            ynhi=ylo=y[0]+tt*log(gen()); //Whenever we look at a vertex, it gets a random thermal fluctuation.
            yhi=y[1]+tt*log(gen());
            if(ylo > yhi) {
                ihi=0;
                ilo=1;
                ynhi=yhi;
                yhi=ylo;
                ylo=ynhi;
            }
            for(i=2; i<mpts; ++i) { //Loop over the points in the simplex.
                yt=y[i]+tt*log(gen()); //More thermal fluctuations.
                if(yt <= ylo) {
                    ilo=i;
                    ylo=yt;
                }
                if(yt > yhi) {
                    ynhi=yhi;
                    ihi=i;
                    yhi=yt;
                } else if(yt > ynhi) {
                    ynhi=yt;
                }
            }
            rtol=2.0f*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo));
            //Compute the fractional range from highest to lowest and return if satisfactory.
            if(rtol < ftol || *iter < 0) {// If returning, put best point and value in slot 1.
                swap=y[0];
                y[0]=y[ilo];
                y[ilo]=swap;
                for(n=0; n<ndim; ++n) {
                    swap=p[0][n];
                    p[0][n]=p[ilo][n];
                    p[ilo][n]=swap;
                }
                break;
            }
            *iter -= 2;
            //Begin a new iteration. First extrapolate by a factor ?1 through the face of the simplex
            //across from the high point, i.e., reflect the simplex from the high point.
            ytry=detail::amotsa(gen,p,bounds,y,psum,ndim,pb,yb,funk,ihi,&yhi,-1.0,tt,ptry,jitter);
            if(ytry <= ylo) {
                //Gives a result better than the best point, so try an additional extrapolation by a factor of 2.
                ytry=detail::amotsa(gen,p,bounds,y,psum,ndim,pb,yb,funk,ihi,&yhi,2.0,tt,ptry,jitter);
            } else if(ytry >= ynhi) {
                //The reflected point is worse than the second-highest, so look for an intermediate
                //lower point, i.e., do a one-dimensional contraction.
                ysave=yhi;
                ytry=detail::amotsa(gen,p,bounds,y,psum,ndim,pb,yb,funk,ihi,&yhi,0.5,tt,ptry,jitter);
                if(ytry >= ysave) { //Can't seem to get rid of that high point. Better contract around the lowest (best) point.
                    for(i=0; i<mpts; ++i) {
                        if(i != ilo) {
                            for(j=0; j<ndim; ++j) {
                                tmpj = 0.5f*(p[i][j]+p[ilo][j]);
                                //if(tmpj < bounds[j].first) tmpj = bounds[j].first + gen()*jitter;
                                //else if(tmpj > bounds[j].second) tmpj = bounds[j].second - gen()*jitter;
                                psum[j]=tmpj;
                                p[i][j]=psum[j];
                            }
                            y[i]=(*funk)(psum, bounds);
                        }
                    }
                    *iter -= ndim;
                    GET_PSUM //Recompute psum.
                }
            } else ++(*iter);// Correct the evaluation count.
        }
        delete [] ptry;
        delete [] psum;
    }

} // end namespace optim

#endif //SIM_ANEAL_NELDER_MEAD_HEADER_INCLUDED_PWEOUPWSBOIJPIIJPHR
