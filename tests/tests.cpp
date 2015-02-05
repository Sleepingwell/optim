#include <limits>
#include <vector>
#include <utility>

#include "../include/optim/sim_aneal_nelder_mead.hpp"
#include "../include/optim/charged_particle_search.hpp"

template<typename Bounds>
double func(double x[], Bounds const& bounds) {
    if(x[0] < bounds[0].first || x[0] > bounds[0].second || x[1] < bounds[1].first || x[1] > bounds[1].second) return std::numeric_limits<double>::infinity();
    //return sin(100*x[0])*exp(x[0]) + sin(100*x[1])*exp(x[1]);
    return x[0]*x[0] + x[1]*x[1];
}

double newTemp(double f) {
    return 0.99f * f;
}

typedef std::pair<size_t, double> Res;

Res testSimAnealNelderMead(int nIters) {
    typedef optim::Bounds<double> Bounds;
    optim::UniformGenerator<double> gen;
    int ndim(2), nloops(0), i, j, iter(nIters), totIters(0);
    double* pbase = new double[(ndim+1) * ndim];
    double** p = new double*[ndim+1];
    double* boundsBase = new double[ndim*3];
    Bounds bounds(ndim);
    for(i=0; i<ndim; ++i) {
        bounds[i].first = -1.0;
        bounds[i].second = 1.0;
        p[i] = &pbase[ndim*i];
    }
    p[i] = &pbase[ndim*i];
    for(i=0; i<=ndim; ++i) for(j=0; j<ndim; ++j) p[i][j] = bounds[j].first + gen()*(bounds[j].second - bounds[j].first);
    std::vector<double> y(ndim+1);
    y[0] = func(p[0], bounds); y[1] = func(p[1], bounds); y[2] = func(p[2], bounds);
    std::vector<double> pb(ndim);
    double yb(std::numeric_limits<double>::max());
    double ftol(1e-9f);
    double temptr(10.0);
    while(temptr > 0.0001f) {
        optim::amebsa(p, bounds, &y[0], ndim, &pb[0], &yb, ftol, &func<Bounds>, &iter, temptr);
        totIters += nIters - iter;
        temptr = newTemp(temptr);
        iter = nIters;
        ++nloops;
    }
    delete [] p;
    delete [] pbase;
    return Res(totIters, yb);
}

double func2(double x[], size_t ndim) {
    //return sin(100*x[0])*exp(x[0]) + sin(100*x[1])*exp(x[1]);
         //sin(100* x1 )*exp( x1 ) + sin(100* x2 )*exp( x2 )
    return x[0]*x[0] + x[1]*x[1];
} 

template<typename FUNC>
Res testChargedParticles(FUNC& func, size_t nParticles, size_t nIters) {
    typedef typename optim::Result<double, double> Result;
    typedef optim::Bounds<double> Bounds;

    int
        i,
        ndim(2),
        seed(42);

    Bounds
        bounds(ndim);

    for(i=0; i<ndim; ++i) {
        bounds[i].first = -1.0;
        bounds[i].second = 1.0;
    }
    Result res(chargedParticleOptimise<double, double, double>(func, bounds, nIters, nParticles, seed));
    return Res(nParticles * nIters, res.y);
}

struct Tester {
    double operator()(double x[], size_t ndim) {
        return func2(x, ndim);
    }
};

int main(int argc, char* argv[]) {
    //Res t1res = testSimAnealNelderMead(50);
    size_t
        nParticles(5);

    //size_t nIters(t1res.first/nParticles);
    Res t2res = testChargedParticles(Tester(), nParticles, 5);
    return 0;
}
