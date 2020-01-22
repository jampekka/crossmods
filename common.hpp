#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <math.h>

namespace crossmods {
using std::vector;
using invec = const std::vector<double>&;

struct Grid1d {
	double _low;
	double dx;
	size_t N;

	Grid1d(double low, double high, size_t N)
		:_low(low), dx((high - low)/N), N(N)
	{}

	size_t bin(double x) const {
		return (x - _low)/dx;
	}
	
	double low() {
		return _low;
	}

	double high() const {
		return _low + N*dx;
	}

	double operator[](size_t i) const {
		return _low + i*dx;
	}
};

struct CrossingPdf {
	Grid1d grid;
	vector<double> ps;
	double uncrossed = 1.0;
	
	CrossingPdf(size_t dur, double dt)
		:grid(0.0, dur*dt, dur), ps(dur, 0.0)
	{}
	
	double operator()(double ct) {
		if(ct < grid.low() || ct > grid.high()) {
			return nan("");
		}
		
		// TODO: Linear interpolation!

		return ps[grid.bin(ct)];
	}

	double loglikelihood(invec cts, double slack=0.0) {
		double loglik = 0.0;
		for(auto ct : cts) {
			if(ct > grid.high()) {
				loglik += std::log(uncrossed + slack);
			} else {
				loglik += std::log((*this)(ct) + slack);
			}
		}

		return loglik;
	}
};

}
