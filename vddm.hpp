#pragma once
#include <cmath>
#include <math.h>
#include <memory>
#include <vector>
#include <algorithm>
#include <functional>
//#include "fastonebigheader.h"

namespace crossmods {
using std::atan;
using std::erf;
using std::erfc;
using std::sqrt;
using std::vector;

using invec = const std::vector<double>&;

/*
inline double fastPow(double a, double b) {
    union {
        double d;
        int x[2];
    } u = { a };
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    return u.d;
}


inline double
fastererfc (double x)
{
  static const double k = 3.3509633149424609;

  return 2.0 / (1.0 + std::pow(2.0, k * x));
}
*/

// Come on C++ people!!
static const double pi = atan(1.0)*4;
inline double stdnormcdf(double x) {
	return (1.0 + erf(x))/2.0;
	//return 1.0 - fastererfc(x)/2.0;
}

inline double normcdf(double x, double m, double v) {
	double scaler = sqrt(2.0*v);
	double z = (x - m)/scaler;
	return stdnormcdf(z);
}

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

	double loglikelihood(invec cts) {
		double loglik = 0.0;
		for(auto ct : cts) {
			if(ct > grid.high()) {
				loglik += uncrossed;
			} else {
				loglik += std::log((*this)(ct));
			}
		}

		return loglik;
	}
};

struct Vddm {
	double dt;
	double std;
	double damping;
	double tau_threshold;
	double pass_threshold=0.0;
	double scale=1.0;
	double act_threshold=1.0;

	Vddm(double dt, double std, double damping, double tau_threshold, double pass_threshold=0.0, double scale=1.0, double act_threshold=1.0)
		:dt(dt), std(std), damping(damping), tau_threshold(tau_threshold), pass_threshold(pass_threshold), scale(scale),
		act_threshold(act_threshold)
		{}
	
	double step(const Grid1d& acts, double tau, const double prev_weights[], double new_weights[], double decision_prob) const {
		auto alpha = 1.0 - exp(-dt*damping);
		auto da = acts.dx;
		auto N = acts.N;

		auto diff_mean_tau = dt*pi/2.0;
		if(tau >= pass_threshold) {
			diff_mean_tau = dt*atan(scale*(tau - tau_threshold));
		}
		
		// TODO: There may be some unutilized symmetry in the normcdf
		#pragma omp parallel for reduction(+:new_weights[:N])
		for(size_t from=0; from < N; ++from) {
			auto diff_mean = diff_mean_tau - alpha*acts[from];
			double too_small = 0.0;
			for(size_t to=0; to < N - 1; ++to) {
				auto diff = acts[to] - acts[from] + da/2.0;
				double small_enough = normcdf(diff, diff_mean, dt*std*std);
				new_weights[to] += (small_enough - too_small)*prev_weights[from];
				too_small = small_enough;
			}
			new_weights[N-1] += (1.0 - too_small)*prev_weights[from];
		}
		
		// TODO: No need to do the whole loop, but one more loop here
    		// doesn't matter so much
		double decided = 0.0;
		double weightsum = 0.0;
		for(size_t i=0; i < N; ++i) {
			auto share_over = (acts[i] - act_threshold)/da + 0.5;
			share_over = share_over < 0.0 ? 0.0 : share_over;
			share_over = share_over > 1.0 ? 1.0 : share_over;
			auto bindec = decision_prob*share_over*new_weights[i];
			new_weights[i] -= bindec;
			decided += bindec;
			weightsum += new_weights[i];
		}

		for(size_t i=0; i < N; ++i) {
			//new_weights[i] /= (1.0 - decided);
			new_weights[i] /= weightsum;
		}

		return decided;
	}
	
	private:
	
	template <typename callback>
	double decisions(const Grid1d& acts, const double taus[], double decidedpdf[], size_t dur, callback cross_prob) const {
		double undecided = 1.0;
		std::vector<double> prev_weights(acts.N, 0.0);
		prev_weights[acts.bin(0.0)] = 1.0;
		std::vector<double> weights(acts.N, 0.0);

		for(size_t t=0; t < dur; ++t) {
			auto cb = cross_prob(t);
			double decided = undecided*step(acts, taus[t], prev_weights.data(), weights.data(), cb);
			prev_weights.swap(weights);
			std::fill(weights.begin(), weights.end(), 0.0);

			undecided -= decided;
			decidedpdf[t] = decided/dt;
		}

		return undecided;
	}

	public:

	double decisions(const Grid1d& acts, const double taus[], double decidedpdf[], size_t dur) const {
		return decisions(acts, taus, decidedpdf, dur, [](size_t t) { return 1.0; });
	}

	double blocker_decisions(const Grid1d& acts, const double taus[], const double taus_b[], double decidedpdf[], size_t dur) const {
		auto has_passed = [&](size_t t) { return taus_b[t] <= pass_threshold ? 1.0 : 0.0; };
		decisions(acts, taus_b, decidedpdf, dur, has_passed);

		auto unblocked = 0.0;
		auto unblocked_share = [&](int t) {
			unblocked += decidedpdf[t]*dt;
			return unblocked;
		};

		auto nolatecross = (*this);
		nolatecross.pass_threshold = 0.0;
		return nolatecross.decisions(acts, taus, decidedpdf, dur, unblocked_share);
	}

	CrossingPdf& decisions(const Grid1d& grid, CrossingPdf& out, invec taus) const {
		out.uncrossed = decisions(grid, taus.data(), out.ps.data(), out.ps.size());
		return out;
	}
	
	std::unique_ptr<CrossingPdf> decisions(const Grid1d& grid, invec taus) const {
		auto out = std::make_unique<CrossingPdf>(taus.size(), dt);
		decisions(grid, *out, taus);
		return out;
	}
	
	CrossingPdf& blocker_decisions(const Grid1d& grid, CrossingPdf& out, invec taus, invec taus_b) const {
		out.uncrossed = blocker_decisions(grid, taus.data(), taus_b.data(), out.ps.data(), out.ps.size());
		return out;
	}
	
	std::unique_ptr<CrossingPdf> blocker_decisions(const Grid1d& grid, invec taus, invec taus_b) const {
		auto out = std::make_unique<CrossingPdf>(taus.size(), dt);
		blocker_decisions(grid, *out, taus, taus_b);
		return out;
	}

};

}
