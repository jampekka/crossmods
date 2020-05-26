#pragma once
#include "common.hpp"
#include <tuple>

namespace crossmods {
using std::atan;
using std::erf;
using std::erfc;
using std::sqrt;


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


struct Vddm {
	double dt;
	double std;
	double damping;
	double tau_threshold;
	double pass_threshold=0.0;
	double scale=1.0;
	double act_threshold=1.0;
	
	Vddm(
			double dt, // Time-step used in the simulation
			double std, // Standard deviation of the accmulation noise
			double damping, // The damping term. TODO: This is currently in weird units
			double tau_threshold, // The threshold when over which the input starts to accumulate positive activation
			double pass_threshold=0.0, // The threshold under which the input starts to accumulate positive activation
			double scale=1.0, // Scaling term for tau - tau_threshold
			double act_threshold=1.0 // Threshold at what activation level the decision is made
		)
		:dt(dt), std(std), damping(damping), tau_threshold(tau_threshold), pass_threshold(pass_threshold), scale(scale),
		act_threshold(act_threshold)
		{}
	
	double step(const Grid1d& acts, double tau, const double prev_weights[], double new_weights[], double decision_prob) const {
		auto alpha = 1.0 - exp(-dt*damping);
		auto da = acts.dx;
		auto N = acts.N;
		
		// pi/2.0 is the maximum for the atan nonlinearity
		auto diff_mean_tau = dt*pi/2.0;
		if(tau >= pass_threshold) {
			// If we're over the pass threshold, run the diffusion
			// "normally"
			diff_mean_tau = dt*atan(scale*(tau - tau_threshold));
		} // Otherwise use the maximum possible activation
		

		// This loop computes the variable drift diffusion. The equation is:
		// act[i] - act[i-1] = D[i] = tau[i]*dt - alpha*act[i-1] + e, e ~ N(0, dt*std**2)
		// For each outgoing bin, we approximate the bin's probability mass as the center
		// of the bin so act[i-1,from] is a scalar, meaning
		// D[i,from] ~ N(tau[i]*dt - alpha*act[i-1,from], dt*std**2)
		//           = N(diff_mean[from], dt*std**2)
		//
		// Thus, the proportion of probability mass from acts[from] to acts[to]
		// is
		// P(D[i,from] < acts[to] - acts[from] + da & D[i,from] > acts[to] - acts[from] - da)
		// = P(D[i,from] < acts[to] - acts[from] + da)*P(D[i,from] > acts[to] - acts[from] - da/2)
		// 
		// In the grid each endpoint is a startpoint of the next bin, so we can reuse the
		// upper edge probability of a previous bin as complement of the lower edge probability
		// of the next bin (too_small becomes small_enough).
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
		

		// After diffusing, we compute how big a share of the probability mass
		// is above the act_threshold, scale by decision_prob sum up these to
		// get share decided and remove the decided share from the distribution
		// leaving only the still undecided.
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

		// And finally normalize the distribution to sum to one.
		for(size_t i=0; i < N; ++i) {
			//new_weights[i] /= (1.0 - decided);
			new_weights[i] /= weightsum;
		}

		return decided;
	}
	
	auto step(const Grid1d& acts, double tau, invec prev_weights, double decision_prob) const {
		std::vector<double> new_weights(acts.N, 0.0);
		double decided = step(acts, tau, prev_weights.data(), new_weights.data(), decision_prob);
		return std::make_tuple(new_weights, decided);
	}

	private:
	
	// This is the actual implementation of the decision process. It's here as a
	// weird private template method, because C++ is a mess (in more detail, it's difficult
	// to pass a function as a parameter in a generic way without templates). This is private
	// because the Python binding generator can not handle templates.
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
	
	// TODO: Refactor these blocker hacks out of here
	double blocker_decisions(const Grid1d& acts, const double taus[], const double taus_b[], double decidedpdf[], size_t dur) const {
		// For the first vehicle, make a parameterization that will never cross before
		// the vehicle by putting the tau_threshold to infinity. This is because the
		// HIKER design where participants were instructed not to pass before the first
		// veicle.
		// NOTE: This does gather quite a bit of negative activation for the first vehicle,
		// which may make the late crosses to be quite late.
		// TODO: Try out a way to do this without accumulating the negative evidence
		auto noearlycross = (*this);
		noearlycross.tau_threshold = std::numeric_limits<double>::infinity();
		noearlycross.decisions(acts, taus_b, decidedpdf, dur);
		
		// Create a function to return share of decisions to cross after
		// the first vehicle.
		auto unblocked = 0.0;
		auto unblocked_share = [&](int t) {
			unblocked += decidedpdf[t]*dt;
			return unblocked;
		};
		
		// For the second vehicle, make a parameterization that never crosses
		// after the vehicle by putting the pass_threshold to infinity.
		// This is due to the HIKER experiment design where the participants
		// were not allowed to cross after the second car.
		// The share of decisions for the first vehicle is passed to the
		// second VDDM so that only the share that decided to cross after the
		// first vehicle can make a crossing decision.
		auto nolatecross = (*this);
		nolatecross.pass_threshold = -std::numeric_limits<double>::infinity();
		return nolatecross.decisions(acts, taus, decidedpdf, dur, unblocked_share);
	}

	// What follows are just overloads of the decisions methods to allow for
	// easier to handle parameter types.

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
