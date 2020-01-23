#include "common.hpp"
#include <functional>

namespace crossmods {

using std::sqrt;
using std::erf;
using std::log;



struct NormalDistribution {
	double mean;
	double std;

	NormalDistribution(double mean=0.0, double std=1.0)
		: mean(mean), std(std)
	{}

	double cdf(double x) {
		double z = (x - mean)/(sqrt(2)*std);
		return (1.0 + erf(z))/2.0;
	}
};

struct LognormalDistribution {
	double mean;
	double std;

	LognormalDistribution(double mean=0.0, double std=1.0)
		: mean(mean), std(std)
	{}
	
	double cdf(double x) {
		if(x <= 0) {
			return 0.0;
		}
		double z = (log(x) - mean)/(sqrt(2)*std);
		return (1.0 + erf(z))/2.0;
	}
};

using ThresholdDist = LognormalDistribution;
using LatencyDist = LognormalDistribution;
using PassedDist = NormalDistribution;

auto always_one = [](size_t t) { return 1.0; };

constexpr double inf = std::numeric_limits<double>::infinity();
struct LognormalTdm {
	private:
	ThresholdDist th_d;
	LatencyDist latency_d;
	PassedDist passed_d;
	public:
	
	LognormalTdm(double thm, double ths, double lagm, double lags, double passm, double passs)
		:th_d(thm, ths), latency_d(lagm, lags), passed_d(passm, passs)
	{}

	CrossingPdf& decisions(CrossingPdf& out, invec taus) {
		decisions_cdf(out, taus, [](size_t t) { return 1.0; });
		double prev = 0.0;
		for(size_t t=0; t < out.grid.N; ++t) {
			double p = out.ps[t];
			out.ps[t] = (p - prev)/out.grid.dx;
			prev = p;
		}

		return out;
	}
	
	// TODO: Refactor this out
	CrossingPdf& blocker_decisions(CrossingPdf& out, invec taus, invec taus_b) {
		auto noearlycross = (*this);
		noearlycross.th_d.mean = inf;
		CrossingPdf unblocked(out.grid.N, out.grid.dx);
		noearlycross.decisions_cdf(unblocked, taus_b, [](size_t t) { return 1.0; });

		auto nolatecross = (*this);
		nolatecross.passed_d.mean = -inf;
		nolatecross.decisions_cdf(out, taus, [&](size_t t) { return unblocked.ps[t]; });

		double prev = 0.0;
		for(size_t t=0; t < out.grid.N; ++t) {
			double p = out.ps[t];
			out.ps[t] = (p - prev)/out.grid.dx;
			prev = p;
		}

		return out;
	}
	
	std::unique_ptr<CrossingPdf> blocker_decisions(invec taus, invec taus_b, double dt) {
		auto out = std::make_unique<CrossingPdf>(taus.size(), dt);
		blocker_decisions(*out, taus, taus_b);
		return out;
	}
	
	std::unique_ptr<CrossingPdf> decisions(invec taus, double dt) {
		auto out = std::make_unique<CrossingPdf>(taus.size(), dt);
		decisions(*out, taus);
		return out;
	}

	private:
	
	template <typename goddamn>
	void decisions_cdf(CrossingPdf& out, invec taus, goddamn unblocked_cdf) {
		// TODO: We don't actually have to grid this, especially
		// for simple trajectories
		auto ts = out.grid;
		auto dt = ts.dx;
		auto N = ts.N;

		double max_seen = -inf;
		double min_seen = inf;
		
		auto max_decided = 0.0;
		auto min_decided = 0.0;

		for(size_t t=0; t < N; ++t) {
			auto tau = taus[t];
			
			if(tau > max_seen) {
				max_seen = tau;
				max_decided = th_d.cdf(max_seen);
			}

			if(tau < min_seen) {
				min_seen = tau;
				min_decided = 1.0 - passed_d.cdf(min_seen);
			}
			
			auto unblocked = unblocked_cdf(t);
			auto total_decided = (max_decided + min_decided - max_decided*min_decided);
			total_decided *= unblocked;
			
			auto prev_p = 0.0; // = latency_d.cdf(0.0);
			for(size_t crossing_t=t; crossing_t < N; ++crossing_t) {
				auto lag = ts[crossing_t] - ts[t];
				auto p = latency_d.cdf(lag + dt);
				out.ps[crossing_t] += total_decided*(p - prev_p);
				prev_p = p;
			}
		}

		out.uncrossed = 1.0 - out.ps[N-1];
	}

};

}
