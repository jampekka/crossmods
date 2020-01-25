#include "common.hpp"
#include <functional>
#include <assert.h>

namespace crossmods {

using std::sqrt;
using std::erf;
using std::log;

const double erfinvhalf = 0.4769362762044698;

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
	
	// TODO: Use the real mean and std!
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
	double pass_th;
	public:
	
	LognormalTdm(double thm, double ths, double lagm, double lags, double pass_th)
		:th_d(thm, ths), latency_d(lagm, lags), pass_th(pass_th)
	{
	}

	CrossingPdf& decisions(CrossingPdf& out, invec taus, size_t start=0) {
		auto ts = out.grid;
		auto dt = ts.dx;
		auto N = ts.N;

		double max_seen = -inf;
		double decided = 0.0;

		for(size_t t=start; t < out.grid.N; ++t) {
			auto tau = taus[t];
			if(tau <= pass_th) tau = inf;
			if(tau <= max_seen) continue;
			max_seen = tau;
			auto deciders = th_d.cdf(max_seen) - decided;
			decided += deciders;
			
			auto prev_p = 0.0; // = latency_d.cdf(0.0);
			for(size_t crossing_t=t; crossing_t < N; ++crossing_t) {
				auto lag = ts[crossing_t] - ts[t];
				auto p = latency_d.cdf(lag + dt);
				
				out.ps[crossing_t] += deciders*(p - prev_p);
				prev_p = p;
			}
		}
		
		out.uncrossed = 1.0;
		for(auto& p : out.ps) {
			out.uncrossed -= p;
			p /= dt;
		}

		return out;
	}

		
	// TODO: Refactor this out
	CrossingPdf& blocker_decisions(CrossingPdf& out, invec taus, invec taus_b) {
		size_t start=0;

		for(; (start < out.grid.N) && (taus_b[start] > pass_th); ++start);
		auto nolatecross = (*this);
		nolatecross.pass_th = -inf;
		nolatecross.decisions(out, taus, start);

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
};

}
