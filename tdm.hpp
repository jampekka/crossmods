#include "common.hpp"

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

	/*GenericTdm(ThresholdDist th_d, LatencyDist latency_d, PassedDist passed_d)
		:th_d(th_d), latency_d(latency_d), passed_d(passed_d)
		{}*/
	
	CrossingPdf& decisions(CrossingPdf& out, invec taus) {
		// TODO: We don't actually have to grid this, especially
		// for simple trajectories
		auto ts = out.grid;
		auto dt = ts.dx;
		auto N = ts.N;
		
		auto threshold_passed = vector<double>(N, 0.0);
		
		double max_seen = -inf;
		double min_seen = inf;
		
		auto max_decided = 0.0;
		auto min_decided = 0.0;

		double total_decided = 0.0;
		for(size_t t=0; t < N; ++t) {
			auto tau = taus[t];
			if(!(tau > max_seen or tau < min_seen)) continue;
			
			if(tau > max_seen) {
				max_seen = tau;
				max_decided = th_d.cdf(max_seen);
			}

			if(tau < min_seen) {
				min_seen = tau;
				min_decided = 1.0 - passed_d.cdf(min_seen);
			}
			
			auto decided = (max_decided + min_decided - max_decided*min_decided) - total_decided;
			total_decided += decided;

			auto prev_p = 0.0; // = latency_d.cdf(0.0);
			for(size_t crossing_t=t; crossing_t < N; ++crossing_t) {
				auto lag = ts[crossing_t] - ts[t];
				auto p = latency_d.cdf(lag + dt);
				out.ps[crossing_t] += (p - prev_p)*decided;
				prev_p = p;
			}
		}

		auto uncrossed = 1.0;
		for(size_t t=0; t < N; ++t) {
			uncrossed -= out.ps[t];
			out.ps[t] /= dt;
		}

		out.uncrossed = uncrossed;
		return out;
	}
	
	std::unique_ptr<CrossingPdf> decisions(invec taus, double dt) {
		auto out = std::make_unique<CrossingPdf>(taus.size(), dt);
		decisions(*out, taus);
		return out;
	}
};

}
