#include "common.hpp"
#include <boost/math/distributions.hpp>

namespace crossmods {

using boost::math::cdf;
using boost::math::normal;
using boost::math::lognormal;

using ThresholdDist = lognormal;
using LatencyDist = lognormal;
using PassedDist = normal;

constexpr double inf = std::numeric_limits<double>::infinity();
struct LognormalTdm {
	private:
	ThresholdDist th_d;
	LatencyDist latency_d;
	PassedDist passed_d;
	double pass_threshold=0.0;
	public:
	
	LognormalTdm(double thm, double ths, double lagm, double lags, double passm, double passs)
		:th_d(thm, ths), latency_d(lagm, lags), passed_d(passm, passs)
	{}

	/*GenericTdm(ThresholdDist th_d, LatencyDist latency_d, PassedDist passed_d)
		:th_d(th_d), latency_d(latency_d), passed_d(passed_d)
		{}*/
	
	CrossingPdf& decisions(CrossingPdf& out, invec taus) const {
		auto ts = out.grid;
		auto dt = ts.dx;
		auto N = ts.N;
		
		auto threshold_passed = vector<double>(N, 0.0);
		
		double max_seen = -inf;
		double min_seen = inf;
		double total_decided = 0.0;
		for(size_t t=0; t < N; ++t) {
			max_seen = std::max(max_seen, taus[t]);
			min_seen = std::min(min_seen, taus[t]);
			
			// TODO: Get rid of the damn boost!
			auto max_decided = cdf(th_d, std::max(0.0, max_seen));
			auto min_decided = 1.0 - cdf(passed_d, min_seen);
			
			auto decided = (max_decided + min_decided - max_decided*min_decided) - total_decided;
			total_decided += decided;

			
			for(size_t crossing_t=t; crossing_t < N; ++crossing_t) {
				auto lag = ts[crossing_t] - ts[t];
				out.ps[crossing_t] += (cdf(latency_d, lag + dt) - cdf(latency_d, lag))*decided;
			}
		}

		auto total_crossed = 0.0;
		for(size_t t=0; t < N; ++t) {
			total_crossed += out.ps[t];
			out.ps[t] /= dt;
		}

		out.uncrossed = 1.0 - total_crossed;
		return out;
	}
	
	std::unique_ptr<CrossingPdf> decisions(invec taus, double dt) const {
		auto out = std::make_unique<CrossingPdf>(taus.size(), dt);
		decisions(*out, taus);
		return out;
	}
};

}
