#include "crossmods.hpp"
#include <iostream>
#include <nlohmann/json.hpp>
#include <chrono>
#include <omp.h>
typedef std::chrono::high_resolution_clock Clock;

using json = nlohmann::json;
using namespace crossmods;

int main() {
	double dt = 1/30.0;
	Vddm vddm = {
		.dt=dt,
		.std=0.75,
		.damping=1.6,
		.tau_threshold=2.3,
		.scale=1.0,
	};

	LognormalTdm tdm(1.0, 1.0, 1.0, 1.0, 0.0);
	
	size_t dur = 10.0/dt;
	Grid1d grid({-3.0, 3.0, 100});
	CrossingPdf pdf(dur, dt);
	
	double tau0 = 3.0;
	//double speed = 20.0;
	vector<double> tau(dur);
	std::generate(tau.begin(), tau.end(),
			[&, t=0]() mutable { return tau0 - t++*dt;}
			);
	
	tdm.decisions(tau, dt);

	size_t niter = 10;
	size_t ntrials = 6;
	// OMG C++!!
	vddm.decisions(grid, pdf, tau);
	
	auto t0 = Clock::now();
	for(size_t i=0; i < niter; ++i) {
		for(size_t t=0; t < ntrials; ++t) {
			CrossingPdf pdf(dur, dt);
			vddm.decisions(grid, pdf, tau);
		}
	}

	std::chrono::duration<double, std::milli> fuckyou(Clock::now() - t0);
	std::cerr << fuckyou.count()/(niter*ntrials) << std::endl;
	std::cout << json(pdf.ps) << std::endl;
}
