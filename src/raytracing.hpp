#ifndef RAYTRACING_H_
#define RAYTRACING_H_

#include <omp.h>
#include "particle.hpp"
#include "physics.hpp"
#include "metric_Kerr_Schild.hpp"
#include "metric_KS_Schwarzschild.hpp"

/**
* Raytracing module
* Requires a metric and a disk model
*/
template <class T, class Model>
class Raytracer
{
	T _M;   // Metric to be specified when instanciating the class
	Model _disk;	
	
public:
	Particles _photons;

        Raytracer(const T &M, const Particles &particles): _M(M), _photons(particles) 
	{}

	Raytracer(const T &M, const Particles &particles, const disk_model & Mod): _M(M), _photons(particles), _disk(Mod) 
	{}

        // Compute derivatives needed for geodesic integration
	Particle derivatives(const Particle particle, double epsilon, double dt);
	
	// Integrate all particles positions by one step using Runge-Kutta 4 integration 
	void auto_integrate_all_once_RK4(double epsilon, double dt, double dmax);
	
	// Integrate each particle position all along their trajectory using simple euler integration
	void auto_integrate_photon_KS_Euler(double epsilon, double dt, double dmax, int steps);

	// Integrate each particle position all along their trajectory using Runge-Kutta 4 integration
	// and taking into account some physical phenomenon: wavelength of emission and doppler effect
	void auto_integrate_photon_KS_RK4(double epsilon, double dt, double dmax, int steps);
	
	Particles get_particles(){return _photons;}

};

// force compiler to compile Raytracer for your metric
template class Raytracer<KerrSchild_Schwarzschild,disk_model>;
template class Raytracer<KerrSchild,disk_model>;

#endif
