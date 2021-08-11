#include "raytracing.hpp"
#include <iostream>



template <class T, class Model>
void Raytracer<T, Model>::auto_integrate_photon_KS_RK4(double epsilon, double dt, double dmax, int steps){
        params par = _M.get_params();
	double M = autodiff::val(par.M);
	double a = autodiff::val(par.a);
	double rs = M + sqrt( M*M - a*a);
	float count = 0;
        int nb_part = _photons._particles.size();
#pragma omp parallel for num_threads(6)
        for(int i=0; i<nb_part; i++)
	{
	        #pragma omp atomic
	        count++;
	        
	        std::cout<<"Progress: "<<100*count/nb_part<<"%"<<std::endl;
	        for(int s=0; s<steps; s++){
	                
	                double z0 = _photons._particles[i].z;
	                Particle p0 = _photons._particles[i];
	                // RK4 integration
	                auto k1 = derivatives(p0, epsilon, _photons._particles[i].t)*dt;
	                auto k2 = derivatives(p0+k1/2., epsilon, _photons._particles[i].t+dt/2.)*dt;
	                auto k3 = derivatives(p0+k2/2., epsilon, _photons._particles[i].t+dt/2.)*dt;
	                auto k4 = derivatives(p0+k3, epsilon, _photons._particles[i].t+dt)*dt;
	                
	                Particle p;
	                p = (k1+k2*2.+k3*2.+k4)*(1./6.);
	                _photons._particles[i] = _photons._particles[i] + p;
	                _photons._particles[i].t += dt;
			
			// check photon position
			double x = _photons._particles[i].x;
			double y = _photons._particles[i].y;
			double z = _photons._particles[i].z;
			double t = _photons._particles[i].t;

			autodiff::dual rc = _M.r(x,y,z,t);
			_photons._particles[i].r = val(rc);
			
			// too far away or captured by BH
			if( _photons._particles[i].r > dmax || _photons._particles[i].r < rs ) break;
			// If collision with disk, change intensity
			if(z*z0 <= 0.){
			if(_photons._particles[i].r > _disk.inner_limit && _photons._particles[i].r < _disk.outer_limit){
				_photons._particles[i].intensity += _disk.emitted_light(_photons._particles[i],p,_disk.inner_limit,_disk.outer_limit);
				if(!_disk.thin) break;
			}
			}
	        }
	}
}

// Integrate each particle position all along their trajectory using simple euler integration
template <class T, class Model>
void Raytracer<T, Model>::auto_integrate_photon_KS_Euler(double epsilon, double dt, double dmax, int steps){
        params par = _M.get_params();
	double M = autodiff::val(par.M);
	double a = autodiff::val(par.a);
	double rs = M + sqrt( M*M - a*a);
	int nb_part = _photons._particles.size();
	int count = 0;
#pragma omp parallel for num_threads(6)
        for(int i=0; i<_photons._particles.size(); i++)
	{
	        #pragma omp atomic
	        count++;
                std::cout<<"Progress: "<<100*count/nb_part<<"%"<<std::endl;
	        for(int s=0; s<steps; s++){
	                // basic Euler integration
	                double z0 = _photons._particles[i].z;
	                auto k = derivatives(_photons._particles[i], epsilon, _photons._particles[i].t)*dt;
	                _photons._particles[i] = _photons._particles[i] + k;
	                _photons._particles[i].t += dt;
			
			autodiff::dual x = _photons._particles[i].x;	
			autodiff::dual y = _photons._particles[i].y;
			autodiff::dual z = _photons._particles[i].z;
			autodiff::dual t = _photons._particles[i].t;

			autodiff::dual rc = _M.r(x,y,z,t);
			_photons._particles[i].r = (val(rc));
			
			// too far away or captured by BH
			if( _photons._particles[i].r > dmax || _photons._particles[i].r < rs ) break;
			// If collision with disk, change intensity
			if(z*z0 < 0.){
			if(_photons._particles[i].r > _disk.inner_limit && _photons._particles[i].r < _disk.outer_limit){
				_photons._particles[i].intensity += _disk.emitted_light(_photons._particles[i],k,_disk.inner_limit,_disk.outer_limit);
				if(!_disk.thin) break;
			}
			}	
	        }
	}
}



template <class T, class Model>
void Raytracer<T, Model>::auto_integrate_all_once_RK4(double epsilon, double dt, double dmax){
        params par = _M.get_params();
	double M = autodiff::val(par.M);
	double a = autodiff::val(par.a);
	double rs = M + sqrt( M*M - a*a);
#pragma omp parallel for
	for(int i=0; i<_photons._particles.size(); i++)
	{
		if(!_photons._particles[i].jump)
		{
			// basic Euler integration
	                double z0 = _photons._particles[i].z;
	                auto k = derivatives(_photons._particles[i], epsilon, dt)*dt;
	                _photons._particles[i] = _photons._particles[i] + k;
	                _photons._particles[i].t += dt;
			
			autodiff::dual x = _photons._particles[i].x;	
			autodiff::dual y = _photons._particles[i].y;
			autodiff::dual z = _photons._particles[i].z;
			autodiff::dual t = _photons._particles[i].t;

			autodiff::dual rc = _M.r(x,y,z,t);
			_photons._particles[i].r = (val(rc));
			
			// too far away or captured by BH
			if( _photons._particles[i].r > dmax || _photons._particles[i].r < rs ) _photons._particles[i].jump = 1;
			// If collision with disk, change intensity
			if(z*z0 <= 0.){
			if(_photons._particles[i].r > _disk.inner_limit && _photons._particles[i].r < _disk.outer_limit){
				_photons._particles[i].intensity += _disk.emitted_light(_photons._particles[i],k,_disk.inner_limit,_disk.outer_limit);
				if(!_disk.thin) _photons._particles[i].jump = 1;
			}
			}
		}
	}
	
}

template <class T, class Model>
Particle Raytracer<T, Model>::derivatives(const Particle particle, double epsilon, double dt)
{
	autodiff::dual x = particle.x;	
	autodiff::dual y = particle.y;
	autodiff::dual z = particle.z;
	autodiff::dual t = particle.t;

	autodiff::dual aa = _M.alpha(x,y,z,t);
	autodiff::dual B_1 = _M.B1(x,y,z,t);
	autodiff::dual B_2 = _M.B2(x,y,z,t);
	autodiff::dual B_3 = _M.B3(x,y,z,t);

	autodiff::dual g_11 = _M.g11(x,y,z,t);
	autodiff::dual g_22 = _M.g22(x,y,z,t);
	autodiff::dual g_33 = _M.g33(x,y,z,t);
	autodiff::dual g_12 = _M.g12(x,y,z,t);
	autodiff::dual g_13 = _M.g13(x,y,z,t);
	autodiff::dual g_23 = _M.g23(x,y,z,t);

	// Alpha derivatives
	double d1_alpha = derivative([&](auto x, auto y,auto z, auto t) { return _M.alpha(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2_alpha = derivative([&](auto x, auto y,auto z, auto t) { return _M.alpha(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3_alpha = derivative([&](auto x, auto y,auto z, auto t) { return _M.alpha(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));
	
	// Beta derivatives
	double d1B1 = derivative([&](auto x, auto y,auto z, auto t) { return _M.B1(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2B1 = derivative([&](auto x, auto y,auto z, auto t) { return _M.B1(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3B1 = derivative([&](auto x, auto y,auto z, auto t) { return _M.B1(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));

	double d1B2 = derivative([&](auto x, auto y,auto z, auto t) { return _M.B2(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2B2 = derivative([&](auto x, auto y,auto z, auto t) { return _M.B2(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3B2 = derivative([&](auto x, auto y,auto z, auto t) { return _M.B2(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));

	double d1B3 = derivative([&](auto x, auto y,auto z, auto t) { return _M.B3(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2B3 = derivative([&](auto x, auto y,auto z, auto t) { return _M.B3(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3B3 = derivative([&](auto x, auto y,auto z, auto t) { return _M.B3(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));

	//gij derivatives
	double d1_g11 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g11(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2_g11 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g11(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3_g11 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g11(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));

	double d1_g12 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g12(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2_g12 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g12(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3_g12 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g12(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));
	
	double d1_g13 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g13(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2_g13 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g13(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3_g13 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g13(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));

	double d1_g22 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g22(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2_g22 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g22(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3_g22 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g22(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));

	double d1_g23 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g23(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2_g23 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g23(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3_g23 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g23(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));

	double d1_g33 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g33(x,y,z,t); }, autodiff::forward::wrt(x), at(x,y,z,t));
	double d2_g33 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g33(x,y,z,t); }, autodiff::forward::wrt(y), at(x,y,z,t));
	double d3_g33 = derivative([&](auto x, auto y,auto z, auto t) { return _M.g33(x,y,z,t); }, autodiff::forward::wrt(z), at(x,y,z,t));

//////////////////////////////////////

	auto u0 = sqrt(g_11*particle.v1*particle.v1 + 2.0*g_12*particle.v1*particle.v2 + 2.0*g_13*particle.v1*particle.v3
			+ g_22*particle.v2*particle.v2 + 2.0*g_23*particle.v2*particle.v3 + g_33*particle.v3*particle.v3 + epsilon)/aa;
	
//////////////////////////////////////

	auto u02 = 2.0*u0;
	/////////////
	auto du1 = -aa*u0*d1_alpha + particle.v1*d1B1 + particle.v2*d1B2 + particle.v3*d1B3 
		     -particle.v1*particle.v1/(u02)*d1_g11 -/*2*/particle.v1*particle.v2/(u0)*d1_g12 -/*2*/particle.v1*particle.v3/(u0)*d1_g13
		     -particle.v2*particle.v2/(u02)*d1_g22 -/*2*/particle.v2*particle.v3/(u0)*d1_g23 -particle.v3*particle.v3/(u02)*d1_g33;
	////////////
	auto du2 = -aa*u0*d2_alpha + particle.v1*d2B1 + particle.v2*d2B2 + particle.v3*d2B3 
		     -particle.v1*particle.v1/(u02)*d2_g11 -/*2*/particle.v1*particle.v2/(u0)*d2_g12 -/*2*/particle.v1*particle.v3/(u0)*d2_g13
		     -particle.v2*particle.v2/(u02)*d2_g22 -/*2*/particle.v2*particle.v3/(u0)*d2_g23 -particle.v3*particle.v3/(u02)*d2_g33;
	/////////////
	auto du3 = -aa*u0*d3_alpha + particle.v1*d3B1 + particle.v2*d3B2 + particle.v3*d3B3 
		     -particle.v1*particle.v1/(u02)*d3_g11 -/*2*/particle.v1*particle.v2/(u02)*d3_g12 -/*2*/particle.v1*particle.v3/(u0)*d3_g13
		     -particle.v2*particle.v2/(u02)*d3_g22 -/*2*/particle.v2*particle.v3/(u02)*d3_g23 -particle.v3*particle.v3/(u02)*d3_g33;

	Particle p;
	p.v1 = val(du1);
	p.v2 = val(du2);
	p.v3 = val(du3);

//////////////////////////////////////
	
	auto d1 = g_11*particle.v1/u0 + g_12*particle.v2/u0 + g_13*particle.v3/u0 - B_1;
	auto d2 = g_12*particle.v1/u0 + g_22*particle.v2/u0 + g_23*particle.v3/u0 - B_2;
	auto d3 = g_13*particle.v1/u0 + g_23*particle.v2/u0 + g_33*particle.v3/u0 - B_3;

//////////////////////////////////////

        // return derivatives
	p.x = val(d1);
	p.y = val(d2);
	p.z = val(d3);
	p.v1 = val(du1);
	p.v2 = val(du2);
	p.v3 = val(du3);

	return p;

}

