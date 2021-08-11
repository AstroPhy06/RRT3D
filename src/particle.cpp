#include "particle.hpp"
#include <math.h>
#include <iostream>

void Particles::init_as_window(int h, int w, double view_field, double view_angle, double camera_distance)
{
    _particles.resize(h*w);
    double pixel_angle = view_field/(w/2.);
        
    int nb = 0;

    for(int i=-h/2; i<h/2; i++)
    {
        for(int j=-w/2; j<w/2; j++)
        {
            //Init cartesian coord
            _particles[nb].x = camera_distance*sin(view_angle);
            _particles[nb].y = 0.;
            _particles[nb].z = camera_distance*cos(view_angle);
            _particles[nb].t = 0;
            
            // Init spherical coord
            _particles[nb].theta = view_angle;
            _particles[nb].phi = 0.;
            
            // Init direction vector
            _particles[nb].xi = pixel_angle*j;
            _particles[nb].xi2 = -pixel_angle*i;
                        
            // Init velocity
            // correct only if phi = 0 but we always can redefine coordinates to get phi = 0
            _particles[nb].v1 = sin(_particles[nb].theta-_particles[nb].xi2)*cos(_particles[nb].xi);
            _particles[nb].v2 = sin(_particles[nb].xi);
            _particles[nb].v3 = cos(_particles[nb].theta-_particles[nb].xi2)*cos(_particles[nb].xi);		

            _particles[nb].jump = 0;
            _particles[nb].intensity = 0.0;
            
            nb++;
	}
    }
    std::cout<<"Initialized "<<nb<<" particles"<<std::endl;
}



void Particles::init_as_diff_disk(double radius, double dr, double thickness, int nb_part)
{
    _particles.resize(nb_part);
    
    // Init random functions
    std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<> rd_radius(radius,dr);
    std::uniform_real_distribution<double> rd_phi(0.,2.*M_PI);
    std::normal_distribution<> rd_thickness(-thickness,thickness);
    double velocity;
    
    for(int i=0; i<nb_part; i++){
	_particles[i].r =  rd_radius(e2);
	_particles[i].phi =  rd_phi(e2);
	
	_particles[i].x =  _particles[i].r*cos(_particles[i].phi);
	_particles[i].y =  _particles[i].r*sin(_particles[i].phi);
	_particles[i].z =  rd_thickness(e2);
	
	_particles[i].xi2 = 0.;
	_particles[i].xi = 0.;
	
	velocity = -1.0/sqrt(_particles[i].r);
		
	_particles[i].v1 = -velocity*sin(_particles[i].phi);
	_particles[i].v2 = velocity*cos(_particles[i].phi);
	_particles[i].v3 = 0.0;

	_particles[i].jump = 0;
    }
    
}


