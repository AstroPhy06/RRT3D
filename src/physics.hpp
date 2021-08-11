#ifndef PHYSICS_H_
#define PHYSICS_H_

#include "particle.hpp"


/**
* A simple class where to implement physics of the disk
*/
class disk_model
{
public:
        disk_model()
        {}
        
        disk_model(double in,double out,bool th): inner_limit(in), outer_limit(out), thin(th)
        {}
        
        /*disk_model(const disk_model &d){
        inner_limit = d.inner_limit;
        outer_limit = d.outer_limit;
        thin = d.thin;}*/
       
        // accretion disk parameters
	double inner_limit;
	double outer_limit;	
        bool thin;
        
        // Relativistic focalization factor
        double focalization_factor(double th, double beta);

        // Luminosity as a function of radius
        double L_radius(double r,double r_in, double r_out);

        // Compute final luminosity as a function of radius
        float emitted_light(const Particle &pos, const Particle &dir, double r_in, double r_out);
        
};

#endif


