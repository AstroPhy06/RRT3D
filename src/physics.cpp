#include "physics.hpp"
#include <iostream>

// Relativistic focalization factor
double disk_model::focalization_factor(double th, double beta){
        return (cos(th)*(cos(th)+beta)+sin(th)*sin(th))/
        (pow(sin(th),2)*sqrt(1-beta*beta)+1/sqrt(1-beta*beta)*pow(beta+cos(th),2));
}

// For an optically thick disk, Luminosity as a function of r is given by
// L(r) ~ 1/r^2 * [1-sqrt(r_in/r)]  [Lacosta eq:60]
// the last part is added to speed-up the dimming of the outer part of the disk
double disk_model::L_radius(double r,double r_in, double r_out){
        //double range = r_out - r_in;
        double y = 1500/(r*r)*(1-sqrt(r_in/r)) - exp(-pow(r_out-r,2.0)/r_out);
        if (y < 0) y = 0.;
        return y;
}

float disk_model::emitted_light(const Particle &pos, const Particle &dir, double r_in, double r_out){

        // disk luminosity
        double fs = L_radius(pos.r, r_in, r_out);
        // velocity (beta)
        double beta = sqrt(1./(pos.r-3)); // assume M=1!
        // gamma
        double gamma = 1./sqrt(1-beta*beta);
        // theta
        double dot_prod = -dir.x*pos.y + dir.y*pos.x;
        double prod_norm = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z)*sqrt(pos.x*pos.x+pos.y*pos.y);
        double ctheta = dot_prod/prod_norm;
        double theta = acos(ctheta);
        
        double shift = 1./(gamma*(1.+beta*ctheta));
        //double factor = focalization_factor(theta,beta);
        return fs*pow(shift,3);
}

