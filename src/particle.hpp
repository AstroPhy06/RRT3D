#ifndef PARTICLE_H_
#define PARTICLE_H_
#include <vector>
#include <random>

/**
* Define particle class
*/
class Particle
{
public:
	double x,y,z,t;                   // cartesian coordinates
	double r,theta,phi;             // spherical coordinates
	double v1,v2,v3;                // cartesian velocities
	double u_r,u_theta,u_phi;       //spherical velocities
	double xi,xi2;                  // initial angle
	int jump;                       // particle escaped or captured by BH
	float intensity;                    // Luminosity of the particle
	
	Particle()
	{}

	
	Particle(double x0, double y0, double z0, double v10, double v20, double v30) : x(x0), y(y0), z(z0), v1(v10), v2(v20), v3(v30)
	{}
	
	// Define a few maths operations: what it means to add/multiply particles?
	Particle &operator+=(Particle &p1)
	{
		x += p1.x;
		y += p1.y;
		z += p1.z;
		v1 += p1.v1;
		v2 += p1.v2;
		v3 += p1.v3;
		return *this;
	}

	
	Particle operator+(const Particle &p1) const
	{
		Particle result;
		result.x = this->x + p1.x;
		result.y = this->y + p1.y;
		result.z = this->z + p1.z;
		result.v1 = this->v1 + p1.v1;
		result.v2 = this->v2 + p1.v2;
		result.v3 = this->v3 + p1.v3;
		return result;
	}
	
	Particle &operator/(double a)
	{
		this->x /= a;
		this->y /= a;
		this->z /= a;
		this->v1 /= a;
		this->v2 /= a;
		this->v3 /= a;
		return *this;
	}
	
	Particle &operator*(double a)
	{
		this->x *= a;
		this->y *= a;
		this->z *= a;
		this->v1 *= a;
		this->v2 *= a;
		this->v3 *= a;
		return *this;
	}
	
	Particle &operator=(Particle p)
	{
	        x=p.x;y=p.y;z=p.z;
	        r=p.r;theta=p.theta,phi=p.phi;
	        v1=p.v1;v2=p.v2;v3=p.v3;
	        u_r=p.u_r;
	        return *this;
	}
};

/**
* Collection of particles
*/
class Particles
{	
protected:
	
public:
        
        std::vector<Particle> _particles;
        
        // for rendering, photons are organized in a window
	void init_as_window(int h, int w, double view_field, double view_angle, double camera_distance);

        // for material particles, create a disk of randomly distributed particles	
	void init_as_diff_disk(double radius, double dr, double thickness, int nb_part);
};

#endif
