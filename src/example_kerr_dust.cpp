#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <cmath>
#include <chrono>

#include "particle.hpp"
#include "draw.hpp"
#include "raytracing.hpp"
#include "physics.hpp"
#include "postprocess.hpp"

using namespace cv;
using namespace std::chrono;


double M = 1.0;
double a = 0.95;

//***************************************************************
// Example of a dust cloud of particles orbiting around a Kerr BH
//***************************************************************

int main(int argc, char** argv )
{
        /////////////////////////////////////////////////
        //    1 - Initialization
        /////////////////////////////////////////////////
        
	// Define window properties
	int slice_hw = 200;
	cv::Mat Ixy(slice_hw,slice_hw,CV_8UC3, cv::Scalar(255,255,255));  // rendering xy slice
	cv::Mat  Ixz(slice_hw,slice_hw,CV_8UC3, cv::Scalar(255,255,255)); // rendering xz slice
	
	int s = 10000;            // Simulation steps
	double dt = -0.5;     // time step in second
	double scale = slice_hw/100.; 
	int nb_part = 1000;      // total number of particles (one for each pixel)
	double camera_distance = 50; // camera distance from BH

	// Define spacetime properties
	params par(a,M);        // a: BH spin; M: BH mass
	KerrSchild KS(par);
        
         
        // Init photons
	Particles dust;
	dust.init_as_diff_disk(25, 15, 1., nb_part);
        	
	// Init raytracer
	Raytracer<KerrSchild, disk_model> rt(KS,dust);

        /////////////////////////////////////////////////
        //    2 - Launch computation
        /////////////////////////////////////////////////
        
        auto start = high_resolution_clock::now();

	for(int j=0; j<s; j++)
	{
	        std::cout<<"Step "<<j<<std::endl;
	        Ixy = cv::Mat(slice_hw,slice_hw,CV_8UC3, cv::Scalar(255,255,255));
	        Ixz = cv::Mat(slice_hw,slice_hw,CV_8UC3, cv::Scalar(255,255,255));
	        // compute all particles position increment
	        rt.auto_integrate_all_once_RK4(1,dt,nb_part); //put epsilon=1 for material particles
	        // Draw particles
	        render_slabxy(Ixy,rt.get_particles(),scale, 1, "dust_xy");
	        render_slabxz(Ixz,rt.get_particles(),scale, 1, "dust_xz");
	}
	
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	std::cout << "Time to advance simulation: "<< duration.count() << " seconds" << std::endl;
        	
	return 0;
}
