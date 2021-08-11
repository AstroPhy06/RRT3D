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


int main(int argc, char** argv )
{
        /////////////////////////////////////////////////
        //    1 - Initialization
        /////////////////////////////////////////////////
        
	// Define rendereing window properties
	int h = 400, w = 600;   // rendering window size
	cv::Mat I(h,w,CV_8UC3, cv::Scalar(255,255,255));        // main rendering window
	int slice_hw = 200;
	cv::Mat Ixy(slice_hw,slice_hw,CV_8UC3, cv::Scalar(255,255,255));  // rendering xy slice
	cv::Mat  Ixz(slice_hw,slice_hw,CV_8UC3, cv::Scalar(255,255,255)); // rendering xz slice
	
	int s = 400;            // Simulation steps
	double dt = -0.35;     // time step in second
	double scale = slice_hw/100.; 
	int nb_part = h*w;      // total number of particles (one for each pixel)
	double camera_distance = 50; // camera distance from BH

	// Define spacetime properties
	params par(a,M);        // a: BH spin; M: BH mass
	KerrSchild_Schwarzschild KS(par);
        
        // Define camera parameters
        double view_angle = M_PI/2.- M_PI/30.; // position above the plan of BH' disk
        double view_field = M_PI/6;             // field of view of the camera
         
        // Init photons
	Particles photons;
	photons.init_as_window(h, w, view_field, view_angle, camera_distance);

        // Init accretion disk model
        disk_model disk;
        disk.thin = 0;
        disk.inner_limit = 7.;
        disk.outer_limit = 25.;
        	
	// Init raytracer
	Raytracer<KerrSchild_Schwarzschild, disk_model> rt(KS,photons,disk);

        /////////////////////////////////////////////////
        //    2 - Launch computation
        /////////////////////////////////////////////////
        
        auto start = high_resolution_clock::now();

	// Choose between slow computation with trajectory rendering or faster computation
	if(0){
	        for(int j=0; j<s; j++)
	        {
		        std::cout<<"Step "<<j<<std::endl;
		        // compute all photons position increment
		        rt.auto_integrate_all_once_RK4(0,dt,camera_distance+1);
		        // Draw trajectories
		        render_slabxy(Ixy,rt.get_particles(),scale, 1, "photons_xy");
		        render_slabxz(Ixz,rt.get_particles(),scale, 1, "photons_xz");
		}
	}
	else{
	        // compute all photons trajectories (no rendering of trajectories)
	        rt.auto_integrate_photon_KS_Euler(0,dt,camera_distance+1,s);
	}
	
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	std::cout << "Time to compute image: "<< duration.count() << " seconds" << std::endl;
	/////////////////////////////////////////////////
	//    3 - Render and post-process image
	/////////////////////////////////////////////////

        // Finally draw the image
        cv::namedWindow("rendering");
        Particles p_rendu = rt.get_particles();
	render_img(I,p_rendu, scale, camera_distance);
	cv::waitKey(0);
	
	// Post process image to smooth results
	PostProcess postp(I);
	postp.post_process_with_scale_up(); // scaling up x2
	cv::Mat imgf = postp.get_img();

        cv::imshow("final",imgf);
        cv::waitKey(0);
        
        /////////////////////////////////////////////////
	//    4 - Save results
	/////////////////////////////////////////////////
        
        // save the resutls
	std::stringstream snum("");
	snum << "../output/SC_" << h<<"x"<< w << ".jpeg";
	std::cout<<" saving "<<snum.str()<<std::endl;
	cv::imwrite(snum.str(),I);
	
	snum.str("");
	snum << "../output/SC_smoothed_" << h<<"x"<< w  << ".jpeg";
	std::cout<<" saving "<<snum.str()<<std::endl;
	cv::imwrite(snum.str(),imgf);

       // }
        	
	return 0;
}
