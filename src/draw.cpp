#include "draw.hpp"

// Draw trajectories of particles in the xy plan
void render_slabxy(cv::Mat & img, const Particles &p, double scale, int step, const std::string &s){
	
	int count = 0;
	for(int i=0; i<p._particles.size(); i+=step){
	if(fabs(p._particles[i].xi2) < 0.001){ 
		// convert to cartesian
		double x = scale * p._particles[i].x + double(img.cols/2);
		double y = scale * p._particles[i].y + double(img.rows/2);
		if(x>=0 && x<img.cols && y>=0 && y<img.rows){
			if(p._particles[i].intensity != 1.0)
			img.at<cv::Vec3b>(img.rows -1 -y,x) = cv::Vec3b(255*p._particles[i].intensity,0.,0.);
		}	
	}
	}
	cv::imshow(s, img);
	cv::waitKey(10);
	

}

// Draw trajectories of particles in the xz plan 
void render_slabxz(cv::Mat & img, const Particles &p, double scale, int step, const std::string &s){
	int count = 0;
	for(int i=0; i<p._particles.size(); i+=step){
		if(fabs(p._particles[i].xi) < 0.001){
		// convert to cartesian
		double x = scale * p._particles[i].x + double(img.cols/2);
		double z = scale * p._particles[i].z + double(img.rows/2);
		if(x>=0 && x<img.cols && z>=0 && z<img.rows){
			img.at<cv::Vec3b>(img.rows -1 -z,x) = cv::Vec3b(255*p._particles[i].intensity,0.,0);	
		}
		}
	}
	cv::imshow(s, img);
	cv::waitKey(10);
}

// Compute final image according to particle color
void render_img(cv::Mat & img, const Particles &p, double scale, double camera_distance){	

        // get min_max intensitys of luminosity
        float min = 100000000;
        float max = 0.;
        for (int i=0; i<p._particles.size(); i++)
        {
                if(p._particles[i].intensity > max) max = p._particles[i].intensity;
                if(p._particles[i].intensity < min) min = p._particles[i].intensity;
        }
        std::cout<<"min and max are "<<min<<" "<<max<<std::endl;
        std::cout<<p._particles.size()<<std::endl;
        int nb = 0;
	int h = img.rows;
	int w = img.cols;
	//std::cout<<h<<"x"<<w<<std::endl;
	for(int i=-h/2; i<h/2; i++){
		for(int j=-w/2; j<w/2; j++){
		        uchar sc = uchar(255*(p._particles[nb].intensity-min)/(max-min));
		        img.at<cv::Vec3b>(h/2-i-1,w/2+j) = cv::Vec3b(sc,sc,sc);
		        nb++;
		}
	}
	std::cout<<"final number of particle drawn "<<nb<<std::endl;
	cv::imshow("rendering", img);
	cv::waitKey(50);
}

