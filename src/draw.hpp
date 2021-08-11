#ifndef DRAW_H_
#define DRAW_H_

#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <string>
#include "particle.hpp"

/**
* A few functions to render final image and monitor particle trajectories
*/

// draw particles trajectories in xz plan	
void render_slabxz(cv::Mat & img, const Particles &p, double scale, int step, const std::string &s);

// draw particles trajectories in xy plan
void render_slabxy(cv::Mat & img, const Particles &p, double scale, int step, const std::string &s);

// Render final image using photon color
void render_img(cv::Mat & img, const Particles &p, double scale, double camera_distance);
	
#endif
