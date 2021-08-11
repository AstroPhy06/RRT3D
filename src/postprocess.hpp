#ifndef POSTPROCESS_H_
#define POSTPROCESS_H_

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>

class PostProcess
{
        cv::Mat _img;
        cv::Mat _output;
        
public:

        PostProcess(const cv::Mat &img){ img.copyTo(_img);}

        void post_process();
        
        void post_process_with_scale_up();
        
        cv::Mat &get_img(){return _output;}
        
};

#endif
