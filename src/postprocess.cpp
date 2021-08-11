#include "postprocess.hpp"


void PostProcess::post_process()
{
        cv::bilateralFilter( _img, _output, 5, 80, 40);
}

void PostProcess::post_process_with_scale_up()
{
        cv::bilateralFilter( _img, _output, 5, 80, 40);
        cv::pyrUp(_output,_output, cv::Size( _img.cols*2.0, _img.rows*2.0 ) );      // inplace modification!
}
