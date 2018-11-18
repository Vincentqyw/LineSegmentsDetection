#include <stdio.h>
#include <fstream>
#include "cv.h"
#include "highgui.h"
#include "CannyLine.h"

using namespace cv;
using namespace std;


void main()
{	
	string fileCur;
	cv::Mat img = imread( fileCur, 0 );

	CannyLine detector;
	std::vector<std::vector<float> > lines;
	detector.cannyLine( img, lines );

	// show
	cv::Mat imgShow( img.rows, img.cols, CV_8UC3, cv::Scalar( 255, 255, 255 ) );
	for ( int m=0; m<lines.size(); ++m )
	{
		cv::line( imgShow, cv::Point( lines[m][0], lines[m][1] ), cv::Point( lines[m][2], lines[m][3] ), cv::Scalar(0,0,0), 1, CV_AA );
	}
	imshow("",imgShow);
	cv::waitKey(0);
}