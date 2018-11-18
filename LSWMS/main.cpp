/*
 * Project:  lineSegments (LSWMS Line Segment using Weighted Mean-Shift)
 *
 * File:     main.cpp
 *
 * Contents: Creation, initialisation and usage of LSWMS object
 *           for the detection of line segments in images or videos
 *
 * Author:   Marcos Nieto <marcos.nieto.doncel@gmail.com>
 *
 * Homepage: www.marcosnieto.net
 */

#ifdef WIN32
	#include <windows.h>
	#include <time.h>
#endif

#ifdef linux
	#include <stdio.h>
	#include <sys/time.h>
	#include <time.h>
#endif

#include <iostream>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "LSWMS.h"

using namespace std;
using namespace cv;

// Timing
#ifdef WIN32
	double t1, t2;	
#else
	int t1, t2;	
	struct timeval ts;
#endif
double t;

void help()
{
	 cout << "/*\n"
         << " **************************************************************************************************\n"
		 << " * Line segment detection using WMS \n"
         << " * ----------------------------------------------------\n"		 
		 << " * \n"
		 << " * Author:Marcos Nieto\n"
		 << " * www.marcosnieto.net\n"
		 << " * marcos.nieto.doncel@gmail.com\n"
		 << " * \n"
		 << " * Date:01/04/2012\n"
		 << " **************************************************************************************************\n"
		 << " * \n"
		 << " * Usage: \n"		 		 
		 << " *		-video		# Specifies video file as input (if not specified, camera is used) \n"
		 << " *		-image		# Specifies image file as input (if not specified, camera is used) \n"
		 << " *		-verbose	# Actives verbose: ON, OFF (default)\n"
		 << " *		-play		# ON: the video runs until the end; OFF: frame by frame (key press event)\n"
		 << " *		-resizedWidth	# Specifies the desired width of the image (the height is computed to keep aspect ratio)\n"
		 << " *		-numMaxLSegs	# Specifies the maximum number of line segments to detected.\n"		 
		 << " *\n"
		 << " * 	-hough		# ON: Applies OpenCV HoughLinesP for comparison (PPHT)\n"
		 << " *\n"
		 << " * Example:\n"
		 << " *		lineSegment -video myVideo.avi -verbose ON -play OFF\n"
		 << " *		lineSegment -image myImage.jpg -resizedWidth 300 -numMaxLSegs 100\n"
		 << " *		lineSegment -image myImage.jpg -hough ON\n"
		 << " * \n"
		 << " * Keys:\n"
		 << " *		Esc: Quit\n"
         << " */\n" << endl;
}
void processPPHT(cv::Mat &img, std::vector<LSEG> &lSegs)
{
	cv::Mat imgGRAY;
	cv::cvtColor(img, imgGRAY, CV_RGB2GRAY);
	cv::Mat dst;
	cv::Canny(imgGRAY, dst, 20, 80, 3);
	std::vector<cv::Vec4i> lines;
	cv::HoughLinesP(dst, lines, 1, CV_PI/180, 80, 30, 10);

	LSEG lSeg;
	cv::Point p1, p2;

	lSegs.clear();

	for(size_t i=0; i<lines.size(); i++)
	{
		lSeg.clear();
		p1.x = lines[i][0];
		p1.y = lines[i][1];

		p2.x = lines[i][2];
		p2.y = lines[i][3];

		lSeg.push_back(p1);
		lSeg.push_back(p2);

		lSegs.push_back(lSeg);			
	}	
}
void drawPPHT(cv::Mat &dst, std::vector<LSEG> &lSegs, cv::Scalar color)
{
	for(size_t i=0; i<lSegs.size(); i++)
	{		
		cv::line(dst, lSegs[i][0], lSegs[i][1], color, 2);	
	}
}
/** Main function*/
int main(int argc, char** argv)
{	
	// Images
	cv::Mat inputImg, imgGRAY;	
	cv::Mat outputImg, outputImgPPHT;
	int procWidth=0, procHeight=0;
	cv::Size procSize;

	// Other variables
	char *videoFileName = 0;
	char *imageFileName = 0;
	cv::VideoCapture video;
	bool useCamera = true;
	
	bool playMode = true;
	bool stillImage = false;
	bool verbose = false;
	int numMaxLSegs = 0;	
	bool usePPHT = false;

	// Line segments (LSWMS and PPHT)
	std::vector<LSEG> lSegs, lSegsPPHT;
	std::vector<double> errors;

	// Start showing help
	help();

	// Parse arguments
	if(argc < 2)
		return -1;	
	for(int i=1; i<argc; i++)
	{
		const char* s = argv[i];

		if(strcmp(s, "-video" ) == 0)
		{
			// Input video is a video file
			videoFileName = argv[++i];
			useCamera = false;
		}
		else if(strcmp(s,"-hough") == 0)
		{
			const char* ss = argv[++i];
			if(strcmp(ss, "ON") == 0 || strcmp(ss, "on") == 0
				|| strcmp(ss, "TRUE") == 0 || strcmp(ss, "true") == 0 
				|| strcmp(ss, "YES") == 0 || strcmp(ss, "yes") == 0 )
				usePPHT = true;	
		}
		else if(strcmp(s,"-image") == 0)
		{
			// Input is a image file
			imageFileName = argv[++i];
			stillImage = true;
			useCamera = false;
		}		
		else if(strcmp(s, "-numMaxLSegs") == 0)
		{
			numMaxLSegs = atoi(argv[++i]);	
		}		
		else if(strcmp(s, "-resizedWidth") == 0)
		{
			procWidth = atoi(argv[++i]);
		}
		else if(strcmp(s, "-verbose" ) == 0)
		{
			const char* ss = argv[++i];
			if(strcmp(ss, "ON") == 0 || strcmp(ss, "on") == 0 
				|| strcmp(ss, "TRUE") == 0 || strcmp(ss, "true") == 0 
				|| strcmp(ss, "YES") == 0 || strcmp(ss, "yes") == 0 )
				verbose = true;			
		}
		else if(strcmp(s, "-play" ) == 0)
		{
			const char* ss = argv[++i];
			if(strcmp(ss, "OFF") == 0 || strcmp(ss, "off") == 0 
				|| strcmp(ss, "FALSE") == 0 || strcmp(ss, "false") == 0 
				|| strcmp(ss, "NO") == 0 || strcmp(ss, "no") == 0 
				|| strcmp(ss, "STEP") == 0 || strcmp(ss, "step") == 0)
				playMode = false;			
		}		
	}

	// Open video input
	if( useCamera )
		video.open(0);
	else
	{
		if(!stillImage)
			video.open(videoFileName);
	}

	// Check video input
	int width = 0, height = 0, fps = 0, fourcc = 0;
	if(!stillImage)
	{
		if( !video.isOpened() )
		{
			printf("ERROR: can not open camera or video file\n");
			return -1;
		}
		else
		{
			// Show video information
			width = (int) video.get(CV_CAP_PROP_FRAME_WIDTH);
			height = (int) video.get(CV_CAP_PROP_FRAME_HEIGHT);
			fps = (int) video.get(CV_CAP_PROP_FPS);
			fourcc = (int) video.get(CV_CAP_PROP_FOURCC);

			if(!useCamera)
				printf("Input video: (%d x %d) at %d fps, fourcc = %d\n", width, height, fps, fourcc);
			else
				printf("Input camera: (%d x %d) at %d fps\n", width, height, fps);
		}
	}
	else
	{
		inputImg = cv::imread(imageFileName);
		if(inputImg.empty())
			return -1;

		width = inputImg.cols;
		height = inputImg.rows;

		printf("Input image: %s, Size (%d x %d)\n", imageFileName, width, height);

		playMode = false;
	}

	// Resize	
	if(procWidth != 0)
	{	
		procHeight = (int)(height*((double)procWidth/width));
		procSize = cv::Size(procWidth, procHeight);

		printf("Resize to: (%d x %d)\n", procWidth, procHeight);	
	}
	else
		procSize = cv::Size(width, height);

	if(numMaxLSegs != 0) printf("NumMaxLSegs=%d\n", numMaxLSegs);
	
	// ---------------------------
	// Create and init LSWMS
	int R = 3;
	LSWMS lswms(procSize, R, numMaxLSegs, verbose);
	if(numMaxLSegs==0)
		printf("LSWMS object created: R=%d\n\n", R);
	else
		printf("LSWMS object created: R=%d, numMaxLSegs=%d\n\n", R, numMaxLSegs);
	// ---------------------------
	
	// MAIN LOOP
	int frameNum=0;
	for( ;; )
	{
		if(!stillImage)
		{
			//if(verbose) printf("\n-------------------------\nFRAME #%6d\n", frameNum);
			frameNum++;

			// Get current image		
			video >> inputImg;
		}	
		else
		{
			printf("-------------------------\n");
		}	

		if( inputImg.empty() )
			break;
		
		// Resize to processing size
		cv::resize(inputImg, inputImg, procSize);		

		// Color Conversion
		if(inputImg.channels() == 3)
		{
			cv::cvtColor(inputImg, imgGRAY, CV_BGR2GRAY);	
			inputImg.copyTo(outputImg);
			if(usePPHT)
				inputImg.copyTo(outputImgPPHT);			
		}
		else
		{
			inputImg.copyTo(imgGRAY);
			cv::cvtColor(inputImg, outputImg, CV_GRAY2BGR);
			if(usePPHT)
				cv::cvtColor(inputImg, outputImgPPHT, CV_GRAY2BGR);			
		}

		// ++++++++++++++++++++++++++++++++++++++++
		// Process LSWMS
		#ifdef WIN32
			t1 = ::GetTickCount();
		#else
			gettimeofday(&ts,0);
			t1 = (ts.tv_sec * 1000 + (ts.tv_usec / 1000));
		#endif
		lswms.run(inputImg, lSegs, errors);				
		#ifdef WIN32
			t2 = ::GetTickCount();
		#else
			gettimeofday(&ts,0);
			t2 = (ts.tv_sec * 1000 + (ts.tv_usec / 1000));
		#endif	

		// process time = t2 - t1		
		t = (double)t2-(double)t1;

		cv::Scalar mean, stddev;
		cv::meanStdDev(errors, mean, stddev);
		if(!stillImage)
			printf("Fr.#%d - LSWMS: %d lines / %.0f ms , Ang.Error: (Mean, Std)=(%.2f, %.2f)(deg)\n", frameNum, lSegs.size(), t, mean.val[0]*180/CV_PI, stddev.val[0]*180/CV_PI);
		else
			printf("LSWMS: %d segments\nAngular Error: Mean = %.2f (deg), Std = %.2f (deg)\nProcess Time = %.0f (ms)\n", lSegs.size(), mean.val[0]*180/CV_PI, stddev.val[0]*180/CV_PI,  t);
		
		//lswms.drawLSegs(outputImg, lSegs,CV_RGB(255,0,0), 2);			// drawing all line segments the same
		lswms.drawLSegs(outputImg, lSegs, errors);				// drawing according to errors
		// ++++++++++++++++++++++++++++++++++++++++				

		// ++++++++++++++++++++++++++++++++++++++++				
		if(usePPHT)
		{
			// Process PPHT
			#ifdef WIN32
				t1 = ::GetTickCount();
			#else
				gettimeofday(&ts,0);
				t1 = (ts.tv_sec * 1000 + (ts.tv_usec / 1000));
			#endif
			processPPHT(inputImg, lSegsPPHT);				
			#ifdef WIN32
				t2 = ::GetTickCount();
			#else
				gettimeofday(&ts,0);
				t2 = (ts.tv_sec * 1000 + (ts.tv_usec / 1000));
			#endif

			// process time = t2 - t1		
			t = (double)t2-(double)t1;

			drawPPHT(outputImgPPHT, lSegsPPHT, CV_RGB(0,0,255));

			if(!stillImage)
				printf("Fr.#%d - PPHT: %d lines / %.0f ms\n", frameNum, lSegsPPHT.size(), t);
			else
				printf("\nPPHT: %d segments\nProcess Time = %.0f (ms)\n", lSegsPPHT.size(), t);
			// ++++++++++++++++++++++++++++++++++++++++	
		}

		// View
		imshow("LSWMS", outputImg);	
		if(usePPHT)
			imshow("PPHT", outputImgPPHT);

		if(stillImage)
		{
			cv::imwrite("lswms.bmp", outputImg);
			if(usePPHT)
				cv::imwrite("ppht.bmp", outputImgPPHT);
		}
		

		if(playMode)
			cv::waitKey(1);
		else
			cv::waitKey(0);

		char q = (char)waitKey(1);
	
		if( q == 27 )
		{
			printf("\nStopped by user request\n");
			break;
		}	

		if(stillImage)
			break;
	} // main while

	if(!stillImage)
		video.release();
	
	return 0;



}
