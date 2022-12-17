#include "CannyLine.h"
#include "MetaLine.h"

CannyLine::CannyLine(void)
{
}

CannyLine::~CannyLine(void)
{
}

void CannyLine::cannyLine(cv::Mat &image,std::vector<std::vector<float> > &lines)
{
	MetaLine deterctor;
	float gausSigma=1.0;
	int gausHalfSize=1;
	deterctor.MetaLineDetection(image,gausSigma,gausHalfSize,lines);
}