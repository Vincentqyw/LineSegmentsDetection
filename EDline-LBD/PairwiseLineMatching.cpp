/*
 * PairwiseLineMatching.cpp
 *
 *  Created on: 2011-8-19
 *      Author: lz
 */

#include "PairwiseLineMatching.hh"
#include <arlsmat.h>
#include <arlssym.h>

using namespace BIAS;
using namespace std;

#define Inf       1e10 //Infinity
//the resolution scale of theta histogram, used when compute angle histogram of lines in the image
#define ResolutionScale  20  //10 degree
/*The following two thresholds are used to decide whether the estimated global rotation is acceptable.
 *Some image pairs don't have a stable global rotation angle, e.g. the image pair of the wide baseline
 *non planar scene. */
#define AcceptableAngleHistogramDifference 0.49
#define AcceptableLengthVectorDifference   0.4


/*The following four thresholds are used to decide whether a line in the left and a line in the right
 *are the possible matched line pair. If they are, then their differences should be smaller than these
 *thresholds.*/
#define LengthDifThreshold                 4
#define AngleDifferenceThreshold           0.7854//45degree
#define DescriptorDifThreshold             0.35//0.35, or o.5 are good values for LBD
/*The following four threshold are used to decide whether the two possible matched line pair agree with each other.
 *They reflect the similarity of pairwise geometric information.
 */
#define RelativeAngleDifferenceThreshold   0.7854//45degree
#define IntersectionRationDifThreshold     1
#define ProjectionRationDifThreshold       1


//this is used when get matching result from principal eigen vector
#define WeightOfMeanEigenVec               0.1


void PairwiseLineMatching::LineMatching(ScaleLines &linesInLeft,ScaleLines &linesInRight,
		std::vector<unsigned int> &matchResult)
{
	//compute the global rotation angle of image pair;
	globalRotationAngle_ = GlobalRotationOfImagePair_(linesInLeft,linesInRight);
	BuildAdjacencyMatrix_(linesInLeft,linesInRight);
	MatchingResultFromPrincipalEigenvector_(linesInLeft,linesInRight,matchResult);
}



double PairwiseLineMatching::GlobalRotationOfImagePair_(ScaleLines &linesInLeft, ScaleLines &linesInRight)
{
  double TwoPI = 2*M_PI;
	double rotationAngle = TwoPI;

  //step 1: compute the angle histogram of lines in the left and right images
  unsigned int dim = 360/ResolutionScale; //number of the bins of histogram
  unsigned int index;//index in the histogram
  double direction;
  double scalar = 180/(ResolutionScale*3.1415927);//used when compute the index
  double angleShift = (ResolutionScale*M_PI)/360;//make sure zero is the middle of the interval

	Vector<double> angleHistLeft(dim);
	Vector<double> angleHistRight(dim);
	Vector<double> lengthLeft(dim);//lengthLeft[i] store the total line length of all the lines in the ith angle bin.
	Vector<double> lengthRight(dim);
  angleHistLeft.SetZero();
  angleHistRight.SetZero();
  lengthLeft.SetZero();
  lengthRight.SetZero();

  for(unsigned int linenum=0; linenum<linesInLeft.size(); linenum++){
  	direction = linesInLeft[linenum][0].direction+M_PI+angleShift;
  	direction = direction<TwoPI?direction:(direction-TwoPI);
  	index = floor(direction*scalar);
  	angleHistLeft[index] ++;
  	lengthLeft[index] += linesInLeft[linenum][0].lineLength;
  }
  for(unsigned int linenum=0; linenum<linesInRight.size(); linenum++){
  	direction = linesInRight[linenum][0].direction+M_PI+angleShift;
  	direction = direction<TwoPI?direction:(direction-TwoPI);
  	index = floor(direction*scalar);
  	angleHistRight[index] ++;
  	lengthRight[index] += linesInRight[linenum][0].lineLength;
  }

  angleHistLeft  = (1/angleHistLeft.NormL2())*angleHistLeft;
  angleHistRight = (1/angleHistRight.NormL2())*angleHistRight;
  lengthLeft     = (1/lengthLeft.NormL2())*lengthLeft;
  lengthRight    = (1/lengthRight.NormL2())*lengthRight;

//  angleHistLeft.Save("histLeft.txt");
//  angleHistRight.Save("histRight.txt");

  //step 2: find shift to decide the approximate global rotation
  Vector<double> difVec(dim);//the difference vector between left histogram and shifted right histogram
  double minDif=10;//the minimal angle histogram difference
  double secondMinDif = 10;//the second minimal histogram difference
  unsigned int minShift;//the shift of right angle histogram when minimal difference achieved
  unsigned int secondMinShift;//the shift of right angle histogram when second minimal difference achieved

  Vector<double> lengthDifVec(dim);//the length difference vector between left and right
  double minLenDif = 10;//the minimal length difference
  double secondMinLenDif = 10;//the second minimal length difference
  unsigned int minLenShift;//the shift of right length vector when minimal length difference achieved
  unsigned int secondMinLenShift;//the shift of right length vector when the second minimal length difference achieved

  double normOfVec;
  for(unsigned int shift=0; shift<dim; shift++){
  	for(unsigned int j=0; j<dim; j++){
  		index = j+shift;
  		index = index<dim?index:(index-dim);
  		difVec[j] = angleHistLeft[j] - angleHistRight[index];
  		lengthDifVec[j] = lengthLeft[j] - lengthRight[index];
  	}
    //find the minShift and secondMinShift for angle histogram
    normOfVec = difVec.NormL2();
    if(normOfVec<secondMinDif){
    	if(normOfVec<minDif){
    		secondMinDif   = minDif;
    		secondMinShift = minShift;
    		minDif   = normOfVec;
    		minShift = shift;
    	}else{
    		secondMinDif    = normOfVec;
    		secondMinShift  = shift;
    	}
    }
    //find the minLenShift and secondMinLenShift of length vector
    normOfVec = lengthDifVec.NormL2();
    if(normOfVec<secondMinLenDif){
    	if(normOfVec<minLenDif){
    		secondMinLenDif    = minLenDif;
    		secondMinLenShift  = minLenShift;
    		minLenDif   = normOfVec;
    		minLenShift = shift;
    	}else{
    		secondMinLenDif    = normOfVec;
    		secondMinLenShift  = shift;
    	}
    }
  }

  //first check whether there exist an approximate global rotation angle between image pair
  if(minDif<AcceptableAngleHistogramDifference && minLenDif<AcceptableLengthVectorDifference){
  	rotationAngle = minShift*ResolutionScale;
  	if(rotationAngle>90 && 360-rotationAngle>90){
  		//In most case we believe the rotation angle between two image pairs should belong to [-Pi/2, Pi/2]
  		rotationAngle = rotationAngle - 180;
  	}
  	rotationAngle = rotationAngle*M_PI/180;
  }

  cout<<"minimal histgram distance = "<<minDif<<", Approximate global rotation angle = "<<rotationAngle<<endl;
	return rotationAngle;
}


void PairwiseLineMatching::BuildAdjacencyMatrix_(ScaleLines &linesInLeft,ScaleLines &linesInRight)
{
	double TwoPI = 2*M_PI;
	unsigned int numLineLeft  = linesInLeft.size();
	unsigned int numLineRight = linesInRight.size();
	/*first step, find nodes which are possible correspondent lines in the left and right images according to
	 *their direction, gray value  and gradient magnitude.
	 */
	nodesList_.clear();
	double angleDif;
	double lengthDif;

	unsigned int dimOfDes = linesInLeft[0][0].descriptor.size();
	Matrix<double> desDisMat(numLineLeft,numLineRight);//store the descriptor distance of lines in left and right images.
	std::vector<float> desLeft;
	std::vector<float> desRight;

// //store descriptor for debug
//	Matrix<double> desCripLeft(numLineLeft,dimOfDes);
//	Matrix<double> desCripRight(numLineRight,dimOfDes);
//	for(unsigned int i=0; i<numLineLeft; i++){
//		for(unsigned int j=0; j<dimOfDes; j++){
//			desCripLeft[i][j] = linesInLeft[i].decriptor[j];
//		}
//	}
//	for(unsigned int i=0; i<numLineRight; i++){
//		for(unsigned int j=0; j<dimOfDes; j++){
//			desCripRight[i][j] = linesInRight[i].decriptor[j];
//		}
//	}
//	desCripLeft.Save("DescriptorLeft.txt");
//	desCripRight.Save("DescriptorRight.txt");

	//first compute descriptor distances

	float *desL, *desR, *desMax, *desOld;

	float minDis,dis,temp;
	for(int idL=0; idL<numLineLeft; idL++){
		short sameLineSize = linesInLeft[idL].size();
		for(int idR=0; idR<numLineRight; idR++){
			minDis = 100;
			short sameLineSizeR = linesInRight[idR].size();
			for(short lineIDInSameLines = 0; lineIDInSameLines<sameLineSize; lineIDInSameLines++){
				desOld = linesInLeft[idL][lineIDInSameLines].descriptor.data();
				for(short lineIDInSameLinesR = 0; lineIDInSameLinesR<sameLineSizeR; lineIDInSameLinesR++){
					desL = desOld;
					desR = linesInRight[idR][lineIDInSameLinesR].descriptor.data();
					desMax = desR+dimOfDes;
					dis = 0;
					while(desR<desMax){
						temp = *(desL++) - *(desR++);
						dis += temp*temp;
					}
					dis = sqrt(dis);
					if(dis<minDis){
						minDis = dis;
					}
				}
			}//end for(short lineIDInSameLines = 0; lineIDInSameLines<sameLineSize; lineIDInSameLines++)
			desDisMat[idL][idR] = minDis;
		}//end for(int idR=0; idR<rightSize; idR++)
	}// end for(int idL=0; idL<leftSize; idL++)


	for(unsigned int i=0; i<numLineLeft; i++){
		for(unsigned int j=0; j<numLineRight; j++){
			if(desDisMat[i][j]>DescriptorDifThreshold){
				continue;//the descriptor difference is too large;
			}

			if(globalRotationAngle_<TwoPI){//there exist a global rotation angle between two image
				lengthDif = fabs(linesInLeft[i][0].lineLength - linesInRight[j][0].lineLength)/MIN(linesInLeft[i][0].lineLength,linesInRight[j][0].lineLength);
				if(lengthDif>LengthDifThreshold){
					continue;//the length difference is too large;
				}
				angleDif = fabs(linesInLeft[i][0].direction+globalRotationAngle_-linesInRight[j][0].direction);
				if(fabs(TwoPI-angleDif)>AngleDifferenceThreshold&&angleDif>AngleDifferenceThreshold){
					continue;//the angle difference is too large;
				}

				Node node;//line i in left image and line j in right image pass the test, (i,j) is a possible matched line pair.
				node.leftLineID  = i;
				node.rightLineID = j;
				nodesList_.push_back(node);
			}else{//there doesn't exist a global rotation angle between two image, so the angle difference test is canceled.
				lengthDif = fabs(linesInLeft[i][0].lineLength - linesInRight[j][0].lineLength)/MIN(linesInLeft[i][0].lineLength,linesInRight[j][0].lineLength);
				if(lengthDif>LengthDifThreshold){
					continue;//the length difference is too large;
				}
				Node node;//line i in left image and line j in right image pass the test, (i,j) is a possible matched line pair.
				node.leftLineID  = i;
				node.rightLineID = j;
				nodesList_.push_back(node);
			}
		}//end inner loop
	}
	cout<<"the number of possible matched line pair = "<<nodesList_.size()<<endl;
//	desDisMat.Save("DescriptorDis.txt");

	/*Second step, build the adjacency matrix which reflect the geometric constraints between nodes.
	 *The matrix is stored in the Compressed Sparse Column(CSC) format.
	 */
	unsigned int dim = nodesList_.size();// Dimension of the problem.
	int nnz = 0;// Number of nonzero elements in adjacenceMat.
	/*adjacenceVec only store the lower part of the adjacency matrix which is a symmetric matrix.
	 *                    | 0  1  0  2  0 |
	 *                    | 1  0  3  0  1 |
	 *eg:  adjMatrix =    | 0  3  0  2  0 |
	 *                    | 2  0  2  0  3 |
	 *                    | 0  1  0  3  0 |
	 *     adjacenceVec = [0,1,0,2,0,0,3,0,1,0,2,0,0,3,0]
	 */
	//	Matrix<double> testMat(dim,dim);
	//	testMat.SetZero();
	Vector<double> adjacenceVec(dim*(dim+1)/2);
	adjacenceVec.SetZero();
	/*In order to save computational time, the following variables are used to store
	 *the pairwise geometric information which has been computed and will be reused many times
	 *latter. The reduction of computational time is at the expenses of memory consumption.
	 */
	Matrix<unsigned int> bComputedLeft(numLineLeft,numLineLeft);//flag to show whether the ith pair of left image has already been computed.
	bComputedLeft.SetZero();
	Matrix<double> intersecRatioLeft(numLineLeft,numLineLeft);//the ratio of intersection point and the line in the left pair
	Matrix<double> projRatioLeft(numLineLeft,numLineLeft);//the point to line distance divided by the projected length of line in the left pair.

	Matrix<unsigned int> bComputedRight(numLineRight,numLineRight);//flag to show whether the ith pair of right image has already been computed.
	bComputedRight.SetZero();
	Matrix<double>       intersecRatioRight(numLineRight,numLineRight);//the ratio of intersection point and the line in the right pair
	Matrix<double>       projRatioRight(numLineRight,numLineRight);//the point to line distance divided by the projected length of line in the right pair.



	unsigned int idLeft1, idLeft2;//the id of lines in the left pair
	unsigned int idRight1, idRight2;//the id of lines in the right pair
	double relativeAngleLeft, relativeAngleRight;//the relative angle of each line pair
	double gradientMagRatioLeft, gradientMagRatioRight;//the ratio of gradient magnitude of lines in each pair

	double iRatio1L,iRatio1R,iRatio2L,iRatio2R;
	double pRatio1L,pRatio1R,pRatio2L,pRatio2R;

	double relativeAngleDif,  gradientMagRatioDif, iRatioDif, pRatioDif;

	double interSectionPointX,interSectionPointY;
	double a1,a2,b1,b2,c1,c2;//line1: a1 x + b1 y + c1 =0; line2: a2 x + b2 y + c2=0
	double a1b2_a2b1;//a1b2-a2b1
	double length1, length2,len;
	double disX,disY;
	double disS,disE;
	double similarity;
	for(unsigned int j=0; j<dim; j++){//column
		idLeft1  = nodesList_[j].leftLineID;
		idRight1 = nodesList_[j].rightLineID;
		for(unsigned int i=j+1; i<dim; i++){//row
			idLeft2  = nodesList_[i].leftLineID;
			idRight2 = nodesList_[i].rightLineID;
			if((idLeft1==idLeft2)||(idRight1==idRight2)){
				continue;//not satisfy the one to one match condition
			}
			//first compute the relative angle between left pair and right pair.
			relativeAngleLeft  = linesInLeft[idLeft1][0].direction - linesInLeft[idLeft2][0].direction;
			relativeAngleLeft  = (relativeAngleLeft<M_PI)?relativeAngleLeft:(relativeAngleLeft-TwoPI);
			relativeAngleLeft  = (relativeAngleLeft>(-M_PI))?relativeAngleLeft:(relativeAngleLeft+TwoPI);
			relativeAngleRight = linesInRight[idRight1][0].direction - linesInRight[idRight2][0].direction;
			relativeAngleRight = (relativeAngleRight<M_PI)?relativeAngleRight:(relativeAngleRight-TwoPI);
			relativeAngleRight = (relativeAngleRight>(-M_PI))?relativeAngleRight:(relativeAngleRight+TwoPI);
			relativeAngleDif   = fabs(relativeAngleLeft - relativeAngleRight);
			if((TwoPI-relativeAngleDif)>RelativeAngleDifferenceThreshold&&relativeAngleDif>RelativeAngleDifferenceThreshold){
				continue;//the relative angle difference is too large;
			}else if((TwoPI-relativeAngleDif)<RelativeAngleDifferenceThreshold){
				relativeAngleDif = TwoPI-relativeAngleDif;
			}

			//at last, check the intersect point ratio and point to line distance ratio
			//check whether the geometric information of pairs (idLeft1,idLeft2) and (idRight1,idRight2) have already been computed.
			if(!bComputedLeft[idLeft1][idLeft2]){//have not been computed yet
				/*compute the intersection point of segment i and j.
				 *a1x + b1y + c1 = 0 and a2x + b2y + c2 = 0.
				 *x = (c2b1 - c1b2)/(a1b2 - a2b1) and
				 *y = (c1a2 - c2a1)/(a1b2 - a2b1)*/
				a1 = linesInLeft[idLeft1][0].endPointY - linesInLeft[idLeft1][0].startPointY;//disY
				b1 = linesInLeft[idLeft1][0].startPointX - linesInLeft[idLeft1][0].endPointX;//-disX
				c1 = (0 - b1*linesInLeft[idLeft1][0].startPointY) - a1 * linesInLeft[idLeft1][0].startPointX;//disX*sy - disY*sx
				length1 = linesInLeft[idLeft1][0].lineLength;

				a2 = linesInLeft[idLeft2][0].endPointY - linesInLeft[idLeft2][0].startPointY;//disY
				b2 = linesInLeft[idLeft2][0].startPointX - linesInLeft[idLeft2][0].endPointX;//-disX
				c2 = (0 - b2*linesInLeft[idLeft2][0].startPointY) - a2 * linesInLeft[idLeft2][0].startPointX;//disX*sy - disY*sx
				length2 = linesInLeft[idLeft2][0].lineLength;

				a1b2_a2b1 = a1 * b2 - a2 * b1;
				if(fabs(a1b2_a2b1)<0.001){//two lines are almost parallel
					iRatio1L = Inf;
					iRatio2L = Inf;
				}else{
					interSectionPointX = (c2 * b1 - c1 * b2)/a1b2_a2b1;
					interSectionPointY = (c1 * a2 - c2 * a1)/a1b2_a2b1;
					//r1 = (s1I*s1e1)/(|s1e1|*|s1e1|)
					disX = interSectionPointX - linesInLeft[idLeft1][0].startPointX;
					disY = interSectionPointY - linesInLeft[idLeft1][0].startPointY;
					len  = disY*a1 - disX*b1;
					iRatio1L = len/(length1*length1);
					//r2 = (s2I*s2e2)/(|s2e2|*|s2e2|)
					disX = interSectionPointX - linesInLeft[idLeft2][0].startPointX;
					disY = interSectionPointY - linesInLeft[idLeft2][0].startPointY;
					len  = disY*a2 - disX*b2;
					iRatio2L = len/(length2*length2);
				}
				intersecRatioLeft[idLeft1][idLeft2] = iRatio1L;
				intersecRatioLeft[idLeft2][idLeft1] = iRatio2L;//line order changed

				/*project the end points of line1 onto line2 and compute their distances to line2;
				 */
				disS = fabs(a2*linesInLeft[idLeft1][0].startPointX + b2*linesInLeft[idLeft1][0].startPointY + c2)/length2;
				disE = fabs(a2*linesInLeft[idLeft1][0].endPointX + b2*linesInLeft[idLeft1][0].endPointY + c2)/length2;
				pRatio1L = (disS+disE)/length1;
				projRatioLeft[idLeft1][idLeft2] = pRatio1L;

				/*project the end points of line2 onto line1 and compute their distances to line1;
				 */
				disS = fabs(a1*linesInLeft[idLeft2][0].startPointX + b1*linesInLeft[idLeft2][0].startPointY + c1)/length1;
				disE = fabs(a1*linesInLeft[idLeft2][0].endPointX + b1*linesInLeft[idLeft2][0].endPointY + c1)/length1;
				pRatio2L = (disS+disE)/length2;
				projRatioLeft[idLeft2][idLeft1] = pRatio2L;

				//mark them as computed
				bComputedLeft[idLeft1][idLeft2] = true;
				bComputedLeft[idLeft2][idLeft1] = true;
			}else{//read these information from matrix;
				iRatio1L = intersecRatioLeft[idLeft1][idLeft2];
				iRatio2L = intersecRatioLeft[idLeft2][idLeft1];
				pRatio1L = projRatioLeft[idLeft1][idLeft2];
				pRatio2L = projRatioLeft[idLeft2][idLeft1];
			}
			if(!bComputedRight[idRight1][idRight2]){//have not been computed yet
				a1 = linesInRight[idRight1][0].endPointY - linesInRight[idRight1][0].startPointY;//disY
				b1 = linesInRight[idRight1][0].startPointX - linesInRight[idRight1][0].endPointX;//-disX
				c1 = (0 - b1*linesInRight[idRight1][0].startPointY) - a1 * linesInRight[idRight1][0].startPointX;//disX*sy - disY*sx
				length1 = linesInRight[idRight1][0].lineLength;

				a2 = linesInRight[idRight2][0].endPointY - linesInRight[idRight2][0].startPointY;//disY
				b2 = linesInRight[idRight2][0].startPointX - linesInRight[idRight2][0].endPointX;//-disX
				c2 = (0 - b2*linesInRight[idRight2][0].startPointY) - a2 * linesInRight[idRight2][0].startPointX;//disX*sy - disY*sx
				length2 = linesInRight[idRight2][0].lineLength;

				a1b2_a2b1 = a1 * b2 - a2 * b1;
				if(fabs(a1b2_a2b1)<0.001){//two lines are almost parallel
					iRatio1R = Inf;
					iRatio2R = Inf;
				}else{
					interSectionPointX = (c2 * b1 - c1 * b2)/a1b2_a2b1;
					interSectionPointY = (c1 * a2 - c2 * a1)/a1b2_a2b1;
					//r1 = (s1I*s1e1)/(|s1e1|*|s1e1|)
					disX = interSectionPointX - linesInRight[idRight1][0].startPointX;
					disY = interSectionPointY - linesInRight[idRight1][0].startPointY;
					len  = disY*a1 - disX*b1;//because b1=-disX
					iRatio1R = len/(length1*length1);
					//r2 = (s2I*s2e2)/(|s2e2|*|s2e2|)
					disX = interSectionPointX - linesInRight[idRight2][0].startPointX;
					disY = interSectionPointY - linesInRight[idRight2][0].startPointY;
					len  = disY*a2 - disX*b2;//because b2=-disX
					iRatio2R = len/(length2*length2);
				}
				intersecRatioRight[idRight1][idRight2] = iRatio1R;
				intersecRatioRight[idRight2][idRight1] = iRatio2R;//line order changed
				/*project the end points of line1 onto line2 and compute their distances to line2;
				 */
				disS = fabs(a2*linesInRight[idRight1][0].startPointX + b2*linesInRight[idRight1][0].startPointY + c2)/length2;
				disE = fabs(a2*linesInRight[idRight1][0].endPointX + b2*linesInRight[idRight1][0].endPointY + c2)/length2;
				pRatio1R = (disS+disE)/length1;
				projRatioRight[idRight1][idRight2] = pRatio1R;

				/*project the end points of line2 onto line1 and compute their distances to line1;
				 */
				disS = fabs(a1*linesInRight[idRight2][0].startPointX + b1*linesInRight[idRight2][0].startPointY + c1)/length1;
				disE = fabs(a1*linesInRight[idRight2][0].endPointX + b1*linesInRight[idRight2][0].endPointY + c1)/length1;
				pRatio2R = (disS+disE)/length2;
				projRatioRight[idRight2][idRight1] = pRatio2R;

				//mark them as computed
				bComputedRight[idRight1][idRight2] = true;
				bComputedRight[idRight2][idRight1] = true;
			}else{//read these information from matrix;
				iRatio1R = intersecRatioRight[idRight1][idRight2];
				iRatio2R = intersecRatioRight[idRight2][idRight1];
				pRatio1R = projRatioRight[idRight1][idRight2];
				pRatio2R = projRatioRight[idRight2][idRight1];
			}
			pRatioDif = MIN(fabs(pRatio1L-pRatio1R), fabs(pRatio2L-pRatio2R));
			if(pRatioDif>ProjectionRationDifThreshold){
				continue;//the projection length ratio difference is too large;
			}
			if((iRatio1L==Inf)||(iRatio2L==Inf)||(iRatio1R==Inf)||(iRatio2R==Inf)){
				//don't consider the intersection length ratio
				similarity = 4 - desDisMat[idLeft1][idRight1]/DescriptorDifThreshold
				- desDisMat[idLeft2][idRight2]/DescriptorDifThreshold
				- pRatioDif/ProjectionRationDifThreshold
				- relativeAngleDif/RelativeAngleDifferenceThreshold;
				adjacenceVec[(2*dim-j-1)*j/2+i] = similarity;
				nnz++;
				//				testMat[i][j] = similarity;
				//				testMat[j][i] = similarity;
			}else{
				iRatioDif = MIN(fabs(iRatio1L-iRatio1R), fabs(iRatio2L-iRatio2R));
				if(iRatioDif>IntersectionRationDifThreshold){
					continue;//the intersection length ratio difference is too large;
				}
				//now compute the similarity score between two line pairs.
				similarity = 5 - desDisMat[idLeft1][idRight1]/DescriptorDifThreshold
				- desDisMat[idLeft2][idRight2]/DescriptorDifThreshold
				- iRatioDif/IntersectionRationDifThreshold - pRatioDif/ProjectionRationDifThreshold
				- relativeAngleDif/RelativeAngleDifferenceThreshold;
				adjacenceVec[(2*dim-j-1)*j/2+i] = similarity;
				nnz++;
				//				testMat[i][j] = similarity;
				//				testMat[j][i] = similarity;
			}
		}
	}

	// pointer to an array that stores the nonzero elements of Adjacency matrix.
	double* adjacenceMat = new double[nnz];
	// pointer to an array that stores the row indices of the non-zeros in adjacenceMat.
	int*    irow = new int[nnz];
	// pointer to an array of pointers to the beginning of each column of adjacenceMat.
	int*    pcol = new int[dim+1];
	int idOfNNZ = 0;//the order of none zero element
	pcol[0] = 0;
	unsigned int tempValue;
	for(unsigned int j=0; j<dim; j++){//column
		for(unsigned int i=j; i<dim; i++){//row
			tempValue = (2*dim-j-1)*j/2+i;
			if(adjacenceVec[tempValue]!=0){
				adjacenceMat[idOfNNZ] = adjacenceVec[tempValue];
				irow[idOfNNZ] = i;
				idOfNNZ++;
			}
		}
		pcol[j+1] = idOfNNZ;
	}

	//	testMat.Save("testmat.txt");

	//	cout<<"CCS Mat"<<endl<<"adjacenceMat= ";
	//	for(int i=0; i<nnz; i++){
	//		cout<<adjacenceMat[i]<<", ";
	//	}
	//	cout<<endl<<"irow = ";
	//	for(int i=0; i<nnz; i++){
	//		cout<<irow[i]<<", ";
	//	}
	//	cout<<endl<<"pcol = ";
	//	for(int i=0; i<dim+1; i++){
	//		cout<<pcol[i]<<", ";
	//	}
	//	cout<<endl;

	/*Third step, solve the principal eigenvector of the adjacency matrix using Arpack lib.
	 */
	ARluSymMatrix<double> arMatrix(dim, nnz, adjacenceMat, irow, pcol);
	ARluSymStdEig<double> dprob(2, arMatrix, "LM");// Defining what we need: the first eigenvector of arMatrix with largest magnitude.
	// Finding eigenvalues and eigenvectors.
	dprob.FindEigenvectors();
	cout<<"Number of 'converged' eigenvalues  : " << dprob.ConvergedEigenvalues() <<endl;
	//  cout<<"eigenvalue is = "<<dprob.Eigenvalue(0)<<", and "<<dprob.Eigenvalue(1)<<endl;
	//  if(dprob.EigenvectorsFound()){
	//  	for(unsigned int j=0; j<dim; j++){
	//  		cout<< dprob.Eigenvector(1,j) <<", ";
	//  	}
	//  	cout<<endl;
	//  }
	eigenMap_.clear();

  double meanEigenVec = 0;
  if(dprob.EigenvectorsFound()){
  	double value;
  	for(unsigned int j=0; j<dim; j++){
  		value = fabs(dprob.Eigenvector(1,j));
  		meanEigenVec += value;
  		eigenMap_.insert(std::make_pair(value,j));
  	}
  }
  minOfEigenVec_ = WeightOfMeanEigenVec*meanEigenVec/dim;
	delete[] adjacenceMat;
	delete[] irow;
	delete[] pcol;
}

void PairwiseLineMatching::MatchingResultFromPrincipalEigenvector_(ScaleLines &linesInLeft,ScaleLines &linesInRight,
		std::vector<unsigned int > &matchResult)
{
	double TwoPI = 2*M_PI;
	std::vector<unsigned int > matchRet1;
	std::vector<unsigned int > matchRet2;
	double matchScore1 = 0;
	double matchScore2 = 0;
	EigenMAP mapCopy = eigenMap_;
	unsigned int dim = nodesList_.size();
	EigenMAP::iterator iter;
	unsigned int id,idLeft2,idRight2;
	double sideValueL, sideValueR;
	double pointX,pointY;
	double relativeAngleLeft, relativeAngleRight;//the relative angle of each line pair
  double relativeAngleDif;


  //store eigenMap for debug
  std::fstream resMap;
  ostringstream fileNameMap;
  fileNameMap<<"eigenVec.txt";
  resMap.open(fileNameMap.str().c_str(), std::ios::out);

  Matrix<double> mat(linesInLeft.size(),linesInRight.size());
  mat.SetZero();
  for(iter = eigenMap_.begin();iter!=eigenMap_.end(); iter++){
  	id = iter->second;
  	resMap<<nodesList_[id].leftLineID<<"    "<<nodesList_[id].rightLineID<<"   "<<iter->first<<endl;
  	mat[nodesList_[id].leftLineID][nodesList_[id].rightLineID] = iter->first;
  }
  mat.Save("eigenMap.txt");
  resMap.flush();
  resMap.close();


	/*first try, start from the top element in eigenmap */
	while(1){
		iter = eigenMap_.begin();
		//if the top element in the map has small value, then there is no need to continue find more matching line pairs;
		if(iter->first < minOfEigenVec_){
			break;
		}
		id = iter->second;
		unsigned int idLeft1 = nodesList_[id].leftLineID;
		unsigned int idRight1= nodesList_[id].rightLineID;
		matchRet1.push_back(idLeft1);
		matchRet1.push_back(idRight1);
		matchScore1 += iter->first;
		eigenMap_.erase(iter++);
		//remove all potential assignments in conflict with top matched line pair
		double xe_xsLeft = linesInLeft[idLeft1][0].endPointX-linesInLeft[idLeft1][0].startPointX;
		double ye_ysLeft = linesInLeft[idLeft1][0].endPointY-linesInLeft[idLeft1][0].startPointY;
		double xe_xsRight = linesInRight[idRight1][0].endPointX-linesInRight[idRight1][0].startPointX;
		double ye_ysRight = linesInRight[idRight1][0].endPointY-linesInRight[idRight1][0].startPointY;
		double coefLeft  = sqrt(xe_xsLeft*xe_xsLeft+ye_ysLeft*ye_ysLeft);
		double coefRight = sqrt(xe_xsRight*xe_xsRight+ye_ysRight*ye_ysRight);
		for( ; iter->first >= minOfEigenVec_; ){
			id = iter->second;
			idLeft2 = nodesList_[id].leftLineID;
			idRight2= nodesList_[id].rightLineID;
			//check one to one match condition
			if((idLeft1==idLeft2)||(idRight1==idRight2)){
				eigenMap_.erase(iter++);
				continue;//not satisfy the one to one match condition
			}
			//check sidedness constraint, the middle point of line2 should lie on the same side of line1.
			//sideValue = (y-ys)*(xe-xs)-(x-xs)*(ye-ys);
			pointX = 0.5*(linesInLeft[idLeft2][0].startPointX+linesInLeft[idLeft2][0].endPointX);
			pointY = 0.5*(linesInLeft[idLeft2][0].startPointY+linesInLeft[idLeft2][0].endPointY);
			sideValueL = (pointY-linesInLeft[idLeft1][0].startPointY)*xe_xsLeft
			 - (pointX-linesInLeft[idLeft1][0].startPointX)*ye_ysLeft;
			sideValueL = sideValueL/coefLeft;
			pointX = 0.5*(linesInRight[idRight2][0].startPointX+linesInRight[idRight2][0].endPointX);
			pointY = 0.5*(linesInRight[idRight2][0].startPointY+linesInRight[idRight2][0].endPointY);
			sideValueR = (pointY-linesInRight[idRight1][0].startPointY)*xe_xsRight
			 - (pointX-linesInRight[idRight1][0].startPointX)*ye_ysRight;
			sideValueR = sideValueR/coefRight;
			if(sideValueL*sideValueR<0&&fabs(sideValueL)>5&&fabs(sideValueR)>5){//have the different sign, conflict happens.
				eigenMap_.erase(iter++);
				continue;
			}
			//check relative angle difference
			relativeAngleLeft  = linesInLeft[idLeft1][0].direction - linesInLeft[idLeft2][0].direction;
			relativeAngleLeft  = (relativeAngleLeft<M_PI)?relativeAngleLeft:(relativeAngleLeft-TwoPI);
			relativeAngleLeft  = (relativeAngleLeft>(-M_PI))?relativeAngleLeft:(relativeAngleLeft+TwoPI);
			relativeAngleRight = linesInRight[idRight1][0].direction - linesInRight[idRight2][0].direction;
			relativeAngleRight = (relativeAngleRight<M_PI)?relativeAngleRight:(relativeAngleRight-TwoPI);
			relativeAngleRight = (relativeAngleRight>(-M_PI))?relativeAngleRight:(relativeAngleRight+TwoPI);
			relativeAngleDif   = fabs(relativeAngleLeft - relativeAngleRight);
			if((TwoPI-relativeAngleDif)>RelativeAngleDifferenceThreshold&&relativeAngleDif>RelativeAngleDifferenceThreshold){
				eigenMap_.erase(iter++);
				continue;//the relative angle difference is too large;
			}
			iter++;
		}
	}//end while(stillLoop)
	matchResult = matchRet1;
	cout<<"matchRet1.size"<<matchRet1.size()<<", minOfEigenVec_= "<<minOfEigenVec_<<endl;
}

