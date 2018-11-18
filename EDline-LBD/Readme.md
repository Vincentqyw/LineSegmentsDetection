This library is an implementation of LBD descriptors developed by Lilian Zhang and Reinhard Koch:
http://www.sciencedirect.com/science/article/pii/S1047320313000874

The original code was written by Lilian Zhang and can be found here:
http://www.mip.informatik.uni-kiel.de/tiki-index.php?page=Lilian+Zhang

To extract edge and lines, this library implements the EDLines Algorithm and the Edge Drawing detector:
http://www.sciencedirect.com/science/article/pii/S0167865511001772
http://www.sciencedirect.com/science/article/pii/S1047320312000831

All original dependencies except Opencv have been removed and the code has been optimized for Opencv 2.4.x
PairWiseLineMatching has not been touched and it still needs original dependencies, because the aim of this porting was to match descriptors with Opencv matchers.

During the 2014 GSoC a binarized version of these descriptors with a dedicated matcher implementation has been included into opencv_contrib:




References:

LBD Descriptors:

[1] Lilian Zhang and Reinhard Koch. 2013. An efficient and robust line segment matching approach based on LBD descriptor and pairwise geometric consistency. J. Vis. Comun. Image Represent. 24, 7 (October 2013), 794-805. DOI=10.1016/j.jvcir.2013.05.006 http://dx.doi.org/10.1016/j.jvcir.2013.05.006 

Edge Drawing (ED):

[1] C. Topal and C. Akinlar, "Edge Drawing: A Combined Real-Time Edge and Segment Detector," Journal of Visual Communication and Image Representation, 23(6), 862-872, DOI: 10.1016/j.jvcir.2012.05.004 (2012). 

[2] C. Topal, C. Akinlar, Y. Genc, Edge Drawing: An Heuristic Approach to Robust Real-Time Edge Detection, Int’l Conf. Pattern Recognition (ICPR), Istanbul, TURKEY, 2010.

EDLines:

[3] C. Akinlar and C. Topal, "EDLines: A Real-time Line Segment Detector with a False Detection Control," Pattern Recognition Letters, 32(13), 1633-1642, DOI: 10.1016/j.patrec.2011.06.001 (2011). 

[4] C. Akinlar, C. Topal, EDLines: Realtime Line Segment Detection by Edge Drawing (ED), IEEE Int’l Conf. Image Processing (ICIP), Brussels, BELGIUM, 2011.
