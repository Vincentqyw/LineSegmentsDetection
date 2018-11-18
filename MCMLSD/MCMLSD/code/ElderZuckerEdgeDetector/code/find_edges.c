
/************************************************************************************************
* 
*  Routine Name: lvedge - Detect luminance edges by selecting triads of significant 
*  2nd derivatives of opposite sign flanking a significant 1st derivative of consistent direction.
* 
*       Purpose: vgtriad
*		 Detect luminance edges by selecting triads of significant
*		 2nd derivatives of opposite sign flanking a directional
*		 maximum in the 1st derivative of consistent direction.
*
*		 Pixels with maximal gradient magnitude in the gradient
*		 direction are considered as candidate edge pixels.  The image
*		 is scanned along the gradient line for pixels with significant
*		 2nd derivative of opposite sign, consistent with the edge
*		 hypothesis.  Scanning is terminated if the line intersects a
*		 pixel with inconsistent or unreliable gradient direction, or
*		 inconsistent 2nd derivative sign if a pixel of consistent sign
*		 has already been detected.  Otherwise, the scan selects the
*		 line segment connecting the edge pixel with the maximum 2nd
*		 derivative pixels.  The local blur width of the blurred
*		 edge is derived as well.  The width is derived from the
*		 distance between G2 terminators of opposite sign, based on an
*		 assumption of Gaussian blur.
*
*         Input: img - Original byte image 
*                g1mag - Multi-scale gradient magnitude map 
*                g1dir - Multi-Scale Gradient direction map
*                g1scale - Gradient scale map
*                g2 - Multi-Scale 2nd Derivative Map
*                g2scale - G2 Scale Map
*                noise - Estimated sensor noise
*                maxwidth - Maximum edge width detected
*        Output: gewidth - Estimate of blur scale based on edgewidth. 
*                          Type is VFF_TYP_FLOAT.
*                dark - Estimate of luminance on dark side of edge. 
*                       Type is VFF_TYP_1_BYTE.
*                light - Estimate of luminance on light side of edge
*       Returns: TRUE (1) on success, FALSE (0) on failure
*  Restrictions: Restricted to real, single-band images.  
*                Does not work on explicit location data.
*    Written By: James Elder
*          Date: Oct 13, 1997
*      Verified: 
*  Side Effects: 
*      Examples: 
* Modifications: 
****************************************************************/

#define  char16_t uint16_T /*JE Sept 2014.  Needed for MATLAB 2012a under OSX19.0 and Xcode 5.1 */
#include "mex.h"
#include <stddef.h>  /* NULL */
#include <math.h>
#include <stdio.h>
#define notDblMtx(it) (mxIsComplex(it))


double 	g2max1,    /* Max G2 response of correct sign */
      	g2max2,    /* Max G2 response of incorrect sign */
      	g2maxscale,  /* Scale of max G2 response of correct sign */
		xfirst, yfirst, zfirst, /* Endpoints for interpolation */
      	/* Sensor noise threshold */
      	k_thresh,  /* Min. contrib. to linear interp. for pixel to be labeled */
      	pixerr[6] = {9.5,8,3.5,0.25,0.04,0.0}, /* Half-range of gradient direction */
                                               /* uncertainty due to pixelation */
                                               /* as function of G1 scale. */
      	g2thresh  = 0.0000596, /* Min on sig. g2 value for interpolation */
      	edgewidth, 	/* Estimated width of edge */
      	xend, yend; /* Estimated gradient ray termination coordinates */

int   	g1len_to_max,maxg1len, 	/* Current and max number of G1 pixels supporting edge assertion */ 
      	g1max,               	/* Is current pixel directional max of gradient? */
      	g2len,                 	/* Number of G2 pixels supporting edge */
      	rows, cols, 			/* Columns in image */
      	g1offsets[1024], 		/* Offsets to G1 pixels in current path */
      	g2offsets[4], 			/* Offsets to G2 pixels in current path */
      	edgescale = 0, 			/* Max scale of G1 or G2 response in path */
      	acuity_problem, 		/* Flag to indicate 1-pixel acuity will not be obtained */
      	zcross,                 /* Flag to indicate zero-crossing found */
      	keepoff1,keepoff2; 		/* Offsets to pixels from which 2nd derivative calculated */
double	G2fact1, G2fact2;		/* Factor in recalculation of G2val: "ki * pow(cos(ddiri),2)" */
int		Xx1,Xx2,Yy1,Yy2;		/* Offsets returned from g2_path for subpixel method 2 */
      
double KPI = 3.141593;
double noise; /* Sensor noise threshold */
FILE *f; 

/* Internal functions */
int round_int(double rd);
int g2_path(double theta,int bx,int by,double *g1mag,double *g1dir,double *g1scale,
            double *g2,double *g2scale,int g2sign);
int compare_dir(double g1mag1,double g1dir1,int g1scale1,double g1mag2,double g1dir2,int g1scale2);
int sign(double r);

int
find_edges(double *img,double *g1mag,double *g1dir,double *g1scale,double *g2,
		double *g2scale,double *g2all,double sensor_noise,int maxwidth,int sublocflag,int maxscale,
		double *xout,double *xout_blur,
		double *xout_dark,double *xout_light,double *xout_x1zero,double *xout_y1zero,
		double *xout_x2zero,double *xout_y2zero)
{

double 	theta; /* Either in gradient direction or against gradient direction */
int 	x,y,i; /* Image dimensions and indices */
double 	*g1magd,*g1dird,*g2d; /* Input G1 and G2 map data */
double 	*imgd;
double 	g2maxp,g2maxn; /* Used to detect false edge pixels */
double 	g2scalep,g2scalen; /* Scale of max pos and min neg G2 responses */
double 	blurvar;
double 	nxend,nyend,pxend,pyend; /* Estimated gradient line terminations */                         
double 	nwidth,pwidth; /*distance from current pixel to g2 extrema */
double 	erf_const,C,D,L,meanl;
double 	l,d;

double 	nxfirst, nyfirst, pxfirst, pyfirst, pzfirst, nzfirst, xinter, yinter;
int		g2sc1, g2sc2, g2sc3, g2sc4, off_side1, off_side2, off_side3, off_side4;
int		gind;
double	G2factor1, G2factor2, G2factor3, G2factor4;
int		X1,X2,X3,X4,Y1,Y2,Y3,Y4,match_flag;
int		newoff1,newoff2,newoff3,newoff4;
double	new_nz,new_pz,bestsnr,bestnz,bestpz;
static double g2norms[6] = {1.873,0.2443,0.0306,0.003817,0.00047715,0.0000596};  /* Thresholds */

int 	count=0; 

g1magd = g1mag;
g1dird = g1dir;

g2d 	= g2;
imgd 	= img; 

noise 		= sensor_noise;
g2thresh 	= 5.2*g2thresh*noise; /* threshold on interpolated g2 values */
maxg1len 	= maxwidth*2;  /* Up to 2 pixels per step */
k_thresh 	= 0.1; /* Min coeff. in linear interp. to include pixel in edge */
erf_const 	= 0.68268949213709;

for (i=0;i<5;i++) /* Convert pixelation uncertainty from deg. to radians */
  pixerr[i] = pixerr[i] * KPI / 180.0;

/**** Order (x,y) reversed from cantata code ****/
for (x=0;x<cols;x++)
for (y=0;y<rows;y++,g1magd++,g1dird++,g2d++,imgd++) /* For each pixel in image */
{

	/* These values are updated in g2_path */
	
   	theta = *g1dird; /* Scan in gradient direction */

   	if (*g1magd > 0.0 && *g2d >= 0.0) /* Non-negative 2nd derivative */
   	{
        
		acuity_problem	= 0;
      	g2max1 			= 0;
		g2max2 			= 0.0;
	  	g2maxscale 		= 0.0;
	  	g1len_to_max 	= 0;
		g2len 			= 0;
	  	zcross			= 0;
	  	edgewidth 		= 0.0;
 
		xfirst = 0.0;
		yfirst = 0.0;
		zfirst = 0.0;
		keepoff1	= 0;
		keepoff2	= 0;
		G2fact1		= 0;
		G2fact2		= 0;

		/* Find g2_path, last parameter '1' makes g2_path selective for negative G2 responses */
	  	g2_path(theta,x,y,g1magd,g1dird,g1scale+count,g2d,g2scale+count,1);

	  	if (g2max1 != 0.0 && g1max) /* Gradient max and found significant negative G2 response */
	  	{

			nxfirst		= xfirst;
			nyfirst		= yfirst;
			nzfirst		= -zfirst;

			off_side1	= keepoff1;
			off_side2	= keepoff2;
			G2factor1	= G2fact1;
			G2factor2	= G2fact2;
			X1			= Xx1;
			X2			= Xx2;
			Y1			= Yy1;
			Y2			= Yy2;

	   		nxend 		= xend;       	/* Save estimated edge ray termination */
	     	nyend 		= yend;
	     	nwidth 		= sqrt(pow(nxend,2) + pow(nyend,2));
	     	g2maxn 		= g2max1;    	/* Save min negative G2 response */
	     	g2maxp 		= g2max2;    	/* Save max positive G2 response */
	     	g2scalen 	= g2maxscale; 	/* Save scale of min negative G2 response */
	     	g2max1 		= 0.0;
			g2max2 		= 0.0; 
	     	g2maxscale 	= 0.0;

	     	if (theta > 0)      /* Now scan opposite to gradient direction */
		   		theta = theta - KPI;
	     	else
		   		theta = theta + KPI;
	    
			xfirst = 0.0;
			yfirst = 0.0;
			zfirst = 0.0;
			keepoff1	= 0;
			keepoff2	= 0;
			G2fact1		= 0;
			G2fact2		= 0;

		 	/* Find g2_path, last parameter '0' makes g2_path selective for positive G2 responses */
	     	g2_path(theta,x,y,g1magd,g1dird,g1scale+count,g2d,g2scale+count,0);

	     	pxend 	= xend;   /* Save estimated edge ray termination */
	     	pyend 	= yend;
	     	pwidth 	= sqrt(pow(pxend,2) + pow(pyend,2));

	     	if ((g2max1>g2maxp) && (g2max2>g2maxn) && g1max)
	     	{
				pxfirst		= xfirst;
				pyfirst		= yfirst;
				pzfirst		= zfirst;

				off_side3	= keepoff1;
				off_side4	= keepoff2;
				G2factor3	= G2fact1;
				G2factor4	= G2fact2;
				X3			= Xx1;
				X4			= Xx2;
				Y3			= Yy1;
				Y4			= Yy2;
	
	        	if ((!acuity_problem && zcross) || (acuity_problem && fabs(pwidth-nwidth)<sqrt(2)))
			    /* Pixel lies between min negative and max positive G2 responses */
		        {
		       		*(xout+count) 	= 255;     
	           		g2scalep 		= g2maxscale; 	/* Save scale of max positive G2 response */
		       		g2scalep 		= pow(2.0,g2scalep)/4.0; 	/* convert integer scale to sigma */
		       		g2scalen 		= pow(2.0,g2scalen)/4.0;
		
	           		blurvar = pow(edgewidth/2.0,2.0) - pow((g2scalen + g2scalep)/2.0,2.0);
		       		if (blurvar > 0.0)
		         		*(xout_blur+count)  = sqrt(blurvar); 
		       		else   
		         		*(xout_blur+count)  = 0.1; 

					l = *(imgd+round_int(nxend)*rows+round_int(nyend));
					d = *(imgd+round_int(pxend)*rows+round_int(pyend));

		      		if (l<=d)
		      		{
	    				meanl = round_int(0.5*(l+d));
		        		*(xout_dark+count)	= meanl;
		        		*(xout_light+count) = meanl;
		      		}
		      		else if (blurvar < 1)
		      		{
		        		*(xout_dark+count)	= round_int(d);
		        		*(xout_light+count)	= round_int(l);
		      		}
		      		else
		     		{
		       			C	= (l - d) / erf_const;
		       			D	= d - 0.5 * C * (1 - erf_const);
		       			L	= C + D;

		       			if (D<0)
			     			*(xout_dark+count) = 0;
		       			else if (D> 255)
			     			*(xout_dark+count)  = 255.0;
		       			else
			     			*(xout_dark+count) = round_int(D);
		       			if (L<0)
			     			*(xout_light+count)= 0;
		       			else if (L> 255)
			     			*(xout_light+count) = 255.0;
		       			else
			     			*(xout_light+count) = round_int(L);
		     		}

					if (sublocflag)
					{

						/* METHOD 1: Interpolate to find better estimate of zero-crossing. 
					   		Store resultant x and y in *xout_xzero and *xout_yzero. */
						xinter	= (nzfirst * (pxfirst+x+1) + pzfirst * (nxfirst+x+1)) / (nzfirst + pzfirst);
						yinter	= (nzfirst * (pyfirst+y+1) + pzfirst * (nyfirst+y+1)) / (nzfirst + pzfirst);

						/* Update in xout_x1zero and xout_y1zero matrices */
						xinter	= (xinter > cols) ? cols : xinter;
						*(xout_x1zero+count) = (xinter < 0) ? 0 : xinter;

						yinter	= (yinter > rows) ? rows : yinter;
						*(xout_y1zero+count) = (yinter < 0) ? 0 : yinter;


						/* METHOD 2: Interpolate as with Method 1, except that scale at the
						  	4 flanking pixels must match. */
						g2sc1 = (int)*(g2scale+count+off_side1);
						g2sc2 = (int)*(g2scale+count+off_side2);
   						g2sc3 = (int)*(g2scale+count+off_side3);
   						g2sc4 = (int)*(g2scale+count+off_side4);

						if ((g2sc1==g2sc2) && (g2sc1==g2sc3) && (g2sc1==g2sc4)) /* Scales match */
						{
							*(xout_x2zero+count) = *(xout_x1zero+count);
							*(xout_y2zero+count) = *(xout_y1zero+count);
						}
						else  /* Scales do not match */
						{
							match_flag  = 0;
							gind		= 1;  /* g2 scale index */

							bestsnr		= 0;

							while (!match_flag && (gind<=maxscale))
							{
								newoff1	= (x+X1) * rows * maxscale + (y+Y1) + rows * (gind-1);
								newoff2	= (x+X2) * rows * maxscale + (y+Y2) + rows * (gind-1);
								newoff3	= (x+X3) * rows * maxscale + (y+Y3) + rows * (gind-1);
								newoff4	= (x+X4) * rows * maxscale + (y+Y4) + rows * (gind-1);

      							new_nz	= -G2factor1 * (*(g2all+newoff1)) - G2factor2 * (*(g2all+newoff2)); 
      							new_pz	=  G2factor3 * (*(g2all+newoff3)) + G2factor4 * (*(g2all+newoff4)); 

								/* If there is no scale for which all 4 are significant, will need to
								 * choose the one with the best SNR, therefore store it. */
								if ((new_nz+new_pz)>bestsnr)
								{
									bestsnr	= new_nz + new_pz;
									bestnz	= new_nz;
									bestpz	= new_pz;
								}

								/* If new estimates exceed the threshold, have found answer */
								if ((new_nz > 5.2*noise*g2norms[gind-1]) && 
											(new_pz > 5.2*noise*g2norms[gind-1]))
								{
									xinter	= (new_nz * (pxfirst+x+1) + new_pz * (nxfirst+x+1)) / 
														(new_nz + new_pz);
									yinter	= (new_nz * (pyfirst+y+1) + new_pz * (nyfirst+y+1)) / 
														(new_nz + new_pz);

									match_flag = 1;
								}

								gind++;
							}

							/* Use scale with best SNR if none match */
							if (!match_flag)
							{
								xinter	= (bestnz * (pxfirst+x+1) + bestpz * (nxfirst+x+1)) / bestsnr;
								yinter	= (bestnz * (pyfirst+y+1) + bestpz * (nyfirst+y+1)) / bestsnr;
							}

							/* Update in xout_x2zero and xout_y2zero matrices */
							xinter	= (xinter > cols) ? cols : xinter;
							*(xout_x2zero+count) = (xinter < 0) ? 0 : xinter;

							yinter	= (yinter > rows) ? rows : yinter;
							*(xout_y2zero+count) = (yinter < 0) ? 0 : yinter;

						}

					}

	      		}
	    	}
	  	}
    }
    count++;
}
return(1);
}


/****************************************************************/
int round_int(double rd)
{ 
double t1, t2;
int	xf, xc;
  
xf = (int)floor(rd);
xc = (int)ceil(rd);
t1 = fabs(rd-xf);
t2 = fabs(rd-xc);

if (t1 >= t2)
	return(xc);
else
	return(xf);

}

/*****************************************************************/
/* g2_path */
int g2_path(double theta,int bx,int by,double *g1mag,double *g1dir,double *g1scale,
					double *g2,double *g2scale,int g2signbit)
{
int 	x1, y1, x2, y2, dx, dy; 
int		off1=0; 
int		off2=0;  /* Offsets and indices for linear interp. */
int 	end_of_path; /* Flag signalling when max (min) found or continuity broken */
int 	same_dir1, same_dir2; 	/* Flags indicating gradient direction compatibility */
int 	g1len, g2halflen; 		/* Num of G2 pixels terminating current ray */
double 	tan_theta, x, y; 		/* Real coords of points on path */
double 	g1magval, g1dirval; 	/* Gradient at centre pixel */
double 	g1mag1, g1dir1, g1mag2, g1dir2; /* Gradients at pixel pair on path */
double 	k1, k2; 				/* Linear interp. coefficients */
double 	g1val, g2val; 			/* Interpolated G1 and G2 responses on path */
int 	g1scaleval; 			/* Scale of gradient at centre pixel */
int 	g1scale1, g1scale2, g2scale1, g2scale2; /* Scales of pixels on path */
double 	ddir1, ddir2; 			/* gradient dir differences along grad line. */
double 	g1_prev_val, acuity, length;
int		found_interp_flag=0;

/*#########################################*/
/* Select the 1st 2 pixels flanking the ray 
	from 0 in the theta direction. */
/*#########################################*/

if (theta > 0)
  	dy = -1;
else
  	dy = 1;

if (fabs(theta) < KPI*0.5)
  	dx = 1;
else
  	dx = -1;
 
x1 = dx;
y1 = 0;
x2 = 0;
y2 = dy;

tan_theta = tan(theta);

/* If theta < 45deg -> (x2,y2) = (dx,dy), else (x1,y1) = (dx,dy) */
if (fabs(tan_theta) < 1)
  	x2 = dx;
else
  	y1 = dy;

/* End of block */
/*##############*/

acuity 		= 20.0;
g2max2 		= 0.0;
g1magval 	= (float)*g1mag;
g1dirval 	= (float)*g1dir;
g1scaleval 	= (int)*g1scale;
g2val 		= (float)*g2;
g1_prev_val = g1magval;

if (g2val != 0.0 && sign(g2val) == g2signbit) /* sig. G2 of right sign */
{
	g2maxscale 	= (float) *g2scale;
    g2max1 		= g2val;
    g2offsets[g2len++] = 0;
    g2halflen 	= 1;
	if (!found_interp_flag)
	{
		G2fact1 = 1;
		G2fact2 = 0;
		Xx1	= 0;
		Xx2	= 0;
		Yy1	= 0;
		Yy2	= 0;
	
		zfirst	= g2max1;
		xfirst	= 0;
		yfirst	= 0;
		keepoff1 = off1;
		keepoff2 = off2;
		found_interp_flag = 1;
	}
}
else /* Initialize G2 variables */
{
      g2max1 		= 0.0;
      g2maxscale 	= 0.0;
      g2halflen 	= 0;
}

end_of_path = 0;
xend = yend = 0.0;
g1len 		= g1len_to_max;
g1max 		= 1;
 

/*#####################################
 Scan for G2 max of the correct sign, 
   until length limit or interrupted 
#####################################*/
while (!end_of_path && g1len < maxg1len) 
{
    
	if ((g2val==0.0) && (g2max1 == 0))
		acuity_problem = 1;

	/* Check if flanking pixels go out of range of image */
    if ((bx+x1<0) || (bx+x1>=cols) || (bx+x2<0) || (bx+x2>=cols) || (by+y1<0) || 
						(by+y1>=rows) || (by+y2<0) || (by+y2>=rows))
      	end_of_path = 1;
    else
    {
      	/* Derive linear interpolation coefficients */
      	if (x1 == x2)
        {
        	y = -x1 * tan_theta;
            k1 = fabs(y-y2);
            k2 = fabs(y-y1);
        }
      	else
        {
            x = -y1 / tan_theta;
            k1 = fabs(x-x2);
            k2 = fabs(x-x1);
        }
      
		/*########################################################*/
      	/* Get G1 and G2 information for the two pixels on path */
		/*########################################################*/

		/**** cols and rows changed from cantata code here as Matlab processes columnwise ****/
      	off1 	= x1 * rows + y1;
      	off2 	= x2 * rows + y2;

      	x 		= k1 * x1 + k2 * x2;
      	y 		= k1 * y1 + k2 * y2;
      	length 	= sqrt(x*x + y*y);


      	g1dir1 	= *(g1dir+off1);
      	g1dir2 	= *(g1dir+off2);
      	g1mag1 	= *(g1mag+off1);
      	g1mag2 	= *(g1mag+off2);

      	g1scale1 = (int)*(g1scale+off1);
      	g1scale2 = (int)*(g1scale+off2);
      	g2scale1 = (int)*(g2scale+off1);
      	g2scale2 = (int)*(g2scale+off2);

      	/* Interpolate the G1 estimate */
      	g1val = k1 * g1mag1 + k2 * g1mag2; 
      
      	if (g1val > g1magval)
        	g1max = 1;
      
      	/* Pixel 1 compatible? */ 
      	same_dir1 = compare_dir(g1magval,g1dirval,(int)g1scaleval,g1mag1,g1dir1,(int)g1scale1); 
      
      	/* Pixel 2 compatible? */                      
      	same_dir2 = compare_dir(g1magval,g1dirval,(int)g1scaleval,g1mag2,g1dir2,(int)g1scale2); 
                                                                                                  
	  	ddir1 = fabs(g1dir1-g1dirval);
	  	ddir2 = fabs(g1dir2-g1dirval);

		/* Interpolate the G2 estimate */
      	g2val = k1 * pow(cos(ddir1),2) * (*(g2+off1)) + k2 * pow(cos(ddir2),2) * (*(g2+off2)); 
            
      	if ((g1len==0) && (g2val<0.0))
	  		zcross = 1;

      	if ((!(g2max1==0.0 && fabs(g2val)>g2thresh && sign(g2val) == g2signbit)   &&
	   			!((same_dir1 && k1>k_thresh) || (same_dir2 && k2>k_thresh)))
          		|| (fabs(g2val) > g2thresh && (sign(g2val) != g2signbit) && zcross)
          		|| ((fabs(g2val) < g2thresh) && (length > acuity)))
		{
        	end_of_path = 1;
      	}
      	else if ((sign(g2val) == g2signbit) && (fabs(g2val) > fabs(g2max1)))
        {
			g2max1 = g2val;
		   	/* Derive scale of G2 response */
           	if (g2scale1 == 0)
        		g2maxscale = (float) g2scale2;
           	else if (g2scale2 == 0)
          		g2maxscale = (float) g2scale1;
    		else
                g2maxscale = k1 * g2scale1 + k2 * g2scale2;

			g1len_to_max = g1len; 		/* Update g1 index */
            g2len 		-= g2halflen; 	/* Erase previous G2 terminators */
            g2halflen 	 = 0;    

			if (k1 > k_thresh) /* If path lies near this pixel, label as G2 */
            {
            	g2offsets[g2len++] = off1;
				g2halflen++;
            }
            if (k2 > k_thresh)
            {
				g2offsets[g2len++] = off2;
                g2halflen++;
            }

			xend = k1 * x1 + k2 * x2; /* Real coordinates of G2 terminator */
            yend = k1 * y1 + k2 * y2;

			if (!found_interp_flag)
			{
				G2fact1	= k1 * pow(cos(ddir1),2);
				G2fact2	= k2 * pow(cos(ddir2),2);
				Xx1	= x1;
				Xx2	= x2;
				Yy1	= y1;
				Yy2	= y2;
	
				zfirst	= g2max1;
				xfirst	= xend;
				yfirst	= yend;
				keepoff1 = off1;
				keepoff2 = off2;
				found_interp_flag = 1;
			}

		}
      	else if (((sign(g2val) != g2signbit) || g2val == 0.0) && g2max1 != 0.0) 
      	{
        	end_of_path = 1; 
      	}
      	else if ((sign(g2val) != g2signbit) && (fabs(g2val) > fabs(g2max2)))
      	{
         	g2max2 = g2val; /* Largest sig. G2 response of incorrect sign */
      	}                    
    
      	g1_prev_val = g1val;

	}

    if (!end_of_path)
    {
    	if (same_dir1) /* 1st pixel in pair is compatible */
        {
       		g1offsets[g1len++] = off1; /* Label as interior pixel */
            if (g1scale1 > edgescale)  /* Update scale to max */
            	edgescale = g1scale1;
		}
        if (same_dir2) /* 2nd pixel in pair is compatible */
        {
	        g1offsets[g1len++] = off2; /* Label as interior pixel */
    		if (g1scale2 > edgescale)  /* Update scale to max */
           		edgescale = g1scale1;
			/**** Is this meant to be g1scale1 or g1scale2? ****/
        }

        /* Step to next intersection of ray with pixel grid */
        if (x1 == x2)
        {
        	x1 += dx;
            if (abs(y2) < fabs(x1*tan_theta))
            	y1+= dy;
            else
                x2+= dx;
      	}
        else
        {
        	y2 += dy;
            if (abs(y2) > fabs(x1*tan_theta))
            	x2+= dx;
            else
                y1+= dy;
        }
	}




}

if (g2maxscale > edgescale) /* Update edgescale to max of G1 and G2 scales */
	edgescale = (int) g2maxscale;
edgewidth += sqrt(xend*xend + yend*yend); /* Add distance from 0 to terminator */

return (1);
}      


/******************************************************************************/
/* compare_dir */
int compare_dir(double g1mag1,double g1dir1,int g1scale1,double g1mag2,double g1dir2,int g1scale2)
{
double a1, a2; /* Pixelation uncertainty at respective scales */
double ddir, tolerance; /* actual and tolerated gradient direction discrepancy */
static double g1noise_norms[6] = {0.765,0.199,0.0499,0.0125,0.00312,0.00078};

if (g1dir2 == 4) /* Gradient not significant */
    return(0);
 
a1 	= pixerr[g1scale1-1];
a2 	= pixerr[g1scale2-1];

ddir = fabs(g1dir1 - g1dir2);
if (ddir > KPI)
  ddir = 2*KPI - ddir;
 
/* tolerance is sum of half-ranges due to pixelation uncertainty 
		plus the range due to sensor noise */
tolerance = a1 + a2 + 5.2 * noise * sqrt(pow(g1noise_norms[g1scale1-1]/g1mag1,2.0)
                                   + pow(g1noise_norms[g1scale2-1]/g1mag2,2.0));
                                                                   
if (ddir <= tolerance)
  	return(1);

return(0);
}


/***************************************************************************/
/* signbit */
int sign(double r)
{
	if (r<0)
    	return (1);
  	else
    	return (0);
}



/******************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

double *img;
double *g1mag;
double *g1dir; 
double *g1scale;
double *g2; 
double *g2scale;
double *g2all;
int x_dim, y_dim;
double sensor_noise;
int maxwidth, maxscale; 
int sublocflag; 
double dmwidth,mxscale;
double subflag;
double *x;
double *x_blur;
double *x_dark;
double *x_light;
double *x_x1zero;
double *x_y1zero;
double *x_x2zero;
double *x_y2zero;

int warnings = 1;
mxArray *arg;
double *mxMat;

/* copy start */


	if (nrhs < 11 ) 
		mexErrMsgTxt("requres  at least 11 args.");
 
  	/* ARG 1: IMAGE  */
  	arg = prhs[0];
  
  	if notDblMtx(arg) mexErrMsgTxt("IMAGE arg must be a real non-sparse matrix.");
  	img = mxGetPr(arg);
  	x_dim = (int) mxGetM(arg); /* X is inner index! */
  	y_dim = (int) mxGetN(arg);
  	rows = x_dim;
  	cols = y_dim;
  
  	/* ARG 2: g1mag */
  	arg = prhs[1];
  	if notDblMtx(arg) mexErrMsgTxt("g1mag arg must be a real non-sparse matrix.");
  	g1mag = mxGetPr(arg);
  

  	/* ARG 3: g1dir */
  	arg = prhs[2];
  	if notDblMtx(arg) mexErrMsgTxt("g1dir arg must be a real matrix");
  	g1dir = mxGetPr(arg);

  	/* ARG 4: g1scale */
  	arg = prhs[3];
  	if (notDblMtx(arg)) mexErrMsgTxt(" g1scale arg must be a int matrix");
  	g1scale = mxGetPr(arg);
  
  	/* ARG 5: g2 */
  	arg = prhs[4];
  	if notDblMtx(arg) mexErrMsgTxt("g2 arg must be a real matrix");
  	g2 = mxGetPr(arg);
  
   	/* ARG 6: g2scale */
  	arg = prhs[5];
  	if (notDblMtx(arg)) mexErrMsgTxt("g2scale arg must be a int matrix");
  	g2scale = mxGetPr(arg);
  
   	/* ARG 7: g2all */
  	arg = prhs[6];
  	if (notDblMtx(arg)) mexErrMsgTxt("g2all arg must be a int matrix");
  	g2all = mxGetPr(arg);
  
   	/* ARG 8: sensor_noise */
  	arg = prhs[7];
  	if notDblMtx(arg) mexErrMsgTxt("sensor_noise arg must be a real matrix");
  	if (mxGetM(arg)*mxGetN(arg) != 1)
      	mexErrMsgTxt("sensor_noise arg must be a scalar");
  	mxMat = mxGetPr(arg);
  	sensor_noise = *mxMat;
  	
  	/* ARG 9: maxwidth */
  	arg = prhs[8];
  	if (notDblMtx(arg)) mexErrMsgTxt("maxwidth arg must be a int matrix");
  	if (mxGetM(arg)*mxGetN(arg) != 1)
      	mexErrMsgTxt("maxwidth arg must be a scalar");
  	mxMat = mxGetPr(arg);
  	dmwidth = *mxMat;
  
  	/* ARG 10: subpixel */
  	arg = prhs[9];
  	if (notDblMtx(arg)) mexErrMsgTxt("subpixelflag arg must be a int matrix");
  	if (mxGetM(arg)*mxGetN(arg) != 1)
      	mexErrMsgTxt("subpixelflag arg must be a scalar");
  	mxMat = mxGetPr(arg);
  	subflag = *mxMat;
  
  	/* ARG 11: maxscale */
  	arg = prhs[10];
  	if (notDblMtx(arg)) mexErrMsgTxt("maxscale arg must be a int matrix");
  	if (mxGetM(arg)*mxGetN(arg) != 1)
      	mexErrMsgTxt("maxscale arg must be a scalar");
  	mxMat = mxGetPr(arg);
  	mxscale = *mxMat;
  
  	/* ARG 6: WARNINGS */
  	if (nrhs>11)
    {
    	arg = prhs[4];
    	if notDblMtx(arg) mexErrMsgTxt("WARINGS arg must be a real scalar.");
    	if (mxGetM(arg) * mxGetN(arg) != 1)
      	mexErrMsgTxt("WARNINGS arg must be a real scalar.");
    	mxMat = mxGetPr(arg);
    	warnings = (int) *mxMat;
    }

  	plhs[0] = mxCreateDoubleMatrix(x_dim,y_dim,mxREAL);
  	if (plhs[0] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
  	x = mxGetPr(plhs[0]); 
  
  	plhs[1] = mxCreateDoubleMatrix(x_dim,y_dim,mxREAL);
  	if (plhs[1] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
  	x_blur = mxGetPr(plhs[1]); 
  
  	plhs[2] = mxCreateDoubleMatrix(x_dim,y_dim,mxREAL);
  	if (plhs[2] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
  	x_dark = mxGetPr(plhs[2]);
  
  	plhs[3] = mxCreateDoubleMatrix(x_dim,y_dim,mxREAL);
  	if (plhs[3] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
  	x_light = mxGetPr(plhs[3]);
  
  	plhs[4] = mxCreateDoubleMatrix(x_dim,y_dim,mxREAL);
  	if (plhs[4] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
  	x_x1zero = mxGetPr(plhs[4]);
  
  	plhs[5] = mxCreateDoubleMatrix(x_dim,y_dim,mxREAL);
  	if (plhs[5] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
  	x_y1zero = mxGetPr(plhs[5]);
  
  	plhs[6] = mxCreateDoubleMatrix(x_dim,y_dim,mxREAL);
  	if (plhs[6] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
  	x_x2zero = mxGetPr(plhs[6]);
  
  	plhs[7] = mxCreateDoubleMatrix(x_dim,y_dim,mxREAL);
  	if (plhs[7] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
  	x_y2zero = mxGetPr(plhs[7]);
  
/* copy ends */

maxwidth = (int)dmwidth;
maxscale = (int)mxscale;
sublocflag = (int)subflag;

find_edges(img,g1mag,g1dir,g1scale,g2,g2scale,g2all,sensor_noise,maxwidth,sublocflag,
			maxscale,x,x_blur,x_dark,x_light,x_x1zero,x_y1zero,x_x2zero,x_y2zero);
}
