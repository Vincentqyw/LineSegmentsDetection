/*
 *  Function:   mexVoteEdges_v2
 *
 *  Input:      Kernel value matrices (regular and flipped), and dimensional
 *              arguments, along with the voting map
 *
 *  Output:     Updated voting map
 *
 *  Author:     Ron Tal
 *
 *
 *  Usage:      Compiled and called through Matlab
 *
 */

/* Include the following header files */
# include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>


#define PI 3.141
#define RMAX 400.0
/* Input Arguments */

#define	EDG_L	prhs[0]
#define	EDG_G	prhs[1]
#define KERN    prhs[2]
#define KERNF   prhs[3]
#define R_RES   prhs[4]
#define T_RES   prhs[5]
#define PP      prhs[6]
#define KMAX    prhs[7]
#define SCALE   prhs[8]


/* Output Arguments */
#define	NMAP	plhs[0]

/* STRUCTS */
typedef struct 
{
    int m;
    int n;
} matPoint;

typedef struct 
{
    double r;
    double th;
} houghPoint;

/* Global variables: */

mwSize mmap, nmap, m_edge_loc, n_edge_loc, m_edge_grad, n_edge_grad, m_kern, m_kern_f;
double* edge_loc, * edge_grad, * pp, r_res, th_res;
double* map,*scale;
mxArray* kernels, *kernels_flip, * k_max;

int round_f (double val)
{
    return floor(val + 0.5);
}
double round2frac(double val, double res)
{
    double lrg = val/res;
    double rnd = (double) round_f(lrg);
    return rnd*res;
}
matPoint hough2mat(double r, double th, double rres, double thres)
{
    matPoint val;
    val.n = round_f((r + (double)RMAX**scale ) / rres  + 1);
    val.m = round_f(th / thres + 1);
    return val;
}
houghPoint mat2hough(int m, int n, double rres, double thres)
{
    houghPoint val;
    val.th = round2frac(m*thres, thres);
    val.r = round2frac(n*rres - (double)RMAX**scale, rres);
    return val;
}
void houghVoteMap(matPoint center, mxArray * k, mxArray * k_f, mxArray * k_m, int x, int y, int r_th, int xd, int yd)
{
    /* Get Kernel properties */
    int i;
    int m = center.m;
    int n = center.n;
    int m_ker, n_max, m_max, n_lim_pos, n_lim_neg, m_lim_pos, m_lim_neg, n_ker;
    double * kernel;
    double * kernel_flip;
    double * k_mx;
    double val;
	double kvalue; 
	int xx, yy;
    m_ker = mxGetM(k);
	n_ker = mxGetN(k);
    kernel = mxGetPr(k);
    
    kernel_flip = mxGetPr(k_f);
    k_mx = mxGetPr(k_m);
    n_max = k_mx[0];
    m_max = k_mx[1];
    n_lim_pos = n + (n_max-1)/2 - 1;
    n_lim_neg = n - (n_max-1)/2 - 1;
    m_lim_pos = m + (m_max-1)/2 - 1;
    m_lim_neg = m - (m_max-1)/2 - 1;
	xx = round_f(x + (double)pp[0]);
	yy = round_f(y + (double)pp[1]);

    /* Vote for each valid kernel value */
    for (i = 0; i < m_ker; i++)
    {
        int m_on_map;
        int n_on_map;
        
        m_on_map = m_lim_neg + (int)kernel[m_ker + i] - 1;
        n_on_map = n_lim_neg + (int)kernel[i] - 1;

        /* Out of R bounds: do nothing */
        if (n_on_map <= 0 ||n_on_map >= nmap) continue; 
        /* Handle wraparound in TH domain: to be added later: */
        else if( m_on_map < 0) 
        {   
            m_on_map = mmap + m_on_map ;
            n_on_map =  nmap - n_on_map ;
            if (n_on_map < 0 ||n_on_map >= nmap) continue;

			kvalue = kernel[2*m_ker + (m_ker - i - 1)]; 
            map[(n_on_map )*mmap + m_on_map] += kvalue;
       
        }else if(m_on_map >= mmap)
        {
            m_on_map = m_on_map - mmap ;
            n_on_map =  nmap - n_on_map ;

			kvalue = kernel[2*m_ker + (m_ker - i)]; 
            map[(n_on_map )*mmap + m_on_map] += kvalue; 
            
        }
        /* No wraparound */
        else
        {
			kvalue = kernel[2*m_ker + i]; 
            map[(n_on_map )*mmap + m_on_map] += kvalue;
        }

    }

    
}


void loopThrough()
{
	int xx, yy;
	int i = 0;
    for (i = 0; i < m_edge_loc;i++ )
    {
         double r, th;
         double * kernel, * kernel_flip, * k_m, * m_l;
		 int r_th;
         matPoint center;
		 int xx, yy, ind;
 		 
         double x = (double)edge_loc[i] - (double)pp[0];
         double y = (double)edge_loc[m_edge_loc + i] - (double)pp[1];
         
		 int xd = round_f((double)edge_loc[i]);
		 int yd = round_f((double)edge_loc[m_edge_loc + i]);

		 th = edge_grad[i];
         if (th == 4.0) continue; 
         if (th < 0) th += PI;
         if (th > PI) th -= PI;
         th = PI - th;
         r = x * cos(th) + y * sin(th);
                  
         th = round2frac(th, th_res);
         r = round2frac(r, r_res);  
         r_th = round_f(- x * sin(th) + y * cos(th));

         if (r_th > (int)(400**scale)) r_th = (int)(400**scale); else if(r_th < - (int)(400**scale)) r_th = -(int)(400**scale);
         
         if (r_th == 0)
         {
             ind = 0;
         } 
         else
         {
             ind = r_th + (double)RMAX**scale;
         }
         kernel = mxGetCell(kernels, ind);
         kernel_flip = mxGetCell(kernels_flip, ind);
         k_m = mxGetCell(k_max, ind);
 
         /* Center of kernel on voting map */
         center = hough2mat(r, th, r_res, th_res);
         houghVoteMap(center,kernel, kernel_flip, k_m, x, y, r_th, xd, yd);        
    }
    
    
}
void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
/*
 *  Function:           mexFunction
 *
 *  Inputs: nrhs ->     Integer representing the number of input arguments
 *                      sent from Matlab.
 *
 *          prhs[] ->   An mxArray structure of size nrhs holding pointers
 *                      to arguments sent from Matlab.
 *
 *          nlhs  ->    Integer representing the number of output arguments
 *                      requested from Matlab.
 *
 *          plhs[] ->   An mxArray structure of size nlhs holding pointers
 *                      to the location of output arguments.
 *
 *  Outputs:            (void)
 *
 *  Functionality:      Analogous to a 'main' function. This function will
 *                      be activated by default when Matlab calls this C
 *                      function. Used to interface between C and Matlab,
 *                      as to ensure output is in the location and format
 *                      Matlab expecs.
 */
{
    /* Check for proper number of arguments */
    if (nrhs != 9) { 
	mexErrMsgTxt("Eight input arguments required.");
    } 
    printf("Inside C function_v2\n");
    /* Get edge locations */
    m_edge_loc = mxGetM(EDG_L); 
    n_edge_loc = mxGetN(EDG_L);
    edge_loc = mxGetPr(EDG_L);
    if (n_edge_loc != 2) { 
        mexErrMsgTxt("Wrong edge locations format.");
    } 
    /* Get edge gradients */
    m_edge_grad = mxGetM(EDG_G); 
    n_edge_grad = mxGetN(EDG_G);
    if (n_edge_grad != 1) { 
        mexErrMsgTxt("Wrong edge gradient format.");
    } 
    if (m_edge_loc != m_edge_grad) { 
        mexErrMsgTxt("Locations and gradients must correspond.");
    } 
    
    edge_grad = mxGetPr(EDG_G);
    /* Get pp */
    pp = mxGetPr(PP);
    r_res = mxGetScalar(R_RES);
    th_res = mxGetScalar(T_RES);
    scale=mxGetPr(SCALE);
    /* Get map*/
    nmap = (int)(800**scale/r_res) + 1;
    mmap = (int)(3.141/th_res) + 1;
    NMAP = mxCreateDoubleMatrix(mmap, nmap, mxREAL);
    map = mxGetPr(NMAP);
    /* Get Kernels*/
    m_kern = mxGetNumberOfElements(KERN); 
    m_kern_f = mxGetNumberOfElements(KERNF);
    kernels = KERN;
    kernels_flip = KERNF;
    k_max = KMAX;
	loopThrough();
    return;
}
/*end mexFunction*/