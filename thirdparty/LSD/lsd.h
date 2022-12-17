/*----------------------------------------------------------------------------

  LSD - Line Segment Detector on digital images

  Copyright 2007-2010 rafael grompone von gioi (grompone@gmail.com)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** @file lsd.h
    LSD module header
    @author rafael grompone von gioi (grompone@gmail.com)
 */
/*----------------------------------------------------------------------------*/
#ifndef LSD_HEADER
#define LSD_HEADER


/*----------------------------------------------------------------------------*/
/*----------------------- 'list of n-tuple' data type ------------------------*/
/*----------------------------------------------------------------------------*/
/** 'list of n-tuple' data type

    The i component, of the n-tuple number j, of an n-tuple list 'ntl'
    is accessed with:

      ntl->values[ i + j * ntl->dim ]

    The dimension of the n-tuple (n) is:

      ntl->dim

    The number of number of n-tuples in the list is:

      ntl->size

    The maximum number of n-tuples that can be stored in the
    list with the allocated memory at a given time is given by:

      ntl->max_size
 */
typedef struct ntuple_list_s
{
  unsigned int size;
  unsigned int max_size;
  unsigned int dim;
  double * values;
} * ntuple_list;

void free_ntuple_list(ntuple_list in);
ntuple_list new_ntuple_list(unsigned int dim);


/*----------------------------------------------------------------------------*/
/*----------------------------- Image Data Types -----------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** char image data type

    The pixel value at (x,y) is accessed by:

      image->data[ x + y * image->xsize ]

    with x and y integer.
 */
typedef struct image_char_s
{
  unsigned char * data;
  unsigned int xsize,ysize;
} * image_char;

void free_image_char(image_char i);
image_char new_image_char(unsigned int xsize, unsigned int ysize);
image_char new_image_char_ini( unsigned int xsize, unsigned int ysize,
                               unsigned char fill_value );

/*----------------------------------------------------------------------------*/
/** int image data type

    The pixel value at (x,y) is accessed by:

      image->data[ x + y * image->xsize ]

    with x and y integer.
 */
typedef struct image_int_s
{
  int * data;
  unsigned int xsize,ysize;
} * image_int;

void free_image_int(image_int i);
image_int new_image_int(unsigned int xsize, unsigned int ysize);
image_int new_image_int_ini( unsigned int xsize, unsigned int ysize,
                             int fill_value );

/*----------------------------------------------------------------------------*/
/** double image data type

    The pixel value at (x,y) is accessed by:

      image->data[ x + y * image->xsize ]

    with x and y integer.
 */
typedef struct image_double_s
{
  double * data;
  unsigned int xsize,ysize;
} * image_double;

void free_image_double(image_double i);
image_double new_image_double(unsigned int xsize, unsigned int ysize);
image_double new_image_double_ini( unsigned int xsize, unsigned int ysize,
                                   double fill_value );


/*----------------------------------------------------------------------------*/
/*-------------------------- Line Segment Detector ---------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* LSD Full Interface                                                         */
/*----------------------------------------------------------------------------*/
/** LSD Full Interface

    @param image       Input image.

    @param scale       When different than 1.0, LSD will scale the image by
                       Gaussian filtering.
                       Example: is scale=0.8, the input image will be subsampled
                       to 80% of its size, and then the line segment detector
                       will be applied.
                       Suggested value: 0.8

    @param sigma_scale When scale!=1.0, the sigma of the Gaussian filter is:
                       sigma = sigma_scale / scale,   if scale <  1.0
                       sigma = sigma_scale,           if scale >= 1.0
                       Suggested value: 0.6

    @param quant       Bound to the quantization error on the gradient norm.
                       Example: if gray level is quantized to integer steps,
                       the gradient (computed by finite differences) error
                       due to quantization will be bounded by 2.0, as the
                       worst case is when the error are 1 and -1, that
                       gives an error of 2.0.
                       Suggested value: 2.0

    @param ang_th      Gradient angle tolerance in the region growing
                       algorithm, in degrees.
                       Suggested value: 22.5

    @param eps         Detection threshold, -log10(NFA).
                       The bigger, the more strict the detector is,
                       and will result in less detections.
                       (Note that the 'minus sign' makes that this
                       behavior is opposite to the one of NFA.)
                       The value -log10(NFA) is equivalent but more
                       intuitive than NFA:
                       - -1.0 corresponds to 10 mean false alarms
                       -  0.0 corresponds to 1 mean false alarm
                       -  1.0 corresponds to 0.1 mean false alarms
                       -  2.0 corresponds to 0.01 mean false alarms
                       .
                       Suggested value: 0.0

    @param density_th  Minimal proportion of region points in a rectangle.
                       Suggested value: 0.7

    @param n_bins      Number of bins used in the pseudo-ordering of gradient
                       modulus.
                       Suggested value: 1024

    @param max_grad    Gradient modulus in the highest bin. For example,
                       for images with integer gray levels in [0,255],
                       the maximum possible gradient value is 255.0.
                       Suggested value: 255.0

    @param region      Optional output: an int image where the pixels used
                       in some line support region are marked. Unused pixels
                       have the value '0' while the used ones have the
                       number of the line segment, numbered 1,2,3,...
                       If desired, a non NULL pointer to an image_int should
                       be used. The resulting image has the size of the image
                       used for the processing, that is, the size of the input
                       image scaled by the given factor 'scale'.
                       Suggested value: NULL

    @return            A 5-tuple list, where each 5-tuple corresponds to a
                       detected line segment. The five values are:
                       - x1,y1,x2,y2,width
                       .
                       for a line segment from (x1,y1) to (x2,y2) and
                       a width 'width'.
 */
ntuple_list LineSegmentDetection( image_double image, double scale,
                                  double sigma_scale, double quant,
                                  double ang_th, double eps, double density_th,
                                  int n_bins, double max_grad,
                                  image_int * region );

/*----------------------------------------------------------------------------*/
/* LSD Simple Interface with Scale                                            */
/*----------------------------------------------------------------------------*/
/** LSD Simple Interface with Scale

    @param image Input image.

    @param scale When different than 1.0, LSD will scale the image by
                 Gaussian filtering.
                 Example: is scale=0.8, the input image will be subsampled
                 to 80% of its size, and then the line segment detector
                 will be applied.
                 Suggested value: 0.8

    @return a 5-tuple list of detected line segments.
 */
ntuple_list lsd_scale(image_double image, double scale);

/*----------------------------------------------------------------------------*/
/* LSD Simple Interface                                                       */
/*----------------------------------------------------------------------------*/
/** LSD Simple Interface

    @param image Input image.

    @return a 5-tuple list of detected line segments.
 */
ntuple_list lsd(image_double image);

#endif /* !LSD_HEADER */
/*----------------------------------------------------------------------------*/
