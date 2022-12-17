#include <stdio.h>
#include "lsd.h"

int main(void)
{
  image_double image;
  ntuple_list out;
  unsigned int x,y,i,j;
  unsigned int X = 128;  /* x image size */
  unsigned int Y = 128;  /* y image size */

  /* create a simple image: left half black, right half gray */
  image = new_image_double(X,Y);
  for(x=0;x<X;x++)
    for(y=0;y<Y;y++)
      image->data[ x + y * image->xsize ] = x<X/2 ? 0.0 : 64.0; /* image(x,y) */

  /* call LSD */
  out = lsd(image);

  /* print output */
  printf("%u line segments found:\n",out->size);
  for(i=0;i<out->size;i++)
    {
      for(j=0;j<out->dim;j++)
        printf("%f ",out->values[ i * out->dim + j ]);
      printf("\n");
    }

  /* free memory */
  free_image_double(image);
  free_ntuple_list(out);

  return 0;
}
