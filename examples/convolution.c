#include "ift.h"

iftKernel *iftReadKernelFile(char *filename)
{
  FILE *fp = fopen(filename,"r");
  iftKernel *K = (iftKernel *)calloc(1,sizeof(iftKernel));
  K->A = (iftAdjRel *)calloc(1,sizeof(iftAdjRel));
  fscanf(fp,"%d",&K->A->n);
  K->A->dx  = iftAllocIntArray(K->A->n);
  K->A->dy  = iftAllocIntArray(K->A->n);
  K->A->dz  = iftAllocIntArray(K->A->n);
  K->weight = iftAllocFloatArray(K->A->n);
  for (int i=0; i < K->A->n; i++)
    fscanf(fp,"%d %d %d %f",&K->A->dx[i],&K->A->dy[i],&K->A->dz[i],&K->weight[i]);

  fclose(fp);
  return(K);
}


int main(int argc, char *argv[]) {
    timer *tstart;

    if (argc!=4)
      iftError("Usage: convolution <...>\n"
	       "[1] input image  (.png, .pgm) \n"
	       "[2] input filter (.txt)\n"
	       "[3] output image (.png, .pgm)\n",
	       "main");

    tstart = iftTic();
    
    iftImage  *img   = iftReadImageByExt(argv[1]);
    iftKernel *K     = iftReadKernelFile(argv[2]);
    iftImage  *conv  = iftLinearFilter(img,K);
    iftImage  *mag   = iftAbs(conv);
    iftImage  *norm  = iftNormalize(mag,0,255);
    
    iftWriteImageByExt(norm,argv[3]);
    iftDestroyImage(&img);
    iftDestroyImage(&conv);
    iftDestroyImage(&norm);
    iftDestroyImage(&mag);
    iftDestroyKernel(&K);
    
    printf("Done ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
    return (0);
}
