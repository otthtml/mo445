#include "ift.h"

iftAdjRel *iftAdjacencyFromFile(char *displ_file)
{
  iftAdjRel *A = (iftAdjRel *)calloc(1,sizeof(iftAdjRel));
  FILE *fp     = fopen(displ_file,"r");
  
  fscanf(fp,"%d",&A->n);
  A->dx = iftAllocIntArray(A->n);
  A->dy = iftAllocIntArray(A->n);
  A->dz = iftAllocIntArray(A->n);
  for (int i=0; i < A->n; i++)
    fscanf(fp,"%d %d %d",&A->dx[i],&A->dy[i],&A->dz[i]);
  fclose(fp);

  return(A);
}

iftImage *Dilation(iftImage *img, iftAdjRel *A)
{
  iftImage *res = iftCopyImage(img);

  for (int p=0; p < img->n; p++) {
    iftVoxel u = iftGetVoxelCoord(img,p);
    for (int i=0; i < A->n; i++){
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){
	int q = iftGetVoxelIndex(img,v);
	if (img->val[q] > res->val[p])
	  res->val[p] = img->val[q];
      }
    }
  }

  return(res);
}

iftImage *Erosion(iftImage *img, iftAdjRel *A)
{
  iftImage *res = iftCopyImage(img);

  for (int p=0; p < img->n; p++) {
    iftVoxel u = iftGetVoxelCoord(img,p);
    for (int i=0; i < A->n; i++){
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){
	int q = iftGetVoxelIndex(img,v);
	if (img->val[q] < res->val[p])
	  res->val[p] = img->val[q];
      }
    }
  }

  return(res);
}

iftImage *Closing(iftImage *img, iftAdjRel *A)
{
  iftImage *dil = Dilation(img,A);
  iftImage *res = Erosion(dil,A);
  iftDestroyImage(&dil);
  return(res);
}

iftImage *Opening(iftImage *img, iftAdjRel *A)
{
  iftImage *ero = Erosion(img,A);
  iftImage *res = Dilation(ero,A);
  iftDestroyImage(&ero);
  return(res);
}

int main(int argc, char *argv[]) {
    timer *tstart;

    if (argc!=5)
      iftError("Usage: adjacency <...>\n"
	       "[1] input image (.png, .pgm) \n"
	       "[2] input displacement file to define adjacency (.txt)\n"
	       "[3] operation -- 0: dilation, 1: erosion, 2: closing, 3:opening \n"
	       "[4] output image (.png, .pgm)\n",
	       "main");

    tstart = iftTic();
    
    iftImage  *img   = iftReadImageByExt(argv[1]);
    iftAdjRel *A     = iftAdjacencyFromFile(argv[2]);
    int        oper  = atoi(argv[3]);
    iftImage  *res   = NULL;
    
    switch(oper) {

    case 0: /* Dilation */
      res = Dilation(img,A);
      break;
    case 1: /* Erosion */
      res = Erosion(img,A);
      break;
    case 2: /* Closing */
      res = Closing(img,A);
      break;
    case 3: /* Opening */
      res = Opening(img,A);
      break;
    default:
      iftError("Invalid operation number %d","main",oper);
    }

    iftWriteImageByExt(res,argv[4]);
    iftDestroyImage(&img);
    iftDestroyImage(&res);
    iftDestroyAdjRel(&A);
    
    printf("Done ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
    return (0);
}
