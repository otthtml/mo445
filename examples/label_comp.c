#include "ift.h"

char *Basename(char *path)
{
  char *basename     = iftBasename(path);
  iftSList *slist    = iftSplitString(basename,"/");
  strcpy(basename,slist->tail->elem);
  iftDestroySList(&slist);
  return(basename);
}

float MaxArcWeight(iftMImage *mimg, iftAdjRel *A)
{
  float max_dist=0.0;
  int dx, dy, dz;
  
  iftMaxAdjShifts(A, &dx, &dy, &dz);

  for (int p=0; p < mimg->n; p++) {
    iftVoxel u = iftMGetVoxelCoord(mimg,p);
    if ((u.x >= dx) && (u.y >= dy) &&
	(u.x < (mimg->xsize-dx))&&(u.y < (mimg->ysize-dy))){
      for (int i=1; i < A->n; i++) {
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	int q = iftMGetVoxelIndex(mimg,v);
	float dist=0.0;
	for (int b=0; b < mimg->m; b++){
	  dist += (mimg->val[q][b]-mimg->val[p][b])*
	    (mimg->val[q][b]-mimg->val[p][b]);
	}
	if (dist > max_dist)
	  max_dist = dist;
      }
    }
  }

  return(max_dist);
}

float Distance(iftMImage *mimg, int p, int q)
{
  float dist=0.0;

  for (int b=0; b < mimg->m; b++)
    dist += (mimg->val[q][b]-mimg->val[p][b])*(mimg->val[q][b]-mimg->val[p][b]);
  return(dist);
}
  
iftImage *LabelComp(iftMImage *mimg, iftAdjRel *A, float thres)
{
  iftImage *label=NULL;
  int l=1;
  iftFIFO *F=NULL;
  float max_dist = MaxArcWeight(mimg,A);

  thres = thres*max_dist;
  
  label  = iftCreateImage(mimg->xsize,mimg->ysize,mimg->zsize);
  F      = iftCreateFIFO(mimg->n);

  for (int r=0; r < label->n; r++){
    if (label->val[r]==0){
      label->val[r]=l;
      iftInsertFIFO(F,r);
      while(!iftEmptyFIFO(F)){
	int p = iftRemoveFIFO(F);
	iftVoxel u = iftGetVoxelCoord(label,p);
	for (int i=1; i < A->n; i++){
	  iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(label,v)){
	    int q = iftGetVoxelIndex(label,v);
	    if ((Distance(mimg,r,q)<=thres)&&(label->val[q] == 0)){
	      label->val[q] = label->val[p];
	      iftInsertFIFO(F,q);
	    }
	  }
	}
      }
      l++;
    }
  }
  iftDestroyFIFO(&F);
  
  return(label);
}

/* It can be considered a watershed transform */

void PropagateLabelByFmax(iftMImage *mimg, iftImage *label, int minarea_thres, iftImage *area)
{
  iftImage   *pathval = NULL;
  iftGQueue  *Q = NULL;
  int            i, p, q, tmp;
  iftVoxel       u, v;
  iftAdjRel  *A      = iftCircular(1.5);
  float max_dist     = MaxArcWeight(mimg,A);
  int   maxarcweight = 65535;
  
  // Initialization
  
  pathval    = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
  Q          = iftCreateGQueue(maxarcweight+1, mimg->n, pathval->val);

  for (p = 0; p < mimg->n; p++)
  {
    pathval->val[p] = IFT_INFINITY_INT;
    if ((label->val[p]>0)&&(area->val[p]>minarea_thres))
      pathval->val[p] = 0;
    iftInsertGQueue(&Q, p);
  }

  // Image Foresting Transform

  while (!iftEmptyGQueue(Q))
  {
    p = iftRemoveGQueue(Q);
    u = iftMGetVoxelCoord(mimg, p);

    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftMValidVoxel(mimg, v))
      {
        q = iftGetVoxelIndex(mimg, v);

	if (Q->L.elem[q].color != IFT_BLACK) {

	  int W = (int)(maxarcweight*Distance(mimg,p,q)/max_dist);
	  
          tmp = iftMax(pathval->val[p], W);

          if (tmp < pathval->val[q])  {
	    if (Q->L.elem[q].color == IFT_GRAY)
	      iftRemoveGQueueElem(Q,q);
            label->val[q]    = label->val[p];
            pathval->val[q]  = tmp;
            iftInsertGQueue(&Q, q);
          }
        }
      }
    }
  }

  iftDestroyAdjRel(&A);
  iftDestroyGQueue(&Q);
  iftDestroyImage(&pathval);

}

int main(int argc, char *argv[]) {
    timer *tstart;

    if (argc!=6)
      iftError("Usage: label_comp <...>\n"
	       "[1] input image .mimg \n"
	       "[2] adjacency radius \n"
	       "[3] feature distance threshold in (0,1) \n"
	       "[4] minimum area threshold in (0,1) \n"
	       "[5] output directory\n",
	       "main");

    tstart = iftTic();

    if (strcmp(iftFileExt(argv[1]),".mimg")!=0)
      iftError("Input image %s does not have .mimg extension","label_comp.c",argv[1]);
    
    iftMImage *mimg     = iftReadMImage(argv[1]);
    char      *basename = Basename(argv[1]);    
    char       filename[200];
    iftMakeDir(argv[5]);

    iftAdjRel *A        = iftCircular(atof(argv[2]));
    iftImage  *label    = LabelComp(mimg,A,atof(argv[3]));    

    /* post-processing to eliminate small components */
    
    iftImage  *area     = iftRegionArea(label);
    int minarea_thres   = iftMax(iftRound(atof(argv[4])*label->n),1);
    printf("%d\n", minarea_thres);
    PropagateLabelByFmax(mimg, label, minarea_thres, area);
    iftDestroyImage(&area);

    /* --------------------------------------------- */
    
    iftImage  *clabel   = iftColorizeComp(label);
    sprintf(filename,"%s/%s_label.png",argv[5],basename);
    iftWriteImageByExt(label,filename);
    sprintf(filename,"%s/%s_colorlabel.png",argv[5],basename);
    iftWriteImageByExt(clabel,filename);
    
    iftDestroyImage(&label);
    iftDestroyImage(&clabel);
    iftDestroyAdjRel(&A);
    iftDestroyMImage(&mimg);
    iftFree(basename);
    
    printf("Done ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
    return (0);
}
