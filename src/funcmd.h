#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include <ctype.h>

int *CalNumTabPerLine(int *VecTab, char *line){
  
  /*
   * CalNumTabPerLine
   * Returns the positions of tabulator characters in a line.
   * Arguments:
   *  VecTab: integer array that will contain the tabulated positions.
   *  line: string to be analyzed.
   */

	int i,t=0;

	for(i=0; i<strlen(line); i++){
		if(line[i]=='\t'){
			VecTab[t]=i;
			t++;
		}
	}

	return VecTab;
}

char *ExtractColFromLine(char *Value, char *line, int *VecTab, int IndCol){
  
  /*
   * ExtractColFromLine
   * Extract an specific field from a line.
   * Arguments:
   *  Value: string that with hold the desired column.
   *  line: string to be analyzed.
   *  VectTab: integer array with the positions of the tabulators in the line.
   *  IndCol: index of the column to be extracted.
   */

	int LenVal,i;

	LenVal=VecTab[IndCol]-VecTab[IndCol-1];

	for(i=0; i<LenVal; i++){
		Value[i]=line[VecTab[IndCol-1]+i];
	}

  Value[LenVal]='\0';

	return Value; 
}

char *ModGeneIDCol(char *GeneIDMod, char *GeneID, int GeneIDSta){
  
  /*
   * ExtractColFromLine
   * Adapt a gene Id string a separate the number.
   * Arguments:
   *  GeneIdMod: string that with hold the modified Id.
   *  GeneId: string to be adapted.
   *  GeneIDSta: Stable Id label.
   */

	int i,i1=0,LenGeneID,FinLen;
		
	LenGeneID=strlen(GeneID);
	for(i=GeneIDSta; i<LenGeneID; i++){
		GeneIDMod[i1]=GeneID[i];
		i1++;
	}
	
	FinLen=LenGeneID-GeneIDSta;
	GeneIDMod[FinLen]='\0';

	return GeneIDMod;
}

char *ModChrString(char *ChrMod, char *Chr, int Thr){

  /*
   * ModChrString
   * Adjusts a compound chromosome name.
   * Arguments:
   *  ChrMod: new chromosome name.
   *  Chr: chromosome name extracted from the annotation file.
   *  Thr: maximum chrosome name size.
   */
  
	int i,SepLoc;

	if(strlen(Chr)>Thr){
		for(i=0; i<strlen(Chr); i++){
			if(Chr[i]=='|'){
				SepLoc=i;
				break;
			}
		}
			
		for(i=0; i<SepLoc; i++){
			ChrMod[i]=Chr[i];
		}
		ChrMod[SepLoc]='\0';
	}

	else{
		for(i=0; i<strlen(Chr); i++){
			ChrMod[i]=Chr[i];
		}
		ChrMod[strlen(Chr)]='\0';
	}

	return ChrMod;
}

int ConvChrXYMToNum(int ValIn, int nAutosomes){

	int ValOut;

	switch(ValIn){

		case 'X':
		ValOut=nAutosomes+1;
		break;

		case 'Y':
		ValOut=nAutosomes+2;
		break;
	
		case 'M':
		ValOut=nAutosomes+3;
		break;
	}

	return ValOut;
}

/*********************/

// Histogram features.
// Based on https://www.gnu.org/software/gsl/doc/html/histogram.html

typedef struct {
  size_t n ;
  double * range ;
  double * bin ;
} gsl_histogram ;



static int find (const size_t n, const double range[], const double x, size_t * i)
{
      size_t i_linear, lower, upper, mid;
      
      if (x < range[0])
      {
        return -1;
      }
      
      if (x >= range[n])
      {
        return +1;
      }
      
      /* optimize for linear case */
      
    #ifdef LINEAR_OPT
    {
      double u =  (x - range[0]) / (range[n] - range[0]);
      i_linear = (size_t) (u * n);
    }
    
    if (x >= range[i_linear] && x < range[i_linear + 1])
    {
      *i = i_linear;
      return 0;
    }
    #endif
  
  /* perform binary search */
  
  upper = n ;
  lower = 0 ;
  
  while (upper - lower > 1)
  {
    mid = (upper + lower) / 2 ;
    
    if (x >= range[mid])
    {
      lower = mid ;
    }
    else
    {
      upper = mid ;
    }
  }
  
  *i = lower ;
  
  /* sanity check the result */
  
  if (x < range[lower] || x >= range[lower + 1])
  {
    Rprintf ("x not found in range");
    exit(-1);
  }
  
  return 0;
}

int gsl_histogram_accumulate (gsl_histogram * h, double x, double weight)
{
  const size_t n = h->n;
  size_t index = 0;

  int status = find (h->n, h->range, x, &index);

  if (status)
    {
    
      if ( x <= h->range[0] )
        h->bin[0] += weight;
      else
        h->bin[n] += weight;
      //return -1;
    }

  if (index >= n)
    {
      Rprintf("index lies outside valid range of 0 .. n - 1" );
      exit(-1);
    }

  h->bin[index] += weight;
  


  return 0;
}

int gsl_histogram_increment (gsl_histogram * h, double x)
{
  int status = gsl_histogram_accumulate (h, x, 1.0);
  return status;
}

double gsl_histogram_get (const gsl_histogram * h, size_t i)
{
  const size_t n = h->n;
  
  if (i >= n)
  {
    Rprintf ("index lies outside valid range of 0 .. n - 1");
    exit(-1);
  }
  
  return h->bin[i];
}

gsl_histogram * gsl_histogram_alloc (size_t n)
{
  gsl_histogram *h;

  if (n == 0)
    {
      Rprintf ("histogram length n must be positive integer");
      exit(-1);
    }

  h = (gsl_histogram *) malloc (sizeof (gsl_histogram));

  if (h == 0)
    {
      Rprintf ("failed to allocate space for histogram struct");
    exit(-1);
    }

  h->range = (double *) malloc ((n + 1) * sizeof (double));

  if (h->range == 0)
    {
      free (h);         /* exception in constructor, avoid memory leak */

      Rprintf ("failed to allocate space for histogram ranges");
      exit(-1);
    }

  h->bin = (double *) malloc (n * sizeof (double));

  if (h->bin == 0)
    {
      free (h->range);
      free (h);         /* exception in constructor, avoid memory leak */

      Rprintf ("failed to allocate space for histogram bins");
      exit(-1);
    }

  h->n = n;

  return h;
}

static void make_uniform (double range[], size_t n, double xmin, double xmax)
{
  size_t i;
  
  for (i = 0; i <= n; i++)
    {
      double f1 = ((double) (n-i) / (double) n);
      double f2 = ((double) i / (double) n);
      range[i] = f1 * xmin +  f2 * xmax;
    }
}

double * gsl_histogram_set_ranges_uniform (gsl_histogram * h, double xmin, double xmax)
{
  size_t i;
  const size_t n = h->n;

  if (xmin >= xmax)
    {
      Rprintf ("xmin must be less than xmax");
    exit(-1);
    }

  /* initialize ranges */

  make_uniform (h->range, n, xmin, xmax);

  /* clear contents */

  for (i = 0; i < n; i++)
    {
      h->bin[i] = 0;
    }

  return h->range;
}

void gsl_histogram_free (gsl_histogram * h)
{
  //RETURN_IF_NULL (h);
  free (h->range);
  free (h->bin);
  free (h);
}

void gsl_histogram_reset (gsl_histogram * h)
{
  size_t i;
  const size_t n = h->n;

  for (i = 0; i < n; i++)
    {
      h->bin[i] = 0;
    }
}


