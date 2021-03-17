#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include <ctype.h>
#include <omp.h> 

char * reverse_seq(char* seq, int seq_len){
  int i;
  char rev_seq[seq_len+1];
  for (i=0; i < seq_len; i++){
    if (tolower(seq[i])=='c') rev_seq[i] = 'g';
    else if (tolower(seq[i])=='g') rev_seq[i] = 'c';
    else if (tolower(seq[i])=='t') rev_seq[i] = 'a';
    else if (tolower(seq[i])=='a') rev_seq[i] = 't';
  }
  rev_seq[seq_len]='\0';
  return(strdup(rev_seq));
}

short contains_wrong_chars(char * seq, int seq_len){
    char c_char;
    int i;
    for (i=0; i < seq_len; i++){
        c_char = tolower(seq[i]);
        if (c_char!='c' && c_char!='g' && c_char!='t' && c_char!='a')
            return(1);
    }
    return(0);
}

/* 
 * General comments:
 * nProt indicates the number of instances being protected from the garbage collector in each function,
 * so as to unprotect them at the end of the execution.
 * 
 */

SEXP SeqDic(SEXP CooMet, SEXP LenDic, SEXP nrow, SEXP SeqChr, SEXP LenChr, SEXP TargetGrowMode, SEXP X, SEXP sense){

  
  /*
   * SeqDic
   * Extracts words from the fasta sequence of the chromosome and returns them
   * along with the methylation frequency of its corresponding target site.
   * CooMet: Vector with coordinates and methylation frequencys of the target sites.
   * LenDic: Length of the words to be extracted.
   * nrow: Length of the coordinates' vector.
   * SeqChr: Sequence of the current chromosome.
   * TargetGrowMode: Growing mode of the words around the start coordinate.
   *  "C" (central), "R" (Lateral right) or "L" (Lateral left).
   * X: Displacement of the words from the start coordinate.
   * 
   */
  
	int w,begin,nrowCooMet,nProt=0,Len,nElem, Xcg, x, MaxLen, sns;
	int TargetLocat;
	double *ColCoo;
	SEXP MatSeq;
	SEXP DicMet;
	SEXP ColMet;
	
	CooMet = PROTECT(coerceVector(CooMet, VECSXP)); nProt++;
	// Separate target coordinates and methylation frequencies.
	ColCoo = REAL(VECTOR_ELT(CooMet,0)); 	
	ColMet = VECTOR_ELT(CooMet,1); 

	LenDic = PROTECT(coerceVector(LenDic, REALSXP)); nProt++;
	w = REAL(LenDic)[0];
	Len = 2*w + 2;

	nElem=Len+1;
	// Auxiliary structure for saving sequences.
	char SubSeq[nElem];

	nrow = PROTECT(coerceVector(nrow, REALSXP)); nProt++;
	nrowCooMet = REAL(nrow)[0];
	// Structure for saving the real beginning coordinate of each word.
	double *Sta = (double*)malloc(nrowCooMet*sizeof(double));

	SeqChr = PROTECT(coerceVector(SeqChr,STRSXP)); nProt++;
	SEXP c = STRING_ELT(SeqChr, 0);
	const char *v = CHAR(c);

	TargetGrowMode = PROTECT(coerceVector(TargetGrowMode,STRSXP)); nProt++;
	SEXP TargetLoc = STRING_ELT(TargetGrowMode, 0);
	const char *Target = CHAR(TargetLoc);
	char TargetMode = Target[0];

	X = PROTECT(coerceVector(X, REALSXP)); nProt++;
	x = REAL(X)[0];
	sense = PROTECT(coerceVector(sense, REALSXP)); nProt++;
	sns = REAL(sense)[0];
	
	LenChr = PROTECT(coerceVector(LenChr, REALSXP)); nProt++;
	MaxLen = REAL(LenChr)[0];

	// Save space for a R structure-like string table.
	MatSeq = PROTECT(allocVector(STRSXP, nrowCooMet)); nProt++;

	// Structured to be returned, will contain the word and the methylation frequency.
	DicMet = PROTECT(allocVector(VECSXP, 2)); nProt++;

	// Extract word from the sequence depending on the growing mode.
	switch( TargetMode)
        {
        	case 'C': // CENTRAL

    // Get a word from each target coordinate.
		for(int i=0; i<nrowCooMet; i++){ 

			Xcg = ColCoo[i];
		  
		  // Check limits.
			Sta[i] = Xcg + x - w - 1;
			if (((Sta[i])<0) || ((Sta[i]+Len)>MaxLen)) {
			  
			  Sta[i]=-1;
			  strcpy(SubSeq, "no");
			  SubSeq[2] = '\0';
			  
			}
			else {
			  
  			begin = Sta[i];
  			strncpy(SubSeq, v+begin, Len);
  			SubSeq[Len] = '\0';
  			if (contains_wrong_chars(SubSeq, Len)){
                            Sta[i]=-1;
                            strcpy(SubSeq, "no");
                            SubSeq[2] = '\0';
                        }
                        else if (sns==0){
  			  strcpy(SubSeq, reverse_seq(SubSeq, Len));
  			}
			}
		
		  // Put the word in its position in the vector, according to the 
		  // position of its target coordinate.
			SET_STRING_ELT(MatSeq, i, mkChar(SubSeq));	
		}
		break;

		case 'R': // LATERAL RIGHT
	
		for(int i=0; i<nrowCooMet; i++){ 
	
			Xcg = ColCoo[i];
			begin = Xcg + x - 1;

			if ((begin<0) || ((begin+Len)>MaxLen)) {
			  
			  //Sta[i]=-1;
			  strcpy(SubSeq, "no");
			  SubSeq[2] = '\0';
			  
			}
			else {
  			strncpy(SubSeq, v+begin, Len);
  			SubSeq[Len] = '\0';
  			 if (contains_wrong_chars(SubSeq, Len)){
                            Sta[i]=-1;
                            strcpy(SubSeq, "no");
                            SubSeq[2] = '\0';
                        }
                        else if (sns==0){
  			  strcpy(SubSeq, reverse_seq(SubSeq, Len));
  			}
  		
  				
			}
			SET_STRING_ELT(MatSeq, i, mkChar(SubSeq));
		}
		break;

		case 'L': // LATERAL LEFT
	
		for(int i=0; i<nrowCooMet; i++){ 
	
			Xcg = ColCoo[i];
			begin = Xcg + x - 2*w - 1;

			if ((begin<0) || ((begin+Len)>MaxLen)) {
			  
			  //Sta[i]=-1;
			  strcpy(SubSeq, "no");
			  SubSeq[2] = '\0';
			  
			  
			}
			else{
  			strncpy(SubSeq, v+begin, Len);
  			SubSeq[Len] = '\0';
  			 if (contains_wrong_chars(SubSeq, Len)){
                            Sta[i]=-1;
                            strcpy(SubSeq, "no");
                            SubSeq[2] = '\0';
                        }
                        else if (sns==0){
  			  strcpy(SubSeq, reverse_seq(SubSeq, Len));
  			}
  			
  		
  			
			}
			SET_STRING_ELT(MatSeq, i, mkChar(SubSeq));
		}
		break;

	}

	// First column of the return, words.
	// Second column of the return, methylation rates.
	SET_VECTOR_ELT(DicMet,0,MatSeq);
	SET_VECTOR_ELT(DicMet,1,ColMet);

	free(Sta);

	UNPROTECT(nProt);

	return DicMet;
}

// TODO: documentar argumentos (w; longitud)
SEXP DissimilarityMatrix(SEXP PomMat, SEXP NumRow, SEXP w, SEXP Metric){

	int nrowPomMat,nVe,l=0,NumInfDiagEle,nRow,nProt=0;
	double *xRans;
	double *VecElem;
	double PeaCor, NorEleMinMea, NorEleMod;
	double EleMinMea, mean, EleMod; //Vector Element Minus Mean 
	SEXP ans;

	PomMat = PROTECT(coerceVector(PomMat, VECSXP)); nProt++;
	nrowPomMat = length(PomMat);

	NumRow = PROTECT(coerceVector(NumRow, REALSXP)); nProt++;
	nRow = REAL(NumRow)[0];

	NumInfDiagEle=(pow(nRow,2)-nRow)/2;

	w = PROTECT(coerceVector(w, REALSXP)); nProt++;
	nVe = REAL(w)[0];

	Metric = PROTECT(coerceVector(Metric,STRSXP)); nProt++;
	SEXP MetricID = STRING_ELT(Metric, 0);
	const char *TypeMetric = CHAR(MetricID);
	char MetricMode = TypeMetric[0];

	ans=PROTECT(allocVector(REALSXP,NumInfDiagEle)); nProt++;
	xRans = REAL(ans);

        //Memory allocation for PomMatStd

        double **PomMatStd = (double**)malloc(nrowPomMat*sizeof(double));
        for(int i =0; i<nrowPomMat; i++){

	         PomMatStd[i]=(double*)malloc(nVe*sizeof(double));
        }
    
        switch(MetricMode)
        {
        	case 'P': // Pearson correlation
 
		// Standarization and normalization of the vectors
		for(int i1=0; i1<nrowPomMat; i1++){
		
			VecElem = REAL(VECTOR_ELT(PomMat, i1));
		
			mean=0.0;
			EleMinMea=0.0;
			
			// Calculate means
			for(int k=0; k<nVe ; k++){
		 		mean += VecElem[k];
		 	}
		 	mean = mean/nVe;
		 
		 	// Calculate modulus
		 	for(int k=0; k<nVe; k++){

		 		VecElem[k]=VecElem[k]-mean;
		 		EleMinMea += pow(VecElem[k],2);
		 	}

		 	// Calculate the norm
		 	NorEleMinMea=sqrt(EleMinMea);
		 
		 	// Normalization
		 	for(int k=0; k<nVe; k++){

			     PomMatStd[i1][k]=VecElem[k]/NorEleMinMea;
	   		} 
		}
		break;

          	case 'C': // Cosine metric

		// Standarization and normalization of the vectors
		for(int i1=0; i1<nrowPomMat; i1++){
		
			VecElem = REAL(VECTOR_ELT(PomMat, i1));
		
			EleMod=0.0;
		 
		 	// Calculate modulus
		 	for(int k=0; k<nVe; k++){

		 		EleMod += pow(VecElem[k],2);
		 	}

		 	// Calculate the norm
		 	NorEleMod=sqrt(EleMod);
		 
		 	// Normalization
		 	for(int k=0; k<nVe; k++){

			     PomMatStd[i1][k]=VecElem[k]/NorEleMod;
	   		} 
		}
		break;

        } // switch(MetricID)
        
        // Calculation of the matrix of distance between the standarized vectors
        
	for(int i1=0; i1<nrowPomMat-1; i1++){

		for(int i2=i1+1; i2<nrowPomMat; i2++){
			
			PeaCor=0.0;
                        // Calculate the scalar product
		        for(int k=0; k<nVe; k++){

		       		PeaCor += PomMatStd[i1][k]*PomMatStd[i2][k];
			}

			xRans[l]=PeaCor;
                 	l += 1;
		}	
	}
	free(PomMatStd);

	UNPROTECT(nProt);
	return ans;
}

SEXP readCooChrFile(SEXP SelRow, SEXP LenArrayRow, SEXP DicLen, SEXP TempDir){

  /*
   * SelRow: 
   * LenArrayRow:
   * DicLen:
   * TempDir:
   * 
   */
  
  //FILE * flogs;
  //flogs = fopen("readCooChrFile15-07-2019-04.logs", "a");
  
	int i=0,i3=0,t,k,Len,LenLine,returnval;
	long int r;
	int nProt=0,w,MaxRow;
	long int PosIni,PosEnd,nElem;
	double *output;
	double *RowArray;
	char *line;
	SEXP FileLine, Out;
	
	SelRow = PROTECT(coerceVector(SelRow,REALSXP)); nProt++;
	RowArray=REAL(SelRow);

	LenArrayRow = PROTECT(coerceVector(LenArrayRow,REALSXP)); nProt++;
	Len = REAL(LenArrayRow)[0];

	//fprintf(flogs, "Protections 1\n");
           
	long int VecIni[Len*2];
	long int VecEnd[Len*2];
	long int VecnElem[Len*2];

	DicLen = PROTECT(coerceVector(DicLen, REALSXP)); nProt++;
	w = REAL(DicLen)[0];

	TempDir = PROTECT(coerceVector(TempDir,STRSXP)); nProt++;
	SEXP temporaldir = STRING_ELT(TempDir, 0);
	const char *tmpdir = CHAR(temporaldir);

	FileLine = PROTECT(allocVector(VECSXP, Len)); nProt++;
	
	//fprintf(flogs, "Protections 2\n");
	
	FILE *pF;
	FILE *pFCol;

	char fileData[1000];
	char fileCol[1000];

	sprintf(fileData,"%s/fileData_%d",tmpdir,w);
	sprintf(fileCol,"%s/fileCol_%d",tmpdir,w);

	pF=fopen(fileData,"r");
	pFCol=fopen(fileCol,"r");
	
	//fprintf(flogs, "Files opened\n");

	/*if(pF == NULL)
	{
	  fprintf(flogs, "Name: %s\tDoes not exist.\n",
           fileData);
	}
	if(pF == NULL)
	{
	  fprintf(flogs, "Name: %s\tDoes not exist.\n",
           fileCol);
	}*/
	
	int j=0;
	int ini=1;
	int fin;
	long line_count=0;
	//fprintf(flogs, "Length for the first for: %d\n", Len);
	for(t=0; t<Len; t++){

		fin=RowArray[t];
	  //fprintf(flogs, "\tFor 1.1 : t = %d\n, fin = %d", t, fin);

		for(r=ini; r<=fin; r++){

		  //fprintf(flogs, "\t\tFor 1.2 : r = %ld, line read= %ld, j=%d\n", r, line_count, j);
			returnval=fscanf(pFCol,"%ld %ld %ld",&PosIni,&PosEnd,&nElem);
			line_count += 1;
			
			
			if(r==fin){

				VecIni[j]=PosIni;
	         		VecEnd[j]=PosEnd;
				VecnElem[j]=nElem;
        	         	j=j+1;
				ini=fin+1;
			}
			//fprintf(flogs, "\t\t\tj=%d, ini=%d\n", j, ini);
		}
	}
	//fprintf(flogs, "First for done\n");
								
	for(j=0; j<Len; j++){

	  //fprintf(flogs, "\tFor 2.1 : j = %d\n", j);
		fseek(pF,VecIni[j],SEEK_SET);
		int n=VecEnd[j]-VecIni[j];
		line = (char*) malloc (sizeof(char)*(n+1)); 
		returnval=fread(line,1,n,pF);
		// Avoid memory corruption issues.
		line[n]='\0';
		
		//fprintf(flogs, "\t\tfseek for %d done\n", j);
		
		char seps[]=" ";
		char *token;
		//double var;
		int var;
		int i=0;
		int nProtLoop=0;

		token=strtok(line,seps); //string separator
		//fprintf(flogs, "\t\tstrtok  done\n");

		LenLine=VecnElem[j];
		int input[LenLine];

		Out = PROTECT(allocVector(REALSXP,LenLine)); nProtLoop++;
		output = REAL(Out);

		while(token!=NULL){
			var=atoi(token);
			input[i++]=var;
			token=strtok(NULL,seps);
		}
    //fprintf(flogs, "\t\twhile  done\n");
		
		for(i=0; i<LenLine; i++){
			output[i]=(double)input[i];
		}

		SET_VECTOR_ELT(FileLine,j,duplicate(Out));
		UNPROTECT(nProtLoop);
		free(line);
	}
	
	//fprintf(flogs, "Second for done\n");

	fclose(pF);
	fclose(pFCol);
	//fclose(flogs);
	
	UNPROTECT(nProt);
	return FileLine;
}

SEXP CooChr(SEXP SeqMetCooChr, SEXP CumSum, SEXP DicLen, SEXP TempDir){

	/*
	* Creates coordinate files for the words related to each target.
	* Created files (being X the length of the words referenced):
	*
	* 	fileData_X : Contains a line for each word of length X.
	*	Each line contains the appearances of the word in the
	*	genome. The length of the lines is variable.
	*	<ChromosomeN> <PositionN> ... <ChromosomeM> <PositionM>
	*	
	*	fileCol_X : Auxiliary file for easing the read of 
	*	fileData_X in R, due to the variable length of its lines.
	*	Each line contains, for its corespondent line in
	*	fileData_X, the start and end characters' indexes and
	*	the number of fields of the lines.
	*	<first_character> <last_character> <number of fields>
	*
	*/

	int nProt=0,i2,i3,fin,nElem,valchr,valcoo;
	double *ColCoo,*ColChr,*CumVec;
	long pos;

	SeqMetCooChr = PROTECT(coerceVector(SeqMetCooChr, VECSXP)); nProt++;
	ColCoo = REAL(VECTOR_ELT(SeqMetCooChr,2)); 
	ColChr = REAL(VECTOR_ELT(SeqMetCooChr,3)); 
	
	CumSum = PROTECT(coerceVector(CumSum,REALSXP)); nProt++;
	CumVec = REAL(CumSum);
	int nCum=length(CumSum);
	
	DicLen = PROTECT(coerceVector(DicLen, REALSXP)); nProt++;
	int w = REAL(DicLen)[0];

	TempDir = PROTECT(coerceVector(TempDir,STRSXP)); nProt++;
	SEXP temporaldir = STRING_ELT(TempDir, 0);
	const char *tmpdir = CHAR(temporaldir);

	FILE *pF;
	FILE *pFCol;

	char fileData[1000];
	char fileCol[1000];

	sprintf(fileData,"%s/fileData_%d",tmpdir,w);
	sprintf(fileCol,"%s/fileCol_%d",tmpdir,w);

	pF=fopen(fileData,"w");
	pFCol=fopen(fileCol,"w");
	
	for(int i1=0; i1<nCum; i1++){
	
		if(i1==0){

			fin=CumVec[i1];
			i2=0; 
			i3=0;
			int chro[fin];
			int coor[fin];
			
			for(i2; i2<fin; i2++){
				if(i2==0){
					pos=ftell(pF);
					fprintf(pFCol,"%ld",pos);
				}
				valchr=ColChr[i2]; i3++;
				valcoo=ColCoo[i2]; i3++;
				fprintf(pF,"%d ",valchr);
				fprintf(pF,"%d ",valcoo);
			}
			fprintf(pF,"\n");
			pos=ftell(pF);
			fprintf(pFCol," %ld",pos);
			fprintf(pFCol," %d\n",i3);
		}
		else{
			fin=CumVec[i1];
			i2=CumVec[i1-1];
			i3=0;
			nElem=fin-i2;
			int chro[nElem];
			int coor[nElem];
			
			for(i2; i2<fin; i2++){
				if(i2==CumVec[i1-1]){
					pos=ftell(pF);
					fprintf(pFCol,"%ld",pos);
				}
				valchr=ColChr[i2]; i3++;
				valcoo=ColCoo[i2]; i3++;
				fprintf(pF,"%d ",valchr);
				fprintf(pF,"%d ",valcoo);
			}
			fprintf(pF,"\n");
			pos=ftell(pF);
			fprintf(pFCol," %ld",pos);
			fprintf(pFCol," %d\n",i3);
		}				
	}

	fclose(pF);
	fclose(pFCol);
			
	UNPROTECT(nProt);

	return R_NilValue;
}

#include "funcmd.h"

SEXP filtmdfile(SEXP NumAutosomes, SEXP dirmd){

  /*
   * filtmdfile
   * Reads an md gene annotation file and selects those that are not 
   * of the type "Un". Creates an R-manageable txt file.
   * Arguments:
   *  NumAutosomes: Total number of autosomes in the species.
   *  dirmd: file where the filtered information will be written.
   *
   */
  
	int nProt=0;
	int nAutosomes;

	NumAutosomes = PROTECT(coerceVector(NumAutosomes, REALSXP)); nProt++;
	nAutosomes = REAL(NumAutosomes)[0];

	dirmd = PROTECT(coerceVector(dirmd,STRSXP)); nProt++;
	SEXP mdfiledir = STRING_ELT(dirmd, 0);
	const char *dir = CHAR(mdfiledir);

	FILE *pFR;
	FILE *pFW;

	char fileRead[1000];
	char fileWrite[1000];

	// Source file.
	sprintf(fileRead,"%s/seq_gene.md",dir);
	// Destination file.
	sprintf(fileWrite,"%s/seq_gene_filt.txt",dir);

	pFR=fopen(fileRead,"r");
	pFW=fopen(fileWrite,"w");

	// Declaration of the columns of interest in the md file.
	int ColChr=1,ColPosIni=2,ColPosEnd=3;
	int ColGeneNam=9,ColGeneID=10,ColFeaTyp=11;
	int GeneIDSta=8,Thr=3,ValIn,ValOut;
	// Possible gene annotation type patterns.
	char *patternGene="GENE",*patternUn="Un",*patternMT="MT";
	char *retval;

	char *line = malloc(sizeof(char) * 500); 
	int *VecTab = malloc(sizeof(int) * 20);
	char *Chr = malloc(sizeof(char) * 100);
	char *ChrMod = malloc(sizeof(char) * 20); 
	char *PosIni = malloc(sizeof(char) * 100); 
	char *PosEnd = malloc(sizeof(char) * 100); 
	char *GeneNam = malloc(sizeof(char) * 100); 
	char *GeneID = malloc(sizeof(char) * 100); 
	char *FeaTyp = malloc(sizeof(char) * 100);
	char *GeneIDMod =malloc(sizeof(char) * 10);

	// Column names' line.
	retval=fgets(line,500,pFR);

	while(fgets(line,500,pFR)){

	  // Tabulated positions in the line.
		VecTab = CalNumTabPerLine(VecTab,line);
		
		// Get the column with the chromosome number.
		Chr=ExtractColFromLine(Chr,line,VecTab,ColChr);
		
		// Adjust the name of the chromosome if it is compound.
		ChrMod=ModChrString(ChrMod,Chr,Thr);
		
		// Adapt the name if the chromosome is mitochondrial.
		if(strstr(ChrMod,patternMT)){
			ChrMod[1]=ChrMod[1];
			ChrMod[2]='\0';
		}

		ValIn=ChrMod[1];

    PosIni=ExtractColFromLine(PosIni,line,VecTab,ColPosIni);
    PosEnd=ExtractColFromLine(PosEnd,line,VecTab,ColPosEnd);
		GeneNam=ExtractColFromLine(GeneNam,line,VecTab,ColGeneNam);

		GeneID=ExtractColFromLine(GeneID,line,VecTab,ColGeneID);
		// Get the number from the Gene Id.
		GeneIDMod=ModGeneIDCol(GeneIDMod,GeneID,GeneIDSta);

		FeaTyp=ExtractColFromLine(FeaTyp,line,VecTab,ColFeaTyp);

		// Extract annotation of type "Gene" and discard those of type "Un".
		if(strstr(FeaTyp,patternGene)!=NULL & strstr(ChrMod,patternUn)==NULL){

		  // 'X', 'Y' and 'M' are coded as numbers.
			if(ValIn=='X' || ValIn=='Y' || ValIn=='M'){
	
				ValOut=ConvChrXYMToNum(ValIn,nAutosomes);
				fprintf(pFW,"\t%d\t%s\t%s\t%s\t%s\n",ValOut,PosIni,PosEnd,GeneNam,GeneIDMod);
			}

			else {
				fprintf(pFW,"%s\t%s\t%s\t%s\t%s\n",ChrMod,PosIni,PosEnd,GeneNam,GeneIDMod);
			}						
		}
	}
	
	// Liberate memory.
	free(line);
	free(VecTab);
	free(Chr);
	free(ChrMod);
	free(PosIni);
	free(PosEnd);
	free(GeneNam);
	free(GeneID);
	free(FeaTyp);	
	free(GeneIDMod);
		
	fclose(pFR);
	fclose(pFW);

	UNPROTECT(nProt);

	return R_NilValue;
}

char Complementarity(char CharIn){

	char CharOut;

	switch(CharIn){

    		case 'a'  :
       		CharOut='t';
       		break; 
    		
		case 'c'  :
      		CharOut='g';
       		break; 

		case 'g' :
		CharOut='c';
		break;

		case 't' :
		CharOut='a';
		break;
		
	  //case 'n' :
	  //CharOut='n';
	  //break;
	}
	return CharOut;
}

char *strrev (char *str){

    char *begin = str;
    char *end = str + strlen (str) - 1;
    char tmp;

    while (end > begin)
    {
        tmp = *end;
        *end-- = *begin;
        *begin++ = tmp;
    }

    return str;
}

SEXP Reverse(SEXP Seq, SEXP LenDic){

	int nProt=0,nSeq,w,Len;
	char CharIn,CharOut;
	SEXP CompVecSeq;
	SEXP c;

	Seq = PROTECT(coerceVector(Seq,STRSXP)); nProt++;
	nSeq = length(Seq);

	LenDic = PROTECT(coerceVector(LenDic, REALSXP)); nProt++;
	w = REAL(LenDic)[0];
	Len = 2*w + 2;

	int nLen=Len+1;
	char CharOutVec[nLen];

	CompVecSeq = PROTECT(allocVector(STRSXP, nSeq)); nProt++;

	for(int i1=0; i1<nSeq; i1++){

		c = STRING_ELT(Seq, i1);
		const char *v = CHAR(c);

		for(int i2=0; i2<Len; i2++){

			CharIn=v[i2];
			CharOut=Complementarity(CharIn);
			CharOutVec[i2]=CharOut;
			CharOutVec[Len] = '\0';
		}

		strrev(CharOutVec);
	 	SET_STRING_ELT(CompVecSeq, i1, mkChar(CharOutVec));
	}

	UNPROTECT(nProt);
	return CompVecSeq;

}

SEXP Scan(SEXP WeiLogPom, SEXP WeiLogPmv, SEXP NumSeq, SEXP w){

  //FILE * flogs;
  //flogs = fopen("Scan15-07-2019-04.logs", "a");
  
  int LenMot,nSeq,IndBas,n=4,nProt=0;
  double Num,Den,s,Sco;
  double *WeiLogPomElem, *WeiLogPmvElem, *NumSeqElem,*xResult;
  SEXP result;
  
  WeiLogPom = PROTECT(coerceVector(WeiLogPom, VECSXP)); nProt++; 
  WeiLogPmv = PROTECT(coerceVector(WeiLogPmv, VECSXP)); nProt++;
  NumSeq = PROTECT(coerceVector(NumSeq, VECSXP)); nProt++;
  w = PROTECT(coerceVector(w, REALSXP)); nProt++;
  
  LenMot = REAL(w)[0];
  nSeq=length(NumSeq);
  
  result = PROTECT(allocVector(REALSXP, nSeq)); nProt++;
  xResult = REAL(result);
  
  WeiLogPomElem = REAL(VECTOR_ELT(WeiLogPom, 0)); 
  WeiLogPmvElem = REAL(VECTOR_ELT(WeiLogPmv, 0));
  
  //fprintf(flogs, "Protections done, nseq=%d\n", nSeq);
  
  for(int i1=0; i1<nSeq; i1++){
    
    
    NumSeqElem = REAL(VECTOR_ELT(NumSeq, i1));
    //fprintf(flogs, "This i1 is i1=%d\n", i1);
    Sco=0.0;
    
    for(int i2=0; i2<LenMot; i2++){
      
      //fprintf(flogs, "NumSeqElem=%f\n", NumSeqElem[i2]);
      IndBas=NumSeqElem[i2];
      Num=WeiLogPomElem[IndBas+i2*n];
      Den=WeiLogPmvElem[i2];
      s=Num-Den;
      Sco=Sco+s;
    }
    //fprintf(flogs, "Score for i1=%d, Sco=%lf\n", i1, Sco);
    
    xResult[i1]=Sco;
  }
  
  //fprintf(flogs, "Scan has ended");
  
  //fclose(flogs);
  UNPROTECT(nProt);
  return result;
}

SEXP scanPOMs (SEXP POMvec, SEXP PMVvec, SEXP NumPOM, SEXP SeqVec, SEXP NumSeq, SEXP Width, SEXP NBins, SEXP NCpu) {
  
  int width, nBins, numSeq, numPOM, nProt, n_cpu;
  int * seqs;
  double * POMs;
  double * PMVs;
  double * breaks;
  double * counts;
  
  nProt = 0;
  
  
  // Initialize data structures
  
  POMvec = PROTECT(coerceVector(POMvec, REALSXP)); nProt++;
  POMs = REAL(POMvec); 
  
  PMVvec = PROTECT(coerceVector(PMVvec, REALSXP)); nProt++;
  PMVs = REAL(PMVvec); 
  
  SeqVec = PROTECT(coerceVector(SeqVec, INTSXP)); nProt++;
  seqs = INTEGER(SeqVec); 
  
  NumPOM = PROTECT(coerceVector(NumPOM, REALSXP)); nProt++;
  numPOM = REAL(NumPOM)[0];
  
  NumSeq = PROTECT(coerceVector(NumSeq, REALSXP)); nProt++;
  numSeq = REAL(NumSeq)[0];
  
  Width = PROTECT(coerceVector(Width, REALSXP)); nProt++;
  width = REAL(Width)[0];
  
  NBins = PROTECT(coerceVector(NBins, REALSXP)); nProt++;
  nBins = REAL(NBins)[0];
  
  NCpu = PROTECT(coerceVector(NCpu, REALSXP)); nProt++;
  n_cpu = REAL(NCpu)[0];
  
  SEXP BreaksCounts = PROTECT(allocVector(REALSXP, ( nBins * numPOM + ( nBins) * numPOM ))); nProt++;
  double * break_counts;
  break_counts = REAL(BreaksCounts);
  
  for ( int j=0; j<( nBins * numPOM + ( nBins) * numPOM ); j++)
    break_counts[j] =1;
  Rprintf("allocated %d break counts\n",  ( nBins * numPOM + ( nBins) * numPOM ) );
  
  int seq_c;
  int pom_c;
  int nt_c;
  int base_c;
  int nucl_counter;
  double score;
  double num;
  double den;
  int i;
  int disp_p_pom;
  int disp_this_pom;
  int disp_this_seq;
  disp_p_pom = width * 4;
  int disp_this_pmv;
  double min;
  double max;
  double h;
  
  int disp_this_pom_breaks;
  int disp_this_pom_counts;
  
  double * current_breaks;
  
  double * scores;
  scores = (double *) malloc(sizeof(double)*numSeq);
  
  gsl_histogram * my_hist;
  
  my_hist = gsl_histogram_alloc( nBins );
  
  Rprintf("Scanning in parallel with %d cpus\n",n_cpu);
  
  for ( pom_c = 0; pom_c < numPOM; pom_c++ ) {
    //for ( pom_c = 0; pom_c < 100; pom_c++ ) {
    
    disp_this_pom = disp_p_pom * pom_c;
    disp_this_pmv = width * pom_c;
    //Rprintf("(pom %d) disp_this_pom %d (%d), disp_this_pmv %d (%d)\n", pom_c, disp_this_pom, numPOM*width*4, disp_this_pmv, numPOM*width);
    
    //Rprintf("Analyzing pom %d\n", pom_c);
    
    /*Rprintf("Scanning POM ");
     for ( int idx=0; idx < width*4; idx += 4 ) {
     Rprintf("%f", POMs[ disp_this_pom + 4 * idx ]);
     Rprintf("-");
     Rprintf("%f", POMs[ disp_this_pom + 4 * idx + 1]);
     Rprintf("-");
     Rprintf("%f", POMs[ disp_this_pom + 4 * idx + 2 ]);
     Rprintf("-");
     Rprintf("%f", POMs[ disp_this_pom + 4 * idx + 3 ]);
     Rprintf(" ");
     
     }
     Rprintf("\n");*/
    
    
      score = 0;
     
     for ( nt_c = 0; nt_c < width; nt_c ++ ){
     num = POMs[ disp_this_pom + 4 * nt_c + seqs[ nt_c ] ];
     den = PMVs[ disp_this_pmv + nt_c ];
     
     score += ( num - den );
     }
     
     //Rprintf("%f ", score);
     scores[0] = score;
     min=max=score;
     
     // Get scores.
     for ( seq_c = 1; seq_c < numSeq; seq_c++ ){
     score = 0;
     disp_this_seq = width * seq_c;
     
     for ( nt_c = 0; nt_c < width; nt_c ++ ){
     num = POMs[ disp_this_pom + 4 * nt_c + seqs[ disp_this_seq + nt_c ] ];
     den = PMVs[ disp_this_pmv + nt_c ];
     
     score += ( num - den );
     }
     
     if( min > score )
     min = score;   
     if( max < score )
     max = score;
     
     //Rprintf("%f ", score);
     scores[seq_c] = score;
     //Rprintf("%f ", scores[seq_c]);
     }
    //Rprintf("Added %d scores\n", numSeq);
    
    /*min=max=my_scores[0];
     
     for(i=1; i<numSeq; i++)
     {
     
     }*/
    //Rprintf("POM %d\n", pom_c); 
    
    
    if ( min >= max ) {
     min = max - 1;
     max = max + 1;
    }
    
    current_breaks = gsl_histogram_set_ranges_uniform(my_hist, min, max);
    
    for ( seq_c = 0; seq_c < numSeq; seq_c++ )
      gsl_histogram_increment( my_hist, scores[seq_c] );
    
    //current_breaks = gsl_histogram_set_ranges_uniform(my_hist, 0, numSeq);
    /*for ( seq_c = 0; seq_c < numSeq; seq_c++ )
     gsl_histogram_increment( my_hist, seq_c );*/
    
    
    disp_this_pom_breaks = nBins * numPOM + pom_c * ( nBins );
    disp_this_pom_counts = pom_c * nBins;
    
    
    //Rprintf("counts [%d, %d] and breaks [%d, %d]\n", 
    //         disp_this_pom_counts, disp_this_pom_counts + nBins -1, 
    //        disp_this_pom_breaks, disp_this_pom_breaks + nBins);
    //for ( i = 0; i < nBins; i++ ){
    //  Rprintf("%d ",  gsl_histogram_get(my_hist, (size_t) i));
    //}
    //Rprintf("\n");
    
    //Rprintf("(%d pom) writing %d - %d and %d - %d (max %d)\n", pom_c, disp_this_pom_counts, nBins+disp_this_pom_counts, disp_this_pom_breaks, nBins+disp_this_pom_breaks, nBins * numPOM + ( nBins) * numPOM );
    for ( i = 0; i < nBins; i++ ) {
     break_counts [ i + disp_this_pom_counts ] = gsl_histogram_get(my_hist, (size_t) i);
     break_counts [ i + disp_this_pom_breaks ] = current_breaks[i];
    }
    
    
    gsl_histogram_reset( my_hist );
  }
  
  UNPROTECT(nProt);
  
  return BreaksCounts;
}

SEXP scanPOMs_par (SEXP POMvec, SEXP PMVvec, SEXP NumPOM, SEXP SeqVec, SEXP NumSeq, SEXP Width, SEXP NBins, SEXP NCpu ){

  
  int width, nBins, numSeq, numPOM, nProt, n_cpu;
  int * seqs;
  double * POMs;
  double * PMVs;
  double * breaks;
  double * counts;
  
  nProt = 0;

  
  // Initialize data structures
  
  POMvec = PROTECT(coerceVector(POMvec, REALSXP)); nProt++;
  POMs = REAL(POMvec); 
  
  PMVvec = PROTECT(coerceVector(PMVvec, REALSXP)); nProt++;
  PMVs = REAL(PMVvec); 
  
  SeqVec = PROTECT(coerceVector(SeqVec, INTSXP)); nProt++;
  seqs = INTEGER(SeqVec); 
  
  NumPOM = PROTECT(coerceVector(NumPOM, REALSXP)); nProt++;
  numPOM = REAL(NumPOM)[0];
  
  NumSeq = PROTECT(coerceVector(NumSeq, REALSXP)); nProt++;
  numSeq = REAL(NumSeq)[0];
  
  Width = PROTECT(coerceVector(Width, REALSXP)); nProt++;
  width = REAL(Width)[0];
  
  NBins = PROTECT(coerceVector(NBins, REALSXP)); nProt++;
  nBins = REAL(NBins)[0];
  
  NCpu = PROTECT(coerceVector(NCpu, REALSXP)); nProt++;
  n_cpu = REAL(NCpu)[0];
  
  SEXP BreaksCounts = PROTECT(allocVector(REALSXP, ( nBins * numPOM + ( nBins) * numPOM ))); nProt++;
  double * break_counts;
  break_counts = REAL(BreaksCounts);
  
  for ( int j=0; j<( nBins * numPOM + ( nBins) * numPOM ); j++)
    break_counts[j] =1;
  Rprintf("allocated %d break counts\n",  ( nBins * numPOM + ( nBins) * numPOM ) );
  
  // Initialization is OK
  
 //omp_set_num_threads(n_cpu);
  
  
  int seq_c;
  int pom_c;
  int nt_c;
  int base_c;
  int nucl_counter;
  double score;
  double num;
  double den;
  int i;
  int disp_p_pom;
  int disp_this_pom;
  int disp_this_seq;
  disp_p_pom = width * 4;
  int disp_this_pmv;
  double min;
  double max;
  double h;
  
  int disp_this_pom_breaks;
  int disp_this_pom_counts;
  
  double * current_breaks;
  
  double * scores;
  double * my_scores;
  scores = (double *) malloc(sizeof(double)*numSeq*n_cpu);

  gsl_histogram * hists[n_cpu];
  gsl_histogram * my_hist;
  for ( i = 0; i < n_cpu; i++ )
    hists[i] = gsl_histogram_alloc( nBins );
  
  //my_hist = gsl_histogram_alloc( nBins );
  
  //omp_set_dynamic(0);         // Explicitly disable dynamic teams
  //omp_set_num_threads(n_cpu);
  
  Rprintf("Scanning in parallel with %d cpus\n",n_cpu);
  
  
  // #pragma omp parallel for private(my_hist, disp_this_pom, disp_this_pmv, disp_this_seq, my_scores, score, num, den, seq_c, nt_c, max, min, current_breaks, disp_this_pom_breaks, disp_this_pom_counts, i, pom_c) shared( width, POMs, PMVs, nBins, numPOM, numSeq, break_counts, scores, hists)
   for ( pom_c = 0; pom_c < numPOM; pom_c++ ) {
  //for ( pom_c = 0; pom_c < 100; pom_c++ ) {

    my_scores = &scores[numSeq*omp_get_thread_num()];
    
    my_hist = hists[omp_get_thread_num()];
    
    disp_this_pom = disp_p_pom * pom_c;
    disp_this_pmv = width * pom_c;
    
    //Rprintf("Analyzing pom %d\n", pom_c);
    
    /*Rprintf("Scanning POM ");
    for ( int idx=0; idx < width*4; idx += 4 ) {
      Rprintf("%f", POMs[ disp_this_pom + 4 * idx ]);
      Rprintf("-");
      Rprintf("%f", POMs[ disp_this_pom + 4 * idx + 1]);
      Rprintf("-");
      Rprintf("%f", POMs[ disp_this_pom + 4 * idx + 2 ]);
      Rprintf("-");
      Rprintf("%f", POMs[ disp_this_pom + 4 * idx + 3 ]);
      Rprintf(" ");

    }
    Rprintf("\n");*/
    
    
    /*score = 0;
    
    for ( nt_c = 0; nt_c < width*4; nt_c += 4 ){
      num = POMs[ disp_this_pom + 4 * nt_c + seqs[ nt_c ] ];
      den = PMVs[ disp_this_pmv + nt_c ];
      
      score += ( num - den );
    }
    
    //Rprintf("%f ", score);
    my_scores[0] = score;
    min=max=score;
    
    // Get scores.
    for ( seq_c = 1; seq_c < numSeq; seq_c++ ){
      score = 0;
      disp_this_seq = width * seq_c;
  
      
      for ( nt_c = 0; nt_c < width*4; nt_c += 4 ){
        num = POMs[ disp_this_pom + 4 * nt_c + seqs[ disp_this_seq + nt_c ] ];
        den = PMVs[ disp_this_pmv + nt_c ];
        
        score += ( num - den );
      }
      
      if( min > score )
        min = score;   
      if( max < score )
        max = score;
      
      //Rprintf("%f ", score);
      my_scores[seq_c] = score;
      //Rprintf("%f ", scores[seq_c]);
    }*/
    //Rprintf("Added %d scores\n", numSeq);
    
    /*min=max=my_scores[0];
    
    for(i=1; i<numSeq; i++)
    {
             
    }*/
    //Rprintf("POM %d\n", pom_c); 
    

   /* if ( min >= max ) {
	      min = max - 1;
        max = max + 1;
   }*/
    
    
    
    
    current_breaks = gsl_histogram_set_ranges_uniform(my_hist, min, max);
    
    for ( seq_c = 0; seq_c < numSeq; seq_c++ )
      gsl_histogram_increment( my_hist, my_scores[seq_c] );
    
    //current_breaks = gsl_histogram_set_ranges_uniform(my_hist, 0, numSeq);
    /*for ( seq_c = 0; seq_c < numSeq; seq_c++ )
      gsl_histogram_increment( my_hist, seq_c );*/
    
    
    disp_this_pom_breaks = nBins * numPOM + pom_c * ( nBins );
    disp_this_pom_counts = pom_c * nBins;
    
    
    //Rprintf("counts [%d, %d] and breaks [%d, %d]\n", 
    //         disp_this_pom_counts, disp_this_pom_counts + nBins -1, 
    //        disp_this_pom_breaks, disp_this_pom_breaks + nBins);
    //for ( i = 0; i < nBins; i++ ){
    //  Rprintf("%d ",  gsl_histogram_get(my_hist, (size_t) i));
    //}
    //Rprintf("\n");
    
    Rprintf("(%d pom) writing %d - %d and %d - %d (max %d)\n", pom_c, disp_this_pom_counts, nBins+disp_this_pom_counts, disp_this_pom_breaks, nBins+disp_this_pom_breaks, nBins * numPOM + ( nBins) * numPOM );
    /*for ( i = 0; i < nBins; i++ ) {
      break_counts [ i + disp_this_pom_counts ] = gsl_histogram_get(my_hist, (size_t) i);
      break_counts [ i + disp_this_pom_breaks ] = current_breaks[i];
    }*/
    
    
    gsl_histogram_reset( my_hist );
  }
  
  for ( i = 0; i < n_cpu; i++ )
    gsl_histogram_free(hists[i]);
  
  UNPROTECT(nProt);
  
  return BreaksCounts;
  
}
