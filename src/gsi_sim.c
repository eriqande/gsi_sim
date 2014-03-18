/*


	GSI_SIM  
	Author:  Eric C Anderson
	Date Initiated: 19 July 2005
	Copyright:  Federal Gov't Work.  No Copyright.
	
	Purpose:  This is a simple program that reads in allele frequencies from 
	a number of different populations (actually it will read allele counts and 
	will also need a dispersion quantity to turn those into dirichlet parameters)
	and then simulates individuals from each of those populations and prints the
	probability of its multilocus genotype given that it came from each of the 
	different populations.  
	
	Now, I have it set up to do and EM algorithm to estimate proportions and a
	few other things.  
	
	The output can easily be processed by awk to calculate posterior probabilities
	that each of the individuals came from any of the populations.




*/

#define UN_EXTERN



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "ECA_MemAlloc.h"
#include "MathStatRand.h"
#include "ranlib.h"
#include "ECA_Opt.h"
#include "MCTypesEtc.h"

#define MAX_POPS 1000
#define MAX_ALLELES 1000
#define MAX_ALLELE_LENGTH 999
#define MAX_LOC_NAME_LEN 500
#define MAX_POP_NAME_LENGTH 5000
#define MAX_FILE_NAME_LENGTH 5000
#define MAX_POST_PRED_RAN_FILES 50
#define MAX_MULTI_FISHERIES 100000
#define MAX_MULTI_MIX_NAME 100

/* a structure to hold information that applies to lots of entities */
typedef struct {
	int *NumAlle;
	int **AlleNames;
	int **LocHash;
	char **LocusNames;
	int *Ploidy;  /* array of ploidy levels for each locus */
	
} alle_info;

/* a structure that will hold a population index, and then a double value.  I use this so that we can quickly sort 
posterior probs, etc, as need be.  In fact, I will also put a dval in there, so that we can average over realizations of 
a markov chain and still have them with us in the end after we have sorted them.  Don't know that I will use that though. 
 
 I will also use these to hold values associated with Reporting Units */
typedef struct {
	int pop;
	double val;
	dval *DV;  
} popindval;

/* a structure to hold individuals */
typedef struct {
	int PopNum;  /* holds the number of the population it was entered with */
	int NumInPop;
	int NumLoc;
	char Name[1000];
	int **y; /* index of allelic types */
	int ***CompactRep; /* for a representation of allele counts that will be computationally more efficient
		                  when it comes to computing CDM densities, etc */
	int *AllelesAppearing;  /* the number of different allelic types appearing in each individual */
	popindval *posts;
	popindval *logls;
} ind_struct;



typedef struct {
	int Idx;
	int IdxInBaseline;  /* this is a weird one---if we are using a holdout set, we need a way of knowing which 
	                       population in the baseline the population in the holdout set corresponds to.  This is 
						   how we do it.  If pop 5 in the holdout set corresponds to pop 8 in the baseline, then
						   IdxInBaseline of pop 5 in the holdout set will be 8 */
	char Name[1000];
	
	double **x;  /* allele counts */
	double **xDP; /* dirichlet pars for simulating given the allele counts, Prior, Factor, and ScalePrior */
	double **xBP;  /* dirichlet pars for computing the probability of an individual using the WORST method */
	double **p;  /* allele frequencies */
	double **simp;  /* for storing the allele frequencies simulated from their posterior.  Used to correct Geneclass Bayesian Fallacy issues. */
	
	double *xsum;  /* the sum of x's at each locus */
	double *xDPsum;  /* the sum of xDP's at each locus */
	
	int NumLoc;
	int *NumAlle;
	int *Ploidy;  /* duplicating the Ploidy in the alle_info */
	
	double Factor;  /* Dirichlet parameters factor */  
	
	double Prior;  /* the prior weight for each allele, or if prior type is Scaled, then the numerator */
	
	int ScalePrior;
	
	double SimPi;  /* relative weight in a simulated mixture when using the -b option to specify baselines */
	int SimCnt;  /*absolute number of fish in a simulated mixture when using the -b option to specify baselines */
	
	/* here are some things for making "constellation" populations.  These are populations with
	   allele frequencies that cluster around the original population according the their Fst values
	   with the original population */
	int NumConstell;
	double *ConstelWts;
	double *ConstelFst;
	
	
	/* and here are some things for storing loss functions and computing expected losses */
	double *Loss;  /*!< vector of losses applied if the fish is assigned to different pops */
	double *ExpScaledLike; /*!< The expected probability of being assigned to any of the other pops scaled to sum to one, assuming equal priors.   */
	
	/* here we record information about the original individuals that may have been 
	summarized into this pop struct.  These two variables get set in IndColl2PopStructs as appropriate.  */
	struct IndCollection *OriginalIndColl;  /* if the pop_struct came from allele freqs and not inds, then this remains null */
	int IdxOfPopInIndColl;  /* default value should be -1.  */
	
	int RepUnitIdx;  /* the index of the reporting unit that this population belongs to */
	
} pop_struct;


typedef struct {
	int NumRepUnits;
	char **Names; /* names of each reporting unit */
	
	int *NumPopsInRepUnit;   /* number of pops in each of the NumRepUnits reporting groups */
	int **PopsInRepUnits; /* the indexes of each of the baseline populations in each of the reporting groups */
	int *RepUnitOfPops;  /* a vector to give the reporting unit of each pop */
	
} reporting_unit_info;


/* This is a structure which is intended to contain all the information that gets passed in,
and eventually computed, on a set of individuals who are passed in using, say, the -b option.
It holds the Individual genotypes themselves and the summaries thereof.  It also will hold the 
allele length hash, and the pop names, etc.  Each different population that gets read into the
program as a group of individuals, will also get summarized into a pop_struct. 
*/
struct IndCollection {

	int NumInds;  /* number of individuals in collection */
	int NumLoc; /* number of loci they are typed at */
	int NumPops;  /* the number of different populations in this collection */
	ind_struct **Inds; /* array of pointers to all the individuals */
	alle_info *AlleInfo;  /* pointer to information about locus names, allele names and the allele name hashes, etc. */
	char **PopNames;
	
	/* now, here are some things that give arrays of the traditional gsi_sim paramters like Factor and Scale
	for each of the populations.  See Factor, Prior, and ScalePrior in pop_struct struct to see what these are */
	double *PopFactors;
	double *PopPriors;
	int *PopScalePriors;
	
	/* We might want to store the indices of individuals that start and end populations,
	as well as the number of individuals in each population. */
	int *PopStarts;
	int *PopEnds;
	int *NumInPop;

	/* and here is a place to keep an array of pop_structs that sumamrize all these inds */
	pop_struct **PopSumms;
	
	popindval *Pi; /* These are mixing proportions of the different populations in a mixture */
	popindval *SimPi;
	
};

/* a silly little structure to hold the options for a GeneClass-like analysis. I will fill out the rest of this later */
typedef struct {
	char ReferenceFile[MAX_FILE_NAME_LENGTH];
	char SampleFile[MAX_FILE_NAME_LENGTH];
} gc_opts;

/* a silly little structure to hold options for EM-based gsi.  Eventually, it will hold 
things like max number of EM iterations, EM-completion criteria.  */
typedef struct {
	double tol;  /* the tolerance to require for a completion of the EM-algorithm */
	int MaxSteps;  /* maximum number of steps for the EM-algorithm. A step is considered to be an iteration of computing 
					the posterior probs of the individuals.  If MaxSteps=1, then the posterior probs will be identical to
					what you get with GeneClass (but the estimates of Pi will no longer be uniform).  If MaxSteps=2 then you 
					get posterior probs that are what you get when computing them using the Pi that you obtained at the end
					of MaxSteps=1, and so forth. */
	int MinSteps;  /* minimum number of steps for the EM-algorithm */
} em_opts;

typedef struct {
	int NumBurnIn;
	int NumSweeps; /* number to do after burn in */

	int PiRecordInterval;  /* how often to accumulate averages on Pi.  Not currently used. */
	int IndivRecordInterval;  /* how often to accumulate averages on individuals. Not currently used. */

	int PiTraceInterval; /* if 0, don't print any trace. If 1 or more, print out at that interval */
	int IndivTraceInterval; /* same as above, but for individuals.  This could be Huge output */
	
	int ZSumTraceInterval;  /* if 0 don't print any trace.  If 1 or more print out at that interval.  This is for the sum of the Zs */
	
} BayesianOpts;

typedef enum {
	RESAMP_LOO_MULTILOCUS,  /* this resamples the vectors of scaled likelihood of individuals to form mixtures. This is equivalent to
							resampling multilocus genotypes and then computing their likelihoods leaving their genes out of the "gene-frequency pile." */
	RESAMP_LOO_SINGLELOCUS,
	RESAMP_LOO_GENECOPY, /* resample gene copies w/o replacement, and compute likelihood using leave-one-out (where you leave out all
							the particular gene copies sampled.  Suggested by Steven Kalinowski */
	NO_LOO_MULTILOCUS,
	NO_LOO_SINGLELOCUS,
	NO_LOO_GENECOPY
} Method;


typedef struct {
	int NumFiles;  /* number of files with posterior samples of sample sizes */
	int *N;   /* number of entries read from each file */
	int **V;  /* V[i] is an array of the values sampled from file i */
	char **Suf;  /* the suffixes to be used for the different files */
} pred_ran_stuff;




/* a global variable or two */
int gVerbosity = 0;
int gPrintMixGenotypes = 0;
char gRepUnitFile[1000];
int gHasRepUnits = 0;
reporting_unit_info* gRepUnitInfo=NULL;
int gBaselinePairLocusLimit = -1;  /* number of non missing loci required to entertain a pair in the FindCloseMatchers function.  If -1 it means don't find close matchers in Baseline */
int gBaselinePairMismatchLimit = -1;
int gMixturePairLocusLimit = -1;  /* number of non missing loci required to entertain a pair in the FindCloseMatchers function.  If -1 it means don't find close matchers in Mixture */
int gMixturePairMismatchLimit = -1;
pred_ran_stuff *gPredRan_p=NULL;


/* prototypes */
int *SlurpPostPredRanFile(char *File, int *n);
void ReadLocFile(pop_struct **P, int L, char *FileName);
pop_struct *InitPopStruct(int I);
void SimInd(int J, int N, pop_struct **Pops, pop_struct **SimPops, double *Outs, Method Meth);
int GetGSI_Options(pop_struct **P, int argc, char *argv[], double **TruePi, int *INDS, int *MIXES, 
	int *MIXSIZE, Method *Meth, struct IndCollection **Blines, struct IndCollection **Mixture, struct IndCollection **Holdout, int *SelfAss, int **FixedNums, int *NoEM, int *ProperBayesWayReps,
	int *NumMixtureLoglSimReps, int *NumBaselinLoglSimReps, int *DumpMixLoglFile, int *DumpBaselineLoglFile,  BayesianOpts **BayesOpts, int **PostPredNs, int *NumPredPostNs,
  char ***PostPredRanFiles, int *NumPredRanFiles, char ***PostPredRanSuffixes, int ***MultiFixMix, char ***MultiFixMixNames, int *NumMultiFixMixes);
void EM_Find_Pi_MLE(double **M, int N, int C, double *SP, double *Pi, double TOL);
double PairwiseFst(pop_struct *A, pop_struct *B);
void StraightUpSim(int INDS, int NumPops, pop_struct **Pops, int NumSimPops, pop_struct **SimPops, Method Meth, double *TruePi);
void MixtureSim(int INDS, int NumPops, pop_struct **Pops, int NumSimPops, pop_struct **SimPops, Method Meth, double *TruePi, 
                int MIXSIZE,int MIXES,  struct IndCollection *Baselines, int *FixedNums,
                const char MultiFixName[], const char *MultiFixPrefix, const char *MultiFixNameHeader);
void PrintLossMatrix(pop_struct **P, int NumPops);
void ExpectedLossLocL(int L, pop_struct **P, double *Pi, int NumPops);
void PrintIndCollection(struct IndCollection *R);
void IndColl2PopStructs(struct IndCollection *R, struct IndCollection *BaseReference);
struct IndCollection *InitIndCollection(void);
long double ProbIndFromPopR_M_SingleLocusDisomic(struct IndCollection *R, ind_struct *T, int p, int LeaveOut, int j, long double *OutNumer);
void ComputePostVector(struct IndCollection *R, ind_struct *T, double *Pi, int LOGLS, int TryToLeaveOut, int CorrectWayReps);
int CmpPopIndByVal(const void *a, const void *b);
void RealizeNewSimP(int NumPops, pop_struct **PopSumms);
int CmpPopIndByPop(const void *a, const void *b);
void ComputeGeneClassSample(struct IndCollection *Baselines, struct IndCollection *TheMixture, gc_opts *GeneClassOpts, int CorrectWayReps, int Sort);
void ComputeGeneClassSelfAssign(struct IndCollection *Baselines, gc_opts *GeneClassOpts);
void ComputeGeneClassSelfAssign_NOSORT(struct IndCollection *Baselines, gc_opts *GeneClassOpts);
void PrintGeneClassCSV(struct IndCollection *R, struct IndCollection *S,  gc_opts *GeneClassOpts, const char *Tag);
void FprintIndivPosteriorFields(struct IndCollection *Baselines, struct IndCollection *TheMixture, const char *PopFileName, reporting_unit_info *RU, const char *RepFileName);
void ComputeGSISample(struct IndCollection *Baselines, struct IndCollection *TheMixture, em_opts *EMOpts, int Sort);
long double ProbIndFromPopR_M_Polysomic(struct IndCollection *R, ind_struct *T, int p, int LeaveOut);
void InitPopPriorsForIndCollection(struct IndCollection *R);
void UpdateBaselinePopPriors(struct IndCollection *R, int Num, pop_struct **Pops);
double *pi_vecFromPopStructPis(struct IndCollection *R, int Num, pop_struct **Pops, int **TrueCnts);
void ResampleBaselines(int NumPops, pop_struct **Pop);
void DrawSingleLocus(int A, int *Y, int Ploidy, int LocIdx, struct IndCollection *Orig, int Idx );
void SimulateAllLoglsForComparison(struct IndCollection *Baselines, struct IndCollection *TheMixture, Method Meth, int NumReps, FILE *loglf, char *outfile_name);
double SimIndForLoglDsnFromSinglePop(int J, pop_struct **Pops, pop_struct **SimPops, Method Meth, int K, int **Yobs, int *NumTyped);
int IdxOfPopName(struct IndCollection *B, char *Name);
int IdxOfString(char *StrArray[], int n, char Name[]);
reporting_unit_info* GatherReportingUnitData(char *RepUFile, struct IndCollection* B);
FILE *OpenTraceFile(const char *FileName, const char *Comment, struct IndCollection *Baselines, BayesianOpts *BO, reporting_unit_info *RU, int interval);
void DoAndPrintPostPredSampling(int NumPredPostNs, int *PredPostNs, double *Pi, int NumPops, reporting_unit_info *RU, int step, FILE **popf, FILE **repf);
void DoAndPrintPostPredRansSampling(pred_ran_stuff *prs, double *Pi, int NumPops, reporting_unit_info *RU, int step, FILE **popf, FILE **repf);
void ComputeGSISample_bayesian(struct IndCollection *Baselines, struct IndCollection *TheMixture, BayesianOpts *BO, reporting_unit_info *RU, int NumPredPostNs, int *PostPredNs);
void FprintPis(struct IndCollection *Baselines, const char *PopFileName, reporting_unit_info *RU, const char *RepFileName);
void FindCloseMatchers(struct IndCollection *C, int S, int NM,  char *outfile_name);



/* slurp the ints out of a PredRan file.  Return them in an array.  Also,
 send back n --- the number of elements in the returned array */
int *SlurpPostPredRanFile(char *File, int *n) 
{
	int i=0, junk;
	FILE *in;
	int *ret;
	
	if( (in=fopen(File,"r"))==NULL) {
		fprintf(stderr, "Failed to open file %s to read pred rans stuff.  Exiting!\n\n", File);
		exit(1);
	}
	while(!feof(in)) {
		fscanf(in," %d", &junk);
		i++;
	}
	ret = (int *)calloc(i, sizeof(int));
	*n = i-1;
	rewind(in);
	i=0;
	while(!feof(in)) {
		fscanf(in," %d", &(ret[i]));
		i++;
	}
	if(i != (*n)+1) {
		fprintf(stderr, "Bad News!! The second time through file %s we read %d instead of %d items. Exiting\n\n", File, i, (*n));
		exit(1);
	}
	return(ret);
	
}

int main(int argc, char *argv[]) 
{
	int i,NumPops=0,NumSimPops=0,j,INDS=0, MIXES=0, MIXSIZE=0;
	pop_struct **Pops=ECA_CALLOC(MAX_POPS, sizeof(pop_struct *));
	pop_struct **SimPops=ECA_CALLOC(MAX_POPS, sizeof(pop_struct *));
	double *TruePi=NULL;
	int *FixedNums=NULL;  /* if we want to simulate mixtures with fixed numbers of indivs, we use this */
	Method Meth=RESAMP_LOO_GENECOPY;  /* make CV-GC the default option */
	struct IndCollection *Baselines=NULL;  /* for the GeneClass/GMA emulation */ 
	struct IndCollection *TheMixture=NULL;
	struct IndCollection *TheHoldout=NULL;
	em_opts EM_Opts;
	int SelfAssign=0;
	int NoEM = 0; /* default is to do the EM if there is a baseline and mixture.  Setting this to 1 makes it
					just assign individuals assuming equal prior probs for each population. */
	int ProperBayesWayReps=0; /* whether to do the geneclass approach properly (>0--integrate the ratio over p using ProperBayesWayReps reps) or the fallacious way (0--rannala-mountain) */ 
	int NumMixtureLoglSimReps=0, NumBaselinLoglSimReps=0;  /* this is a flag which gives the number of LogL Simulations to do to copmare the raw log like values.  Default is zero.  Can be set through command line */
	int DumpMixLoglFile=0, DumpBaselineLoglFile=0;  /*  flags as to whether the histograms of the the LogL simulations should be dumped to a text file */
	BayesianOpts *BayesOpts;
	int *PredPostNs, NumPredPostNs;
	char **PostPredRanFiles;
	char **PostPredRanSuffixes;
	int NumPredRanFiles=0;
	pred_ran_stuff PRS;
  int **MultiFixMix = NULL;
  char **MultiFixMixNames = NULL;
  int NumMultiFixMixes = 0;
	
	
	
		
	/* some em_opts defaults for now */
	EM_Opts.MaxSteps=1000;
	EM_Opts.MinSteps=50;
	EM_Opts.tol=.000001;
	PRS.NumFiles=0;
	
	gPredRan_p = &(PRS);  /* make this global pointer so I don't have to pass it into functions way down the road */
	
	
	/* get all the inputs */
	NumPops = GetGSI_Options(Pops, argc, argv, &TruePi, &INDS, &MIXES, &MIXSIZE, &Meth, &Baselines, &TheMixture,&TheHoldout, &SelfAssign, &FixedNums, &NoEM, &ProperBayesWayReps,&NumMixtureLoglSimReps,
							 &NumBaselinLoglSimReps, &DumpMixLoglFile, &DumpBaselineLoglFile, &BayesOpts, &PredPostNs, &NumPredPostNs, &PostPredRanFiles, &NumPredRanFiles, &PostPredRanSuffixes, 
                           &MultiFixMix, &MultiFixMixNames, &NumMultiFixMixes);

	for(i=0;i<NumPredPostNs;i++)  {
		printf("PRED POST: %d  %d\n",i,PredPostNs[i]);
	}
	
	/* down here we process the Pred Rans Stuff */
	if(NumPredRanFiles) {
		PRS.NumFiles = NumPredRanFiles;
		PRS.V = (int **)calloc(NumPredRanFiles, sizeof(int *));
		PRS.N = (int *)calloc(NumPredRanFiles, sizeof(int));
		PRS.Suf = PostPredRanSuffixes;
		for(i=0;i<NumPredRanFiles;i++)  {
			printf("PRED_POST_RANS: Preparing to read file %s tied to suffix %s.  File %d of %d\n",PostPredRanFiles[i], PostPredRanSuffixes[i], i+1, NumPredRanFiles);
			PRS.V[i] = SlurpPostPredRanFile(PostPredRanFiles[i],  &(PRS.N[i]) );
			printf("PRED_POST_RANS: Extracted %d integers [%d, %d, %d, ..., %d, %d] from that file!\n", PRS.N[i], PRS.V[i][0], PRS.V[i][1], PRS.V[i][2], PRS.V[i][PRS.N[i]-2], PRS.V[i][PRS.N[i]-1]);
		}
	}
	

	/* convert the baseline into a summary */
	if(Baselines != NULL) {
		InitPopPriorsForIndCollection(Baselines);
		UpdateBaselinePopPriors(Baselines, NumPops, Pops);
		IndColl2PopStructs(Baselines,NULL);
		printf("IND_COLL_SUMMARY_FOR: BASELINE\n");
		PrintIndCollection(Baselines);
		TruePi = pi_vecFromPopStructPis(Baselines,NumPops,Pops,&FixedNums);
		
		if(gBaselinePairLocusLimit>=0) {
			FindCloseMatchers(Baselines, gBaselinePairMismatchLimit, gBaselinePairLocusLimit,  "close_matches_baseline.txt");
		}
		
	}
	
	
	
	if(gHasRepUnits==1) {
		printf("REPORTING_UNITS:  Preparing to collect from File : %s\n",gRepUnitFile);
		gRepUnitInfo = GatherReportingUnitData(gRepUnitFile, Baselines);
	}
	

	
	
	
	/* also, convert the mixture into a summary */
	if(TheMixture != NULL) {
		InitPopPriorsForIndCollection(TheMixture);
		IndColl2PopStructs(TheMixture,NULL);
		printf("IND_COLL_SUMMARY_FOR: MIXTURE\n");
		PrintIndCollection(TheMixture);
		
		if(gMixturePairLocusLimit>=0) {
			FindCloseMatchers(TheMixture, gMixturePairMismatchLimit, gMixturePairLocusLimit,  "close_matches_mixture.txt");
		}
	}
	
	
	if(gMixturePairLocusLimit>=0 || gBaselinePairLocusLimit>=0 ) {
		printf("Terminating Execution After Searching For Close Matchers\n\n");
		return(0);
	}
	
	/* And convert the holdout into a summary as well */
	if(TheHoldout != NULL) {
		InitPopPriorsForIndCollection(TheHoldout);
		IndColl2PopStructs(TheHoldout,Baselines);
		printf("IND_COLL_SUMMARY_FOR: HOLDOUT\n");
		PrintIndCollection(TheHoldout);
	}

	
	/* set the RNG */
	SeedFromFile("gsi_sim_seeds");
	
	
	/* If we have both a baseline and a mixture sample, then analyze the mixture! */
	if(TheMixture != NULL && Baselines != NULL) { 
		
		/* testing the bayesian stuff */
		if(BayesOpts->NumSweeps>0) {
			ComputeGSISample_bayesian(Baselines,TheMixture,BayesOpts,gRepUnitInfo, NumPredPostNs, PredPostNs);
			
			/* here we do the LogL simulation based on the population assignment in the posterior means if requested */
			if(NumMixtureLoglSimReps>0) { FILE *loglf;
				loglf=NULL;
				if(DumpMixLoglFile) {
					loglf=fopen("bayes_mixture_logl_compare_dump_file.txt","w");
				}
				if(loglf) {
					fprintf(loglf,"FishID\tLogL\n");
				}
				SimulateAllLoglsForComparison(Baselines, TheMixture, RESAMP_LOO_GENECOPY, NumMixtureLoglSimReps, loglf, "bayes_mixture_logl_summary.txt");
				fclose(loglf);
			}
		}
		
		/* do the GeneClass Bit */
		ComputeGeneClassSample(Baselines,TheMixture,NULL,ProperBayesWayReps,0);
		PrintGeneClassCSV(Baselines,TheMixture,NULL,"UNSORTED_GENE_CLASS_CSV:");
		FprintIndivPosteriorFields(Baselines, TheMixture, "pop_pofz_scaled_likelihoods.txt", gRepUnitInfo, "rep_unit_pofz_scaled_likelihoods.txt");
		
		ComputeGeneClassSample(Baselines,TheMixture,NULL,ProperBayesWayReps,1);
		PrintGeneClassCSV(Baselines,TheMixture,NULL,"GENECLASS_CSV:");
		
		
		if( !(NoEM) ) {	
			/* do the Full-EM conditional GSI */
			ComputeGSISample(Baselines, TheMixture, &EM_Opts,0);
			FprintPis(Baselines, "pop_pi_full_em_mle.txt", gRepUnitInfo, "rep_unit_pi_full_em_mle.txt");
			PrintGeneClassCSV(Baselines,TheMixture,NULL,"UNSORTED_GMA_FULL_EM_IND_CSV:");  /* give us the output in canonical baseline population order instead of sorting for each individual */
			FprintIndivPosteriorFields(Baselines, TheMixture, "pop_pofz_full_em_mle.txt", gRepUnitInfo, "rep_unit_pofz_full_em_mle.txt");
			
			
			ComputeGSISample(Baselines, TheMixture, &EM_Opts,1);
			qsort( (void *)(Baselines->Pi), (size_t)(Baselines->NumPops), sizeof(popindval), CmpPopIndByVal);
			for(i=0;i<Baselines->NumPops;i++)  {
				printf("GMA_FULL_EM_MLE_OF_PI: %s   %f\n",Baselines->PopNames[Baselines->Pi[i].pop], Baselines->Pi[i].val);
			}
			PrintGeneClassCSV(Baselines,TheMixture,NULL,"GMA_FULL_EM_INDIVS_CSV:");
			
			/* here we do the LogL simulation if requested */
			if(NumMixtureLoglSimReps>0) { FILE *loglf;
				loglf=NULL;
				if(DumpMixLoglFile) {
					loglf=fopen("em_mixture_logl_compare_dump_file.txt","w");
				}
				if(loglf) {
					fprintf(loglf,"FishID\tLogL\n");
				}
				SimulateAllLoglsForComparison(Baselines, TheMixture, RESAMP_LOO_GENECOPY, NumMixtureLoglSimReps, loglf, "em_mixture_logl_summary.txt");
				fclose(loglf);
			}
			
			/* do the single step of EM conditional GSI */
			EM_Opts.MaxSteps=2;
			ComputeGSISample(Baselines, TheMixture, &EM_Opts,0);
			PrintGeneClassCSV(Baselines,TheMixture,NULL,"UNSORTED_GMA_ONE_STEP_EM_IND_CSV:");
			FprintIndivPosteriorFields(Baselines, TheMixture, "pop_pofz_one_step_em_mle.txt", gRepUnitInfo, "rep_unit_pofz_one_step_em_mle.txt");
			
			ComputeGSISample(Baselines, TheMixture, &EM_Opts,1);
			qsort( (void *)(Baselines->Pi), (size_t)(Baselines->NumPops), sizeof(popindval), CmpPopIndByVal);
			for(i=0;i<Baselines->NumPops;i++)  {
				printf("GMA_ONE_STEP_EM_MLE_OF_PI: %s   %f\n",Baselines->PopNames[Baselines->Pi[i].pop], Baselines->Pi[i].val);
			}
			PrintGeneClassCSV(Baselines,TheMixture,NULL,"GMA_ONE_STEP_EM_INDIVS_CSV:");
		}
		
	}
	if(SelfAssign==1) {
		if(Baselines != NULL) {  
			ComputeGeneClassSelfAssign(Baselines,NULL);
			PrintGeneClassCSV(Baselines,Baselines,NULL,"SELF_ASSIGN_A_LA_GC_CSV:");
			
			/* here we do the LogL simulation if requested */
			if(NumBaselinLoglSimReps>0) { FILE* loglf;
				loglf=NULL;
				if(DumpBaselineLoglFile) {
					loglf=fopen("baseline_logl_compare_dump_file.txt","w");
				}
				if(loglf) {
					fprintf(loglf,"FishID\tLogL\n");
				}
				SimulateAllLoglsForComparison(Baselines, Baselines, RESAMP_LOO_GENECOPY, NumBaselinLoglSimReps, loglf, "baseline_logl_summary.txt");
				fclose(loglf);
			}
			
			/* and then, as long as we are here we will do the self assignment but not sort things */
			ComputeGeneClassSelfAssign_NOSORT(Baselines,NULL);
			PrintGeneClassCSV(Baselines,Baselines,NULL,"UNSORTED_SELF_ASS_LIKE_GC_CSV:");
		}
	}
	
	
	
	
	
///ExpectedLossLocL(0, Pops, TruePi,NumPops);
//return(0);


	/* down here, if we were given data in a baseline format (-b) then we have to set Pops and NumPops
		before sending everything along, and we will set NumSimPops and SimPops as well, though if there
		is a holdout set, then we will reassign new values to those variables. */
	if(Baselines != NULL)  {
		NumPops = Baselines->NumPops; 
		Pops = Baselines->PopSumms;
		NumSimPops=NumPops;
		SimPops = Baselines->PopSumms;
	}
	
	/* if there is a holdout set, then that is the one that we are going to be using for simulating individuals */
	if(TheHoldout != NULL) {
		NumSimPops=TheHoldout->NumPops;
		SimPops = TheHoldout->PopSumms;
	}
	
	/*************************
	Here we print out the haploid pairwise Fst's between the populations
	**************************/
	for(i=0;i<NumPops;i++)  {
		for(j=i+1;j<NumPops;j++)  {
			printf("PAIRWISE_FST :  %s  :  %s  :  %f\n",Pops[i]->Name,Pops[j]->Name,PairwiseFst(Pops[i], Pops[j]));
		}
	}

	/*******
	print the loss matrix only if it is actually there
	*****/
	if(Pops[0]->Loss != NULL)  {
		PrintLossMatrix(Pops, NumPops);
	}
	
	/*************
	This makes INDS straight-up simmed individuals
	*************/
	StraightUpSim(INDS,NumPops,Pops,NumSimPops,SimPops,Meth,TruePi);
	
	
	
	/*********
	This makes MIXES mixtures of size MIXSIZE 
	*******/
	if(TruePi!=NULL) {
		if(Meth==RESAMP_LOO_MULTILOCUS) {  /* if we are doing the shrewd method, then we need to compute all the scaled likelihoods !!but don't sort them!!
							before we enter the mixture sim */
			ComputeGeneClassSelfAssign_NOSORT(Baselines, NULL);
		}
		MixtureSim(INDS,NumPops, Pops, NumSimPops, SimPops, Meth, TruePi, MIXSIZE, MIXES,Baselines,FixedNums, NULL, NULL, " ");
	}


  /********
  Finally, if requested, bang out the MultiFixMix simulations
  *********/
  if(Meth==RESAMP_LOO_MULTILOCUS) {  /* if we are doing the shrewd method, then we need to compute all the scaled likelihoods !!but don't sort them!!
                                      before we enter the mixture sim. Just sort of threw this down here. */
    ComputeGeneClassSelfAssign_NOSORT(Baselines, NULL);
  }
  for(i=0; i<NumMultiFixMixes; i++)  {
      MixtureSim(INDS,NumPops, Pops, NumSimPops, SimPops, Meth, TruePi, 1, 1,Baselines, MultiFixMix[i], MultiFixMixNames[i], "MULTI_FIX_MIX_", "MultiFixMixName");
  }
	
	
	SeedToFile("gsi_sim_seeds");

	return(0);
}






/* simple function to do all pairwise comparisons within an IndCollection C and find 
 pairs that mismatch at X or fewer genotypes out of at least NM diploid loci that are not
 missing in either member of the pair.
 
 The results get written into the file outfile_name
 
 Currently this is implemented only for the diploid loci
 
 */
void FindCloseMatchers(struct IndCollection *C, int X, int NM,  char *outfile_name)
{
	int i,j,k;
	int N = C->NumInds;
	int L = C->NumLoc;
	int non_missers;
	int mismatches;
	int **x, **y; 
	int a,b,c,d;
	int *TotNonMiss;
	int NumDiploid=0;
	int TossIt;
	FILE *out;
	
	/* make a header for the file */
	if( (out=fopen(outfile_name,"w"))==NULL) {
		fprintf(stderr,"Failed to open file %s to write to it.  Exiting to system\n",outfile_name);
		exit(1);
	}
	
	fprintf(out,"# Pairs of fish having %d or fewer mismatching genotypes and at least %d loci that are non-missing in each member of the pair \n",X,NM);
	fprintf(out,"# Each pair should be listed at most once. Indiv1 appears before Indiv2 in the data set.\n");
	fprintf(out,"Indiv1\tIndiv2\tNumMismatchLoci\tNumNonMissingLoci\n"); 
	
	/* first, cycle over all the individuals and record how much non-missing data each has */
	TotNonMiss = (int *)calloc(N,sizeof(int));
	for(i=0;i<N;i++) {
		non_missers = 0;
		for(j=0;j<L;j++)  {
			if(C->AlleInfo->Ploidy[j]==2 && C->Inds[i]->y[j][0]>=0 && C->Inds[i]->y[j][1]>=0) {
				non_missers++;
			}
		}
		TotNonMiss[i] = non_missers;
		/*printf("CLOSE_MATCH: TotNonMiss of Fish %d is %d\n",i,TotNonMiss[i]); */
	}
	
	/* count the number of diploids loci */
	for(j=0;j<L;j++)  {
		if(C->AlleInfo->Ploidy[j]==2)  NumDiploid++;
	}
	
	
	
	for(i=0;i<N;i++) {
		/*printf("CLOSE_MATCH: Starting Matching of Fish %d\n",i); */
				
		if(TotNonMiss[i]>=NM) for(k=i+1;k<N;k++)  {
			non_missers=0;
			mismatches=0;
			x = C->Inds[i]->y;
			y = C->Inds[k]->y;
			TossIt = 0;
			
			if(TotNonMiss[j]>=NM) for(j=0;j<L;j++)  {
				if(C->AlleInfo->Ploidy[j]==2) {
					if(x[j][0]>=0 && x[j][1]>=0 && y[j][0]>=0 && y[j][1]>=0) {
						non_missers++;
						if(x[j][0]<=x[j][1])  { 
							a = x[j][0]; 
							b = x[j][1];
						}
						else {
							a = x[j][1]; 
							b = x[j][0];
						}
						if(y[j][0]<=y[j][1])  { 
							c = y[j][0]; 
							d = y[j][1];
						}
						else {
							c = y[j][1]; 
							d = y[j][0];
						}
						if(!(a==c && b==d)) {
							mismatches++;
						}
					}
				}
				if(mismatches > X) {
					TossIt=1;
					j=L;
					/*printf("CLOSE_MATCH: Just tossed %d %d.  mismatches=%d  j=%d  NumDiploid= %d  NM= %d\n",i,k, mismatches, j, NumDiploid, NM); */
					
				}
			}
			if(TossIt==0 && non_missers>=NM && mismatches <= X) {
				printf("CLOSE_MATCH_GENOTYPES: %s and %s differ at %d genotypes out of %d non-missing loci\n",C->Inds[i]->Name, C->Inds[k]->Name, mismatches, non_missers);
				fprintf(out,"%s\t%s\t%d\t%d\n",C->Inds[i]->Name, C->Inds[k]->Name, mismatches, non_missers);
			}
		}  /* closes loop over k */
	}
	
	fclose(out);
	free(TotNonMiss);
}



/* just a function to do the posterior predictive sampling over allocations of other fish with no genotypes.  
 This prints the results to the files whose handles are in the arrays of file handles.  
 
 If RU is non-null then it takes each realization with the pops and summarizes them into reporting units, too.
 */
void DoAndPrintPostPredSampling(int NumPredPostNs, int *PredPostNs, double *Pi, int NumPops, reporting_unit_info *RU, int step, FILE **popf, FILE **repf)
{
	int i,j;
	int PredZSum_pop[MAX_POPS], PredZSum_rep[MAX_POPS];
	
	for(i=0;i<NumPredPostNs;i++)  {
		/* first simulate the numbers in each population */
		D_MultinomialRV(PredPostNs[i], Pi ,NumPops, PredZSum_pop);
		
		/* then print those to the appropriate output file */
		fprintf(popf[i],"%d",step);
		for(j=0;j<NumPops;j++)  {
			fprintf(popf[i],"\t%d", PredZSum_pop[j]);
		}
		fprintf(popf[i],"\n");
		
		/* then, if RU is non-null, agglomerate them into reporting units and print those too */
		if(RU) {
			for(j=0;j<RU->NumRepUnits;j++)  {  /* initialize to accumulate a sum */
				PredZSum_rep[j] = 0;
			}
			for(j=0;j<NumPops;j++)  {
				PredZSum_rep[RU->RepUnitOfPops[j]] += PredZSum_pop[j];
			}
			fprintf(repf[i],"%d",step);
			for(j=0;j<RU->NumRepUnits;j++)  {  /* initialize to accumulate a sum */
				fprintf(repf[i],"\t%d",PredZSum_rep[j]);
			}
			fprintf(repf[i],"\n");
		}
		
	}
	
}





/* just a function to do the posterior predictive sampling over allocations of other fish with no genotypes.  
 This prints the results to the files whose handles are in the arrays of file handles.  
 
 If RU is non-null then it takes each realization with the pops and summarizes them into reporting units, too.
 
 This version uses the files of sample sizes so that we can propagate uncertainty in sample sizes too.
 */
void DoAndPrintPostPredRansSampling(pred_ran_stuff *prs, double *Pi, int NumPops, reporting_unit_info *RU, int step, FILE **popf, FILE **repf)
{
	int i,j;
	int PredZSum_pop[MAX_POPS], PredZSum_rep[MAX_POPS];
	int TheExN;
	
	for(i=0;i<prs->NumFiles;i++)  {
		/* first we sample from the number that we want to expand to: */
		TheExN = prs->V[i][ UniformRV(0, prs->N[i]-1) ];
		
		/* second simulate the numbers in each population */
		D_MultinomialRV(TheExN, Pi ,NumPops, PredZSum_pop);
		
		/* then print those to the appropriate output file */
		fprintf(popf[i],"%d",step);
		for(j=0;j<NumPops;j++)  {
			fprintf(popf[i],"\t%d", PredZSum_pop[j]);
		}
		fprintf(popf[i],"\n");
		
		/* then, if RU is non-null, agglomerate them into reporting units and print those too */
		if(RU) {
			for(j=0;j<RU->NumRepUnits;j++)  {  /* initialize to accumulate a sum */
				PredZSum_rep[j] = 0;
			}
			for(j=0;j<NumPops;j++)  {
				PredZSum_rep[RU->RepUnitOfPops[j]] += PredZSum_pop[j];
			}
			fprintf(repf[i],"%d",step);
			for(j=0;j<RU->NumRepUnits;j++)  {  /* initialize to accumulate a sum */
				fprintf(repf[i],"\t%d",PredZSum_rep[j]);
			}
			fprintf(repf[i],"\n");
		}
		
	}
	
}






/* a little function to open files and write a few things in them for the traces.  If RU is not NULL then it 
 sets the headers up for Reporting units.  Otherwise it does it for the populations.  */
FILE *OpenTraceFile(const char *FileName, const char *Comment, struct IndCollection *Baselines, BayesianOpts *BO, reporting_unit_info *RU, int interval)
{
	int j;
	FILE *outf;
	
	if(RU) {
		if( (outf = fopen(FileName,"w"))==NULL) {
			fprintf(stderr,"Failed to open file %s to write to it.  Fatal Error.  Exiting\n",FileName);
			exit(1);
		}
		fprintf(outf,"## %s  NumBurnIn= %d  NumSweeps= %d ThinningInterval= %d\n",Comment,BO->NumBurnIn, BO->NumSweeps, interval);
		fprintf(outf,"SweepNumber");  /* printing our header */
		for(j=0;j<RU->NumRepUnits;j++) {
			fprintf(outf,"\t%s",RU->Names[j]);
		}
		fprintf(outf,"\n");
	}
	else {
		if( (outf = fopen(FileName,"w"))==NULL) {
			fprintf(stderr,"Failed to open file %s to write to it.  Fatal Error.  Exiting\n",FileName);
			exit(1);
		}
		fprintf(outf,"## %s  NumBurnIn= %d  NumSweeps= %d ThinningInterval= %d\n",Comment,BO->NumBurnIn, BO->NumSweeps, interval);
		fprintf(outf,"SweepNumber");  /* printing our header */
		for(j=0;j<Baselines->NumPops;j++) {
			fprintf(outf,"\t%s",Baselines->PopSumms[j]->Name);
		}
		fprintf(outf,"\n");
	}
	return(outf);
}

/* this function simply computes GSI-related output (only Rannala-Mountain implemented so far)
 for TheMixture given reference samples of Baselines.  Note that the data have to be in the
 structures, all proper and ready to go.  This function also sorts the Posterior probs into 
 descending order.  If Sort is 0, then the posterior vector of each individual is NOT sorted
 from high to low posterior probability.  Otherwise it is. 
 
 This version does a simple bayesian mcmc for the mixing proportions (I need to add mcmc structs to the mixing proportions and
 the individuals still, this is exploratory. ).   Also, various parameters are totally hard-wired at the moment.  This is
 just exploratory to see if it works.  
 
 Check it out---apparently the EM thing recalculated the posteriors vector at each step.  Totally bogus...
 I am going to precompute all of those and then it will scream.  */
void ComputeGSISample_bayesian(struct IndCollection *Baselines, struct IndCollection *TheMixture, BayesianOpts *BO, reporting_unit_info *RU, int NumPredPostNs, int *PredPostNs)
{
	int j,i,steps;
	double init;
	double *Pi = (double *)ECA_CALLOC(Baselines->NumPops,sizeof(double));
	int dummy=0;
	double **M; /* matrix indexed by mix indiv and then by pop */
	int N=TheMixture->NumInds;  /* number of fish in the mixture */
	int P=Baselines->NumPops;  /* number of pops in baseline */
 	int *Z;  /* Z[i] is the population to which indiv i is allocated */
	double *PZ = (double *)ECA_CALLOC(Baselines->NumPops,sizeof(double)); /* the number of individuals currently allocated to pop i */
	int *ZSum = (int *)ECA_CALLOC(Baselines->NumPops,sizeof(int));
	double *PofZ = (double *)ECA_CALLOC(Baselines->NumPops,sizeof(double)); /* variable to be used temporarily for each fish in the mixture */
	double normo;
	dval **RepUnitPis;
	dval ***RepUnitPofZs;
	dval **PopPis;
	int *RepUnitZSum;  /* sum of Zs */
	
	FILE *outf;
	FILE *popZSumtracef;  /* for traces of the Sum of Zs */
	FILE *repZSumtracef;
	FILE *popPi_tracef;
	FILE *repPi_tracef;
	FILE **PostPredFiles_pop;
	FILE **PostPredFiles_rep;
	FILE **PostPredRanFiles_pop;
	FILE **PostPredRanFiles_rep;	
	
	PopPis = DvalVector(0,P,0.0,1.0,.01);
	if(RU!=NULL) {
		RepUnitZSum = (int *)calloc(RU->NumRepUnits,sizeof(int));
		RepUnitPis = DvalVector(0,RU->NumRepUnits, 0.0, 1.0, .01);
		RepUnitPofZs = (dval ***)calloc(N,sizeof(dval **));
		for(i=0;i<N;i++) {
			RepUnitPofZs[i] = DvalVector(0,RU->NumRepUnits-1,  0.0,1.0,   -1);  /* no histograms for these */
		}
		if(BO->PiTraceInterval>0) {
			repPi_tracef = OpenTraceFile("rep_unit_pi_trace.txt", "trace file of reporting unit Pi values from gsi_sim.", Baselines, BO, RU, BO->PiTraceInterval);
		}
		if(BO->ZSumTraceInterval>0) {
			repZSumtracef = OpenTraceFile("rep_unit_zsum_trace.txt", "trace file of reporting unit ZSum values from gsi_sim.", Baselines, BO, RU, BO->ZSumTraceInterval);
		}
	}
	
	/* now we also prepare files to print the population Pi and ZSum traces */
	if(BO->PiTraceInterval>0) {
		popPi_tracef = OpenTraceFile("pop_pi_trace.txt", "trace file of population Pi values from gsi_sim.", Baselines, BO, NULL, BO->PiTraceInterval);
	}
	if(BO->ZSumTraceInterval>0) {
		popZSumtracef =  OpenTraceFile("pop_zsum_trace.txt", "trace file of population ZSum values from gsi_sim.", Baselines, BO, NULL, BO->ZSumTraceInterval);
	}
	
	/* now, we also open the trace files (if needed) for the draws from the posterior predictive of the allocations */
	if(BO->ZSumTraceInterval>0 && NumPredPostNs>0) {
		
		/* memory allocation for file handles */
		PostPredFiles_pop = (FILE **)calloc(NumPredPostNs,sizeof(FILE *));
		if(RU) {
			PostPredFiles_rep = (FILE **)calloc(NumPredPostNs,sizeof(FILE *));
		}
		for(i=0;i<NumPredPostNs;i++)  { char temp[3000];
			sprintf(temp,"pop_pred_post_draws_%d.txt",PredPostNs[i]);
			PostPredFiles_pop[i] = OpenTraceFile(temp, "trace file of population allocation draws from the posterior predictive.", Baselines, BO, NULL, BO->ZSumTraceInterval);
			if(RU) {
				sprintf(temp,"rep_unit_pred_post_draws_%d.txt",PredPostNs[i]);
				PostPredFiles_rep[i] = OpenTraceFile(temp, "trace file of reporting unit allocation draws from the posterior predictive.", Baselines, BO, RU, BO->ZSumTraceInterval);
			}
		}
	}
	
	if(BO->ZSumTraceInterval>0 && gPredRan_p->NumFiles>0) {   /* open up files for output if we have the Post Pred Rans samples from the posteriors to put in there */
		
		/* memory allocation for file handles */
		PostPredRanFiles_pop = (FILE **)calloc(gPredRan_p->NumFiles, sizeof(FILE *));
		if(RU) {
			PostPredRanFiles_rep = (FILE **)calloc(gPredRan_p->NumFiles, sizeof(FILE *));
		}
		for(i=0;i<gPredRan_p->NumFiles;i++)  { char temp[3000];
			sprintf(temp,"pop_pred_post_draws_%s.txt",gPredRan_p->Suf[i]);
			PostPredRanFiles_pop[i] = OpenTraceFile(temp, "trace file of population allocation draws from the posterior predictive.", Baselines, BO, NULL, BO->ZSumTraceInterval);
			if(RU) {
				sprintf(temp,"rep_unit_pred_post_draws_%s.txt",gPredRan_p->Suf[i]);
				PostPredRanFiles_rep[i] = OpenTraceFile(temp, "trace file of reporting unit allocation draws from the posterior predictive.", Baselines, BO, RU, BO->ZSumTraceInterval);
			}
		}
	}
	
	
	
	
	/* first, we allocate memory to the Pi vector.  This is silly here, I am setting myself up
	 for a memory leak, but I am going to just allocate to it, whether it needs it or not.  We also will initialize it
	 to uniform. */
	Baselines->Pi = (popindval *)ECA_CALLOC(Baselines->NumPops,sizeof(popindval));
	init = 1.0/Baselines->NumPops;  /* initialize to uniform */
	for(i=0;i<Baselines->NumPops;i++)  {
		Pi[i] = init;
		Baselines->Pi[i].pop = i;
	}
	
	/* now, we allocate memory to the matrix of posteriors we need */
	M = (double **)ECA_CALLOC(N,sizeof(double *));
	for(i=0;i<N;i++) {
		M[i] = (double *)ECA_CALLOC(P,sizeof(double));
	}
	/* and allocate to the Z and PZ */
	Z = (int *)ECA_CALLOC(N,sizeof(int));
	
	/* now, fill the M matrix */
	for(i=0;i<N;i++) {
		ComputePostVector(Baselines, TheMixture->Inds[i], Pi, 1, 0, dummy);
		for(j=0;j<P;j++)  {
			M[i][j] = TheMixture->Inds[i]->posts[j].val;
		}
	}
	
	
	/* now, everything is in place and we just have to do our sweeps */
	for(steps=0;steps<BO->NumBurnIn + BO->NumSweeps;steps++) {  /* MUST ADD NUMSWEEPS HERE */
		
		/*printf("SIMPLE_MCMC_MIXTURE_DEBUG: PiValues  Step=%d  ",steps);
		for(j=0;j<P;j++) { 
			printf(" %f",Pi[j]);
		}
		printf("\n");*/
		
		if(steps%500==0) {
			printf("DOING_MCMC:  Step %d of %d\n",steps,BO->NumBurnIn + BO->NumSweeps);
		}
		
		
		for(j=0;j<P;j++) { 
			PZ[j] = 1.0/P; /* now, we add the prior on the mixing proportions.  Currently it is just hard wired at 1/10. */
			ZSum[j] = 0.0;  /* this is redundant, but I want to have one that doesn't include the Pi Prior in it */
		} 
		for(i=0;i<N;i++) {
			/* first, compute the posterior given current value of Pi */
			for(j=0,normo=0.0;j<P;j++) { 
				PofZ[j] = Pi[j] * M[i][j];
				normo += PofZ[j];
			} 
			for(j=0;j<P;j++) {
				PofZ[j]/=normo;
			}
			
			/*printf("SIMPLE_MCMC_MIXTURE_DEBUG: M Step=%d_Fish=%d  ",steps,i);
			for(j=0;j<P;j++) { 
				printf(" %.3f",M[i][j]);
			}
			printf("\n");
			
			
			printf("SIMPLE_MCMC_MIXTURE_DEBUG: PoZ Step=%d_Fish=%d  ",steps,i);
			for(j=0;j<P;j++) { 
				printf(" %.3f",PofZ[j]);
			}
			printf("\n");
			*/
			
			/**** RECORDING AND PRINTING PofZ VALUES, BOTH RAW AND FOR REPORTING UNITS *****/
			if(steps>BO->NumBurnIn) {
				if(RU) for(j=0;j<RU->NumRepUnits;j++)  {
					RepUnitPofZs[i][j]->v = 0.0; /* initialize to accumulate a sum */
				}
				for(j=0;j<P;j++) {
					TheMixture->Inds[i]->posts[j].DV->v = PofZ[j];
					IncrementDval(TheMixture->Inds[i]->posts[j].DV);
					if(RU) RepUnitPofZs[i][ RU->RepUnitOfPops[j] ]->v += PofZ[j];
				}
				if(RU) for(j=0;j<RU->NumRepUnits;j++)  {
					IncrementDval(RepUnitPofZs[i][j]);
				}
			}
			
			 
			
			/* now, randomly draw a Z for fish i (i.e. allocate it to a population) */
			Z[i] = IntFromProbsRV(PofZ, 0, P); /* note, we will want to store the PoZ for a Rao-Blackwellized estimator for each fish */
			
			/* and increment the appropriate element of PZ */
			PZ[Z[i]]+=1.0;
			ZSum[Z[i]]+=1;
		}
		
		/*printf("SIMPLE_MCMC_MIXTURE_DEBUG: PZ Step=%d  ",steps);
		for(j=0;j<P;j++) { 
			printf(" %.0f",PZ[j]);
		}
		printf("\n"); */
		
		/* now we just choose a new Pi */
		DirichletRV(PZ,P,Pi);
		
		
		/***** RECORDING THE Pi and the PZ VALUES (AND PRINTING TRACES) ******/
		if(steps>BO->NumBurnIn) for(j=0;j<P;j++) { 
			PopPis[j]->v = Pi[j];
			IncrementDval(PopPis[j]);
		}
		
		/* print pop traces, if indicated */
		if(BO->ZSumTraceInterval>0 && steps % BO->ZSumTraceInterval==0) {
			fprintf(popZSumtracef,"%d",steps);
			for(j=0;j<P;j++) {
				fprintf(popZSumtracef,"\t%d",ZSum[j]);
			}
			fprintf(popZSumtracef,"\n");
		}
		if(BO->PiTraceInterval>0 && steps % BO->PiTraceInterval==0) {
			fprintf(popPi_tracef,"%d",steps);
			for(j=0;j<P;j++) {
				fprintf(popPi_tracef,"\t%.8f",Pi[j]);
			}
			fprintf(popPi_tracef,"\n");
		}
		
		if(RU!=NULL) { 
			/* sum those up into Rep Units */
			for(j=0;j<RU->NumRepUnits;j++)  {
				RepUnitPis[j]->v = 0.0;  /* initialized to accumulate a sum */
				RepUnitZSum[j] = 0;
			}
			for(j=0;j<P;j++) { 
				RepUnitPis[ RU->RepUnitOfPops[j] ]->v += Pi[j];
				RepUnitZSum[ RU->RepUnitOfPops[j] ] += ZSum[j];
			}
			if(steps>BO->NumBurnIn) for(j=0;j<RU->NumRepUnits;j++)  {
				IncrementDval(RepUnitPis[j]);  /* accumulate the sum */
			}
			if(BO->ZSumTraceInterval>0 && steps % BO->ZSumTraceInterval==0) {
				fprintf(repZSumtracef,"%d",steps);
				for(j=0;j<RU->NumRepUnits;j++) {
					fprintf(repZSumtracef,"\t%d",RepUnitZSum[j]);
				}
				fprintf(repZSumtracef,"\n");
			}
			if(BO->PiTraceInterval>0 && steps % BO->PiTraceInterval==0) {
				fprintf(repPi_tracef,"%d",steps);
				for(j=0;j<RU->NumRepUnits;j++) {
					fprintf(repPi_tracef,"\t%.8f",RepUnitPis[j]->v);
				}
				fprintf(repPi_tracef,"\n");
			}
			/*printf("SIMPLE_MCMC_REP_UNIT_PI: PZ Step=%d  ",steps);
			for(j=0;j<RU->NumRepUnits;j++)  {
				printf("  %.5f",RepUnitPis[j]->v);
			}
			printf("\n"); */
		}
		
		/* and here is a little function that does the sampling of zsums from the posterior predictive for allocations, and then prints those to the files for those */
		if(BO->ZSumTraceInterval>0 && NumPredPostNs>0 && steps % BO->ZSumTraceInterval==0) {
			DoAndPrintPostPredSampling(NumPredPostNs, PredPostNs, Pi, Baselines->NumPops, RU, steps, PostPredFiles_pop, PostPredFiles_rep);
		}
		/* and here is the same for the Pred Rans Samp version */
		if(BO->ZSumTraceInterval>0 && gPredRan_p->NumFiles>0 && steps % BO->ZSumTraceInterval==0) {
			DoAndPrintPostPredRansSampling(gPredRan_p, Pi, Baselines->NumPops, RU, steps, PostPredRanFiles_pop, PostPredRanFiles_rep);
		}
	}

	/* closing all the trace files here */
	if(BO->ZSumTraceInterval>0) { 
		fclose(popZSumtracef);
		if(RU) fclose(repZSumtracef);
		if(NumPredPostNs>0) {
			for(i=0;i<NumPredPostNs;i++)  { 
				fclose(PostPredFiles_pop[i]);
				if(RU) {
					fclose(PostPredFiles_rep[i]);
				}
			}
		}
	}
	
	
	
	if(BO->PiTraceInterval>0) {
		fclose(popPi_tracef);
		if(RU )fclose(repPi_tracef);
	}
	
	
	
	
	
	
	/***** PRINTING THE AVERAGES OF THINGS *****/
	/* print out the population Pis. We will include a column that has their reporting unit if RU's are available */
	if( (outf = fopen("pop_pi_posterior_means.txt","w"))==NULL) {
		fprintf(stderr,"Failed to open file pop_pi_posterior_means.txt to write to it.  Fatal Error.  Exiting\n");
		exit(1);
	}
	fprintf(outf,"## population Pi values from gsi_sim.  NumBurnIn= %d  NumSweeps= %d  NumValuesInAverage= %d\n",BO->NumBurnIn, BO->NumSweeps, PopPis[0]->NumAved);
	if(RU) {
		fprintf(outf,"Pop\tRepUnit\tMean.Pi\tSt.Dev.of.Pis\n");  /* this is our header line if we have reporting units */
	}
	else {
		fprintf(outf,"Pop\tMean.Pi\tSt.Dev.of.Pis\n");  /* this is our header line */
	}
	for(j=0;j<P;j++)  {
		if(RU) {
			fprintf(outf,"%s\t%s\t%.8f\t%.8f\n",Baselines->PopSumms[j]->Name, RU->Names[ RU->RepUnitOfPops[j]], PopPis[j]->Ave, sqrt(PopPis[j]->Var * PopPis[j]->NumAved) );
		}
		else {
			fprintf(outf,"%s\t%.8f\t%.8f\n",Baselines->PopSumms[j]->Name, PopPis[j]->Ave, sqrt(PopPis[j]->Var * PopPis[j]->NumAved) );
		}
	}
	fclose(outf);
	
	/* then print the population Indiv PofZs */
	if( (outf = fopen("pop_pofz_posterior_means.txt","w"))==NULL) {
		fprintf(stderr,"Failed to open file pop_pofz_posterior_means.txt to write to it.  Fatal Error.  Exiting\n");
		exit(1);
	}
	fprintf(outf,"## population PofZ values from gsi_sim.  NumBurnIn= %d  NumSweeps= %d  NumValuesInAverage= %d\n",BO->NumBurnIn, BO->NumSweeps, TheMixture->Inds[0]->posts[0].DV->NumAved);
	fprintf(outf,"IndivName\tValueType");  /* this is our header line */
	for(j=0;j<P;j++)  {
		fprintf(outf,"\t%s",Baselines->PopSumms[j]->Name);
	}
	fprintf(outf,"\n");
	for(i=0;i<N;i++) {
		fprintf(outf,"%s\tPosterior.Mean",TheMixture->Inds[i]->Name);
		for(j=0;j<P;j++)  {
			fprintf(outf,"\t%.8f",TheMixture->Inds[i]->posts[j].DV->Ave);
		}
		fprintf(outf,"\n");
		fprintf(outf,"%s\tSD",TheMixture->Inds[i]->Name);
		for(j=0;j<P;j++)  {
			fprintf(outf,"\t%.8f",sqrt(TheMixture->Inds[i]->posts[j].DV->Var * TheMixture->Inds[i]->posts[j].DV->NumAved));
		}
		fprintf(outf,"\n");
	}
	fclose(outf);
	
	
	/* here print out all the reporting unit stuff */
	if(RU) {
		
		/* first print the reporting unit Pis */
		if( (outf = fopen("rep_unit_pi_posterior_means.txt","w"))==NULL) {
			fprintf(stderr,"Failed to open file rep_unit_pi_posterior_means.txt to write to it.  Fatal Error.  Exiting\n");
			exit(1);
		}
		fprintf(outf,"## reporting unit Pi values from gsi_sim.  NumBurnIn= %d  NumSweeps= %d  NumValuesInAverage= %d\n",BO->NumBurnIn, BO->NumSweeps, RepUnitPis[0]->NumAved);
		fprintf(outf,"RepUnit\tMean.Pi\tSt.Dev.of.Pis\n");  /* this is our header line */
		for(j=0;j<RU->NumRepUnits;j++)  {
			fprintf(outf,"%s\t%.8f\t%.8f\n",RU->Names[j], RepUnitPis[j]->Ave, sqrt(RepUnitPis[j]->Var * RepUnitPis[j]->NumAved) );
		}
		fclose(outf);
		
		
		
		/* then print the reporting unit Indiv PofZs */
		if( (outf = fopen("rep_unit_pofz_posterior_means.txt","w"))==NULL) {
			fprintf(stderr,"Failed to open file rep_unit_pofz_posterior_means.txt to write to it.  Fatal Error.  Exiting\n");
			exit(1);
		}
		fprintf(outf,"## reporting unit PofZ values from gsi_sim.  NumBurnIn= %d  NumSweeps= %d  NumValuesInAverage= %d\n",BO->NumBurnIn, BO->NumSweeps, RepUnitPofZs[0][0]->NumAved);
		fprintf(outf,"IndivName\tValueType");  /* this is our header line */
		for(j=0;j<RU->NumRepUnits;j++)  {
			fprintf(outf,"\t%s",RU->Names[j]);
		}
		fprintf(outf,"\n");
		for(i=0;i<N;i++) {
			fprintf(outf,"%s\tPosterior.Mean",TheMixture->Inds[i]->Name);
			for(j=0;j<RU->NumRepUnits;j++)  {
				fprintf(outf,"\t%.8f",RepUnitPofZs[i][j]->Ave);
			}
			fprintf(outf,"\n");
			fprintf(outf,"%s\tSD",TheMixture->Inds[i]->Name);
			for(j=0;j<RU->NumRepUnits;j++)  {
				fprintf(outf,"\t%.8f",sqrt(RepUnitPofZs[i][j]->Var * RepUnitPofZs[i][j]->NumAved));
			}
			fprintf(outf,"\n");
		}
		fclose(outf);
		
	}
	
	/* and for a final hoorah, we sort the popindvals on the basis of posterior mean so that we can use that
	 sorted array for the LoglIndSims */
	/* first copy the posterior mean into the val field in the popindval */
	for(i=0;i<N;i++) {
		for(j=0;j<P;j++)  {
			TheMixture->Inds[i]->posts[j].val = TheMixture->Inds[i]->posts[j].DV->Ave;
		}
	}
	/* then cycle over individuals and sort by that val field */
	for(i=0;i<N;i++) {
		qsort( (void *)(TheMixture->Inds[i]->posts), (size_t)(Baselines->NumPops), sizeof(popindval), CmpPopIndByVal);
	}
	
	/* now, we should be ready to do some LoglIndSims with them */
	
	/* must do more memory freeing here */
	free(Pi);
}




/* this function simply computes GSI-related output (only Rannala-Mountain implemented so far)
for TheMixture given reference samples of Baselines.  Note that the data have to be in the
structures, all proper and ready to go.  This function also sorts the Posterior probs into 
descending order.  If Sort is 0, then the posterior vector of each individual is NOT sorted
from high to low posterior probability.  Otherwise it is. */
void ComputeGSISample(struct IndCollection *Baselines, struct IndCollection *TheMixture, em_opts *EMOpts, int Sort)
{
	int j,i,steps,k;
	double init;
	double *PiSum = (double *)ECA_CALLOC(Baselines->NumPops,sizeof(double));
	double *PiTemp = (double *)ECA_CALLOC(Baselines->NumPops,sizeof(double));
	double **M; /* matrix indexed by mix indiv and then by pop */
	double dist;
	int dummy=0;
	int N=TheMixture->NumInds;  /* number of fish in the mixture */
	int P=Baselines->NumPops;  /* number of pops in baseline */
	double *PofZ = (double *)ECA_CALLOC(Baselines->NumPops,sizeof(double)); /* variable to be used temporarily for each fish in the mixture */
	double PiNormo;
	
	/* first, we allocate memory to the Pi vector.  This is silly here, I am setting myself up
	for a memory leak, but I am going to just allocate to it, whether it needs it or not.  We also will initialize it
	to uniform. */
	Baselines->Pi = (popindval *)ECA_CALLOC(Baselines->NumPops,sizeof(popindval));
	init = 1.0/Baselines->NumPops;
	for(i=0;i<Baselines->NumPops;i++)  {
		PiTemp[i] = init;
		Baselines->Pi[i].pop = i;
	}
	
	
	/* for some crazy reason I was recomputing the likelihood each iteration through (I don't know why...
	 I was pretty certain in the very first iteration I precomputed and stored them).  Whatever.  I do that
	 here now */
	/* now, we allocate memory to the matrix of posteriors we need */
	M = (double **)ECA_CALLOC(N,sizeof(double *));
	for(i=0;i<N;i++) {
		M[i] = (double *)ECA_CALLOC(P,sizeof(double));
	}
	/* now, fill the M matrix */
	for(i=0;i<N;i++) {
		ComputePostVector(Baselines, TheMixture->Inds[i], PiTemp, 1, 0, dummy);
		for(j=0;j<P;j++)  {
			M[i][j] = TheMixture->Inds[i]->posts[j].val;
		}
	}
	
	/* now, the matrix M contains posteriors given a uniform prior, so they can 
	 behave like scaled likelihoods */
	for(steps=0;steps<EMOpts->MaxSteps;steps++) {
		for(j=0;j<Baselines->NumPops;j++) {
			PiSum[j] = 0.0;
		} 
		for(i=0;i<TheMixture->NumInds;i++) {
			PiNormo = 0.0;
			for(j=0;j<Baselines->NumPops;j++) { double temp;
				temp = M[i][j] * PiTemp[j];
				PofZ[j] = temp;  /* this is the part where we are re-estimating Pi */
				PiNormo += temp;
			}
			for(j=0;j<Baselines->NumPops;j++) {
				PiSum[j] += PofZ[j] / PiNormo;  /* normalize the contribution of each fish to sum to one */
			}
			
		}
		/* here we normalize PiSum and compute distance from previous value */
		dist = 0.0;
		for(j=0;j<Baselines->NumPops;j++) {
			PiSum[j] /= (double)N;
			dist += ECA_ABS(PiSum[j] - PiTemp[j]);
			PiTemp[j] = PiSum[j];
			/*printf("EM_STEP %d: %s  %f\n",steps,Baselines->PopNames[j],PiTemp[j]);*/
		}
		printf("EM_ALGORITHM_PROGRESS: step= %d  difference= %.20f : ",steps+1,dist);
		for(k=0;k<Baselines->NumPops;k++)  {
			printf(" %f",PiSum[k]);
		}
		printf("\n");
		if(dist<EMOpts->tol && steps>=EMOpts->MinSteps-1)  {/* here we can set a tolerance in EMOpts */
			printf("EM_COMPLETION_AT_TOL:  step=%d  dist=%.20f\n",steps+1,dist);
			break;
		}
		if(steps==EMOpts->MaxSteps-1) {
			printf("EM_COMPLETION_AT_MAX_STEP:  step=%d  dist=%.20f\n",steps+1,dist);
		}
	}
	
	/* here we copy those into the Pi field in Baselines */
	for(j=0;j<Baselines->NumPops;j++) {
		Baselines->Pi[j].val = PiTemp[j];
	}
	
	/* and here we compute the posterior probs for each individual using M[i][j] and weighting by what we have in PiTemp[j] and we put the normalized result
	 back over to Inds[i]->posts[j].val so that we can print them.  This is just the plug-in estimator using Bayes Theorem.  */
	for(i=0;i<TheMixture->NumInds;i++) { double normo;
		normo = 0.0;
		for(j=0;j<Baselines->NumPops;j++)  {
			TheMixture->Inds[i]->posts[j].val = M[i][j] * PiTemp[j];
			normo += TheMixture->Inds[i]->posts[j].val;
		}
		for(j=0;j<Baselines->NumPops;j++)  {
			TheMixture->Inds[i]->posts[j].val /= normo;
		}
	}
	
	/* then, cycle over the individuals one last time and sort them according to their posterior probs */
	if(Sort) {
			for(i=0;i<TheMixture->NumInds;i++) {
				qsort( (void *)(TheMixture->Inds[i]->posts), (size_t)(Baselines->NumPops), sizeof(popindval), CmpPopIndByVal);
			}
	} 
	
	free(PiSum);
	free(PiTemp);
}

/* given the population summaries in NumPops pop structs in P, use the xDP variable as the parameter of
a Dirichlet distribution to simulate new values for the "simp" field of each locus in each pop struct. 
*/
void RealizeNewSimP(int NumPops, pop_struct **PopSumms)
{
	int i,j;
	pop_struct *P;  /* for a shorthand */
	
	for(i=0;i<NumPops;i++)  {
		P = PopSumms[i];
		
		for(j=0;j<P->NumLoc;j++)  {
			DirichletRV(P->xDP[j], P->NumAlle[j], P->simp[j]);
		}
	}
}

/* this function simply computes the GeneClass output (only Rannala-Mountain implemented so far)
for TheMixture given reference samples of Baselines.  Note that the data have to be in the
structures, all proper and ready to go.  This function also sorts the Posterior probs into 
descending order. 

CorrectWay is the number of Monte Carlo reps to use if we want to sample over allele freqs from their posterior dsn given the 
reference samples (i.e. we wish to do it correctly).  If it is <=zero we just do Rannala-Mountain on it. 
 
 If Sort is zero, the populations are not sorted according to posterior prob for each indiv.  Otherwise, they are. */
void ComputeGeneClassSample(struct IndCollection *Baselines, struct IndCollection *TheMixture, gc_opts *GeneClassOpts, int CorrectWayReps, int Sort)
{
	int i,k,j;
	
	if(CorrectWayReps<=0) {
		for(i=0;i<TheMixture->NumInds;i++) {
			ComputePostVector(Baselines, TheMixture->Inds[i], NULL, 1, 0, CorrectWayReps);
			if(Sort) {
				qsort( (void *)(TheMixture->Inds[i]->posts), (size_t)(Baselines->NumPops), sizeof(popindval), CmpPopIndByVal);
			}
		}
	}
	else { /* if doing the proper-bayes thing, we must first do some initialization stuff, then cycle over the CorrectWayReps, and for each rep we must:
	             1. Simulate a new value for simp given the baselines 
				 2. Cycle over indivs in the mixture and compute the post vector for each one.
				 3. Add the value of the post vectors to a running average.
		   At the end, we can stick the averages back into TheMixture-->Inds[i] for all i so it can be sorted and printed.
		   At some point we may wish to output the other parts (like the Monte Carlo standard deviation, etc.  */
		
		/* first, initialize all the dvals to zero to accumulate some averages and deal with allocation if necessary (doing that here is a major Kluge.  Oh well. */
		for(i=0;i<TheMixture->NumInds;i++) {
			/* deal with allocation if necessary */
			if(TheMixture->Inds[i]->posts==NULL) {
				TheMixture->Inds[i]->posts = (popindval *)ECA_CALLOC(Baselines->NumPops,sizeof(popindval));
				/* and allocate to the Dval in that struct too */
				for(j=0;j<Baselines->NumPops;j++)  {
					TheMixture->Inds[i]->posts[j].DV = AllocDval(0,1,.01);
				}
			}
			/* even if it doesn't need allocation, make sure you initialize each dval to zero */
			for(j=0;j<Baselines->NumPops;j++)  {
				InitDvalToZero(TheMixture->Inds[i]->posts[j].DV);
			}

		}
		
		for(k=0;k<CorrectWayReps;k++) {
		
			/* simulate new values for the simp's */
			RealizeNewSimP(Baselines->NumPops, Baselines->PopSumms);
			
			for(i=0;i<TheMixture->NumInds;i++) {
				/* compute this individual's posterior prob vector given the current simp's */
				ComputePostVector(Baselines, TheMixture->Inds[i], NULL, 1, 0, CorrectWayReps);
				
				/* then put those values in the "v" field of the individuals' posts dval and increment! */
				for(j=0;j<Baselines->NumPops;j++)  {
					TheMixture->Inds[i]->posts[j].DV->v = TheMixture->Inds[i]->posts[j].val;
					IncrementDval(TheMixture->Inds[i]->posts[j].DV);
				}
			}
		}
		
		/* finally, we have to transfer the averages back to posts[j].val so that they can be sorted and printed. Here we also sort them! */
		for(i=0;i<TheMixture->NumInds;i++) {
			for(j=0;j<Baselines->NumPops;j++)  {
				TheMixture->Inds[i]->posts[j].val = TheMixture->Inds[i]->posts[j].DV->Ave;
			}
			if(Sort) {
				qsort( (void *)(TheMixture->Inds[i]->posts), (size_t)(Baselines->NumPops), sizeof(popindval), CmpPopIndByVal);
			}
		}
		
	} 
}


/* this function simply computes the GeneClass output (only Rannala-Mountain implemented so far)
for self-assignment of Baselines.  Note that the data have to be in the
structures, all proper and ready to go.  This function also sorts the Posterior probs into 
descending order.  */
void ComputeGeneClassSelfAssign(struct IndCollection *Baselines, gc_opts *GeneClassOpts)
{
	int i;
	int dummy=0;
	
	for(i=0;i<Baselines->NumInds;i++) {
		ComputePostVector(Baselines, Baselines->Inds[i], NULL, 1, 1, dummy);
		qsort( (void *)(Baselines->Inds[i]->posts), (size_t)(Baselines->NumPops), sizeof(popindval), CmpPopIndByVal);
	}
}



/* this function simply computes the GeneClass output (only Rannala-Mountain implemented so far)
for self-assignment of Baselines.  Note that the data have to be in the
structures, all proper and ready to go.  This version doesn't sort stuff into descending order  */
void ComputeGeneClassSelfAssign_NOSORT(struct IndCollection *Baselines, gc_opts *GeneClassOpts)
{
	int i;
	int dummy = 0;
	
	for(i=0;i<Baselines->NumInds;i++) {
		ComputePostVector(Baselines, Baselines->Inds[i], NULL, 1, 1, dummy);
	}
}



/* quick function to print the values in posts[j].val after doing EM or gene-class or OneStep.  This just prints them in whatever
 order they are in.  So you have to be careful not to sort them by value first!  They should be sorted by population index!
 The purpose of this is to create output files that can easily be read by R and are in the same format as those that get printed
 after the Bayesian analysis.  If RU is not NULL then it prints out the same for the reporting units into another file named
 as given in RepFileName
*/
void FprintIndivPosteriorFields(struct IndCollection *Baselines, struct IndCollection *TheMixture, const char *PopFileName, reporting_unit_info *RU, const char *RepFileName)
{
	int i,j;
	FILE *outf;
	int N=TheMixture->NumInds;  /* number of fish in the mixture */
	int P=Baselines->NumPops;  /* number of pops in baseline */
	int NumLocUsed;
	double RepUnitPosts[MAX_POPS];
	
	if( (outf = fopen(PopFileName,"w"))==NULL) {
		fprintf(stderr,"Failed to open file %s to write to it.  Fatal Error.  Exiting\n",PopFileName);
		exit(1);
	}
	
	fprintf(outf,"## population PofZ values from gsi_sim.  These are plug-in values of type indicated by the file name\n");
	fprintf(outf,"IndivName\tNumLoc");  /* this is our header line */
	for(j=0;j<P;j++)  {
		fprintf(outf,"\t%s",Baselines->PopSumms[j]->Name);
	}
	fprintf(outf,"\n");
	for(i=0;i<N;i++) {
		
		/* now, we count how many loci were used */
		NumLocUsed = 0;
		for(j=0;j<Baselines->NumLoc;j++)  {
			if(TheMixture->Inds[i]->y[j][0]>=0 && TheMixture->Inds[i]->y[j][1]>=0) {
				NumLocUsed++;
			}
		}
		
		fprintf(outf,"%s\t%d",TheMixture->Inds[i]->Name,NumLocUsed);
		for(j=0;j<P;j++)  {
			fprintf(outf,"\t%.8f",TheMixture->Inds[i]->posts[j].val);
		}
		fprintf(outf,"\n");
	}
	fclose(outf);
	
	
	if(RU)  {
		if( (outf = fopen(RepFileName,"w"))==NULL) {
			fprintf(stderr,"Failed to open file %s to write to it.  Fatal Error.  Exiting\n",RepFileName);
			exit(1);
		}
		fprintf(outf,"## reporting unit PofZ values from gsi_sim.  These are plug-in values of type indicated by the file name\n");
		fprintf(outf,"IndivName\tLocNum");  /* this is our header line */
		for(j=0;j<RU->NumRepUnits;j++)  {
			fprintf(outf,"\t%s",RU->Names[j]);
		}
		fprintf(outf,"\n");
		for(i=0;i<N;i++) {
			
			/* now, we count how many loci were used */
			NumLocUsed = 0;
			for(j=0;j<Baselines->NumLoc;j++)  {
				if(TheMixture->Inds[i]->y[j][0]>=0 && TheMixture->Inds[i]->y[j][1]>=0) {
					NumLocUsed++;
				}
			}
			for(j=0;j<RU->NumRepUnits;j++)  {
				RepUnitPosts[j]=0.0; /* to accumulate a sum */
			}
			for(j=0;j<P;j++)  {
				RepUnitPosts[ RU->RepUnitOfPops[j] ] += TheMixture->Inds[i]->posts[j].val;
			}
			
			fprintf(outf,"%s\t%d",TheMixture->Inds[i]->Name,NumLocUsed);
			for(j=0;j<RU->NumRepUnits;j++)  {
				fprintf(outf,"\t%.8f",RepUnitPosts[j]);
			}
			fprintf(outf,"\n");
		}
		fclose(outf);		
	}
}
 

/* Here is a quick function to open a file and print to it the Pi values
 that are in the Baselines struct.  This is useful for printing the MLE's.
 
 It also sums values up into reporting groups.
 
*/
void FprintPis(struct IndCollection *Baselines, const char *PopFileName, reporting_unit_info *RU, const char *RepFileName) 
{
	int j;
	FILE *outf;
	int P=Baselines->NumPops;  /* number of pops in baseline */
	double RepUnitMLEs[MAX_POPS];
	
	
	/* do the populations ones */
	if( (outf = fopen(PopFileName,"w"))==NULL) {
		fprintf(stderr,"Failed to open file %s to write to it.  Fatal Error.  Exiting\n",PopFileName);
		exit(1);
	}
	fprintf(outf,"## population Pi values from gsi_sim.  Likelihood estimates of type indicated by file name: %s\n", PopFileName);
	if(RU) {
		fprintf(outf,"Pop\tRepUnit\tEst.Pi\n");  /* this is our header line if we have reporting units */
	}
	else {
		fprintf(outf,"Pop\tEst.Pi\n");  /* this is our header line */
	}
	for(j=0;j<P;j++)  {
		if(RU) {
			fprintf(outf,"%s\t%s\t%.8f\n",Baselines->PopSumms[j]->Name, RU->Names[RU->RepUnitOfPops[j]], Baselines->Pi[j].val);
		}
		else {
			fprintf(outf,"%s\t%.8f\n",Baselines->PopSumms[j]->Name,Baselines->Pi[j].val);
		}
	}
	fclose(outf);
	
	
	
	
	/* do the rep-unit ones if indicated */
	if(RU) {
		if( (outf = fopen(RepFileName,"w"))==NULL) {
			fprintf(stderr,"Failed to open file %s to write to it.  Fatal Error.  Exiting\n",RepFileName);
			exit(1);
		}
		fprintf(outf,"## population Pi values from gsi_sim.  Likelihood estimates of type indicated by file name: %s\n", RepFileName);
		fprintf(outf,"RepUnit\tEst.Pi\n");  /* this is our header line if we have reporting units */
	
		for(j=0;j<RU->NumRepUnits;j++)  {
			RepUnitMLEs[j] = 0.0;
		}
		for(j=0;j<P;j++)  {
			RepUnitMLEs[RU->RepUnitOfPops[j]]+=Baselines->Pi[j].val;
		}
		for(j=0;j<RU->NumRepUnits;j++)  {
			fprintf(outf,"%s\t%.8f\n",RU->Names[j],RepUnitMLEs[j]);
		}
		fclose(outf);
	}
}



/* incomplete function.  currently it just prints out the posterior probs in S (the mixture), that correspond to
the populations in R (the baselines).  It prints out the top X posterior probs.  You can use Tag to specify how 
the lines should be tagged. Include the colon. */
void PrintGeneClassCSV(struct IndCollection *R, struct IndCollection *S,  gc_opts *GeneClassOpts, const char *Tag) 
{
	int i,j;
	int NumLocUsed;
	
	for(i=0;i<S->NumInds;i++)  {
		printf("%s/%s",Tag,S->Inds[i]->Name);
		for(j=0;j<R->NumPops;j++)  {
			printf(";%s;%.3f;",R->PopNames[S->Inds[i]->posts[j].pop],100.0 * S->Inds[i]->posts[j].val);
		}
		/* here we include the logls */
		for(j=0;j<R->NumPops;j++)  {
			printf(";%.3f",-S->Inds[i]->logls[j].val);
		}
		/* now, we count how many loci were used */
		NumLocUsed = 0;
		for(j=0;j<S->NumLoc;j++)  {
			if(S->Inds[i]->y[j][0]>=0 && S->Inds[i]->y[j][0]>=0) {
				NumLocUsed++;
			}
		}
		printf(";%d ;",NumLocUsed);
		/* now print the names of the ones used */
		for(j=0;j<S->NumLoc;j++)  {
			if(S->Inds[i]->y[j][0]>=0 && S->Inds[i]->y[j][0]>=0) {
				printf(" %s ",S->AlleInfo->LocusNames[j]);
			}
		}
		printf(" ;");
		/* then print the ones not used */
		for(j=0;j<S->NumLoc;j++)  {
			if( !(S->Inds[i]->y[j][0]>=0 && S->Inds[i]->y[j][0]>=0) ) {
				printf(" %s ",S->AlleInfo->LocusNames[j]);
			}
		}
		printf("\n");
	}
}


/* This is a simple function to compute the probability of an indiv's genotype given simp, the allele freqs, using a simple binomial model
and only for disomic critters */
long double ProbIndFromPop_P_SingleLocusDisomic(struct IndCollection *R, ind_struct *T, int p, int LeaveOut, int j)
{
	int *y;  /* shorthand for the genotype */
	double *simp;  /* temporary, shorthand to get to the alle freq vec */
	long double ret = 1.0;
	
	/* if not doing leave out */
	if(LeaveOut==0) {
		y = T->y[j];
		if( (y[0] >= 0) && (y[1] >= 0) ) {  /* if the data aren't missing */
			simp = R->PopSumms[p]->simp[j];
			
			if(y[0]==y[1]) {  /* if it is a homozygote */
				ret = simp[y[0]] * simp[y[0]];  /* just the frequency of the homozygote */
			}
			else {  /* otherwise it is a heterozygote */
				ret = 2.0 * simp[y[0]] * simp[y[1]];
			}
		}
	}
	
	else {
		fprintf(stderr, "Error! Currently leave-one-out is not available in ProbIndFromPop_P_SingleLocusDisomic\n");
		exit(1);
	}

	return(ret);
	
}




/* this just returns an indiv's genotype prob given the current value of simp */
long double ProbIndFromPop_P_Disomic(struct IndCollection *R, ind_struct *T, int p, int LeaveOut)
{
	int j;
	double ret=1.0; /* set to 1 to accumulate a product */
	
	/* if not doing leave out */
	for(j=0;j<R->NumLoc;j++)  {/* cycle over loci */
		if(R->AlleInfo->Ploidy[j]==2) {
				 ret *= ProbIndFromPop_P_SingleLocusDisomic(R,T,p,LeaveOut,j);
		}
		else {
			fprintf(stderr,"Sorry! Can't do Ploidy>2 with ProbIndFromPop_P_Disomic(). Exiting\n");
			exit(1);
		}
	}
	
	return(ret);
}



/* This returns the probability of indiv T's genotype given it came from population p.
It gets called with the whole struct IndCollection so all the other necessary variables are accessible. 
It uses the Rannala-Mountain
method using the xDP field of each pop_struct.  It requires that xDP sum be appropriately computed ahead of time.

If flag LeaveOut != 0, then the individual's genotype is subtracted off the xDP value ahead of time.

This calculation is based on the equation 35.151 repentation of the CDM in Johnson, Kotz, and Balakrishnan.

This version is now changed to a single locus.  I have written a new version that is relevant to fully penetrant 
polysomic loci.  But then I modified this one to more quickly do disomic loci. 

 */
long double ProbIndFromPopR_M_SingleLocusDisomic(struct IndCollection *R, ind_struct *T, int p, int LeaveOut, int j, long double *OutNumer)
{
	int *y;  /* shorthand for the genotype */
	double *a;  /* temporary, shorthand to get to the dirichlet pars */
	long double numer = 1.0, denom=1.0;
	
	/* if not doing leave out */
	if(LeaveOut==0) {
		y = T->y[j];
		if( (y[0] >= 0) && (y[1] >= 0) ) {  /* if the data aren't missing */
			a = R->PopSumms[p]->xDP[j];
			
			numer *= 2.0;
			denom *= R->PopSumms[p]->xDPsum[j] * (R->PopSumms[p]->xDPsum[j] + 1.0);
			
			if(y[0]==y[1]) {  /* if it is a homozygote */
				numer *= a[y[0]] * (a[y[0]] + 1.0);
				denom *= 2.0;
			}
			else {  /* otherwise it is a heterozygote */
				numer *= a[y[0]] * a[y[1]];
			}
		}
	}
	else {
		y = T->y[j];
		if( (y[0] >= 0) && (y[1] >= 0) ) {  /* if the data aren't missing */
			a = R->PopSumms[p]->xDP[j];
			
			numer *= 2.0;
			/* when doing LeaveOneOut, xDPsum will be 2 less than as recorded */
			denom *= (R->PopSumms[p]->xDPsum[j]-2.0) * (R->PopSumms[p]->xDPsum[j] - 1.0);
			
			if(y[0]==y[1]) {  /* if it is a homozygote */
				numer *= (a[y[0]]-2.0) * (a[y[0]] - 1.0); /* here we subtract 2 off of the single a */
				denom *= 2.0;
			}
			else {  /* otherwise it is a heterozygote */
				numer *= (a[y[0]]-1.0) * (a[y[1]]-1.0);  /* here we subtract 1 off of each of the two a's */
			}
		}
	}
	
	*OutNumer = numer;
	return(denom);
}

/* 
Note, I had originally set this up so that the numerator and the denominator were computed 
separately, and then divided at the end (thus speeding thing up becuase it eliminated a lot of
floating point divides) but that leads to under/over-flow with lots of loci, so I had
to rework that.

*/
long double ProbIndFromPopR_M_Polysomic(struct IndCollection *R, ind_struct *T, int p, int LeaveOut)
{
	int i,j,k;
	int **C;  /* shorthand for the genotype, compactified */
	double *a;  /* temporary, shorthand to get to the dirichlet pars */
	double adot,ai;  /* temporary shorthand for the sum of the dirichlet pars */
	double A; /* temporary shorthand for the number of alleles appearing */
	long double denom,numer,result=1.0;
	int Pl, Xi; /* the Ploidy of each locus */
	double Pld; /* ploidy as a double */
	
	/* if not doing leave out */
	for(j=0;j<R->NumLoc;j++)  {/* cycle over loci */
		A = T->AllelesAppearing[j];
		if(A>0) {  /* if the data aren't missing */
			C = T->CompactRep[j];
			Pl = R->AlleInfo->Ploidy[j];
			
			/* if it is a diploid locus, then do it a little faster */
			if(Pl==2) {
				denom = ProbIndFromPopR_M_SingleLocusDisomic(R,T,p,LeaveOut,j,&numer);
				result *= numer/denom;
			}
			else {
				numer=1.0;
				denom=1.0;
				Pld = (double)Pl;
				a = R->PopSumms[p]->xDP[j];
				adot = R->PopSumms[p]->xDPsum[j];
				
				if(LeaveOut) {
					adot -= Pld;  /* if LeaveOut, then these alleles don't appear in the xsum */
				}
				
				/* get the "normalizing constant" part of the CDM */ 
				numer *= Factorial(Pl);
				for(k=Pl-1;k>=0;k--)  {
					denom *= ((double)k + adot);  /* this is faster than computing a gamma function for small Pl */
				}
				
				/* now, get the terms for all alleles present */
				for(i=0;i<A;i++)  {
					Xi = C[i][1];
					ai = a[C[i][0]];
					if(LeaveOut) {
						ai -= (double)Xi;
					}
					for(k=Xi-1;k>=0;k--) {
						numer *= (double)k + ai;
					}
					denom *= Factorial(Xi);
				}
				result *= numer/denom;
			}
		}
	}

	
	return(result);
}





/* this function fills an individuals posts array for all the populations in the struct IndCollection.
   If LOGLS is not zero, then it also fills the logls array.  It takes mixing poportion prior Pi.
   If Pi==NULL, then uniform mixing priors are assumed.  If posts or logls are NULL, then it allocates
   memory to them, as necessary. 
   
   If TryToLeaveOut is not zero, then this function will assume that the population index carried by T refers
   to the indices of the populations in R, and it will leave the individual out as appropriate.  Otherwise not.  
   
	
*/
void ComputePostVector(struct IndCollection *R, ind_struct *T, double *Pi, int LOGLS, int TryToLeaveOut, int CorrectWayReps) 
{
	int i;
	double normo=0.0;
	int LeaveOut;
	
	/* deal with allocation if necessary */
	if(T->posts==NULL) {
		T->posts = (popindval *)ECA_CALLOC(R->NumPops,sizeof(popindval));
		/* and allocate to the Dval in that struct too */
		for(i=0;i<R->NumPops;i++)  {
			T->posts[i].DV = AllocDval(0,1,-1); /* don't bother with the histograms on these, hence the -1 for the d parameter */
		}
	}
	if(LOGLS && T->logls==NULL) {
		T->logls = (popindval *)ECA_CALLOC(R->NumPops,sizeof(popindval));
	}
	
	/* then, cycle over populations */
	for(i=0;i<R->NumPops;i++)  {
		LeaveOut = TryToLeaveOut && T->PopNum==i;
		if(CorrectWayReps==0) {
			T->posts[i].val =  Pi ?  Pi[i] *  ProbIndFromPopR_M_Polysomic(R,T,i, LeaveOut)  : ProbIndFromPopR_M_Polysomic(R,T,i,LeaveOut);
		}
		else {
			if(Pi!=NULL) {
				fprintf(stderr,"Sorry! Can't do --proper-bayes > 0 when Pi is non null in ComputePostVector. Exiting\n");
				exit(1);
			}
			if(TryToLeaveOut) {
				fprintf(stderr,"Sorry! Currently can't do --proper-bayes > 0 and also TryToLeaveOut in ComputePostVector. Exiting\n");
				exit(1);
			}
			T->posts[i].val = ProbIndFromPop_P_Disomic(R,T,i,0);
		}
		T->posts[i].pop = i;
		normo += T->posts[i].val;
		if(LOGLS) {
			T->logls[i].val =  Pi ? log10(T->posts[i].val/Pi[i])   :   log10(T->posts[i].val);  /* Here, I  divide Pi[i] out of it if Pi is not NULL, so I can get the raw log likelihood! */
			T->logls[i].pop = i;
		}
	}
	
	/* in the end, we normalize all the posts */
	for(i=0;i<R->NumPops;i++)  {
		T->posts[i].val /= normo;
	}
		
}

int CmpPopIndByVal(const void *a, const void *b)
{
	popindval aa = *(popindval *)a;
	popindval bb = *(popindval *)b;
	
	if(aa.val > bb.val)
		return(-1);
	else if(aa.val < bb.val) 
		return(1);
	
	return(0);
}

int CmpPopIndByPop(const void *a, const void *b)
{
	popindval aa = *(popindval *)a;
	popindval bb = *(popindval *)b;
	
	if(aa.pop > bb.pop)
		return(1);
	else if(aa.pop < bb.pop) 
		return(-1);
	
	return(0);
}



struct IndCollection *InitIndCollection(void) 
{
	struct IndCollection *Baselines;
	
	Baselines=(struct IndCollection *)ECA_MALLOC(sizeof(struct IndCollection));
	Baselines->AlleInfo = (alle_info *)ECA_MALLOC(sizeof(alle_info));
	Baselines->AlleInfo->LocHash = NULL;
	Baselines->AlleInfo->AlleNames = NULL;
	Baselines->AlleInfo->NumAlle = NULL;
	Baselines->PopFactors=NULL;
	Baselines->PopScalePriors=NULL;
	Baselines->PopPriors=NULL;
	Baselines->SimPi = NULL;
	
	Baselines->PopStarts=NULL;
	Baselines->PopEnds=NULL;
	Baselines->NumInPop=NULL;
	Baselines->PopSumms=NULL;

	return(Baselines);
}

/* allocate stuff to a pop_struct.  For now we aren't doing any of the constell stuff */
pop_struct *AllocPopStruct(int NumLoc, int *NumAlle, int NumPops)
{
	pop_struct *ret;
	int j;
	
	ret = (pop_struct *)ECA_MALLOC(sizeof(pop_struct));
	
	ret->NumAlle = (int *)ECA_CALLOC(NumLoc, sizeof(int));
	ret->Ploidy = (int *)ECA_CALLOC(NumLoc,sizeof(int));
	
	ret->x = (double **)ECA_CALLOC(NumLoc,sizeof(double *));
	ret->xDP = (double **)ECA_CALLOC(NumLoc,sizeof(double *));
	ret->xBP = (double **)ECA_CALLOC(NumLoc,sizeof(double *));
	ret->p = (double **)ECA_CALLOC(NumLoc,sizeof(double *));
	ret->simp = (double **)ECA_CALLOC(NumLoc,sizeof(double *));

	ret->xsum = (double *)ECA_CALLOC(NumLoc,sizeof(double));	
	ret->xDPsum = (double *)ECA_CALLOC(NumLoc,sizeof(double));
	ret->Loss = NULL; /*(double *)ECA_CALLOC(NumPops,sizeof(double *));*/
	ret->ExpScaledLike = (double *)ECA_CALLOC(NumPops,sizeof(double *));
	
	for(j=0;j<NumLoc;j++)  {
		ret->x[j] = (double *)ECA_CALLOC(NumAlle[j],sizeof(double));
		ret->xDP[j] = (double *)ECA_CALLOC(NumAlle[j],sizeof(double));
		ret->xBP[j] = (double *)ECA_CALLOC(NumAlle[j],sizeof(double));
		ret->p[j] = (double *)ECA_CALLOC(NumAlle[j],sizeof(double));
		ret->simp[j] = (double *)ECA_CALLOC(NumAlle[j],sizeof(double));
	}
	
	
	/* here we record information about the original individuals that may have been 
	summarized into this pop struct.  These two variables get set in IndColl2PopStructs as appropriate.  */
	ret->OriginalIndColl = NULL;  /* if the pop_struct came from allele freqs and not inds, then this remains null */
	ret->IdxOfPopInIndColl = -1;  /* default value.  */
	ret->IdxInBaseline = -1;
	ret->RepUnitIdx = -1;
	
	return(ret);
}



/*
	If you use the "-n name -f factor -s -p p --end-pop -b Baseline.txt" type of specification 
	for the baseline allele freqs, you have to get that Pop Prior information into the struct IndCollection
	priors before making pop_structs out of the struct IndCollection.
	
	 Num is the number of pops specified with --end-pop.  
	 Pops is the array of pop_structs that came out of getting the options.
	  
*/
void UpdateBaselinePopPriors(struct IndCollection *R, int Num, pop_struct **Pops) 
{
	int i,k,idx;
	int Throw=0;
	
	for(i=0;i<Num;i++)  {
		/* do a really slow thing to find the name---would be better to have a hash table, but we only have to
		   do this once */
		idx = -1;
		for(k=0;k<R->NumPops;k++)  {
			if(strcmp(Pops[i]->Name,R->PopNames[k])==0) {
				idx = k;
				break;
			}
		}
		if(idx<0) {
			fprintf(stderr,"Error! Searched for name %s in baseline and didn't find it.  This will cause a failure...\n", Pops[i]->Name);
			Throw++;
		}
		else {
			printf("SETTING_SOME_VALUES: Setting f, s, p, and Pi options for %s\n",Pops[i]->Name);
			R->PopScalePriors[idx] = Pops[i]->ScalePrior;
			R->PopFactors[idx] = Pops[i]->Factor;
			R->PopPriors[idx] = Pops[i]->Prior;
			R->SimPi[idx].val = Pops[i]->SimPi; 
			
		}
	}
	if(Throw) {
		fprintf(stderr,"Error!  %d of the names associated with modified population priors were not names of populations found in the baseline. Exiting...\n",Throw);
	}
}



/*
	If you use the "-n name --Pi pi --end-pop -b Baseline.txt" type of specification 
	for the mixing proportions of populations from the Baseline in the mixtures simulated under
	the -x option, you have to get that Pi information into the TruePi vector.  This function
	returns that vector so it may be assigned to TruePi.
	
	It also handles the setup of everything for the TrueCounts if the --fixed-pi option was used.
	
	 Num is the number of pops specified with --end-pop.  
	 Pops is the array of pop_structs that came out of getting the options.
	  
*/
double *pi_vecFromPopStructPis(struct IndCollection *R, int Num, pop_struct **Pops, int **TrueCnts) 
{
	int i,k,idx;
	int Throw=0;
	double *ret = (double *)ECA_CALLOC(R->NumPops,sizeof(double));
	double normo = 0.0;
	int HasAbsCnts=0;
	
	/* first we are going to see if we have any absolute counts */
	for(i=0;i<Num;i++)  {
		if(Pops[i]->SimCnt>0) {
			HasAbsCnts=1;
		}
	}
	if(HasAbsCnts==1) {
		(*TrueCnts) = (int *)ECA_CALLOC(R->NumPops,sizeof(int));
	}
	
	for(i=0;i<Num;i++)  {
		/* do a really slow thing to find the name---would be better to have a hash table, but we only have to
		   do this once */
		idx = -1;
		for(k=0;k<R->NumPops;k++)  {
			if(strcmp(Pops[i]->Name,R->PopNames[k])==0) {
				idx = k;
				break;
			}
		}
		if(idx<0) {
			fprintf(stderr,"Error! Searched for name %s in baseline and didn't find it.  This will cause a failure...\n", Pops[i]->Name);
			Throw++;
		}
		else {
			if(HasAbsCnts) {
				(*TrueCnts)[idx] = Pops[i]->SimCnt;
				ret[idx] = (double)(Pops[i]->SimCnt);
				if(Pops[i]->SimCnt>0) {
					printf("SETTING_TRUE_CNTS_VALUES_FROM_FixedPi_OPTION: for %s.  Value is %d\n",Pops[i]->Name,Pops[i]->SimCnt);
				}
			}
			else {
				ret[idx] = Pops[i]->SimPi;
				if(Pops[i]->SimPi>0) {
					printf("SETTING_TRUE_PI_VALUES_FROM_Pi_OPTION: for %s.  Unscaled value is %f\n",Pops[i]->Name,Pops[i]->SimPi);
				}
			}
		}
	}
	if(Throw) {
		fprintf(stderr,"Error!  %d of the names associated with modified population priors were not names of populations found in the baseline. Exiting...\n",Throw);
		exit(1);
	}
	/* here we normalize, and also see if there was any information in there */
	for(i=0;i<R->NumPops;i++) {
		normo += ret[i];
	}
	if(normo==0.0) { /* if no information was given about relative weights, just return NULL */
		free(ret);
		return(NULL);  /* just sent it back NULL */
	}
	else { /* otherwise normalize */
		for(i=0;i<R->NumPops;i++) {
			ret[i] /= normo;
		}
	}
	return(ret);
	
}



/* if Factors, Priors or ScalePriors are null in an struct IndCollection, this sets them
 to default values.  The default values give you Rannala-Mountain with unit-info prior. */
void InitPopPriorsForIndCollection(struct IndCollection *R) 
{
	int i;
	
	
	if(R->PopScalePriors==NULL)  {
		R->PopScalePriors = (int *)ECA_CALLOC(R->NumPops,sizeof(int));
		for(i=0;i<R->NumPops;i++) {
			R->PopScalePriors[i] = 1;
		}
	}
	if(R->PopFactors==NULL)  {
		R->PopFactors = (double *)ECA_CALLOC(R->NumPops,sizeof(double)); 
		for(i=0;i<R->NumPops;i++) {
			R->PopFactors[i] = 1.0;
		}
	}
	if(R->PopPriors==NULL)  {
		R->PopPriors = (double *)ECA_CALLOC(R->NumPops,sizeof(double)); 
		for(i=0;i<R->NumPops;i++) {
			R->PopPriors[i] = 1.0;
		}
	}
	if(R->SimPi==NULL)  {
		R->SimPi = (popindval *)ECA_CALLOC(R->NumPops,sizeof(popindval)); 
		for(i=0;i<R->NumPops;i++) {
			R->SimPi[i].pop = i;
			R->SimPi[i].val = 0.0;
		}
	}
}

/* this takes an struct IndCollection and it uses the information in it to fill up 
some pop_structs, which it also allocates memory to.  

If you are filling up a Holdout Set, you need to pass the Baseline IndCollection
in as BaseReference so that it can figure out which pops in the holdout set correspond
to which pops in the Baseline.  This will be set in the IdxInBaseline field in the 
popsummary.  

Otherwise, just pass it as NULL and IdxInBaseline will remain set at -1.

 */
void IndColl2PopStructs(struct IndCollection *R, struct IndCollection *BaseReference) 
{
	int i,j,k, InPop, NumInPop=0;
	
	
	
	/* Now we will allocate memory and/or initialize values as necessary */
	if(R->PopSumms==NULL)  {
		R->PopSumms = (pop_struct **)ECA_CALLOC(R->NumPops,sizeof(pop_struct *));
	}
	for(i=0;i<R->NumPops;i++)  {
		if(R->PopSumms[i]==NULL) {
			R->PopSumms[i] = AllocPopStruct(R->NumLoc,R->AlleInfo->NumAlle,R->NumPops);
		}
		else { /* otherwise, just initialize to zero what we must to collect a sum */
			for(j=0;j<R->NumLoc;j++)  {
				R->PopSumms[i]->xsum[j] = 0.0;
				for(k=0;k<R->AlleInfo->NumAlle[j];k++)  {
					R->PopSumms[i]->x[j][k] = 0.0;
				}
			}
		}
	}
	if(R->PopStarts==NULL) {
		R->PopStarts = (int *)ECA_CALLOC(R->NumPops,sizeof(int));
	}
	if(R->PopEnds==NULL) {
		R->PopEnds = (int *)ECA_CALLOC(R->NumPops,sizeof(int));
	}
	if(R->NumInPop==NULL) {
		R->NumInPop = (int *)ECA_CALLOC(R->NumPops,sizeof(int));
	}
	
		
	/* and now we cycle over the individuals and summarize into the x and xsum fields */
	InPop = -1;
	for(i=0;i<R->NumInds;i++)  { int PN;
		PN = R->Inds[i]->PopNum;
		if(PN != InPop)  {
			R->PopStarts[PN] = i;
			if(InPop>=0) {
				R->PopEnds[InPop] = i-1;
				R->NumInPop[InPop] = NumInPop;
			}
			InPop = R->Inds[i]->PopNum;
			NumInPop = 0;
		}
		for(j=0;j<R->NumLoc;j++)  { int HasAll;
			/* if data are not missing at any allele */
			HasAll=1;
			for(k=0;k<R->AlleInfo->Ploidy[j];k++)  {
				if(R->Inds[i]->y[j][k] < 0) {
					HasAll=0;
					break;
				}
			}
			if(HasAll)  { 
				R->PopSumms[PN]->xsum[j] += (double)R->AlleInfo->Ploidy[j];
				for(k=0;k<R->AlleInfo->Ploidy[j];k++) {
					R->PopSumms[PN]->x[j][ R->Inds[i]->y[j][k] ] += 1.0;
				}
			}
		}
		NumInPop++;
		
		/* we have to assign PopEnds and the NumInPop for the very last individual */
		if(i==R->NumInds-1) {
			R->NumInPop[PN] = NumInPop;
			R->PopEnds[PN] = i;
		}
	}
	
	/* now that that is done, we can fill the xDP field---the Dirichlet Pars for the population,
		and, at the same time, compute and store the xDPsum's.  And we'll fill out the Ploidy
		while we are at it too. */
	for(i=0;i<R->NumPops;i++)  {
	
			
		/* put a pointer in the PopSumms back to this original struct IndCollection,
		and record which population in this struct IndCollection it is */
		R->PopSumms[i]->OriginalIndColl = R;
		R->PopSumms[i]->IdxOfPopInIndColl = i;
		
		/* copy the factor, prior, scaleprior info */
		R->PopSumms[i]->Factor = R->PopFactors[i];
		R->PopSumms[i]->Prior = R->PopPriors[i];
		R->PopSumms[i]->ScalePrior = R->PopScalePriors[i];
		
		/* copy over the NumLoc info */
		R->PopSumms[i]->NumLoc = R->NumLoc;
		R->PopSumms[i]->NumAlle = (int *)ECA_CALLOC(R->NumLoc,sizeof(int));
		
		for(j=0;j<R->NumLoc;j++)  { double add;
		
			R->PopSumms[i]->NumAlle[j] = R->AlleInfo->NumAlle[j];
			
			if(R->PopSumms[i]->ScalePrior) {
				add = R->PopSumms[i]->Prior / R->AlleInfo->NumAlle[j];
			}
			else {
				add = R->PopSumms[i]->Prior;
			}
			
			R->PopSumms[i]->xDPsum[j] = 0.0;  /* initialize to accumulate a sum */
			for(k=0;k<R->AlleInfo->NumAlle[j];k++) {
				R->PopSumms[i]->xDP[j][k] = R->PopSumms[i]->Factor * R->PopSumms[i]->x[j][k] + add;
				R->PopSumms[i]->xDPsum[j] += R->PopSumms[i]->xDP[j][k];
				
				R->PopSumms[i]->xBP[j][k] = R->PopSumms[i]->x[j][k] + add;
				
				R->PopSumms[i]->p[j][k] = R->PopSumms[i]->x[j][k] / R->PopSumms[i]->xsum[j];
				
			}
			
			/* copy the Ploidy info too */
			R->PopSumms[i]->Ploidy[j] = R->AlleInfo->Ploidy[j];
			
		}
	}
	
	
	/* in the end, we ought to copy the names over.  And, while we do that, we will fill out the IdxInBaseline field */
	for(i=0;i<R->NumPops;i++)  {
		sprintf(R->PopSumms[i]->Name,"%s",R->PopNames[i]);
		
		/* if this is for the holdout sample */
		if(BaseReference != NULL)  { int k; int foundit=0;
			for(k=0;k<BaseReference->NumPops;k++)  {
				if(strcmp(R->PopSumms[i]->Name,BaseReference->PopNames[k])==0) {
					R->PopSumms[i]->IdxInBaseline = k;
					foundit=1;
					printf("HOLDOUT_POPULATION %s is number %d and Corresponds to BASELINE POPULATION %d\n",R->PopSumms[i]->Name,i,R->PopSumms[i]->IdxInBaseline);
					break;
				}
			}
			if(!foundit) {
				printf("HOLDOUT_POPULATION %s is number %d and Corresponds to NO BASELINE POPULATIONS\n",R->PopSumms[i]->Name,i);
			}
		}
	}
}


/* Print a Summary of the information in an struct IndCollection */
void PrintIndCollection(struct IndCollection *R) 
{
	int i,j,k;
	printf("IND_COLL_SUMMARY_basic: NumInds= %d  NumLoc= %d  NumPops= %d\n",R->NumInds, R->NumLoc, R->NumPops);
	printf("IND_COLL_SUMMARY_popnames:");
	for(i=0;i<R->NumPops;i++)  {
		printf(" %s",R->PopNames[i]);
	}
	printf("\n");
	
	/* then print out which individuals are in which populations, etc */
	for(i=0;i<R->NumPops;i++)  {
		printf("IND_COLL_SUMMARY_popsummaries: Population %d is %s with %d indivs with indices between %d and %d inclusive\n",
				i,R->PopNames[i],R->NumInPop[i],R->PopStarts[i],R->PopEnds[i]);
	}
	
	/* now, summarize the locus names and allele hash stuff */
	for(j=0;j<R->NumLoc;j++) {
		printf("IND_COLL_SUMMARY_loci: Locus %s has %d alleles named as follows: ",R->AlleInfo->LocusNames[j],R->AlleInfo->NumAlle[j]);
		for(k=0;k<MAX_ALLELE_LENGTH+1;k++)  {
			if(R->AlleInfo->LocHash[j][k]>0) {
				printf("%d => %d   ",k,R->AlleInfo->LocHash[j][k]-1); 
			}
		}
		printf("\n");
	}
	
	/* and now we print the summarized the population information if PopSumm is not NULL */
	if(R->PopSumms != NULL) {
		for(j=0;j<R->NumLoc;j++) {
			printf("IND_COLL_SUMMARY_popfreqsum: LOCUS=%s",R->AlleInfo->LocusNames[j]);
			for(k=0;k<MAX_ALLELE_LENGTH+1;k++) {
				if(R->AlleInfo->LocHash[j][k]>0) {
					printf("\t%d",k);
				}
			}
			printf("\n");
			for(i=0;i<R->NumPops;i++)  {
				printf("IND_COLL_SUMMARY_popfreqsum: Pop=%s",R->PopSumms[i]->Name);
				for(k=0;k<MAX_ALLELE_LENGTH+1;k++) {
					if(R->AlleInfo->LocHash[j][k]>0) {
						printf("\t%.0f (%.2f)",R->PopSumms[i]->x[j][ R->AlleInfo->LocHash[j][k] - 1 ],
											R->PopSumms[i]->xDP[j][ R->AlleInfo->LocHash[j][k] - 1 ]);
					}
				}
				printf("\n");
			}
			for(i=0;i<R->NumPops;i++)  {
				printf("IND_COLL_SUMMARY_missing_data_report: Pop= %s   PopN= %d  Locus= %s   NonMissingGeneCopies= %.1f   FractionMissing= %.4f    \n",
						R->PopSumms[i]->Name, R->NumInPop[i], R->AlleInfo->LocusNames[j],R->PopSumms[i]->xsum[j],
					    (R->NumInPop[i]*2.0 - R->PopSumms[i]->xsum[j])/(R->NumInPop[i]*2.0));
			}
		}
	}
}


/* given an array of Ploidy allelic types, this makes an array of the number of occurrences 
of each type.  K is the number */
int **Compactify(int Ploidy, int *y, int *Appearing)
{
	int i,j,A=0, InIt=0;
	int **ret;
	int cur;
	int Missing = 0;
		
	/* first, scan for missing values.  Any missing values means you flush this down the toilet */
	for(i=0;i<Ploidy;i++)  {
		if(y[i]<0) {
			Missing++;
		}
	}
	if(Missing>0) {
		*Appearing=0;
		return(NULL);
	}
	
	ret = (int **)ECA_CALLOC(Ploidy,sizeof(int *));
	for(i=0;i<Ploidy;i++)  {
		cur = y[i];  /* get the allelic type of the i-th gene copy */
		
		/* see if the individual has another copy of this */
		InIt=0;
		for(j=0;j<A;j++) {
			if(cur==ret[j][0]) {
				InIt=1;
				break;
			}
		}
		if(InIt) { /* if he does, count another copy of it */
			ret[j][1]++;
		}
		else {  /* if not, add it to the list and increment A */
			ret[j]=(int *)ECA_CALLOC(2,sizeof(int));
			ret[j][0] = cur;
			ret[j][1]++; 
			A++;
		}
	}
	*Appearing = A;
	return(ret);
}  

/* to read individuals from a file.  This returns a pointer to 
   all the indivs */
void ReadIndsFromFile(char *FileName, struct IndCollection *TheInds) 
{
	int i,j;
	FILE *in;
	ind_struct **I;
	int NumInds, NumLoc;
	char **LN;  /* Local LocNames */
	int *Ploidies; /* local ploidies */
	char **PN=(char **)ECA_CALLOC(MAX_POPS,sizeof(char*)); /* local popnames */
	char temp[MAX_POP_NAME_LENGTH];	
	int NumPops=0, NumInPop;
	register int Verbose = gVerbosity;
	
	if( (in=fopen(FileName,"r"))==NULL ) {
		fprintf(stderr,"\n\nCouldn't open file \"%s\" to get individuals.\nExiting...\n\n",FileName);
		exit(1);
	}
	
	/* read the number of inds and the number of loci */
	fscanf(in," %d %d", &NumInds,&NumLoc);
	printf("READ_INDS: NumInds= %d NumLoc= %d\n",NumInds,NumLoc);
	
	/* then read the locus names, allocating memory first */
	LN = (char **)ECA_CALLOC(NumLoc, sizeof(char *));
	Ploidies = (int *)ECA_CALLOC(NumLoc, sizeof(int));
	for(i=0;i<NumLoc;i++) { char temp[MAX_LOC_NAME_LEN]; int ploid;
		LN[i] = (char *)ECA_CALLOC(MAX_LOC_NAME_LEN,sizeof(char));
		fscanf(in," %s",temp);
		if(strcmp(temp,"PLOIDY")==0) {
			fscanf(in," %d",&ploid);
			Ploidies[i] = ploid;
			fscanf(in," %s",temp); /* then scan for locus name */
			if(strcmp(temp,"PLOIDY")==0) {
				fprintf(stderr,"Read PLOIDY specifier twice without reading a locus name.  Locus i=%d.  Exiting...\n",i);
				exit(1);
			}
		}
		else {
			Ploidies[i] = 2;
		}
		sprintf(LN[i],"%s",temp);
		printf("READ_INDS_LocusName: %s  PLOIDY = %d\n",LN[i],Ploidies[i]);
	}
	
	/* now, check for POP specifier and a name */
	fscanf(in," %s",temp);
	if(strcmp(temp,"POP")==0 || strcmp(temp,"Pop")==0 || strcmp(temp,"pop")==0) {
		PN[NumPops]=(char *)ECA_CALLOC(MAX_POP_NAME_LENGTH,sizeof(char));
		fscanf(in," %s",PN[NumPops]);
		printf("READ_INDS_Reached_Population: %s\n",PN[NumPops]);
		NumPops++;
		NumInPop=0;
	}
	else {
		fprintf(stderr,"Error! Expecting POP specifier and pop name in file %s. Exiting...\n",FileName);
		exit(1);
	}
	
	/* now allocate memory to the indiv array, and the lochash if necessary, then read all the indivs */
	I = (ind_struct **)ECA_CALLOC(NumInds,sizeof(ind_struct *));
	if(TheInds->AlleInfo->LocHash==NULL)  {
		TheInds->AlleInfo->NumAlle = (int *)ECA_CALLOC(NumLoc,sizeof(int));
		TheInds->AlleInfo->LocHash = (int **)ECA_CALLOC(NumLoc,sizeof(int*));
		TheInds->AlleInfo->AlleNames = (int **)ECA_CALLOC(NumLoc,sizeof(int*));
		for(i=0;i<NumLoc;i++)  {
			TheInds->AlleInfo->LocHash[i] = (int *)ECA_CALLOC(MAX_ALLELE_LENGTH+1,sizeof(int));
			TheInds->AlleInfo->AlleNames[i] = (int *)ECA_CALLOC(MAX_ALLELES,sizeof(int));
		}
	}
	
	printf("About to cycle over indivs\n");

	for(i=0;i<NumInds;i++)  {
		/* first, check to see if we are starting into a new population */
		fscanf(in," %s",temp);
		if(strcmp(temp,"POP")==0 || strcmp(temp,"Pop")==0 || strcmp(temp,"pop")==0) {
			PN[NumPops]=(char *)ECA_CALLOC(MAX_POP_NAME_LENGTH,sizeof(char));
			fscanf(in," %s",PN[NumPops]);
			printf("READ_INDS_Reached_Population: %s\n",PN[NumPops]);
			NumPops++;
			NumInPop=0;
			i--; /* take one back, because this didn't read an indivdidual */
		}
		else { /* otherwise we get the individual */
			/* first, allocate memory to it */
			I[i] = (ind_struct *)ECA_MALLOC(sizeof(ind_struct));
			I[i]->PopNum = NumPops-1;
			I[i]->NumInPop = NumInPop++;
			I[i]->NumLoc = NumLoc;
			sprintf(I[i]->Name,"%s",temp);
			I[i]->y = (int **)ECA_CALLOC(NumLoc,sizeof(int *));
			I[i]->CompactRep = (int ***)ECA_CALLOC(NumLoc,sizeof(int *));
			I[i]->AllelesAppearing = (int *)ECA_CALLOC(NumLoc,sizeof(int));
			I[i]->posts = NULL;
			I[i]->logls = NULL;
			if(Verbose) {
				printf("READ_INDS_Indiv: %s has index %d overall and index %d in Pop %d : ",I[i]->Name,i,I[i]->NumInPop,I[i]->PopNum);
			}
			
			/* then get the alleles that they carry */
			for(j=0;j<NumLoc;j++)  {  int a; int k;
				I[i]->y[j] = (int *)ECA_CALLOC(Ploidies[j], sizeof(int));
				for(k=0;k<Ploidies[j];k++) {
					fscanf(in," %d",&a);
					/*printf( " [a=%d] ",a);*/
					if(a<=0) {  /* so MISSING DATA IS REPRESENTED BY A ZERO OR A NEGATIVE NUMBER!! */
						I[i]->y[j][k] = -1;
					}
					else {
						if(TheInds->AlleInfo->LocHash[j][a]==0) {  /* if we've not already seen this allele */
							TheInds->AlleInfo->LocHash[j][a] = ++(TheInds->AlleInfo->NumAlle[j]);
							TheInds->AlleInfo->AlleNames[j][TheInds->AlleInfo->NumAlle[j] - 1] = a;
						}
						I[i]->y[j][k] = (TheInds->AlleInfo->LocHash[j][a]) - 1;
					}
				}	
				
				/* now, store that in our "Compact" (which isn't so compact in terms of storage!) representation */
				I[i]->CompactRep[j] = Compactify(Ploidies[j], I[i]->y[j], &(I[i]->AllelesAppearing[j]));
				if(Verbose) { printf( "  %d %d",I[i]->y[j][0],I[i]->y[j][1]);}
				for(k=1;k<Ploidies[j];k++)  {
					if(Verbose) { printf( "/%d",I[i]->y[j][k]);}
				}
				if(Verbose) { printf(" ("); }
				for(k=0;k<Ploidies[j];k++)  {
					if(I[i]->y[j][0]>=0) {
						if(Verbose) { printf("%d",TheInds->AlleInfo->AlleNames[j][ I[i]->y[j][k] ]);}
					}
					else {
						if(Verbose) { printf("0");}
					}
					if(k<Ploidies[j]-1) {if(Verbose) { printf("/"); }}
					else if(Verbose) { printf(")");}
				}
				/* and in square brackets we print out the "Compact" representation */
				if(Verbose) { printf(" [ "); }
				for(k=0;k<I[i]->AllelesAppearing[j];k++)  {
					if(Verbose) { printf("%d->%d ",I[i]->CompactRep[j][k][0],I[i]->CompactRep[j][k][1]); }
				}
				if(Verbose) { printf("]"); }
				
			}
		}
		if(Verbose) { printf("\n");}
	}
	
	
	/* down here at the end, we assign all the local variables to the appropriate spots in TheInds */
	TheInds->NumInds = NumInds;
	TheInds->NumLoc = NumLoc;
	TheInds->NumPops = NumPops;
	TheInds->Inds = I;
	TheInds->AlleInfo->LocusNames = LN;
	TheInds->AlleInfo->Ploidy = Ploidies;
	TheInds->PopNames = PN;
	
}



/*!
Here we have a tidy little block of code to do 
the simulations of mixtures.  We cycle over the
number of mixtures to simulate, and for each one
we simulate the indivdiuals, then we do an EM-algorithm
to get the mle's of the mixing proportions. Once we 
have those, we print out the plug-in estimates of the
posterior probabilities of membership for each of the 
individuals.

if doing methods RESAMP_LOO_MULTILOCUS, you must have computed the posterior vectors in the 
struct IndCollection outside of this function. 
\todo Free memory...
*/
void MixtureSim(int INDS, int NumPops, pop_struct **Pops, int NumSimPops, pop_struct **SimPops,
                Method Meth, double *TruePi, int MIXSIZE,int MIXES,  struct IndCollection *Baselines, int *FixedNums,
                const char *MultiFixName, const char *MultiFixPrefix, const char *MultiFixNameHeader)
{
	int inds,i,j,k,mixes,FixedCntSum=0;;
	double **MixContainer;
	int Counts[MAX_POPS];
	
	/* count up the fixed nums if we have them */
	if(FixedNums != NULL)  {
		for(i=0;i<NumSimPops;i++) {
			Counts[i] = FixedNums[i];
			FixedCntSum += Counts[i];
		}
		MIXSIZE = FixedCntSum;
	}
	
	
	/* allocate to MixContainer */
	MixContainer = (double **)ECA_CALLOC(MIXSIZE, sizeof(double *));
	for(inds=0;inds<MIXSIZE;inds++) {
		MixContainer[inds] = (double *)ECA_CALLOC(NumPops,sizeof(double));
	}
	
	for(mixes=1;mixes<=MIXES;mixes++)  { double StartPoint[MAX_POPS]; double Pi[MAX_POPS]; int *TheirPops;
		
		if(FixedNums == NULL)  {
			for(i=0;i<NumSimPops;i++)  {  /* get set to count the actual numbers from each population.  If doing FixedNums then just input those here */			
				Counts[i] = 0;
			}
		}
		
		
		
		
		TheirPops = (int *)ECA_CALLOC(MIXSIZE,sizeof(int));  /* this is for recording the population they came from */
		
		
		if(FixedNums==NULL) {
			D_MultinomialRV(MIXSIZE, TruePi, NumSimPops, Counts);  /* draw the numbers from each pop */
		}

		
		for(j=0,i=0;i<NumSimPops;i++)  {
			for(k=0;k<Counts[i];k++)  {
				TheirPops[j++] = i;
			}
		}
				
		for(inds=0;inds<MIXSIZE;inds++) { int the_pop; /* simulate the individuals and store their likelihoods in MixContainer */
			 
			the_pop = TheirPops[inds];
			
			
			
			if(gPrintMixGenotypes) {
				printf("MIX_GENOS: Mixture= %d Pop= %d  MixInd= %d : Geno_Mix%d_Ind%d_pop%d ",mixes,the_pop,inds,mixes,inds,the_pop);
			}
			
			/* now, right here, if the we are doing the LOO approach inspired by conversations with Robin Waples, we
			don't simulate individuals, but rather we just put their pre-computed Leave-One-Out scaled likelihoods
			into the MixContainer */
			if(Meth==RESAMP_LOO_MULTILOCUS)  { int the_ind; int j;
				/* choose the index of the randomly chosen individual from the_pop */
				the_ind = UniformRV(Baselines->PopStarts[the_pop],Baselines->PopEnds[the_pop]);
				for(j=0;j<NumPops;j++)  {
					MixContainer[inds][j] = Baselines->Inds[the_ind]->posts[j].val;
				}
			}
			else {
				SimInd(the_pop,NumPops,Pops,SimPops,MixContainer[inds],Meth);
			}
		}
		
		/* print the TruePi's */
    if(TruePi != NULL) {
      printf("%sMIXFISH_TRUE_PIS: ", MultiFixPrefix);
      for(i=0;i<NumPops;i++)  {
        printf("%f ",TruePi[i]);
      }
      printf("\n");
    }
		
		/* print the actual numbers in the mixture */
		printf("%sMIXFISH_ACTUAL_NUMBERS: ", MultiFixPrefix);
		for(i=0;i<NumPops;i++)  {
			printf("%d ",Counts[i]);
		}
		printf("\n");
		
		/* now, use the results in MixContainer to find the MLE of the mixing proportions */
		/* first make the starting point uniform */
		for(i=0;i<NumPops;i++)  {
			StartPoint[i] = 1.0/(double)NumPops;
		}
		EM_Find_Pi_MLE(MixContainer, MIXSIZE, NumPops, StartPoint, Pi, .000001);
		
		printf("%sMIXFISH_PI_MLES: ", MultiFixPrefix);
		for(i=0;i<NumPops;i++)  {
			printf("%f ",Pi[i]);
		}
		printf("\n");
		
		
		/* print out headers for each mixture */
		printf("%sMIXFISH_NAMES_HEADER: %s  MixRep MixIndNum  PopNo  PopName ", MultiFixPrefix, MultiFixNameHeader);
		for(i=0;i<NumPops;i++)  {
			printf("%s ",Pops[i]->Name);
		}
		printf("\n");
		printf("%sMIX_FISH_NUMBERS_HEADER: %s  MixRep MixIndNum PopNo  PopName ", MultiFixPrefix, MultiFixNameHeader);
		for(i=0;i<NumPops;i++)  {
			printf("Post.%d ",i+1);
		}
		for(i=0;i<NumPops;i++)  {
			printf("ScLik.%d ",i+1);
		}
		printf("\n");
		
		/* then print out the fish in the mixture, and record the maximum while we are at it. */
		for(inds=0;inds<MIXSIZE;inds++) { double normo; double max; int maxi; double temp; double randoloss; double ScLikNorm;
			printf("%sMIXED_FISH_INDIVS:  %s  %d  %d  %d  %s  ",MultiFixPrefix, MultiFixName, mixes,inds+1,TheirPops[inds]+1, Pops[TheirPops[inds]]->Name);
			normo = 0.0;
			ScLikNorm = 0.0;
			max = -99999.9;
			for(i=0;i<NumPops;i++)  {
				temp = MixContainer[inds][i] * Pi[i];
				normo += temp;
				ScLikNorm += MixContainer[inds][i];
				if(temp>max) {
					max = temp;
					maxi = i;
				}
			}
			randoloss = 0.0;
			for(i=0;i<NumPops;i++)  {  /* here we print out the posterior prob */
				temp = MixContainer[inds][i] * Pi[i]/normo;
				printf("%e ",temp);
				if(Pops[TheirPops[inds]]->Loss != NULL)  {
					randoloss += temp * Pops[TheirPops[inds]]->Loss[i];
				}
			}
			for(i=0;i<NumPops;i++)  {  /* here we print out the scaled likelihood */
				printf("%e ",MixContainer[inds][i]/ScLikNorm);
			}
			printf("\n");
			
			/* and down here we compute and print the loss for that, if there is a Loss vector for the population from which it came from */
			if(Pops[TheirPops[inds]]->Loss != NULL)  {
				printf("MIXED_FISH_LOSS: %d  %d  %d  %s : AsgnToNum= %d   AsgnToName= %s  Loss= %.4f\n",
					mixes,inds+1,TheirPops[inds]+1, Pops[TheirPops[inds]]->Name, maxi+1, Pops[maxi]->Name, Pops[TheirPops[inds]]->Loss[maxi]);
				printf("MIXED_FISH_RANDO_LOSS: %d  %d  %d  %s :  %f\n",mixes,inds+1,TheirPops[inds]+1, Pops[TheirPops[inds]]->Name, randoloss);
			}
		}

	}

  /* Free the MixContainer */
  for(inds=0;inds<MIXSIZE;inds++) {
		free(MixContainer[inds]);
	}
	free(MixContainer);
}



/* this function resamples alleles from the baselines and then puts the results
in xBP along with the required priors, etc */
void ResampleBaselines(int NumPops, pop_struct **Pop)
{
	int i,j,k;
	int A,N;
	int Counts[MAX_ALLELES];
	double dA;
	
	/* cycle over loci */
	for(j=0;j<Pop[0]->NumLoc;j++)  {
		A = Pop[0]->NumAlle[j];  /* number of alleles at the locus */
		dA = (double)A;
		/* now, cycle over the pops and fill their xBP fields */
		for(i=0;i<NumPops;i++)  {
			N = (int)Pop[i]->xsum[j];  /* this is the number of gene copies sampled from pop i at loc j */
			
			/* here we resample with replacement from the alleles in the baseline */
			D_MultinomialRV(N,Pop[i]->p[j],A,Counts);

			/* now, we combine that with the priors to make xBP */
			for(k=0;k<A;k++)  {
				Pop[i]->xBP[j][k] = Counts[k];
				Pop[i]->xBP[j][k] += Pop[i]->ScalePrior  ?  Pop[i]->Prior/dA :  Pop[i]->Prior;
			}
		}
	}
}


/*!
Here we have a block of code to do the
"Straight-Up Simmed" ones.  To get a feel 
for the distribution of likelihood values
and also for posterior distributions given
the true mixing proportions.  We do one individual
from each population INDS times.  

\param INDS The number of indivdiduals from each population to simulate
\param NumPops The number of populations
\param Pops The array holding all the population information
*/

void StraightUpSim(int INDS, int NumPops, pop_struct **Pops, int NumSimPops, pop_struct **SimPops, Method Meth, double *TruePi)
{
	int inds,i,j;
	double *Boing = (double *)ECA_CALLOC(NumPops, sizeof(double));  /* this is just to hold the output variable from SimInd */
	double max; 
	int maxi;
	double randoloss;
	
	
	for(inds=1;inds<=INDS;inds++) {
		if(inds==1) {
			printf("SIMPLE_IND_NAMES_HEADER: Rep  PopNo  PopName ");
			for(i=0;i<NumPops;i++)  {
				printf("%s ",Pops[i]->Name);
			}
			printf("\n");
			printf("SIMPLE_IND_NUMBERS_HEADER: Rep  PopNo  PopName ");
			for(i=0;i<NumPops;i++)  {
				printf("%d ",i+1);
			}
			printf("\n");
		}
		
		
		for(i=0;i<NumSimPops;i++)  { 
			SimInd(i,NumPops,Pops,SimPops,Boing,Meth);
			
			
			printf("SIMPLE_IND_SIM:  %d   %d  %s ",inds,i+1,SimPops[i]->Name);
			max = -99999.9;
			randoloss = 0.0;  /* initialize to accumulate a sum */
			for(j=0;j<NumPops;j++)  {
				printf("%.4e  ",Boing[j]);
				if(Boing[j]>max) {  /* keep track of the maximum for computing the loss */
					max=Boing[j];
					maxi = j;
				}
				if(Pops[i]->Loss != NULL) {
					randoloss += Pops[i]->Loss[j] * Boing[j];
				}
			}
			printf("\n");
			
			if(Pops[i]->Loss != NULL) {
				printf("SIMPLE_LOSS_IND:  %d   %d  %s : AsgnToNum= %d AsgnToName= %s Loss= %.4f\n",inds+1,i+1,Pops[i]->Name, maxi+1, Pops[maxi]->Name, Pops[i]->Loss[maxi]);
				printf("SIMPLE_RANDO_LOSS:  %d   %d  %s : %f\n",inds+1,i+1,Pops[i]->Name,randoloss);
			}
			
			if(TruePi) { double normo=0.0; double temp;
				for(j=0;j<NumPops;j++)  {
					normo += Boing[j] * TruePi[j];
				}
				printf("TRUE_WEIGHTED_IND %d  %d  %s ",inds+1,i+1,Pops[i]->Name);
				
				max = -9999.9;
				randoloss = 0.0;
				for(j=0;j<NumPops;j++)  {
					temp = Boing[j]*TruePi[j]/normo;
					printf("%.4e  ",temp);
					if(temp > max) {
						max = temp;
						maxi = j;
					}
					if(Pops[i]->Loss != NULL) {
						randoloss += temp * Pops[i]->Loss[j];
					}
				}
				printf("\n");
				
				/* finally print the losses */
				if(Pops[i]->Loss != NULL) {
					printf("TRUE_WEIGHTED_LOSS_IND:  %d   %d  %s : AsgnToNum= %d AsgnToName= %s Loss= %.4f\n",inds+1,i+1,Pops[i]->Name,maxi+1, Pops[maxi]->Name,Pops[i]->Loss[maxi]);
					printf("TRUE_WEIGHTED_RANDO_LOSS:  %d   %d  %s : %f\n",inds+1,i+1,Pops[i]->Name,randoloss);
				}
			}
		}
	}
}





/* compute the haploid Fst between pops A and B, using the method for estimating
theta given in Weir's Genetic Data Analysis II (pp. 172-174). */
double PairwiseFst(pop_struct *A, pop_struct *B)
{
	int j,k;
	double MSP,MSG,nc, pave;
	double T1 = 0.0, T2 = 0.0;
	
	for(j=0;j<A->NumLoc;j++)  {
		for(k=0;k<A->NumAlle[j];k++)  {
			pave = (A->x[j][k] + B->x[j][k]) / (A->xsum[j] + B->xsum[j]);
			MSP =  A->xsum[j] * (A->p[j][k] - pave) * (A->p[j][k] - pave) +
					B->xsum[j] * (B->p[j][k] - pave) * (B->p[j][k] - pave);
			MSG = ( A->xsum[j] * A->p[j][k] * (1.0 - A->p[j][k]) +
					B->xsum[j] * B->p[j][k] * (1.0 - B->p[j][k]) ) / 
					 ( A->xsum[j] + B->xsum[j] - 2.0 );
			nc = (A->xsum[j] + B->xsum[j]) - ( A->xsum[j]*A->xsum[j] + B->xsum[j]*B->xsum[j]) / (A->xsum[j] + B->xsum[j]);
			
			T1 += (MSP-MSG);
			T2 += MSP + (nc-1.0)*MSG; 
		}
	}
	
	return(T1/T2);
}


/* given a (N x C) matrix M of likelihoods for N individuals possibly having membership in 
   C different sources, and an initial starting condition SP for the mixing proportions,
   this function runs a simple EM algorithm to find the MLE of the mixing proportions.  
   This MLE is output in the output variable Pi which must have memory allocated to it 
   already.  The EM algorithm terminates when the sum of absolute difference between 
   components is less then TOL, which should be 1e-05 or so.
   
   Note that for this function it doesn't matter if SP sums to one or not. 
*/
void EM_Find_Pi_MLE(double **M, int N, int C, double *SP, double *Pi, double TOL)
{
	int i,j, steps=0;
	double *old = (double *)ECA_CALLOC(C,sizeof(double));
	double *cur = (double *)ECA_CALLOC(C,sizeof(double));
	double diff=0.0, normo;
	
	/* copy SP to old */
	for(i=0;i<C;i++)  {
		old[i] = SP[i];
	}
	
	/* then start the loop */
	do {
		/* set cur to zeroes */
		for(i=0;i<C;i++)  {
			cur[i] = 0.0;
		}
		
		/* cycle over individuals and set cur to the average of their posterior probs */
		for(j=0;j<N;j++)  {
			normo=0.0;
			for(i=0;i<C;i++)  {
				normo += old[i] * M[j][i];
			}
			for(i=0;i<C;i++)  {
				cur[i] += old[i] * M[j][i] / (normo * (double)N);
			}
		}
		
		/* compute the difference between old and cur, and, at the same time, set old equal to cur */
		diff = 0.0;
		for(i=0;i<C;i++)  {
			diff += ECA_ABS(cur[i] - old[i]);
			old[i] = cur[i];
		}
		printf("EM_ALGORITHM_PROGRESS:  %d  :  %f  : ",++steps,diff);
		for(i=0;i<C;i++)  {
			printf("%f ",cur[i]);
		}
		printf("\n");
		
	} while (diff>TOL);
	
	/* copy cur over to Pi and we are done */
	for(i=0;i<C;i++)  {
		Pi[i] = cur[i];
	}
	
}



void PrintLossMatrix(pop_struct **P, int NumPops)
{
	int i,j;
	
	/* first print the header */
	printf("LOSSMATRIX : Name");
	for(i=0;i<NumPops;i++)  {
		printf("\t%s",P[i]->Name);
	}
	printf("\n");
	
	/* then print the rest */
	for(j=0;j<NumPops;j++)  {
		printf("LOSSMATRIX : %s",P[j]->Name);
		for(i=0;i<NumPops;i++)  {
			printf("\t%.4f",P[j]->Loss[i]);
		}
		printf("\n");
	}
	
}


/* in order to draw a single locus from an individuals in a certain range in an IndCollection.
This is used in SimInd for the Meth==RESAMP_LOO_SINGLELOCUS. 

A is the number of alleles at the locus.
Y is for returning the genotype of that locus.
Ploidy is the ploidy at the locus.
LocIdx is the index of the locus that we wish to sample a realization for.
Orig is a pointer to the IndCollection from which the pop was summarized.
Idx is the index of the pop within the IndCollection pointed to by OrigIndColl.

If OrigIndColl is NULL, this exits with an error.


  */
void DrawSingleLocus(int A, int *Y, int Ploidy, int LocIdx, struct IndCollection *Orig, int Idx )
{
	int i,rando,keep_at_it;
	int *Loc;
	
	/* initialize Y */
	for(i=0;i<A;i++)  {
		Y[i] = 0;
	}
	
	
	
/* for debugging 
printf("Entered DrawSingleLocus A=%d, LocIdx=%d,  Start=%d,  Stop=%d, ",A,LocIdx,Orig->PopStarts[Idx],Orig->PopEnds[Idx]); */

	
	/* now, as it currently stands, we are going to sample loci, until we find one that is not missing. 
	In the future, we can change this to allow for missing data of various sorts */
	do {
		keep_at_it=0;
		rando = UniformRV(Orig->PopStarts[Idx],Orig->PopEnds[Idx]);
		Loc = Orig->Inds[rando]->y[LocIdx];
		for(i=0;i<Ploidy;i++)  {
			if(Loc[i]<0) {
				keep_at_it++;
			}
		}
/* printf("rando=%d  types: ",rando); */
		if(!keep_at_it)  {  /* if there was no missing data at the locus */
			for(i=0;i<Ploidy;i++)  {
				Y[Loc[i]]++;
/* printf("%d ",Loc[i]); */
			}
		}
/* printf("\n"); */
	} while(keep_at_it);
}


/* simulate an individual from population J with characteristics given by SimPops[J],
 and compute the likelihood that it came from each of the N different populations whose characteristics are given
in the Pops array. 

In some instances Pops and SimPops will be the same (i.e., the population being simulated
from will be the same as the population being used to do the inference).  This will be apparent
because Pops==SimPops (i.e. the pointers will have the same value).  In such a case, it will
always be appropriate to use a leave-one-out procedure of some sort.  In some cases, however,
SimPops might refer to a holdout set---it may include indivdiuals that were not included when 
locus selection was done, for example.  In such a case, it would be typical not to do any sort of 
leave-one-out thing.  

However, occasionally you may want to include the holdout individuals in the baseline (after the
original baseline was used to select the loci, for example), in which case you DO want to
leave-them-out of the baseline.  So, for that case we have to be a little bit tricky and 
use the IdxOfBaseline field in the popstruct.

The resulting likelihoods are put into the output
variable Outs.  They are scaled so that the largest of them equals 1.000  

MIXSIZE is the number of individuals in the mixture and CountsInPop are the counts of
individuals from each population in the mixture.  

Note that there are N populations in Pops, but there may be a different number of pop ins SimPops.

*/
void SimInd(int J, int N, pop_struct **Pops, pop_struct **SimPops, double *Outs, Method Meth) 
{
	int i,j,k;
	int L=Pops[0]->NumLoc; /* number of loci */
	int A;  /* temp to hold number of alleles */
	int Ploidy; /* temp to hold the ploidy of a locus */
	int *Y = (int *)ECA_CALLOC(MAX_ALLELES,sizeof(int));  /* to gather the results for individuals */
	double Max = -9999999999.9, normo=0.0;
	double Phony_xBP[MAX_ALLELES];
	int PopsIsNotSimPops=0; /* 0 if SimPops and Pops are the same, 1 otherwise */
	int BaseIdxOfJ=J;  /* if Pops==SimPops, then this is correct.  Otherwise we must change it */
	
	
	if(Pops!=SimPops) {
		PopsIsNotSimPops = 1;
		BaseIdxOfJ = SimPops[J]->IdxInBaseline;
	}
	
	/* first, initialize Outs to 0.0 to accumulate a sum of logs */
	for(i=0;i<N;i++) {
		Outs[i] = 0.0;
	}
	
	/* cycle over loci */
	for(j=0;j<L;j++)  {
	
		A = Pops[0]->NumAlle[j];
		Ploidy = Pops[0]->Ploidy[j];
		
		/* draw the genotype for the simulated fish in the appropriate manner */
		if(Meth==RESAMP_LOO_GENECOPY) {
			DrawWithoutReplacement_d_Squash(A, SimPops[J]->x[j], Y, Ploidy);
		}
		else if(Meth==NO_LOO_GENECOPY) {
			D_MultinomialRV(Ploidy,SimPops[J]->x[j], A ,Y);
		}
		else if(Meth==RESAMP_LOO_SINGLELOCUS || Meth==NO_LOO_SINGLELOCUS)  {  /* in this case, we have to sample from amongst the loci in the 
													  individuals in the IndCollection from which the pop_struct was summarized */
			DrawSingleLocus(A, Y, Ploidy, j, SimPops[J]->OriginalIndColl, SimPops[J]->IdxOfPopInIndColl );
		}
		else {
			fprintf(stderr,"Error! Unknown Method in SimInd.  Exiting!\n");
			exit(1);
		}
		
		
		/* here we print out the simulated genotypes */
		if(gPrintMixGenotypes) { int a; int b;
			printf("    ");
			for(a=0;a<A;a++) {  /* cycle over all possible alleles */
				for(b=0;b<Y[a];b++) {
					if(Pops[J]->OriginalIndColl != NULL) {  /* if the alleles came from a genotype list then the alleles have real
					                                           names and we ought to use them */
						printf("%d ",Pops[J]->OriginalIndColl->AlleInfo->AlleNames[j][a]);
					}
					else {
						printf("%d ",a);
					}
				}
			}
			if(j==L-1) {
				printf("\n");
			}
		}
		
		/* now, compute the probability of that genotype coming up in each of the N populations in Pops.  Note that
		using this function, we are taking care of the factor of 2 in heterozygotes automatically */
		for(i=0;i<N;i++)  {
			/* here is the drill.  If i != BaseIdxOfJ, then we never do any sort of leave one out thing.  However,
			if i==BaseIdxOfJ, then we do leave one out for all he LOO methods and we don't if Meth is any of the NO_LOO
			methods.  Pretty simple, but not very elegant */
			if(i != BaseIdxOfJ) {
				Outs[i] += log(CMDPMF_small_N(Y, Pops[i]->xBP[j], A, Ploidy, -1));
			}
			else {
				if(Meth==NO_LOO_GENECOPY || Meth==NO_LOO_SINGLELOCUS || Meth==NO_LOO_MULTILOCUS) {
					Outs[i] += log(CMDPMF_small_N(Y, Pops[i]->xBP[j], A, Ploidy, -1));
				}
				else if(Meth==RESAMP_LOO_GENECOPY || Meth==RESAMP_LOO_SINGLELOCUS || Meth==RESAMP_LOO_MULTILOCUS) { double Xsum;
					/* in this case, we have to effect a leave-one-out procedure.  This is easily done by copying xBP to a new
					variable, called Phony_xBP, then subtracting from that the genotype values Y.  While we are at it, we can
					compute the sum of Phony_xBP and include that in the call to CMDPMF_small_N, so it doesn't take much more time. */
					for(Xsum=0.0,k=0;k<A;k++)  {
						Phony_xBP[k] = Pops[i]->xBP[j][k] - Y[k];
						Xsum += Phony_xBP[k];
					}
					Outs[i] += log(CMDPMF_small_N(Y, Phony_xBP, A, Ploidy, Xsum));
				}
			}
		}
		
		
		/* then zero out the Y vector */
		for(k=0;k<A;k++) {
			Y[k] = 0;
		}
	}
	
	/* then find the maximum of the Outs and normalize and exponentiate */
	for(i=0;i<N;i++)  {
		if(Outs[i] > Max) 
			Max = Outs[i];
	}
	for(i=0;i<N;i++)  {
		Outs[i] = exp(Outs[i]-Max);
		normo += Outs[i];
	}
	
	/* then scale them all so they sum to one */
	for(i=0;i<N;i++)  {
		Outs[i] /= normo;
	}
	
	free(Y);
}


/*
	This is almost exactly like SimInd, but it just computes the log10 likelihood that the simulated
	individual came from population K.  *and* in which the simulated individual has missing data pattern
	as given in the array Yobs (values <0 in Y indicate missing data ).  It returns this Log10 value.
 
	NumTyped is the number of non-missing loci.
 
 */
double SimIndForLoglDsnFromSinglePop(int J, pop_struct **Pops, pop_struct **SimPops, Method Meth, int K, int **Yobs, int *NumTyped) 
{
	int i,j,k,xx;
	int L=Pops[0]->NumLoc; /* number of loci */
	int A;  /* temp to hold number of alleles */
	int Ploidy; /* temp to hold the ploidy of a locus */
	int *Y = (int *)ECA_CALLOC(MAX_ALLELES,sizeof(int));  /* to gather the results for individuals */
	double Phony_xBP[MAX_ALLELES],AlleCountSum;
	int PopsIsNotSimPops=0; /* 0 if SimPops and Pops are the same, 1 otherwise */
	int BaseIdxOfJ=J;  /* if Pops==SimPops, then this is correct.  Otherwise we must change it */
	double LoglSum;
	int NotMissing=0;
	
	if(Pops!=SimPops) {
		PopsIsNotSimPops = 1;
		BaseIdxOfJ = SimPops[J]->IdxInBaseline;
	}
	
	/* first, initialize LoglSum to 0.0 to accumulate a sum of logs */
	LoglSum = 0.0;
	
	/* cycle over loci */
	for(j=0;j<L;j++)  {
		
		
		if(Yobs[j][0]>=0 && Yobs[j][1]>=0   ) {  /* only simulate and compute a likelihood for loci not missing data in the individuals Yobs's genotype */
			
			NotMissing++;
			A = Pops[0]->NumAlle[j];
			Ploidy = Pops[0]->Ploidy[j];
			
			/* now, we have to check to see what the sum of the SimPops[J]->x[j] vector is.  
			 If the sum is less than Ploidy it will crap out when we try RESAMP_LOO drawing from it.  If it is 
			 zero it will crap out if we try drawing a multinomial from it.  This sort of 
			 situation should only occur for populations in which nearly everyone is missing data at a particular locus in the reference
			 sample, so that locus ought to be skipped anyway.  Except of course, we are trying to use this to compare LogLs with the 
			 observed individual.  So, if the reference sample doesn't have any observations for this locus, then we should just
			 choose the locus from the observed individual.  */
			AlleCountSum=0.0;
			for(xx=0;xx<A;xx++) {
				AlleCountSum+=SimPops[J]->x[j][xx];
				Y[xx]=0;  /* go ahead and zero out Y, in case we need to fill in the observed value by just adding some 1's in there */
			}
			
			/* draw the genotype for the simulated fish in the appropriate manner */
			if(Meth==RESAMP_LOO_GENECOPY) {
				if(AlleCountSum>Ploidy) {
					DrawWithoutReplacement_d_Squash(A, SimPops[J]->x[j], Y, Ploidy);
				}
				else {
					Y[Yobs[j][0]]++;
					Y[Yobs[j][1]]++;
				}
			}
			else if(Meth==NO_LOO_GENECOPY) {
				if(AlleCountSum>0) {
					D_MultinomialRV(Ploidy,SimPops[J]->x[j], A ,Y);
				}
				else {
					Y[Yobs[j][0]]++;
					Y[Yobs[j][1]]++;
				}
			}
			else if(Meth==RESAMP_LOO_SINGLELOCUS || Meth==NO_LOO_SINGLELOCUS)  {  /* in this case, we have to sample from amongst the loci in the 
			 individuals in the IndCollection from which the pop_struct was summarized */
				if(AlleCountSum>0) {
					DrawSingleLocus(A, Y, Ploidy, j, SimPops[J]->OriginalIndColl, SimPops[J]->IdxOfPopInIndColl );
				}
				else {
					Y[Yobs[j][0]]++;
					Y[Yobs[j][1]]++;
				}
			}
			else {
				fprintf(stderr,"Error! Unknown Method in SimInd.  Exiting!\n");
				exit(1);
			}
			
			
			/* here we print out the simulated genotypes */
			if(gPrintMixGenotypes) { int a; int b;
				printf("    ");
				for(a=0;a<A;a++) {  /* cycle over all possible alleles */
					for(b=0;b<Y[a];b++) {
						if(Pops[J]->OriginalIndColl != NULL) {  /* if the alleles came from a genotype list then the alleles have real
						 names and we ought to use them */
							printf("%d ",Pops[J]->OriginalIndColl->AlleInfo->AlleNames[j][a]);
						}
						else {
							printf("%d ",a);
						}
					}
				}
				if(j==L-1) {
					printf("\n");
				}
			}
			
			/* now, compute the probability of that genotype coming up in Population K.  Note that
			 using this function, we are taking care of the factor of 2 in heterozygotes automatically */
			for(i=K;i==K;i++)  {    /* this is silly and lazy.  Originally this function cycled over the N pops, but I just changed this so it 
										picks out just Population K */
				/* here is the drill.  If i != BaseIdxOfJ, then we never do any sort of leave one out thing.  However,
				 if i==BaseIdxOfJ, then we do leave one out for all he LOO methods and we don't if Meth is any of the NO_LOO
				 methods.  Pretty simple, but not very elegant */
				if(i != BaseIdxOfJ) {
					LoglSum += log10(CMDPMF_small_N(Y, Pops[i]->xBP[j], A, Ploidy, -1));
				}
				else {
					if(Meth==NO_LOO_GENECOPY || Meth==NO_LOO_SINGLELOCUS || Meth==NO_LOO_MULTILOCUS) {
						LoglSum += log10(CMDPMF_small_N(Y, Pops[i]->xBP[j], A, Ploidy, -1));
					}
					else if(Meth==RESAMP_LOO_GENECOPY || Meth==RESAMP_LOO_SINGLELOCUS || Meth==RESAMP_LOO_MULTILOCUS) { double Xsum;
						/* in this case, we have to effect a leave-one-out procedure.  This is easily done by copying xBP to a new
						 variable, called Phony_xBP, then subtracting from that the genotype values Y.  While we are at it, we can
						 compute the sum of Phony_xBP and include that in the call to CMDPMF_small_N, so it doesn't take much more time.
						 
						 Note that there is an extra weirdness here.  Fishery samples assigned to the coho population might have non-missing data at
						 a locus which is missing in the entire "coho-within-chinook" baseline.  This causes things to bomb with LOO because I have just assigned
						 Yobs to Y instead of simulating something, but there is nothing to pull of out of the baseline to look like LOO---however, you don't need to
						 do LOO anyway!  So, I have a little thing to catch that two lines down */
						for(Xsum=0.0,k=0;k<A;k++)  {
							if(AlleCountSum>Ploidy) {
								Phony_xBP[k] = Pops[i]->xBP[j][k] - Y[k];
							}
							else {
								Phony_xBP[k] = Pops[i]->xBP[j][k];
							}
							Xsum += Phony_xBP[k];
						}
						LoglSum += log10(CMDPMF_small_N(Y, Phony_xBP, A, Ploidy, Xsum));
					}
				}
			}
			
			
			/* then zero out the Y vector */
			for(k=0;k<A;k++) {
				Y[k] = 0;
			}
		}  
	}
	free(Y);
	*NumTyped = NotMissing;
	return(LoglSum);
	
}




/* this should be run after an EM run.  It does simulations for each individual in the mixture
 to compare the LogL of that individual to what you expect from the population it is most likely from
 
 Note that if Baselines and TheMixture point to the same thing, this function will set the AssignedTo
 variable to the PopNum field because it means that we are going to be simulating this distribution for 
 a baseline individual.  
 
 */
void SimulateAllLoglsForComparison(struct IndCollection *Baselines, struct IndCollection *TheMixture, Method Meth, int NumReps, FILE *loglf, char *outfile_name)
{
	int i,j,AssignedTo;
	ind_struct *Ind;
	double ObsLogL,SimLogL;
	dval *LogVals;
	double z, sd;
	double NumSimGreater;
	int NumTyped;
	FILE *outf;
	
	if( (outf=fopen(outfile_name,"w"))==NULL) {
		fprintf(stderr,"Failed to open file %s to write to it.  Exiting.\n",outfile_name);
		exit(1);
	}
	
	LogVals = AllocDval(0,0,-1);
	InitDvalSummaryToZero(LogVals);
	
	
	fprintf(outf,"FishId\tAssignedTo\tNumLoci\tObsLogL\tzScore\tFractionSimmedGreater\n");
	/* cycle over the individuals */
	for(i=0;i<TheMixture->NumInds;i++)  {
		Ind =  TheMixture->Inds[i];
		if(Baselines==TheMixture) {
			AssignedTo = Ind->PopNum;
		}
		else {
			AssignedTo = Ind->posts[0].pop;
		}
		ObsLogL = Ind->logls[AssignedTo].val;
		NumSimGreater=0.0;
		
		if(loglf) {
			fprintf(loglf,"%s\t%f\n",Ind->Name,ObsLogL);
		}
		
		InitDvalSummaryToZero(LogVals);
		for(j=0;j<NumReps;j++) {
			SimLogL = SimIndForLoglDsnFromSinglePop(AssignedTo,  Baselines->PopSumms, Baselines->PopSumms, Meth, AssignedTo, Ind->y, &NumTyped);
			if(loglf) {
				fprintf(loglf,"%s\t%f\n",Ind->Name,SimLogL);
			}
			LogVals->v = SimLogL;
			IncrementDval(LogVals);
			if(ObsLogL > SimLogL)  {
				NumSimGreater += 1.0;
			}
		}
		sd = sqrt(LogVals->Var * LogVals->NumAved);
		z = (ObsLogL - LogVals->Ave)/sd;
		fprintf(outf,"%s\t%s\t%d\t%f\t%f\t%f\n",Ind->Name,Baselines->PopNames[AssignedTo],NumTyped,ObsLogL,z,NumSimGreater/NumReps);

		
	}
}








/* allocate memory to and initialize a pop_struct */
pop_struct *InitPopStruct(int I) 
{
	pop_struct *temp;
	
	temp = (pop_struct *)ECA_MALLOC(sizeof(pop_struct));
	
	temp->Idx = I;
	sprintf(temp->Name,"NoName");
	temp->NumLoc = 0;
	temp->x = NULL;
	temp->NumAlle = NULL;
	temp->Factor = 1.0;
	temp->ScalePrior = 1;
	temp->Prior = 1.0;
	
	temp->NumConstell = 0;
	temp->ConstelWts = NULL;
	temp->ConstelFst = NULL;
	
	temp->Loss = NULL;
	temp->ExpScaledLike = NULL;
	
	
	temp->SimPi = 0.0;
	temp->SimCnt = 0;
	
	
	return(temp);
}


/* Compute and print the expected loss for locus L.  If Pi==NULL, then do the calculation only under the assumption of equal priors.  */
void ExpectedLossLocL(int L, pop_struct **P, double *Pi, int NumPops)
{
	int i,j,k,a,b;
	double res=0.0, reseq=0.0;
	int A = P[0]->NumAlle[L]; /* number of alleles */
	int *Y;  /* to hold the genotypes */
	double EqPri = 1.0/NumPops;
	double eqnormo,normo;
	double Likes[MAX_POPS];
	double EqPriPosts[MAX_POPS];
	double PiPosts[MAX_POPS];
	
	/* allocate memory */
	Y = (int *)ECA_CALLOC(A,sizeof(int));
	
		
		
	/* cycle over the first gene copy */
	for(a=0;a<A;a++)  {
		
		Y[a]++;  /* make the genotype variable Y hold a 1 in the a-position */
		/* cycle over the second gene copy */
		for(b=a;b<A;b++)  {
			Y[b]++; /* add the b gene copy to the Y variable */
			
			/***** here the Y variable holds the current genotype ******/
			
			/* cycle over populations to compute the likelihoods */
			for(i=0;i<NumPops;i++)  {
				Likes[i] = exp(LogCMDPMF(Y, P[i]->xDP[j], A)); /* NOTE: This is where we use xDP!! */
			}
			
			/* now, make normalized versions */
			eqnormo = 0.0;
			normo = 0.0;
			for(k=0;k<NumPops;k++)  {
				EqPriPosts[k] = Likes[k];
				eqnormo += EqPriPosts[k];
				
				if(Pi != NULL) {
					PiPosts[k] = Likes[k] * Pi[k];
					normo += PiPosts[k];
				}
			}
			for(k=0;k<NumPops;k++)  {
				EqPriPosts[k] /= eqnormo;
				if(Pi != NULL) {
					PiPosts[k] /= normo;
				}
			}
			
			/* now, we have the posteriors, and we just need to sum over all the population origins and over all the losses
			possible with those, to get the expected losses.  Note that we are going to want to do both a randomized loss and a max-assignment loss here, eventually */
			for(k=0;k<NumPops;k++)  {  double smallsum; double smallsumeq; /* cycle over the population of origins */

				smallsum=0.0;
				smallsumeq=0.0;
				
				/* sum up the losses */
				for(j=0;j<NumPops;j++) {
					smallsumeq += EqPriPosts[j] * P[k]->Loss[j];
					if(Pi != NULL) {
						smallsum += PiPosts[j] * P[k]->Loss[j];
					}
				}
				
				/* then weight that by the probability that the fish came from pop k and had genotype Y */
				smallsumeq *= (EqPri * Likes[k]);
				if(Pi != NULL) {
					smallsum *= (Pi[k] * Likes[k]);
				}
				
				/* and go ahead and add that to the overall loss */
				res += smallsum;
				reseq += smallsumeq;
			
			}
			
			
			printf("\n");
			
			
			
			Y[b]--; /* remove the b gene copy from the Y variable */
		}
		
		Y[a]--; /* remove this gene copy from the a-position */
	}
		
		
		
	
	
	
	
	
	/* free memory */
	free(Y);
	
}


/* read data out of the named locus parameters file */
void ReadLocFile(pop_struct **P, int L, char *FileName)
{
	int j,k;
	FILE *in;
	
	/* now get all the locus pars */
	if( (in=fopen(FileName,"r"))==NULL ) {
		fprintf(stderr,"\n\nCouldn't open file \"%s\" for locus information.\nExiting...\n\n",FileName);
		exit(1);
	}
	
	printf("ALLECOUNTS : opened file \"%s\" to get genetic data for population %d\n",FileName,L+1);
	while(eat_comments(in,'&')) ;
	/* get the number of loci */
	fscanf(in,"%d",&(P[L]->NumLoc) );
	
	if(L>0) {
		if(P[L]->NumLoc != P[L-1]->NumLoc) {
			fprintf(stderr,"Number of loci in file %s is %d, but for previous population was %d.  Exiting...\n",FileName,P[L]->NumLoc,P[L-1]->NumLoc);
			exit(1);
		}
	}
	
	printf("ALLEFREQS : number of loci in file %s is %d\n",FileName,P[L]->NumLoc);
	
	while(eat_comments(in,'&')) ;
	/* then get all the rest of the locus parameters */
	P[L]->NumAlle = (int *)ECA_CALLOC(P[L]->NumLoc, sizeof(int));
	
	while(eat_comments(in,'&')) ;
	P[L]->x = (double **)ECA_CALLOC(P[L]->NumLoc, sizeof(double *));
	P[L]->xDP = (double **)ECA_CALLOC(P[L]->NumLoc, sizeof(double *));
	P[L]->xBP = (double **)ECA_CALLOC(P[L]->NumLoc, sizeof(double *));
	P[L]->p = (double **)ECA_CALLOC(P[L]->NumLoc, sizeof(double *));
	P[L]->simp = (double **)ECA_CALLOC(P[L]->NumLoc, sizeof(double *));
	P[L]->xsum = (double *)ECA_CALLOC(P[L]->NumLoc, sizeof(double));
	P[L]->xDPsum = (double *)ECA_CALLOC(P[L]->NumLoc, sizeof(double));
	
	for(j=0;j<P[L]->NumLoc;j++)  {
		while(eat_comments(in,'&')) ;
		fscanf(in,"%d",&(P[L]->NumAlle[j]));
		
		if(L>0) {
			if(P[L]->NumAlle[j] != P[L-1]->NumAlle[j]) {
				fprintf(stderr,"Number of alleles in file %s at locus %d is %d, but for previous population was %d.  Exiting...\n",FileName,j+1,P[L]->NumAlle[j],P[L-1]->NumAlle[j]);
				exit(1);
			}
		}
		
		P[L]->x[j] = (double *)ECA_CALLOC(P[L]->NumAlle[j], sizeof(double));
		P[L]->xDP[j] = (double *)ECA_CALLOC(P[L]->NumAlle[j], sizeof(double));
		P[L]->xBP[j] = (double *)ECA_CALLOC(P[L]->NumAlle[j], sizeof(double));
		P[L]->p[j] = (double *)ECA_CALLOC(P[L]->NumAlle[j], sizeof(double));
		P[L]->simp[j] = (double *)ECA_CALLOC(P[L]->NumAlle[j], sizeof(double));
		
		
		for(k=0;k<P[L]->NumAlle[j];k++)  { /* in this loop, we collect all of the allele counts, and we also make the Dirichlet priors up for 
											them based on the Factor, ScalePrior, and Prior variables. */
			while(eat_comments(in,'&')) ;
			fscanf(in,"%lf",&(P[L]->x[j][k]));
			
			/* add the x to the xsum */
			P[L]->xsum[j] += P[L]->x[j][k];
			
			/* set the Dirichlet Pars for simulating */
			P[L]->xDP[j][k] = P[L]->ScalePrior  ?  P[L]->Prior/(double)P[L]->NumAlle[j] + P[L]->Factor * P[L]->x[j][k] : 
												P[L]->Prior + P[L]->Factor * P[L]->x[j][k] ;
			P[L]->xDPsum[j] += P[L]->xDP[j][k];
												
			/* set the baseline Dirichlet pars for computing using the WORST method */
			P[L]->xBP[j][k] = P[L]->ScalePrior  ?  P[L]->Prior/(double)P[L]->NumAlle[j] + P[L]->x[j][k] : 
													P[L]->Prior + P[L]->x[j][k];
												
		}
		printf("ALLE_COUNTS : Locus %d : %d Alleles : ", j+1,P[L]->NumAlle[j]);
		for(k=0;k<P[L]->NumAlle[j];k++)  {
			printf("%f ",P[L]->x[j][k]);
		}
		printf("\n");
		
		printf("ALLE_FREQUENCIES : Locus %d : %d Alleles : ", j+1,P[L]->NumAlle[j]);
		for(k=0;k<P[L]->NumAlle[j];k++)  {
			P[L]->p[j][k] = P[L]->x[j][k] / P[L]->xsum[j];  /* compute the allele frequencies */
			printf("%f ",P[L]->p[j][k]);
		}
		printf("\n");
	
		printf("SIMULATION_DIRICHLET_PARS : Locus %d : %d Alleles : ", j+1,P[L]->NumAlle[j]);
		for(k=0;k<P[L]->NumAlle[j];k++)  {
			printf("%f ",P[L]->xDP[j][k]);
		}
		printf("\n");
		
		printf("BASELINE_DIRICHLET_PARS : Locus %d : %d Alleles : ", j+1,P[L]->NumAlle[j]);
		for(k=0;k<P[L]->NumAlle[j];k++)  {
			printf("%f ",P[L]->xBP[j][k]);
		}
		printf("\n");
	
	}
	fclose(in);
}




/*
Simple function that returns the index of a given population name within
and IndCollection
or returns -1 if that name does not exist
*/
int IdxOfPopName(struct IndCollection *B, char *Name)
{
	int ret = -1;
	int k;
	if(B==NULL) {
		fprintf(stderr,"Null pointer B while checking for name %s in function IdxOfPopName\nExiting....\n\n",Name);
		exit(1);
	}
	for(k=0;k<B->NumPops;k++)  {
/*printf("       Seeking: \"%s\" from amongst  \"%s\"  strcmp value = %d\n",Name, B->PopNames[k],strcmp(Name,B->PopNames[k])); */
		if(strcmp(Name,B->PopNames[k])==0) {
			ret = k;
			break;
		}
	}
	return(ret);
}


/* here is a simple function to return the index of name from amongst an array of n strings. */
int IdxOfString(char *StrArray[], int n, char Name[])
{
	int ret = -1;
	int k;
	if(StrArray==NULL) {
		fprintf(stderr,"Null pointer StrArray while checking for name %s in function IdxOfPopName\nExiting....\n\n",Name);
		exit(1);
	}
	for(k=0;k<n;k++)  {
		if(strcmp(Name,StrArray[k])==0) {
			ret = k;
			break;
		}
	}
	return(ret);
	
}


/* Once the populations have been summarized in the Baselines File you can call this function
 with the name of a file that holds the reporting unit info and it will parse it, check to make
 sure the repunit names are unique for each one and that the populations are found.  Also it
 checks to make sure that each population is in a repunit. 
 
 The RepUFile should look like this:
 
 REPUNIT  Repunit 1
    Pop_1
	Pop_2
 
 REPUNIT  Repunit 2
   Pop3
   Pop8
 
 etc.
 
REPUNIT is a reserved word in that file.
 
*/

reporting_unit_info* GatherReportingUnitData(char *RepUFile, struct IndCollection* B)
{
	FILE *file;
	int n=-1,i,j;
	int grab_name=0;
	char tempstr[1000];
	reporting_unit_info *ret = (reporting_unit_info *)malloc(sizeof(reporting_unit_info));
	int Orphans = 0;
	
	ret->NumPopsInRepUnit = (int *)calloc(MAX_POPS,sizeof(int));
	
	
	if( (file=fopen(RepUFile,"r"))==NULL) {
		fprintf(stderr,"Failed to open file %s to read reporting units.\nExiting....\n\n",RepUFile);
		exit(1);
	}
	
	/* take a first pass just to count stuff up */
	while(!feof(file)) {
		fscanf(file," %s",tempstr);
		if(n==-1 && strcmp(tempstr,"REPUNIT")!=0) {
			fprintf(stderr,"Expecting file %s to start with REPUNIT identifier, not %s. \nExiting\n\n",RepUFile,tempstr);
			exit(1);
		}
		if(strcmp(tempstr,"REPUNIT")==0) {
			n++; 
			grab_name=1;
		}
		else {
			if(grab_name==1) {
				/*printf("%d  %s\n",n,tempstr);*/
				grab_name=0;
			}
			else {
				if(!feof(file)) {
					ret->NumPopsInRepUnit[n]++;
				}
			}
		}
	}
	
    /* record number of rep units read and do some memory allocation */
	ret->NumRepUnits = n+1;
	ret->Names = (char **)calloc(ret->NumRepUnits,sizeof(char *));
	ret->PopsInRepUnits = (int **)calloc(ret->NumRepUnits,sizeof(int *));
	ret->RepUnitOfPops = (int *)calloc(B->NumPops,sizeof(int));

	
	printf("REPORTING_UNITS:  NumRepUnits : %d\n",ret->NumRepUnits);
	
	/* now, rewind the file and do a second pass to grab all the stuff and check it */
	rewind(file);
	
	
	for(i=0;i<ret->NumRepUnits;i++)  {
		ret->Names[i] = (char *)calloc(1000,sizeof(char));  /* allocate to hold name */
		ret->PopsInRepUnits[i] = (int *)calloc(ret->NumPopsInRepUnit[i],sizeof(int));
		
		/* eat REPUNIT */
		fscanf(file," %s",tempstr);
		if(strcmp(tempstr,"REPUNIT")!=0) {
			fprintf(stderr,"Whoa! Expecting to Read REPUNIT, but read %s parsing reporting unit file.\nExiting\n\n",tempstr);
			exit(1);
		}
		
		/* get the repunit name */
		fscanf(file," %s",tempstr);
		printf("REPORTING_UNITS:  RepUnit with Index %d is %s\n",i,tempstr);
		
		/* then check to make sure it is not the same as one given before */
		if(i>0) { int check;
			if( (check = IdxOfString(ret->Names, i,tempstr)) > -1) {
				fprintf(stderr,"Error Reading Reporting units file.  Repunit %s with index %d appears to have the same name as the one with index %d.\nExiting\n\n",tempstr,i,check);
				exit(1);
			}
		}
		
		/* record the repunit name */
		sprintf(ret->Names[i],"%s",tempstr);
		
		/* then cycle over the populations and get each of them */
		for(j=0;j<ret->NumPopsInRepUnit[i];j++)  { int tempidx;
			fscanf(file," %s",tempstr);
			
			
			/* check to make sure this population exists in the baseline */
			tempidx = IdxOfPopName(B, tempstr);
			if(tempidx==-1) {
				fprintf(stderr,"Error! Bad news.  You requested population \"%s\" to be in repunit \"%s\", but that population is not found in the baseline!\nExiting\n\n",tempstr,ret->Names[i]);
				exit(1);
			}
			/* if we got here, then record that index */
			ret->PopsInRepUnits[i][j] = tempidx;
			
			
			printf("REPORTING_UNITS:            Population %s with index %d is population index %d in repunit %d\n",tempstr,tempidx,j,i);
			
			/* and also record the repunit of the population */
			B->PopSumms[tempidx]->RepUnitIdx = i;
			ret->RepUnitOfPops[tempidx] = i;
		}
	}
	
	/* now that this is done, we just have to look at all the different populations and make sure
	 that they have been assigned a reporting unit  */
	for(i=0;i<B->NumPops;i++)  {
		if(B->PopSumms[i]->RepUnitIdx == -1) {
			fprintf(stderr,"Bad News Dude! Population %s is not assigned to a reporting unit in the file %s.  This is a fatal error.\n",B->PopSumms[i]->Name, RepUFile);
			Orphans++;
		}
	}
	if(Orphans) {
		fprintf(stderr,"ERROR!  %d populations are not in Reporting Groups.\nExiting\n\n",Orphans);
		exit(1);
	}
	
	printf("REPORTING_UNITS:  Apparently every population in the baseline was found in a reporting unit.  Good.\n");
	
	return(ret);
}

int GetGSI_Options(pop_struct **P, int argc, char *argv[], double **TruePi, int *INDS, int *MIXES, 
	int *MIXSIZE, Method *Meth, struct IndCollection **Blines, struct IndCollection **Mixture, struct IndCollection **Holdout, int *SelfAss, int **FixedNums, int *NoEM, int *ProperBayesWayReps,
	int *NumMixtureLoglSimReps, int *NumBaselinLoglSimReps, int *DumpMixLoglFile, int *DumpBaselineLoglFile, BayesianOpts **BayesOpts, int **PostPredNs, int *NumPredPostNs, char ***PostPredRanFiles, int *NumPredRanFiles, char ***PostPredRanSuffixes,
                   int ***MultiFixMix, char ***MultiFixMixNames, int *NumMultFixMixes)
{
	int i;
	char locfilename[10000];
	char baselinefilename[10000];
	char mixturefilename[10000];
	char holdoutfilename[10000];
	struct IndCollection *Baselines;  /* this is a local variable---an address we may send back in Blines */
	struct IndCollection *TheMixture;  /* this is a local variable---an address we may send back in Mixture */
	struct IndCollection *TheHoldout;  /* this is a local variable---an address we may send back in Holdout */
	BayesianOpts *BO; /* this is a local variable.  An address we can send back in BayesOpts */  
	int *LocalPloidies = NULL;
	int LOO_style; /* the leave one out style.  0 = don't leave out and 1 = do so */
	int Samp_Unit; /* for level of resampling. 0 = genecopies; 1 = single locus genotypes; 2 = multilocus genotypes */
	int LocFileF=0,
		ploidyF = 0,
		PiF=0,
		endpopF=0,
		NameF=0,
		FactorF=0,
		NoScaleF=0,
		PriorF=0,
		indsim_f=0,
		mixwts_f=0,
		mixfishery_f=0,
		quasiuncond_f=0,
		constelfst_f=0,
		resample_baselines_f=0,
		resample_adfg_f=0,
		sampunit_f=0,
		loostyle_f=0,
		loss_mat_f = 0,
		baseline_f = 0,
		repunit_f = 0,
		mixture_f = 0,
		holdout_f = 0,
		selfass_f = 0,
		loss_vec_f = 0,
		fixednums_f=0,
		fixedpiF=0,
		resamp_like_vecs_f=0,
		resamp_genes_loo_f=0,
		resamp_loci_loo_f=0,
		noem_f=0,
		nomixcdf_f=0,
		print_genos_f=0,
		num_mix_logl_reps_f=0,
		num_base_logl_reps_f=0,
		num_mcmc_sweeps_f=0,
		num_mcmc_burn_in_f=0,
		pi_trace_interval_f=0,
		zsum_trace_interval_f=0,
		indiv_trace_interval_f=0,
		post_pred_ns_f=0,
		post_pred_rando_samp_ns_f=0,
		baseline_close_matchers_f=0,
    multi_fix_f=0,
		mixture_close_matchers_f=0;
	int CurrentPop=0; 
	DECLARE_ECA_OPT_VARS;


	SET_OPT_WIDTH(28);
	SET_ARG_WIDTH(17);
	SET_VERSION("\ngsi_sim -- simulate genetic stock ID --\n\nVERSION: 1.0\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 13 August 2005\nCOPYRIGHT: None -- US Govt Employee Work\n\n");
	SET_DESCRIBE("\ngsi_sim -- simulate genetic stock ID --\n");

	/* allocate memory to the first population, and do other initializations that we must do  */
	P[0] = InitPopStruct(CurrentPop);
	BO = (BayesianOpts *)malloc(sizeof(BayesianOpts));
	BO->NumBurnIn = 50000;
	BO->NumSweeps = 0;
	BO->PiTraceInterval = 0;
	BO->IndivTraceInterval = 0;
	BO->ZSumTraceInterval = 0;
	*BayesOpts = BO;
	(*PostPredNs)=NULL;
	*NumPredPostNs=0;
	
	
	BEGIN_OPT_LOOP 	 
		
		if(__OptCommentLevel > 0) {
			printf("\n   ****  Command Line Switches Controlling gsi_sim AlleFreq Input / Population Characteristics  ****\n\n");
		}
		
		if(MULT_USE_OPTION(NameF,
			n,
			name,
			[S],
			[S]=the name of the current population,
			[S]=the name of the current population.  [S] should have no spaces.
				The default name is NoName,
			MAX_POPS)) {
			
			if(ARGS_EQ(1)) {
				GET_STR(P[CurrentPop]->Name);
				printf("POPNAME: Population number %d named %s\n",CurrentPop+1,P[CurrentPop]->Name);
			}
		}
		
		
		if(MULT_USE_OPTION(PriorF,
			p,
			prior,
			[R],
			Make the prior for allele freqs [R] or [R]/K,
			Issuing this option makes the prior for allele freqs [R] for the current 
				population if the -s or --scale-prior option is not given.  If the 
				-s or --scale-prior option is given then the prior for allele freqs
				will be [R]/K where K is the number
				 of alleles at the locus.  The default is 1.0,
			MAX_POPS)) {
			
			if(ARGS_EQ(1)) {
				P[CurrentPop]->Prior = GET_DUB;
				printf("FREQ_PRIOR: Population number %d gets alle freq prior %f\n",CurrentPop+1,P[CurrentPop]->Prior);
			}
		}

		
		if(MULT_USE_OPTION(FactorF,
			f,
			sample-size-factor,
			[R],
			Set Dirichlet pars to allele counts times [R],
			Allele frequencies for each population are simulated from something
				akin to the predictive posterior for them given the baselines.  Using this
				option causes the Dirichlet parameters defining that predictive-posterior-like
				distribution to be multiplied by [R] for the current population.  Default is 1.0.,
			MAX_POPS)) {
			
			if(ARGS_EQ(1)) {
				P[CurrentPop]->Factor = GET_DUB;
				printf("SAMPLE_SIZE_FACTOR: Population number %d gets sample size factor %f\n",CurrentPop+1,P[CurrentPop]->Factor);
			}
		}
		
		
		if(MULT_USE_OPTION(NoScaleF,
			s,
			no-scale-prior,
			,
			Make the prior for allele freqs p instead of default p/K,
			Issuing this option makes the prior for allele freqs p for the current 
				population where
				 p is the argument to the --prior or -p option.  The default is for the prior to be p/K for
				 each allele where K is the number
				 of alleles at the locus,
			MAX_POPS)) {
			
			if(ARGS_EQ(0)) {
				P[CurrentPop]->ScalePrior = 0;
				printf("SCALE_PRIOR: Population number %d given the scale-prior option.\n",CurrentPop+1);
			}
		}
		
		if(MULT_USE_OPTION(
			constelfst_f,
			c,
			constellation-fst,
			[R1] [R2] ...,
			NOT IMPLEMENTED defines  constellation pops with Fst = R1 R2 ...,
			NOT IMPLEMENTED. Use this option to define new populations that are all part of a cluster that the current
			population belongs to.  R1 is the Fst value between the current population and the first constellation
			population.  R2 is the Fst value between the current population and the second constellation 
			population.  This can be used only once for each population.  Each time a mixed fishery sample 
			is simulated the allele frequencies in these constellation populations are set by simulating allele
			frequencies q from their posterior distribution given the baseline sample for the current population.  Then 
			the frequencies for the constellation populations are simulated from a Dirichlet distribution with parameters
			of q times the quantity 1/Fst  - 1.  Individuals simulated from the current population then have a chance that 
			they actually came from one of the constellation populations.  The default is for any individual from the current
			population to have come from any of the constellation populations or the original one with equal probability.  This
			behavior can be changed with the -w or --constellation-wts option.  Note that for any population the -c or the --constellation-wts option
			MUST APPEAR BEFORE THE -w or --constellation-wts option., 
			MAX_POPS)) { 
				if(P[CurrentPop]->NumConstell > 0) {
					fprintf(stderr,"Error! It appears that -c or --constellation-fst has been issued more than once for the same population (population %d)\n",CurrentPop+1);
					OPT_ERROR;
				}
				else if(ARGS_GT(0)) { int i;
					P[CurrentPop]->NumConstell = COUNT_ARGS;
					P[CurrentPop]->ConstelFst = (double *)ECA_CALLOC(P[CurrentPop]->NumConstell, sizeof(double));
					for(i=0;i<P[CurrentPop]->NumConstell;i++)  {
						P[CurrentPop]->ConstelFst[i] = GET_DUB;
					}
				}
		}
		
		if ( MULT_USE_OPTION(
			LocFileF,
			a,
			allele-freq-file,
			[S],
			S=pathname to file with locus information,
			S=pathname to file with locus information.  The format of this allele-frequency information
					file is simple.  It includes only integers or real numbers (or comments enclosed by a pair of & characters).
					The first integer is the number of loci in the file.  Then for each locus you give the
					number of alleles at that locus followed by the counts (possibly real-valued) of each allele at the locus.
					All these things should be separated by whitespace.  No punctuation!  If the number of
					allele counts listed does not match the number of alleles given for the locus then the program
					may fail without warning or give otherwise unexpected results.  Each time this command is executed
					the population counter is incremented by one.  So any commands previous to the current --freq-file option
					but after the last --freq-file option apply to the population with allele counts given in the current
					--freq-file option., 
					MAX_POPS) ) {
			if(  ARGS_EQ(1)  && !HAS_BUT_SHOULD_NOT(mixwts_f,"-m or the --mixture-wts")  && 
								!HAS_BUT_SHOULD_NOT(loss_mat_f,"-L or the --loss-matrix")) {
				if(endpopF) {
					fprintf(stderr,"Error! You have mixed the -a/--allele-freq-file and the --end-pop option.  This is not allowed.  Exiting...\n");
					exit(1); 
				}
				GET_STR(locfilename);
				ReadLocFile(P,CurrentPop,locfilename);
				
				/* increment the population counter, and allocate to the next population */
				CurrentPop++;
				P[CurrentPop] = InitPopStruct(CurrentPop);
			}
		}
		if(OPTION(
			ploidyF,
			,
			ploidies,
			J1 ... JL,
			List of ploidies for loci 1 through L.  For -a baseline freq spec only,
			When you are just going to be doing simulation from allele frequencies specified in files included 
			with the -a option you can set the ploidy of each locus with this option.  This option must be issued
			after at least one -a option.  It takes L arguments where L is the number of loci in each of the 
			allele freq files.  This is optional.  The default is to consider all loci as diploid. ))
		{
			if(ALREADY_HAS(LocFileF,-a or --allele-freq-file) ) {
				if(ARGS_EQ(P[0]->NumLoc))  { int j;
					LocalPloidies = (int *)ECA_CALLOC(P[0]->NumLoc, sizeof(int));
					for(j=0;j<P[0]->NumLoc;j++)  {
						LocalPloidies[j] = GET_INT;
					}
				}
			}
		}
			
		if ( MULT_USE_OPTION(
			endpopF,
			,
			end-pop,
			,
			advance focus to another population,
			Advances focus to another population.  After issuing this command the options -n -p -f -s and --Pi will apply
					to another population in memory.  This option should also be issued after the last set of -n -p -f and 
					-s options.  This option cannot be mixed with the -a or --allele-freq-file option.  The purpose of this option
					is to allow the user to enter population specific parameters that will be linked to the name of the
					population.  This is useful if the baseline is read from a genotypes file and you want to add more parameters
					to each separate population.  Note that the argument of the -n option---i.e. the name---must match a name
					given in the baseline file or the program will exit when it searches for that name.  Example: imagine that 
					the baseline file has among others the populations Klamath_F and Feather_H in it.  To give those both weight 5 
					in the simulated mixture and to set their 
					scale priors -s to .7 you would do -n Klamath_F -s .7 --Pi 5 --end-pop -n Feather_H -s .7 --Pi 5. --end-pop , 
					MAX_POPS) ) {
			if(  ARGS_EQ(0)) {
				if(LocFileF) {
					fprintf(stderr,"Error! You have mixed the -a/--allele-freq-file and the --end-pop option.  This is not allowed.  Exiting...\n");
					exit(1); 
				}
				CurrentPop++;
				P[CurrentPop] = InitPopStruct(CurrentPop);
			}
		}
		
		if(__OptCommentLevel > 0) {
			printf("\n   ****  Command Line Switches Controlling gsi_sim Indiv-Data Output / Mixture Characteristics  ****\n\n");
		}
		if(OPTION(
			baseline_f,
			b,
			baseline-genotypes,
			[S],
			path to file with genotypes of baseline,
			)) {
			if(ARGS_EQ(1)) {
				/* do the memory allocations for the baselines collection of individuals */
				Baselines=InitIndCollection();
				
				/* if the mixture has already been read in, then transfer the AlleInfo from there */
				if(mixture_f) {
					Baselines->AlleInfo = (*Mixture)->AlleInfo;
				}
				GET_STR(baselinefilename);
				ReadIndsFromFile(baselinefilename,Baselines);
				
				/* then assign the value of Baselines to *Blines */
				*Blines = Baselines;
				
			}
		}
	if(OPTION(
			  repunit_f,
			  r,
			  rep-unit-file,
			  S,
			  path to a file with reporting units,
			  S is the path to a file to reporting units which has this format:\n\n
REPUNIT NameOfRepUnitOne\n
Pop10\n
Pop2\n
REPUNIT NameOfRepUnitTwo\n
Pop3\n
Pop4\n
etc.\n\n
Every population must be put into a reporting unit. The program will exit with an error if
	   you specify a population that does not appear in the baseline. This can only be given after the 
	   -b option has been given.)) {
		if(ALREADY_HAS(baseline_f,"-b/--baseline-genotypes") && ARGS_EQ(1)) {
			GET_STR(gRepUnitFile);
			gHasRepUnits=1;
			
		}
			
	}
		if(OPTION(
			mixture_f,
			t,
			mixture-genotypes,
			[S],
			path to file with genotypes of the mixture to do inference about,
			)) {
			if(ARGS_EQ(1)) {
				/* do the memory allocations for the baselines collection of individuals */
				TheMixture=InitIndCollection();
				
				/* if the baselines have already been read in, then transfer the AlleInfo from there */
				if(baseline_f) {
					TheMixture->AlleInfo = (*Blines)->AlleInfo;
				}
				GET_STR(mixturefilename);
				ReadIndsFromFile(mixturefilename,TheMixture);
				
				/* then assign the value of Baselines to *Blines */
				*Mixture = TheMixture;
			}
		}
		
		
		if(OPTION(
			holdout_f,
			h,
			holdout-genotypes,
			F,
			path to file with genotypes in a holdout set,
			This should be the path to a file that has the genotypes of individuals that are in a holdout set.  The format for this file
			is identical to that of the Baseline file.  Note that if you have included these individuals also in the baseline and so 
			want to leave each of their genotypes out of the baseline in turn when computing genotype probs then you will need to ensure that
			the population names used in the holdout file are exactly the same as those used in the baseline file (when the holdout file includes
			any populations that are also in the baseline---as it ought to if it is a holdout file!)  This option can only be issued after the baseline
			file and is currently incompatible with the --mixture-genotypes option. 
			)) {
			if(ARGS_EQ(1)) {
				if(ALREADY_HAS(baseline_f,-b/--baseline-genotypes)) {
					/* do the memory allocations for the baselines collection of individuals */
					TheHoldout=InitIndCollection();
					
					TheHoldout->AlleInfo = (*Blines)->AlleInfo;
					
					GET_STR(holdoutfilename);
					ReadIndsFromFile(holdoutfilename,TheHoldout);
					
					/* then assign the value of TheHoldout to the output variable */
					*Holdout = TheHoldout;
				}
			}
		}
		
		
		if(OPTION(
			selfass_f,
			,
			self-assign,
			,
			do self-assignment on the baseline genotypes using Leave-One-Out,
			Issuing this option when the -b option has been invoked will cause the program to 
			do self-assignment of the individuals in the baseline much the way that GeneClass does.
			If this option is given but the -b option was not given then the program will abort.  Output
			from this option is on lines tagged by: SELF_ASSIGN_A_LA_GC_CSV))  {
		
			if(ARGS_EQ(0)) {
				*SelfAss = 1;
			}
		}
	
		if(OPTION(
			num_mix_logl_reps_f,
			,
			mix-logl-sims,
			J1 J2,
			J1 is number of reps for simulating log-likelihoods of mixture fish. J2 nonzero means print histograms to file,
			Default is zero. If J1 > 0 then\054 for every fish in the mixture\054 J1 new genotypes will be simulated from 
		      the population to which the fish is assigned (via the Full-EM method) and the loglikelihood of these simulated
			  fish will be compared to the log likelihood of the fish in the mixture.  Z scores and quantiles are reported in 
		      standard output on tagged lines.  If J2 is nonzero\054 then the full histogram output is dumped to the file GSI_MixtureLoglCompare_DumpFile.txt. 
			  Beware! These files can be quite large.  Note: when the genotypes are simulated\054 the pattern of missing data in each fish is duplicated. )) {	
			if(ARGS_EQ(2)) {
				*NumMixtureLoglSimReps = GET_INT;
				*DumpMixLoglFile = GET_INT;
			}
		}
	
	if(OPTION(
			  num_base_logl_reps_f,
			  ,
			  base-logl-sims,
			  J1 J2,
			  J1 is number of reps for simulating log-likelihoods of baseline fish. J2 nonzero means print histograms to file,
			  Default is zero. If J1 > 0 then\054 for every fish in the baseline\054 J1 new genotypes will be simulated from 
		      the reference population that fish comes from and the loglikelihood of these simulated
			  fish will be compared to the log likelihood of the fish in the baseline.  Z scores and quantiles are reported in 
		      standard output on tagged lines.  If J2 is nonzero\054 then the full histogram output is dumped to the file GSI_BaselineLoglCompare_DumpFile.txt. 
			  Beware! These files can be quite large.  Note: when the genotypes are simulated\054 the pattern of missing data in each fish is duplicated. )) {	
		if(ARGS_EQ(2)) {
			*NumBaselinLoglSimReps= GET_INT;
			*DumpBaselineLoglFile = GET_INT;
		}
	}
	if(OPTION(baseline_close_matchers_f,
			  ,
			  close-match-base,
			  J1 J2,
			  Find closely matching pairs of genotypes in the baseline THEN EXIT,
			  Find pairs of individuals in the baseline that have only J1 (or fewer) mismatching genotypes out of at least J2 loci at which they both have genotypes that are not missing data.
			  Invoking this option will cause the 
			  program to quit execution before carrying out any of the GSI analyses.  The output goes to a file called close_matches_baseline.txt.  Requires invoking the -b/--baseline-genotypes option first.  )) {
		if(ALREADY_HAS(baseline_f, "-b/--baseline-genotypes") && ARGS_EQ(2)) {
			gBaselinePairMismatchLimit = GET_INT
			gBaselinePairLocusLimit = GET_INT; 
		}
	}
	if(OPTION(mixture_close_matchers_f,
			  ,
			  close-match-mix,
			  J1 J2,
			  Find closely matching pairs of genotypes in the mixture THEN EXIT,
			  Find pairs of individuals in the mixture that have only J1 (or fewer) mismatching genotypes out of at least J2 loci at which they both have genotypes that are not missing data.
			  Invoking this option will cause the 
			  program to quit execution before carrying out any of the GSI analyses.  The output goes to a file called close_matches_mixture.txt.  Requires invoking the -t/--mixture-genotypes option first.  )) {
		if(ALREADY_HAS(baseline_f, "-t/--mixture-genotypes") && ARGS_EQ(2)) {
			gMixturePairMismatchLimit = GET_INT
			gMixturePairLocusLimit = GET_INT; 
		}
	}
	
	if(__OptCommentLevel > 0) {
		printf("\n   ****  Command Line Switches Controlling gsi_sim Bayesian method  ****\n\n");
	}
		
		if(OPTION(
			num_mcmc_sweeps_f,
			,
			mcmc-sweeps,
			J,
			Number of sweeps to do after Burn-In in Bayesian analysis,
			Number of sweeps to do after Burn-In in Bayesian analysis. Default is 0 --- i.e. 
		   do no Bayesian analysis. ))  {
		
			if(ARGS_EQ(1)) {
				BO->NumSweeps = GET_INT;
			}
		}
		if(OPTION(
			  num_mcmc_burn_in_f,
			  ,
			  mcmc-burnin,
			  J,
			  Number of sweeps of Burn-In to do in Bayesian analysis,
			  Number of sweeps of Burn-In to do in Bayesian analysis. Default is 50 thousand. However\054 if 
			    mcmc-sweeps is not given\054 or is zero\054 then no mcmc will be done\054 regardless of the
		        burn in value.))  {
		
			if(ARGS_EQ(1)) {
				BO->NumBurnIn = GET_INT;
			}
		}
	if(OPTION(
			  pi_trace_interval_f,
			  ,
			  pi-trace-interval,
			  J,
			  Thinning interval at which to print trace of mixing proportions from mcmc,
			  Thinning interval at which to print trace of mixing proportions from mcmc.  If this is 0\054
	           then there will be no trace printed.  Otherwise\054 the current values of the MCMC will be 
			   printed every J iterations\054 from the beginning (i.e. including the burn-in). Default is zero. ))  {
		
		if(ARGS_EQ(1)) {
			BO->PiTraceInterval = GET_INT;
		}
	}
	if(OPTION(
			  zsum_trace_interval_f,
			  ,
			  zsum-trace-interval,
			  J,
			  Thinning interval at which to print trace of ZSums from mcmc,
			  Thinning interval at which to print trace of ZSums from mcmc.  If this is 0\054
			  then there will be no trace printed.  Otherwise\054 the current values of the MCMC will be 
			  printed every J iterations\054 from the beginning (i.e. including the burn-in). Default is zero.
			  The ZSum is the sum over all fish of the current allocation of the fish to each population or reporting
			  group.  It is the number of fish currently allocated to each population (or reporting group) done by randomly 
			  assigning them according to their posterior probabilities. ))  {
		
		if(ARGS_EQ(1)) {
			BO->ZSumTraceInterval = GET_INT;
		}
	}
	if(OPTION(
			  indiv_trace_interval_f,
			  ,
			  ind-trace-int,
			  J,
			  Thinning interval at which to print trace of individual PofZs,
			  Thinning interval at which to print trace of posterior probabilities from individual fish from mcmc.  If this is 0\054
			  then there will be no trace printed.  Otherwise\054 the current values of the MCMC will be 
			  printed every J iterations\054 from the beginning (i.e. including the burn-in).  Note that this option could
			  produce massive amounts of output.  Default is 0. APPARENTLY NOT YET IMPLEMENTED.))  {
		
		if(ARGS_EQ(1)) {
			BO->IndivTraceInterval = GET_INT;
		}
	}
	if(OPTION(
			  post_pred_ns_f,
			  ,
			  post-pred-ns,
			  J1 J2 ...,
			  Size of samples to draw population allocations from the posterior predictive,
			  This is a little something to generate some realizations that should be helpful for propagating uncertainty
			  in fisheries management.  The idea is that we might have 400 fish genotyped\054 but we know that 500 were sampled
	          from the fishery.  If we want to have an idea of the posterior distribution of all 500 of those fish we can simulate
			  the origins for 100 new fish from the posterior predictive distribution and add those to the current allocations of the 
	          first 400.  With this option we tell gsi_sim to print these out at the same frequency as the zsum-trace-interval.  You can 
	          specify up to 50 Js (i.e. different sample sizes to draw from the posterior predictive).  The output will be exactly parallel
	          to the output in the zsum files.  The will be named pop_pred_post_draws_XX.txt or rep_unit_post_pred_draws_XX.txt where XX is the 
	          size of the sample (i.e. one of J1 or J2 or J3\054 etc.)  Note that the Js must be unique or there will be an error.  ) ) {  
	   int NumArgs; int i; int j;
	   
		NumArgs = COUNT_ARGS;
		printf("POSTPREDDEBUG: NumArgs is %d\n",NumArgs);
		
		if(NumArgs>50) {
			fprintf(stderr,"Error! No more than 50 arguments are allowed to --post-pred-ns.\nExiting\n\n");
			exit(1);
		}
		if(ARGS_GEQ(1)) { 
			(*PostPredNs) = (int *)calloc(NumArgs,sizeof(int));
			for(i=0;i<NumArgs;i++)  {
				(*PostPredNs)[i] = GET_INT;
				
				printf("POSTPREDDEBUG: i= %d  Ns= %d\n",i,(*PostPredNs)[i]);
				for(j=0;j<i;j++) {  /* make sure none are duplicated */
					if((*PostPredNs)[i]==(*PostPredNs)[j]) {
						fprintf(stderr,"Error!  post-pred-ns argument %d is identical in value (%d) to argument %d\nExiting\n\n",i+1,(*PostPredNs)[i],j+1);
						exit(1);
					}
				}
			}
			*NumPredPostNs = NumArgs;
		}	
	}
	if(MULT_USE_OPTION(
			  post_pred_rando_samp_ns_f,
			  ,
			  post-pred-ran-ns,
			  S1 S2,
			  Variable size of samples to draw population allocations from the posterior predictive,
			  This is an extension of the post-pred-ns option.  In this case\054 we imagine that we are uncertain
			  about how many fish were actually sampled in the fishery in the first place. For example\054 maybe we
			  have a posterior sample of the number sampled in the fishery and want to incorporate that uncertainty.
			  With this option you can include a sample from the posterior in a file that contains nothing but integers
			  that are that sample.  S1 is the name you want for the YYY part int the output file names pop_pred_post_draws_YYY.txt
			  and rep_unit_post_pred_draws_YYY.txt.   S2 is the path to the file that includes the sample from the
			  posterior.  This option can be used MAX_POST_PRED_RAN_FILES times so that you can expand up to different sample sizes
			  that might reflect different posteriors etc.  MAX_POST_PRED_RAN_FILES is set to 50 at the moment.,
			  MAX_POST_PRED_RAN_FILES
  ) ) {  
	   int i;
	   if(ARGS_EQ(2)) {
		   i = (*NumPredRanFiles)++;
		   if(i==0) {
			   (*PostPredRanFiles) = (char **)calloc(MAX_POST_PRED_RAN_FILES, sizeof(char *)); 
			   (*PostPredRanSuffixes) = (char **)calloc(MAX_POST_PRED_RAN_FILES, sizeof(char *));
		   }
		   (*PostPredRanFiles)[i] = (char *)calloc(MAX_FILE_NAME_LENGTH, sizeof(char));
		   (*PostPredRanSuffixes)[i] = (char *)calloc(MAX_FILE_NAME_LENGTH, sizeof(char));
		   GET_STR((*PostPredRanSuffixes)[i]);
		   GET_STR((*PostPredRanFiles)[i]);
	   }
	}
	

		
		
		if(__OptCommentLevel > 0) {
			printf("\n   ****  Command Line Switches Controlling gsi_sim Data Output / Mixture Characteristics  ****\n\n");
		}
		
		if(OPTION(
			indsim_f,
			i,
			ind-sim,
			[J],
			simulate [J] individuals from each population,
			This option controls how many individuals should be simulated.  In this case [J] individuals
			will be simulated from each of the populations.  Output is on lines tagged by the phrase SIMPLE_IND_SIM which contain the
			likeihood for each simulated fish of coming from each of the different populations.  These likelihoods have been scaled so
			so that they sum to one.   If the -m or --mixture-wts option then there is also a line tagged with the phrase
			TRUE_WEIGHTED_IND which gives the posterior probability of each fish coming from each of the populations given
			the true mixing proportions.  All ind-simmed individuals appear in the output before any mix-sim individuals.) ) {
				if(ARGS_EQ(1)) {
					*INDS = GET_INT;
				}
			
		}
		if ( MULT_USE_OPTION(
			PiF,
			,
			Pi,
			[R],
			Relative weight of focal population in simulated mixtures,
			Sets the relative weight of the population having name given by -n in the current focal population to Pi in 
			the simulated mixtures of the populations.  This option clashes with -m/--mixture-weights.  You can set the 
			mixing proportions using -m only if you are using the -a options for allele frequencies. It is so cumbersome 
			when using the -b option for baselines that it is not allowed.  You can use the --Pi option to set mixing 
			proportions under either the -a or the -b approach to allele freq input., 
					MAX_POPS) ) {
			if(  ARGS_EQ(1)) {
				if(mixwts_f || fixednums_f) {
					fprintf(stderr,"Error! You have mixed either the -m/--mixture-weights or the --fixed-nums options and the --Pi option.  This is not allowed.  Exiting...\n");
					exit(1); 
				}
				if(fixedpiF) {
					fprintf(stderr,"Error! You have mixed either the --Pi option and the --fixed-Pi option.  This is not allowed.  Exiting...\n");
					exit(1); 
				}
				P[CurrentPop]->SimPi = GET_DUB;
			}
		}
		if ( MULT_USE_OPTION(
			fixedpiF,
			,
			fixed-Pi,
			[J],
			Absolute count of the focal population in simulated mixtures,
			Sets the absolute number of individuals from the population having name given by -n in the current focal population to J in 
			the simulated mixtures of the populations.  This option clashes with -m/--mixture-weights and with --fixed-nums.  You can set the 
			mixing proportions/counts using -m or --fixed-nums only if you are using the -a options for allele frequencies. It is so cumbersome 
			when using the -b option for baselines that it is not allowed.  You can use the --fixed-Pi option to set fixed mixing 
			counts under either the -a or the -b approach to allele freq input., 
					MAX_POPS) ) {
			if(  ARGS_EQ(1)) {
				if(mixwts_f || fixednums_f) {
					fprintf(stderr,"Error! You have mixed either the -m/--mixture-weights or the --fixed-nums options and the --fixed-Pi option.  This is not allowed.  Exiting...\n");
					exit(1); 
				}
				if(PiF) {
					fprintf(stderr,"Error! You have mixed either the --Pi option and the --fixed-Pi option.  This is not allowed.  Exiting...\n");
					exit(1); 
				}

				P[CurrentPop]->SimCnt = GET_INT;
				printf("P[CurrentPop]->SimCnt set to %d.  CurrentPop= %d\n",P[CurrentPop]->SimCnt,CurrentPop);
			}
		}
    if(MULT_USE_OPTION(
      multi_fix_f,
      ,
      multi-fix-mix,
      S J0 ... J_NumPops,
      define another mixed fishery sample with fixed numbers from each population,
      This is a fully dodgy option that I wrote to be able to do lots of mixed fishery samples of different
         compositions and different sizes while reading in the baseline only once.  It is dodgy because it does not check
         that you have the right number of populations while it is reading this in.  S is a unique string tag you want for 
         the fishery and J0 is the number from population 0\054 J1 is the number from pop 1\054 and so forth.  It is super 
                       dodgy because this does no checking to make sure you have the right number in here. If you do not have a 
                       number for every pop in the baseline then you will likely overrun an array bound., 
      MAX_MULTI_FISHERIES
      )) {
      int rem_args;
      int cc;
      if(ALREADY_HAS(baseline_f, -b or baseline-genotypes)) {
        /* allocate memory if necesssary */
        if((*MultiFixMix) == NULL) (*MultiFixMix) = (int **)calloc(MAX_MULTI_FISHERIES, sizeof(int *));
        if((*MultiFixMixNames) == NULL) (*MultiFixMixNames) = (char**)calloc(MAX_MULTI_FISHERIES, sizeof(char *));
        (*MultiFixMixNames)[*NumMultFixMixes] = (char *)calloc(MAX_MULTI_MIX_NAME, sizeof(char));
        GET_STR((*MultiFixMixNames)[*NumMultFixMixes]) /* get the name of the fishery */
        printf("Reading multi-fix-mix number %d,  Name= %s,  Numbers=  ", *NumMultFixMixes, (*MultiFixMixNames)[*NumMultFixMixes]);
        rem_args = COUNT_ARGS;
        if(rem_args != Baselines->NumPops) {
          fprintf(stderr, "Error reading arguments to --multi-fix-mix invocation number %d with Name= %s.  There should be %d remaining arguments (one for each pop in the baseline.  But there are %d\n",
                  *NumMultFixMixes, (*MultiFixMixNames)[*NumMultFixMixes], Baselines->NumPops, rem_args);
          OPT_ERROR
        }
        (*MultiFixMix)[*NumMultFixMixes] = (int *)calloc(rem_args, sizeof(int));
        for(cc=0; cc<rem_args; cc++) {
          (*MultiFixMix)[*NumMultFixMixes][cc] = GET_INT;
          printf(" %d", (*MultiFixMix)[*NumMultFixMixes][cc]);
        }
        printf("\n");
        (*NumMultFixMixes) = (*NumMultFixMixes) + 1;
      }
    }
		if(COND_REQ_OPTION(
			mixwts_f,
			m,
			mixture-wts,
			[R1] [R2] ...,
			proportions from different populations in the mixture,
			Use this option to give the true proportion of individuals from each of the populations
			in the mixture.  The number of arguments given to this function should be the number of times that
			the -a or --allele-freq-file option has been used---i.e. the number of populations.  Once this option
			has been used the -a option can no longer be used.  So use this option only after you are done inputting all
			other information about the populations.  R1 is the proportion of population 1 in the mixture.  R2 is the proportion
			of population 2 and so forth.  If the Rs do not sum to 1 they will be scaled so that they do sum to one.,
			mixfishery_f > 0  && LocFileF > 0 && !(fixednums_f) ,
			the -x or --mixed-fishery option in conjunction with the -a/--alle-freq-file option and you are not using --fixed-nums instead.
			)) {
				if(ARGS_EQ(CurrentPop) ) { int i; double normo=0.0;
				
					if(PiF) {
					fprintf(stderr,"Error! You have mixed the -m/--mixture-weights and the --Pi option.  This is not allowed.  Exiting...\n");
						exit(1); 
					}
					(*TruePi) = (double *)ECA_CALLOC(CurrentPop,sizeof(double));
					for(i=0;i<CurrentPop;i++)  { /* get the numbers off the command line */
						(*TruePi)[i] = GET_DUB;
						normo += (*TruePi)[i];
					}
					for(i=0;i<CurrentPop;i++)  { /* then normalize them so they sum to one */
						(*TruePi)[i] /= normo;
					}
				}
		}
		
		
		if(COND_REQ_OPTION(
			fixednums_f,
			,
			fixed-nums,
			J1 J2 ...,
			absolute numbers from different populations in the mixture,
			Use this option to dictate the number of individuals from each of the populations
			in the mixture.  The number of arguments given to this function should be the number of times that
			the -a or --allele-freq-file option has been used---i.e. the number of populations.  Once this option
			has been used the -a option can no longer be used.  So use this option only after you are done inputting all
			other information about the populations.  J1 is the number of population 1 individuals in the mixture.  R2 is the number
			of population 2 individuals and so forth.  The sum of these numbers will be used in place of the J2 argument to the
			-x/--mixed-fishery option.,
			mixfishery_f > 0  && LocFileF > 0 && !(mixwts_f) ,
			the -x or --mixed-fishery option in conjunction with the -a/--alle-freq-file option and you are not using -m/--mixture-wts option instead.
			)) {
				if(ARGS_EQ(CurrentPop) ) { int i; double normo=0.0;
				
					if(PiF) {
					fprintf(stderr,"Error! You have mixed the --fixed-nums and the --Pi option.  This is not allowed.  Exiting...\n");
						exit(1); 
					}
					(*FixedNums) = (int *)ECA_CALLOC(CurrentPop,sizeof(int));
					(*TruePi) = (double *)ECA_CALLOC(CurrentPop,sizeof(double));  /* note that we set TruePi here as well */
					for(i=0;i<CurrentPop;i++)  { /* get the numbers off the command line */
						(*FixedNums)[i] = GET_INT;
						(*TruePi)[i] = (double)(*FixedNums)[i];
						normo += (*TruePi)[i];
					}
					for(i=0;i<CurrentPop;i++)  { /* then normalize them so they sum to one */
						(*TruePi)[i] /= normo;
					}

				}
		}
	
		if(COND_REQ_OPTION(
			mixfishery_f,
			x,
			mixed-fishery,
			[J1] [J2],
			simulate J1 mixed fisheries of J2 individuals,
			When this option is invoked the program will simulate J1 mixed fishery samples of size
			J2 individuals each.  The number of fish in the mixture will be drawn from a mutlinomial
			distribution with cell probabilities as given in the arguments to the -m or --mixture-wts 
			option.  So if you use this option you must also issue the -m or --mixture-wts option. 
			Alternatively you can issue the --fixed-nums option.  In that case the mixtures will be composed according
			to the arguments to --fixed-nums.  It is still necessary to enter the argument J2 to the -x/--mixed-fishery
			option but its value will be ignored and the total number of fish in each mixture will be determined from
			the argument to the --fixed-nums option.
			The simulated individuals from each data set will be output on lines tagged with the phrase MIXED_FISH_INDIVS.
			The output on these lines is the plug-in estimate of the posterior probability that each fish 
			came from each of the different populations.  These are posterior probabilities based on the maximum
			likelihood estimate of the mixing proportions of each contributing stock.  The MLE of the mixing proportions
			is found using a simple EM algorithm.  For each simulated mixture the progress of the EM algorithm 
			gets reported on lines tagged by 
			EM_ALGORITHM_PROGRESS and the actual MLEs of the mixture proportions are reported on lines tagged
			with the phrase MIXFISH_PI_MLES.,
			quasiuncond_f > 0,
			"--quasi-unconditional" ))  {
				if(ARGS_EQ(2)) {
					*MIXES = GET_INT;
					*MIXSIZE = GET_INT;
				}
		}
		
		if(OPTION(
			sampunit_f,
			,
			samp-unit,
			S,
			the unit to resample from baseline or holdouts,
			The unit to resample from baselines or holdout sets.  Allowed values are: 
			genes; loci; or multilocus)) {
				if(ARGS_EQ(1)) { char temp[1000];
					 GET_STR(temp);
					 if(strcmp(temp,"genes")==0) {
						Samp_Unit = 0;
					 }
					 else if(strcmp(temp,"loci")==0) {
						Samp_Unit = 1;
					 }
					 else if(strcmp(temp,"multilocus")==0) {
						Samp_Unit = 2;
					 }
					 else {
						fprintf(stderr,"Error, unknown string argument %s to option --samp-unit.  Exiting...\n",temp);
					 }
					 
				}
		}
		
		if(OPTION(
			loostyle_f,
			,
			leave-one-out,
			S,
			set S to yes to ensure leave one out and no otherwise,
			If S is yes then the program will do leave-one-out routine if possible.  If S is no then it
			will not to a leave one out procedure\051 even if it is possible.  If this option is not given the default
			depends on whether there is holdout set.  If no holdout set then default is S=yes.  If there is a holdout set
			then S=no. )) {
				if(ARGS_EQ(1)) { char temp[1000];
					 GET_STR(temp);
					 if(strcmp(temp,"no")==0) {
						LOO_style = 0;
					 }
					 else if(strcmp(temp,"yes")==0) {
						LOO_style = 1;
					 }
					 else {
						fprintf(stderr,"Error, unknown string argument %s to option --leave-one-out.  Exiting...\n",temp);
					 }
					 
				}
		}

				
		if(OPTION(
			noem_f,
			,
			no-em,
			,
			inhibit the default behavior of doing full EM mixture analysis,
			By default when gsi_sim has received a baseline file and a mixture file the EM algorithm will be 
			applied to find the MLE of the mixing proportions in the mixture.  This can be time consuming.  Sometimes
			all that you want is to assign individuals in the mixture to different baseline populations assuming an equal
			prior weight for each population---like GeneClass does.  This is achieved by invoking this --no-em option. )) {
			if(ARGS_EQ(0)) {
				*NoEM = 1;
			}
		}

	
		if(OPTION(
			print_genos_f,
			,
			print-mix-genos,
			,
			print out the genotypes of individuals that appear in the simulated mixtures,
			If you issue this option\054 then the genotypes of individuals that get simulated in the mixture will be printed to stdout on lines
			tagged with the phrase: MIX_GENOS.)) {
				if(ARGS_EQ(0)) {
					gPrintMixGenotypes=1;
				}
		}
	
		if(OPTION(
			loss_mat_f,
			L,
			loss-matrix,
			R1 R2 ...,
			Define a matrix of loss values,
			This option lets the user input a matrix of losses.  This is a matrix with N rows and N columns where
			N is the number of populations for which genetic data are available.  The loss is the amount that it costs if you 
			get the assignment incorrect.  i.e. the ij-th element of the matrix is the price you pay for assigning a fish from
			population i to population j.  Note that the ii-th element is typically zero and the others are positive. Once this option
			has been used the -a option can no longer be used.  So use this option only after you are done inputting all
			other information about the populations.)) {
			if(ARGS_EQ(CurrentPop*CurrentPop)) { int i; int j;
				for(i=0;i<CurrentPop;i++)  {
					P[i]->Loss = (double *)ECA_CALLOC(CurrentPop, sizeof(double));
					P[i]->ExpScaledLike = (double *)ECA_CALLOC(CurrentPop, sizeof(double));
					printf("SETTING_LOSS_MATRIX_ROW : %d  :  %s  : ",i+1,P[i]->Name);
					for(j=0;j<CurrentPop;j++)  {
						P[i]->Loss[j] = GET_DUB;
						printf("%.4f ",P[i]->Loss[j]);
					}
					printf("\n");
					
				}
				
			}
		
		}
		
		if(OPTION(
			loss_vec_f,
			l,
			loss-vector,
			J1 R1 R2 ...,
			Define a vector of loss values,
			This option lets the user input a vector of losses for population J1.  This is a vector where
			N is the number of populations for which genetic data are available.  The loss is the amount that it costs if you 
			get the assignment incorrect.  i.e. the i-th element of the matrix is the price you pay for assigning a fish from
			population J1 to population i.  Note that the iJ1-th element is typically zero and the others are positive.   J1 is the index of the 
			population you want to set the loss vector for and the index starts at 1.  Once this option
			has been used the -a option can no longer be used.  So use this option only after you are done inputting all
			other information about the populations.  If this option is used more than once on a single population the last one used is the
			one that takes effect.  If this option is not called on all populations---i.e. all possible values of J1---and the --loss-matrix option
			was not given then the populations which did
			not get loss vectors explicitly assigned to them will get simple 0-1 losses---0 for the right population and 1 for the others. )) {
			if(ARGS_EQ(1+CurrentPop)) { int i; int j; int J1;
				J1 = GET_INT;
				J1--;
				for(i=0;i<CurrentPop;i++)  {
					if(P[i]->Loss == NULL)  {
						P[i]->Loss = (double *)ECA_CALLOC(CurrentPop, sizeof(double));
						if(i != J1) {
							for(j=0;j<CurrentPop;j++)  {  /* this puts in default values */
								P[i]->Loss[j] = (double)(i!=j);
							}
						}
					}
					if(P[i]->ExpScaledLike==NULL) {
						P[i]->ExpScaledLike = (double *)ECA_CALLOC(CurrentPop, sizeof(double));
					}
					
					if(i == J1) {
						printf("SETTING_LOSS_VECTOR : %d  :  %s  : ",i+1,P[i]->Name);
						for(j=0;j<CurrentPop;j++)  {
							P[i]->Loss[j] = GET_DUB;
							printf("%.4f ",P[i]->Loss[j]);
						}
						printf("\n");
					}
				}
				
			}
		
		}

		if(__OptCommentLevel > 0) {
			printf("\n\nOUTPUT OF GSI_SIM:  If...\n");
			printf("\t-b and -t are given, then the mixture will automatically be analyzed conditional on the baselines using:\n");
			printf("\t\t1. The GeneClass method (assumes equal missing proportions for all populations\n");
			printf("\t\t2. The standard GSI method (estimates mixing proportions via an EM-algorithm\n");
			printf("\t\t3. Single Pass EM (like 2, but only a single iteration of the EM is done\n");
			printf("\t\tOutput from the above come out on separate lines in a GeneClass-like CSV format\n");
			printf("\n");
			printf("\t-b is given (with or without -t) then the --self-assign option can be invoked to do\n");
			printf("\t\tself-assignment of individuals in the baselines a la GeneClass (using Leave-One-Out)\n");
			printf("\n");
			printf("\tOne of -b and -a are given, then simulated individuals may be analyzed using the -x or -i options.\n");
			printf("\tThis can be done whether or not a mixture has been specified with the -t option.\n");
			printf("\tNotes:\n");
			printf("\t\t1. If you specified baseline frequencies using -a options, then you set mixing proportions for\n");
			printf("\t\t   mixtures simulated under the -x option using the -m/--mixture-weights option.\n");
			printf("\t\t2. If you specify baseline frequencies using the -b option, then the mixing proportions for\n");
			printf("\t\t   mixtures simulated under the -x option must be specified for each named population using the\n");
			printf("\t\t   --Pi option followed by an --end-pop option.  Don't forget the --end-pop at the end!!\n");
			printf("\n");
			printf("NOTE:\n");
			printf("\t1. You cannot specify a baseline using -a options for analyzing a mixture under the -t option\n");
			printf("\t2. You cannot mix a -b option with -a options.  You must use one or the other.\n");
			printf("\n");
		}

	
	
	END_OPT_LOOP

	/* some final error checking */
	if(fixednums_f && mixwts_f) {
		fprintf(stderr,"Error! you have issued both the --fixed-nums and the -m/--mixture-wts options.  You may only issue one of them.  Exitin...\n");
		exit(1);
	}
	if(quasiuncond_f + indsim_f > 1) {
		fprintf(stderr,"Error! You have issued both the --quasi-unconditional and the -i or --ind-sim options.  The -i or --ind-sim option is not allowed in conjunction with the --quasi-unconditional option.\n Exiting...\n");
		exit(1);
	}
	if(LocFileF + baseline_f == 0)  {
		fprintf(stderr,"Error! You must issue at least one of the options -a or -b\n Exiting...\n");
		exit(1);
	}
	if(LocFileF && baseline_f)  {
		fprintf(stderr,"Error! You must issue exactly one of the options -a or -b.  Not both!!\n Exiting...\n");
		exit(1);
	}
	if(baseline_f && mixwts_f) {
		fprintf(stderr,"Error! You can't use the -m/--mixture-wts option with baselines given with the -b option.  You have to use the --Pi option.\n Exiting...\n");
		exit(1);
	}
	if(selfass_f && !baseline_f) {
		fprintf(stderr,"Error! You can't use the --self-assign option without baselines given with the -b option.\n Exiting...\n");
		exit(1);
	}
	if(baseline_f && mixfishery_f && !(PiF || fixedpiF) ) {
		fprintf(stderr,"Error! You can't use the -x/--mixfishery option in conjunction with a baseline given with the -b option without giving at least one recognized --Pi or --fixed-Pi option.  By \"recognized\" we mean that it is correctly followed by an --end-pop option.  \n Exiting...\n");
		exit(1);

	}
	if(resamp_like_vecs_f && !mixfishery_f)  {
		fprintf(stderr,"Error! The method invoked follows the --resamp-like-vecs option, but there is no mixed fishery to simulate.  To use the --resamp-like-vecs method, you must be simulating mixed fisheries.  Exiting...\n");
		exit(1);
	}
	if(resamp_like_vecs_f && !baseline_f)  {
		fprintf(stderr,"Error! The method invoked follows the --resamp-like-vecs option which requires that the -b/--baseline-genotypes option be used to specify baseline genotypes. Exiting...\n");
		exit(1);
	}
	if(resamp_loci_loo_f && !baseline_f)  {
		fprintf(stderr,"Error! The method invoked follows the --resamp-loci-loo option which requires that the -b/--baseline-genotypes option be used to specify baseline genotypes. Exiting...\n");
		exit(1);
	}
	if(resamp_loci_loo_f + resamp_like_vecs_f + resamp_genes_loo_f + nomixcdf_f + resample_baselines_f + resample_adfg_f > 1 ) {
		fprintf(stderr,"Error! You can issue only one of:  --resamp-base, --resamp-adfg, --resamp-like-vecs, --resamp-genes-loo, --resamp-loci-loo, or --no-cdf-mix. Exiting...\n");
		exit(1);
	}
	if(mixture_f + holdout_f > 1) {
		fprintf(stderr,"Error! You can only issue one of: -t/--mixture-genotypes or -h/--holdout-genotypes. Exiting...\n");
		exit(1);
	}
	
	/* and here we transfer across the LocalPloidies to all the pop_structs, unless LocalPloidies is NULL, in which
	case we transfer over diploid for all those Ploidies */
	for(i=0;i<CurrentPop;i++) {	int j;
		P[i]->Ploidy = (int *)ECA_CALLOC(P[0]->NumLoc,sizeof(int));
		
		for(j=0;j<P[0]->NumLoc;j++) {
			if(LocalPloidies) {
				P[i]->Ploidy[j] = LocalPloidies[j];
			}
			else {
				P[i]->Ploidy[j] = 2;
			}
		}
	}
	
	
	/* and finally we figure out the method to use */
	/* start with setting defaults if they didn't get set here */
	if(!sampunit_f) Samp_Unit=0;
	if(!loostyle_f) {
		if(holdout_f) {
			LOO_style = 0;
		}
		else {
			LOO_style = 1;
		}
	}
	/* then go ahead and set the rest */
	if(LOO_style==0) {
		switch (Samp_Unit) {
			case(0):
				*Meth=NO_LOO_GENECOPY;
				break;
			case(1):
				*Meth=NO_LOO_SINGLELOCUS;
				break;
			case(2):
				*Meth=NO_LOO_MULTILOCUS;
				break;
			default:
				fprintf(stderr,"Error.  Unrecognized Samp_Unit value %d in GetGSI_Options.  Exiting...\n",Samp_Unit);
				exit(1);
				break;
		}
	}
	else if(LOO_style==1) {
		switch (Samp_Unit) {
			case(0):
				*Meth=RESAMP_LOO_GENECOPY;
				break;
			case(1):
				*Meth=RESAMP_LOO_SINGLELOCUS;
				break;
			case(2):
				*Meth=RESAMP_LOO_MULTILOCUS;
				break;
			default:
				fprintf(stderr,"Error.  Unrecognized Samp_Unit value %d in GetGSI_Options.  Exiting...\n",Samp_Unit);
				exit(1);
				break;
		}

	}
	
	return(CurrentPop);
	
}


