/* **********************************************************************************************************
 * William T. Jarratt
 * Root v6.16
 * Version: 1.00
 * Last Updated: 14 FEB 2021
 *
 * The purpose of this script is to create histograms
 * using the tree data provided by TUNL from the summer 2018 Ar run. Please note that
 * if the result file already exist you will recieve warnings for every duplicate object
 * you attempt to write. hist_2D() was seperated from the TOF generation in order to reduce
 * script run time. You only need to generate the complete histograms once (unless you
 * accidently delete them.) filepaths.h must be inculded, and the #define sections must
 * be completely filled out. No errors are returned for incorrect file paths, the script
 * just creates the file.
 *
 * *******************************************************************************************************/

#include <dirent.h>		// Needed for non root file system operations.
#include <time.h>		// Needed for time operations.
#include <unistd.h>
#include <TMath.h>
#include <stdio.h>		// Standard input and output
#include <string.h>		// String operations
#include "filepaths.h"	// Header file with all the important file path and names
#include "fitfunctions.c"	// Code that takes care of fitting functions.

int bcs11mev(void){

	int 	error 	= 0;	// Indicates errors from functions
	int 	angle	= 0;	// The angle the user chooses to analyze
	int 	i		= 0;

	char	file_name[100];	// String to hold user given file name
	char	root_path[100];	// Path for root file
	char	log_path[100];	// Path for log text file
	char	menChoi;		// Indicates which operation to perform

	char	savepSTR[]		= SAVEAREA; // From filepaths.h
	char	tree1STR[]		= GI_TREE0;	// Variables that are all caps
	char	tree2STR[]		= GI_TREE1; // are definitions from the header
	char	tree3STR[]		= GI_TREE2;	// See filepaths.h for more information
	char	tree4STR[]		= GI_TREE3;
	char	tree5STR[]		= GO_TREE0;
	char	cutSTRs[][21] 	= TCUTGS;	// Imports cut file names.

	TFile*	fRoot 	= NULL;	// File pointer for the root file
	FILE*	fLog	= NULL;	// File pointer for the log text file

	menu(&angle, &menChoi);	// Allows the user to select what they want to do.

	// The user chose to quite the program
	if((menChoi == '5') || (angle == 21)){
		printf("Quitting Application, Goodbye.\n\n");
		return 0;
	}

	// User chooses the file name. Entering a file name that already
	// exist will generate an error.
	printf("Enter a file name, do not include the file extension: ");
	scanf("%s", file_name);

	sprintf(root_path, "%s%s%s", SAVEAREA, file_name, ".root");
	sprintf(log_path, "%s%s%s", SAVEAREA, file_name, "_log.txt");
	fRoot = new TFile(root_path, "CREATE");
	fLog = fopen(log_path, "w");

	// Check to make sure that the filename is valid
	if((fRoot->IsZombie()) || (fLog == NULL)){
		printf("There was an error creating the ROOT or log file, check if it already exist\n");
		printf("or if the file path is incorrect. Quitting application.\n\n");
		return 1;
	}

	printf("\nA file has been created that contains the results:\n\t%s", root_path);
	printf("\n\nA text file has been created that logs important information:\n\t%s\n", log_path);
	printf("\nYou have chosen %d degree for your analysis.\n\n", angle);

	// Starts the log file.
	fprintf(fLog, "\nThis is a log file created for the purpose of recording the analysis of %d degree data.\n", angle);
	fprintf(fLog, "\nPath where the results are written: %s", savepSTR);
	fprintf(fLog, "\nPath for gas in tree 3672: %s", tree1STR);
	fprintf(fLog, "\nPath for gas in tree 3673: %s", tree2STR);
	fprintf(fLog, "\nPath for gas in tree 3674: %s", tree3STR);
	fprintf(fLog, "\nPath for gas in tree 3675: %s", tree4STR);
	fprintf(fLog, "\nPath for gas out tree 3702: %s", tree5STR);
	fprintf(fLog, "\n\n");

	fprintf(fLog, "The cut files used for this script cycle are: \n");

	fprintf(fLog, "%s", cutSTRs[0]);
	for (i=1; i<11; i++){
		fprintf(fLog, ", %s", cutSTRs[i]);
	}
	fprintf(fLog, "\n%s", cutSTRs[11]);
	for (i=12; i<21; i++){
		fprintf(fLog, ", %s", cutSTRs[i]);
	}
	fprintf(fLog, "\n\n");

	// A switch statement is used along with menChoi to determine which function to call.
	switch(menChoi){

		case '1':
		printf("\nWe are in the process of implementing this option (1).\n");
		break;

		case '2':
		printf("\nWe are in the process of implementing this option (2).\n");
		break;

		case '3':
		printf("\nBeginning TOF generation. The plot names are set in filepaths.h.\n");
		error = Generate_TOF(angle, fRoot, fLog);
		break;

		case '4':
		printf("\nBeginning curve fit. The TOF names and source file name are specified in filepaths.h.\n\n");
		error = Apply_Fits(angle, fRoot, fLog);
		break;

		default:
		return 1;
	}


	fRoot->Close();
	fclose(fLog);

	printf("\nThe script is finished, goodbye.\n");

	if(error == 1){ return 1; }
	return 0;

}


/* ******************************************************************************************************
 * This function takes the tree data and applies the cuts for the gas in/out runs and generates
 * TOFs. Replaces iandiixCs_Test1.C and iandiixCs_Test12.C
 *
 * Inputs:
 *
 * Returns:
 *
 * ******************************************************************************************************/
int Generate_TOF(int angle, TFile* f1, FILE* f3){

	/* Here we declare all the important pointers and variables. */
	char	leafPH[][10] = { "0TPH", "0MPH", "0BPH", "3TPH", "3MPH", "3BPH", "6TPH", "6MPH", "6BPH",\
							"9TPH", "9MPH", "9BPH", "12TPH", "12MPH", "12BPH", "15TPH", "15MPH", "15BPH",\
							"18TPH", "18MPH", "18BPH" };
	char	leafPSD[][10] = { "0TPSD", "0MPSD", "0BPSD", "3TPSD", "3MPSD", "3BPSD", "6TPSD", "6MPSD", "6BPSD",\
							"9TPSD", "9MPSD", "9BPSD", "12TPSD", "12MPSD", "12BPSD", "15TPSD", "15MPSD",\
							"15BPSD", "18TPSD", "18MPSD", "18BPSD"};
	char	leafTDC[][10] = { "0TTDC", "0MTDC", "0BTDC", "3TTDC", "3MTDC", "3BTDC", "6TTDC", "6MTDC", "6BTDC",\
							"9TTDC", "9MTDC", "9BTDC", "12TTDC", "12MTDC", "12BTDC", "15TTDC", "15MTDC",\
							"15BTDC", "18TTDC", "18MTDC", "18BTDC"};
	char	cutSTR[][21] = 	TCUTGS;

	TCutG*	tCutG[21];
	TH1F* 	cut1cs[21];		// Gas in TOF gated by 1xCs with cut applied.
	TH1F* 	cut1cs2[21];	// Gas out TOF gated by 1xCs with cut applied.
	TH1F*	cut2cs[21];		// Gas in TOF gated by 2xCs with cut applied.
	TH1F*	cut2cs2[21];	// Gas out TOF gated by 2xCs with cut applied.
	TH1F*	uncut1cs[21];	// Gas in TOF gated by 1xCs with no cut applied.
	TH1F*	uncut1cs2[21];	// Gas out TOF gated by 1xCs with no cut applied.
	TH1F*	uncut2cs[21];	// Gas in TOF gated by 2xCs with no cut applied.
	TH1F*	uncut2cs2[21];	// Gas out TOF gated by 2xCs with no cut applied.
	TH1F* 	gasin_nocut[21];	// Gas in TOF no Cs no cut for top.
	TH1F* 	gasout_nocut[21];	// Gas out TOF no Cs no cut for top.


	TCanvas *c_gasin_comb = new TCanvas("gasin_combined");
	TCanvas *c_gasin_comb_nocut = new TCanvas("gasin_combined_nocut");
  TCanvas *c_gasout_comb = new TCanvas("gasout_combined");

	int		phLimB1[21] = PHLIMB1C;
	int		phLimB2[21] = PHLIMB2C;
	int		psdLimB1[21] = PSDLIMB1C;
	int		psdLimB2[21] = PSDLIMB2C;
	int		rfOffIn[21] = RFINOFF;
	int		rfOffOut[21] = RFOUTOFF;
	int		i = 0;
	int   rebin_h = 4;
	float 	ph1[21];
	float	psd1[21];
	float	ph2[21];
	float	psd2[21];
	float	tdc1[21];
	float	tdc2[21];
	float	rfpick1;
	float	rfpick2;

	TChain 	cch1("InTree");
	TChain	cch2("OutTree");

	gROOT->SetBatch(kTRUE);

	TFile*	f2 = NULL;
	f2 = new TFile(CUTPATH, "READ");

	if(f2->IsZombie())
	{
		printf("\nFile path for cuts does not exist. Update filepaths.h\n");
		return 1;
	}

	/* Here we merge all the trees together because we want to read the data from all four runs. We
	 * use SetBranchAddress() to link each leaf to a variable. The variable will recieve the
	 * data. */
	cch1.Add(GI_TREE0);
	cch1.Add(GI_TREE1);
	cch1.Add(GI_TREE2);
	cch1.Add(GI_TREE3);
	cch2.Add(GO_TREE0);

	for(i = angle; i <= (angle+2); i++){
		cch1.SetBranchAddress(leafPH[i], &ph1[i]);
		cch1.SetBranchAddress(leafPSD[i], &psd1[i]);
		cch1.SetBranchAddress(leafTDC[i], &tdc1[i]);
		cch2.SetBranchAddress(leafPH[i], &ph2[i]);
		cch2.SetBranchAddress(leafPSD[i], &psd2[i]);
		cch2.SetBranchAddress(leafTDC[i], &tdc2[i]);
	}

	cch1.SetBranchAddress("RFTDC", &rfpick1);
	cch2.SetBranchAddress("RFTDC", &rfpick2);

	for(i = angle; i <= (angle+2); i++){
		f2->GetObject(cutSTR[i], tCutG[i]);
	}

	int bin_number = (2048*2);

	// Gas in TOF gated by 1xCs with cut applied. 1
	cut1cs[0] = new TH1F("TOF_1Cs_0T_CUT_IN", "0 Degree Top Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[1] = new TH1F("TOF_1Cs_0M_CUT_IN", "0 Degree Mid Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[2] = new TH1F("TOF_1Cs_0B_CUT_IN", "0 Degree Bot Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[3] = new TH1F("TOF_1Cs_3T_CUT_IN", "3 Degree Top Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[4] = new TH1F("TOF_1Cs_3M_CUT_IN", "3 Degree Mid Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[5] = new TH1F("TOF_1Cs_3B_CUT_IN", "3 Degree Bot Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[6] = new TH1F("TOF_1Cs_6T_CUT_IN", "6 Degree Top Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[7] = new TH1F("TOF_1Cs_6M_CUT_IN", "6 Degree Mid Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[8] = new TH1F("TOF_1Cs_6B_CUT_IN", "6 Degree Bot Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[9] = new TH1F("TOF_1Cs_9T_CUT_IN", "9 Degree Top Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[10] = new TH1F("TOF_1Cs_9M_CUT_IN", "9 Degree Mid Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[11] = new TH1F("TOF_1Cs_9B_CUT_IN", "9 Degree Bot Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[12] = new TH1F("TOF_1Cs_12T_CUT_IN", "12 Degree Top Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[13] = new TH1F("TOF_1Cs_12M_CUT_IN", "12 Degree Mid Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[14] = new TH1F("TOF_1Cs_12B_CUT_IN", "12 Degree Bot Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[15] = new TH1F("TOF_1Cs_15T_CUT_IN", "15 Degree Top Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[16] = new TH1F("TOF_1Cs_15M_CUT_IN", "15 Degree Mid Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[17] = new TH1F("TOF_1Cs_15B_CUT_IN", "15 Degree Bot Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[18] = new TH1F("TOF_1Cs_18T_CUT_IN", "18 Degree Top Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[19] = new TH1F("TOF_1Cs_18M_CUT_IN", "18 Degree Mid Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs[20] = new TH1F("TOF_1Cs_18B_CUT_IN", "18 Degree Bot Gas In TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);

	// Gas out TOF gated by 1xCs with cut applied. 2
	cut1cs2[0] = new TH1F("TOF_1Cs_0T_CUT_OUT", "0 Degree Top Gas Out TOF Gated By 1xCs w/cut", (bin_number), 0., 4095.);
	cut1cs2[1] = new TH1F("TOF_1Cs_0M_CUT_OUT", "0 Degree Mid Gas Out TOF Gated By 1xCs w/cut", (bin_number), 0., 4095.);
	cut1cs2[2] = new TH1F("TOF_1Cs_0B_CUT_OUT", "0 Degree Bot Gas Out TOF Gated By 1xCs w/cut", (bin_number), 0., 4095.);
	cut1cs2[3] = new TH1F("TOF_1Cs_3T_CUT_OUT", "3 Degree Top Gas Out TOF Gated By 1xCs w/cut", (bin_number), 0., 4095.);
	cut1cs2[4] = new TH1F("TOF_1Cs_3M_CUT_OUT", "3 Degree Mid Gas Out TOF Gated By 1xCs w/cut", (bin_number), 0., 4095.);
	cut1cs2[5] = new TH1F("TOF_1Cs_3B_CUT_OUT", "3 Degree Bot Gas Out TOF Gated By 1xCs w/cut", (bin_number), 0., 4095.);
	cut1cs2[6] = new TH1F("TOF_1Cs_6T_CUT_OUT", "6 Degree Top Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[7] = new TH1F("TOF_1Cs_6M_CUT_OUT", "6 Degree Mid Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[8] = new TH1F("TOF_1Cs_6B_CUT_OUT", "6 Degree Bot Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[9] = new TH1F("TOF_1Cs_9T_CUT_OUT", "9 Degree Top Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[10] = new TH1F("TOF_1Cs_9M_CUT_OUT", "9 Degree Mid Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[11] = new TH1F("TOF_1Cs_9B_CUT_OUT", "9 Degree Bot Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[12] = new TH1F("TOF_1Cs_12T_CUT_OUT", "12 Degree Top Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[13] = new TH1F("TOF_1Cs_12M_CUT_OUT", "12 Degree Mid Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[14] = new TH1F("TOF_1Cs_12B_CUT_OUT", "12 Degree Bot Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[15] = new TH1F("TOF_1Cs_15T_CUT_OUT", "15 Degree Top Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[16] = new TH1F("TOF_1Cs_15M_CUT_OUT", "15 Degree Mid Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[17] = new TH1F("TOF_1Cs_15B_CUT_OUT", "15 Degree Bot Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[18] = new TH1F("TOF_1Cs_18T_CUT_OUT", "18 Degree Top Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[19] = new TH1F("TOF_1Cs_18M_CUT_OUT", "18 Degree Mid Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);
	cut1cs2[20] = new TH1F("TOF_1Cs_18B_CUT_OUT", "18 Degree Bot Gas Out TOF Gated By 1xCs w/cut", bin_number, 0., 4095.);

	// Gas in TOF gated by 2xCs with cut applied. 3
	cut2cs[0] = new TH1F("TOF_2Cs_0T_CUT_IN", "0 Degree Top Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[1] = new TH1F("TOF_2Cs_0M_CUT_IN", "0 Degree Mid Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[2] = new TH1F("TOF_2Cs_0B_CUT_IN", "0 Degree Bot Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[3] = new TH1F("TOF_2Cs_3T_CUT_IN", "3 Degree Top Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[4] = new TH1F("TOF_2Cs_3M_CUT_IN", "3 Degree Mid Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[5] = new TH1F("TOF_2Cs_3B_CUT_IN", "3 Degree Bot Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[6] = new TH1F("TOF_2Cs_6T_CUT_IN", "6 Degree Top Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[7] = new TH1F("TOF_2Cs_6M_CUT_IN", "6 Degree Mid Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[8] = new TH1F("TOF_2Cs_6B_CUT_IN", "6 Degree Bot Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[9] = new TH1F("TOF_2Cs_9T_CUT_IN", "9 Degree Top Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[10] = new TH1F("TOF_2Cs_9M_CUT_IN", "9 Degree Mid Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[11] = new TH1F("TOF_2Cs_9B_CUT_IN", "9 Degree Bot Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[12] = new TH1F("TOF_2Cs_12T_CUT_IN", "12 Degree Top Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[13] = new TH1F("TOF_2Cs_12M_CUT_IN", "12 Degree Mid Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[14] = new TH1F("TOF_2Cs_12B_CUT_IN", "12 Degree Bot Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[15] = new TH1F("TOF_2Cs_15T_CUT_IN", "15 Degree Top Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[16] = new TH1F("TOF_2Cs_15M_CUT_IN", "15 Degree Mid Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[17] = new TH1F("TOF_2Cs_15B_CUT_IN", "15 Degree Bot Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[18] = new TH1F("TOF_2Cs_18T_CUT_IN", "18 Degree Top Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[19] = new TH1F("TOF_2Cs_18M_CUT_IN", "18 Degree Mid Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs[20] = new TH1F("TOF_2Cs_18B_CUT_IN", "18 Degree Bot Gas In TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);

	// Gas out TOF gated by 2xCs with cut applied. 4
	cut2cs2[0] = new TH1F("TOF_2Cs_0T_CUT_OUT", "0 Degree Top Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[1] = new TH1F("TOF_2Cs_0M_CUT_OUT", "0 Degree Mid Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[2] = new TH1F("TOF_2Cs_0B_CUT_OUT", "0 Degree Bot Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[3] = new TH1F("TOF_2Cs_3T_CUT_OUT", "3 Degree Top Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[4] = new TH1F("TOF_2Cs_3M_CUT_OUT", "3 Degree Mid Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[5] = new TH1F("TOF_2Cs_3B_CUT_OUT", "3 Degree Bot Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[6] = new TH1F("TOF_2Cs_6T_CUT_OUT", "6 Degree Top Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[7] = new TH1F("TOF_2Cs_6M_CUT_OUT", "6 Degree Mid Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[8] = new TH1F("TOF_2Cs_6B_CUT_OUT", "6 Degree Bot Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[9] = new TH1F("TOF_2Cs_9T_CUT_OUT", "9 Degree Top Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[10] = new TH1F("TOF_2Cs_9M_CUT_OUT", "9 Degree Mid Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[11] = new TH1F("TOF_2Cs_9B_CUT_OUT", "9 Degree Bot Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[12] = new TH1F("TOF_2Cs_12T_CUT_OUT", "12 Degree Top Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[13] = new TH1F("TOF_2Cs_12M_CUT_OUT", "12 Degree Mid Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[14] = new TH1F("TOF_2Cs_12B_CUT_OUT", "12 Degree Bot Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[15] = new TH1F("TOF_2Cs_15T_CUT_OUT", "15 Degree Top Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[16] = new TH1F("TOF_2Cs_15M_CUT_OUT", "15 Degree Mid Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[17] = new TH1F("TOF_2Cs_15B_CUT_OUT", "15 Degree Bot Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[18] = new TH1F("TOF_2Cs_18T_CUT_OUT", "18 Degree Top Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[19] = new TH1F("TOF_2Cs_18M_CUT_OUT", "18 Degree Mid Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);
	cut2cs2[20] = new TH1F("TOF_2Cs_18B_CUT_OUT", "18 Degree Bot Gas Out TOF Gated By 2xCs w/cut", bin_number, 0., 4095.);

	// Gas in TOF gated by 1xCs with no cut applied. 5
	uncut1cs[0] = new TH1F("TOF_1Cs_0T_NOCUT_IN", "0 Degree Top Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[1] = new TH1F("TOF_1Cs_0M_NOCUT_IN", "0 Degree Mid Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[2] = new TH1F("TOF_1Cs_0B_NOCUT_IN", "0 Degree Bot Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[3] = new TH1F("TOF_1Cs_3T_NOCUT_IN", "3 Degree Top Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[4] = new TH1F("TOF_1Cs_3M_NOCUT_IN", "3 Degree Mid Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[5] = new TH1F("TOF_1Cs_3B_NOCUT_IN", "3 Degree Bot Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[6] = new TH1F("TOF_1Cs_6T_NOCUT_IN", "6 Degree Top Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[7] = new TH1F("TOF_1Cs_6M_NOCUT_IN", "6 Degree Mid Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[8] = new TH1F("TOF_1Cs_6B_NOCUT_IN", "6 Degree Bot Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[9] = new TH1F("TOF_1Cs_9T_NOCUT_IN", "9 Degree Top Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[10] = new TH1F("TOF_1Cs_9M_NOCUT_IN", "9 Degree Mid Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[11] = new TH1F("TOF_1Cs_9B_NOCUT_IN", "9 Degree Bot Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[12] = new TH1F("TOF_1Cs_12T_NOCUT_IN", "12 Degree Top Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[13] = new TH1F("TOF_1Cs_12M_NOCUT_IN", "12 Degree Mid Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[14] = new TH1F("TOF_1Cs_12B_NOCUT_IN", "12 Degree Bot Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[15] = new TH1F("TOF_1Cs_15T_NOCUT_IN", "15 Degree Top Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[16] = new TH1F("TOF_1Cs_15M_NOCUT_IN", "15 Degree Mid Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[17] = new TH1F("TOF_1Cs_15B_NOCUT_IN", "15 Degree Bot Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[18] = new TH1F("TOF_1Cs_18T_NOCUT_IN", "18 Degree Top Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[19] = new TH1F("TOF_1Cs_18M_NOCUT_IN", "18 Degree Mid Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs[20] = new TH1F("TOF_1Cs_18B_NOCUT_IN", "18 Degree Bot Gas In TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);

	// Gas out TOF gated by 1xCs with no cut applied. 6
	uncut1cs2[0] = new TH1F("TOF_1Cs_0T_NOCUT_OUT", "0 Degree Top Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[1] = new TH1F("TOF_1Cs_0M_NOCUT_OUT", "0 Degree Mid Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[2] = new TH1F("TOF_1Cs_0B_NOCUT_OUT", "0 Degree Bot Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[3] = new TH1F("TOF_1Cs_3T_NOCUT_OUT", "3 Degree Top Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[4] = new TH1F("TOF_1Cs_3M_NOCUT_OUT", "3 Degree Mid Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[5] = new TH1F("TOF_1Cs_3B_NOCUT_OUT", "3 Degree Bot Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[6] = new TH1F("TOF_1Cs_6T_NOCUT_OUT", "6 Degree Top Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[7] = new TH1F("TOF_1Cs_6M_NOCUT_OUT", "6 Degree Mid Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[8] = new TH1F("TOF_1Cs_6B_NOCUT_OUT", "6 Degree Bot Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[9] = new TH1F("TOF_1Cs_9T_NOCUT_OUT", "9 Degree Top Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[10] = new TH1F("TOF_1Cs_9M_NOCUT_OUT", "9 Degree Mid Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[11] = new TH1F("TOF_1Cs_9B_NOCUT_OUT", "9 Degree Bot Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[12] = new TH1F("TOF_1Cs_12T_NOCUT_OUT", "12 Degree Top Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[13] = new TH1F("TOF_1Cs_12M_NOCUT_OUT", "12 Degree Mid Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[14] = new TH1F("TOF_1Cs_12B_NOCUT_OUT", "12 Degree Bot Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[15] = new TH1F("TOF_1Cs_15T_NOCUT_OUT", "15 Degree Top Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[16] = new TH1F("TOF_1Cs_15M_NOCUT_OUT", "15 Degree Mid Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[17] = new TH1F("TOF_1Cs_15B_NOCUT_OUT", "15 Degree Bot Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[18] = new TH1F("TOF_1Cs_18T_NOCUT_OUT", "18 Degree Top Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[19] = new TH1F("TOF_1Cs_18M_NOCUT_OUT", "18 Degree Mid Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);
	uncut1cs2[20] = new TH1F("TOF_1Cs_18B_NOCUT_OUT", "18 Degree Bot Gas Out TOF Gated By 1xCs w/o cut", bin_number, 0., 4095.);

	// Gas in TOF gated by 2xCs with no cut applied. 7
	uncut2cs[0] = new TH1F("TOF_2Cs_0T_NOCUT_IN", "0 Degree Top Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[1] = new TH1F("TOF_2Cs_0M_NOCUT_IN", "0 Degree Mid Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[2] = new TH1F("TOF_2Cs_0B_NOCUT_IN", "0 Degree Bot Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[3] = new TH1F("TOF_2Cs_3T_NOCUT_IN", "3 Degree Top Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[4] = new TH1F("TOF_2Cs_3M_NOCUT_IN", "3 Degree Mid Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[5] = new TH1F("TOF_2Cs_3B_NOCUT_IN", "3 Degree Bot Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[6] = new TH1F("TOF_2Cs_6T_NOCUT_IN", "6 Degree Top Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[7] = new TH1F("TOF_2Cs_6M_NOCUT_IN", "6 Degree Mid Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[8] = new TH1F("TOF_2Cs_6B_NOCUT_IN", "6 Degree Bot Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[9] = new TH1F("TOF_2Cs_9T_NOCUT_IN", "9 Degree Top Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[10] = new TH1F("TOF_2Cs_9M_NOCUT_IN", "9 Degree Mid Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[11] = new TH1F("TOF_2Cs_9B_NOCUT_IN", "9 Degree Bot Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[12] = new TH1F("TOF_2Cs_12T_NOCUT_IN", "12 Degree Top Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[13] = new TH1F("TOF_2Cs_12M_NOCUT_IN", "12 Degree Mid Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[14] = new TH1F("TOF_2Cs_12B_NOCUT_IN", "12 Degree Bot Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[15] = new TH1F("TOF_2Cs_15T_NOCUT_IN", "15 Degree Top Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[16] = new TH1F("TOF_2Cs_15M_NOCUT_IN", "15 Degree Mid Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[17] = new TH1F("TOF_2Cs_15B_NOCUT_IN", "15 Degree Bot Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[18] = new TH1F("TOF_2Cs_18T_NOCUT_IN", "18 Degree Top Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[19] = new TH1F("TOF_2Cs_18M_NOCUT_IN", "18 Degree Mid Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs[20] = new TH1F("TOF_2Cs_18B_NOCUT_IN", "18 Degree Bot Gas In TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);

	// Gas out TOF gated by 2xCs with no cut applied. 8
	uncut2cs2[0] = new TH1F("TOF_2Cs_0T_NOCUT_OUT", "0 Degree Top Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[1] = new TH1F("TOF_2Cs_0M_NOCUT_OUT", "0 Degree Mid Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[2] = new TH1F("TOF_2Cs_0B_NOCUT_OUT", "0 Degree Bot Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[3] = new TH1F("TOF_2Cs_3T_NOCUT_OUT", "3 Degree Top Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[4] = new TH1F("TOF_2Cs_3M_NOCUT_OUT", "3 Degree Mid Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[5] = new TH1F("TOF_2Cs_3B_NOCUT_OUT", "3 Degree Bot Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[6] = new TH1F("TOF_2Cs_6T_NOCUT_OUT", "6 Degree Top Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[7] = new TH1F("TOF_2Cs_6M_NOCUT_OUT", "6 Degree Mid Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[8] = new TH1F("TOF_2Cs_6B_NOCUT_OUT", "6 Degree Bot Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[9] = new TH1F("TOF_2Cs_9T_NOCUT_OUT", "9 Degree Top Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[10] = new TH1F("TOF_2Cs_9M_NOCUT_OUT", "9 Degree Mid Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[11] = new TH1F("TOF_2Cs_9B_NOCUT_OUT", "9 Degree Bot Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[12] = new TH1F("TOF_2Cs_12T_NOCUT_OUT", "12 Degree Top Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[13] = new TH1F("TOF_2Cs_12M_NOCUT_OUT", "12 Degree Mid Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[14] = new TH1F("TOF_2Cs_12B_NOCUT_OUT", "12 Degree Bot Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[15] = new TH1F("TOF_2Cs_15T_NOCUT_OUT", "15 Degree Top Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[16] = new TH1F("TOF_2Cs_15M_NOCUT_OUT", "15 Degree Mid Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[17] = new TH1F("TOF_2Cs_15B_NOCUT_OUT", "15 Degree Bot Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[18] = new TH1F("TOF_2Cs_18T_NOCUT_OUT", "18 Degree Top Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[19] = new TH1F("TOF_2Cs_18M_NOCUT_OUT", "18 Degree Mid Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);
	uncut2cs2[20] = new TH1F("TOF_2Cs_18B_NOCUT_OUT", "18 Degree Bot Gas Out TOF Gated By 2xCs w/o cut", bin_number, 0., 4095.);

	// Gas in TOF no Cs no cut.
	gasin_nocut[0] = new TH1F("TOF_0T_IN_RAW", "0 Degree Top Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[1] = new TH1F("TOF_0M_IN_RAW", "0 Degree Mid Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[2] = new TH1F("TOF_0B_IN_RAW", "0 Degree Bot Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[3] = new TH1F("TOF_3T_IN_RAW", "3 Degree Top Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[4] = new TH1F("TOF_3M_IN_RAW", "3 Degree Mid Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[5] = new TH1F("TOF_3B_IN_RAW", "3 Degree Bot Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[6] = new TH1F("TOF_6T_IN_RAW", "6 Degree Top Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[7] = new TH1F("TOF_6M_IN_RAW", "6 Degree Mid Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[8] = new TH1F("TOF_6B_IN_RAW", "6 Degree Bot Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[9] = new TH1F("TOF_9T_IN_RAW", "9 Degree Top Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[10] = new TH1F("TOF_9M_IN_RAW", "9 Degree Mid Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[11] = new TH1F("TOF_9B_IN_RAW", "9 Degree Bot Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[12] = new TH1F("TOF_12T_IN_RAW", "12 Degree Top Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[13] = new TH1F("TOF_12M_IN_RAW", "12 Degree Mid Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[14] = new TH1F("TOF_12B_IN_RAW", "12 Degree Bot Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[15] = new TH1F("TOF_15T_IN_RAW", "15 Degree Top Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[16] = new TH1F("TOF_15M_IN_RAW", "15 Degree Mid Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[17] = new TH1F("TOF_15B_IN_RAW", "15 Degree Bot Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[18] = new TH1F("TOF_18T_IN_RAW", "18 Degree Top Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[19] = new TH1F("TOF_18M_IN_RAW", "18 Degree Mid Gas In TOF with no cuts", bin_number, 0., 4095.);
	gasin_nocut[20] = new TH1F("TOF_18B_IN_RAW", "18 Degree Bot Gas In TOF with no cuts", bin_number, 0., 4095.);

	// Gas out TOF no Cs no cut.
	gasout_nocut[0] = new TH1F("TOF_0T_OUT_RAW", "0 Degree Top Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[1] = new TH1F("TOF_0M_OUT_RAW", "0 Degree Mid Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[2] = new TH1F("TOF_0B_OUT_RAW", "0 Degree Bot Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[3] = new TH1F("TOF_3T_OUT_RAW", "3 Degree Top Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[4] = new TH1F("TOF_3M_OUT_RAW", "3 Degree Mid Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[5] = new TH1F("TOF_3B_OUT_RAW", "3 Degree Bot Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[6] = new TH1F("TOF_6T_OUT_RAW", "6 Degree Top Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[7] = new TH1F("TOF_6M_OUT_RAW", "6 Degree Mid Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[8] = new TH1F("TOF_6B_OUT_RAW", "6 Degree Bot Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[9] = new TH1F("TOF_9T_OUT_RAW", "9 Degree Top Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[10] = new TH1F("TOF_9M_OUT_RAW", "9 Degree Mid Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[11] = new TH1F("TOF_9B_OUT_RAW", "9 Degree Bot Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[12] = new TH1F("TOF_12T_OUT_RAW", "12 Degree Top Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[13] = new TH1F("TOF_12M_OUT_RAW", "12 Degree Mid Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[14] = new TH1F("TOF_12B_OUT_RAW", "12 Degree Bot Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[15] = new TH1F("TOF_15T_OUT_RAW", "15 Degree Top Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[16] = new TH1F("TOF_15M_OUT_RAW", "15 Degree Mid Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[17] = new TH1F("TOF_15B_OUT_RAW", "15 Degree Bot Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[18] = new TH1F("TOF_18T_OUT_RAW", "18 Degree Top Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[19] = new TH1F("TOF_18M_OUT_RAW", "18 Degree Mid Gas Out TOF with no cuts", bin_number, 0., 4095.);
	gasout_nocut[20] = new TH1F("TOF_18B_OUT_RAW", "18 Degree Bot Gas Out TOF with no cuts", bin_number, 0., 4095.);


	printf("\nNow beginning TOF generation...\n\n");

	printf("%lld gas in entries will be read.\n", cch1.GetEntries());

	/* First loop that creates the gas in TOF. Note that each angle TOF has
	 * an offset. This is required to line up TOF. Unsure of reasoning as of
	 * right now. It is either to line up degrees or gas in and gas out. */
	for(i = 0; i <= cch1.GetEntries(); i++){
		cch1.GetEvent(i);
		printf("\r%d", i);


		// top
		// Gas in, 1xCs, w cut
		if (tCutG[angle]->IsInside(psd1[angle], ph1[angle]) && (tdc1[angle] > 400) && (ph1[angle] > phLimB1[angle]) && (ph1[angle] < 3800)){
			cut1cs[angle]->Fill(rfpick1 + rfOffIn[angle]);
		}
		// Gas in, 2xCs, w cut
		if (tCutG[angle]->IsInside(psd1[angle], ph1[angle]) && (tdc1[angle] > 400) && (ph1[angle] > phLimB2[angle]) && (ph1[angle] < 3800)){
			cut2cs[angle]->Fill(rfpick1 + rfOffIn[angle]);
		}
		// Gas in, 1xCs, wo cut
		//if ((tdc1[angle] > 400) && (ph1[angle] > phLimB1[angle]) && (ph1[angle] < 3800) && (psd1[angle] > psdLimB1[angle])){
			//uncut1cs[angle]->Fill(rfpick1 + rfOffIn[angle]);
		//}



		//no cut
		if ((tdc1[angle] > 400) && (ph1[angle] > phLimB1[angle]) && (ph1[angle] < 3800)){
			uncut1cs[angle]->Fill(rfpick1 + rfOffIn[angle]);
		}



		// Gas in, 2xCs, wo cut
		if ((tdc1[angle] > 400) && (ph1[angle] > phLimB2[angle]) && (ph1[angle] < 3800) && (psd1[angle] > psdLimB2[angle])){
			uncut2cs[angle]->Fill(rfpick1 + rfOffIn[angle]);
		}
		gasin_nocut[angle]->Fill(rfpick1 + rfOffIn[angle]);


		//middle
		// Gas in, 1xCs, w cut
		if (tCutG[angle+1]->IsInside(psd1[angle+1], ph1[angle+1]) && (tdc1[angle+1] > 400) && (ph1[angle+1] > phLimB1[angle+1]) && (ph1[angle+1] < 3800)){
			cut1cs[angle+1]->Fill(rfpick1 + rfOffIn[angle+1]);
		}
		// Gas in, 2xCs, w cut
		if (tCutG[angle+1]->IsInside(psd1[angle+1], ph1[angle+1]) && (tdc1[angle+1] > 400) && (ph1[angle+1] > phLimB2[angle+1]) && (ph1[angle+1] < 3800)){
			cut2cs[angle+1]->Fill(rfpick1 + rfOffIn[angle+1]);
		}
		// Gas in, 1xCs, wo cut
		//if ((tdc1[angle+1] > 400) && (ph1[angle+1] > phLimB1[angle+1]) && (ph1[angle+1] < 3800) && (psd1[angle+1] > psdLimB1[angle+1])){
			//uncut1cs[angle+1]->Fill(rfpick1 + rfOffIn[angle+1]);
		//}

			//no cut
			if ((tdc1[angle+1] > 400) && (ph1[angle+1] > phLimB1[angle+1]) && (ph1[angle+1] < 3800)){
				uncut1cs[angle+1]->Fill(rfpick1 + rfOffIn[angle+1]);
			}



		// Gas in, 2xCs, wo cut
		if ((tdc1[angle+1] > 400) && (ph1[angle+1] > phLimB2[angle+1]) && (ph1[angle+1] < 3800) && (psd1[angle+1] > psdLimB2[angle+1])){
			uncut2cs[angle+1]->Fill(rfpick1 + rfOffIn[angle+1]);
		}
		gasin_nocut[angle+1]->Fill(rfpick1 + rfOffIn[angle+1]);


		//bottom
		// Gas in, 1xCs, w cut
		if (tCutG[angle+2]->IsInside(psd1[angle+2], ph1[angle+2]) && (tdc1[angle+2] > 400) && (ph1[angle+2] > phLimB1[angle+2]) && (ph1[angle+2] < 3800)){
			cut1cs[angle+2]->Fill(rfpick1 + rfOffIn[angle+2]);
		}
		// Gas in, 2xCs, w cut
		if (tCutG[angle+2]->IsInside(psd1[angle+2], ph1[angle+2]) && (tdc1[angle+2] > 400) && (ph1[angle+2] > phLimB2[angle+2]) && (ph1[angle+2] < 3800)){
			cut2cs[angle+2]->Fill(rfpick1 + rfOffIn[angle+2]);
		}
		// Gas in, 1xCs, wo cut
		//if ((tdc1[angle+2] > 400) && (ph1[angle+2] > phLimB1[angle+2]) && (ph1[angle+2] < 3800) && (psd1[angle+2] > psdLimB1[angle+2])){
		//	uncut1cs[angle+2]->Fill(rfpick1 + rfOffIn[angle+2]);
		//}


		//no cut
		if ((tdc1[angle+2] > 400) && (ph1[angle+2] > phLimB1[angle+2]) && (ph1[angle+2] < 3800)){
			uncut1cs[angle+2]->Fill(rfpick1 + rfOffIn[angle+2]);
		}




		// Gas in, 2xCs, wo cut
		if ((tdc1[angle+2] > 400) && (ph1[angle+2] > phLimB2[angle+2]) && (ph1[angle+2] < 3800) && (psd1[angle+2] > psdLimB2[angle+2])){
			uncut2cs[angle+2]->Fill(rfpick1 + rfOffIn[angle+2]);
		}
		gasin_nocut[angle+2]->Fill(rfpick1 + rfOffIn[angle+2]);
	}


	printf("\n\n%lld gas out entries will be read.\n", cch2.GetEntries());
	/* Second loop that creates the gas out TOF. */

	for(i = 0; i <= cch2.GetEntries(); i++)
	{
		cch2.GetEvent(i);
		printf("\r%d", i);

		//top
		// Gas out, 1xCs, w cut
		if (tCutG[angle]->IsInside(psd2[angle], ph2[angle]) && (tdc2[angle] > 400) && (ph2[angle] > phLimB1[angle]) && (ph2[angle] < 3800)){
			cut1cs2[angle]->Fill(rfpick2 + rfOffIn[angle]);
		}
		// Gas out, 2xCs, w cut
		if (tCutG[angle]->IsInside(psd2[angle], ph2[angle]) && (tdc2[angle] > 400) && (ph2[angle] > phLimB2[angle]) && (ph2[angle] < 3800)){
			cut2cs2[angle]->Fill(rfpick2 + rfOffIn[angle]);
		}
		// Gas out, 1xCs, wo cut
		if ((tdc2[angle] > 400) && (ph2[angle] > phLimB1[angle]) && (ph2[angle] < 3800) && (psd2[angle] > psdLimB1[angle])){
			uncut1cs2[angle]->Fill(rfpick2 + rfOffIn[angle]);
		}
		// Gas out, 2xCs, wo cut
		if ((tdc2[angle] > 400) && (ph2[angle] > phLimB2[angle]) && (ph2[angle] < 3800) && (psd2[angle] > psdLimB2[angle])){
			uncut2cs2[angle]->Fill(rfpick2 + rfOffIn[angle]);
		}
		gasout_nocut[angle]->Fill(rfpick2 + rfOffIn[angle]);


		//middle
		// Gas out, 1xCs, w cut
		if (tCutG[angle+1]->IsInside(psd2[angle+1], ph2[angle+1]) && (tdc2[angle+1] > 400) && (ph2[angle+1] > phLimB1[angle+1]) && (ph2[angle+1] < 3800)){
			cut1cs2[angle+1]->Fill(rfpick2 + rfOffIn[angle+1]);
		}
		// Gas out, 2xCs, w cut
		if (tCutG[angle+1]->IsInside(psd2[angle+1], ph2[angle+1]) && (tdc2[angle+1] > 400) && (ph2[angle+1] > phLimB2[angle+1]) && (ph2[angle+1] < 3800)){
			cut2cs2[angle+1]->Fill(rfpick2 + rfOffIn[angle+1]);
		}
		// Gas out, 1xCs, wo cut
		if ((tdc2[angle+1] > 400) && (ph2[angle+1] > phLimB1[angle+1]) && (ph2[angle+2] < 3800) && (psd2[angle+1] > psdLimB1[angle+1])){
			uncut1cs2[angle+1]->Fill(rfpick2 + rfOffIn[angle+1]);
		}
		// Gas out, 2xCs, wo cut
		if ((tdc2[angle+1] > 400) && (ph2[angle+1] > phLimB2[angle+1]) && (ph2[angle+1] < 3800) && (psd2[angle+1] > psdLimB2[angle+1])){
			uncut2cs2[angle+1]->Fill(rfpick2 + rfOffIn[angle+1]);
		}
		gasout_nocut[angle+1]->Fill(rfpick2 + rfOffIn[angle+1]);


		//bottom
		// Gas out, 1xCs, w cut
		if (tCutG[angle+2]->IsInside(psd2[angle+2], ph2[angle+2]) && (tdc2[angle+2] > 400) && (ph2[angle+2] > phLimB1[angle+2]) && (ph2[angle+2] < 3800)){
			cut1cs2[angle+2]->Fill(rfpick2 + rfOffIn[angle+2]);
		}
		// Gas out, 2xCs, w cut
		if (tCutG[angle+2]->IsInside(psd2[angle+2], ph2[angle+2]) && (tdc2[angle+2] > 400) && (ph2[angle+2] > phLimB2[angle+2]) && (ph2[angle+2] < 3800)){
			cut2cs2[angle+2]->Fill(rfpick2 + rfOffIn[angle+2]);
		}
		// Gas out, 1xCs, wo cut
		if ((tdc2[angle+2] > 400) && (ph2[angle+2] > phLimB1[angle+2]) && (ph2[angle+2] < 3800) && (psd2[angle+2] > psdLimB1[angle+2])){
			uncut1cs2[angle+2]->Fill(rfpick2 + rfOffIn[angle+2]);
		}
		// Gas out, 2xCs, wo cut
		if ((tdc2[angle+2] > 400) && (ph2[angle+2] > phLimB2[angle+2]) && (ph2[angle+2] < 3800) && (psd2[angle+2] > psdLimB2[angle+2])){
			uncut2cs2[angle+2]->Fill(rfpick2 + rfOffIn[angle+2]);
		}
		gasout_nocut[angle+2]->Fill(rfpick2 + rfOffIn[angle+2]);
	}

	/* Writes all the TOF spectra into one .root file. */
	f1->cd();

	for(i = angle; i <= (angle+2); i++){
		cut1cs[i]->Write();
		cut1cs2[i]->Write();
		cut2cs[i]->Write();
		uncut1cs[i]->Write();
		uncut2cs[i]->Write();
		cut2cs2[i]->Write();
		uncut1cs2[i]->Write();
		uncut2cs2[i]->Write();
		gasin_nocut[i]->Write();
		gasout_nocut[i]->Write();

	}


	c_gasin_comb->cd();

	cut1cs[angle]->Rebin(rebin_h);
	cut1cs[angle]->SetTitle("Gas In Top");
	cut1cs[angle]->SetLineColor(kRed);
	cut1cs[angle]->Draw();

	cut1cs[angle+1]->Rebin(rebin_h);
	cut1cs[angle+1]->SetTitle("Gas In Middle");
	cut1cs[angle+1]->SetLineColor(kBlue);
	cut1cs[angle+1]->Draw("SAME");

	cut1cs[angle+2]->Rebin(rebin_h);
	cut1cs[angle+2]->SetTitle("Gas In Bottom");
	cut1cs[angle+2]->SetLineColor(kGreen);
	cut1cs[angle+2]->Draw("SAME");


	c_gasin_comb_nocut->cd();

	uncut1cs[angle]->Rebin(rebin_h);
	uncut1cs[angle]->SetTitle("Gas In Top");
	uncut1cs[angle]->SetLineColor(kRed);
	uncut1cs[angle]->Draw();

	uncut1cs[angle+1]->Rebin(rebin_h);
	uncut1cs[angle+1]->SetTitle("Gas In Middle");
	uncut1cs[angle+1]->SetLineColor(kBlue);
	uncut1cs[angle+1]->Draw("SAME");

	uncut1cs[angle+2]->Rebin(rebin_h);
	uncut1cs[angle+2]->SetTitle("Gas In Bottom");
	uncut1cs[angle+2]->SetLineColor(kGreen);
	uncut1cs[angle+2]->Draw("SAME");


	c_gasout_comb->cd();

	cut1cs2[angle]->Rebin(rebin_h);
	cut1cs2[angle]->SetTitle("Gas Out Top");
	cut1cs2[angle]->SetLineColor(kRed);
	cut1cs2[angle]->Draw();

	cut1cs2[angle+1]->Rebin(rebin_h);
	cut1cs2[angle+1]->SetTitle("Gas Out Middle");
	cut1cs2[angle+1]->SetLineColor(kBlue);
	cut1cs2[angle+1]->Draw("SAME");

	cut1cs2[angle+2]->Rebin(rebin_h);
	cut1cs2[angle+2]->SetTitle("Gas Out Bottom");
	cut1cs2[angle+2]->SetLineColor(kGreen);
	cut1cs2[angle+2]->Draw("SAME");

	c_gasin_comb->BuildLegend();
	c_gasin_comb_nocut->BuildLegend();
	c_gasout_comb->BuildLegend();
	f1->cd();
	c_gasin_comb->Write();
	c_gasin_comb_nocut->Write();
	c_gasout_comb->Write();

	gROOT->SetBatch(kFALSE);

	cch1.ResetBranchAddresses();
	cch2.ResetBranchAddresses();

	f2->Close();

	return 0;
}

/* ******************************************************************************************************
 * This function applies fits to the TOF. Replaces hardtest_full.C. It accesses the fit functions
 * through fitfunctions.c, which was kept seperate due to its length.
 *
 * Inputs:
 *
 * Returns:
 *
 * ******************************************************************************************************/
int Apply_Fits(int angle, TFile* f2, FILE* f3){

	TFile*	fFit	= NULL;

	switch(angle){
		case 0:
		fFit = new TFile(TOFSPECTRA0, "READ");
		break;

		case 3:
		fFit = new TFile(TOFSPECTRA3, "READ");
		break;

		case 6:
		fFit = new TFile(TOFSPECTRA6, "READ");
		break;

		case 9:
		fFit = new TFile(TOFSPECTRA9, "READ");
		break;

		case 12:
		fFit = new TFile(TOFSPECTRA12, "READ");
		break;

		case 15:
		fFit = new TFile(TOFSPECTRA15, "READ");
		break;

		case 18:
		fFit = new TFile(TOFSPECTRA18, "READ");
		break;
		default: break;
	}

	// Check to make sure that the filename is valid
	if(fFit->IsZombie()){
		printf("There was an error opening the root file with the TOF spectra.\n");
		printf("Check if the file path is incorrect. Quitting application.\n\n");
		return 1;
	}

	switch(angle){
		case 0:
			fit0D(angle, f2, f3, fFit);
			break;
		case 3:
			fit3D(angle, f2, f3, fFit);
			break;
		case 6:
			fit6D(angle, f2, f3, fFit);
			break;
		case 9:
			fit9D(angle, f2, f3, fFit);
			break;
		case 12:
			fit12D(angle, f2, f3, fFit);
			break;
		case 15:
			fit15D(angle, f2, f3, fFit);
			break;
		case 18:
			fit18D(angle, f2, f3, fFit);
			break;
		default:
			return 1;
	}
	return 0;
}


/* ******************************************************************************************************
 * This function generates a user menu and allows them to choose what operation they want,
 * along with the angle they would like to work with.
 *
 * Inputs: The address for the angle input (angle) and the address for the menu choice (menCho).
 *
 * Returns: Nothing useful. The addresses are given to the function and the variables are altered
 * directly. Forced to return zero. No error reports needed at this time.
 *
 * *****************************************************************************************************/
int menu(int* angle, char* menCho){
	int 	menuFlag = 0;
	int 	badC = 1;
	int		i = 0;
	char 	angSTR[5];		// Buffer for checking angle input
	char 	input[25];
	char	menuVal = '6';


	/* Creates user menu. */
	printf("\n\nWelcome to the 11 MeV TUNL analysis of the Summer 2018 Argon run.\n");
	printf("The filled sample cell runs are 3672, 3673, 3674, and 3675. The empty run\n");
	printf("is 3702. Select an option below and then select the angle to be analyzed.\n");
	printf("Make sure you have set up all the paths in the script header file correctly.");
	printf("\nIF YOU HAVE NOT PROPERLY SET-UP filepaths.h, THEN QUIT SCRIPT AND DO SO.");

	printf("\n\n\
	[1] Generate PH versus PSD histogram\n\
	[2] Generate PH vs. PSD with cuts shown\n\
	[3] Generate TOF spectra\n\
	[4] Fit functions to combined top, middle, and bottom TOF spectra\n\
	[5] Quit Script\n\n");

	/* Handles the user input. Provides error if user doesn't input
	 * allowable choice. Crashes if input exceeds 24 characters.
	 */
	while ((menuVal < '1') || (menuVal > '5'))
	{
		printf("Enter choice (1-5): ");
		scanf("%s", input);
		sscanf(input, "%c", &menuVal);
		if ((menuVal < '1') || (menuVal > '5'))
		{
			printf("Sorry, invalid user input.\n");
		}
	}

	*menCho = menuVal;
	if(*menCho == '5'){
		return 0;
	}
	printf("You have chosen option %c, next select the angle.\n", *menCho);

	// This while loop gets the angle from the user. It handles incorrect inputs.
	while(badC){
		printf("\nPlease enter the angle (0,3,6,9,12,15,18). Type 21 to abort: ");
		scanf("%s", angSTR);
		if(sscanf(angSTR, "%d", angle)){
			for (i = 0; i <= 7; i++){
				if(*angle == 3*i){
					badC = 0;
				}
			}
		}
		if (badC == 1){
			printf("Invalid input, please try again.");
		}
	}

	printf("\n");

	return 0;
}

// End of script.
