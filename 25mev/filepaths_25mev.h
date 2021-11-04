
/* Here are all the paths you will need to set up. They are user specific.
 * SAVEAREA is where you want to write your results. CUTPATH is where you
 * have your cuts stored. Gx_TREEx are the paths to the trees. You may have
 * a variable number of trees. This code is only prepared to handle 11MeV data
 * at this particular moment.
 *
 */
#define XBIN 2

#define SAVEAREA	"/Users/joannawuko/bcs_25mev/"
#define CUTPATH 	"/Users/joannawuko/bcs_25mev/25mev_cuts_0802.root"
//#define TOFSPECTRA	"/Users/joannawuko/bcs_11mev/0degtofnew.root"


#define TOFSPECTRA0   "/Users/joannawuko/bcs_25mev/deg0_tof/deg0tof_0718_25mev.root"
#define TOFSPECTRA3   "/Users/joannawuko/bcs_25mev/deg3_tof/deg3tof_0718_25mev.root"
#define TOFSPECTRA6   "/Users/joannawuko/bcs_25mev/deg6_tof/deg6tof_0718_25mev.root"
#define TOFSPECTRA9   "/Users/joannawuko/bcs_25mev/deg9_tof/deg9tof_0718_25mev.root"
#define TOFSPECTRA12  "/Users/joannawuko/bcs_25mev/deg12_tof/deg12tof_0718_25mev.root"
#define TOFSPECTRA15  "/Users/joannawuko/bcs_25mev/deg15_tof/deg15tof_0718_25mev.root"
#define TOFSPECTRA18  "/Users/joannawuko/bcs_25mev/deg18_tof/deg18tof_0718_25mev.root"


#define GI_TREE0 	"/Users/joannawuko/bcs_25mev/trees/nme_3681_nt.root/ntuple;14"
#define GI_TREE1	"/Users/joannawuko/bcs_25mev/trees/nme_3682_nt.root/ntuple;8"
#define GI_TREE2	"/Users/joannawuko/bcs_25mev/trees/nme_3687_nt.root/ntuple;14"
#define GI_TREE3	"/Users/joannawuko/bcs_25mev/trees/nme_3688_nt.root/ntuple;6"
#define GO_TREE0	"/Users/joannawuko/bcs_25mev/trees/nme_3690_nt.root/ntuple;12"
#define GO_TREE1	"/Users/joannawuko/bcs_25mev/trees/nme_3691_nt.root/ntuple;9"
#define GO_TREE2	"/Users/joannawuko/bcs_25mev/trees/nme_3692_nt.root/ntuple;11"
#define GO_TREE3	"/Users/joannawuko/bcs_25mev/trees/nme_3693_nt.root/ntuple;10"
#define GO_TREE4	"/Users/joannawuko/bcs_25mev/trees/nme_3694_nt.root/ntuple;8"

/* These paths are where you want to save your ph versus psd 2D plots. */
#define PHPSD		"/Users/joannawuko/bcs_25mev/phpsd_25mev_0718.root"
#define PHPSDCUT	"/Users/joannawuko/bcs_25mev/25mev_phpsdwcut_0802.root"

/* Here you must put your cut names, in the correct order.
 *
 */
#define TCUTGS		{"deg0Tcut_25mev_0802", "deg0Mcut_25mev_0802", "deg0Bcut_25mev_0802", "deg3Tcut_25mev_0718",\
                     "deg3Mcut_25mev_0718", "deg3Bcut_25mev_0718", "deg6Tcut_25mev_0718", "deg6Mcut_25mev_0718",\
                     "deg6Bcut_25mev_0718", "deg9Tcut_25mev_0718", "deg9Mcut_25mev_0718", "deg9Bcut_25mev_0718",\
                     "deg12Tcut_25mev_0718", "deg12Mcut_25mev_0718", "deg12Bcut_25mev_0718", "deg15Tcut_25mev_0718",\
                     "deg15Mcut_25mev_0718","deg15Bcut_25mev_0718", "deg18Tcut_25mev_0718", "deg18Mcut_25mev_0718",\
                     "deg18Bcut_25mev_0718"}


/* This section of code does not need to be changed in order to run the script on your
 * own workspace. They can be changed to modify the program if need be. If you are not well
 * versed in programming, I would not recommend changing these without supervision.
 *
 */

#define HIST2D		{"phPSD0Tin", "phPSD0Min", "phPSD0Bin", "phPSD3Tin", "phPSD3Min", "phPSD3Bin",\
					 "phPSD6Tin", "phPSD6Min", "phPSD6Bin", "phPSD9Tin", "phPSD9Min", "phPSD9Bin",\
					 "phPSD12Tin", "phPSD12Min", "phPSD12Bin", "phPSD15Tin", "phPSD15Min", "phPSD15Bin",\
					 "phPSD18Tin", "phPSD18Min", "phPSD18Bin"}

#define CUTIN		{"TOF_1Cs_0T_CUT_IN", "TOF_1Cs_0M_CUT_IN", "TOF_1Cs_0B_CUT_IN", "TOF_1Cs_3T_CUT_IN", "TOF_1Cs_3M_CUT_IN",\
					 "TOF_1Cs_3B_CUT_IN", "TOF_1Cs_6T_CUT_IN", "TOF_1Cs_6M_CUT_IN", "TOF_1Cs_6B_CUT_IN", "TOF_1Cs_9T_CUT_IN",\
					 "TOF_1Cs_9M_CUT_IN", "TOF_1Cs_9B_CUT_IN", "TOF_1Cs_12T_CUT_IN", "TOF_1Cs_12M_CUT_IN", "TOF_1Cs_12B_CUT_IN",\
					 "TOF_1Cs_15T_CUT_IN", "TOF_1Cs_15M_CUT_IN", "TOF_1Cs_15B_CUT_IN", "TOF_1Cs_18T_CUT_IN", "TOF_1Cs_18M_CUT_IN",\
					 "TOF_1Cs_18B_CUT_IN"}


#define CUTOUT		{"TOF_1Cs_0T_CUT_OUT", "TOF_1Cs_0M_CUT_OUT", "TOF_1Cs_0B_CUT_OUT", "TOF_1Cs_3T_CUT_OUT", "TOF_1Cs_3M_CUT_OUT",\
					 "TOF_1Cs_3B_CUT_OUT", "TOF_1Cs_6T_CUT_OUT", "TOF_1Cs_6M_CUT_OUT", "TOF_1Cs_6B_CUT_OUT", "TOF_1Cs_9T_CUT_OUT",\
					 "TOF_1Cs_9M_CUT_OUT", "TOF_1Cs_9B_CUT_OUT", "TOF_1Cs_12T_CUT_OUT", "TOF_1Cs_12M_CUT_OUT", "TOF_1Cs_12B_CUT_OUT",\
					 "TOF_1Cs_15T_CUT_OUT", "TOF_1Cs_15M_CUT_OUT", "TOF_1Cs_15B_CUT_OUT", "TOF_1Cs_18T_CUT_OUT", "TOF_1Cs_18M_CUT_OUT",\
					 "TOF_1Cs_18B_CUT_OUT"}

/* Here we have the limits given by TUNL labratory. Can be found in the PAW
 * configuration file. It goes 0T, 0M, 0B, 3T, 3M, 3B, etc.
 * The PSD top limit was left out and not used.
 *
 */
#define TDCLIM		400
#define	PHLIMT		3800

//1xcs
#define	PHLIMB1C	{280, 315, 293, 328, 431, 361, 384, 353, 243, 237, 386, 220,\
					 290, 266, 316, 329, 221, 331, 268, 268, 268}

#define	PSDLIMB1C	{1700, 1785, 1280, 1038, 1125, 1075, 1125, 1160, 2090,\
					 2120, 2090, 2190, 1780, 1730, 1565, 1650, 1305, 1341,\
					 1350, 1350, 1350}
//2xcs
#define	PHLIMB2C	{432, 449, 422, 475, 692, 563, 593, 517, 327, 334, 501, 276,\
					 406, 356, 453, 458, 348, 552, 410, 402, 432}

#define	PSDLIMB2C	{1700, 1785, 1280, 1038, 1125, 1075, 1125, 1160, 2090,\
					 2120, 2090, 2190, 1780, 1730, 1565, 1650, 1305, 1341,\
					 1350, 1350, 1350}

/* Note that there is an anomaly with the 0B rfpick offset for gas in.
 * There is a -14 in place of one of the -10 for 1Cs uncut. This is how
 * Ryans script is different than the current script. */
//#define	RFINOFF		{0, 26, -10, (-1+10), -24, 18, -22, 8, -1, 0, -26, 0, 1, 0, -47,\
					 -21, 28, -3, 0, 0, 0}

/* Note that the offset is not applied to the uncut TOF for 3M and 3B.
 * This is how ryans script is different than the current script. */
//#define RFOUTOFF	{0, 26, -10, (-1+10), -24, 18, -22, 8, -1, 0, -26, 0, 1, 0, -47,\
					 -21, 28, -3, 0, 0, 0}

//***********************************************************************

//offsets for gamma peak alignment


//#define	RFINOFF		{0, 26, -10, 9, -17, 25, 0, -2, -1, 0, -26, 0, 1, 0, -47,\
            -21, 28, -3, 0, 0, 0}


//#define RFOUTOFF	{0, 26, -10, 9, -17, 25, 0, -2, -1, 0, -26, 0, 1, 0, -47,\
            -21, 28, -3, 0, 0, 0}


//***********************************************************************

//updated test offsets


#define	RFINOFF		{0, 26, -10, 0, -25, 17, -22, 9, 0, -12, 0, 27, 1, 0, -45,\
              -20, 29, 0, 0, 0, 0}


#define RFOUTOFF	{0, 26, -10, 0, -25, 17, -22, 9, 0, -12, 0, 27, 1, 0, -45,\
              -20, 29, 0, 0, 0, 0}



//***********************************************************************

//0 offsets


//#define	RFINOFF		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\
              0, 0, 0, 0, 0, 0}


//#define RFOUTOFF	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\
              0, 0, 0, 0, 0, 0}




/* Do not change this section of code. These are function definitions.
 * They are required for the script to run. There should be no reason to need to change this code.
 */
int hist_2D(void);
int Paint_Cuts(void);
int Generate_TOF(int angle, TFile* f1, FILE* f3);
int Apply_Fits(int angle, TFile* f2, FILE* f3);
int menu(int* angle, char* menCho);
