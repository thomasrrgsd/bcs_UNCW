/*
 *
 * fitfunctions.c
 *
 * William T. Jarratt
 * Root V6.16
 * Version 1.00
 * Last Updated:  14 FEB 2021
 *
 * This script contains the code used for applying fits to the combined TOF spectra.
 * The code is directly ported from Ryan's framework and is still being modified
 * to work out a few bugs.
 *
 */

#include <TMath.h>

// Quadratic background function.
Double_t background(Double_t *x, Double_t *par) {
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Gaussian Peak function.
Double_t gaussianpeak(Double_t *x, Double_t *par) {
	Double_t arg1 = 0;
	Double_t arg2 = 0;
	Double_t arg3 = 0;
	Double_t arg4 = 0;
	Double_t arg5 = 0;
	Double_t arg6 = 0;
	Double_t arg7 = 0;
	if (par[0]) {
	  arg1 = (x[0] - par[1])/par[0];
	  arg2 = (x[0] - par[1] + 41.4)/par[0];
	  arg3 = (x[0] - par[1] + 50.5)/par[0];
	  arg4 = (x[0] - par[1] + 68.2)/par[0];
	  arg5 = (x[0] - par[1] + 96.1)/par[0];
		arg6 = (x[0] - par[1] + 108)/par[0];
		arg7 = (x[0] - par[1] + 132)/par[0];
	}
	//return 	par[2]*TMath::Exp(-0.5*arg1*arg1) +
		//par[3]*TMath::Exp(-0.5*arg2*arg2) +
		//par[4]*TMath::Exp(-0.5*arg3*arg3) +
		//par[5]*TMath::Exp(-0.5*arg4*arg4) +
		//par[6]*TMath::Exp(-0.5*arg5*arg5);


		return 	par[2]*TMath::Exp(-0.5*arg1*arg1) + 	// ground state amplitude
			par[2]*par[3]*TMath::Exp(-0.5*arg2*arg2) +
			par[2]*par[4]*TMath::Exp(-0.5*arg3*arg3) +
			par[2]*par[5]*TMath::Exp(-0.5*arg4*arg4) +
			par[2]*par[6]*TMath::Exp(-0.5*arg5*arg5) +
			par[2]*par[7]*TMath::Exp(-0.5*arg6*arg6) +
			par[2]*par[8]*TMath::Exp(-0.5*arg7*arg7);
}

// Gaussian Peak function.
Double_t gaussianpeak1(Double_t *x, Double_t *par) {
	Double_t arg1 = 0;
	if (par[0])
	  arg1 = (x[0] - par[1])/par[0];
	return par[2]*TMath::Exp(-0.5*arg1*arg1);
}


// Combination of gaussian peak and background function.
Double_t fitfuncc(Double_t *x, Double_t *par) {
	return gaussianpeak(x,&par[3])+ background(x,par);
}


/*
 * Fit function for 0 degrees.
 *
 */
int fit0D(int angle, TFile* f2, FILE* f3, TFile* fFit){

	TH1F	*tofTop = NULL;
	TH1F	*tofMid = NULL;
	TH1F	*tofBot = NULL;
	TH1F	*tofTopOut = NULL;
	TH1F	*tofMidOut = NULL;
	TH1F	*tofBotOut = NULL;
	char 	cutinSTR[][21] = CUTIN;
	char 	cutoutSTR[][21] = CUTOUT;
	char	txtSTR[100];

	// rebinned (div by x)
	float h = 4;

	gROOT->SetBatch(kTRUE);

	fprintf(f3, "\nYou have chosen to apply a fit to the combined TOF Spectra.\n");

	/* This portion of code grabs the TOF objects that will be fitted. */
	fFit->GetObject(cutinSTR[angle],tofTop);
	fFit->GetObject(cutinSTR[angle+1],tofMid);
	fFit->GetObject(cutinSTR[angle+2],tofBot);
	fFit->GetObject(cutoutSTR[angle],tofTopOut);
	fFit->GetObject(cutoutSTR[angle+1],tofMidOut);
	fFit->GetObject(cutoutSTR[angle+2],tofBotOut);

	TH1F *mycutzu1_clone = (TH1F*)tofTop->Clone();
	mycutzu1_clone->SetName("mycutzu1_clone");
	mycutzu1_clone->SetLineColor(kMagenta);
	mycutzu1_clone->Add(tofMid);
	mycutzu1_clone->Add(tofBot);
	mycutzu1_clone->Rebin(h);
	sprintf(txtSTR, "Combined %d Degree TOF with Applied Fits", angle);
	mycutzu1_clone->SetTitle(txtSTR);

	TH1F *gasout1_1clone = (TH1F*)tofTopOut->Clone();
	gasout1_1clone->SetName("gasout1_1clone");
	gasout1_1clone->SetLineColor(kAzure);
	gasout1_1clone->Add(tofMidOut);
	gasout1_1clone->Add(tofBotOut);
	gasout1_1clone->Rebin(h);


	gasout1_1clone->Scale(3.62); //product of ratio of gas in to gas out bcis and ratio of gas in to gas out live times

	f2->cd();

	/* Fit range parameters */

	Int_t lu1 = 2406;
	Int_t hu1 = 2600;


	sprintf(txtSTR, "1xCs_degree_%d_fitted_data", angle);
	TCanvas *c1 = new TCanvas(txtSTR);
	c1->cd();

	TH1F *gasin_subtract = (TH1F*)mycutzu1_clone->Clone();
	gasin_subtract->Sumw2();
	gasout1_1clone->Sumw2();
	gasin_subtract->Add(gasout1_1clone,-1);

//*******new things that were added**********************************

	//sprintf(txtSTR, "1xCs_degree_%d_fitted_data", angle);
	//TCanvas *c1temp = new TCanvas("temp");
	//c1temp->cd();

	//TH1F *gasin_subtract = (TH1F*)mycutzu1_clone->Clone();

	//gasin_subtract->Sumw2();
	//gasout1_1clone->Sumw2();
	//gasin_subtract->Add(gasout1_1clone,-1);

	//gasout1_1clone->Draw();
	//gasin_subtract->Draw("SAME");
	//mycutzu1_clone->Draw("SAME");
	//c1temp->Write();

	//TCanvas *c1 = new TCanvas(txtSTR);
	//c1->cd();

//********************************************************************

	/*
	 * Creation of fit parameters
	 */

	Double_t params[12] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

	TF1 *fitfunc2 = new TF1("fitfunc2", fitfuncc, lu1, hu1, 12);
	fitfunc2->SetParameter(0,(5.1));  		//constant
	fitfunc2->FixParameter(1,0);  		//linear
	fitfunc2->FixParameter(2,0);  		//quad
	fitfunc2->SetParameter(3,11.7);  	//sigma
	fitfunc2->SetParameter(4,2513.0);  //mean
	fitfunc2->SetParameter(5,206.8);   //amp1

	//***************************************

	//parameter 6, first 2+ state
	//uncomment to include

	fitfunc2->SetParameter(6,0.001);     //amp2   2+ states
	fitfunc2->SetParLimits(6,0,1000);

	//uncomment to exclude

	//fitfunc2->FixParameter(6,0.0);

	//***************************************

	fitfunc2->SetParameter(7,0.13);    	//amp3  50->0.1
	fitfunc2->SetParameter(8,0.001);     //amp4   2+ states
	fitfunc2->SetParameter(9,0.11);    	//amp5  50->0.1
	fitfunc2->FixParameter(10,0.0);			//amp6
	fitfunc2->FixParameter(11,0.0);			//amp7


	fitfunc2->SetParLimits(0,0,40);
	fitfunc2->SetParLimits(5,0,1000);

	fitfunc2->SetParLimits(7,0,1000);
	fitfunc2->SetParLimits(8,0,1000);
	fitfunc2->SetParLimits(9,0,1000);
	//fitfunc2->SetParLimits(10,0,1000);
	//fitfunc2->SetParLimits(11,0,1000);

	TF1 *gausfuncu1 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu11 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu12 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu13 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu14 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu15 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu16 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);

	TF1 *bkdfuncu1 = new TF1("bkdfunc",background, lu1, hu1, 3);

	printf("\n\n");
	TFitResultPtr fit_result = mycutzu1_clone->Fit("fitfunc2","RS");

	fitfunc2->GetParameters(params);

	bkdfuncu1->SetParameter(0,params[0]);   //constant
	bkdfuncu1->SetParameter(1,params[1]);   //linear
	bkdfuncu1->SetParameter(2,params[2]);   //quad

	gausfuncu1->SetParameter(0,params[3]);  //sigma
	gausfuncu1->SetParameter(1,params[4]);  //mean
	gausfuncu1->SetParameter(2,params[5]);  //amp


	gausfuncu11->SetParameter(0,params[3]);  //Sigma
	gausfuncu11->SetParameter(1,params[4] - 41.4);  //amp2
	gausfuncu11->SetParameter(2,params[5]*params[6]);  //amp3  //*p5

	gausfuncu12->SetParameter(0,params[3]);  //Sigma
	gausfuncu12->SetParameter(1,params[4] - 50.5);  //Mean
	gausfuncu12->SetParameter(2,params[5]*params[7]);  //amp5

	gausfuncu13->SetParameter(0,params[3]);  //Sigma
	gausfuncu13->SetParameter(1,params[4] - 68.2);  //Mean
	gausfuncu13->SetParameter(2,params[5]*params[8]);  //amp5

	gausfuncu14->SetParameter(0,params[3]);  //Sigma
	gausfuncu14->SetParameter(1,params[4] - 96.1);  //Mean
	gausfuncu14->SetParameter(2,params[5]*params[9]);   //amp5

	gausfuncu15->SetParameter(0,params[3]);  //Sigma
	gausfuncu15->SetParameter(1,params[4] - 108);  //Mean
	gausfuncu15->SetParameter(2,params[5]*params[10]);   //amp6

	gausfuncu16->SetParameter(0,params[3]);  //Sigma
	gausfuncu16->SetParameter(1,params[4] - 132);  //Mean
	gausfuncu16->SetParameter(2,params[5]*params[11]);   //amp6


	//legend titles

	mycutzu1_clone->SetBit(TH1::kNoTitle);


	mycutzu1_clone->SetTitle("Combined Gas In");
	mycutzu1_clone->SetLineColor(kMagenta);
	mycutzu1_clone->Draw();
	gStyle->SetOptFit(1);

	gasout1_1clone->SetTitle("Combined Gas Out");
	gasout1_1clone->SetLineColor(kYellow + 1);
	gasout1_1clone->Draw("same");

	gasin_subtract->SetTitle("Gas In Minus Gas Out");
	gasin_subtract->SetLineColor(kBlue - 7);
	gasin_subtract->Draw("SAME");					// Using "SAMES" changes the info box to this fit.


	fitfunc2->SetTitle("Complete Fit Function");
	fitfunc2->SetLineColor(kRed);
	fitfunc2->Draw("same");

	bkdfuncu1->SetTitle("Background");
	bkdfuncu1->SetLineColor(kGreen - 3);
	bkdfuncu1->Draw("same");

	gausfuncu1->SetTitle("Ground State");
	gausfuncu1->SetLineColor(kAzure);
	gausfuncu1->Draw("same");

	gausfuncu11->SetTitle("2+ First Excited State 1.52 MeV");
	gausfuncu11->SetLineColor(kCyan);
	gausfuncu11->Draw("same");

	gausfuncu12->SetTitle("0+ First Excited State 1.84 MeV");
	gausfuncu12->SetLineColor(kGreen);
	gausfuncu12->Draw("same");

	gausfuncu13->SetTitle("2+ Second Excited State 2.42 MeV");
	gausfuncu13->SetLineColor(kOrange - 3);
	gausfuncu13->Draw("same");

	gausfuncu14->SetTitle("0+ Second Excited State 3.30 MeV");
	gausfuncu14->SetLineColor(kAzure+9);
	gausfuncu14->Draw("same");

	gausfuncu15->SetTitle("3.65 MeV");
	gausfuncu15->SetLineColor(kRed+4);
	gausfuncu15->Draw("same");

	gausfuncu16->SetTitle("??? MeV/126 TDC offset");
	gausfuncu16->SetLineColor(kCyan-2);
	gausfuncu16->Draw("same");


	//title: formatting needs to be fixed
	TPaveText *t = new TPaveText(0.0, 0.9, 0.3, 1.0, "brNDC"); // left-up
	t->AddText("0 Degree Fit");
	t->SetBorderSize(0);
	t->SetFillColor(gStyle->GetTitleFillColor());
	t->Draw();

	fit_result->Print("V");

	c1->BuildLegend();
	c1->Update();
	c1->Write();

	gROOT->SetBatch(kFALSE);

	double chisk = fit_result->Chi2();
	unsigned int ndf = fit_result->Ndf();
	double edm = fit_result->Edm();
	unsigned int nCall = fit_result->NCalls();
	int i,j;

	double para10[10];
	double erro10[10];
	for(i=0; i<10; i++){
		para10[i] = fit_result->Parameter(i);
		erro10[i] = fit_result->Error(i);
	}

	double cor66[10][10];
	double cov66[10][10];

	for(i=0; i<10; i++){
		for(j=0; j<10; j++){
			cor66[i][j] = fit_result->Correlation(i, j);
			cov66[i][j] = fit_result->CovMatrix(i, j);
		}
	}

	fprintf(f3, "\nChi2\t=\t%.3e", chisk);
	fprintf(f3, "\nNDf\t=\t%u", ndf);
	fprintf(f3, "\nEdm\t=\t%.3e", edm);
	fprintf(f3, "\nNCalls\t=\t%u", nCall);
	fprintf(f3, "\np0\t=\t%.3e +/- %.3e", para10[0], erro10[0]);
	fprintf(f3, "\np1\t=\t%.3e +/- %.3e", para10[1], erro10[1]);
	fprintf(f3, "\np2\t=\t%.3e +/- %.3e", para10[2], erro10[2]);
	fprintf(f3, "\np3\t=\t%.3e +/- %.3e", para10[3], erro10[3]);
	fprintf(f3, "\np4\t=\t%.3e +/- %.3e", para10[4], erro10[4]);
	fprintf(f3, "\np5\t=\t%.3e +/- %.3e", para10[5], erro10[5]);
	fprintf(f3, "\np6\t=\t%.3e +/- %.3e", para10[6], erro10[6]);
	fprintf(f3, "\np7\t=\t%.3e +/- %.3e", para10[7], erro10[7]);
	fprintf(f3, "\np8\t=\t%.3e +/- %.3e", para10[8], erro10[8]);
	fprintf(f3, "\np9\t=\t%.3e +/- %.3e", para10[9], erro10[9]);


	/* Propagation of errors begins here */
	printf("\n\nPropagation of Errors: \n\n");
	printf("p3 --> %f \n", para10[3]);
	printf("p5 --> %f \n", para10[5]);
	printf("p7 --> %f \n", para10[7]);
	printf("p8 --> %f \n", para10[8]);
	printf("p9 --> %f \n", para10[9]);

	printf("p3u --> %f \n", cov66[3][3]);
	printf("p5u --> %f \n", cov66[5][5]);
	printf("p7u --> %f \n", cov66[7][7]);
	printf("p8u --> %f \n", cov66[8][8]);
	printf("p9u --> %f \n", cov66[9][9]);

	printf("p3p5 --> %f \n", cov66[3][5]);
	printf("p3p7 --> %f \n", cov66[3][7]);
	printf("p3p8 --> %f \n", cov66[3][8]);
	printf("p3p9 --> %f \n", cov66[3][9]);
	printf("p5p7 --> %f \n", cov66[5][7]);
	printf("p5p8 --> %f \n", cov66[5][8]);
	printf("p5p9 --> %f \n", cov66[5][9]);

	double b = TMath::Sqrt(2*TMath::Pi());
	double a1, a1u, a2, a2u, a3, a3u, a4, a4u, a5, a5u;

	a1 = para10[3]*para10[5]*b;
	a1u = TMath::Power((b*para10[5]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]), 2)*cov66[5][5] +
						2*b*b*para10[5]*para10[3]*cov66[3][5];
  a1 = a1/4.0;
  a1u = (TMath::Sqrt(a1u))/4.0;
	printf("\n\nGround State Area: %.4f +/- %.4f", a1, a1u);


	a2 = para10[3]*para10[5]*para10[7]*b;;
	a2u = TMath::Power((b*para10[5]*para10[7]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[7]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[7][7] +
					2*b*b*para10[5]*para10[3]*para10[7]*para10[7]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[7]*cov66[5][7] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[7]*cov66[3][7];
	a2 = a2/4.0;
	a2u = (TMath::Sqrt(a2u))/4.0;
	printf("\n\n0+ First Excited State Area: %.4f +/- %.4f", a2, a2u);

	a3 = para10[3]*para10[5]*para10[8]*b;;
	a3u = TMath::Power((b*para10[5]*para10[8]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[8]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[8][8] +
					2*b*b*para10[5]*para10[3]*para10[8]*para10[8]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[8]*cov66[5][8] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[8]*cov66[3][8];
	a3 = a3/4.0;
	a3u = (TMath::Sqrt(a3u))/4.0;
	printf("\n\n2+ Second Excited State Area: %.4f +/- %.4f", a3, a3u);

	a4 = para10[3]*para10[5]*para10[9]*b;;
	a4u = TMath::Power((b*para10[5]*para10[9]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[9]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[9][9] +
					2*b*b*para10[5]*para10[3]*para10[9]*para10[9]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[9]*cov66[5][9] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[9]*cov66[3][9];
	a4 = a4/4.0;
	a4u = (TMath::Sqrt(a4u))/4.0;
	printf("\n\n0+ Second Excited State Area: %.4f +/- %.4f", a4, a4u);

	a5 = para10[3]*para10[5]*para10[6]*b;;
	a5u = TMath::Power((b*para10[5]*para10[6]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[6]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[6][6] +
					2*b*b*para10[5]*para10[3]*para10[6]*para10[6]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[6]*cov66[5][6] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[6]*cov66[3][6];
	a5 = a5/4.0;
	a5u = (TMath::Sqrt(a5u))/4.0;
	printf("\n\n2+ First Excited State Area: %.4f +/- %.4f", a5, a5u);



	printf("\n\n");
	fprintf(f3, "\n\n");
	fclose(f3);
	return 0;
}


/*
 * Fit function for 3 degrees.
 *
 */
int fit3D(int angle, TFile* f2, FILE* f3, TFile* fFit){

	TH1F	*tofTop = NULL;
	TH1F	*tofMid = NULL;
	TH1F	*tofBot = NULL;
	TH1F	*tofTopOut = NULL;
	TH1F	*tofMidOut = NULL;
	TH1F	*tofBotOut = NULL;
	char 	cutinSTR[][21] = CUTIN;
	char 	cutoutSTR[][21] = CUTOUT;
	char	txtSTR[100];

	// rebinned (div by x)
	float h = 4;

	gROOT->SetBatch(kTRUE);

	fprintf(f3, "\nYou have chosen to apply a fit to the combined TOF Spectra.\n");

	/* This portion of code grabs the TOF objects that will be fitted. */
	fFit->GetObject(cutinSTR[angle],tofTop);
	fFit->GetObject(cutinSTR[angle+1],tofMid);
	fFit->GetObject(cutinSTR[angle+2],tofBot);
	fFit->GetObject(cutoutSTR[angle],tofTopOut);
	fFit->GetObject(cutoutSTR[angle+1],tofMidOut);
	fFit->GetObject(cutoutSTR[angle+2],tofBotOut);

	TH1F *mycutzu1_clone = (TH1F*)tofTop->Clone();
	mycutzu1_clone->SetName("mycutzu1_clone");
	mycutzu1_clone->SetLineColor(kMagenta);
	mycutzu1_clone->Add(tofMid);
	mycutzu1_clone->Add(tofBot);
	mycutzu1_clone->Rebin(h);
	sprintf(txtSTR, "Combined %d Degree TOF with Applied Fits", angle);
	mycutzu1_clone->SetTitle(txtSTR);

	TH1F *gasout1_1clone = (TH1F*)tofTopOut->Clone();
	gasout1_1clone->SetName("gasout1_1clone");
	gasout1_1clone->SetLineColor(kAzure);
	gasout1_1clone->Add(tofMidOut);
	gasout1_1clone->Add(tofBotOut);
	gasout1_1clone->Rebin(h);


	gasout1_1clone->Scale(3.62); //product of ratio of gas in to gas out bcis and ratio of gas in to gas out live times

	f2->cd();

	/* Fit range parameters */

	Int_t lu1 = 2353;
	Int_t hu1 = 2558;


	sprintf(txtSTR, "1xCs_degree_%d_fitted_data", angle);
	TCanvas *c1 = new TCanvas(txtSTR);
	c1->cd();

	TH1F *gasin_subtract = (TH1F*)mycutzu1_clone->Clone();
	gasin_subtract->Sumw2();
	gasout1_1clone->Sumw2();
	gasin_subtract->Add(gasout1_1clone,-1);

	/*
	 * Creation of fit parameters
	 */

	Double_t params[12] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

	TF1 *fitfunc2 = new TF1("fitfunc2", fitfuncc, lu1, hu1, 12);
	fitfunc2->SetParameter(0,(8.8));  		//constant
	fitfunc2->FixParameter(1,0);  		//linear
	fitfunc2->FixParameter(2,0);  		//quad
	fitfunc2->SetParameter(3,11.6);  	//sigma
	fitfunc2->SetParameter(4,2461.0);  //mean
	fitfunc2->SetParameter(5,165.9);   //amp1
	fitfunc2->SetParameter(6,0.1);     //amp2   2+ states
	fitfunc2->SetParameter(7,0.15);    	//amp3  50->0.1
	fitfunc2->SetParameter(8,0.001);     //amp4   2+ states
	fitfunc2->SetParameter(9,0.09);    	//amp5  50->0.1
	fitfunc2->FixParameter(10,0.0);			//amp6
	fitfunc2->FixParameter(11,0.0);			//amp7

	fitfunc2->SetParLimits(0,0,40);
	fitfunc2->SetParLimits(5,0,1000);
	fitfunc2->SetParLimits(6,0,1000);
	fitfunc2->SetParLimits(7,0,1000);
	fitfunc2->SetParLimits(8,0,1000);
	fitfunc2->SetParLimits(9,0,1000);
	//fitfunc2->SetParLimits(10,0,1000);
	//fitfunc2->SetParLimits(11,0,1000);

	TF1 *gausfuncu1 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu11 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu12 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu13 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu14 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu15 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu16 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);

	TF1 *bkdfuncu1 = new TF1("bkdfunc",background, lu1, hu1, 3);

	printf("\n\n");
	TFitResultPtr fit_result = mycutzu1_clone->Fit("fitfunc2","RS");

	fitfunc2->GetParameters(params);

	bkdfuncu1->SetParameter(0,params[0]);   //constant
	bkdfuncu1->SetParameter(1,params[1]);   //linear
	bkdfuncu1->SetParameter(2,params[2]);   //quad

	gausfuncu1->SetParameter(0,params[3]);  //sigma
	gausfuncu1->SetParameter(1,params[4]);  //mean
	gausfuncu1->SetParameter(2,params[5]);  //amp


	gausfuncu11->SetParameter(0,params[3]);  //Sigma
	gausfuncu11->SetParameter(1,params[4] - 41.4);  //amp2
	gausfuncu11->SetParameter(2,params[5]*params[6]);  //amp3  //*p5

	gausfuncu12->SetParameter(0,params[3]);  //Sigma
	gausfuncu12->SetParameter(1,params[4] - 50.5);  //Mean
	gausfuncu12->SetParameter(2,params[5]*params[7]);  //amp5

	gausfuncu13->SetParameter(0,params[3]);  //Sigma
	gausfuncu13->SetParameter(1,params[4] - 68.2);  //Mean
	gausfuncu13->SetParameter(2,params[5]*params[8]);  //amp5

	gausfuncu14->SetParameter(0,params[3]);  //Sigma
	gausfuncu14->SetParameter(1,params[4] - 96.1);  //Mean
	gausfuncu14->SetParameter(2,params[5]*params[9]);   //amp5

	gausfuncu15->SetParameter(0,params[3]);  //Sigma
	gausfuncu15->SetParameter(1,params[4] - 108);  //Mean
	gausfuncu15->SetParameter(2,params[5]*params[10]);   //amp6

	gausfuncu16->SetParameter(0,params[3]);  //Sigma
	gausfuncu16->SetParameter(1,params[4] - 126);  //Mean
	gausfuncu16->SetParameter(2,params[5]*params[11]);   //amp6


	//legend titles

	mycutzu1_clone->SetBit(TH1::kNoTitle);


	mycutzu1_clone->SetTitle("Combined Gas In");
	mycutzu1_clone->SetLineColor(kMagenta);
	mycutzu1_clone->Draw();
	gStyle->SetOptFit(1);

	gasout1_1clone->SetTitle("Combined Gas Out");
	gasout1_1clone->SetLineColor(kYellow + 1);
	gasout1_1clone->Draw("same");

	gasin_subtract->SetTitle("Gas In Minus Gas Out");
	gasin_subtract->SetLineColor(kBlue - 7);
	gasin_subtract->Draw("SAME");					// Using "SAMES" changes the info box to this fit.


	fitfunc2->SetTitle("Complete Fit Function");
	fitfunc2->SetLineColor(kRed);
	fitfunc2->Draw("same");

	bkdfuncu1->SetTitle("Background");
	bkdfuncu1->SetLineColor(kGreen - 3);
	bkdfuncu1->Draw("same");

	gausfuncu1->SetTitle("Ground State");
	gausfuncu1->SetLineColor(kAzure);
	gausfuncu1->Draw("same");

	gausfuncu11->SetTitle("2+ First Excited State 1.52 MeV");
	gausfuncu11->SetLineColor(kCyan);
	gausfuncu11->Draw("same");

	gausfuncu12->SetTitle("0+ First Excited State 1.84 MeV");
	gausfuncu12->SetLineColor(kGreen);
	gausfuncu12->Draw("same");

	gausfuncu13->SetTitle("2+ Second Excited State 2.42 MeV");
	gausfuncu13->SetLineColor(kOrange - 3);
	gausfuncu13->Draw("same");

	gausfuncu14->SetTitle("0+ Second Excited State 3.30 MeV");
	gausfuncu14->SetLineColor(kAzure+9);
	gausfuncu14->Draw("same");

	gausfuncu15->SetTitle("3.65 MeV");
	gausfuncu15->SetLineColor(kRed+4);
	gausfuncu15->Draw("same");

	gausfuncu16->SetTitle("??? MeV/126 TDC offset");
	gausfuncu16->SetLineColor(kCyan-2);
	gausfuncu16->Draw("same");

	//title: formatting needs to be fixed
	TPaveText *t = new TPaveText(0.0, 0.9, 0.3, 1.0, "brNDC"); // left-up
	t->AddText("3 Degree Fit");
	t->SetBorderSize(0);
	t->SetFillColor(gStyle->GetTitleFillColor());
	t->Draw();

	fit_result->Print("V");

	c1->BuildLegend();
	c1->Update();
	c1->Write();

	gROOT->SetBatch(kFALSE);

	double chisk = fit_result->Chi2();
	unsigned int ndf = fit_result->Ndf();
	double edm = fit_result->Edm();
	unsigned int nCall = fit_result->NCalls();
	int i,j;

	double para10[10];
	double erro10[10];
	for(i=0; i<10; i++){
		para10[i] = fit_result->Parameter(i);
		erro10[i] = fit_result->Error(i);
	}

	double cor66[10][10];
	double cov66[10][10];

	for(i=0; i<10; i++){
		for(j=0; j<10; j++){
			cor66[i][j] = fit_result->Correlation(i, j);
			cov66[i][j] = fit_result->CovMatrix(i, j);
		}
	}

	fprintf(f3, "\nChi2\t=\t%.3e", chisk);
	fprintf(f3, "\nNDf\t=\t%u", ndf);
	fprintf(f3, "\nEdm\t=\t%.3e", edm);
	fprintf(f3, "\nNCalls\t=\t%u", nCall);
	fprintf(f3, "\np0\t=\t%.3e +/- %.3e", para10[0], erro10[0]);
	fprintf(f3, "\np1\t=\t%.3e +/- %.3e", para10[1], erro10[1]);
	fprintf(f3, "\np2\t=\t%.3e +/- %.3e", para10[2], erro10[2]);
	fprintf(f3, "\np3\t=\t%.3e +/- %.3e", para10[3], erro10[3]);
	fprintf(f3, "\np4\t=\t%.3e +/- %.3e", para10[4], erro10[4]);
	fprintf(f3, "\np5\t=\t%.3e +/- %.3e", para10[5], erro10[5]);
	fprintf(f3, "\np6\t=\t%.3e +/- %.3e", para10[6], erro10[6]);
	fprintf(f3, "\np7\t=\t%.3e +/- %.3e", para10[7], erro10[7]);
	fprintf(f3, "\np8\t=\t%.3e +/- %.3e", para10[8], erro10[8]);
	fprintf(f3, "\np9\t=\t%.3e +/- %.3e", para10[9], erro10[9]);


	/* Propagation of errors begins here */
	printf("\n\nPropagation of Errors: \n\n");
	printf("p3 --> %f \n", para10[3]);
	printf("p5 --> %f \n", para10[5]);
	printf("p7 --> %f \n", para10[7]);
	printf("p8 --> %f \n", para10[8]);
	printf("p9 --> %f \n", para10[9]);

	printf("p3u --> %f \n", cov66[3][3]);
	printf("p5u --> %f \n", cov66[5][5]);
	printf("p7u --> %f \n", cov66[7][7]);
	printf("p8u --> %f \n", cov66[8][8]);
	printf("p9u --> %f \n", cov66[9][9]);

	printf("p3p5 --> %f \n", cov66[3][5]);
	printf("p3p7 --> %f \n", cov66[3][7]);
	printf("p3p8 --> %f \n", cov66[3][8]);
	printf("p3p9 --> %f \n", cov66[3][9]);
	printf("p5p7 --> %f \n", cov66[5][7]);
	printf("p5p8 --> %f \n", cov66[5][8]);
	printf("p5p9 --> %f \n", cov66[5][9]);

	double b = TMath::Sqrt(2*TMath::Pi());
	double a1, a1u, a2, a2u, a3, a3u, a4, a4u, a5, a5u;

	a1 = para10[3]*para10[5]*b;
	a1u = TMath::Power((b*para10[5]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]), 2)*cov66[5][5] +
						2*b*b*para10[5]*para10[3]*cov66[3][5];
  a1 = a1/4.0;
	a1u = (TMath::Sqrt(a1u))/4.0;
	printf("\n\nGround State Area: %.4f +/- %.4f", a1, a1u);

	a2 = para10[3]*para10[5]*para10[7]*b;;
	a2u = TMath::Power((b*para10[5]*para10[7]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[7]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[7][7] +
					2*b*b*para10[5]*para10[3]*para10[7]*para10[7]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[7]*cov66[5][7] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[7]*cov66[3][7];
	a2 = a2/4.0;
	a2u = (TMath::Sqrt(a2u))/4.0;
	printf("\n\n0+ First Excited State Area: %.4f +/- %.4f", a2, a2u);

	a3 = para10[3]*para10[5]*para10[8]*b;;
	a3u = TMath::Power((b*para10[5]*para10[8]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[8]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[8][8] +
					2*b*b*para10[5]*para10[3]*para10[8]*para10[8]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[8]*cov66[5][8] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[8]*cov66[3][8];
	a3 = a3/4.0;
	a3u = (TMath::Sqrt(a3u))/4.0;
	printf("\n\n2+ Second Excited State Area: %.4f +/- %.4f", a3, a3u);

	a4 = para10[3]*para10[5]*para10[9]*b;;
	a4u = TMath::Power((b*para10[5]*para10[9]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[9]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[9][9] +
					2*b*b*para10[5]*para10[3]*para10[9]*para10[9]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[9]*cov66[5][9] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[9]*cov66[3][9];
	a4 = a4/4.0;
	a4u = (TMath::Sqrt(a4u))/4.0;


	printf("\n\n0+ Second Excited State Area: %.4f +/- %.4f", a4, a4u);

		a5 = para10[3]*para10[5]*para10[6]*b;;
		a5u = TMath::Power((b*para10[5]*para10[6]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[6]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[6][6] +
						2*b*b*para10[5]*para10[3]*para10[6]*para10[6]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[6]*cov66[5][6] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[6]*cov66[3][6];
		a5 = a5/4.0;
		a5u = (TMath::Sqrt(a5u))/4.0;
		printf("\n\n2+ First Excited State Area: %.4f +/- %.4f", a5, a5u);



	printf("\n\n");
	fprintf(f3, "\n\n");
	fclose(f3);
	return 0;

}


/*
 * Fit function for 6 degrees.
 *
 */
int fit6D(int angle, TFile* f2, FILE* f3, TFile* fFit){

	TH1F	*tofTop = NULL;
	TH1F	*tofMid = NULL;
	TH1F	*tofBot = NULL;
	TH1F	*tofTopOut = NULL;
	TH1F	*tofMidOut = NULL;
	TH1F	*tofBotOut = NULL;
	char 	cutinSTR[][21] = CUTIN;
	char 	cutoutSTR[][21] = CUTOUT;
	char	txtSTR[100];

	// rebinned (div by x)
	float h = 4;

	gROOT->SetBatch(kTRUE);

	fprintf(f3, "\nYou have chosen to apply a fit to the combined TOF Spectra.\n");

	/* This portion of code grabs the TOF objects that will be fitted. */
	fFit->GetObject(cutinSTR[angle],tofTop);
	fFit->GetObject(cutinSTR[angle+1],tofMid);
	fFit->GetObject(cutinSTR[angle+2],tofBot);
	fFit->GetObject(cutoutSTR[angle],tofTopOut);
	fFit->GetObject(cutoutSTR[angle+1],tofMidOut);
	fFit->GetObject(cutoutSTR[angle+2],tofBotOut);

	TH1F *mycutzu1_clone = (TH1F*)tofTop->Clone();
	mycutzu1_clone->SetName("mycutzu1_clone");
	mycutzu1_clone->SetLineColor(kMagenta);
	mycutzu1_clone->Add(tofMid);
	mycutzu1_clone->Add(tofBot);
	mycutzu1_clone->Rebin(h);
	sprintf(txtSTR, "Combined %d Degree TOF with Applied Fits", angle);
	mycutzu1_clone->SetTitle(txtSTR);

	TH1F *gasout1_1clone = (TH1F*)tofTopOut->Clone();
	gasout1_1clone->SetName("gasout1_1clone");
	gasout1_1clone->SetLineColor(kAzure);
	gasout1_1clone->Add(tofMidOut);
	gasout1_1clone->Add(tofBotOut);
	gasout1_1clone->Rebin(h);


	gasout1_1clone->Scale(3.62); //product of ratio of gas in to gas out bcis and ratio of gas in to gas out live times

	f2->cd();

	/* Fit range parameters */

	Int_t lu1 = 2364;
	Int_t hu1 = 2558;


	sprintf(txtSTR, "1xCs_degree_%d_fitted_data", angle);
	TCanvas *c1 = new TCanvas(txtSTR);
	c1->cd();

	TH1F *gasin_subtract = (TH1F*)mycutzu1_clone->Clone();
	gasin_subtract->Sumw2();
	gasout1_1clone->Sumw2();
	gasin_subtract->Add(gasout1_1clone,-1);

	/*
	 * Creation of fit parameters
	 */

	Double_t params[12] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

	TF1 *fitfunc2 = new TF1("fitfunc2", fitfuncc, lu1, hu1, 12);
	fitfunc2->SetParameter(0,(4.89));  		//constant
	fitfunc2->FixParameter(1,0);  		//linear
	fitfunc2->FixParameter(2,0);  		//quad
	fitfunc2->SetParameter(3,11.2);  	//sigma
	fitfunc2->SetParameter(4,2471.0);  //mean
	fitfunc2->SetParameter(5,163.5);   //amp1
	fitfunc2->SetParameter(6,0.1);     //amp2   2+ states
	fitfunc2->SetParameter(7,0.16);    	//amp3  50->0.1
	fitfunc2->SetParameter(8,0.001);     //amp4   2+ states
	fitfunc2->SetParameter(9,0.11);    	//amp5  50->0.1
	fitfunc2->FixParameter(10,0.0);			//amp6
	fitfunc2->FixParameter(11,0.0);			//amp7

	fitfunc2->SetParLimits(0,0,60);
	fitfunc2->SetParLimits(5,0,1000);
	fitfunc2->SetParLimits(6,0,1000);
	fitfunc2->SetParLimits(7,0,1000);
	fitfunc2->SetParLimits(8,0,1000);
	fitfunc2->SetParLimits(9,0,1000);
	//fitfunc2->SetParLimits(10,0,1000);
	//fitfunc2->SetParLimits(11,0,1000);

	TF1 *gausfuncu1 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu11 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu12 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu13 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu14 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu15 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
	TF1 *gausfuncu16 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);

	TF1 *bkdfuncu1 = new TF1("bkdfunc",background, lu1, hu1, 3);

	printf("\n\n");
	TFitResultPtr fit_result = mycutzu1_clone->Fit("fitfunc2","RS");

	fitfunc2->GetParameters(params);

	bkdfuncu1->SetParameter(0,params[0]);   //constant
	bkdfuncu1->SetParameter(1,params[1]);   //linear
	bkdfuncu1->SetParameter(2,params[2]);   //quad

	gausfuncu1->SetParameter(0,params[3]);  //sigma
	gausfuncu1->SetParameter(1,params[4]);  //mean
	gausfuncu1->SetParameter(2,params[5]);  //amp


	gausfuncu11->SetParameter(0,params[3]);  //Sigma
	gausfuncu11->SetParameter(1,params[4] - 41.4);  //amp2
	gausfuncu11->SetParameter(2,params[5]*params[6]);  //amp3  //*p5

	gausfuncu12->SetParameter(0,params[3]);  //Sigma
	gausfuncu12->SetParameter(1,params[4] - 50.5);  //Mean
	gausfuncu12->SetParameter(2,params[5]*params[7]);  //amp5

	gausfuncu13->SetParameter(0,params[3]);  //Sigma
	gausfuncu13->SetParameter(1,params[4] - 68.2);  //Mean
	gausfuncu13->SetParameter(2,params[5]*params[8]);  //amp5

	gausfuncu14->SetParameter(0,params[3]);  //Sigma
	gausfuncu14->SetParameter(1,params[4] - 96.1);  //Mean
	gausfuncu14->SetParameter(2,params[5]*params[9]);   //amp5

	gausfuncu15->SetParameter(0,params[3]);  //Sigma
	gausfuncu15->SetParameter(1,params[4] - 108);  //Mean
	gausfuncu15->SetParameter(2,params[5]*params[10]);   //amp6

	gausfuncu16->SetParameter(0,params[3]);  //Sigma
	gausfuncu16->SetParameter(1,params[4] - 126);  //Mean
	gausfuncu16->SetParameter(2,params[5]*params[11]);   //amp6


	//legend titles

	mycutzu1_clone->SetBit(TH1::kNoTitle);

	mycutzu1_clone->SetTitle("Combined Gas In");
	mycutzu1_clone->SetLineColor(kMagenta);
	mycutzu1_clone->Draw();
	gStyle->SetOptFit(1);

	gasout1_1clone->SetTitle("Combined Gas Out");
	gasout1_1clone->SetLineColor(kYellow + 1);
	gasout1_1clone->Draw("same");

	gasin_subtract->SetTitle("Gas In Minus Gas Out");
	gasin_subtract->SetLineColor(kBlue - 7);
	gasin_subtract->Draw("SAME");					// Using "SAMES" changes the info box to this fit.


	fitfunc2->SetTitle("Complete Fit Function");
	fitfunc2->SetLineColor(kRed);
	fitfunc2->Draw("same");

	bkdfuncu1->SetTitle("Background");
	bkdfuncu1->SetLineColor(kGreen - 3);
	bkdfuncu1->Draw("same");

	gausfuncu1->SetTitle("Ground State");
	gausfuncu1->SetLineColor(kAzure);
	gausfuncu1->Draw("same");

	gausfuncu11->SetTitle("2+ First Excited State 1.52 MeV");
	gausfuncu11->SetLineColor(kCyan);
	gausfuncu11->Draw("same");

	gausfuncu12->SetTitle("0+ First Excited State 1.84 MeV");
	gausfuncu12->SetLineColor(kGreen);
	gausfuncu12->Draw("same");

	gausfuncu13->SetTitle("2+ Second Excited State 2.42 MeV");
	gausfuncu13->SetLineColor(kOrange - 3);
	gausfuncu13->Draw("same");

	gausfuncu14->SetTitle("0+ Second Excited State 3.30 MeV");
	gausfuncu14->SetLineColor(kAzure+9);
	gausfuncu14->Draw("same");

	gausfuncu15->SetTitle("3.65 MeV");
	gausfuncu15->SetLineColor(kRed+4);
	gausfuncu15->Draw("same");

	gausfuncu16->SetTitle("??? MeV/126 TDC offset");
	gausfuncu16->SetLineColor(kCyan-2);
	gausfuncu16->Draw("same");

	//title: formatting needs to be fixed
		TPaveText *t = new TPaveText(0.0, 0.9, 0.3, 1.0, "brNDC"); // left-up
		t->AddText("6 Degree Fit");
		t->SetBorderSize(0);
		t->SetFillColor(gStyle->GetTitleFillColor());
		t->Draw();

	fit_result->Print("V");

	c1->BuildLegend();
	c1->Update();
	c1->Write();

	gROOT->SetBatch(kFALSE);

	double chisk = fit_result->Chi2();
	unsigned int ndf = fit_result->Ndf();
	double edm = fit_result->Edm();
	unsigned int nCall = fit_result->NCalls();
	int i,j;

	double para10[10];
	double erro10[10];
	for(i=0; i<10; i++){
		para10[i] = fit_result->Parameter(i);
		erro10[i] = fit_result->Error(i);
	}

	double cor66[10][10];
	double cov66[10][10];

	for(i=0; i<10; i++){
		for(j=0; j<10; j++){
			cor66[i][j] = fit_result->Correlation(i, j);
			cov66[i][j] = fit_result->CovMatrix(i, j);
		}
	}

	fprintf(f3, "\nChi2\t=\t%.3e", chisk);
	fprintf(f3, "\nNDf\t=\t%u", ndf);
	fprintf(f3, "\nEdm\t=\t%.3e", edm);
	fprintf(f3, "\nNCalls\t=\t%u", nCall);
	fprintf(f3, "\np0\t=\t%.3e +/- %.3e", para10[0], erro10[0]);
	fprintf(f3, "\np1\t=\t%.3e +/- %.3e", para10[1], erro10[1]);
	fprintf(f3, "\np2\t=\t%.3e +/- %.3e", para10[2], erro10[2]);
	fprintf(f3, "\np3\t=\t%.3e +/- %.3e", para10[3], erro10[3]);
	fprintf(f3, "\np4\t=\t%.3e +/- %.3e", para10[4], erro10[4]);
	fprintf(f3, "\np5\t=\t%.3e +/- %.3e", para10[5], erro10[5]);
	fprintf(f3, "\np6\t=\t%.3e +/- %.3e", para10[6], erro10[6]);
	fprintf(f3, "\np7\t=\t%.3e +/- %.3e", para10[7], erro10[7]);
	fprintf(f3, "\np8\t=\t%.3e +/- %.3e", para10[8], erro10[8]);
	fprintf(f3, "\np9\t=\t%.3e +/- %.3e", para10[9], erro10[9]);


	/* Propagation of errors begins here */
	printf("\n\nPropagation of Errors: \n\n");
	printf("p3 --> %f \n", para10[3]);
	printf("p5 --> %f \n", para10[5]);
	printf("p7 --> %f \n", para10[7]);
	printf("p8 --> %f \n", para10[8]);
	printf("p9 --> %f \n", para10[9]);

	printf("p3u --> %f \n", cov66[3][3]);
	printf("p5u --> %f \n", cov66[5][5]);
	printf("p7u --> %f \n", cov66[7][7]);
	printf("p8u --> %f \n", cov66[8][8]);
	printf("p9u --> %f \n", cov66[9][9]);

	printf("p3p5 --> %f \n", cov66[3][5]);
	printf("p3p7 --> %f \n", cov66[3][7]);
	printf("p3p8 --> %f \n", cov66[3][8]);
	printf("p3p9 --> %f \n", cov66[3][9]);
	printf("p5p7 --> %f \n", cov66[5][7]);
	printf("p5p8 --> %f \n", cov66[5][8]);
	printf("p5p9 --> %f \n", cov66[5][9]);

	double b = TMath::Sqrt(2*TMath::Pi());
	double a1, a1u, a2, a2u, a3, a3u, a4, a4u, a5, a5u;

	a1 = para10[3]*para10[5]*b;
	a1u = TMath::Power((b*para10[5]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]), 2)*cov66[5][5] +
						2*b*b*para10[5]*para10[3]*cov66[3][5];
	a1 = a1/4.0;
  a1u = (TMath::Sqrt(a1u))/4.0;
	printf("\n\nGround State Area: %.4f +/- %.4f", a1, a1u);

	a2 = para10[3]*para10[5]*para10[7]*b;;
	a2u = TMath::Power((b*para10[5]*para10[7]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[7]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[7][7] +
					2*b*b*para10[5]*para10[3]*para10[7]*para10[7]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[7]*cov66[5][7] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[7]*cov66[3][7];
	a2 = a2/4.0;
	a2u = (TMath::Sqrt(a2u))/4.0;
	printf("\n\n0+ First Excited State Area: %.4f +/- %.4f", a2, a2u);

	a3 = para10[3]*para10[5]*para10[8]*b;;
	a3u = TMath::Power((b*para10[5]*para10[8]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[8]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[8][8] +
					2*b*b*para10[5]*para10[3]*para10[8]*para10[8]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[8]*cov66[5][8] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[8]*cov66[3][8];
	a3 = a3/4.0;
	a3u = (TMath::Sqrt(a3u))/4.0;
	printf("\n\n2+ Second Excited State Area: %.4f +/- %.4f", a3, a3u);

	a4 = para10[3]*para10[5]*para10[9]*b;;
	a4u = TMath::Power((b*para10[5]*para10[9]), 2)*cov66[3][3] +
			TMath::Power((b*para10[3]*para10[9]), 2)*cov66[5][5] +
				TMath::Power((b*para10[3]*para10[5]), 2)*cov66[9][9] +
					2*b*b*para10[5]*para10[3]*para10[9]*para10[9]*cov66[3][5] +
						2*b*b*para10[3]*para10[3]*para10[5]*para10[9]*cov66[5][9] +
							2*b*b*para10[3]*para10[5]*para10[5]*para10[9]*cov66[3][9];
	a4 = a4/4.0;
	a4u = (TMath::Sqrt(a4u))/4.0;
	printf("\n\n0+ Second Excited State Area: %.4f +/- %.4f", a4, a4u);



	a5 = para10[3]*para10[5]*para10[6]*b;;
	a5u = TMath::Power((b*para10[5]*para10[6]), 2)*cov66[3][3] +
					TMath::Power((b*para10[3]*para10[6]), 2)*cov66[5][5] +
						TMath::Power((b*para10[3]*para10[5]), 2)*cov66[6][6] +
							2*b*b*para10[5]*para10[3]*para10[6]*para10[6]*cov66[3][5] +
								2*b*b*para10[3]*para10[3]*para10[5]*para10[6]*cov66[5][6] +
									2*b*b*para10[3]*para10[5]*para10[5]*para10[6]*cov66[3][6];
	a5 = a5/4.0;
	a5u = (TMath::Sqrt(a5u))/4.0;
	printf("\n\n2+ First Excited State Area: %.4f +/- %.4f", a5, a5u);



	printf("\n\n");
	fprintf(f3, "\n\n");
	fclose(f3);
	return 0;

}



int fit9D(int angle, TFile* f2, FILE* f3, TFile* fFit){

		//TH1F	*tofTop = NULL;
		TH1F	*tofMid = NULL;
		TH1F	*tofBot = NULL;
		//TH1F	*tofTopOut = NULL;
		TH1F	*tofMidOut = NULL;
		TH1F	*tofBotOut = NULL;
		char 	cutinSTR[][21] = CUTIN;
		char 	cutoutSTR[][21] = CUTOUT;
		char	txtSTR[100];

		// rebinned (div by x)
		float h = 4;

		gROOT->SetBatch(kTRUE);

		fprintf(f3, "\nYou have chosen to apply a fit to the combined TOF Spectra.\n");

		/* This portion of code grabs the TOF objects that will be fitted. */
		//fFit->GetObject(cutinSTR[angle],tofTop);
		fFit->GetObject(cutinSTR[angle+1],tofMid);
		fFit->GetObject(cutinSTR[angle+2],tofBot);
		//fFit->GetObject(cutoutSTR[angle],tofTopOut);
		fFit->GetObject(cutoutSTR[angle+1],tofMidOut);
		fFit->GetObject(cutoutSTR[angle+2],tofBotOut);

		//TH1F *mycutzu1_clone = (TH1F*)tofTop->Clone();
		TH1F *mycutzu1_clone = (TH1F*)tofMid->Clone();
		mycutzu1_clone->SetName("mycutzu1_clone");
		mycutzu1_clone->SetLineColor(kMagenta);
		//mycutzu1_clone->Add(tofMid);
		mycutzu1_clone->Add(tofBot);
		mycutzu1_clone->Rebin(h);
		sprintf(txtSTR, "Combined %d Degree TOF with Applied Fits", angle);
		mycutzu1_clone->SetTitle(txtSTR);

		//TH1F *gasout1_1clone = (TH1F*)tofTopOut->Clone();
		TH1F *gasout1_1clone = (TH1F*)tofMidOut->Clone();
		gasout1_1clone->SetName("gasout1_1clone");
		gasout1_1clone->SetLineColor(kAzure);
		//gasout1_1clone->Add(tofMidOut);
		gasout1_1clone->Add(tofBotOut);
		gasout1_1clone->Rebin(h);


		gasout1_1clone->Scale(3.62); //product of ratio of gas in to gas out bcis and ratio of gas in to gas out live times

		f2->cd();

		/* Fit range parameters */

		Int_t lu1 = 2338;
		Int_t hu1 = 2532;


		sprintf(txtSTR, "1xCs_degree_%d_fitted_data", angle);
		TCanvas *c1 = new TCanvas(txtSTR);
		c1->cd();

		TH1F *gasin_subtract = (TH1F*)mycutzu1_clone->Clone();
		gasin_subtract->Sumw2();
		gasout1_1clone->Sumw2();
		gasin_subtract->Add(gasout1_1clone,-1);

		/*
		 * Creation of fit parameters
		 */

		Double_t params[12] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

		TF1 *fitfunc2 = new TF1("fitfunc2", fitfuncc, lu1, hu1, 12);
		fitfunc2->SetParameter(0,(9.1));  		//constant
		fitfunc2->FixParameter(1,0);  		//linear
		fitfunc2->FixParameter(2,0);  		//quad
		fitfunc2->SetParameter(3,11.5);  	//sigma
		fitfunc2->SetParameter(4,2445.0);  //mean
		fitfunc2->SetParameter(5,88.0);   //amp1
		fitfunc2->SetParameter(6,0.1);     //amp2   2+ states
		fitfunc2->SetParameter(7,0.19);    	//amp3  50->0.1
		fitfunc2->SetParameter(8,0.001);     //amp4   2+ states
		fitfunc2->SetParameter(9,0.15);    	//amp5  50->0.1
		fitfunc2->FixParameter(10,0.0);			//amp6
		fitfunc2->FixParameter(11,0.0);			//amp7

		fitfunc2->SetParLimits(0,0,40);
		fitfunc2->SetParLimits(5,0,1000);
		fitfunc2->SetParLimits(6,0,1000);
		fitfunc2->SetParLimits(7,0,1000);
		fitfunc2->SetParLimits(8,0,1000);
		fitfunc2->SetParLimits(9,0,1000);
		//fitfunc2->SetParLimits(10,0,1000);
		//fitfunc2->SetParLimits(11,0,1000);

		TF1 *gausfuncu1 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu11 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu12 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu13 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu14 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu15 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu16 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);

		TF1 *bkdfuncu1 = new TF1("bkdfunc",background, lu1, hu1, 3);

		printf("\n\n");
		TFitResultPtr fit_result = mycutzu1_clone->Fit("fitfunc2","RS");

		fitfunc2->GetParameters(params);

		bkdfuncu1->SetParameter(0,params[0]);   //constant
		bkdfuncu1->SetParameter(1,params[1]);   //linear
		bkdfuncu1->SetParameter(2,params[2]);   //quad

		gausfuncu1->SetParameter(0,params[3]);  //sigma
		gausfuncu1->SetParameter(1,params[4]);  //mean
		gausfuncu1->SetParameter(2,params[5]);  //amp


		gausfuncu11->SetParameter(0,params[3]);  //Sigma
		gausfuncu11->SetParameter(1,params[4] - 41.4);  //amp2
		gausfuncu11->SetParameter(2,params[5]*params[6]);  //amp3  //*p5

		gausfuncu12->SetParameter(0,params[3]);  //Sigma
		gausfuncu12->SetParameter(1,params[4] - 50.5);  //Mean
		gausfuncu12->SetParameter(2,params[5]*params[7]);  //amp5

		gausfuncu13->SetParameter(0,params[3]);  //Sigma
		gausfuncu13->SetParameter(1,params[4] - 68.2);  //Mean
		gausfuncu13->SetParameter(2,params[5]*params[8]);  //amp5

		gausfuncu14->SetParameter(0,params[3]);  //Sigma
		gausfuncu14->SetParameter(1,params[4] - 96.1);  //Mean
		gausfuncu14->SetParameter(2,params[5]*params[9]);   //amp5

		gausfuncu15->SetParameter(0,params[3]);  //Sigma
		gausfuncu15->SetParameter(1,params[4] - 108);  //Mean
		gausfuncu15->SetParameter(2,params[5]*params[10]);   //amp6

		gausfuncu16->SetParameter(0,params[3]);  //Sigma
		gausfuncu16->SetParameter(1,params[4] - 126);  //Mean
		gausfuncu16->SetParameter(2,params[5]*params[11]);   //amp6


		//legend titles

		mycutzu1_clone->SetBit(TH1::kNoTitle);

		mycutzu1_clone->SetTitle("Combined Gas In");
		mycutzu1_clone->SetLineColor(kMagenta);
		mycutzu1_clone->Draw();
		gStyle->SetOptFit(1);

		gasout1_1clone->SetTitle("Combined Gas Out");
		gasout1_1clone->SetLineColor(kYellow + 1);
		gasout1_1clone->Draw("same");

		gasin_subtract->SetTitle("Gas In Minus Gas Out");
		gasin_subtract->SetLineColor(kBlue - 7);
		gasin_subtract->Draw("SAME");					// Using "SAMES" changes the info box to this fit.


		fitfunc2->SetTitle("Complete Fit Function");
		fitfunc2->SetLineColor(kRed);
		fitfunc2->Draw("same");

		bkdfuncu1->SetTitle("Background");
		bkdfuncu1->SetLineColor(kGreen - 3);
		bkdfuncu1->Draw("same");

		gausfuncu1->SetTitle("Ground State");
		gausfuncu1->SetLineColor(kAzure);
		gausfuncu1->Draw("same");

		gausfuncu11->SetTitle("2+ First Excited State 1.52 MeV");
		gausfuncu11->SetLineColor(kCyan);
		gausfuncu11->Draw("same");

		gausfuncu12->SetTitle("0+ First Excited State 1.84 MeV");
		gausfuncu12->SetLineColor(kGreen);
		gausfuncu12->Draw("same");

		gausfuncu13->SetTitle("2+ Second Excited State 2.42 MeV");
		gausfuncu13->SetLineColor(kOrange - 3);
		gausfuncu13->Draw("same");

		gausfuncu14->SetTitle("0+ Second Excited State 3.30 MeV");
		gausfuncu14->SetLineColor(kAzure+9);
		gausfuncu14->Draw("same");

		gausfuncu15->SetTitle("3.65 MeV");
		gausfuncu15->SetLineColor(kRed+4);
		gausfuncu15->Draw("same");

		gausfuncu16->SetTitle("??? MeV/126 TDC offset");
		gausfuncu16->SetLineColor(kCyan-2);
		gausfuncu16->Draw("same");

		//title: formatting needs to be fixed
		TPaveText *t = new TPaveText(0.0, 0.9, 0.3, 1.0, "brNDC"); // left-up
		t->AddText("9 Degree Fit");
		t->SetBorderSize(0);
		t->SetFillColor(gStyle->GetTitleFillColor());
		t->Draw();

		fit_result->Print("V");

		c1->BuildLegend();
		c1->Update();
		c1->Write();

		gROOT->SetBatch(kFALSE);

		double chisk = fit_result->Chi2();
		unsigned int ndf = fit_result->Ndf();
		double edm = fit_result->Edm();
		unsigned int nCall = fit_result->NCalls();
		int i,j;

		double para10[10];
		double erro10[10];
		for(i=0; i<10; i++){
			para10[i] = fit_result->Parameter(i);
			erro10[i] = fit_result->Error(i);
		}

		double cor66[10][10];
		double cov66[10][10];

		for(i=0; i<10; i++){
			for(j=0; j<10; j++){
				cor66[i][j] = fit_result->Correlation(i, j);
				cov66[i][j] = fit_result->CovMatrix(i, j);
			}
		}

		fprintf(f3, "\nChi2\t=\t%.3e", chisk);
		fprintf(f3, "\nNDf\t=\t%u", ndf);
		fprintf(f3, "\nEdm\t=\t%.3e", edm);
		fprintf(f3, "\nNCalls\t=\t%u", nCall);
		fprintf(f3, "\np0\t=\t%.3e +/- %.3e", para10[0], erro10[0]);
		fprintf(f3, "\np1\t=\t%.3e +/- %.3e", para10[1], erro10[1]);
		fprintf(f3, "\np2\t=\t%.3e +/- %.3e", para10[2], erro10[2]);
		fprintf(f3, "\np3\t=\t%.3e +/- %.3e", para10[3], erro10[3]);
		fprintf(f3, "\np4\t=\t%.3e +/- %.3e", para10[4], erro10[4]);
		fprintf(f3, "\np5\t=\t%.3e +/- %.3e", para10[5], erro10[5]);
		fprintf(f3, "\np6\t=\t%.3e +/- %.3e", para10[6], erro10[6]);
		fprintf(f3, "\np7\t=\t%.3e +/- %.3e", para10[7], erro10[7]);
		fprintf(f3, "\np8\t=\t%.3e +/- %.3e", para10[8], erro10[8]);
		fprintf(f3, "\np9\t=\t%.3e +/- %.3e", para10[9], erro10[9]);


		/* Propagation of errors begins here */
		printf("\n\nPropagation of Errors: \n\n");
		printf("p3 --> %f \n", para10[3]);
		printf("p5 --> %f \n", para10[5]);
		printf("p7 --> %f \n", para10[7]);
		printf("p8 --> %f \n", para10[8]);
		printf("p9 --> %f \n", para10[9]);

		printf("p3u --> %f \n", cov66[3][3]);
		printf("p5u --> %f \n", cov66[5][5]);
		printf("p7u --> %f \n", cov66[7][7]);
		printf("p8u --> %f \n", cov66[8][8]);
		printf("p9u --> %f \n", cov66[9][9]);

		printf("p3p5 --> %f \n", cov66[3][5]);
		printf("p3p7 --> %f \n", cov66[3][7]);
		printf("p3p8 --> %f \n", cov66[3][8]);
		printf("p3p9 --> %f \n", cov66[3][9]);
		printf("p5p7 --> %f \n", cov66[5][7]);
		printf("p5p8 --> %f \n", cov66[5][8]);
		printf("p5p9 --> %f \n", cov66[5][9]);

		double b = TMath::Sqrt(2*TMath::Pi());
		double a1, a1u, a2, a2u, a3, a3u, a4, a4u, a5, a5u;

		a1 = para10[3]*para10[5]*b;
		a1u = TMath::Power((b*para10[5]), 2)*cov66[3][3] +
					TMath::Power((b*para10[3]), 2)*cov66[5][5] +
							2*b*b*para10[5]*para10[3]*cov66[3][5];
		a1 = a1/4.0;
		a1u = (TMath::Sqrt(a1u))/4.0;
		printf("\n\nGround State Area: %.4f +/- %.4f", a1, a1u);

		a2 = para10[3]*para10[5]*para10[7]*b;;
		a2u = TMath::Power((b*para10[5]*para10[7]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[7]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[7][7] +
						2*b*b*para10[5]*para10[3]*para10[7]*para10[7]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[7]*cov66[5][7] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[7]*cov66[3][7];
		a2 = a2/4.0;
		a2u = (TMath::Sqrt(a2u))/4.0;
		printf("\n\n0+ First Excited State Area: %.4f +/- %.4f", a2, a2u);

		a3 = para10[3]*para10[5]*para10[8]*b;;
		a3u = TMath::Power((b*para10[5]*para10[8]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[8]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[8][8] +
						2*b*b*para10[5]*para10[3]*para10[8]*para10[8]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[8]*cov66[5][8] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[8]*cov66[3][8];
		a3 = a3/4.0;
		a3u = (TMath::Sqrt(a3u))/4.0;
		printf("\n\n2+ Second Excited State Area: %.4f +/- %.4f", a3, a3u);

		a4 = para10[3]*para10[5]*para10[9]*b;;
		a4u = TMath::Power((b*para10[5]*para10[9]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[9]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[9][9] +
						2*b*b*para10[5]*para10[3]*para10[9]*para10[9]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[9]*cov66[5][9] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[9]*cov66[3][9];
		a4 = a4/4.0;
		a4u = (TMath::Sqrt(a4u))/4.0;
		printf("\n\n0+ Second Excited State Area: %.4f +/- %.4f", a4, a4u);

		printf("\n\n0+ Second Excited State Area: %.4f +/- %.4f", a4, a4u);

		a5 = para10[3]*para10[5]*para10[6]*b;;
		a5u = TMath::Power((b*para10[5]*para10[6]), 2)*cov66[3][3] +
						TMath::Power((b*para10[3]*para10[6]), 2)*cov66[5][5] +
							TMath::Power((b*para10[3]*para10[5]), 2)*cov66[6][6] +
								2*b*b*para10[5]*para10[3]*para10[6]*para10[6]*cov66[3][5] +
									2*b*b*para10[3]*para10[3]*para10[5]*para10[6]*cov66[5][6] +
										2*b*b*para10[3]*para10[5]*para10[5]*para10[6]*cov66[3][6];
		a5 = a5/4.0;
		a5u = (TMath::Sqrt(a5u))/4.0;
		printf("\n\n2+ First Excited State Area: %.4f +/- %.4f", a5, a5u);


		printf("\n\n");
		fprintf(f3, "\n\n");
		fclose(f3);
	return 0;
}



int fit12D(int angle, TFile* f2, FILE* f3, TFile* fFit){

		TH1F	*tofTop = NULL;
		TH1F	*tofMid = NULL;
		TH1F	*tofBot = NULL;
		TH1F	*tofTopOut = NULL;
		TH1F	*tofMidOut = NULL;
		TH1F	*tofBotOut = NULL;
		char 	cutinSTR[][21] = CUTIN;
		char 	cutoutSTR[][21] = CUTOUT;
		char	txtSTR[100];

		// rebinned (div by x)
		float h = 4;

		gROOT->SetBatch(kTRUE);

		fprintf(f3, "\nYou have chosen to apply a fit to the combined TOF Spectra.\n");

		/* This portion of code grabs the TOF objects that will be fitted. */
		fFit->GetObject(cutinSTR[angle],tofTop);
		fFit->GetObject(cutinSTR[angle+1],tofMid);
		fFit->GetObject(cutinSTR[angle+2],tofBot);
		fFit->GetObject(cutoutSTR[angle],tofTopOut);
		fFit->GetObject(cutoutSTR[angle+1],tofMidOut);
		fFit->GetObject(cutoutSTR[angle+2],tofBotOut);

		TH1F *mycutzu1_clone = (TH1F*)tofTop->Clone();
		mycutzu1_clone->SetName("mycutzu1_clone");
		mycutzu1_clone->SetLineColor(kMagenta);
		mycutzu1_clone->Add(tofMid);
		mycutzu1_clone->Add(tofBot);
		mycutzu1_clone->Rebin(h);
		sprintf(txtSTR, "Combined %d Degree TOF with Applied Fits", angle);
		mycutzu1_clone->SetTitle(txtSTR);

		TH1F *gasout1_1clone = (TH1F*)tofTopOut->Clone();
		gasout1_1clone->SetName("gasout1_1clone");
		gasout1_1clone->SetLineColor(kAzure);
		gasout1_1clone->Add(tofMidOut);
		gasout1_1clone->Add(tofBotOut);
		gasout1_1clone->Rebin(h);


		gasout1_1clone->Scale(3.62); //product of ratio of gas in to gas out bcis and ratio of gas in to gas out live times

		f2->cd();

		/* Fit range parameters */

		Int_t lu1 = 2399;
		Int_t hu1 = 2594;


		sprintf(txtSTR, "1xCs_degree_%d_fitted_data", angle);
		TCanvas *c1 = new TCanvas(txtSTR);
		c1->cd();

		TH1F *gasin_subtract = (TH1F*)mycutzu1_clone->Clone();
		gasin_subtract->Sumw2();
		gasout1_1clone->Sumw2();
		gasin_subtract->Add(gasout1_1clone,-1);

		/*
		 * Creation of fit parameters
		 */

		Double_t params[12] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

		TF1 *fitfunc2 = new TF1("fitfunc2", fitfuncc, lu1, hu1, 12);
		fitfunc2->SetParameter(0,(13.3));  		//constant
		fitfunc2->FixParameter(1,0);  		//linear
		fitfunc2->FixParameter(2,0);  		//quad
		fitfunc2->SetParameter(3,9.8);  	//sigma
		fitfunc2->SetParameter(4,2507.0);  //mean
		fitfunc2->SetParameter(5,96.2);   //amp1
		fitfunc2->SetParameter(6,0.1);     //amp2   2+ states
		fitfunc2->SetParameter(7,0.001);    	//amp3  50->0.1
		fitfunc2->SetParameter(8,0.001);     //amp4   2+ states
		fitfunc2->SetParameter(9,0.05);    	//amp5  50->0.1
		fitfunc2->FixParameter(10,0.0);			//amp6
		fitfunc2->FixParameter(11,0.0);			//amp7

		fitfunc2->SetParLimits(0,0,40);
		fitfunc2->SetParLimits(5,0,1000);
		fitfunc2->SetParLimits(6,0,1000);
		fitfunc2->SetParLimits(7,0,1000);
		fitfunc2->SetParLimits(8,0,1000);
		fitfunc2->SetParLimits(9,0,1000);
		//fitfunc2->SetParLimits(10,0,1000);
		//fitfunc2->SetParLimits(11,0,1000);

		TF1 *gausfuncu1 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu11 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu12 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu13 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu14 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu15 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu16 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);

		TF1 *bkdfuncu1 = new TF1("bkdfunc",background, lu1, hu1, 3);

		printf("\n\n");
		TFitResultPtr fit_result = mycutzu1_clone->Fit("fitfunc2","RS");

		fitfunc2->GetParameters(params);

		bkdfuncu1->SetParameter(0,params[0]);   //constant
		bkdfuncu1->SetParameter(1,params[1]);   //linear
		bkdfuncu1->SetParameter(2,params[2]);   //quad

		gausfuncu1->SetParameter(0,params[3]);  //sigma
		gausfuncu1->SetParameter(1,params[4]);  //mean
		gausfuncu1->SetParameter(2,params[5]);  //amp


		gausfuncu11->SetParameter(0,params[3]);  //Sigma
		gausfuncu11->SetParameter(1,params[4] - 41.4);  //amp2
		gausfuncu11->SetParameter(2,params[5]*params[6]);  //amp3  //*p5

		gausfuncu12->SetParameter(0,params[3]);  //Sigma
		gausfuncu12->SetParameter(1,params[4] - 50.5);  //Mean
		gausfuncu12->SetParameter(2,params[5]*params[7]);  //amp5

		gausfuncu13->SetParameter(0,params[3]);  //Sigma
		gausfuncu13->SetParameter(1,params[4] - 68.2);  //Mean
		gausfuncu13->SetParameter(2,params[5]*params[8]);  //amp5

		gausfuncu14->SetParameter(0,params[3]);  //Sigma
		gausfuncu14->SetParameter(1,params[4] - 96.1);  //Mean
		gausfuncu14->SetParameter(2,params[5]*params[9]);   //amp5

		gausfuncu15->SetParameter(0,params[3]);  //Sigma
		gausfuncu15->SetParameter(1,params[4] - 108);  //Mean
		gausfuncu15->SetParameter(2,params[5]*params[10]);   //amp6

		gausfuncu16->SetParameter(0,params[3]);  //Sigma
		gausfuncu16->SetParameter(1,params[4] - 126);  //Mean
		gausfuncu16->SetParameter(2,params[5]*params[11]);   //amp6


		//legend titles

		mycutzu1_clone->SetBit(TH1::kNoTitle);

		mycutzu1_clone->SetTitle("Combined Gas In");
		mycutzu1_clone->SetLineColor(kMagenta);
		mycutzu1_clone->Draw();
		gStyle->SetOptFit(1);

		gasout1_1clone->SetTitle("Combined Gas Out");
		gasout1_1clone->SetLineColor(kYellow + 1);
		gasout1_1clone->Draw("same");

		gasin_subtract->SetTitle("Gas In Minus Gas Out");
		gasin_subtract->SetLineColor(kBlue - 7);
		gasin_subtract->Draw("SAME");					// Using "SAMES" changes the info box to this fit.


		fitfunc2->SetTitle("Complete Fit Function");
		fitfunc2->SetLineColor(kRed);
		fitfunc2->Draw("same");

		bkdfuncu1->SetTitle("Background");
		bkdfuncu1->SetLineColor(kGreen - 3);
		bkdfuncu1->Draw("same");

		gausfuncu1->SetTitle("Ground State");
		gausfuncu1->SetLineColor(kAzure);
		gausfuncu1->Draw("same");

		gausfuncu11->SetTitle("2+ First Excited State 1.52 MeV");
		gausfuncu11->SetLineColor(kCyan);
		gausfuncu11->Draw("same");

		gausfuncu12->SetTitle("0+ First Excited State 1.84 MeV");
		gausfuncu12->SetLineColor(kGreen);
		gausfuncu12->Draw("same");

		gausfuncu13->SetTitle("2+ Second Excited State 2.42 MeV");
		gausfuncu13->SetLineColor(kOrange - 3);
		gausfuncu13->Draw("same");

		gausfuncu14->SetTitle("0+ Second Excited State 3.30 MeV");
		gausfuncu14->SetLineColor(kAzure+9);
		gausfuncu14->Draw("same");

		gausfuncu15->SetTitle("3.65 MeV");
		gausfuncu15->SetLineColor(kRed+4);
		gausfuncu15->Draw("same");

		gausfuncu16->SetTitle("??? MeV/126 TDC offset");
		gausfuncu16->SetLineColor(kCyan-2);
		gausfuncu16->Draw("same");

		//title: formatting needs to be fixed
		TPaveText *t = new TPaveText(0.0, 0.9, 0.3, 1.0, "brNDC"); // left-up
		t->AddText("12 Degree Fit");
		t->SetBorderSize(0);
		t->SetFillColor(gStyle->GetTitleFillColor());
		t->Draw();

		fit_result->Print("V");

		c1->BuildLegend();
		c1->Update();
		c1->Write();

		gROOT->SetBatch(kFALSE);

		double chisk = fit_result->Chi2();
		unsigned int ndf = fit_result->Ndf();
		double edm = fit_result->Edm();
		unsigned int nCall = fit_result->NCalls();
		int i,j;

		double para10[10];
		double erro10[10];
		for(i=0; i<10; i++){
			para10[i] = fit_result->Parameter(i);
			erro10[i] = fit_result->Error(i);
		}

		double cor66[10][10];
		double cov66[10][10];

		for(i=0; i<10; i++){
			for(j=0; j<10; j++){
				cor66[i][j] = fit_result->Correlation(i, j);
				cov66[i][j] = fit_result->CovMatrix(i, j);
			}
		}

		fprintf(f3, "\nChi2\t=\t%.3e", chisk);
		fprintf(f3, "\nNDf\t=\t%u", ndf);
		fprintf(f3, "\nEdm\t=\t%.3e", edm);
		fprintf(f3, "\nNCalls\t=\t%u", nCall);
		fprintf(f3, "\np0\t=\t%.3e +/- %.3e", para10[0], erro10[0]);
		fprintf(f3, "\np1\t=\t%.3e +/- %.3e", para10[1], erro10[1]);
		fprintf(f3, "\np2\t=\t%.3e +/- %.3e", para10[2], erro10[2]);
		fprintf(f3, "\np3\t=\t%.3e +/- %.3e", para10[3], erro10[3]);
		fprintf(f3, "\np4\t=\t%.3e +/- %.3e", para10[4], erro10[4]);
		fprintf(f3, "\np5\t=\t%.3e +/- %.3e", para10[5], erro10[5]);
		fprintf(f3, "\np6\t=\t%.3e +/- %.3e", para10[6], erro10[6]);
		fprintf(f3, "\np7\t=\t%.3e +/- %.3e", para10[7], erro10[7]);
		fprintf(f3, "\np8\t=\t%.3e +/- %.3e", para10[8], erro10[8]);
		fprintf(f3, "\np9\t=\t%.3e +/- %.3e", para10[9], erro10[9]);


		/* Propagation of errors begins here */
		printf("\n\nPropagation of Errors: \n\n");
		printf("p3 --> %f \n", para10[3]);
		printf("p5 --> %f \n", para10[5]);
		printf("p7 --> %f \n", para10[7]);
		printf("p8 --> %f \n", para10[8]);
		printf("p9 --> %f \n", para10[9]);

		printf("p3u --> %f \n", cov66[3][3]);
		printf("p5u --> %f \n", cov66[5][5]);
		printf("p7u --> %f \n", cov66[7][7]);
		printf("p8u --> %f \n", cov66[8][8]);
		printf("p9u --> %f \n", cov66[9][9]);

		printf("p3p5 --> %f \n", cov66[3][5]);
		printf("p3p7 --> %f \n", cov66[3][7]);
		printf("p3p8 --> %f \n", cov66[3][8]);
		printf("p3p9 --> %f \n", cov66[3][9]);
		printf("p5p7 --> %f \n", cov66[5][7]);
		printf("p5p8 --> %f \n", cov66[5][8]);
		printf("p5p9 --> %f \n", cov66[5][9]);

		double b = TMath::Sqrt(2*TMath::Pi());
		double a1, a1u, a2, a2u, a3, a3u, a4, a4u, a5, a5u;

		a1 = para10[3]*para10[5]*b;
		a1u = TMath::Power((b*para10[5]), 2)*cov66[3][3] +
					TMath::Power((b*para10[3]), 2)*cov66[5][5] +
							2*b*b*para10[5]*para10[3]*cov66[3][5];
		a1 = a1/4.0;
		a1u = (TMath::Sqrt(a1u))/4.0;
		printf("\n\nGround State Area: %.4f +/- %.4f", a1, a1u);

		a2 = para10[3]*para10[5]*para10[7]*b;;
		a2u = TMath::Power((b*para10[5]*para10[7]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[7]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[7][7] +
						2*b*b*para10[5]*para10[3]*para10[7]*para10[7]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[7]*cov66[5][7] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[7]*cov66[3][7];
		a2 = a2/4.0;
		a2u = (TMath::Sqrt(a2u))/4.0;
		printf("\n\n0+ First Excited State Area: %.4f +/- %.4f", a2, a2u);

		a3 = para10[3]*para10[5]*para10[8]*b;;
		a3u = TMath::Power((b*para10[5]*para10[8]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[8]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[8][8] +
						2*b*b*para10[5]*para10[3]*para10[8]*para10[8]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[8]*cov66[5][8] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[8]*cov66[3][8];
		a3 = a3/4.0;
		a3u = (TMath::Sqrt(a3u))/4.0;
		printf("\n\n2+ Second Excited State Area: %.4f +/- %.4f", a3, a3u);

		a4 = para10[3]*para10[5]*para10[9]*b;;
		a4u = TMath::Power((b*para10[5]*para10[9]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[9]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[9][9] +
						2*b*b*para10[5]*para10[3]*para10[9]*para10[9]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[9]*cov66[5][9] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[9]*cov66[3][9];
		a4 = a4/4.0;
		a4u = (TMath::Sqrt(a4u))/4.0;
		printf("\n\n0+ Second Excited State Area: %.4f +/- %.4f", a4, a4u);


		a5 = para10[3]*para10[5]*para10[6]*b;;
		a5u = TMath::Power((b*para10[5]*para10[6]), 2)*cov66[3][3] +
						TMath::Power((b*para10[3]*para10[6]), 2)*cov66[5][5] +
							TMath::Power((b*para10[3]*para10[5]), 2)*cov66[6][6] +
								2*b*b*para10[5]*para10[3]*para10[6]*para10[6]*cov66[3][5] +
									2*b*b*para10[3]*para10[3]*para10[5]*para10[6]*cov66[5][6] +
										2*b*b*para10[3]*para10[5]*para10[5]*para10[6]*cov66[3][6];
		a5 = a5/4.0;
		a5u = (TMath::Sqrt(a5u))/4.0;
		printf("\n\n2+ First Excited State Area: %.4f +/- %.4f", a5, a5u);



		printf("\n\n");
		fprintf(f3, "\n\n");
		fclose(f3);
	return 0;
}



int fit15D(int angle, TFile* f2, FILE* f3, TFile* fFit){

		TH1F	*tofTop = NULL;
		TH1F	*tofMid = NULL;
		TH1F	*tofBot = NULL;
		TH1F	*tofTopOut = NULL;
		TH1F	*tofMidOut = NULL;
		TH1F	*tofBotOut = NULL;
		char 	cutinSTR[][21] = CUTIN;
		char 	cutoutSTR[][21] = CUTOUT;
		char	txtSTR[100];

		// rebinned (div by x)
		float h = 4;

		gROOT->SetBatch(kTRUE);

		fprintf(f3, "\nYou have chosen to apply a fit to the combined TOF Spectra.\n");

		/* This portion of code grabs the TOF objects that will be fitted. */
		fFit->GetObject(cutinSTR[angle],tofTop);
		fFit->GetObject(cutinSTR[angle+1],tofMid);
		fFit->GetObject(cutinSTR[angle+2],tofBot);
		fFit->GetObject(cutoutSTR[angle],tofTopOut);
		fFit->GetObject(cutoutSTR[angle+1],tofMidOut);
		fFit->GetObject(cutoutSTR[angle+2],tofBotOut);

		TH1F *mycutzu1_clone = (TH1F*)tofTop->Clone();
		mycutzu1_clone->SetName("mycutzu1_clone");
		mycutzu1_clone->SetLineColor(kMagenta);
		mycutzu1_clone->Add(tofMid);
		mycutzu1_clone->Add(tofBot);
		mycutzu1_clone->Rebin(h);
		sprintf(txtSTR, "Combined %d Degree TOF with Applied Fits", angle);
		mycutzu1_clone->SetTitle(txtSTR);

		TH1F *gasout1_1clone = (TH1F*)tofTopOut->Clone();
		gasout1_1clone->SetName("gasout1_1clone");
		gasout1_1clone->SetLineColor(kAzure);
		gasout1_1clone->Add(tofMidOut);
		gasout1_1clone->Add(tofBotOut);
		gasout1_1clone->Rebin(h);


		gasout1_1clone->Scale(3.62); //product of ratio of gas in to gas out bcis and ratio of gas in to gas out live times

		f2->cd();

		/* Fit range parameters */

		Int_t lu1 = 2382;
		Int_t hu1 = 2600;


		sprintf(txtSTR, "1xCs_degree_%d_fitted_data", angle);
		TCanvas *c1 = new TCanvas(txtSTR);
		c1->cd();

		TH1F *gasin_subtract = (TH1F*)mycutzu1_clone->Clone();
		gasin_subtract->Sumw2();
		gasout1_1clone->Sumw2();
		gasin_subtract->Add(gasout1_1clone,-1);

		/*
		 * Creation of fit parameters
		 */

		Double_t params[12] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

		TF1 *fitfunc2 = new TF1("fitfunc2", fitfuncc, lu1, hu1, 12);
		fitfunc2->SetParameter(0,(8.1));  		//constant
		fitfunc2->FixParameter(1,0);  		//linear
		fitfunc2->FixParameter(2,0);  		//quad
		fitfunc2->SetParameter(3,13.8);  	//sigma
		fitfunc2->SetParameter(4,2491.0);  //mean
		fitfunc2->SetParameter(5,43.80);   //amp1
		fitfunc2->SetParameter(6,0.1);     //amp2   2+ states
		fitfunc2->SetParameter(7,0.35);    	//amp3  50->0.1
		fitfunc2->SetParameter(8,0.001);     //amp4   2+ states
		fitfunc2->SetParameter(9,0.3);    	//amp5  50->0.1
		fitfunc2->FixParameter(10,0.0);			//amp6
		fitfunc2->FixParameter(11,0.0);			//amp7

		fitfunc2->SetParLimits(0,0,60);
		fitfunc2->SetParLimits(5,0,1000);
		fitfunc2->SetParLimits(6,0,1000);
		fitfunc2->SetParLimits(7,0,1000);
		fitfunc2->SetParLimits(8,0,1000);
		fitfunc2->SetParLimits(9,0,1000);
		//fitfunc2->SetParLimits(10,0,1000);
		//fitfunc2->SetParLimits(11,0,1000);

		TF1 *gausfuncu1 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu11 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu12 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu13 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu14 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu15 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu16 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);

		TF1 *bkdfuncu1 = new TF1("bkdfunc",background, lu1, hu1, 3);

		printf("\n\n");
		TFitResultPtr fit_result = mycutzu1_clone->Fit("fitfunc2","RS");

		fitfunc2->GetParameters(params);

		bkdfuncu1->SetParameter(0,params[0]);   //constant
		bkdfuncu1->SetParameter(1,params[1]);   //linear
		bkdfuncu1->SetParameter(2,params[2]);   //quad

		gausfuncu1->SetParameter(0,params[3]);  //sigma
		gausfuncu1->SetParameter(1,params[4]);  //mean
		gausfuncu1->SetParameter(2,params[5]);  //amp


		gausfuncu11->SetParameter(0,params[3]);  //Sigma
		gausfuncu11->SetParameter(1,params[4] - 41.4);  //amp2
		gausfuncu11->SetParameter(2,params[5]*params[6]);  //amp3  //*p5

		gausfuncu12->SetParameter(0,params[3]);  //Sigma
		gausfuncu12->SetParameter(1,params[4] - 50.5);  //Mean
		gausfuncu12->SetParameter(2,params[5]*params[7]);  //amp5

		gausfuncu13->SetParameter(0,params[3]);  //Sigma
		gausfuncu13->SetParameter(1,params[4] - 68.2);  //Mean
		gausfuncu13->SetParameter(2,params[5]*params[8]);  //amp5

		gausfuncu14->SetParameter(0,params[3]);  //Sigma
		gausfuncu14->SetParameter(1,params[4] - 96.1);  //Mean
		gausfuncu14->SetParameter(2,params[5]*params[9]);   //amp5

		gausfuncu15->SetParameter(0,params[3]);  //Sigma
		gausfuncu15->SetParameter(1,params[4] - 108);  //Mean
		gausfuncu15->SetParameter(2,params[5]*params[10]);   //amp6

		gausfuncu16->SetParameter(0,params[3]);  //Sigma
		gausfuncu16->SetParameter(1,params[4] - 126);  //Mean
		gausfuncu16->SetParameter(2,params[5]*params[11]);   //amp6


		//legend titles

		mycutzu1_clone->SetBit(TH1::kNoTitle);

		mycutzu1_clone->SetTitle("Combined Gas In");
		mycutzu1_clone->SetLineColor(kMagenta);
		mycutzu1_clone->Draw();
		gStyle->SetOptFit(1);

		gasout1_1clone->SetTitle("Combined Gas Out");
		gasout1_1clone->SetLineColor(kYellow + 1);
		gasout1_1clone->Draw("same");

		gasin_subtract->SetTitle("Gas In Minus Gas Out");
		gasin_subtract->SetLineColor(kBlue - 7);
		gasin_subtract->Draw("SAME");					// Using "SAMES" changes the info box to this fit.


		fitfunc2->SetTitle("Complete Fit Function");
		fitfunc2->SetLineColor(kRed);
		fitfunc2->Draw("same");

		bkdfuncu1->SetTitle("Background");
		bkdfuncu1->SetLineColor(kGreen - 3);
		bkdfuncu1->Draw("same");

		gausfuncu1->SetTitle("Ground State");
		gausfuncu1->SetLineColor(kAzure);
		gausfuncu1->Draw("same");

		gausfuncu11->SetTitle("2+ First Excited State 1.52 MeV");
		gausfuncu11->SetLineColor(kCyan);
		gausfuncu11->Draw("same");

		gausfuncu12->SetTitle("0+ First Excited State 1.84 MeV");
		gausfuncu12->SetLineColor(kGreen);
		gausfuncu12->Draw("same");

		gausfuncu13->SetTitle("2+ Second Excited State 2.42 MeV");
		gausfuncu13->SetLineColor(kOrange - 3);
		gausfuncu13->Draw("same");

		gausfuncu14->SetTitle("0+ Second Excited State 3.30 MeV");
		gausfuncu14->SetLineColor(kAzure+9);
		gausfuncu14->Draw("same");

		gausfuncu15->SetTitle("3.65 MeV");
		gausfuncu15->SetLineColor(kRed+4);
		gausfuncu15->Draw("same");

		gausfuncu16->SetTitle("??? MeV/126 TDC offset");
		gausfuncu16->SetLineColor(kCyan-2);
		gausfuncu16->Draw("same");

		//title: formatting needs to be fixed
		TPaveText *t = new TPaveText(0.0, 0.9, 0.3, 1.0, "brNDC"); // left-up
		t->AddText("15 Degree Fit");
		t->SetBorderSize(0);
		t->SetFillColor(gStyle->GetTitleFillColor());
		t->Draw();

		fit_result->Print("V");

		c1->BuildLegend();
		c1->Update();
		c1->Write();

		gROOT->SetBatch(kFALSE);

		double chisk = fit_result->Chi2();
		unsigned int ndf = fit_result->Ndf();
		double edm = fit_result->Edm();
		unsigned int nCall = fit_result->NCalls();
		int i,j;

		double para10[10];
		double erro10[10];
		for(i=0; i<10; i++){
			para10[i] = fit_result->Parameter(i);
			erro10[i] = fit_result->Error(i);
		}

		double cor66[10][10];
		double cov66[10][10];

		for(i=0; i<10; i++){
			for(j=0; j<10; j++){
				cor66[i][j] = fit_result->Correlation(i, j);
				cov66[i][j] = fit_result->CovMatrix(i, j);
			}
		}

		fprintf(f3, "\nChi2\t=\t%.3e", chisk);
		fprintf(f3, "\nNDf\t=\t%u", ndf);
		fprintf(f3, "\nEdm\t=\t%.3e", edm);
		fprintf(f3, "\nNCalls\t=\t%u", nCall);
		fprintf(f3, "\np0\t=\t%.3e +/- %.3e", para10[0], erro10[0]);
		fprintf(f3, "\np1\t=\t%.3e +/- %.3e", para10[1], erro10[1]);
		fprintf(f3, "\np2\t=\t%.3e +/- %.3e", para10[2], erro10[2]);
		fprintf(f3, "\np3\t=\t%.3e +/- %.3e", para10[3], erro10[3]);
		fprintf(f3, "\np4\t=\t%.3e +/- %.3e", para10[4], erro10[4]);
		fprintf(f3, "\np5\t=\t%.3e +/- %.3e", para10[5], erro10[5]);
		fprintf(f3, "\np6\t=\t%.3e +/- %.3e", para10[6], erro10[6]);
		fprintf(f3, "\np7\t=\t%.3e +/- %.3e", para10[7], erro10[7]);
		fprintf(f3, "\np8\t=\t%.3e +/- %.3e", para10[8], erro10[8]);
		fprintf(f3, "\np9\t=\t%.3e +/- %.3e", para10[9], erro10[9]);


		/* Propagation of errors begins here */
		printf("\n\nPropagation of Errors: \n\n");
		printf("p3 --> %f \n", para10[3]);
		printf("p5 --> %f \n", para10[5]);
		printf("p7 --> %f \n", para10[7]);
		printf("p8 --> %f \n", para10[8]);
		printf("p9 --> %f \n", para10[9]);

		printf("p3u --> %f \n", cov66[3][3]);
		printf("p5u --> %f \n", cov66[5][5]);
		printf("p7u --> %f \n", cov66[7][7]);
		printf("p8u --> %f \n", cov66[8][8]);
		printf("p9u --> %f \n", cov66[9][9]);

		printf("p3p5 --> %f \n", cov66[3][5]);
		printf("p3p7 --> %f \n", cov66[3][7]);
		printf("p3p8 --> %f \n", cov66[3][8]);
		printf("p3p9 --> %f \n", cov66[3][9]);
		printf("p5p7 --> %f \n", cov66[5][7]);
		printf("p5p8 --> %f \n", cov66[5][8]);
		printf("p5p9 --> %f \n", cov66[5][9]);

		double b = TMath::Sqrt(2*TMath::Pi());
		double a1, a1u, a2, a2u, a3, a3u, a4, a4u, a5, a5u;

		a1 = para10[3]*para10[5]*b;
		a1u = TMath::Power((b*para10[5]), 2)*cov66[3][3] +
					TMath::Power((b*para10[3]), 2)*cov66[5][5] +
							2*b*b*para10[5]*para10[3]*cov66[3][5];
		a1 = a1/4.0;
		a1u = (TMath::Sqrt(a1u))/4.0;
		printf("\n\nGround State Area: %.4f +/- %.4f", a1, a1u);

		a2 = para10[3]*para10[5]*para10[7]*b;;
		a2u = TMath::Power((b*para10[5]*para10[7]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[7]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[7][7] +
						2*b*b*para10[5]*para10[3]*para10[7]*para10[7]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[7]*cov66[5][7] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[7]*cov66[3][7];
		a2 = a2/4.0;
		a2u = (TMath::Sqrt(a2u))/4.0;
		printf("\n\n0+ First Excited State Area: %.4f +/- %.4f", a2, a2u);

		a3 = para10[3]*para10[5]*para10[8]*b;;
		a3u = TMath::Power((b*para10[5]*para10[8]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[8]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[8][8] +
						2*b*b*para10[5]*para10[3]*para10[8]*para10[8]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[8]*cov66[5][8] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[8]*cov66[3][8];
		a3 = a3/4.0;
		a3u = (TMath::Sqrt(a3u))/4.0;
		printf("\n\n2+ Second Excited State Area: %.4f +/- %.4f", a3, a3u);

		a4 = para10[3]*para10[5]*para10[9]*b;;
		a4u = TMath::Power((b*para10[5]*para10[9]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[9]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[9][9] +
						2*b*b*para10[5]*para10[3]*para10[9]*para10[9]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[9]*cov66[5][9] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[9]*cov66[3][9];
		a4 = a4/4.0;
		a4u = (TMath::Sqrt(a4u))/4.0;
		printf("\n\n0+ Second Excited State Area: %.4f +/- %.4f", a4, a4u);


		a5 = para10[3]*para10[5]*para10[6]*b;;
		a5u = TMath::Power((b*para10[5]*para10[6]), 2)*cov66[3][3] +
						TMath::Power((b*para10[3]*para10[6]), 2)*cov66[5][5] +
							TMath::Power((b*para10[3]*para10[5]), 2)*cov66[6][6] +
								2*b*b*para10[5]*para10[3]*para10[6]*para10[6]*cov66[3][5] +
									2*b*b*para10[3]*para10[3]*para10[5]*para10[6]*cov66[5][6] +
										2*b*b*para10[3]*para10[5]*para10[5]*para10[6]*cov66[3][6];
		a5 = a5/4.0;
		a5u = (TMath::Sqrt(a5u))/4.0;
		printf("\n\n2+ First Excited State Area: %.4f +/- %.4f", a5, a5u);



		printf("\n\n");
		fprintf(f3, "\n\n");
		fclose(f3);
	return 0;
}


//18 degrees
int fit18D(int angle, TFile* f2, FILE* f3, TFile* fFit){

		TH1F	*tofTop = NULL;
		TH1F	*tofMid = NULL;
		TH1F	*tofBot = NULL;
		TH1F	*tofTopOut = NULL;
		TH1F	*tofMidOut = NULL;
		TH1F	*tofBotOut = NULL;
		char 	cutinSTR[][21] = CUTIN;
		char 	cutoutSTR[][21] = CUTOUT;
		char	txtSTR[100];

		// rebinned (div by x)
		float h = 4;

		gROOT->SetBatch(kTRUE);

		fprintf(f3, "\nYou have chosen to apply a fit to the combined TOF Spectra.\n");

		/* This portion of code grabs the TOF objects that will be fitted. */
		fFit->GetObject(cutinSTR[angle],tofTop);
		fFit->GetObject(cutinSTR[angle+1],tofMid);
		fFit->GetObject(cutinSTR[angle+2],tofBot);
		fFit->GetObject(cutoutSTR[angle],tofTopOut);
		fFit->GetObject(cutoutSTR[angle+1],tofMidOut);
		fFit->GetObject(cutoutSTR[angle+2],tofBotOut);

		TH1F *mycutzu1_clone = (TH1F*)tofTop->Clone();
		mycutzu1_clone->SetName("mycutzu1_clone");
		mycutzu1_clone->SetLineColor(kMagenta);
		mycutzu1_clone->Add(tofMid);
		mycutzu1_clone->Add(tofBot);
		mycutzu1_clone->Rebin(h);
		sprintf(txtSTR, "Combined %d Degree TOF with Applied Fits", angle);
		mycutzu1_clone->SetTitle(txtSTR);

		TH1F *gasout1_1clone = (TH1F*)tofTopOut->Clone();
		gasout1_1clone->SetName("gasout1_1clone");
		gasout1_1clone->SetLineColor(kAzure);
		gasout1_1clone->Add(tofMidOut);
		gasout1_1clone->Add(tofBotOut);
		gasout1_1clone->Rebin(h);


		gasout1_1clone->Scale(3.62); //product of ratio of gas in to gas out bcis and ratio of gas in to gas out live times

		f2->cd();

		/* Fit range parameters */

		Int_t lu1 = 2350;
		Int_t hu1 = 2600;


		sprintf(txtSTR, "1xCs_degree_%d_fitted_data", angle);
		TCanvas *c1 = new TCanvas(txtSTR);
		c1->cd();

		TH1F *gasin_subtract = (TH1F*)mycutzu1_clone->Clone();
		gasin_subtract->Sumw2();
		gasout1_1clone->Sumw2();
		gasin_subtract->Add(gasout1_1clone,-1);

		/*
		 * Creation of fit parameters
		 */

		Double_t params[12] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

		TF1 *fitfunc2 = new TF1("fitfunc2", fitfuncc, lu1, hu1, 12);
		fitfunc2->SetParameter(0,(5.1));  		//constant
		fitfunc2->FixParameter(1,0);  		//linear
		fitfunc2->FixParameter(2,0);  		//quad
		fitfunc2->SetParameter(3,11.7);  	//sigma
		fitfunc2->SetParameter(4,2485.0);  //mean
		fitfunc2->SetParameter(5,20);   //amp1

		//***************************************

		//parameter 6, first 2+ state
		//uncomment to include

		fitfunc2->SetParameter(6,0.001);     //amp2   2+ states
		fitfunc2->SetParLimits(6,0,1000);

		//uncomment to exclude

		//fitfunc2->FixParameter(6,0.0);

		//***************************************

		fitfunc2->SetParameter(7,0.13);    	//amp3  50->0.1
		fitfunc2->SetParameter(8,0.001);     //amp4   2+ states
		fitfunc2->SetParameter(9,0.11);    	//amp5  50->0.1
		fitfunc2->FixParameter(10,0.0);			//amp6
		fitfunc2->FixParameter(11,0.0);			//amp7


		fitfunc2->SetParLimits(0,0,40);
		fitfunc2->SetParLimits(5,0,1000);

		fitfunc2->SetParLimits(7,0,1000);
		fitfunc2->SetParLimits(8,0,1000);
		fitfunc2->SetParLimits(9,0,1000);
		//fitfunc2->SetParLimits(10,0,1000);
		//fitfunc2->SetParLimits(11,0,1000);

		TF1 *gausfuncu1 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu11 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu12 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu13 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu14 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu15 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);
		TF1 *gausfuncu16 = new TF1("gausfunc", gaussianpeak1, lu1, hu1, 3);

		TF1 *bkdfuncu1 = new TF1("bkdfunc",background, lu1, hu1, 3);

		printf("\n\n");
		TFitResultPtr fit_result = mycutzu1_clone->Fit("fitfunc2","RS");

		fitfunc2->GetParameters(params);

		bkdfuncu1->SetParameter(0,params[0]);   //constant
		bkdfuncu1->SetParameter(1,params[1]);   //linear
		bkdfuncu1->SetParameter(2,params[2]);   //quad

		gausfuncu1->SetParameter(0,params[3]);  //sigma
		gausfuncu1->SetParameter(1,params[4]);  //mean
		gausfuncu1->SetParameter(2,params[5]);  //amp


		gausfuncu11->SetParameter(0,params[3]);  //Sigma
		gausfuncu11->SetParameter(1,params[4] - 41.4);  //amp2
		gausfuncu11->SetParameter(2,params[5]*params[6]);  //amp3  //*p5

		gausfuncu12->SetParameter(0,params[3]);  //Sigma
		gausfuncu12->SetParameter(1,params[4] - 50.5);  //Mean
		gausfuncu12->SetParameter(2,params[5]*params[7]);  //amp5

		gausfuncu13->SetParameter(0,params[3]);  //Sigma
		gausfuncu13->SetParameter(1,params[4] - 68.2);  //Mean
		gausfuncu13->SetParameter(2,params[5]*params[8]);  //amp5

		gausfuncu14->SetParameter(0,params[3]);  //Sigma
		gausfuncu14->SetParameter(1,params[4] - 96.1);  //Mean
		gausfuncu14->SetParameter(2,params[5]*params[9]);   //amp5

		gausfuncu15->SetParameter(0,params[3]);  //Sigma
		gausfuncu15->SetParameter(1,params[4] - 108);  //Mean
		gausfuncu15->SetParameter(2,params[5]*params[10]);   //amp6

		gausfuncu16->SetParameter(0,params[3]);  //Sigma
		gausfuncu16->SetParameter(1,params[4] - 132);  //Mean
		gausfuncu16->SetParameter(2,params[5]*params[11]);   //amp6


		//legend titles

		mycutzu1_clone->SetBit(TH1::kNoTitle);


		mycutzu1_clone->SetTitle("Combined Gas In");
		mycutzu1_clone->SetLineColor(kMagenta);
		mycutzu1_clone->Draw();
		gStyle->SetOptFit(1);

		gasout1_1clone->SetTitle("Combined Gas Out");
		gasout1_1clone->SetLineColor(kYellow + 1);
		gasout1_1clone->Draw("same");

		gasin_subtract->SetTitle("Gas In Minus Gas Out");
		gasin_subtract->SetLineColor(kBlue - 7);
		gasin_subtract->Draw("SAME");					// Using "SAMES" changes the info box to this fit.


		fitfunc2->SetTitle("Complete Fit Function");
		fitfunc2->SetLineColor(kRed);
		fitfunc2->Draw("same");

		bkdfuncu1->SetTitle("Background");
		bkdfuncu1->SetLineColor(kGreen - 3);
		bkdfuncu1->Draw("same");

		gausfuncu1->SetTitle("Ground State");
		gausfuncu1->SetLineColor(kAzure);
		gausfuncu1->Draw("same");

		gausfuncu11->SetTitle("2+ First Excited State 1.52 MeV");
		gausfuncu11->SetLineColor(kCyan);
		gausfuncu11->Draw("same");

		gausfuncu12->SetTitle("0+ First Excited State 1.84 MeV");
		gausfuncu12->SetLineColor(kGreen);
		gausfuncu12->Draw("same");

		gausfuncu13->SetTitle("2+ Second Excited State 2.42 MeV");
		gausfuncu13->SetLineColor(kOrange - 3);
		gausfuncu13->Draw("same");

		gausfuncu14->SetTitle("0+ Second Excited State 3.30 MeV");
		gausfuncu14->SetLineColor(kAzure+9);
		gausfuncu14->Draw("same");

		gausfuncu15->SetTitle("3.65 MeV");
		gausfuncu15->SetLineColor(kRed+4);
		gausfuncu15->Draw("same");

		gausfuncu16->SetTitle("??? MeV/126 TDC offset");
		gausfuncu16->SetLineColor(kCyan-2);
		gausfuncu16->Draw("same");


		//title: formatting needs to be fixed
		TPaveText *t = new TPaveText(0.0, 0.9, 0.3, 1.0, "brNDC"); // left-up
		t->AddText("18 Degree Fit");
		t->SetBorderSize(0);
		t->SetFillColor(gStyle->GetTitleFillColor());
		t->Draw();

		fit_result->Print("V");

		c1->BuildLegend();
		c1->Update();
		c1->Write();

		gROOT->SetBatch(kFALSE);

		double chisk = fit_result->Chi2();
		unsigned int ndf = fit_result->Ndf();
		double edm = fit_result->Edm();
		unsigned int nCall = fit_result->NCalls();
		int i,j;

		double para10[10];
		double erro10[10];
		for(i=0; i<10; i++){
			para10[i] = fit_result->Parameter(i);
			erro10[i] = fit_result->Error(i);
		}

		double cor66[10][10];
		double cov66[10][10];

		for(i=0; i<10; i++){
			for(j=0; j<10; j++){
				cor66[i][j] = fit_result->Correlation(i, j);
				cov66[i][j] = fit_result->CovMatrix(i, j);
			}
		}

		fprintf(f3, "\nChi2\t=\t%.3e", chisk);
		fprintf(f3, "\nNDf\t=\t%u", ndf);
		fprintf(f3, "\nEdm\t=\t%.3e", edm);
		fprintf(f3, "\nNCalls\t=\t%u", nCall);
		fprintf(f3, "\np0\t=\t%.3e +/- %.3e", para10[0], erro10[0]);
		fprintf(f3, "\np1\t=\t%.3e +/- %.3e", para10[1], erro10[1]);
		fprintf(f3, "\np2\t=\t%.3e +/- %.3e", para10[2], erro10[2]);
		fprintf(f3, "\np3\t=\t%.3e +/- %.3e", para10[3], erro10[3]);
		fprintf(f3, "\np4\t=\t%.3e +/- %.3e", para10[4], erro10[4]);
		fprintf(f3, "\np5\t=\t%.3e +/- %.3e", para10[5], erro10[5]);
		fprintf(f3, "\np6\t=\t%.3e +/- %.3e", para10[6], erro10[6]);
		fprintf(f3, "\np7\t=\t%.3e +/- %.3e", para10[7], erro10[7]);
		fprintf(f3, "\np8\t=\t%.3e +/- %.3e", para10[8], erro10[8]);
		fprintf(f3, "\np9\t=\t%.3e +/- %.3e", para10[9], erro10[9]);


		/* Propagation of errors begins here */
		printf("\n\nPropagation of Errors: \n\n");
		printf("p3 --> %f \n", para10[3]);
		printf("p5 --> %f \n", para10[5]);
		printf("p7 --> %f \n", para10[7]);
		printf("p8 --> %f \n", para10[8]);
		printf("p9 --> %f \n", para10[9]);

		printf("p3u --> %f \n", cov66[3][3]);
		printf("p5u --> %f \n", cov66[5][5]);
		printf("p7u --> %f \n", cov66[7][7]);
		printf("p8u --> %f \n", cov66[8][8]);
		printf("p9u --> %f \n", cov66[9][9]);

		printf("p3p5 --> %f \n", cov66[3][5]);
		printf("p3p7 --> %f \n", cov66[3][7]);
		printf("p3p8 --> %f \n", cov66[3][8]);
		printf("p3p9 --> %f \n", cov66[3][9]);
		printf("p5p7 --> %f \n", cov66[5][7]);
		printf("p5p8 --> %f \n", cov66[5][8]);
		printf("p5p9 --> %f \n", cov66[5][9]);

		double b = TMath::Sqrt(2*TMath::Pi());
		double a1, a1u, a2, a2u, a3, a3u, a4, a4u, a5, a5u;

		a1 = para10[3]*para10[5]*b;
		a1u = TMath::Power((b*para10[5]), 2)*cov66[3][3] +
					TMath::Power((b*para10[3]), 2)*cov66[5][5] +
							2*b*b*para10[5]*para10[3]*cov66[3][5];
		a1 = a1/4.0;
		a1u = (TMath::Sqrt(a1u))/4.0;
		printf("\n\nGround State Area: %.4f +/- %.4f", a1, a1u);

		a2 = para10[3]*para10[5]*para10[7]*b;
		a2u = TMath::Power((b*para10[5]*para10[7]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[7]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[7][7] +
						2*b*b*para10[5]*para10[3]*para10[7]*para10[7]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[7]*cov66[5][7] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[7]*cov66[3][7];
		a2 = a2/4.0;
		a2u = (TMath::Sqrt(a2u))/4.0;
		printf("\n\n0+ First Excited State Area: %.4f +/- %.4f", a2, a2u);

		a3 = para10[3]*para10[5]*para10[8]*b;
		a3u = TMath::Power((b*para10[5]*para10[8]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[8]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[8][8] +
						2*b*b*para10[5]*para10[3]*para10[8]*para10[8]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[8]*cov66[5][8] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[8]*cov66[3][8];
		a3 = a3/4.0;
		a3u = (TMath::Sqrt(a3u))/4.0;
		printf("\n\n2+ Second Excited State Area: %.4f +/- %.4f", a3, a3u);

		a4 = para10[3]*para10[5]*para10[9]*b;
		a4u = TMath::Power((b*para10[5]*para10[9]), 2)*cov66[3][3] +
				TMath::Power((b*para10[3]*para10[9]), 2)*cov66[5][5] +
					TMath::Power((b*para10[3]*para10[5]), 2)*cov66[9][9] +
						2*b*b*para10[5]*para10[3]*para10[9]*para10[9]*cov66[3][5] +
							2*b*b*para10[3]*para10[3]*para10[5]*para10[9]*cov66[5][9] +
								2*b*b*para10[3]*para10[5]*para10[5]*para10[9]*cov66[3][9];
		a4 = a4/4.0;
		a4u = (TMath::Sqrt(a4u))/4.0;
		printf("\n\n0+ Second Excited State Area: %.4f +/- %.4f", a4, a4u);


		a5 = para10[3]*para10[5]*para10[6]*b;;
		a5u = TMath::Power((b*para10[5]*para10[6]), 2)*cov66[3][3] +
						TMath::Power((b*para10[3]*para10[6]), 2)*cov66[5][5] +
							TMath::Power((b*para10[3]*para10[5]), 2)*cov66[6][6] +
								2*b*b*para10[5]*para10[3]*para10[6]*para10[6]*cov66[3][5] +
									2*b*b*para10[3]*para10[3]*para10[5]*para10[6]*cov66[5][6] +
										2*b*b*para10[3]*para10[5]*para10[5]*para10[6]*cov66[3][6];
		a5 = a5/4.0;
		a5u = (TMath::Sqrt(a5u))/4.0;
		printf("\n\n2+ First Excited State Area: %.4f +/- %.4f", a5, a5u);



		printf("\n\n");
		fprintf(f3, "\n\n");
		fclose(f3);

	return 0;
}
