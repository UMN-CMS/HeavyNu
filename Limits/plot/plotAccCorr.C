//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////root -b -q -l upperLimit.C+
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TText.h"
#include "TLatex.h"

#include "tdrstyle.C"
#include "../acceptance/acceptance_db.C"

#include <cstdio>
#include <map>
#include <utility>

int mwrs[] = {1000, 1500, 2000, 3000};
const int nmwrs = sizeof(mwrs) / sizeof(int);

TF1 * fits[nmwrs];
std::map<int, std::map<int, double> > uncertMap;

double getAccRatio(double mWr, double mNu);

void plotAccCorr(int channel, int mWr = 1500)
{
    setTDRStyle();

    gStyle->SetPalette(1);

    AcceptanceDB *db = 0;

    std::map<int, std::map<int, double> > recoAccMap;
    if(channel == 0)
    {
        db = new AcceptanceDB("2012_muon");
                 
        recoAccMap[1000][62]   = 0.007;
        recoAccMap[1000][125]  = 0.086;
        recoAccMap[1000][187]  = 0.208;
        recoAccMap[1000][250]  = 0.304;
        recoAccMap[1000][500]  = 0.483;
        recoAccMap[1000][833]  = 0.443;
        recoAccMap[1500][93]   = 0.017;
        recoAccMap[1500][187]  = 0.140;
        recoAccMap[1500][281]  = 0.296;
        recoAccMap[1500][375]  = 0.411;
        recoAccMap[1500][750]  = 0.601;
        recoAccMap[1500][1250] = 0.592;
        recoAccMap[2000][125]  = 0.027;
        recoAccMap[2000][250]  = 0.171;
        recoAccMap[2000][375]  = 0.355;
        recoAccMap[2000][500]  = 0.471;
        recoAccMap[2000][1000] = 0.655;
        recoAccMap[2000][1666] = 0.649;
        recoAccMap[3000][187]  = 0.075;
        recoAccMap[3000][375]  = 0.258;
        recoAccMap[3000][562]  = 0.435;
        recoAccMap[3000][750]  = 0.554;
        recoAccMap[3000][1500] = 0.686;
        recoAccMap[3000][2500] = 0.692;

        recoAccMap[700][350] = 0.30173;
        recoAccMap[800][400] = 0.38558;
        recoAccMap[900][450] = 0.44423;
        recoAccMap[1100][550] = 0.52121;
        recoAccMap[1200][600] = 0.54649;
        recoAccMap[1300][650] = 0.56841;
        recoAccMap[1400][700] = 0.58167;
        recoAccMap[1600][800] = 0.61871;
        recoAccMap[1700][850] = 0.63428;
        recoAccMap[1800][900] = 0.63401;
        recoAccMap[1900][950] = 0.64491;
        recoAccMap[2100][1050] = 0.65548;
        recoAccMap[2200][1100] = 0.66099;
        recoAccMap[2300][1150] = 0.66346;
        recoAccMap[2400][1200] = 0.67203;
        recoAccMap[2500][1250] = 0.67863;
        recoAccMap[2600][1300] = 0.68455;
        recoAccMap[2700][1350] = 0.67886;
        recoAccMap[2800][1400] = 0.68715;
        recoAccMap[2900][1450] = 0.69472;
        recoAccMap[3100][1550] = 0.69171;
        recoAccMap[3200][1600] = 0.69157;
        recoAccMap[3300][1650] = 0.69888;
        recoAccMap[3400][1700] = 0.69309;
        recoAccMap[3500][1750] = 0.69731;
        recoAccMap[3600][1800] = 0.69738;
        recoAccMap[3700][1850] = 0.69302;
        recoAccMap[3800][1900] = 0.70030;
        recoAccMap[3900][1950] = 0.69767;
        recoAccMap[4000][2000] = 0.70181;


        uncertMap[1000][62]   = 1.095891;
        uncertMap[1000][125]  = 1.028128;
        uncertMap[1000][187]  = 1.018112;
        uncertMap[1000][250]  = 1.015022;
        uncertMap[1000][500]  = 1.007297;
        uncertMap[1000][833]  = 1.012223;
        uncertMap[1500][93]   = 1.060321;
        uncertMap[1500][187]  = 1.021925;
        uncertMap[1500][281]  = 1.015186;
        uncertMap[1500][375]  = 1.012939;
        uncertMap[1500][750]  = 1.00651;
        uncertMap[1500][1250] = 1.010912;
        uncertMap[2000][125]  = 1.048666;
        uncertMap[2000][250]  = 1.019827;
        uncertMap[2000][375]  = 1.01622;
        uncertMap[2000][500]  = 1.011948;
        uncertMap[2000][1000] = 1.009828;
        uncertMap[2000][1666] = 1.010292;
        uncertMap[3000][187]  = 1.028286;
        uncertMap[3000][375]  = 1.016066;
        uncertMap[3000][562]  = 1.01258;
        uncertMap[3000][750]  = 1.011183;
        uncertMap[3000][1500] = 1.007214;
        uncertMap[3000][2500] = 1.009941;
    }
    else if(channel == 1)
    {
        db = new AcceptanceDB("2012_elec");
        
        recoAccMap[1000][62]   = 0.005863;
        recoAccMap[1000][125]  = 0.062181;
        recoAccMap[1000][187]  = 0.170653;
        recoAccMap[1000][250]  = 0.250767;
        recoAccMap[1000][500]  = 0.409907;
        recoAccMap[1000][833]  = 0.377353;
        recoAccMap[1500][93]   = 0.011335;
        recoAccMap[1500][187]  = 0.106364;
        recoAccMap[1500][281]  = 0.240461;
        recoAccMap[1500][375]  = 0.345594;
        recoAccMap[1500][750]  = 0.508572;
        recoAccMap[1500][1250] = 0.502404;
        recoAccMap[2000][125]  = 0.023142;
        recoAccMap[2000][250]  = 0.141413;
        recoAccMap[2000][375]  = 0.289359;
        recoAccMap[2000][500]  = 0.3893;
        recoAccMap[2000][1000] = 0.551397;
        recoAccMap[2000][1666] = 0.535799;
        recoAccMap[3000][187]  = 0.063246;
        recoAccMap[3000][375]  = 0.211899;
        recoAccMap[3000][562]  = 0.360551;
        recoAccMap[3000][750]  = 0.457733;
        recoAccMap[3000][1500] = 0.580682;
        recoAccMap[3000][2500] = 0.580486;

        recoAccMap[1100][550] = 0.43299;
        recoAccMap[1200][600] = 0.46256;
        recoAccMap[1300][650] = 0.47850;
        recoAccMap[1400][700] = 0.49682;
        recoAccMap[1600][800] = 0.51957;
        recoAccMap[1700][850] = 0.53121;
        recoAccMap[1800][900] = 0.53734;
        recoAccMap[1900][950] = 0.54517;
        recoAccMap[2100][1050] = 0.54497;
        recoAccMap[2200][1100] = 0.55318;
        recoAccMap[2300][1150] = 0.56327;
        recoAccMap[2400][1200] = 0.56309;
        recoAccMap[2500][1250] = 0.56904;
        recoAccMap[2600][1300] = 0.56930;
        recoAccMap[2700][1350] = 0.57698;
        recoAccMap[2800][1400] = 0.57210;
        recoAccMap[2900][1450] = 0.57619;
        recoAccMap[3100][1550] = 0.57916;
        recoAccMap[3200][1600] = 0.57855;
        recoAccMap[3300][1650] = 0.57927;
        recoAccMap[3400][1700] = 0.58021;
        recoAccMap[3500][1750] = 0.58921;
        recoAccMap[3600][1800] = 0.58235;
        recoAccMap[3700][1850] = 0.58422;
        recoAccMap[3800][1900] = 0.59192;
        recoAccMap[3900][1950] = 0.58760;
        recoAccMap[4000][2000] = 0.58752;


        uncertMap[1000][62]   = 1.108679;
        uncertMap[1000][125]  = 1.033026;
        uncertMap[1000][187]  = 1.020082;
        uncertMap[1000][250]  = 1.016411;
        uncertMap[1000][500]  = 1.008008;
        uncertMap[1000][833]  = 1.013416;
        uncertMap[1500][93]   = 1.074489;
        uncertMap[1500][187]  = 1.025311;
        uncertMap[1500][281]  = 1.017055;
        uncertMap[1500][375]  = 1.014136;
        uncertMap[1500][750]  = 1.007175;
        uncertMap[1500][1250] = 1.011933;
        uncertMap[2000][125]  = 1.052579;
        uncertMap[2000][250]  = 1.022152;
        uncertMap[2000][375]  = 1.018111;
        uncertMap[2000][500]  = 1.013424;
        uncertMap[2000][1000] = 1.010762;
        uncertMap[2000][1666] = 1.011407;
        uncertMap[3000][187]  = 1.031147;
        uncertMap[3000][375]  = 1.017772;
        uncertMap[3000][562]  = 1.013851;
        uncertMap[3000][750]  = 1.012354;
        uncertMap[3000][1500] = 1.007975;
        uncertMap[3000][2500] = 1.01108;

    }
    
    if(db == 0) return;

    TCanvas * c1[nmwrs];
    //c1->SetMargin(0.18, 0.13, 0.13, 0.08);

    TH2 * dummy[nmwrs];

    TGraph * tgGen[nmwrs];
    TGraphErrors * tgReco[nmwrs], *tgRatio[nmwrs];

    TLegend * leg[nmwrs];
    
    TGraphErrors *tgRatio2 = new TGraphErrors();
    tgRatio2->SetMarkerStyle(24);
    tgRatio2->SetMarkerColor(kRed);

    TGraph2DErrors *tg2D = new TGraph2DErrors();
    int ctr2D = 0;
    //int ctr = 0;

    for(int iWr = 0; iWr < nmwrs; iWr++)
    {
        char cname[32];
        sprintf(cname, "c%d", iWr);
        c1[iWr] = new TCanvas(cname, "c1", 800, 800);
        c1[iWr]->SetMargin(0.15, 0.06, 0.15, 0.10);

//        dummy[iWr] = new TH2D("dummy", "dummy", 1000, 0, mwrs[iWr], 1000, 0, 1.15);
        dummy[iWr] = new TH2D("dummy", "dummy", 1000, 0, 1.0, 1000, 0, 1.15);
        if(channel == 0) dummy[iWr]->GetXaxis()->SetTitle("M_{N_{#mu}} / M_{W_{R}}");
        else if(channel == 1) dummy[iWr]->GetXaxis()->SetTitle("M_{N_{e}} / M_{W_{R}}");
        dummy[iWr]->GetXaxis()->SetTitleOffset(1.05);
        dummy[iWr]->GetYaxis()->SetTitle("Efficiency");
        dummy[iWr]->GetYaxis()->SetTitleOffset(1.00);
        dummy[iWr]->GetYaxis()->SetNdivisions(3, 5, 0);
        dummy[iWr]->SetStats(0);
        dummy[iWr]->SetTitle(0);
        dummy[iWr]->GetXaxis()->SetNdivisions(6, 5, 0);

        tgGen[iWr] = new TGraph();
        tgReco[iWr] = new TGraphErrors();
        tgRatio[iWr] = new TGraphErrors();
        tgRatio[iWr]->SetMarkerStyle(24);
        tgRatio[iWr]->SetMarkerColor(kRed);

        fits[iWr] = new TF1("fit", "1-expo(0)", 0, mwrs[iWr]);
        //fits[iWr] = new TF1("fit", "1-exp([0]*x)", 0, mwrs[iWr]);
        fits[iWr]->SetParLimits(0, -10.0, 0.0);

        leg[iWr] = new TLegend(0.5, 0.17, 0.83, 0.4);
        leg[iWr]->SetFillStyle(0);
        leg[iWr]->SetBorderSize(0);
        leg[iWr]->AddEntry(tgGen[iWr], "Gen Eff", "L");
        leg[iWr]->AddEntry(tgReco[iWr], "Reco Eff", "P");
        leg[iWr]->AddEntry(tgRatio[iWr], "Eff Ratio", "P");
        leg[iWr]->AddEntry(fits[iWr], "Eff Ratio Fit", "L");

        //double extrapoX[] = {125.0, 250.0, 375.0, 500.0, 1000.0, 1666.0};
        //double extrapoY[] = {0.918369, 0.986134, 0.998650, 1.000348, 1.000022, 1.000000};
        //double extrapoY[] = {0.873972, 0.952278, 0.973968, 0.996863, 0.999263, 0.999982};
        //double extrapoX[] = {93.0, 187.0, 281.0, 375.0, 750.0, 1250.0};
        //double extrapoY[] = {0.850008, 0.881803, 0.936398, 0.974038, 0.997375, 0.999784};
        //TGraph *extra = new TGraph(6, extrapoX, extrapoY);

        double eff = recoAccMap[mwrs[iWr]][mwrs[iWr] / 2];

        for(int iNu = 1; iNu <= 30; iNu++)
        {
            int mNu = iNu * 100;
            double accCor = db->getBestEstimate(mwrs[iWr], mNu, 2012) / db->getBestEstimate(mwrs[iWr], mwrs[iWr] / 2, 2012);
            if(accCor > 0) tgGen[iWr]->SetPoint(iNu - 1, (double)mNu/mwrs[iWr], eff * accCor);
        }

        int ctr = 0;
        for(std::map<int, double>::const_iterator mNu = recoAccMap[mwrs[iWr]].begin(); mNu != recoAccMap[mwrs[iWr]].end(); ++mNu)
        {
            tgReco[iWr]->SetPoint(ctr, double(mNu->first)/mwrs[iWr], mNu->second);
            double ratio = mNu->second / (eff * db->getBestEstimate(mwrs[iWr], mNu->first, 2012) / db->getBestEstimate(mwrs[iWr], mwrs[iWr] / 2, 2012));
            std::cout << ratio << std::endl;
            double error = (uncertMap[mwrs[iWr]][mNu->first] - 1) * ratio;
            tgRatio[iWr]->SetPoint(ctr, double(mNu->first)/mwrs[iWr], ratio);
            tgReco[iWr]->SetPointError(ctr, 0.0, (uncertMap[mwrs[iWr]][mNu->first] - 1) * mNu->second);
            tgRatio[iWr]->SetPointError(ctr, 0.0, error);
            
            tgRatio2->SetPoint(ctr2D, double(mNu->first) / mwrs[iWr], ratio);
            tgRatio2->SetPointError(ctr2D, 0.0, error);
            
            tg2D->SetPoint(ctr2D, mwrs[iWr], double(mNu->first) / mwrs[iWr], ratio);
            tg2D->SetPointError(ctr2D, 0.0, 0.0, error);
            ctr++;
            ctr2D++;
        }
        printf("\n");

        tgRatio[iWr]->Fit(fits[iWr], "QN");
        tgRatio[iWr]->Fit(fits[iWr], "QN");
        tgRatio[iWr]->Fit(fits[iWr], "QN");
        tgRatio[iWr]->Fit(fits[iWr], "QN");
        
        printf("%d\t%f\t%f\n", mwrs[iWr], fits[iWr]->GetParameter(0), fits[iWr]->GetParameter(1));

        dummy[iWr]->Draw();
        leg[iWr]->Draw("same");

        tgGen[iWr]->Draw("SAME L");
        tgReco[iWr]->Draw("SAME P");
        tgRatio[iWr]->Draw("SAME P");
        fits[iWr]->Draw("same");

        char fname[128] = "";
        if(channel == 0) sprintf(fname, "eff_gen_vs_reco_mWr_mm_%d.png", mwrs[iWr]);
        else if(channel == 1) sprintf(fname, "eff_gen_vs_reco_mWr_ee_%d.png", mwrs[iWr]);
        c1[iWr]->Print(fname);
        if(channel == 0) sprintf(fname, "eff_gen_vs_reco_mWr_mm_%d.pdf", mwrs[iWr]);
        else if(channel == 1) sprintf(fname, "eff_gen_vs_reco_mWr_ee_%d.pdf", mwrs[iWr]);
        c1[iWr]->Print(fname);
    }

    printf("%f\t%f\t%f\t%f\t%f\t%f\n",
           getAccRatio(2000.0, 125.0),
           getAccRatio(2000.0, 250.0),
           getAccRatio(2000.0, 375.0),
           getAccRatio(2000.0, 500.0),
           getAccRatio(2000.0, 1000.0),
           getAccRatio(2000.0, 1666.0));
    
    TCanvas *c44 = new TCanvas("c44", "c1", 800, 800);
    c44->SetMargin(0.15, 0.06, 0.15, 0.10);
 
    TH1 *dummy2 = new TH2D("dummy2", "dummy2", 1000, 0, 1.0, 1000, 0, 1.15);
    if(channel == 0) dummy2->GetXaxis()->SetTitle("M_{N_{#mu}} / M_{W_{R}}");
    else if(channel == 1) dummy2->GetXaxis()->SetTitle("M_{N_{e}} / M_{W_{R}}");
    dummy2->GetXaxis()->SetTitleOffset(1.05);
    dummy2->GetYaxis()->SetTitle("Efficiency");
    dummy2->GetYaxis()->SetTitleOffset(1.00);
    dummy2->GetYaxis()->SetNdivisions(3, 5, 0);
    dummy2->SetStats(0);
    dummy2->SetTitle(0);
    dummy2->GetXaxis()->SetNdivisions(6, 5, 0);
 
    TF1 *fit2 = new TF1("fit2", "1-expo(0)", 0, 1);
    fit2->SetParLimits(0, -10, 0);
 
    TLegend *leg2 = new TLegend(0.3, 0.17, 0.83, 0.4);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->AddEntry(tgRatio2, "Eff Ratio", "P");
    
    tgRatio2->Fit(fit2, "NQ");
    tgRatio2->Fit(fit2, "NQ");
    tgRatio2->Fit(fit2, "NQ");
    tgRatio2->Fit(fit2, "NQ");    
    c44->cd();
    dummy2->Draw();
    tgRatio2->Draw("same P");
    fit2->Draw("same L");
    
    char leglabel[128];
    sprintf(leglabel, "Fit (a = %f, b = %f, #chi^{2} = %f / %d)", fit2->GetParameter(0), fit2->GetParameter(1), fit2->GetChisquare(), fit2->GetNDF());
    leg2->AddEntry(fit2, leglabel, "L");
    
    printf("Fit parameters: a = %f, b = %f, #chi^{2} = %f / %d\n", fit2->GetParameter(0), fit2->GetParameter(1), fit2->GetChisquare(), fit2->GetNDF());
    
    double mins[6], maxs[6], merrs[6];
    for(int i = 0; i < 6; i++)
    {
        mins[i] = 9.9e99;
        maxs[i] = 0.0;
        merrs[i] = 0.0;
    }
    
    for(int iWr = 0; iWr < nmwrs; iWr++)
    {
        double eff = recoAccMap[mwrs[iWr]][mwrs[iWr] / 2];
        
        int iNu = 0;
        for(std::map<int, double>::const_iterator mNu = recoAccMap[mwrs[iWr]].begin(); mNu != recoAccMap[mwrs[iWr]].end(); ++mNu)
        {
            double ratio = mNu->second / (eff * db->getBestEstimate(mwrs[iWr], mNu->first, 2012) / db->getBestEstimate(mwrs[iWr], mwrs[iWr] / 2, 2012));
            double error = (uncertMap[mwrs[iWr]][mNu->first] - 1) * ratio;
            
            mins[iNu] = std::min(mins[iNu], ratio);
            maxs[iNu] = std::max(maxs[iNu], ratio);
            merrs[iNu] = std::max(merrs[iNu], error);
            iNu++;
        }
    }
    
    double errX[] = {1.0/16, 1.0/8, 3.0/16, 1.0/4, 1.0/2, 5.0/6};
    double errY[6];
    for(int i = 0; i < 6; i++) errY[i] = std::max((maxs[i] - mins[i]) / 2, merrs[i]);
    TGraph *tgerr2 = new TGraph(6, errX, errY);
    
    TF1 *tef = new TF1("tef", "expo", 0, 1);
    tgerr2->Fit(tef, "NQ");
    TF1 *errfit = new TF1("rf","expo(0) + [2]", 0, 1);
    errfit->SetLineColor(kBlue);
    errfit->SetParameter(0, tef->GetParameter(0));
    errfit->SetParameter(1, tef->GetParameter(1));
    errfit->SetParameter(2, 0.05);
    tgerr2->Fit(errfit, "NQ");
    tgerr2->Fit(errfit, "NQ");
    tgerr2->Fit(errfit, "NQ");
    tgerr2->Fit(errfit, "NQ");
    tgerr2->Draw("same P");
    errfit->Draw("same L");
    
    leg2->AddEntry(tgerr2, "Systematic Estimate", "P");
    leg2->AddEntry(errfit, "Systematic Fit", "L");
    leg2->Draw("same");

    char fname[128] = "";
    if(channel == 0)      sprintf(fname, "eff_gen_vs_reco_totalFit_mm.png");
    else if(channel == 1) sprintf(fname, "eff_gen_vs_reco_totalFit_ee.png");
    c44->Print(fname);    
    if(channel == 0)      sprintf(fname, "eff_gen_vs_reco_totalFit_mm.pdf");
    else if(channel == 1) sprintf(fname, "eff_gen_vs_reco_totalFit_ee.pdf");
    c44->Print(fname);

    
    printf("error fit parameters: a = %f, b = %f, c = %f, #chi^{2} = %f / %d\n", errfit->GetParameter(0), errfit->GetParameter(1), errfit->GetParameter(2), errfit->GetChisquare(), errfit->GetNDF());
    for(double rNu = 0.05; rNu < 1.0; rNu+=0.05) printf("%f\t", errfit->Eval(rNu));
    printf("\n");
}

double getAccRatio(double mWr, double mNu)
{
    double rNu = mNu / mWr;

    TGraphErrors *tg = new TGraphErrors();

    for(int iWr = 0; iWr < nmwrs; iWr++)
    {
        tg->SetPoint(iWr, mwrs[iWr], fits[iWr]->Eval(mwrs[iWr] * rNu));
        double uncert1 = -1, uncert2 = -1;
        double mNu1 = 0, mNu2 = 0;
        std::map<int, double> &mum = uncertMap[mwrs[iWr]];
        std::map<int, double>::const_iterator itLast = --mum.end();
        for(std::map<int, double>::const_iterator it = mum.begin(); it != mum.end(); ++it)
        {
            if(uncert1 > 0)
            {
                mNu2 = it->first;
                uncert2 = it->second;
                break;
            }
            else if(rNu < (double(it->first) / mwrs[iWr]) || it == itLast)
            {
                mNu1 = it->first;
                uncert1 = it->second;
            }
        }
        if(uncert1 > 0 && uncert2 > 0)
        {
            double m = (uncert2 - uncert1) / (mNu2 - mNu1);
            tg->SetPointError(iWr, 0.0, uncert1 + (rNu * mwrs[iWr] - mNu1) * m - 1);
        }
        else
        {
            tg->SetPointError(iWr, 0.0, uncert1-1);
        }
    }

    //TF1 *fit = new TF1("fit", "[0]-exp([1]-x*[2])", 1000, 3000);
    TF1 *fit = new TF1("fit", "pol2", 1000, 3000);
    tg->Fit(fit, "NQ");
    tg->Fit(fit, "NQ");

    //TCanvas *c24 = new TCanvas("c24", "c24", 800, 800);
    //TH1 *dummy = new TH1D("there", "there", 1000, 0, 4000);
    //dummy->Draw();
    //tg->Draw("SAME P");
    //fit->Draw("same");

    return fit->Eval(mWr);
}
