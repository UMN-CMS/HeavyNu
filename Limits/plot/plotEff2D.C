//root -b -q -l upperLimit.C+
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TText.h"
#include "TLatex.h"

#include "tdrstyle.C"
#include "../acceptance/acceptance_db.C"

#include <cstdio>

void plotEff2D(int channel)
{
    setTDRStyle();

    gStyle->SetPalette(1);

    AcceptanceDB *db = new AcceptanceDB("2012");

    double elecEff[] = {0.40991, 0.43299, 0.46256, 0.47850, 0.49682, 0.50857, 0.51957, 0.53121, 0.53734, 0.54517, 0.55140, 0.54497, 0.55318, 0.56327, 0.56309, 0.56904, 0.56930, 0.57698, 0.57210, 0.57619, 0.58068, 0.57916, 0.57855, 0.57927, 0.58021, 0.58921, 0.58235, 0.58422, 0.59192, 0.58760, 0.58752 };
    double muonEff[] = {0.48311, 0.52121, 0.54649, 0.56841, 0.58167, 0.60076, 0.61871, 0.63428, 0.63401, 0.64491, 0.65453, 0.65548, 0.66099, 0.66346, 0.67203, 0.67863, 0.68455, 0.67886, 0.68715, 0.69472, 0.68600, 0.69171, 0.69157, 0.69888, 0.69309, 0.69731, 0.69738, 0.69302, 0.70030, 0.69767, 0.70181 };

    TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);
    c1->SetMargin(0.18, 0.13, 0.13, 0.08);
    TH2 * h = new TH2D("h", "h", 21, 950, 3050, 30, 50, 3050);

    for(int iWr = 1; iWr <= 21; iWr++)
    {
        for(int iNu = 1; iNu <= 30; iNu++)
        {
            h->SetBinContent(iWr, iNu, -10000);
        }
    }

    for(int iWr = 1; iWr <= 21; iWr++)
    {
        for(int iNu = 1; iNu <= 30; iNu++)
        {
            int mWR = 900 + iWr * 100, mNu = iNu * 100;
            double accCor = db->getBestEstimate(mWR, mNu, 2012) / db->getBestEstimate(mWR, mWR / 2, 2012);
            double eff = 1.0;
            if(channel == 0) eff = muonEff[iWr];
            else if(channel == 1) eff = elecEff[iWr];
            h->SetBinContent(iWr, iNu, eff * accCor);
        }
    }

    h->GetZaxis()->SetRangeUser(0, 0.70);
    h->GetXaxis()->SetTitle("M_{W_{R}} [GeV]");
    if(channel == 0) h->GetYaxis()->SetTitle("M_{N_{#mu}} [GeV]");
    else if(channel == 1) h->GetYaxis()->SetTitle("M_{N_{e}} [GeV]");
    h->GetYaxis()->SetTitleOffset(1.50);
    
    char buffer[64];
    sprintf(buffer, "CMS Preliminary    #sqrt{s} = %d TeV    %.1f fb^{-1}", 8, 19.7);
    TLatex *tex = new TLatex(0.180, 0.94, buffer);
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);

    h->Draw("colz");
    tex->Draw("same");
    
    if(channel == 0) 
    {
        c1->Print("eff2D_muon.pdf");
        c1->Print("eff2D_muon.png");
    }
    else if(channel == 1) 
    {
        c1->Print("eff2D_elec.pdf");
        c1->Print("eff2D_elec.png");
    }
}
