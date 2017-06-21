#ifndef SLIDINGWINDOWFIT_H
#define SLIDINGWINDOWFIT_H

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include "RooBinning.h"
#include "TMath.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPad.h"
#include "TString.h"
//#include "RooTFnBinding.h"
#include <vector>
#include <iostream>
#include <string>

using namespace RooFit;

class SlidingWindowFit {

private:
    double emin;
    double emax;
    double ecenter;
    double fscale;
    double sigma_w;
    double index[3];
    double chi2;
    int nbins;
    int nentries;
    int iter;
    int seed;
    int numCPU;
    TString *pdfname_fit;
    RooWorkspace *ws;
    RooDataSet *r_data;
    RooBinning *custom_binning;
    bool use_custom_binning;
    TFile *fOutput;
    TDirectory *glob;
    bool savePlot;
    bool silent;
    bool is_binned;
    bool include_proton_bkg;
    bool include_eff;

public:
    SlidingWindowFit(double v1, double v2){
        emin = v1;
        emax = v2;
        nbins = 50;
        ws = new RooWorkspace("ws");
        sigma_w = 0.01;
        savePlot = true;
        silent = false;
        include_proton_bkg = false;
        include_eff = false;
        custom_binning = new RooBinning(emin,emax);
        seed = 1111;
        use_custom_binning = false;
        pdfname_fit = new TString("bmodel");
        is_binned = false;
    }
    void addSignalEfficiency(char fname[128], char histname[64]);

    void addProtonBkg(char fname[128], char histname[64]);

    void setNbins(int val){
        nbins = val;
    }

    void setSeed(int val){
        seed = val;
    }

    void setBinnedData(bool val){
        is_binned = val;
    }

    void setSignalPdf(TString val){
        pdfname_fit = new TString(val.Data());
    }

    void setPlotSave(bool val){
        savePlot = val;
    }


    RooWorkspace *getWorkspace(){
        return ws;
    }

    void setNumCPU(int val){
        numCPU = val;
    }

    void setSilent(){
        silent = true;
        if (silent) RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    }

    void setOutputData(TFile *fData){
        fOutput = fData;
        fOutput->cd();
        glob = gDirectory;
    }

    double getResult(int idx){
        // returns index, parameter & chi2
        if (idx > 3){
            std::cout << "result is array of length 4, not allowed query" << std::endl;
            throw;
        }
        std::vector<double> data;
        data.push_back(index[0]);
        data.push_back(index[1]);
        data.push_back(index[2]);
        data.push_back(chi2);
        return data.at(idx);
    }

    virtual ~SlidingWindowFit(){

    };

    int getEntries(){ return nentries; }

    void setLogBinning(int ival){
        nbins = ival;
        use_custom_binning = true;
        double logEmin = TMath::Log10(emin);
        double logEmax = TMath::Log10(emax);
        double delta = (logEmax - logEmin) / float(nbins-1);
        double val = logEmin;
        double valE = pow(10.,val);
        while (val < logEmax){
            val+=delta;
            valE=pow(10.,val);
            if (!silent) std::cout << "valE: " << valE <<  " logE: " << val << std::endl;
            custom_binning->addBoundary(valE);
        }
    }

    void setCustomBinning(RooBinning *binning){
        use_custom_binning = true;
        custom_binning = binning;
        nbins = custom_binning->numBins();
    }

    double getSignalFraction(int idx=0){
        RooRealVar *frac = ws->var("frac");
        if (idx == 0) return frac->getVal();
        else return frac->getError();
    }

    void setMin(TString *var, double val){
        ws->var(var->Data())->setMin(val);
    }

    void setMax(TString *var, double val){
        ws->var(var->Data())->setMax(val);
    }

    void setVal(TString *var, double val){
        ws->var(var->Data())->setVal(val);
    }

    void setFixed(TString *var, bool val){
        ws->var(var->Data())->setConstant(val);
    }

    void setIteration(int val){ iter = val; }

    void setData(char fname[128], char tname[24]);

    void setEnergyResolution(double val){ sigma_w = val; }

    void setWindowDef(double v1, double v2, double v3){
        emin   = v1;
        ecenter= v2;
        emax   = v3;
    }

    void setWindowFscale(double val){
        fscale = val;
    }

    void buildModel();

    void toyMC(char signal_pdf[64], int ntoys = 1000);
    void fit(bool doToyMC);

private :

    ClassDef(SlidingWindowFit,1) // SlidingWindowFit

};



ClassImp(SlidingWindowFit);

#endif /* SlidingWindowFit */
