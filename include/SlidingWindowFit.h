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
#include "TMath.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"

#include <vector>
#include <iostream>

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
    int nentries;
    int iter;
    int seed;
    int numCPU;
    RooWorkspace *ws;
    RooDataSet *r_data;
    TFile *fOutput;
    TDirectory *glob;
    bool savePlot;
    bool silent;
    RooPlot *rplot;

public:
    SlidingWindowFit(double v1, double v2){
        emin = v1;
        emax = v2;
        ws = new RooWorkspace("ws");
        sigma_w = 0.01;
        savePlot = true;
        silent = false;
        seed = 1111;
        RooPlot *rplot = new RooPlot();
    }
    void setSeed(int val){
        seed = val;
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
        delete rplot;
        delete ws;
    };

    int getEntries(){ return nentries; }

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

    void fit(bool doToyMC = false);

private :

    ClassDef(SlidingWindowFit,1) // SlidingWindowFit

};



ClassImp(SlidingWindowFit);

#endif /* SlidingWindowFit */
