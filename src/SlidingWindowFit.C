#include "SlidingWindowFit.h"

void SlidingWindowFit::setData(char fname[128], char tname[24]){
    if (!silent){
        std::cout << "file name: " << fname << std::endl;
        std::cout << "tree name: " << tname << std::endl;
    }
    // fname be the name of the file, tname the name of the tree.
    TFile *fIn = TFile::Open(fname);
    if (!silent) {
        fIn->Print();
        fIn->ls();
    }
    if ( !fIn->GetListOfKeys()->Contains(tname) ) {
        std::cout << "ERROR: could not find tree: " << tname << std::endl;
        throw;
    }
    TTree *ttree_data = (TTree*)fIn->Get(tname);
    if (!silent) ttree_data->Print();
    //RooRealVar *w = ws->var("weight");
    RooRealVar *E = ws->var("E");
    RooDataSet *data = new RooDataSet("data", "data", RooArgSet(*E), Import(*ttree_data),Cut("E <= 2e3 && E > 25."));
    ws->import(*data);
}

void SlidingWindowFit::buildModel(){
    RooRealVar E("E","E",emin,emax);
    ws->import(E);
    double f_nobs = 1e8;
    RooRealVar gamma("gamma", "#gamma", 3.07, 1, 6);
    RooRealVar scale("scale", "scale" , 1.0, 1e-4,1e4);
    //RooRealVar weight("weight","weight",1.0,0.,1.0);
    RooRealVar norm("norm","norm",f_nobs, f_nobs - TMath::Sqrt(f_nobs), f_nobs + TMath::Sqrt(f_nobs));
    ws->import(gamma);
    ws->import(scale);
    ws->import(norm);
    //ws->import(weight);
    if (!silent) ws->Print();
    // build powerlaw
    //RooGenericPdf pwl("pwl","powerlaw","scale*pow(E,-gamma)",RooArgSet(E,gamma,scale));
    RooFormulaVar pwl("pwl","(scale*E)**(-gamma)",RooArgList(E,gamma,scale));
    //ws->import(pwl);
    RooGenericPdf pwl_pdf("pwl_pdf","pwl",RooArgSet(pwl));
    RooAddPdf bmodel("bmodel","bmodel",
                       RooArgList(pwl_pdf),
                       RooArgList(norm));
    ws->import(bmodel);
    if (!silent) ws->Print();
}
void SlidingWindowFit::toyMC(char signal_pdf[64], int ntoys){
    gROOT->SetBatch(true);
    RooRandom::randomGenerator()->SetSeed(seed);
    // use pdf to generate events.
    silent = true;
    savePlot = false;
    RooRealVar *E = ws->var("E");
    int n;
    char buffer[128];
    n = sprintf(buffer, "E > %1.4f && E <= %1.4f", emin, emax);
    r_data = (RooDataSet *) ws->data("data")->reduce(RooArgSet(*E), buffer);
    n = sprintf(buffer,"h1_window%d_emin_%d_emax_%d", iter, int(emin), int(emax));
    int nobs = r_data->sumEntries();

    RooAbsPdf *pdf = ws->pdf(signal_pdf);
    double g = ws->var("gamma")->getVal();
    ws->var("gamma")->setVal(3.10);
    fOutput->cd("toyMC");
    TH1D *h1_window = new TH1D("h1_window","toyMC",100,1.8,4.4);
    TH1D *h1_chi2 = new TH1D("h1_chi2","toyMC chi2 log10",60,-1,1);
    // loop over all toys.
    for (int i = 0; i < ntoys; i++){
        int n_obs_mock = gRandom->Poisson(nobs);
        r_data = pdf->generate(*E,n_obs_mock);
        this->fit(true);
        h1_window->Fill(index[0]);
        h1_chi2->Fill(TMath::Log10(chi2));
    }
    h1_window->SetName(buffer);
    n = sprintf(buffer,"h1_chi2_%d_emin_%d_emax_%d", iter, int(emin), int(emax));
    h1_chi2->SetName(buffer);
    h1_window->Write();
    h1_chi2->Write();
    fOutput->cd("");
    delete h1_window;
    delete h1_chi2;
}

void SlidingWindowFit::fit(bool doToyMC){
    int n;
    char buffer[128];
    n = sprintf(buffer, "E > %1.4f && E <= %1.4f", emin, emax);
    if (!silent) std::cout << "cut: " << buffer << std::endl;

    RooRealVar *E = ws->var("E");
    E->setMax(emax);
    E->setMin(emin);
    E->setVal(ecenter);

    RooRealVar *gamma = ws->var("gamma");
    RooRealVar *scale = ws->var("scale");
    RooRealVar *norm = ws->var("norm");
    gamma->setConstant(false);
    scale->setConstant(false);
    //RooRealVar *w=ws->var("weight");
    if (!doToyMC) {
        r_data = (RooDataSet *) ws->data("data")->reduce(RooArgSet(*E), buffer);
    }

    if (!silent) r_data->Print();
    int nobs = r_data->sumEntries();
    nentries = nobs;
    if (!nobs){
        std::cout << " found null entries in window, do nothing. " << std::endl;
        return;
    }
    float f_nobs = float(nobs);
    norm->setVal(f_nobs);
    norm->setMin(f_nobs - TMath::Sqrt(f_nobs));
    norm->setMax(f_nobs + TMath::Sqrt(f_nobs));
    //norm->setConstant(true);
    RooAbsPdf *pdf = ws->pdf("bmodel");
    if (!silent) pdf->Print();
    RooAbsReal* nll = pdf->createNLL(*r_data,NumCPU(4)) ;
    RooMinimizer *minuit= new RooMinimizer(*nll);
    minuit->setMinimizerType("Minuit");
    minuit->setVerbose(false);
    minuit->setPrintLevel(-1);
    RooFitResult *r = new RooFitResult();
    int calls, max_calls, status;
    calls = 1;
    max_calls = 6;
    double nm, gam;
    double min_norm_old, max_norm_old, min_gamma_old, max_gamma_old;
    min_norm_old = norm->getMin();
    max_norm_old = norm->getMax();
    min_gamma_old= gamma->getMin();
    max_gamma_old= gamma->getMax();
    status = -1;
    while (status != 0) {
        if (calls > max_calls) {
            std::cout << "could not find good fit for MINUIT, giving up" << std::endl;
            norm->setMin(min_norm_old);
            norm->setMax(max_norm_old);
            gamma->setMax(max_gamma_old);
            gamma->setMin(min_gamma_old);
            break;
        }
        if (calls > 1){
            std::cout << "fitting - attempt: " << calls << std::endl;
            // if it's pegging, increase
            nm = norm->getVal();
            gam = gamma->getVal();
            if (nm == norm->getMax() || nm == norm->getMin()){
                norm->setVal(norm->getMax());
                norm->setMin(.1 * norm->getMin());
                norm->setMax(10. * norm->getMax());
            }
            if (gam == gamma->getMax() || gam == norm->getMin()){
                gamma->setVal(gamma->getMax());
                gamma->setMin(1.5 * gamma->getMin());
                gamma->setMax(1.5 * gamma->getMax());
            }
        }
        if (calls > 2){
            std::cout << "couldn't find good fit after 3 iterations, changing strategy!" << std::endl;
            norm->setConstant(true);
            minuit->simplex();
            minuit->migrad();
            norm->setConstant(false);
            gamma->setConstant(true);
            minuit->simplex();
            minuit->migrad();
            gamma->setConstant(false);
            minuit->migrad();
        }
        else minuit->migrad();
        r = minuit->save();
        status = r->status();
        //minuit->minos(*gamma);
        calls++;
    }
    minuit->hesse();
    r = minuit->save();

    //RooFitResult *r = pdf->fitTo(*r_data, RooFit::Save(true),RooFit::Minimizer("Minuit", "Migrad"),
    //                             PrintLevel(-1),Verbose(false),Extended(true));
    //RooArgSet* floatPars = pdf->getParameters(*r_data)->selectByAttribute("Constant",kFALSE)
    int ndof = r->floatParsFinal().getSize();
    if (!silent) r->Print();
    index[0] = gamma->getVal();
    //index[1] = gamma->getPropagatedError(*r);
    //index[1] = gamma->getAsymErrorLo();
    //index[2] = gamma->getAsymErrorHi();
    index[1] = -gamma->getPropagatedError(*r);
    index[2] = gamma->getPropagatedError(*r);
    chi2 = 0.;
    //r->Print();
    TCanvas *c = new TCanvas("c",buffer,800,600);
    TPad *pad1 = new TPad("pad1","pad 60%", 0., 0.4, 1., 1.);
    TPad *pad2 = new TPad("pad2","pad 40%", 0., 0., 1., .4);
    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.);
    pad2->SetTopMargin(0.);
    pad2->SetBottomMargin(0.15);
    pad1->SetGridx();
    pad1->SetGridy();
    pad2->SetGridx();
    pad2->SetGridy();
    pad1->Draw();
    pad2->Draw();
    // let's add some plotting.
    pad1->cd();
    RooPlot *plot = E->frame(emin,emax,30);//int(TMath::Sqrt(nobs))));
    pdf->paramOn(plot,Format("NEU",AutoPrecision(1)));
    r_data->plotOn(plot,DataError(RooAbsData::SumW2));

    n = sprintf(buffer,"E > %1.4f && E <= %1.4f (nObs=%1.1e)", emin, emax, float(nobs));

    plot->SetTitle(buffer);
    plot->GetXaxis()->SetTitle("BgoTotalEcorr");
    pdf->plotOn(plot,VisualizeError(*r,1,true),FillColor(kYellow));
    pdf->plotOn(plot,LineStyle(kDashed),LineColor(kBlue),DrawOption("Lsame"));
/*        TCanvas *c = new TCanvas("c",buffer,800,600);
    c->Divide(1,2);
    c->cd(1); */
    //plot->SetAxisRange(emin,emax,"X");
    plot->Draw("same");
    //gPad->SetLogx();
    //gPad->SetLogy();
    plot->GetXaxis()->SetMoreLogLabels();
    plot->GetXaxis()->SetNoExponent();
    plot->GetYaxis()->SetMoreLogLabels();
    plot->GetYaxis()->SetNoExponent();
    if (savePlot) fOutput->cd("fit_panels");
    n = sprintf(buffer,"pwl_fit_window%d_emin_%d_emax_%d", iter, int(emin), int(emax));
    pad2->cd();
    //c->cd(2);
    RooHist *h1 = (RooHist*)plot->residHist();
    h1->GetXaxis()->SetRangeUser(emin,emax);
    h1->Draw();
    //c->Update();
    c->SetName(buffer);
    if (savePlot) c->Write();
    c->Close();
    chi2 = plot->chiSquare(ndof);
    if (!silent) std::cout << "E: " << ecenter << " gamma: " << index[0] << "" << index[1] << "+" << index[2] << " chi2/ndof: " << chi2 << std::endl;
/*    delete c;
    //delete plot;
    delete pad1;
    delete pad2;
    delete minuit;
    delete r;*/
    glob->cd();
}
