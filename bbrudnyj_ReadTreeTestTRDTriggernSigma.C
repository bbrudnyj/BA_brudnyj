/*--------------------------------------------------------------------------------------------------*\
|                                                                                                    |
|    Macro that makes histograms with branches of "nmartin_TestTriggerTRDTree.root"                  |
|                                                                                                    |
|    - TPCnSigma<seperated particle> vs. Rigidity & TOFnSigma<seperated particle> vs. Rigidity       |
|      - without preselection and with preselection of TPC, TOF, TPC && TOF                          |
|   ( <seperated particle> = Electron, Pion, Kaon, Proton, Deuteron, Triton, Helium3, Alpha )        |
|                                                                                                    |
| Author: Benjamin Brudnyj (2016)                                                                    | 
\*--------------------------------------------------------------------------------------------------*/
void bbrudnyj_ReadTreeTestTRDTriggernSigma() {

  static const Int_t nFiles = 2;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Define Variables
  ///__________________
  Int_t Sigma      = 3;

  //char *Collisions = "PbPb_11h_tracklets";
  char *Collisions = "PbPb_11h"; char *Tracklets = "";
  //char *Collisions = "PbPb_11hStd";
  //char *Collisions = "pPb_13b";
  //char *Collisions = "pPb_13c";
  //char *Collisions = "pp_15f";

  //char *Tracklets                  = "4trkl"; Double_t fnTrkl = 4; // exists only in PbPb_11h_tracklets
  //char *Tracklets                  = "5trkl"; Double_t fnTrkl = 5;
  //char *Tracklets                  = "6trkl"; Double_t fnTrkl = 6;

  //Double_t pTCut                   = 2.; char *pT2 = "pT2";
  Double_t pTCut                   = 1000.; char *pT2 = "";  // equivalent to no pT Cut

  char *File[nFiles]                  = {"pdf","root"};

  if(Sigma == 2) {
    Double_t ClearCutDeuteronTPC    = 1.5; // Rigidity cuts on TPCnSigma histograms
    Double_t ClearCutDeuteronTPCTOF = 1.7;
    Double_t ClearCutTritonTPC      = 1.8;
    Double_t ClearCutTritonTPCTOF   = 2.;
  }
  else if(Sigma == 3) {
    Double_t ClearCutDeuteronTPC    = 1.4;
    Double_t ClearCutDeuteronTPCTOF = 1.7;
    Double_t ClearCutTritonTPC      = 1.7;
    Double_t ClearCutTritonTPCTOF   = 2.;
  }
  else {
    Double_t ClearCutDeuteronTPC    = 20.;
    Double_t ClearCutDeuteronTPCTOF = 20.;
    Double_t ClearCutTritonTPC      = 20.;
    Double_t ClearCutTritonTPCTOF   = 20.;
  }


  printf(Form("\nCreating 10 canvases with Sigma = %i from %s\n",Sigma,Collisions));
  if(Tracklets != "" && pTCut != 2.) printf(Form("with %.f tracklet tracks\n\n",fnTrkl));
  else if(Tracklets != "" && pTCut == 2.) printf(Form("with %.f tracklet tracks and pT < %.f\n\n",fnTrkl,pTCut));
  else if(Tracklets == "" && pTCut == 2.) printf(Form("and pT < %.f\n\n",pTCut));

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Get data from file /lustre/nyx/alice/users/bbrudnyj/trunk/nmartin_nuclei/AliAnalysisTaskTestTriggerTRD.cxx
  /// Get histograms from /lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger3.C
  ///_________________________________________________________________________________________________
  TFile *cinput = TFile::Open(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/bbrudnyj_TestTRDTrigger3%s%iSigma%s%s.root",Collisions,Tracklets,Collisions,Sigma,Tracklets,pT2));
  TList *inlist = (TList*)cinput->Get(Form("bbrudnyj_TestTRDTrigger3%s%iSigma%s%s;1",Collisions,Sigma,Tracklets,pT2));

  // TPCnSigma
  TH2D *fHistTPCnSigmavRigvElectron             = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvElectron");
  TH2D *fHistTPCnSigmavRigvPion                 = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvPion");
  TH2D *fHistTPCnSigmavRigvKaon                 = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvKaon");
  TH2D *fHistTPCnSigmavRigvProton               = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvProton");
  TH2D *fHistTPCnSigmavRigvDeuteron             = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvDeuteron");
  TH2D *fHistTPCnSigmavRigvTriton               = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvTriton");
  TH2D *fHistTPCnSigmavRigvHelium3              = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvHelium3");
  TH2D *fHistTPCnSigmavRigvAlpha                = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvAlpha");

  TH2D *fHistTPCnSigmavRigvElectronTPC          = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvElectronTPC");
  TH2D *fHistTPCnSigmavRigvPionTPC              = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvPionTPC");
  TH2D *fHistTPCnSigmavRigvKaonTPC              = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvKaonTPC");
  TH2D *fHistTPCnSigmavRigvProtonTPC            = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvProtonTPC");
  TH2D *fHistTPCnSigmavRigvDeuteronTPC          = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvDeuteronTPC");
  TH2D *fHistTPCnSigmavRigvDeuteronTPC_clear    = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvDeuteronTPC_clear");
  TH2D *fHistTPCnSigmavRigvTritonTPC            = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvTritonTPC");
  TH2D *fHistTPCnSigmavRigvTritonTPC_clear      = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvTritonTPC_clear");
  TH2D *fHistTPCnSigmavRigvHelium3TPC           = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvHelium3TPC");
  TH2D *fHistTPCnSigmavRigvAlphaTPC             = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvAlphaTPC");

  TH2D *fHistTPCnSigmavRigvElectronTOF          = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvElectronTOF");
  TH2D *fHistTPCnSigmavRigvPionTOF              = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvPionTOF");
  TH2D *fHistTPCnSigmavRigvKaonTOF              = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvKaonTOF");
  TH2D *fHistTPCnSigmavRigvProtonTOF            = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvProtonTOF");
  TH2D *fHistTPCnSigmavRigvDeuteronTOF          = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvDeuteronTOF");
  TH2D *fHistTPCnSigmavRigvTritonTOF            = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvTritonTOF");
  TH2D *fHistTPCnSigmavRigvHelium3TOF           = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvHelium3TOF");
  TH2D *fHistTPCnSigmavRigvAlphaTOF             = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvAlphaTOF");

  TH2D *fHistTPCnSigmavRigvElectronTPCTOF       = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvElectronTPCTOF");
  TH2D *fHistTPCnSigmavRigvPionTPCTOF           = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvPionTPCTOF");
  TH2D *fHistTPCnSigmavRigvKaonTPCTOF           = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvKaonTPCTOF");
  TH2D *fHistTPCnSigmavRigvProtonTPCTOF         = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvProtonTPCTOF");
  TH2D *fHistTPCnSigmavRigvDeuteronTPCTOF       = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvDeuteronTPCTOF");
  TH2D *fHistTPCnSigmavRigvDeuteronTPCTOF_clear = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvDeuteronTPCTOF_clear");
  TH2D *fHistTPCnSigmavRigvTritonTPCTOF         = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvTritonTPCTOF");
  TH2D *fHistTPCnSigmavRigvTritonTPCTOF_clear   = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvTritonTPCTOF_clear");
  TH2D *fHistTPCnSigmavRigvHelium3TPCTOF        = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvHelium3TPCTOF");
  TH2D *fHistTPCnSigmavRigvAlphaTPCTOF          = (TH2D*)inlist->FindObject("fHistTPCnSigmavRigvAlphaTPCTOF");

  // TOFnSigma
  TH2D *fHistTOFnSigmavRigvElectron             = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvElectron");
  TH2D *fHistTOFnSigmavRigvPion                 = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvPion");
  TH2D *fHistTOFnSigmavRigvKaon                 = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvKaon");
  TH2D *fHistTOFnSigmavRigvProton               = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvProton");
  TH2D *fHistTOFnSigmavRigvDeuteron             = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvDeuteron");
  TH2D *fHistTOFnSigmavRigvTriton               = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvTriton");
  TH2D *fHistTOFnSigmavRigvHelium3              = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvHelium3");
  TH2D *fHistTOFnSigmavRigvAlpha                = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvAlpha");

  TH2D *fHistTOFnSigmavRigvElectronTPC          = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvElectronTPC");
  TH2D *fHistTOFnSigmavRigvPionTPC              = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvPionTPC");
  TH2D *fHistTOFnSigmavRigvKaonTPC              = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvKaonTPC");
  TH2D *fHistTOFnSigmavRigvProtonTPC            = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvProtonTPC");
  TH2D *fHistTOFnSigmavRigvDeuteronTPC          = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvDeuteronTPC");
  TH2D *fHistTOFnSigmavRigvDeuteronTPC_clear    = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvDeuteronTPC_clear");
  TH2D *fHistTOFnSigmavRigvTritonTPC            = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvTritonTPC");
  TH2D *fHistTOFnSigmavRigvTritonTPC_clear      = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvTritonTPC_clear");
  TH2D *fHistTOFnSigmavRigvHelium3TPC           = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvHelium3TPC");
  TH2D *fHistTOFnSigmavRigvAlphaTPC             = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvAlphaTPC");

  TH2D *fHistTOFnSigmavRigvElectronTOF          = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvElectronTOF");
  TH2D *fHistTOFnSigmavRigvPionTOF              = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvPionTOF");
  TH2D *fHistTOFnSigmavRigvKaonTOF              = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvKaonTOF");
  TH2D *fHistTOFnSigmavRigvProtonTOF            = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvProtonTOF");
  TH2D *fHistTOFnSigmavRigvDeuteronTOF          = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvDeuteronTOF");
  TH2D *fHistTOFnSigmavRigvTritonTOF            = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvTritonTOF");
  TH2D *fHistTOFnSigmavRigvHelium3TOF           = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvHelium3TOF");
  TH2D *fHistTOFnSigmavRigvAlphaTOF             = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvAlphaTOF");

  TH2D *fHistTOFnSigmavRigvElectronTPCTOF       = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvElectronTPCTOF");
  TH2D *fHistTOFnSigmavRigvPionTPCTOF           = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvPionTPCTOF");
  TH2D *fHistTOFnSigmavRigvKaonTPCTOF           = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvKaonTPCTOF");
  TH2D *fHistTOFnSigmavRigvProtonTPCTOF         = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvProtonTPCTOF");
  TH2D *fHistTOFnSigmavRigvDeuteronTPCTOF       = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvDeuteronTPCTOF");
  TH2D *fHistTOFnSigmavRigvDeuteronTPCTOF_clear = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvDeuteronTPCTOF_clear");
  TH2D *fHistTOFnSigmavRigvTritonTPCTOF         = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvTritonTPCTOF");
  TH2D *fHistTOFnSigmavRigvTritonTPCTOF_clear   = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvTritonTPCTOF_clear");
  TH2D *fHistTOFnSigmavRigvHelium3TPCTOF        = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvHelium3TPCTOF");
  TH2D *fHistTOFnSigmavRigvAlphaTPCTOF          = (TH2D*)inlist->FindObject("fHistTOFnSigmavRigvAlphaTPCTOF");


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create legends
  ///__________________________________________
  TLegend *legTPC = new TLegend(.58,.15,.9,.25);
  legTPC->SetFillStyle(0);
  legTPC->SetBorderSize(0);
  legTPC->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  TLegend *legTPCwCutDeuteron = new TLegend(.6,.15,.9,.3);
  legTPCwCutDeuteron->SetFillStyle(0);
  legTPCwCutDeuteron->SetBorderSize(0);
  legTPCwCutDeuteron->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legTPCwCutDeuteron->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
  TLegend *legTPCwCutTriton = new TLegend(.6,.15,.9,.3);
  legTPCwCutTriton->SetFillStyle(0);
  legTPCwCutTriton->SetBorderSize(0);
  legTPCwCutTriton->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legTPCwCutTriton->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");

  TLegend *legTOF = new TLegend(.58,.15,.9,.25);
  legTOF->SetFillStyle(0);
  legTOF->SetBorderSize(0);
  legTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legTOF2 = new TLegend(.6,.55,.9,.65);
  legTOF2->SetFillStyle(0);
  legTOF2->SetBorderSize(0);
  legTOF2->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");

  TLegend *legTPCTOF = new TLegend(.47,.15,.9,.25);
  legTPCTOF->SetFillStyle(0);
  legTPCTOF->SetBorderSize(0);
  legTPCTOF->AddEntry((TObject*)0,Form("|(TPC&TOF)nSigma| < %i#sigma",Sigma),"");
  TLegend *legTPCTOFwCutDeuteron = new TLegend(.45,.15,.88,.3);
  legTPCTOFwCutDeuteron->SetFillStyle(0);
  legTPCTOFwCutDeuteron->SetBorderSize(0);
  legTPCTOFwCutDeuteron->AddEntry((TObject*)0,Form("|(TPC&TOF)nSigma| < %i#sigma",Sigma),"");
  legTPCTOFwCutDeuteron->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
  TLegend *legTPCTOFwCutTriton = new TLegend(.45,.15,.88,.3);
  legTPCTOFwCutTriton->SetFillStyle(0);
  legTPCTOFwCutTriton->SetBorderSize(0);
  legTPCTOFwCutTriton->AddEntry((TObject*)0,Form("|(TPC&TOF)nSigma| < %i#sigma",Sigma),"");
  legTPCTOFwCutTriton->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");

    TLegend *legTracklets = new TLegend(.57,.52,.95,.6);
  legTracklets->SetFillStyle(0);
  legTracklets->SetBorderSize(0);
  if(Tracklets != "") legTracklets->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");

  TLegend *legpTCut = new TLegend(.57,.6,.95,.66);
  legpTCut->SetFillStyle(0);
  legpTCut->SetBorderSize(0);
  if(pTCut == 2.) legpTCut->AddEntry((TObject*)0,Form("#it{p}_{T} < %.f GeV/#it{c}",pTCut),"");

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create and fill canvases
  ///______________________
  // TPCnSigma vs. Rigidity
    TCanvas *c10 = new TCanvas("c10","TPCnSigma vs. Rigidity without preselection for separated particles");
    c10->Divide(3,3);
    c10->cd(2)->SetLogz();
      fHistTPCnSigmavRigvElectron->Draw("colz");
    c10->cd(3)->SetLogz();
      fHistTPCnSigmavRigvPion->Draw("colz");
    c10->cd(4)->SetLogz();
      fHistTPCnSigmavRigvKaon->Draw("colz");
    c10->cd(5)->SetLogz();
      fHistTPCnSigmavRigvProton->Draw("colz");
    c10->cd(6)->SetLogz();
      fHistTPCnSigmavRigvDeuteron->Draw("colz");
    c10->cd(7)->SetLogz();
      fHistTPCnSigmavRigvTriton->Draw("colz");
    c10->cd(8)->SetLogz();
      fHistTPCnSigmavRigvHelium3->Draw("colz");
    c10->cd(9)->SetLogz();
      fHistTPCnSigmavRigvAlpha->Draw("colz");
    for(Int_t i = 2; i <= 9; i++) {
      c10->cd(i)->SetLogx();
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

    TCanvas *c11 = new TCanvas("c11",Form("TPCnSigma vs. Rigidity w TPC preselection (Sigma=%i) for separated particles",Sigma));
    c11->Divide(3,3);
    c11->cd(2)->SetLogz();
      fHistTPCnSigmavRigvElectronTPC->Draw("colz");
      legTPC->Draw("same");
    c11->cd(3)->SetLogz();
      fHistTPCnSigmavRigvPionTPC->Draw("colz");
      legTPC->Draw("same");
    c11->cd(4)->SetLogz();
      fHistTPCnSigmavRigvKaonTPC->Draw("colz");
      legTPC->Draw("same");
    c11->cd(5)->SetLogz();
      fHistTPCnSigmavRigvProtonTPC->Draw("colz");
      legTPC->Draw("same");
    c11->cd(6)->SetLogz();
      fHistTPCnSigmavRigvDeuteronTPC->Draw("colz");
      legTPC->Draw("same");
    c11->cd(7)->SetLogz();
      fHistTPCnSigmavRigvTritonTPC->Draw("colz");
      legTPC->Draw("same");
    c11->cd(8)->SetLogz();
      fHistTPCnSigmavRigvHelium3TPC->Draw("colz");
      legTPC->Draw("same");
    c11->cd(9)->SetLogz();
      fHistTPCnSigmavRigvAlphaTPC->Draw("colz");
      legTPC->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c11->cd(i)->SetLogx();
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

    TCanvas *c12 = new TCanvas("c12",Form("TPCnSigma vs. Rigidity (with Cut)  w TPC preselection (Sigma=%i) for separated particles",Sigma));
    c12->Divide(3,3);
    c12->cd(6)->SetLogz();
      fHistTPCnSigmavRigvDeuteronTPC_clear->Draw("colz");
      legTPCwCutDeuteron->Draw("same");
    c12->cd(7)->SetLogz();
      fHistTPCnSigmavRigvTritonTPC_clear->Draw("colz");
      legTPCwCutTriton->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c12->cd(i)->SetLogx();
      if(i==6 || i==7) {
        legTracklets->Draw("same");
        legpTCut->Draw("same");
      }
    }

    TCanvas *c13 = new TCanvas("c13",Form("TPCnSigma vs. Rigidity w TOF preselection (Sigma=%i) for separated particles",Sigma));
    c13->Divide(3,3);
    c13->cd(2)->SetLogz();
      fHistTPCnSigmavRigvElectronTOF->Draw("colz");
      legTOF->Draw("same");
    c13->cd(3)->SetLogz();
      fHistTPCnSigmavRigvPionTOF->Draw("colz");
      legTOF->Draw("same");
    c13->cd(4)->SetLogz();
      fHistTPCnSigmavRigvKaonTOF->Draw("colz");
      legTOF->Draw("same");
    c13->cd(5)->SetLogz();
      fHistTPCnSigmavRigvProtonTOF->Draw("colz");
      legTOF->Draw("same");
    c13->cd(6)->SetLogz();
      fHistTPCnSigmavRigvDeuteronTOF->Draw("colz");
      legTOF->Draw("same");
    c13->cd(7)->SetLogz();
      fHistTPCnSigmavRigvTritonTOF->Draw("colz");
      legTOF->Draw("same");
    c13->cd(8)->SetLogz();
      fHistTPCnSigmavRigvHelium3TOF->Draw("colz");
      legTOF->Draw("same");
    c13->cd(9)->SetLogz();
      fHistTPCnSigmavRigvAlphaTOF->Draw("colz");
      legTOF->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c13->cd(i)->SetLogx();
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

    TCanvas *c14 = new TCanvas("c14",Form("TPCnSigma vs. Rigidity w TPC & TOF preselection (Sigma=%i) for separated particles",Sigma));
    c14->Divide(3,3);
    c14->cd(2)->SetLogz();
      fHistTPCnSigmavRigvElectronTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c14->cd(3)->SetLogz();
      fHistTPCnSigmavRigvPionTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c14->cd(4)->SetLogz();
      fHistTPCnSigmavRigvKaonTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c14->cd(5)->SetLogz();
      fHistTPCnSigmavRigvProtonTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c14->cd(6)->SetLogz();
      fHistTPCnSigmavRigvDeuteronTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c14->cd(7)->SetLogz();
      fHistTPCnSigmavRigvTritonTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c14->cd(8)->SetLogz();
      fHistTPCnSigmavRigvHelium3TPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c14->cd(9)->SetLogz();
      fHistTPCnSigmavRigvAlphaTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c14->cd(i)->SetLogx();
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

    TCanvas *c15 = new TCanvas("c15",Form("TPCnSigma vs. Rigidity (w Cut)  w TPC & TOF preselection (Sigma=%i) for separated particles",Sigma));
    c15->Divide(3,3);
    c15->cd(6)->SetLogz();
      fHistTPCnSigmavRigvDeuteronTPCTOF_clear->Draw("colz");
      legTPCTOFwCutDeuteron->Draw("same");
    c15->cd(7)->SetLogz();
      fHistTPCnSigmavRigvTritonTPCTOF_clear->Draw("colz");
      legTPCTOFwCutTriton->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c15->cd(i)->SetLogx();
      if(i==6 || i==7) {
        legTracklets->Draw("same");
        legpTCut->Draw("same");
      }
    }


  // TOFnSigma vs. Rigidity
    TCanvas *c16 = new TCanvas("c16","TOFnSigma vs. Rigidity without preselection for separated particles");
    c16->Divide(3,3);
    c16->cd(2)->SetLogz();
      fHistTOFnSigmavRigvElectron->Draw("colz");
    c16->cd(3)->SetLogz();
      fHistTOFnSigmavRigvPion->Draw("colz");
    c16->cd(4)->SetLogz();
      fHistTOFnSigmavRigvKaon->Draw("colz");
    c16->cd(5)->SetLogz();
      fHistTOFnSigmavRigvProton->Draw("colz");
    c16->cd(6)->SetLogz();
      fHistTOFnSigmavRigvDeuteron->Draw("colz");
    c16->cd(7)->SetLogz();
      fHistTOFnSigmavRigvTriton->Draw("colz");
    c16->cd(8)->SetLogz();
      fHistTOFnSigmavRigvHelium3->Draw("colz");
    c16->cd(9)->SetLogz();
      fHistTOFnSigmavRigvAlpha->Draw("colz");
    for(Int_t i = 2; i <= 9; i++) {
      c16->cd(i)->SetLogx();
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

    TCanvas *c17 = new TCanvas("c17",Form("TOFnSigma vs. Rigidity w TPC preselection (Sigma=%i) for separated particles",Sigma));
    c17->Divide(3,3);
    c17->cd(2)->SetLogz();
      fHistTOFnSigmavRigvElectronTPC->Draw("colz");
      legTPC->Draw("same");
    c17->cd(3)->SetLogz();
      fHistTOFnSigmavRigvPionTPC->Draw("colz");
      legTPC->Draw("same");
    c17->cd(4)->SetLogz();
      fHistTOFnSigmavRigvKaonTPC->Draw("colz");
      legTPC->Draw("same");
    c17->cd(5)->SetLogz();
      fHistTOFnSigmavRigvProtonTPC->Draw("colz");
      legTPC->Draw("same");
    c17->cd(6)->SetLogz();
      fHistTOFnSigmavRigvDeuteronTPC->Draw("colz");
      legTPC->Draw("same");
    c17->cd(7)->SetLogz();
      fHistTOFnSigmavRigvTritonTPC->Draw("colz");
      legTPC->Draw("same");
    c17->cd(8)->SetLogz();
      fHistTOFnSigmavRigvHelium3TPC->Draw("colz");
      legTPC->Draw("same");
    c17->cd(9)->SetLogz();
      fHistTOFnSigmavRigvAlphaTPC->Draw("colz");
      legTPC->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c17->cd(i)->SetLogx();
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

    TCanvas *c18 = new TCanvas("c18",Form("TOFnSigma vs. Rigidity w TOF preselection (Sigma=%i) for separated particles",Sigma));
    c18->Divide(3,3);
    c18->cd(2)->SetLogz();
      fHistTOFnSigmavRigvElectronTOF->Draw("colz");
      legTOF->Draw("same");
    c18->cd(3)->SetLogz();
      fHistTOFnSigmavRigvPionTOF->Draw("colz");
      legTOF->Draw("same");
    c18->cd(4)->SetLogz();
      fHistTOFnSigmavRigvKaonTOF->Draw("colz");
      legTOF->Draw("same");
    c18->cd(5)->SetLogz();
      fHistTOFnSigmavRigvProtonTOF->Draw("colz");
      legTOF->Draw("same");
    c18->cd(6)->SetLogz();
      fHistTOFnSigmavRigvDeuteronTOF->Draw("colz");
      legTOF->Draw("same");
    c18->cd(7)->SetLogz();
      fHistTOFnSigmavRigvTritonTOF->Draw("colz");
      legTOF->Draw("same");
    c18->cd(8)->SetLogz();
      fHistTOFnSigmavRigvHelium3TOF->Draw("colz");
      legTOF->Draw("same");
    c18->cd(9)->SetLogz();
      fHistTOFnSigmavRigvAlphaTOF->Draw("colz");
      legTOF->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c18->cd(i)->SetLogx();
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

    TCanvas *c19 = new TCanvas("c19",Form("TOFnSigma vs. Rigidity w TPC & TOF preselection (Sigma=%i) for separated particles",Sigma));
    c19->Divide(3,3);
    c19->cd(2)->SetLogz();
      fHistTOFnSigmavRigvElectronTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c19->cd(3)->SetLogz();
      fHistTOFnSigmavRigvPionTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c19->cd(4)->SetLogz();
      fHistTOFnSigmavRigvKaonTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c19->cd(5)->SetLogz();
      fHistTOFnSigmavRigvProtonTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c19->cd(6)->SetLogz();
      fHistTOFnSigmavRigvDeuteronTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c19->cd(7)->SetLogz();
      fHistTOFnSigmavRigvTritonTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c19->cd(8)->SetLogz();
      fHistTOFnSigmavRigvHelium3TPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    c19->cd(9)->SetLogz();
      fHistTOFnSigmavRigvAlphaTPCTOF->Draw("colz");
      legTPCTOF->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c19->cd(i)->SetLogx();
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Save canvases
  ///_________________________
  for(Int_t f = 0; f < nFiles; f++) {

  c10->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TPCnSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Collisions,Tracklets,pT2,File[f]));
  c11->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TPCnSigmaTPC%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Sigma,Collisions,Tracklets,pT2,File[f]));
  c12->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TPCnSigmaTPC%iSigma%s_clear%s%s.%s",Collisions,Tracklets,Sigma,pT2,Sigma,Collisions,Tracklets,pT2,File[f]));
  c13->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TPCnSigmaTOF%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Sigma,Collisions,Tracklets,pT2,File[f]));
  c14->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TPCnSigmaTPCTOF%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Sigma,Collisions,Tracklets,pT2,File[f]));
  c15->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TPCnSigmaTPCTOF%iSigma%s_clear%s%s.%s",Collisions,Tracklets,Sigma,pT2,Sigma,Collisions,Tracklets,pT2,File[f]));
  c16->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TOFnSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Collisions,Tracklets,pT2,File[f]));
  c17->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TOFnSigmaTPC%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Sigma,Collisions,Tracklets,pT2,File[f]));
  c18->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TOFnSigmaTOF%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Sigma,Collisions,Tracklets,pT2,File[f]));
  c19->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/plots%iSigma%s/bbrudnyj_TOFnSigmaTPCTOF%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Sigma,Collisions,Tracklets,pT2,File[f]));

  }

  printf("\ndone!\n\n");
}
