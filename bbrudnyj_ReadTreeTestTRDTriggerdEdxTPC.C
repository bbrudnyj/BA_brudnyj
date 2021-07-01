 /*--------------------------------------------------------------------------------------------------*\
|                                                                                                    |
|    Macro that makes projection histograms with branches of "nmartin_TestTriggerTRDTree.root"       |
|                                                                                                    |
|    - TPC dE/dx vs. rigidity for: all particles                                                     |
|      - with preselection of TPC, TOF, TPC && TOF for Electrons, Pions, Kaons, Protons, Deuterons,  |
|                                                      Tritons, Helium3, Alphas                      |
|                                                                                                    |
| Author: Benjamin Brudnyj (2016)                                                                    | 
\*--------------------------------------------------------------------------------------------------*/
void bbrudnyj_ReadTreeTestTRDTriggerdEdxTPC() {

  static const Int_t nFiles = 2;

  gStyle->SetTitleW(.8f);
  gStyle->SetTitleH(.08f);
  gStyle->SetTitleSize(.07,"xy");
  //gStyle->SetTitleOffset(.8,"x");
  gROOT->ForceStyle();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Define Variables
  ///__________________
  Int_t Sigma      = 3; // 2 or 3

  //char *Collisions = "PbPb_11h_tracklets";
  char *Collisions = "PbPb_11h"; char *Tracklets = "";
  //char *Collisions = "PbPb_11hStd";
  //char *Collisions = "pPb_13b";
  //char *Collisions = "pPb_13c";
  //char *Collisions = "pp_15f";

  //char *Tracklets                   = "4trkl"; Double_t fnTrkl = 4; // exists only in PbPb_11h_tracklets
  //char *Tracklets                   = "5trkl"; Double_t fnTrkl = 5;
  //char *Tracklets                   = "6trkl"; Double_t fnTrkl = 6;

  //Double_t pTCut                    = 2.; char *pT2 = "pT2";
  Double_t pTCut                    = 1000.; char *pT2 = "";  // equivalent to no pT Cut

  char *File[nFiles]                = {"pdf","root"};

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

  char *Cut = "yes"; char *RigidityCut = "with"; // Rigidity cuts on d & t for a clear signal
  //char *Cut = "no"; char *RigidityCut = "without";

  printf(Form("\nCreating canvases with Sigma = %i from %s %s rigidity-cuts on d & t\n",Sigma,Collisions,RigidityCut));
  if(Tracklets != "" && pTCut != 2.) printf(Form("with %.f tracklet tracks\n\n",fnTrkl));
  else if(Tracklets != "" && pTCut == 2.) printf(Form("with %.f tracklet tracks and pT < %.f\n\n",fnTrkl,pTCut));
  else if(Tracklets == "" && pTCut == 2.) printf(Form("and pT < %.f\n\n",pTCut));


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Get data from file /lustre/nyx/alice/users/bbrudnyj/trunk/nmartin_nuclei/AliAnalysisTaskTestTriggerTRD.cxx
  /// Get histograms from /lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2dEdx.C
  ///_________________________________________________________________________________________________
  TFile *cinput       = TFile::Open(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/bbrudnyj_TestTRDTrigger2%s%iSigma%sdEdx%s.root",Collisions,Tracklets,Collisions,Sigma,Tracklets,pT2));
  TList *inlist       = (TList*)cinput->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigma%sdEdx%s;1",Collisions,Sigma,Tracklets,pT2));
  TList *inlistTPC    = (TList*)cinput->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigma%sdEdxTPC%s;1",Collisions,Sigma,Tracklets,pT2));
  TList *inlistTOF    = (TList*)cinput->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigma%sdEdxTOF%s;1",Collisions,Sigma,Tracklets,pT2));
  TList *inlistTPCTOF = (TList*)cinput->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigma%sdEdxTPCTOF%s;1",Collisions,Sigma,Tracklets,pT2));

  TH2D *fHistdEdxvRig                                 = (TH2D*)inlist->FindObject("fHistdEdxvRig");

  TH2D *fHistdEdxvRigvElectronTPC                     = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvElectronTPC");
  TH2D *fHistdEdxvRigvPionTPC                         = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvPionTPC");
  TH2D *fHistdEdxvRigvKaonTPC                         = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvKaonTPC");
  TH2D *fHistdEdxvRigvProtonTPC                       = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvProtonTPC");
  if(Cut == "no")  TH2D *fHistdEdxvRigvDeuteronTPC    = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvDeuteronTPC");
  if(Cut == "yes") TH2D *fHistdEdxvRigvDeuteronTPC    = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvDeuteronTPC_clear");
  if(Cut == "no")  TH2D *fHistdEdxvRigvTritonTPC      = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvTritonTPC");
  if(Cut == "yes") TH2D *fHistdEdxvRigvTritonTPC      = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvTritonTPC_clear");
  TH2D *fHistdEdxvRigvHelium3TPC                      = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvHelium3TPC");
  TH2D *fHistdEdxvRigvHelium3TPC_bb                   = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvHelium3TPC_bb");
  TH2D *fHistdEdxvRigvAlphaTPC                        = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvAlphaTPC");
  TH2D *fHistdEdxvRigvAlphaTPC_bb                     = (TH2D*)inlistTPC->FindObject("fHistdEdxvRigvAlphaTPC_bb");

  TH2D *fHistdEdxvRigvElectronTOF                     = (TH2D*)inlistTOF->FindObject("fHistdEdxvRigvElectronTOF");
  TH2D *fHistdEdxvRigvPionTOF                         = (TH2D*)inlistTOF->FindObject("fHistdEdxvRigvPionTOF");
  TH2D *fHistdEdxvRigvKaonTOF                         = (TH2D*)inlistTOF->FindObject("fHistdEdxvRigvKaonTOF");
  TH2D *fHistdEdxvRigvProtonTOF                       = (TH2D*)inlistTOF->FindObject("fHistdEdxvRigvProtonTOF");
  TH2D *fHistdEdxvRigvDeuteronTOF                     = (TH2D*)inlistTOF->FindObject("fHistdEdxvRigvDeuteronTOF");
  TH2D *fHistdEdxvRigvTritonTOF                       = (TH2D*)inlistTOF->FindObject("fHistdEdxvRigvTritonTOF");
  TH2D *fHistdEdxvRigvHelium3TOF                      = (TH2D*)inlistTOF->FindObject("fHistdEdxvRigvHelium3TOF");
  TH2D *fHistdEdxvRigvAlphaTOF                        = (TH2D*)inlistTOF->FindObject("fHistdEdxvRigvAlphaTOF");

  TH2D *fHistdEdxvRigvElectronTPCTOF                  = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvElectronTPCTOF");
  TH2D *fHistdEdxvRigvPionTPCTOF                      = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvPionTPCTOF");
  TH2D *fHistdEdxvRigvKaonTPCTOF                      = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvKaonTPCTOF");
  TH2D *fHistdEdxvRigvProtonTPCTOF                    = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvProtonTPCTOF");
  if(Cut == "no")  TH2D *fHistdEdxvRigvDeuteronTPCTOF = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvDeuteronTPCTOF");
  if(Cut == "yes") TH2D *fHistdEdxvRigvDeuteronTPCTOF = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvDeuteronTPCTOF_clear");
  if(Cut == "no")  TH2D *fHistdEdxvRigvTritonTPCTOF   = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvTritonTPCTOF");
  if(Cut == "yes") TH2D *fHistdEdxvRigvTritonTPCTOF   = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvTritonTPCTOF_clear");
  TH2D *fHistdEdxvRigvHelium3TPCTOF                   = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvHelium3TPCTOF");
  TH2D *fHistdEdxvRigvHelium3TPCTOF_bb                = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvHelium3TPCTOF_bb");
  TH2D *fHistdEdxvRigvAlphaTPCTOF                     = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvAlphaTPCTOF");
  TH2D *fHistdEdxvRigvAlphaTPCTOF_bb                  = (TH2D*)inlistTPCTOF->FindObject("fHistdEdxvRigvAlphaTPCTOF_bb");


  fHistdEdxvRig->SetStats(0);
  fHistdEdxvRig->SetTitle("");
  fHistdEdxvRig->GetYaxis()->SetTitle("TPC  dE/dx (arb. units)");
  fHistdEdxvRig->GetXaxis()->SetTitle("#frac{#it{p}}{z} (GeV/#it{c})   ");

  fHistdEdxvRigvHelium3TPC_bb->SetStats(0);
  fHistdEdxvRigvHelium3TPC_bb->SetTitle("");
  fHistdEdxvRigvHelium3TPC_bb->GetYaxis()->SetTitle("TPC  dE/dx (arb. units)");

  fHistdEdxvRigvHelium3TPCTOF_bb->SetStats(0);
  fHistdEdxvRigvHelium3TPCTOF_bb->SetTitle("");
  fHistdEdxvRigvHelium3TPCTOF_bb->GetYaxis()->SetTitle("TPC  dE/dx (arb. units)");
  fHistdEdxvRigvHelium3TPCTOF_bb->GetZaxis()->SetRangeUser(0.5,10000);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create legends
  ///_________________________________________
  TLegend *leg = new TLegend(.05,.63,.42,.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if(Collisions == "PbPb_11h" || Collisions == "PbPb_11hStd" || Collisions == "PbPb_11h_tracklets") {
    leg->AddEntry((TObject*)0,"Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV","");
    //leg->AddEntry((TObject*)0,"run LHC11h pass2","");
    leg->AddEntry((TObject*)0," ","");
  }
  else if(Collisions == "pPb_13b") {
    leg->AddEntry((TObject*)0,"p-Pb, #sqrt{s_{NN}} = 5.023 TeV","");
    leg->AddEntry((TObject*)0,"run LHC13b pass3","");
  }
  else if(Collisions == "pPb_13c") {
    leg->AddEntry((TObject*)0,"p-Pb, #sqrt{s_{NN}} = 5.023 TeV","");
    leg->AddEntry((TObject*)0,"run LHC13c pass2","");
  }
  else { //(Collisions == "pp_15f")
    leg->AddEntry((TObject*)0,"p-p, #sqrt{s} = 13 TeV","");
    leg->AddEntry((TObject*)0,"run LHC15f pass2","");
  }

  TLegend *legTracklets = new TLegend(.57,.52,.95,.6);
  legTracklets->SetFillStyle(0);
  legTracklets->SetBorderSize(0);
  if(Tracklets != "") legTracklets->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");

  TLegend *legpTCut = new TLegend(.57,.6,.95,.66);
  legpTCut->SetFillStyle(0);
  legpTCut->SetBorderSize(0);
  if(pTCut == 2.) legpTCut->AddEntry((TObject*)0,Form("#it{p}_{T} < %.f GeV/#it{c}",pTCut),"");

  TLegend *legTPC = new TLegend(.05,.53,.4,.65);
  legTPC->SetFillStyle(0);
  legTPC->SetBorderSize(0);
  legTPC->AddEntry((TObject*)0,"with preselection:","");
  legTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  TLegend *legTPCwCutDeuteron = new TLegend(.05,.47,.4,.65);
  legTPCwCutDeuteron->SetFillStyle(0);
  legTPCwCutDeuteron->SetBorderSize(0);
  legTPCwCutDeuteron->AddEntry((TObject*)0,"with preselection:","");
  legTPCwCutDeuteron->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legTPCwCutDeuteron->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
  TLegend *legTPCwCutTriton = new TLegend(.05,.47,.4,.65);
  legTPCwCutTriton->SetFillStyle(0);
  legTPCwCutTriton->SetBorderSize(0);
  legTPCwCutTriton->AddEntry((TObject*)0,"with preselection:","");
  legTPCwCutTriton->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legTPCwCutTriton->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");

  TLegend *legTOF = new TLegend(.05,.53,.4,.65);
  legTOF->SetFillStyle(0);
  legTOF->SetBorderSize(0);
  legTOF->AddEntry((TObject*)0,"with preselection:","");
  legTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");

  TLegend *legTPCTOF = new TLegend(.05,.55,.4,.73);
  legTPCTOF->SetFillStyle(0);
  legTPCTOF->SetBorderSize(0);
  legTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  TLegend *legTPCTOFwCutDeuteron = new TLegend(.05,.41,.4,.65);
  legTPCTOFwCutDeuteron->SetFillStyle(0);
  legTPCTOFwCutDeuteron->SetBorderSize(0);
  legTPCTOFwCutDeuteron->AddEntry((TObject*)0,"with preselection:","");
  legTPCTOFwCutDeuteron->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legTPCTOFwCutDeuteron->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  legTPCTOFwCutDeuteron->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
  TLegend *legTPCTOFwCutTriton = new TLegend(.05,.41,.4,.65);
  legTPCTOFwCutTriton->SetFillStyle(0);
  legTPCTOFwCutTriton->SetBorderSize(0);
  legTPCTOFwCutTriton->AddEntry((TObject*)0,"with preselection:","");
  legTPCTOFwCutTriton->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legTPCTOFwCutTriton->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  legTPCTOFwCutTriton->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");


  TLegend *legHe3 = new TLegend(.05,.73,.4,.85);
  legHe3->SetFillStyle(0);
  legHe3->SetBorderSize(0);
  legHe3->AddEntry((TObject*)0,"^{3}He","");


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create and fill canvases
  ///_________________________________________________________________________________________________
  TCanvas *c7 = new TCanvas("c7",Form("dEdx vs. Rigidity w TPC preselection (Sigma=%i) for separated particles",Sigma));
  c7->Divide(3,3);
  c7->cd(1)->SetLogz();
    fHistdEdxvRig->Draw("colz");
    leg->Draw("same");
  c7->cd(2)->SetLogz();
    fHistdEdxvRigvElectronTPC->Draw("colz");
    legTPC->Draw("same");
  c7->cd(3)->SetLogz();
    fHistdEdxvRigvPionTPC->Draw("colz");
    legTPC->Draw("same");
  c7->cd(4)->SetLogz();
    fHistdEdxvRigvKaonTPC->Draw("colz");
    legTPC->Draw("same");
  c7->cd(5)->SetLogz();
    fHistdEdxvRigvProtonTPC->Draw("colz");
    legTPC->Draw("same");
  c7->cd(6)->SetLogz();
    fHistdEdxvRigvDeuteronTPC->Draw("colz");
    if(Cut == "no")  legTPC->Draw("same");
    if(Cut == "yes") legTPCwCutDeuteron->Draw("same");
  c7->cd(7)->SetLogz();
    fHistdEdxvRigvTritonTPC->Draw("colz");
    if(Cut == "no")  legTPC->Draw("same");
    if(Cut == "yes") legTPCwCutTriton->Draw("same");
  c7->cd(8)->SetLogz();
    //fHistdEdxvRigvHelium3TPC->Draw("colz");
    fHistdEdxvRigvHelium3TPC_bb->Draw("colz");
    legTPC->Draw("same");
  c7->cd(9)->SetLogz();
    //fHistdEdxvRigvAlphaTPC->Draw("colz");
    fHistdEdxvRigvAlphaTPC_bb->Draw("colz");
    legTPC->Draw("same");
  for(Int_t i = 1; i <= 9 ; i++) {
    c7->cd(i)->SetTicks();
    //c7->cd(i)->SetLogx();
    legTracklets->Draw("same");
    legpTCut->Draw("same");
  }

  TCanvas *c8 = new TCanvas("c8",Form("dEdx vs. Rigidity w TOF preselection (Sigma=%i) for separated particles",Sigma));
  c8->Divide(3,3);
  c8->cd(1)->SetLogz();
    fHistdEdxvRig->Draw("colz");
    leg->Draw("same");
  c8->cd(2)->SetLogz();
    fHistdEdxvRigvElectronTOF->Draw("colz");
    legTOF->Draw("same");
  c8->cd(3)->SetLogz();
    fHistdEdxvRigvPionTOF->Draw("colz");
    legTOF->Draw("same");
  c8->cd(4)->SetLogz();
    fHistdEdxvRigvKaonTOF->Draw("colz");
    legTOF->Draw("same");
  c8->cd(5)->SetLogz();
    fHistdEdxvRigvProtonTOF->Draw("colz");
    legTOF->Draw("same");
  c8->cd(6)->SetLogz();
    fHistdEdxvRigvDeuteronTOF->Draw("colz");
    legTOF->Draw("same");
  c8->cd(7)->SetLogz();
    fHistdEdxvRigvTritonTOF->Draw("colz");
    legTOF->Draw("same");
  c8->cd(8)->SetLogz();
    fHistdEdxvRigvHelium3TOF->Draw("colz");
    legTOF->Draw("same");
  c8->cd(9)->SetLogz();
    fHistdEdxvRigvAlphaTOF->Draw("colz");
    legTOF->Draw("same");
  for(Int_t i = 1; i <= 9 ; i++) {
    c8->cd(i)->SetTicks();
    //c8->cd(i)->SetLogx();
    legTracklets->Draw("same");
    legpTCut->Draw("same");
  }

  TCanvas *c9 = new TCanvas("c9",Form("dEdx vs. Rigidity w TPC & TOF preselection (Sigma=%i) for separated particles",Sigma));
  c9->Divide(3,3);
  c9->cd(1)->SetLogz();
    fHistdEdxvRig->Draw("colz");
    leg->Draw("same");
  c9->cd(2)->SetLogz();
    fHistdEdxvRigvElectronTPCTOF->Draw("colz");
    legTPCTOF->Draw("same");
  c9->cd(3)->SetLogz();
    fHistdEdxvRigvPionTPCTOF->Draw("colz");
    legTPCTOF->Draw("same");
  c9->cd(4)->SetLogz();
    fHistdEdxvRigvKaonTPCTOF->Draw("colz");
    legTPCTOF->Draw("same");
  c9->cd(5)->SetLogz();
    fHistdEdxvRigvProtonTPCTOF->Draw("colz");
    legTPCTOF->Draw("same");
  c9->cd(6)->SetLogz();
    fHistdEdxvRigvDeuteronTPCTOF->Draw("colz");
    if(Cut == "no")  legTPCTOF->Draw("same");
    if(Cut == "yes") legTPCTOFwCutDeuteron->Draw("same");
  c9->cd(7)->SetLogz();
    fHistdEdxvRigvTritonTPCTOF->Draw("colz");
    if(Cut == "no")  legTPCTOF->Draw("same");
    if(Cut == "yes") legTPCTOFwCutTriton->Draw("same");
  c9->cd(8)->SetLogz();
    //fHistdEdxvRigvHelium3TPCTOF->Draw("colz");
    fHistdEdxvRigvHelium3TPCTOF_bb->Draw("colz");
    legTPCTOF->Draw("same");
  c9->cd(9)->SetLogz();
    //fHistdEdxvRigvAlphaTPCTOF->Draw("colz");
    fHistdEdxvRigvAlphaTPCTOF_bb->Draw("colz");
    legTPCTOF->Draw("same");
  for(Int_t i = 1; i <= 9 ; i++) {
    c9->cd(i)->SetTicks();
    //c9->cd(i)->SetLogx();
    legTracklets->Draw("same");
    legpTCut->Draw("same");
  }

  TCanvas *TPCdEdx = new TCanvas("TPCdEdx",Form(""));
  //TPCdEdx->SetCanvasSize(1600,600);
  TPCdEdx->SetLogz();
    fHistdEdxvRig->Draw("colz");
    leg ->Draw("same");
  TPCdEdx->SaveAs(Form("~/TPCdEdx.jpg"));

  TCanvas *TPCdEdxHe3 = new TCanvas("TPCdEdxHe3",Form(""));
  TPCdEdxHe3->SetLogz();
    fHistdEdxvRigvHelium3TPC_bb->Draw("colz");
    legTPC->Draw("same");
    legHe3->Draw("same");
  TPCdEdxHe3->SaveAs(Form("~/TPCdEdxHe3.pdf"));

  TCanvas *TPCdEdxHe3TOF = new TCanvas("TPCdEdxHe3TOF",Form(""));
  TPCdEdxHe3TOF->SetLogz();
    fHistdEdxvRigvHelium3TPCTOF_bb->Draw("colz");
    legTPCTOF->Draw("same");
    legHe3->Draw("same");
  TPCdEdxHe3TOF->SaveAs(Form("~/TPCdEdxHe3TOF.pdf"));


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Save canvases
  ///_________________________________________________________________________________________________
  for(Int_t f = 0; f < nFiles; f++) {
    c8->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_TPCdEdx%sTOF%iSigma%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,File[f]));
    if(Cut == "yes") {
      c7->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_TPCdEdx%sTPC%iSigma%swC.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,File[f]));
      c9->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_TPCdEdx%sTPCTOF%iSigma%swC.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,File[f]));
    }
    else if(Cut == "no") {
      c7->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_TPCdEdx%sTPC%iSigma%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,File[f]));
      c9->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_TPCdEdx%sTPCTOF%iSigma%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,File[f]));
    }
  }

  printf("\ndone!\n\n");
}
