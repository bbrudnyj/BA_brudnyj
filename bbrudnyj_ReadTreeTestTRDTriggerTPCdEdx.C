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
void bbrudnyj_ReadTreeTestTRDTriggerTPCdEdx() {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Define Variables
  ///__________________
  Int_t Sigma      = 2;

  char *Collisions = "PbPb_11h";
  //char *Collisions = "pPb_13b";
  //char *Collisions = "pPb_13c";
  //char *Collisions = "pp_15f";


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Get data from file /lustre/nyx/alice/users/bbrudnyj/trunk/nmartin_nuclei/AliAnalysisTaskTestTriggerTRD.cxx
  /// Get histograms from /lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2.C
  ///_________________________________________________________________________________________________
  TFile *cinput = TFile::Open(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTrigger2/bbrudnyj_TestTRDTrigger2%s%iSigma.root",Collisions,Collisions,Sigma));
  TList *inlist = (TList*)cinput->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigma;1");

  TH2D *fHistdEdxvRigvElectronTPC      = (TH2D*)inlist->FindObject("fHistdEdxvRigvElectronTPC");
  TH2D *fHistdEdxvRigvPionTPC          = (TH2D*)inlist->FindObject("fHistdEdxvRigvPionTPC");
  TH2D *fHistdEdxvRigvKaonTPC          = (TH2D*)inlist->FindObject("fHistdEdxvRigvKaonTPC");
  TH2D *fHistdEdxvRigvProtonTPC        = (TH2D*)inlist->FindObject("fHistdEdxvRigvProtonTPC");
  TH2D *fHistdEdxvRigvDeuteronTPC      = (TH2D*)inlist->FindObject("fHistdEdxvRigvDeuteronTPC");
  TH2D *fHistdEdxvRigvTritonTPC        = (TH2D*)inlist->FindObject("fHistdEdxvRigvTritonTPC");
  TH2D *fHistdEdxvRigvHelium3TPC       = (TH2D*)inlist->FindObject("fHistdEdxvRigvHelium3TPC");
  TH2D *fHistdEdxvRigvHelium3TPC_bb    = (TH2D*)inlist->FindObject("fHistdEdxvRigvHelium3TPC_bb");
  TH2D *fHistdEdxvRigvAlphaTPC         = (TH2D*)inlist->FindObject("fHistdEdxvRigvAlphaTPC");
  TH2D *fHistdEdxvRigvAlphaTPC_bb      = (TH2D*)inlist->FindObject("fHistdEdxvRigvAlphaTPC_bb");

  TH2D *fHistdEdxvRigvElectronTOF      = (TH2D*)inlist->FindObject("fHistdEdxvRigvElectronTOF");
  TH2D *fHistdEdxvRigvPionTOF          = (TH2D*)inlist->FindObject("fHistdEdxvRigvPionTOF");
  TH2D *fHistdEdxvRigvKaonTOF          = (TH2D*)inlist->FindObject("fHistdEdxvRigvKaonTOF");
  TH2D *fHistdEdxvRigvProtonTOF        = (TH2D*)inlist->FindObject("fHistdEdxvRigvProtonTOF");
  TH2D *fHistdEdxvRigvDeuteronTOF      = (TH2D*)inlist->FindObject("fHistdEdxvRigvDeuteronTOF");
  TH2D *fHistdEdxvRigvTritonTOF        = (TH2D*)inlist->FindObject("fHistdEdxvRigvTritonTOF");
  TH2D *fHistdEdxvRigvHelium3TOF       = (TH2D*)inlist->FindObject("fHistdEdxvRigvHelium3TOF");
  TH2D *fHistdEdxvRigvAlphaTOF         = (TH2D*)inlist->FindObject("fHistdEdxvRigvAlphaTOF");

  TH2D *fHistdEdxvRigvElectronTPCTOF   = (TH2D*)inlist->FindObject("fHistdEdxvRigvElectronTPCTOF");
  TH2D *fHistdEdxvRigvPionTPCTOF       = (TH2D*)inlist->FindObject("fHistdEdxvRigvPionTPCTOF");
  TH2D *fHistdEdxvRigvKaonTPCTOF       = (TH2D*)inlist->FindObject("fHistdEdxvRigvKaonTPCTOF");
  TH2D *fHistdEdxvRigvProtonTPCTOF     = (TH2D*)inlist->FindObject("fHistdEdxvRigvProtonTPCTOF");
  TH2D *fHistdEdxvRigvDeuteronTPCTOF   = (TH2D*)inlist->FindObject("fHistdEdxvRigvDeuteronTPCTOF");
  TH2D *fHistdEdxvRigvTritonTPCTOF     = (TH2D*)inlist->FindObject("fHistdEdxvRigvTritonTPCTOF");
  TH2D *fHistdEdxvRigvHelium3TPCTOF    = (TH2D*)inlist->FindObject("fHistdEdxvRigvHelium3TPCTOF");
  TH2D *fHistdEdxvRigvHelium3TPCTOF_bb = (TH2D*)inlist->FindObject("fHistdEdxvRigvHelium3TPCTOF_bb");
  TH2D *fHistdEdxvRigvAlphaTPCTOF      = (TH2D*)inlist->FindObject("fHistdEdxvRigvAlphaTPCTOF");
  TH2D *fHistdEdxvRigvAlphaTPCTOF_bb   = (TH2D*)inlist->FindObject("fHistdEdxvRigvAlphaTPCTOF_bb");


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create legends
  ///_________________________________________________________________________________________________
  TLegend *leg = new TLegend(.3,.6,.8,.9);
  if(Collisions == "PbPb_11h") {
    leg->AddEntry((TObject*)0,"Pb-Pb @ #sqrt{s_{NN}} = 2.76 TeV","");
    leg->AddEntry((TObject*)0,"run LHC11h pass2","");
  }
  else if(Collisions == "pPb_13b") {
    leg->AddEntry((TObject*)0,"p-Pb @ #sqrt{s_{NN}} = 5.023 TeV","");
    leg->AddEntry((TObject*)0,"run LHC13b pass3","");
  }
  else if(Collisions == "pPb_13c") {
    leg->AddEntry((TObject*)0,"p-Pb @ #sqrt{s_{NN}} = 5.023 TeV","");
    leg->AddEntry((TObject*)0,"run LHC13c pass2","");
  }
  else { //(Collisions == "pp_15f")
    leg->AddEntry((TObject*)0,"p-p @ #sqrt{s}} = 13 TeV","");
    leg->AddEntry((TObject*)0,"run LHC15f pass2","");
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create and fill canvases
  ///_________________________________________________________________________________________________
  TCanvas *c7 = new TCanvas("c7",Form("dEdx vs. Rigidity w TPC preselection (Sigma=%i) for separated particles",Sigma));
  c7->Divide(3,3);
  for(Int_t i = 1; i <= 9 ; i++) {
    c7->cd(i)->SetTicks();
  }
  c7->cd(1)->SetLogx();
  c7->cd(1)->SetLogz();
    fHistdEdxvRig->Draw("colz");
  c7->cd(2)->SetLogx();
  c7->cd(2)->SetLogz();
    fHistdEdxvRigvElectronTPC->Draw("colz");
  c7->cd(3)->SetLogx();
  c7->cd(3)->SetLogz();
    fHistdEdxvRigvPionTPC->Draw("colz");
  c7->cd(4)->SetLogx();
  c7->cd(4)->SetLogz();
    fHistdEdxvRigvKaonTPC->Draw("colz");
  c7->cd(5)->SetLogx();
  c7->cd(5)->SetLogz();
    fHistdEdxvRigvProtonTPC->Draw("colz");
  c7->cd(6)->SetLogx();
  c7->cd(6)->SetLogz();
    fHistdEdxvRigvDeuteronTPC->Draw("colz");
  c7->cd(7)->SetLogx();
  c7->cd(7)->SetLogz();
    fHistdEdxvRigvTritonTPC->Draw("colz");
  c7->cd(8)->SetLogx();
  c7->cd(8)->SetLogz();
    //fHistdEdxvRigvHelium3TPC->Draw("colz");
    fHistdEdxvRigvHelium3TPC_bb->Draw("colz");
  c7->cd(9)->SetLogx();
  c7->cd(9)->SetLogz();
    //fHistdEdxvRigvAlphaTPC->Draw("colz");
    fHistdEdxvRigvAlphaTPC_bb->Draw("colz");

  TCanvas *c8 = new TCanvas("c8",Form("dEdx vs. Rigidity w TOF preselection (Sigma=%i) for separated particles",Sigma));
  c8->Divide(3,3);
  for(Int_t i = 1; i <= 9 ; i++) {
    c8->cd(i)->SetTicks();
  }
  c8->cd(1)->SetLogx();
  c8->cd(1)->SetLogz();
    fHistdEdxvRig->Draw("colz");
  c8->cd(2)->SetLogx();
  c8->cd(2)->SetLogz();
    fHistdEdxvRigvElectronTOF->Draw("colz");
  c8->cd(3)->SetLogx();
  c8->cd(3)->SetLogz();
    fHistdEdxvRigvPionTOF->Draw("colz");
  c8->cd(4)->SetLogx();
  c8->cd(4)->SetLogz();
    fHistdEdxvRigvKaonTOF->Draw("colz");
  c8->cd(5)->SetLogx();
  c8->cd(5)->SetLogz();
    fHistdEdxvRigvProtonTOF->Draw("colz");
  c8->cd(6)->SetLogx();
  c8->cd(6)->SetLogz();
    fHistdEdxvRigvDeuteronTOF->Draw("colz");
  c8->cd(7)->SetLogx();
  c8->cd(7)->SetLogz();
    fHistdEdxvRigvTritonTOF->Draw("colz");
  c8->cd(8)->SetLogx();
  c8->cd(8)->SetLogz();
    fHistdEdxvRigvHelium3TOF->Draw("colz");
  c8->cd(9)->SetLogx();
  c8->cd(9)->SetLogz();
    fHistdEdxvRigvAlphaTOF->Draw("colz");

  TCanvas *c9 = new TCanvas("c9",Form("dEdx vs. Rigidity w TPC & TOF preselection (Sigma=%i) for separated particles",Sigma));
  c9->Divide(3,3);
  for(Int_t i = 1; i <= 9 ; i++) {
    c9->cd(i)->SetTicks();
  }
  c9->cd(1)->SetLogx();
  c9->cd(1)->SetLogz();
    fHistdEdxvRig->Draw("colz");
  c9->cd(2)->SetLogx();
  c9->cd(2)->SetLogz();
    fHistdEdxvRigvElectronTPCTOF->Draw("colz");
  c9->cd(3)->SetLogx();
  c9->cd(3)->SetLogz();
    fHistdEdxvRigvPionTPCTOF->Draw("colz");
  c9->cd(4)->SetLogx();
  c9->cd(4)->SetLogz();
    fHistdEdxvRigvKaonTPCTOF->Draw("colz");
  c9->cd(5)->SetLogx();
  c9->cd(5)->SetLogz();
    fHistdEdxvRigvProtonTPCTOF->Draw("colz");
  c9->cd(6)->SetLogx();
  c9->cd(6)->SetLogz();
    fHistdEdxvRigvDeuteronTPCTOF->Draw("colz");
  c9->cd(7)->SetLogx();
  c9->cd(7)->SetLogz();
    fHistdEdxvRigvTritonTPCTOF->Draw("colz");
  c9->cd(8)->SetLogx();
  c9->cd(8)->SetLogz();
    //fHistdEdxvRigvHelium3TPCTOF->Draw("colz");
    fHistdEdxvRigvHelium3TPCTOF_bb->Draw("colz");
  c9->cd(9)->SetLogx();
  c9->cd(9)->SetLogz();
    //fHistdEdxvRigvAlphaTPCTOF->Draw("colz");
    fHistdEdxvRigvAlphaTPCTOF_bb->Draw("colz");


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Save canvases
  ///_________________________________________________________________________________________________
  c7 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_dEdxTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
  c8 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_dEdxTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
  c9 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_dEdxTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));


}
