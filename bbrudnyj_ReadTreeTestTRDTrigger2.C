/*--------------------------------------------------------------------------------------------------*\
|                                                                                                    |
|    Macro that makes projection histograms with branches of "nmartin_TestTriggerTRDTree.root"       |
|    - TRD PID vs. pT for all particles                                                              |
|      - with preselection of TPC, TOF, TPC && TOF for:<seperated particles> = Electrons,Pions,Kaons,|
|                                                      Protons, Deuterons, Tritons, Helium3, Alphas  |
|    - TRD dE/dx vs. rigidity for: all particles                                                     |
|      - with preselection of TPC, TOF, TPC && TOF for:<seperated particles>                         |
|                                                                                                    |
|                                                                                                    |
| Author: Benjamin Brudnyj (2016)                                                                    | 
\*--------------------------------------------------------------------------------------------------*/
void bbrudnyj_ReadTreeTestTRDTrigger2() {

  static const Int_t binDD = 100000;
  Int_t N = 0; // for Event printing in loop

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Get data from file /lustre/nyx/alice/users/bbrudnyj/trunk/nmartin_nuclei/AliAnalysisTaskTestTriggerTRD.cxx
  ///_________________________________________________________________________________________________

  TFile *ftree= TFile::Open("/u/brudnyj/Desktop/Bachelor/codes/PbPb_11h/bbrudnyj_ReadTreeTestTRDTrigger2/nmartin_TestTriggerTRDTree.root"); char *Collisions ="PbPb_11h";N=1;
  //TFile *ftree= TFile::Open("/lustre/nyx/alice/users/bbrudnyj/train/V012.pPb/                                                             /nmartin_TestTriggerTRDTree.root"); char *Collisions ="pPb_13b";N=10;
  //TFile *ftree= TFile::Open("/lustre/nyx/alice/users/bbrudnyj/train/V012.pPb/                                                             /nmartin_TestTriggerTRDTree.root"); char *Collisions ="pPb_13c";N=10;
  //TFile *ftree= TFile::Open("/lustre/nyx/alice/users/bbrudnyj/train/V012.pp/                                                              /nmartin_TestTriggerTRDTree.root"); char *Collisions ="pp_15f"; N=10;


  TTree *t = (TTree *)ftree->Get("tree;1");
  //t->Print();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Define Variables
  /// Declaration of leaf types
  ///_________________________________________________
  Int_t      TotEvent               = t->GetEntries();
  Int_t      mod                    = 100000*N;

  Int_t      Sigma                  = 3;
  Int_t      LessStatisticFactor    = 10000;

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

  Int_t      fItrk;

  //Double_t   *fEta  = new Double_t[binDD];
  Double_t   *fdEdx = new Double_t[binDD];
  Double_t   *fpT = new Double_t[binDD];
  //Float_t    *fCentrality = new Float_t[binDD];
  Double_t   *fRigidity = new Double_t[binDD];
  //Double_t   *fRigidityTRD = new Double_t[binDD];
  Int_t      *fTRDPID = new Int_t[binDD];
  //Float_t    *fTOFMass = new Float_t[binDD];
  //Bool_t     *fLabelGTU = new Bool_t[binDD];
  //Bool_t     *fLabelHW = new Bool_t[binDD];
  //Bool_t     *fElectron = new Bool_t[binDD];
  //Bool_t     *fPion = new Bool_t[binDD];
  //Bool_t     *fProton = new Bool_t[binDD];
  //Bool_t     *fHelium3 = new Bool_t[binDD];
  Float_t    *fTPCnSigmaElectron = new Float_t[binDD];
  Float_t    *fTPCnSigmaPion = new Float_t[binDD];
  Float_t    *fTPCnSigmaKaon = new Float_t[binDD];
  Float_t    *fTPCnSigmaProton = new Float_t[binDD];
  Float_t    *fTPCnSigmaDeuteron = new Float_t[binDD];
  Float_t    *fTPCnSigmaTriton = new Float_t[binDD];
  Float_t    *fTPCnSigmaHelium3 = new Float_t[binDD];
  Float_t    *fTPCnSigmaAlpha = new Float_t[binDD];
  Float_t    *fTOFnSigmaElectron = new Float_t[binDD];
  Float_t    *fTOFnSigmaPion = new Float_t[binDD];
  Float_t    *fTOFnSigmaKaon = new Float_t[binDD];
  Float_t    *fTOFnSigmaProton = new Float_t[binDD];
  Float_t    *fTOFnSigmaDeuteron = new Float_t[binDD];
  Float_t    *fTOFnSigmaTriton = new Float_t[binDD];
  Float_t    *fTOFnSigmaHelium3 = new Float_t[binDD];
  Float_t    *fTOFnSigmaAlpha = new Float_t[binDD];

  t->SetBranchAddress("fItrk",&fItrk);
  t->SetBranchAddress("fdEdx",fdEdx);
  t->SetBranchAddress("fpT",fpT);
  //t->SetBranchAddress("fEta",fEta);
  //t->SetBranchAddress("fCentrality",fCentrality);
  t->SetBranchAddress("fRigidity",fRigidity);
  //t->SetBranchAddress("fRigidityTRD",fRigidityTRD);
  t->SetBranchAddress("fTRDPID",fTRDPID);
  //t->SetBranchAddress("fTOFMass",fTOFMass);
  //t->SetBranchAddress("fLabelGTU",fLabelGTU);
  //t->SetBranchAddress("fLabelHW",fLabelHW);
  //t->SetBranchAddress("fElectron",fElectron);
  //t->SetBranchAddress("fPion",fPion);
  //t->SetBranchAddress("fProton",fProton);
  //t->SetBranchAddress("fHelium3",fHelium3);
  t->SetBranchAddress("fTPCnSigmaElectron",fTPCnSigmaElectron);
  t->SetBranchAddress("fTPCnSigmaPion",fTPCnSigmaPion);
  t->SetBranchAddress("fTPCnSigmaKaon",fTPCnSigmaKaon);
  t->SetBranchAddress("fTPCnSigmaProton",fTPCnSigmaProton);
  t->SetBranchAddress("fTPCnSigmaDeuteron",fTPCnSigmaDeuteron);
  t->SetBranchAddress("fTPCnSigmaTriton",fTPCnSigmaTriton);
  t->SetBranchAddress("fTPCnSigmaHelium3",fTPCnSigmaHelium3);
  t->SetBranchAddress("fTPCnSigmaAlpha",fTPCnSigmaAlpha);
  t->SetBranchAddress("fTOFnSigmaElectron",fTOFnSigmaElectron);
  t->SetBranchAddress("fTOFnSigmaPion",fTOFnSigmaPion);
  t->SetBranchAddress("fTOFnSigmaKaon",fTOFnSigmaKaon);
  t->SetBranchAddress("fTOFnSigmaProton",fTOFnSigmaProton);
  t->SetBranchAddress("fTOFnSigmaDeuteron",fTOFnSigmaDeuteron);
  t->SetBranchAddress("fTOFnSigmaTriton",fTOFnSigmaTriton);
  t->SetBranchAddress("fTOFnSigmaHelium3",fTOFnSigmaHelium3);
  t->SetBranchAddress("fTOFnSigmaAlpha",fTOFnSigmaAlpha);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create and fill histograms
  ///_______________________________________________________
  printf("\nParticle preselection with Sigma = %i\n",Sigma);
  printf("Creating 117 histograms and 9 canvases...\n");
  TH2D *fHistTRDpTvPID = new TH2D("fHistTRDpTvPID","TRD PID vs. #it{p}_{T}",500,-25,25,260,1,260);
  t->Project("fHistTRDpTvPID","fTRDPID:fpT");
  fHistTRDpTvPID->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHistTRDpTvPID->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
  printf("1/117 histograms done\n");

  TH2D *fHistdEdxvRig = new TH2D("fHistdEdxvRig","TPC dEdx vs. Rigidity",2000,0.3,20,4000,5,800);
  t->Project("fHistdEdxvRig","fdEdx:fRigidity");
  fHistdEdxvRig->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
  fHistdEdxvRig->GetYaxis()->SetTitle("dE/dx (a.u.)");
  printf("2/117 histograms done\n");

  // PID vs. pT

    // TPC preselection histograms
    TH2D *fHistTRDpTvPIDvElectronTPC = new TH2D("fHistTRDpTvPIDvElectronTPC","TRD PID vs. #it{p}_{T} for Electron",500,-25,25,260,1,260);
    fHistTRDpTvPIDvElectronTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvElectronTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvPionTPC = new TH2D("fHistTRDpTvPIDvPionTPC","TRD PID vs. #it{p}_{T} for Pion",500,-25,25,260,1,260);
    fHistTRDpTvPIDvPionTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvPionTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvKaonTPC = new TH2D("fHistTRDpTvPIDvKaonTPC","TRD PID vs. #it{p}_{T} for Kaon",500,-25,25,260,1,260);
    fHistTRDpTvPIDvKaonTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvKaonTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvProtonTPC = new TH2D("fHistTRDpTvPIDvProtonTPC","TRD PID vs. #it{p}_{T} for Proton",500,-25,25,260,1,260);
    fHistTRDpTvPIDvProtonTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvProtonTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvDeuteronTPC = new TH2D("fHistTRDpTvPIDvDeuteronTPC","TRD PID vs. #it{p}_{T} for Deuteron",500,-25,25,260,1,260);
    fHistTRDpTvPIDvDeuteronTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvDeuteronTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvDeuteronTPC_clear = new TH2D("fHistTRDpTvPIDvDeuteronTPC_clear",Form("TRD PID vs. #it{p}_{T} (w #it{p}/#it{z} < %.1f GeV/it#{c}) for Deuteron",ClearCutDeuteronTPC),500,-25,25,260,1,260);
    fHistTRDpTvPIDvDeuteronTPC_clear->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvDeuteronTPC_clear->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvTritonTPC = new TH2D("fHistTRDpTvPIDvTritonTPC","TRD PID vs. #it{p}_{T} for Triton (TPC ps)",500,-25,25,260,1,260);
    fHistTRDpTvPIDvTritonTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvTritonTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvTritonTPC_clear = new TH2D("fHistTRDpTvPIDvTritonTPC_clear",Form("TRD PID vs. #it{p}_{T} (w #it{p}/#it{z} < %.1f GeV/it#{c}) for Triton",ClearCutTritonTPC),500,-25,25,260,1,260);
    fHistTRDpTvPIDvTritonTPC_clear->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvTritonTPC_clear->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvHelium3TPC = new TH2D("fHistTRDpTvPIDvHelium3TPC","TRD PID vs. #it{p}_{T} for Helium3",160,-8,8,260,1,260);
    fHistTRDpTvPIDvHelium3TPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvHelium3TPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAlphaTPC = new TH2D("fHistTRDpTvPIDvAlphaTPC","TRD PID vs. #it{p}_{T} for Alpha",160,-8,8,260,1,260);
    fHistTRDpTvPIDvAlphaTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAlphaTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvHelium3AlphaTPC = new TH2D("fHistTRDpTvPIDvHelium3AlphaTPC","TRD PID vs. #it{p}_{T} for Helium3 & Alpha",160,-8,8,260,1,260);
    fHistTRDpTvPIDvHelium3AlphaTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvHelium3AlphaTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");

    TH2D *fHistTRDpTvPIDvAllexceptDeuteronTPC_clear = new TH2D("fHistTRDpTvPIDvAllexceptDeuteronTPC_clear","TRD PID vs. #it{p}_{T} for All except Deuteron",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptDeuteronTPC_clear->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptDeuteronTPC_clear->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAllexceptTritonTPC_clear = new TH2D("fHistTRDpTvPIDvAllexceptTritonTPC_clear","TRD PID vs. #it{p}_{T} for All except Triton",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptTritonTPC_clear->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptTritonTPC_clear->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAllexceptHelium3TPC = new TH2D("fHistTRDpTvPIDvAllexceptHelium3TPC","TRD PID vs. #it{p}_{T} for All except Helium3",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptHelium3TPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptHelium3TPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAllexceptAlphaTPC = new TH2D("fHistTRDpTvPIDvAllexceptAlphaTPC","TRD PID vs. #it{p}_{T} for All except Alpha",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptAlphaTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptAlphaTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAllexceptHelium3AlphaTPC = new TH2D("fHistTRDpTvPIDvAllexceptHelium3AlphaTPC","TRD PID vs. #it{p}_{T} for All except Helium3 & Alpha",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptHelium3AlphaTPC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptHelium3AlphaTPC->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");


    // TOF preselection histograms
    TH2D *fHistTRDpTvPIDvElectronTOF = new TH2D("fHistTRDpTvPIDvElectronTOF","TRD PID vs. #it{p}_{T} for Electron",500,-25,25,260,1,260);
    fHistTRDpTvPIDvElectronTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvElectronTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvPionTOF = new TH2D("fHistTRDpTvPIDvPionTOF","TRD PID vs. #it{p}_{T} for Pion",500,-25,25,260,1,260);
    fHistTRDpTvPIDvPionTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvPionTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvKaonTOF = new TH2D("fHistTRDpTvPIDvKaonTOF","TRD PID vs. #it{p}_{T} for Kaon",500,-25,25,260,1,260);
    fHistTRDpTvPIDvKaonTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvKaonTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvProtonTOF = new TH2D("fHistTRDpTvPIDvProtonTOF","TRD PID vs. #it{p}_{T} for Proton",500,-25,25,260,1,260);
    fHistTRDpTvPIDvProtonTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvProtonTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvDeuteronTOF = new TH2D("fHistTRDpTvPIDvDeuteronTOF","TRD PID vs. #it{p}_{T} for Deuteron",500,-25,25,260,1,260);
    fHistTRDpTvPIDvDeuteronTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvDeuteronTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvTritonTOF = new TH2D("fHistTRDpTvPIDvTritonTOF","TRD PID vs. #it{p}_{T} for Triton",500,-25,25,260,1,260);
    fHistTRDpTvPIDvTritonTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvTritonTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvHelium3TOF = new TH2D("fHistTRDpTvPIDvHelium3TOF","TRD PID vs. #it{p}_{T} for Helium3",500,-25,25,260,1,260);
    fHistTRDpTvPIDvHelium3TOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvHelium3TOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAlphaTOF = new TH2D("fHistTRDpTvPIDvAlphaTOF","TRD PID vs. #it{p}_{T} for Alpha",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAlphaTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAlphaTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");

    // TPC & TOF preselection histograms
    TH2D *fHistTRDpTvPIDvElectronTPCTOF = new TH2D("fHistTRDpTvPIDvElectronTPCTOF","TRD PID vs. #it{p}_{T} for Electron",500,-25,25,260,1,260);
    fHistTRDpTvPIDvElectronTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvElectronTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvPionTPCTOF = new TH2D("fHistTRDpTvPIDvPionTPCTOF","TRD PID vs. #it{p}_{T} for Pion",500,-25,25,260,1,260);
    fHistTRDpTvPIDvPionTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvPionTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvKaonTPCTOF = new TH2D("fHistTRDpTvPIDvKaonTPCTOF","TRD PID vs. #it{p}_{T} for Kaon",500,-25,25,260,1,260);
    fHistTRDpTvPIDvKaonTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvKaonTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvProtonTPCTOF = new TH2D("fHistTRDpTvPIDvProtonTPCTOF","TRD PID vs. #it{p}_{T} for Proton",500,-25,25,260,1,260);
    fHistTRDpTvPIDvProtonTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvProtonTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvDeuteronTPCTOF = new TH2D("fHistTRDpTvPIDvDeuteronTPCTOF","TRD PID vs. #it{p}_{T} for Deuteron",500,-25,25,260,1,260);
    fHistTRDpTvPIDvDeuteronTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvDeuteronTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvDeuteronTPCTOF_clear = new TH2D("fHistTRDpTvPIDvDeuteronTPCTOF_clear",Form("TRD PID vs. #it{p}_{T} (w #it{p}/#it{z} < %.1f GeV/it#{c}) for Deuteron",ClearCutDeuteronTPCTOF),500,-25,25,260,1,260);
    fHistTRDpTvPIDvDeuteronTPCTOF_clear->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvDeuteronTPCTOF_clear->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvTritonTPCTOF = new TH2D("fHistTRDpTvPIDvTritonTPCTOF","TRD PID vs. #it{p}_{T} for Triton (TPC & TOF ps)",500,-25,25,260,1,260);
    fHistTRDpTvPIDvTritonTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvTritonTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvTritonTPCTOF_clear = new TH2D("fHistTRDpTvPIDvTritonTPCTOF_clear",Form("TRD PID vs. #it{p}_{T} (w #it{p}/#it{z} < %.1f GeV/it#{c}) for Triton",ClearCutTritonTPCTOF),500,-25,25,260,1,260);
    fHistTRDpTvPIDvTritonTPCTOF_clear->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvTritonTPCTOF_clear->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvHelium3TPCTOF = new TH2D("fHistTRDpTvPIDvHelium3TPCTOF","TRD PID vs. #it{p}_{T} for Helium3",160,-8,8,260,1,260);
    fHistTRDpTvPIDvHelium3TPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvHelium3TPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAlphaTPCTOF = new TH2D("fHistTRDpTvPIDvAlphaTPCTOF","TRD PID vs. #it{p}_{T} for Alpha",160,-8,8,260,1,260);
    fHistTRDpTvPIDvAlphaTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAlphaTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvHelium3AlphaTPCTOF = new TH2D("fHistTRDpTvPIDvHelium3AlphaTPCTOF","TRD PID vs. #it{p}_{T} for Helium3 & Alpha",160,-8,8,260,1,260);
    fHistTRDpTvPIDvHelium3AlphaTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvHelium3AlphaTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");

    TH2D *fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear = new TH2D("fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear","TRD PID vs. #it{p}_{T} for All except Deuteron",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear = new TH2D("fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear","TRD PID vs. #it{p}_{T} for All except Triton",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAllexceptHelium3TPCTOF = new TH2D("fHistTRDpTvPIDvAllexceptHelium3TPCTOF","TRD PID vs. #it{p}_{T} for All except Helium3",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptHelium3TPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptHelium3TPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAllexceptAlphaTPCTOF = new TH2D("fHistTRDpTvPIDvAllexceptAlphaTPCTOF","TRD PID vs. #it{p}_{T} for All except Alpha",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptAlphaTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptAlphaTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");
    TH2D *fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF = new TH2D("fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF","TRD PID vs. #it{p}_{T} for All except Helium3 & Alpha",500,-25,25,260,1,260);
    fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF->GetYaxis()->SetTitle("<Q_{s}> (a.u.)");


  // dE/dx vs. Rigidity

    // TPC preselection histograms
    TH2D *fHistdEdxvRigvElectronTPC = new TH2D("fHistdEdxvRigvElectronTPC","TPC dEdx vs. Rigidity for Electron",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvElectronTPC->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvElectronTPC->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvPionTPC = new TH2D("fHistdEdxvRigvPionTPC","TPC dEdx vs. Rigidity for Pion",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvPionTPC->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvPionTPC->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvKaonTPC = new TH2D("fHistdEdxvRigvKaonTPC","TPC dEdx vs. Rigidity for Kaon",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvKaonTPC->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvKaonTPC->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvProtonTPC = new TH2D("fHistdEdxvRigvProtonTPC","TPC dEdx vs. Rigidity for Proton",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvProtonTPC->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvProtonTPC->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvDeuteronTPC = new TH2D("fHistdEdxvRigvDeuteronTPC","TPC dEdx vs. Rigidity for Deuteron",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvDeuteronTPC->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvDeuteronTPC->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvTritonTPC = new TH2D("fHistdEdxvRigvTritonTPC","TPC dEdx vs. Rigidity for Triton",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvTritonTPC->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvTritonTPC->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvHelium3TPC = new TH2D("fHistdEdxvRigvHelium3TPC","TPC dEdx vs. Rigidity for Helium3",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvHelium3TPC->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvHelium3TPC->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvHelium3TPC_bb = new TH2D("fHistdEdxvRigvHelium3TPC_bb","TPC dEdx vs. Rigidity for Helium3",1000,0.3,20,795,5,800);
    fHistdEdxvRigvHelium3TPC_bb->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvHelium3TPC_bb->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvAlphaTPC = new TH2D("fHistdEdxvRigvAlphaTPC","TPC dEdx vs. Rigidity for Alpha",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvAlphaTPC->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvAlphaTPC->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvAlphaTPC_bb = new TH2D("fHistdEdxvRigvAlphaTPC_bb","TPC dEdx vs. Rigidity for Alpha",1000,0.3,20,795,5,800);
    fHistdEdxvRigvAlphaTPC_bb->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvAlphaTPC_bb->GetYaxis()->SetTitle("dE/dx (a.u.)");

    // TOF preselection histograms
    TH2D *fHistdEdxvRigvElectronTOF = new TH2D("fHistdEdxvRigvElectronTOF","TPC dEdx vs. Rigidity for Electron",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvElectronTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvElectronTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvPionTOF = new TH2D("fHistdEdxvRigvPionTOF","TPC dEdx vs. Rigidity for Pion",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvPionTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvPionTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvKaonTOF = new TH2D("fHistdEdxvRigvKaonTOF","TPC dEdx vs. Rigidity for Kaon",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvKaonTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvKaonTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvProtonTOF = new TH2D("fHistdEdxvRigvProtonTOF","TPC dEdx vs. Rigidity for Proton",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvProtonTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvProtonTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvDeuteronTOF = new TH2D("fHistdEdxvRigvDeuteronTOF","TPC dEdx vs. Rigidity for Deuteron",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvDeuteronTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvDeuteronTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvTritonTOF = new TH2D("fHistdEdxvRigvTritonTOF","TPC dEdx vs. Rigidity for Triton",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvTritonTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvTritonTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvHelium3TOF = new TH2D("fHistdEdxvRigvHelium3TOF","TPC dEdx vs. Rigidity for Helium3",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvHelium3TOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvHelium3TOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvAlphaTOF = new TH2D("fHistdEdxvRigvAlphaTOF","TPC dEdx vs. Rigidity for Alpha",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvAlphaTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvAlphaTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");

    // TPC & TOF preselection histograms
    TH2D *fHistdEdxvRigvElectronTPCTOF = new TH2D("fHistdEdxvRigvElectronTPCTOF","TPC dEdx vs. Rigidity for Electron",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvElectronTPCTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvElectronTPCTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvPionTPCTOF = new TH2D("fHistdEdxvRigvPionTPCTOF","TPC dEdx vs. Rigidity for Pion",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvPionTPCTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvPionTPCTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvKaonTPCTOF = new TH2D("fHistdEdxvRigvKaonTPCTOF","TPC dEdx vs. Rigidity for Kaon",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvKaonTPCTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvKaonTPCTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvProtonTPCTOF = new TH2D("fHistdEdxvRigvProtonTPCTOF","TPC dEdx vs. Rigidity for Proton",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvProtonTPCTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvProtonTPCTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvDeuteronTPCTOF = new TH2D("fHistdEdxvRigvDeuteronTPCTOF","TPC dEdx vs. Rigidity for Deuteron",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvDeuteronTPCTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvDeuteronTPCTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvTritonTPCTOF = new TH2D("fHistdEdxvRigvTritonTPCTOF","TPC dEdx vs. Rigidity for Triton",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvTritonTPCTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvTritonTPCTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvHelium3TPCTOF = new TH2D("fHistdEdxvRigvHelium3TPCTOF","TPC dEdx vs. Rigidity for Helium3",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvHelium3TPCTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvHelium3TPCTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvHelium3TPCTOF_bb = new TH2D("fHistdEdxvRigvHelium3TPCTOF_bb","TPC dEdx vs. Rigidity for Helium3",1000,0.3,20,795,5,800);
    fHistdEdxvRigvHelium3TPCTOF_bb->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvHelium3TPCTOF_bb->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvAlphaTPCTOF = new TH2D("fHistdEdxvRigvAlphaTPCTOF","TPC dEdx vs. Rigidity for Alpha",2000,0.3,20,4000,5,800);
    fHistdEdxvRigvAlphaTPCTOF->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvAlphaTPCTOF->GetYaxis()->SetTitle("dE/dx (a.u.)");
    TH2D *fHistdEdxvRigvAlphaTPCTOF_bb = new TH2D("fHistdEdxvRigvAlphaTPCTOF_bb","TPC dEdx vs. Rigidity for Alpha",1000,0.3,20,795,5,800);
    fHistdEdxvRigvAlphaTPCTOF_bb->GetXaxis()->SetTitle("#frac{#it{p}}{#it{z}} (GeV/#it{c})");
    fHistdEdxvRigvAlphaTPCTOF_bb->GetYaxis()->SetTitle("dE/dx (a.u.)");


  // loop to fill histograms
  printf("loop over %i Events... Fill 3.-64. histogram...\n",TotEvent);
  if(LessStatisticFactor==1)  printf("");
  else if(LessStatisticFactor==2) printf("go only through every %ind Event\n",LessStatisticFactor);
  else if(LessStatisticFactor==3) printf("go only through every %ird Event\n",LessStatisticFactor);
  else printf("go only through every %ith Event\n",LessStatisticFactor);

  for(Int_t i=1; i < TotEvent; i++) {
    if(i%mod==0) printf("%i Events done\n",i);

    t->GetEntry(i);

    if(fItrk >= binDD) {
      printf("fItrk to big: %i",fItrk);
      return;
    }
    if(i%LessStatisticFactor==0) {
      for(Int_t j=0; j < fItrk; j++) {
        // Electron
        if(abs(fTOFnSigmaElectron[j]) < Sigma) {
     	  fHistTRDpTvPIDvElectronTOF->Fill(fpT[j],fTRDPID[j]);
  	  fHistdEdxvRigvElectronTOF ->Fill(fRigidity[j],fdEdx[j]);
        }
        if(abs(fTPCnSigmaElectron[j]) < Sigma) {
          fHistTRDpTvPIDvElectronTPC->Fill(fpT[j],fTRDPID[j]);
          fHistdEdxvRigvElectronTPC ->Fill(fRigidity[j],fdEdx[j]);

          if(abs(fTOFnSigmaElectron[j]) < Sigma) {
	    fHistTRDpTvPIDvElectronTPCTOF->Fill(fpT[j],fTRDPID[j]);
	    fHistdEdxvRigvElectronTPCTOF ->Fill(fRigidity[j],fdEdx[j]);
          }
        }


        // Pion
        if(abs(fTOFnSigmaPion[j]) < Sigma) {
	  fHistTRDpTvPIDvPionTOF->Fill(fpT[j],fTRDPID[j]);
	  fHistdEdxvRigvPionTOF ->Fill(fRigidity[j],fdEdx[j]);
        }
        if(abs(fTPCnSigmaPion[j]) < Sigma) {
          fHistTRDpTvPIDvPionTPC->Fill(fpT[j],fTRDPID[j]);
          fHistdEdxvRigvPionTPC ->Fill(fRigidity[j],fdEdx[j]);

          if(abs(fTOFnSigmaPion[j]) < Sigma) {
	    fHistTRDpTvPIDvPionTPCTOF->Fill(fpT[j],fTRDPID[j]);
	    fHistdEdxvRigvPionTPCTOF ->Fill(fRigidity[j],fdEdx[j]);
	  }
        }


        // Kaon
        if(abs(fTOFnSigmaKaon[j]) < Sigma) {
	  fHistTRDpTvPIDvKaonTOF->Fill(fpT[j],fTRDPID[j]);
	  fHistdEdxvRigvKaonTOF ->Fill(fRigidity[j],fdEdx[j]);
        }
        if(abs(fTPCnSigmaKaon[j]) < Sigma) {
          fHistTRDpTvPIDvKaonTPC->Fill(fpT[j],fTRDPID[j]);
          fHistdEdxvRigvKaonTPC ->Fill(fRigidity[j],fdEdx[j]);

          if(abs(fTOFnSigmaKaon[j]) < Sigma) {
	    fHistTRDpTvPIDvKaonTPCTOF->Fill(fpT[j],fTRDPID[j]);
	    fHistdEdxvRigvKaonTPCTOF ->Fill(fRigidity[j],fdEdx[j]);
          }
        }


        // Proton
        if(abs(fTOFnSigmaProton[j]) < Sigma) {
	  fHistTRDpTvPIDvProtonTOF->Fill(fpT[j],fTRDPID[j]);
	  fHistdEdxvRigvProtonTOF ->Fill(fRigidity[j],fdEdx[j]);
        }
        if(abs(fTPCnSigmaProton[j]) < Sigma) {
          fHistTRDpTvPIDvProtonTPC->Fill(fpT[j],fTRDPID[j]);
          fHistdEdxvRigvProtonTPC ->Fill(fRigidity[j],fdEdx[j]);

          if(abs(fTOFnSigmaProton[j]) < Sigma) {
	    fHistTRDpTvPIDvProtonTPCTOF->Fill(fpT[j],fTRDPID[j]);
	    fHistdEdxvRigvProtonTPCTOF ->Fill(fRigidity[j],fdEdx[j]);
          }
        }


        // Deuteron
        if(abs(fTOFnSigmaDeuteron[j]) < Sigma) {
	  fHistTRDpTvPIDvDeuteronTOF->Fill(fpT[j],fTRDPID[j]);
	  fHistdEdxvRigvDeuteronTOF ->Fill(fRigidity[j],fdEdx[j]);
        }
        if(abs(fTPCnSigmaDeuteron[j]) < Sigma) {
          fHistTRDpTvPIDvDeuteronTPC->Fill(fpT[j],fTRDPID[j]);
          fHistdEdxvRigvDeuteronTPC ->Fill(fRigidity[j],fdEdx[j]);

          if(TMath::Abs(fRigidity[j]) < ClearCutDeuteronTPC) fHistTRDpTvPIDvDeuteronTPC_clear->Fill(fpT[j],fTRDPID[j]);

          if(abs(fTOFnSigmaDeuteron[j]) < Sigma) {
	    fHistTRDpTvPIDvDeuteronTPCTOF->Fill(fpT[j],fTRDPID[j]);
	    fHistdEdxvRigvDeuteronTPCTOF ->Fill(fRigidity[j],fdEdx[j]);

            if(TMath::Abs(fRigidity[j]) < ClearCutDeuteronTPCTOF) fHistTRDpTvPIDvDeuteronTPCTOF_clear->Fill(fpT[j],fTRDPID[j]);
	  }
        }
        //if(!(abs(fTPCnSigmaDeuteron[j]) < Sigma && TMath::Abs(fRigidity[j]) < ClearCutDeuteronTPC))                                          fHistTRDpTvPIDvAllexceptDeuteronTPC_clear   ->Fill(fpT[j],fTRDPID[j]);
        //if(!(abs(fTPCnSigmaDeuteron[j]) < Sigma && TMath::Abs(fRigidity[j]) < ClearCutDeuteronTPCTOF && abs(fTOFnSigmaDeuteron[j]) < Sigma)) fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear->Fill(fpT[j],fTRDPID[j]);
        if(!(abs(fTPCnSigmaDeuteron[j]) < Sigma)) {
	  if(!(TMath::Abs(fRigidity[j]) < ClearCutDeuteronTPC))      fHistTRDpTvPIDvAllexceptDeuteronTPC_clear   ->Fill(fpT[j],fTRDPID[j]);
          if(!(abs(fTOFnSigmaDeuteron[j]) < Sigma)) {
	    if(!(TMath::Abs(fRigidity[j]) < ClearCutDeuteronTPCTOF)) fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear->Fill(fpT[j],fTRDPID[j]);
	  }
	}



        // Triton
        if(abs(fTOFnSigmaTriton[j]) < Sigma) {
	  fHistTRDpTvPIDvTritonTOF->Fill(fpT[j],fTRDPID[j]);
	  fHistdEdxvRigvTritonTOF ->Fill(fRigidity[j],fdEdx[j]);
        }
        if(abs(fTPCnSigmaTriton[j]) < Sigma) {
          fHistTRDpTvPIDvTritonTPC->Fill(fpT[j],fTRDPID[j]);
          fHistdEdxvRigvTritonTPC ->Fill(fRigidity[j],fdEdx[j]);

          if(TMath::Abs(fRigidity[j]) < ClearCutTritonTPC) fHistTRDpTvPIDvTritonTPC_clear->Fill(fpT[j],fTRDPID[j]);

          if(abs(fTOFnSigmaTriton[j]) < Sigma) {
	    fHistTRDpTvPIDvTritonTPCTOF->Fill(fpT[j],fTRDPID[j]);
	    fHistdEdxvRigvTritonTPCTOF ->Fill(fRigidity[j],fdEdx[j]);

            if(TMath::Abs(fRigidity[j]) < ClearCutTritonTPCTOF) fHistTRDpTvPIDvTritonTPCTOF_clear->Fill(fpT[j],fTRDPID[j]);
	  }
        }
	  //if(!(abs(fTPCnSigmaTriton[j]) < Sigma && TMath::Abs(fRigidity[j]) < ClearCutTritonTPC))                                        fHistTRDpTvPIDvAllexceptTritonTPC_clear   ->Fill(fpT[j],fTRDPID[j]);
	  //if(!(abs(fTPCnSigmaTriton[j]) < Sigma && TMath::Abs(fRigidity[j]) < ClearCutTritonTPCTOF && abs(fTOFnSigmaTriton[j]) < Sigma)) fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear->Fill(fpT[j],fTRDPID[j]);

        if(!(abs(fTPCnSigmaTriton[j]) < Sigma)) {
	  if(!(TMath::Abs(fRigidity[j]) < ClearCutTritonTPC))     fHistTRDpTvPIDvAllexceptTritonTPC_clear   ->Fill(fpT[j],fTRDPID[j]);
          if(!(abs(fTOFnSigmaTriton[j]) < Sigma)) {
            if(!(TMath::Abs(fRigidity[j]) < ClearCutTritonTPCTOF)) fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear->Fill(fpT[j],fTRDPID[j]);
	  }
	}

        // Helium3
        if(abs(fTOFnSigmaHelium3[j]) < Sigma) {
	  fHistTRDpTvPIDvHelium3TOF->Fill(fpT[j],fTRDPID[j]);
	  fHistdEdxvRigvHelium3TOF ->Fill(fRigidity[j],fdEdx[j]);
        }
        if(abs(fTPCnSigmaHelium3[j]) < Sigma) {
          fHistTRDpTvPIDvHelium3TPC  ->Fill(fpT[j],fTRDPID[j]);
          fHistdEdxvRigvHelium3TPC   ->Fill(fRigidity[j],fdEdx[j]);
          fHistdEdxvRigvHelium3TPC_bb->Fill(fRigidity[j],fdEdx[j]);

          if(abs(fTOFnSigmaHelium3[j]) < Sigma) {
	    fHistTRDpTvPIDvHelium3TPCTOF  ->Fill(fpT[j],fTRDPID[j]);
	    fHistdEdxvRigvHelium3TPCTOF   ->Fill(fRigidity[j],fdEdx[j]);
	    fHistdEdxvRigvHelium3TPCTOF_bb->Fill(fRigidity[j],fdEdx[j]);
          }
        }

        if(!(abs(fTPCnSigmaHelium3[j]) < Sigma)) {
	  fHistTRDpTvPIDvAllexceptHelium3TPC->Fill(fpT[j],fTRDPID[j]);
          if(!(abs(fTOFnSigmaHelium3[j]) < Sigma)) fHistTRDpTvPIDvAllexceptHelium3TPCTOF->Fill(fpT[j],fTRDPID[j]);
	}


        // Alpha
        if(abs(fTOFnSigmaAlpha[j]) < Sigma) {
	  fHistTRDpTvPIDvAlphaTOF->Fill(fpT[j],fTRDPID[j]);
	  fHistdEdxvRigvAlphaTOF ->Fill(fRigidity[j],fdEdx[j]);
        }
        if(abs(fTPCnSigmaAlpha[j]) < Sigma) {
          fHistTRDpTvPIDvAlphaTPC  ->Fill(fpT[j],fTRDPID[j]);
          fHistdEdxvRigvAlphaTPC   ->Fill(fRigidity[j],fdEdx[j]);
          fHistdEdxvRigvAlphaTPC_bb->Fill(fRigidity[j],fdEdx[j]);

          if(abs(fTOFnSigmaAlpha[j]) < Sigma) {
	    fHistTRDpTvPIDvAlphaTPCTOF  ->Fill(fpT[j],fTRDPID[j]);
	    fHistdEdxvRigvAlphaTPCTOF   ->Fill(fRigidity[j],fdEdx[j]);
	    fHistdEdxvRigvAlphaTPCTOF_bb->Fill(fRigidity[j],fdEdx[j]);
          }
        }

        if(!(abs(fTPCnSigmaAlpha[j]) < Sigma)) {
	  fHistTRDpTvPIDvAllexceptAlphaTPC->Fill(fpT[j],fTRDPID[j]);
          if(!(abs(fTOFnSigmaAlpha[j]) < Sigma)) fHistTRDpTvPIDvAllexceptAlphaTPCTOF->Fill(fpT[j],fTRDPID[j]);
	}


        // Helium3 & Alpha
        if(abs(fTPCnSigmaHelium3[j]) < Sigma || abs(fTPCnSigmaAlpha[j]) < Sigma) {
          fHistTRDpTvPIDvHelium3AlphaTPC->Fill(fpT[j],fTRDPID[j]);

          if(abs(fTOFnSigmaHelium3[j]) < Sigma || abs(fTOFnSigmaAlpha[j]) < Sigma) {
	    fHistTRDpTvPIDvHelium3AlphaTPCTOF->Fill(fpT[j],fTRDPID[j]);
          }
        }

        if(!(abs(fTPCnSigmaHelium3[j]) < Sigma || abs(fTPCnSigmaAlpha[j]) < Sigma)) {
          fHistTRDpTvPIDvAllexceptHelium3AlphaTPC->Fill(fpT[j],fTRDPID[j]);

          if(!(abs(fTOFnSigmaHelium3[j]) < Sigma || abs(fTOFnSigmaAlpha[j]) < Sigma)) {
	  fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF->Fill(fpT[j],fTRDPID[j]);
          }
        }
      }
    }
  }
  printf("64/117 histograms done\n");

  // projection histograms TPC
  printf("53 projection histograms...\n");
    TH2D *fHistTRDpTvPIDvElectronTPC_Clone = (TH2D *)fHistTRDpTvPIDvElectronTPC->Clone("fHistTRDpTvPIDvElectronTPC_Clone");
    TH1D *fHistTRDPIDElectronTPC = fHistTRDpTvPIDvElectronTPC_Clone->ProjectionY("fHistTRDPIDElectronTPC",1,260);
    fHistTRDPIDElectronTPC->SetTitle("TRD PID distribution for Electron;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDElectronTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDElectronTPC->SetLineColor(4);
    Double_t intElectronTPC = fHistTRDPIDElectronTPC->Integral(1,260);
    fHistTRDPIDElectronTPC->Scale(1/intElectronTPC);

    TH2D *fHistTRDpTvPIDvPionTPC_Clone = (TH2D *)fHistTRDpTvPIDvPionTPC->Clone("fHistTRDpTvPIDvPionTPC_Clone");
    TH1D *fHistTRDPIDPionTPC = fHistTRDpTvPIDvPionTPC_Clone->ProjectionY("fHistTRDPIDPionTPC",1,260);
    fHistTRDPIDPionTPC->SetTitle("TRD PID distribution for Pion;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDPionTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDPionTPC->SetLineColor(2);
    Double_t intPionTPC = fHistTRDPIDPionTPC->Integral(1,260);
    fHistTRDPIDPionTPC->Scale(1/intPionTPC);

    TH2D *fHistTRDpTvPIDvKaonTPC_Clone = (TH2D *)fHistTRDpTvPIDvKaonTPC->Clone("fHistTRDpTvPIDvKaonTPC_Clone");
    TH1D *fHistTRDPIDKaonTPC = fHistTRDpTvPIDvKaonTPC_Clone->ProjectionY("fHistTRDPIDKaonTPC",1,260);
    fHistTRDPIDKaonTPC->SetTitle("TRD PID distribution for Kaon;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDKaonTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDKaonTPC->SetLineColor(5);
    Double_t intKaonTPC = fHistTRDPIDKaonTPC->Integral(1,260);
    fHistTRDPIDKaonTPC->Scale(1/intKaonTPC);

    TH2D *fHistTRDpTvPIDvProtonTPC_Clone = (TH2D *)fHistTRDpTvPIDvProtonTPC->Clone("fHistTRDpTvPIDvProtonTPC_Clone");
    TH1D *fHistTRDPIDProtonTPC = fHistTRDpTvPIDvProtonTPC_Clone->ProjectionY("fHistTRDPIDProtonTPC",1,260);
    fHistTRDPIDProtonTPC->SetTitle("TRD PID distribution for Proton;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDProtonTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDProtonTPC->SetLineColor(3);
    Double_t intProtonTPC = fHistTRDPIDProtonTPC->Integral(1,260);
    fHistTRDPIDProtonTPC->Scale(1/intProtonTPC);

    TH2D *fHistTRDpTvPIDvDeuteronTPC_Clone = (TH2D *)fHistTRDpTvPIDvDeuteronTPC->Clone("fHistTRDpTvPIDvDeuteronTPC_Clone");
    TH1D *fHistTRDPIDDeuteronTPC = fHistTRDpTvPIDvDeuteronTPC_Clone->ProjectionY("fHistTRDPIDDeuteronTPC",1,260);
    fHistTRDPIDDeuteronTPC->SetTitle("TRD PID distribution for Deuteron;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDDeuteronTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDDeuteronTPC->SetLineColor(7);
    Double_t intDeuteronTPC = fHistTRDPIDDeuteronTPC->Integral(1,260);
    fHistTRDPIDDeuteronTPC->Scale(1/intDeuteronTPC);

    TH2D *fHistTRDpTvPIDvDeuteronTPC_clear_Clone = (TH2D *)fHistTRDpTvPIDvDeuteronTPC_clear->Clone("fHistTRDpTvPIDvDeuteronTPC_clear_Clone");
    TH1D *fHistTRDPIDDeuteronTPC_clear = fHistTRDpTvPIDvDeuteronTPC_clear_Clone->ProjectionY("fHistTRDPIDDeuteronTPC_clear",1,260);
    fHistTRDPIDDeuteronTPC_clear->SetTitle(Form("TRD PID distribution for Deuteron (Rigidity < %.1f GeV/#it{cz});<Q_{s}> (a.u.);Normalized Yield",ClearCutDeuteronTPC));
    fHistTRDPIDDeuteronTPC_clear->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDDeuteronTPC_clear->SetLineColor(7);
    Double_t intDeuteronTPC_clear = fHistTRDPIDDeuteronTPC_clear->Integral(1,260);
    fHistTRDPIDDeuteronTPC_clear->Scale(1/intDeuteronTPC_clear);

    TH2D *fHistTRDpTvPIDvTritonTPC_Clone = (TH2D *)fHistTRDpTvPIDvTritonTPC->Clone("fHistTRDpTvPIDvTritonTPC_Clone");
    TH1D *fHistTRDPIDTritonTPC = fHistTRDpTvPIDvTritonTPC_Clone->ProjectionY("fHistTRDPIDTritonTPC",1,260);
    fHistTRDPIDTritonTPC->SetTitle("TRD PID distribution for Triton;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDTritonTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDTritonTPC->SetLineColor(28);
    Double_t intTritonTPC = fHistTRDPIDTritonTPC->Integral(1,260);
    fHistTRDPIDTritonTPC->Scale(1/intTritonTPC);

    TH2D *fHistTRDpTvPIDvTritonTPC_clear_Clone = (TH2D *)fHistTRDpTvPIDvTritonTPC_clear->Clone("fHistTRDpTvPIDvTritonTPC_clear_Clone");
    TH1D *fHistTRDPIDTritonTPC_clear = fHistTRDpTvPIDvTritonTPC_clear_Clone->ProjectionY("fHistTRDPIDTritonTPC_clear",1,260);
    fHistTRDPIDTritonTPC_clear->SetTitle(Form("TRD PID distribution for Triton (Rigidity < %.1f GeV/#it{cz});<Q_{s}> (a.u.);Normalized Yield",ClearCutTritonTPC));
    fHistTRDPIDTritonTPC_clear->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDTritonTPC_clear->SetLineColor(28);
    Double_t intTritonTPC_clear = fHistTRDPIDTritonTPC_clear->Integral(1,260);
    fHistTRDPIDTritonTPC_clear->Scale(1/intTritonTPC_clear);

    TH2D *fHistTRDpTvPIDvHelium3TPC_Clone = (TH2D *)fHistTRDpTvPIDvHelium3TPC->Clone("fHistTRDpTvPIDvHelium3TPC_Clone");
    TH1D *fHistTRDPIDHelium3TPC = fHistTRDpTvPIDvHelium3TPC_Clone->ProjectionY("fHistTRDPIDHelium3TPC",1,260);
    fHistTRDPIDHelium3TPC->SetTitle("TRD PID distribution for Helium3;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDHelium3TPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDHelium3TPC->SetLineColor(1);
    Double_t intHelium3TPC = fHistTRDPIDHelium3TPC->Integral(1,260);
    fHistTRDPIDHelium3TPC->Scale(1/intHelium3TPC);

    TH2D *fHistTRDpTvPIDvAlphaTPC_Clone = (TH2D *)fHistTRDpTvPIDvAlphaTPC->Clone("fHistTRDpTvPIDvAlphaTPC_Clone");
    TH1D *fHistTRDPIDAlphaTPC = fHistTRDpTvPIDvAlphaTPC_Clone->ProjectionY("fHistTRDPIDAlphaTPC",1,260);
    fHistTRDPIDAlphaTPC->SetTitle("TRD PID distribution for Alpha;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAlphaTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAlphaTPC->SetLineColor(17);
    Double_t intAlphaTPC = fHistTRDPIDAlphaTPC->Integral(1,260);
    fHistTRDPIDAlphaTPC->Scale(1/intAlphaTPC);

    TH2D *fHistTRDpTvPIDvHelium3AlphaTPC_Clone = (TH2D *)fHistTRDpTvPIDvHelium3AlphaTPC->Clone("fHistTRDpTvPIDvHelium3AlphaTPC_Clone");
    TH1D *fHistTRDPIDHelium3AlphaTPC = fHistTRDpTvPIDvHelium3AlphaTPC_Clone->ProjectionY("fHistTRDPIDHelium3AlphaTPC",1,260);
    fHistTRDPIDHelium3AlphaTPC->SetTitle("TRD PID distribution for Helium3 & Alpha;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDHelium3AlphaTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDHelium3AlphaTPC->SetLineColor(4);
    Double_t intHelium3AlphaTPC = fHistTRDPIDHelium3AlphaTPC->Integral(1,260);
    fHistTRDPIDHelium3AlphaTPC->Scale(1/intHelium3AlphaTPC);

    // Allexcept
    TH2D *fHistTRDpTvPIDvAllexceptDeuteronTPC_clear_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptDeuteronTPC_clear->Clone("fHistTRDpTvPIDvAllexceptDeuteronTPC_clear_Clone");
    TH1D *fHistTRDPIDAllexceptDeuteronTPC_clear = fHistTRDpTvPIDvAllexceptDeuteronTPC_clear_Clone->ProjectionY("fHistTRDPIDAllexceptDeuteronTPC_clear",1,260);
    fHistTRDPIDAllexceptDeuteronTPC_clear->SetTitle("TRD PID distribution for All except Deuteron;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptDeuteronTPC_clear->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptDeuteronTPC_clear->SetLineColor(2);
    Double_t intAllexceptDeuteronTPC_clear = fHistTRDPIDAllexceptDeuteronTPC_clear->Integral(1,260);
    fHistTRDPIDAllexceptDeuteronTPC_clear->Scale(1/intAllexceptDeuteronTPC_clear);

    TH2D *fHistTRDpTvPIDvAllexceptTritonTPC_clear_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptTritonTPC_clear->Clone("fHistTRDpTvPIDvAllexceptTritonTPC_clear_Clone");
    TH1D *fHistTRDPIDAllexceptTritonTPC_clear = fHistTRDpTvPIDvAllexceptTritonTPC_clear_Clone->ProjectionY("fHistTRDPIDAllexceptTritonTPC_clear",1,260);
    fHistTRDPIDAllexceptTritonTPC_clear->SetTitle("TRD PID distribution for All except Triton;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptTritonTPC_clear->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptTritonTPC_clear->SetLineColor(2);
    Double_t intAllexceptTritonTPC_clear = fHistTRDPIDAllexceptTritonTPC_clear->Integral(1,260);
    fHistTRDPIDAllexceptTritonTPC_clear->Scale(1/intAllexceptTritonTPC_clear);

    TH2D *fHistTRDpTvPIDvAllexceptHelium3TPC_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptHelium3TPC->Clone("fHistTRDpTvPIDvAllexceptHelium3TPC_Clone");
    TH1D *fHistTRDPIDAllexceptHelium3TPC = fHistTRDpTvPIDvAllexceptHelium3TPC_Clone->ProjectionY("fHistTRDPIDAllexceptHelium3TPC",1,260);
    fHistTRDPIDAllexceptHelium3TPC->SetTitle("TRD PID distribution for All except Helium3;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptHelium3TPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptHelium3TPC->SetLineColor(2);
    Double_t intAllexceptHelium3TPC = fHistTRDPIDAllexceptHelium3TPC->Integral(1,260);
    fHistTRDPIDAllexceptHelium3TPC->Scale(1/intAllexceptHelium3TPC);

    TH2D *fHistTRDpTvPIDvAllexceptAlphaTPC_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptAlphaTPC->Clone("fHistTRDpTvPIDvAllexceptAlphaTPC_Clone");
    TH1D *fHistTRDPIDAllexceptAlphaTPC = fHistTRDpTvPIDvAllexceptAlphaTPC_Clone->ProjectionY("fHistTRDPIDAllexceptAlphaTPC",1,260);
    fHistTRDPIDAllexceptAlphaTPC->SetTitle("TRD PID distribution for All except Alpha;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptAlphaTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptAlphaTPC->SetLineColor(2);
    Double_t intAllexceptAlphaTPC = fHistTRDPIDAllexceptAlphaTPC->Integral(1,260);
    fHistTRDPIDAllexceptAlphaTPC->Scale(1/intAllexceptAlphaTPC);

    TH2D *fHistTRDpTvPIDvAllexceptHelium3AlphaTPC_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptHelium3AlphaTPC->Clone("fHistTRDpTvPIDvAllexceptHelium3AlphaTPC_Clone");
    TH1D *fHistTRDPIDAllexceptHelium3AlphaTPC = fHistTRDpTvPIDvAllexceptHelium3AlphaTPC_Clone->ProjectionY("fHistTRDPIDAllexceptHelium3AlphaTPC",1,260);
    fHistTRDPIDAllexceptHelium3AlphaTPC->SetTitle("TRD PID distribution for All except Helium3 & Alpha;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptHelium3AlphaTPC->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptHelium3AlphaTPC->SetLineColor(2);
    Double_t intAllexceptHelium3AlphaTPC = fHistTRDPIDAllexceptHelium3AlphaTPC->Integral(1,260);
    fHistTRDPIDAllexceptHelium3AlphaTPC->Scale(1/intAllexceptHelium3AlphaTPC);

    // proj_Clone
    TH1D *fHistTRDPIDallTPC = (TH1D *)fHistTRDPIDPionTPC->Clone("fHistTRDPIDallTPC");
    fHistTRDPIDallTPC->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallDeuteronTPC = (TH1D *)fHistTRDPIDAllexceptAlphaTPC->Clone("fHistTRDPIDallDeuteronTPC");
    fHistTRDPIDallDeuteronTPC->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallTritonTPC = (TH1D *)fHistTRDPIDAllexceptAlphaTPC->Clone("fHistTRDPIDallTritonTPC");
    fHistTRDPIDallTritonTPC->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallHelium3TPC = (TH1D *)fHistTRDPIDAllexceptAlphaTPC->Clone("fHistTRDPIDallHelium3TPC");
    fHistTRDPIDallHelium3TPC->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallAlphaTPC = (TH1D *)fHistTRDPIDAllexceptAlphaTPC->Clone("fHistTRDPIDallAlphaTPC");
    fHistTRDPIDallAlphaTPC->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallHelium3AlphaTPC = (TH1D *)fHistTRDPIDAllexceptAlphaTPC->Clone("fHistTRDPIDallHelium3AlphaTPC");
    fHistTRDPIDallHelium3AlphaTPC->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");


    TH1D *fHistTRDPIDElectronTPC_Clone = (TH1D *)fHistTRDPIDElectronTPC->Clone("fHistTRDPIDElectronTPC_Clone");
    fHistTRDPIDElectronTPC_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDPionTPC_Clone = (TH1D *)fHistTRDPIDPionTPC->Clone("fHistTRDPIDPionTPC_Clone");
    fHistTRDPIDPionTPC_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDKaonTPC_Clone = (TH1D *)fHistTRDPIDKaonTPC->Clone("fHistTRDPIDKaonTPC_Clone");
    fHistTRDPIDKaonTPC_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDProtonTPC_Clone = (TH1D *)fHistTRDPIDProtonTPC->Clone("fHistTRDPIDProtonTPC_Clone");
    fHistTRDPIDProtonTPC_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDDeuteronTPC_Clone = (TH1D *)fHistTRDPIDDeuteronTPC->Clone("fHistTRDPIDDeuteronTPC_Clone");
    fHistTRDPIDDeuteronTPC_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDDeuteronTPC_clear_Clone = (TH1D *)fHistTRDPIDDeuteronTPC_clear->Clone("fHistTRDPIDDeuteronTPC_clear_Clone");
    fHistTRDPIDDeuteronTPC_clear_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDDeuteronTPC_clear_Clone2 = (TH1D *)fHistTRDPIDDeuteronTPC_clear->Clone("fHistTRDPIDDeuteronTPC_clear_Clone2");
    fHistTRDPIDDeuteronTPC_clear_Clone2->SetLineColor(4);
    TH1D *fHistTRDPIDTritonTPC_Clone = (TH1D *)fHistTRDPIDTritonTPC->Clone("fHistTRDPIDTritonTPC_Clone");
    fHistTRDPIDTritonTPC_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDTritonTPC_clear_Clone = (TH1D *)fHistTRDPIDTritonTPC_clear->Clone("fHistTRDPIDTritonTPC_clear_Clone");
    fHistTRDPIDTritonTPC_clear_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDTritonTPC_clear_Clone2 = (TH1D *)fHistTRDPIDTritonTPC_clear->Clone("fHistTRDPIDTritonTPC_clear_Clone2");
    fHistTRDPIDTritonTPC_clear_Clone2->SetLineColor(4);
    TH1D *fHistTRDPIDHelium3TPC_Clone = (TH1D *)fHistTRDPIDHelium3TPC->Clone("fHistTRDPIDHelium3TPC_Clone");
    fHistTRDPIDHelium3TPC_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDHelium3TPC_Clone2 = (TH1D *)fHistTRDPIDHelium3TPC->Clone("fHistTRDPIDHelium3TPC_Clone2");
    fHistTRDPIDHelium3TPC_Clone2->SetLineColor(4);
    TH1D *fHistTRDPIDAlphaTPC_Clone = (TH1D *)fHistTRDPIDAlphaTPC->Clone("fHistTRDPIDAlphaTPC_Clone");
    fHistTRDPIDAlphaTPC_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDAlphaTPC_Clone2 = (TH1D *)fHistTRDPIDAlphaTPC->Clone("fHistTRDPIDAlphaTPC_Clone2");
    fHistTRDPIDAlphaTPC_Clone2->SetLineColor(4);


  // projection histograms TOF
    TH2D *fHistTRDpTvPIDvElectronTOF_Clone = (TH2D *)fHistTRDpTvPIDvElectronTOF->Clone("fHistTRDpTvPIDvElectronTOF_Clone");
    TH1D *fHistTRDPIDElectronTOF = fHistTRDpTvPIDvElectronTOF_Clone->ProjectionY("fHistTRDPIDElectronTOF",1,260);
    fHistTRDPIDElectronTOF->SetTitle("TRD PID distribution for Electron;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDElectronTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDElectronTOF->SetLineColor(4);
    Double_t intElectronTOF = fHistTRDPIDElectronTOF->Integral(1,260);
    fHistTRDPIDElectronTOF->Scale(1/intElectronTOF);

    TH2D *fHistTRDpTvPIDvPionTOF_Clone = (TH2D *)fHistTRDpTvPIDvPionTOF->Clone("fHistTRDpTvPIDvPionTOF_Clone");
    TH1D *fHistTRDPIDPionTOF = fHistTRDpTvPIDvPionTOF_Clone->ProjectionY("fHistTRDPIDPionTOF",1,260);
    fHistTRDPIDPionTOF->SetTitle("TRD PID distribution for Pion;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDPionTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDPionTOF->SetLineColor(2);
    Double_t intPionTOF = fHistTRDPIDPionTOF->Integral(1,260);
    fHistTRDPIDPionTOF->Scale(1/intPionTOF);

    TH2D *fHistTRDpTvPIDvKaonTOF_Clone = (TH2D *)fHistTRDpTvPIDvKaonTOF->Clone("fHistTRDpTvPIDvKaonTOF_Clone");
    TH1D *fHistTRDPIDKaonTOF = fHistTRDpTvPIDvKaonTOF_Clone->ProjectionY("fHistTRDPIDKaonTOF",1,260);
    fHistTRDPIDKaonTOF->SetTitle("TRD PID distribution for Kaon;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDKaonTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDKaonTOF->SetLineColor(5);
    Double_t intKaonTOF = fHistTRDPIDKaonTOF->Integral(1,260);
    fHistTRDPIDKaonTOF->Scale(1/intKaonTOF);

    TH2D *fHistTRDpTvPIDvProtonTOF_Clone = (TH2D *)fHistTRDpTvPIDvProtonTOF->Clone("fHistTRDpTvPIDvProtonTOF_Clone");
    TH1D *fHistTRDPIDProtonTOF = fHistTRDpTvPIDvProtonTOF_Clone->ProjectionY("fHistTRDPIDProtonTOF",1,260);
    fHistTRDPIDProtonTOF->SetTitle("TRD PID distribution for Proton;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDProtonTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDProtonTOF->SetLineColor(3);
    Double_t intProtonTOF = fHistTRDPIDProtonTOF->Integral(1,260);
    fHistTRDPIDProtonTOF->Scale(1/intProtonTOF);

    TH2D *fHistTRDpTvPIDvDeuteronTOF_Clone = (TH2D *)fHistTRDpTvPIDvDeuteronTOF->Clone("fHistTRDpTvPIDvDeuteronTOF_Clone");
    TH1D *fHistTRDPIDDeuteronTOF = fHistTRDpTvPIDvDeuteronTOF_Clone->ProjectionY("fHistTRDPIDDeuteronTOF",1,260);
    fHistTRDPIDDeuteronTOF->SetTitle("TRD PID distribution for Deuteron;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDDeuteronTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDDeuteronTOF->SetLineColor(7);
    Double_t intDeuteronTOF = fHistTRDPIDDeuteronTOF->Integral(1,260);
    fHistTRDPIDDeuteronTOF->Scale(1/intDeuteronTOF);

    TH2D *fHistTRDpTvPIDvTritonTOF_Clone = (TH2D *)fHistTRDpTvPIDvTritonTOF->Clone("fHistTRDpTvPIDvTritonTOF_Clone");
    TH1D *fHistTRDPIDTritonTOF = fHistTRDpTvPIDvTritonTOF_Clone->ProjectionY("fHistTRDPIDTritonTOF",1,260);
    fHistTRDPIDTritonTOF->SetTitle("TRD PID distribution for Triton;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDTritonTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDTritonTOF->SetLineColor(28);
    Double_t intTritonTOF = fHistTRDPIDTritonTOF->Integral(1,260);
    fHistTRDPIDTritonTOF->Scale(1/intTritonTOF);

    TH2D *fHistTRDpTvPIDvHelium3TOF_Clone = (TH2D *)fHistTRDpTvPIDvHelium3TOF->Clone("fHistTRDpTvPIDvHelium3TOF_Clone");
    TH1D *fHistTRDPIDHelium3TOF = fHistTRDpTvPIDvHelium3TOF_Clone->ProjectionY("fHistTRDPIDHelium3TOF",1,260);
    fHistTRDPIDHelium3TOF->SetTitle("TRD PID distribution for Helium3;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDHelium3TOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDHelium3TOF->SetLineColor(1);
    Double_t intHelium3TOF = fHistTRDPIDHelium3TOF->Integral(1,260);
    fHistTRDPIDHelium3TOF->Scale(1/intHelium3TOF);

    TH2D *fHistTRDpTvPIDvAlphaTOF_Clone = (TH2D *)fHistTRDpTvPIDvAlphaTOF->Clone("fHistTRDpTvPIDvAlphaTOF_Clone");
    TH1D *fHistTRDPIDAlphaTOF = fHistTRDpTvPIDvAlphaTOF_Clone->ProjectionY("fHistTRDPIDAlphaTOF",1,260);
    fHistTRDPIDAlphaTOF->SetTitle("TRD PID distribution for Alpha;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAlphaTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAlphaTOF->SetLineColor(17);
    Double_t intAlphaTOF = fHistTRDPIDAlphaTOF->Integral(1,260);
    fHistTRDPIDAlphaTOF->Scale(1/intAlphaTOF);

    // proj_Clone
    TH1D *fHistTRDPIDallTOF = (TH1D *)fHistTRDPIDPionTOF->Clone("fHistTRDPIDallTOF");
    fHistTRDPIDallTOF->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");

    TH1D *fHistTRDPIDElectronTOF_Clone = (TH1D *)fHistTRDPIDElectronTOF->Clone("fHistTRDPIDElectronTOF_Clone");
    fHistTRDPIDElectronTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDPionTOF_Clone = (TH1D *)fHistTRDPIDPionTOF->Clone("fHistTRDPIDPionTOF_Clone");
    fHistTRDPIDPionTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDKaonTOF_Clone = (TH1D *)fHistTRDPIDKaonTOF->Clone("fHistTRDPIDKaonTOF_Clone");
    fHistTRDPIDKaonTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDProtonTOF_Clone = (TH1D *)fHistTRDPIDProtonTOF->Clone("fHistTRDPIDProtonTOF_Clone");
    fHistTRDPIDProtonTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDDeuteronTOF_Clone = (TH1D *)fHistTRDPIDDeuteronTOF->Clone("fHistTRDPIDDeuteronTOF_Clone");
    fHistTRDPIDDeuteronTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDTritonTOF_Clone = (TH1D *)fHistTRDPIDTritonTOF->Clone("fHistTRDPIDTritonTOF_Clone");
    fHistTRDPIDTritonTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDHelium3TOF_Clone = (TH1D *)fHistTRDPIDHelium3TOF->Clone("fHistTRDPIDHelium3TOF_Clone");
    fHistTRDPIDHelium3TOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDAlphaTOF_Clone = (TH1D *)fHistTRDPIDAlphaTOF->Clone("fHistTRDPIDAlphaTOF_Clone");
    fHistTRDPIDAlphaTOF_Clone->SetLineColor(1);


  // projection histograms TPC & TOF
    TH2D *fHistTRDpTvPIDvElectronTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvElectronTPCTOF->Clone("fHistTRDpTvPIDvElectronTPCTOF_Clone");
    TH1D *fHistTRDPIDElectronTPCTOF = fHistTRDpTvPIDvElectronTPCTOF_Clone->ProjectionY("fHistTRDPIDElectronTPCTOF",1,260);
    fHistTRDPIDElectronTPCTOF->SetTitle("TRD PID distribution for Electron;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDElectronTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDElectronTPCTOF->SetLineColor(4);
    Double_t intElectronTPCTOF = fHistTRDPIDElectronTPCTOF->Integral(1,260);
    fHistTRDPIDElectronTPCTOF->Scale(1/intElectronTPCTOF);

    TH2D *fHistTRDpTvPIDvPionTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvPionTPCTOF->Clone("fHistTRDpTvPIDvPionTPCTOF_Clone");
    TH1D *fHistTRDPIDPionTPCTOF = fHistTRDpTvPIDvPionTPCTOF_Clone->ProjectionY("fHistTRDPIDPionTPCTOF",1,260);
    fHistTRDPIDPionTPCTOF->SetTitle("TRD PID distribution for Pion;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDPionTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDPionTPCTOF->SetLineColor(2);
    Double_t intPionTPCTOF = fHistTRDPIDPionTPCTOF->Integral(1,260);
    fHistTRDPIDPionTPCTOF->Scale(1/intPionTPCTOF);

    TH2D *fHistTRDpTvPIDvKaonTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvKaonTPCTOF->Clone("fHistTRDpTvPIDvKaonTPCTOF_Clone");
    TH1D *fHistTRDPIDKaonTPCTOF = fHistTRDpTvPIDvKaonTPCTOF_Clone->ProjectionY("fHistTRDPIDKaonTPCTOF",1,260);
    fHistTRDPIDKaonTPCTOF->SetTitle("TRD PID distribution for Kaon;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDKaonTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDKaonTPCTOF->SetLineColor(5);
    Double_t intKaonTPCTOF = fHistTRDPIDKaonTPCTOF->Integral(1,260);
    fHistTRDPIDKaonTPCTOF->Scale(1/intKaonTPCTOF);

    TH2D *fHistTRDpTvPIDvProtonTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvProtonTPCTOF->Clone("fHistTRDpTvPIDvProtonTPCTOF_Clone");
    TH1D *fHistTRDPIDProtonTPCTOF = fHistTRDpTvPIDvProtonTPCTOF_Clone->ProjectionY("fHistTRDPIDProtonTPCTOF",1,260);
    fHistTRDPIDProtonTPCTOF->SetTitle("TRD PID distribution for Proton;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDProtonTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDProtonTPCTOF->SetLineColor(3);
    Double_t intProtonTPCTOF = fHistTRDPIDProtonTPCTOF->Integral(1,260);
    fHistTRDPIDProtonTPCTOF->Scale(1/intProtonTPCTOF);

    TH2D *fHistTRDpTvPIDvDeuteronTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvDeuteronTPCTOF->Clone("fHistTRDpTvPIDvDeuteronTPCTOF_Clone");
    TH1D *fHistTRDPIDDeuteronTPCTOF = fHistTRDpTvPIDvDeuteronTPCTOF_Clone->ProjectionY("fHistTRDPIDDeuteronTPCTOF",1,260);
    fHistTRDPIDDeuteronTPCTOF->SetTitle("TRD PID distribution for Deuteron;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDDeuteronTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDDeuteronTPCTOF->SetLineColor(7);
    Double_t intDeuteronTPCTOF = fHistTRDPIDDeuteronTPCTOF->Integral(1,260);
    fHistTRDPIDDeuteronTPCTOF->Scale(1/intDeuteronTPCTOF);

    TH2D *fHistTRDpTvPIDvDeuteronTPCTOF_clear_Clone = (TH2D *)fHistTRDpTvPIDvDeuteronTPCTOF_clear->Clone("fHistTRDpTvPIDvDeuteronTPCTOF_clear_Clone");
    TH1D *fHistTRDPIDDeuteronTPCTOF_clear = fHistTRDpTvPIDvDeuteronTPCTOF_clear_Clone->ProjectionY("fHistTRDPIDDeuteronTPCTOF_clear",1,260);
    fHistTRDPIDDeuteronTPCTOF_clear->SetTitle(Form("TRD PID distribution for Deuteron (Rigidity < %.1f GeV/#it{cz});<Q_{s}> (a.u.);Normalized Yield",ClearCutDeuteronTPCTOF));
    fHistTRDPIDDeuteronTPCTOF_clear->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDDeuteronTPCTOF_clear->SetLineColor(7);
    Double_t intDeuteronTPCTOF_clear = fHistTRDPIDDeuteronTPCTOF_clear->Integral(1,260);
    fHistTRDPIDDeuteronTPCTOF_clear->Scale(1/intDeuteronTPCTOF_clear);

    TH2D *fHistTRDpTvPIDvTritonTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvTritonTPCTOF->Clone("fHistTRDpTvPIDvTritonTPCTOF_Clone");
    TH1D *fHistTRDPIDTritonTPCTOF = fHistTRDpTvPIDvTritonTPCTOF_Clone->ProjectionY("fHistTRDPIDTritonTPCTOF",1,260);
    fHistTRDPIDTritonTPCTOF->SetTitle("TRD PID distribution for Triton;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDTritonTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDTritonTPCTOF->SetLineColor(28);
    Double_t intTritonTPCTOF = fHistTRDPIDTritonTPCTOF->Integral(1,260);
    fHistTRDPIDTritonTPCTOF->Scale(1/intTritonTPCTOF);

    TH2D *fHistTRDpTvPIDvTritonTPCTOF_clear_Clone = (TH2D *)fHistTRDpTvPIDvTritonTPCTOF_clear->Clone("fHistTRDpTvPIDvTritonTPCTOF_clear_Clone");
    TH1D *fHistTRDPIDTritonTPCTOF_clear = fHistTRDpTvPIDvTritonTPCTOF_clear_Clone->ProjectionY("fHistTRDPIDTritonTPCTOF_clear",1,260);
    fHistTRDPIDTritonTPCTOF_clear->SetTitle(Form("TRD PID distribution for Triton (Rigidity < %.1f GeV/#it{cz});<Q_{s}> (a.u.);Normalized Yield",ClearCutTritonTPCTOF));
    fHistTRDPIDTritonTPCTOF_clear->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDTritonTPCTOF_clear->SetLineColor(28);
    Double_t intTritonTPCTOF_clear = fHistTRDPIDTritonTPCTOF_clear->Integral(1,260);
    fHistTRDPIDTritonTPCTOF_clear->Scale(1/intTritonTPCTOF_clear);

    TH2D *fHistTRDpTvPIDvHelium3TPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvHelium3TPCTOF->Clone("fHistTRDpTvPIDvHelium3TPCTOF_Clone");
    TH1D *fHistTRDPIDHelium3TPCTOF = fHistTRDpTvPIDvHelium3TPCTOF_Clone->ProjectionY("fHistTRDPIDHelium3TPCTOF",1,260);
    fHistTRDPIDHelium3TPCTOF->SetTitle("TRD PID distribution for Helium3;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDHelium3TPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDHelium3TPCTOF->SetLineColor(1);
    Double_t intHelium3TPCTOF = fHistTRDPIDHelium3TPCTOF->Integral(1,260);
    fHistTRDPIDHelium3TPCTOF->Scale(1/intHelium3TPCTOF);

    TH2D *fHistTRDpTvPIDvAlphaTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvAlphaTPCTOF->Clone("fHistTRDpTvPIDvAlphaTPCTOF_Clone");
    TH1D *fHistTRDPIDAlphaTPCTOF = fHistTRDpTvPIDvAlphaTPCTOF_Clone->ProjectionY("fHistTRDPIDAlphaTPCTOF",1,260);
    fHistTRDPIDAlphaTPCTOF->SetTitle("TRD PID distribution for Alpha;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAlphaTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAlphaTPCTOF->SetLineColor(17);
    Double_t intAlphaTPCTOF = fHistTRDPIDAlphaTPCTOF->Integral(1,260);
    fHistTRDPIDAlphaTPCTOF->Scale(1/intAlphaTPCTOF);

    TH2D *fHistTRDpTvPIDvHelium3AlphaTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvHelium3AlphaTPCTOF->Clone("fHistTRDpTvPIDvHelium3AlphaTPCTOF_Clone");
    TH1D *fHistTRDPIDHelium3AlphaTPCTOF = fHistTRDpTvPIDvHelium3AlphaTPCTOF_Clone->ProjectionY("fHistTRDPIDHelium3AlphaTPCTOF",1,260);
    fHistTRDPIDHelium3AlphaTPCTOF->SetTitle("TRD PID distribution for Helium3 & Alpha;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDHelium3AlphaTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDHelium3AlphaTPCTOF->SetLineColor(4);
    Double_t intHelium3AlphaTPCTOF = fHistTRDPIDHelium3AlphaTPCTOF->Integral(1,260);
    fHistTRDPIDHelium3AlphaTPCTOF->Scale(1/intHelium3AlphaTPCTOF);

    // Allexcept
    TH2D *fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear->Clone("fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear_Clone");
    TH1D *fHistTRDPIDAllexceptDeuteronTPCTOF_clear = fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear_Clone->ProjectionY("fHistTRDPIDAllexceptDeuteronTPCTOF_clear",1,260);
    fHistTRDPIDAllexceptDeuteronTPCTOF_clear->SetTitle("TRD PID distribution for All except Deuteron;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptDeuteronTPCTOF_clear->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptDeuteronTPCTOF_clear->SetLineColor(2);
    Double_t intAllexceptDeuteronTPCTOF_clear = fHistTRDPIDAllexceptDeuteronTPCTOF_clear->Integral(1,260);
    fHistTRDPIDAllexceptDeuteronTPCTOF_clear->Scale(1/intAllexceptDeuteronTPCTOF_clear);

    TH2D *fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear->Clone("fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear_Clone");
    TH1D *fHistTRDPIDAllexceptTritonTPCTOF_clear = fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear_Clone->ProjectionY("fHistTRDPIDAllexceptTritonTPCTOF_clear",1,260);
    fHistTRDPIDAllexceptTritonTPCTOF_clear->SetTitle("TRD PID distribution for All except Triton;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptTritonTPCTOF_clear->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptTritonTPCTOF_clear->SetLineColor(2);
    Double_t intAllexceptTritonTPCTOF_clear = fHistTRDPIDAllexceptTritonTPCTOF_clear->Integral(1,260);
    fHistTRDPIDAllexceptTritonTPCTOF_clear->Scale(1/intAllexceptTritonTPCTOF_clear);

    TH2D *fHistTRDpTvPIDvAllexceptHelium3TPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptHelium3TPCTOF->Clone("fHistTRDpTvPIDvAllexceptHelium3TPCTOF_Clone");
    TH1D *fHistTRDPIDAllexceptHelium3TPCTOF = fHistTRDpTvPIDvAllexceptHelium3TPCTOF_Clone->ProjectionY("fHistTRDPIDAllexceptHelium3TPCTOF",1,260);
    fHistTRDPIDAllexceptHelium3TPCTOF->SetTitle("TRD PID distribution for All except Helium3;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptHelium3TPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptHelium3TPCTOF->SetLineColor(2);
    Double_t intAllexceptHelium3TPCTOF = fHistTRDPIDAllexceptHelium3TPCTOF->Integral(1,260);
    fHistTRDPIDAllexceptHelium3TPCTOF->Scale(1/intAllexceptHelium3TPCTOF);

    TH2D *fHistTRDpTvPIDvAllexceptAlphaTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptAlphaTPCTOF->Clone("fHistTRDpTvPIDvAllexceptAlphaTPCTOF_Clone");
    TH1D *fHistTRDPIDAllexceptAlphaTPCTOF = fHistTRDpTvPIDvAllexceptAlphaTPCTOF_Clone->ProjectionY("fHistTRDPIDAllexceptAlphaTPCTOF",1,260);
    fHistTRDPIDAllexceptAlphaTPCTOF->SetTitle("TRD PID distribution for All except Alpha;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptAlphaTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptAlphaTPCTOF->SetLineColor(2);
    Double_t intAllexceptAlphaTPCTOF = fHistTRDPIDAllexceptAlphaTPCTOF->Integral(1,260);
    fHistTRDPIDAllexceptAlphaTPCTOF->Scale(1/intAllexceptAlphaTPCTOF);

    TH2D *fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF_Clone = (TH2D *)fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF->Clone("fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF_Clone");
    TH1D *fHistTRDPIDAllexceptHelium3AlphaTPCTOF = fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF_Clone->ProjectionY("fHistTRDPIDAllexceptHelium3AlphaTPCTOF",1,260);
    fHistTRDPIDAllexceptHelium3AlphaTPCTOF->SetTitle("TRD PID distribution for All except Helium3 & Alpha;<Q_{s}> (a.u.);Normalized Yield");
    fHistTRDPIDAllexceptHelium3AlphaTPCTOF->GetYaxis()->SetTitleOffset(1.5);
    fHistTRDPIDAllexceptHelium3AlphaTPCTOF->SetLineColor(2);
    Double_t intAllexceptHelium3AlphaTPCTOF = fHistTRDPIDAllexceptHelium3AlphaTPCTOF->Integral(1,260);
    fHistTRDPIDAllexceptHelium3AlphaTPCTOF->Scale(1/intAllexceptHelium3AlphaTPCTOF);

    // proj_Clone
    TH1D *fHistTRDPIDallTPCTOF = (TH1D *)fHistTRDPIDPionTPCTOF->Clone("fHistTRDPIDallTPCTOF");
    fHistTRDPIDallTPCTOF->SetTitle("TRD PID distribution;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallDeuteronTPCTOF = (TH1D *)fHistTRDPIDAllexceptAlphaTPCTOF->Clone("fHistTRDPIDallDeuteronTPCTOF");
    fHistTRDPIDallDeuteronTPCTOF->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallTritonTPCTOF = (TH1D *)fHistTRDPIDAllexceptAlphaTPCTOF->Clone("fHistTRDPIDallTritonTPCTOF");
    fHistTRDPIDallTritonTPCTOF->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallHelium3TPCTOF = (TH1D *)fHistTRDPIDAllexceptAlphaTPCTOF->Clone("fHistTRDPIDallHelium3TPCTOF");
    fHistTRDPIDallHelium3TPCTOF->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallAlphaTPCTOF = (TH1D *)fHistTRDPIDAllexceptAlphaTPCTOF->Clone("fHistTRDPIDallAlphaTPCTOF");
    fHistTRDPIDallAlphaTPCTOF->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");
    TH1D *fHistTRDPIDallHelium3AlphaTPCTOF = (TH1D *)fHistTRDPIDAllexceptAlphaTPCTOF->Clone("fHistTRDPIDallHelium3AlphaTPCTOF");
    fHistTRDPIDallHelium3AlphaTPCTOF->SetTitle("TRD PID distributions;<Q_{s}> (a.u.);Normalized Yield");

    TH1D *fHistTRDPIDElectronTPCTOF_Clone = (TH1D *)fHistTRDPIDElectronTPCTOF->Clone("fHistTRDPIDElectronTPCTOF_Clone");
    fHistTRDPIDElectronTPCTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDPionTPCTOF_Clone = (TH1D *)fHistTRDPIDPionTPCTOF->Clone("fHistTRDPIDPionTPCTOF_Clone");
    fHistTRDPIDPionTPCTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDKaonTPCTOF_Clone = (TH1D *)fHistTRDPIDKaonTPCTOF->Clone("fHistTRDPIDKaonTPCTOF_Clone");
    fHistTRDPIDKaonTPCTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDProtonTPCTOF_Clone = (TH1D *)fHistTRDPIDProtonTPCTOF->Clone("fHistTRDPIDProtonTPCTOF_Clone");
    fHistTRDPIDProtonTPCTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDDeuteronTPCTOF_Clone = (TH1D *)fHistTRDPIDDeuteronTPCTOF->Clone("fHistTRDPIDDeuteronTPCTOF_Clone");
    fHistTRDPIDDeuteronTPCTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone = (TH1D *)fHistTRDPIDDeuteronTPCTOF_clear->Clone("fHistTRDPIDDeuteronTPCTOF_clear_Clone");
    fHistTRDPIDDeuteronTPCTOF_clear_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone2 = (TH1D *)fHistTRDPIDDeuteronTPCTOF_clear->Clone("fHistTRDPIDDeuteronTPCTOF_clear_Clone2");
    fHistTRDPIDDeuteronTPCTOF_clear_Clone2->SetLineColor(4);
    TH1D *fHistTRDPIDTritonTPCTOF_Clone = (TH1D *)fHistTRDPIDTritonTPCTOF->Clone("fHistTRDPIDTritonTPCTOF_Clone");
    fHistTRDPIDTritonTPCTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone = (TH1D *)fHistTRDPIDTritonTPCTOF_clear->Clone("fHistTRDPIDTritonTPCTOF_clear_Clone");
    fHistTRDPIDTritonTPCTOF_clear_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone2 = (TH1D *)fHistTRDPIDTritonTPCTOF_clear->Clone("fHistTRDPIDTritonTPCTOF_clear_Clone2");
    fHistTRDPIDTritonTPCTOF_clear_Clone2->SetLineColor(4);
    TH1D *fHistTRDPIDHelium3TPCTOF_Clone = (TH1D *)fHistTRDPIDHelium3TPCTOF->Clone("fHistTRDPIDHelium3TPCTOF_Clone");
    fHistTRDPIDHelium3TPCTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDHelium3TPCTOF_Clone2 = (TH1D *)fHistTRDPIDHelium3TPCTOF->Clone("fHistTRDPIDHelium3TPCTOF_Clone2");
    fHistTRDPIDHelium3TPCTOF_Clone2->SetLineColor(4);
    TH1D *fHistTRDPIDAlphaTPCTOF_Clone = (TH1D *)fHistTRDPIDAlphaTPCTOF->Clone("fHistTRDPIDAlphaTPCTOF_Clone");
    fHistTRDPIDAlphaTPCTOF_Clone->SetLineColor(1);
    TH1D *fHistTRDPIDAlphaTPCTOF_Clone2 = (TH1D *)fHistTRDPIDAlphaTPCTOF->Clone("fHistTRDPIDAlphaTPCTOF_Clone2");
    fHistTRDPIDAlphaTPCTOF_Clone2->SetLineColor(4);

  printf("117/117 histograms done\n");


  // Create legends
  // TPC legends
  TLegend *legTPC = new TLegend(.6,.45,.8,.85);
  legTPC->SetFillStyle(0);
  legTPC->SetBorderSize(0);
  legTPC->AddEntry(fHistTRDPIDElectronTPC,"Electron","l");
  legTPC->AddEntry(fHistTRDPIDPionTPC,"Pion","l");
  legTPC->AddEntry(fHistTRDPIDKaonTPC,"Kaon","l");
  legTPC->AddEntry(fHistTRDPIDProtonTPC,"Proton","l");
  legTPC->AddEntry(fHistTRDPIDDeuteronTPC,"Deuteron","l");
  legTPC->AddEntry(fHistTRDPIDTritonTPC,"Triton","l");
  legTPC->AddEntry(fHistTRDPIDHelium3TPC,"Helium3","l");
  legTPC->AddEntry(fHistTRDPIDAlphaTPC,"Alpha","l");
  TLegend *legElectronTPC = new TLegend(.5,.75,.7,.85);
  legElectronTPC->SetFillStyle(0);
  legElectronTPC->SetBorderSize(0);
  legElectronTPC->AddEntry(fHistTRDPIDElectronTPC_Clone,"Electron","l");
  TLegend *legPionTPC = new TLegend(.5,.75,.7,.85);
  legPionTPC->SetFillStyle(0);
  legPionTPC->SetBorderSize(0);
  legPionTPC->AddEntry(fHistTRDPIDPionTPC_Clone,"Pion","l");
  TLegend *legKaonTPC = new TLegend(.5,.75,.7,.85);
  legKaonTPC->SetFillStyle(0);
  legKaonTPC->SetBorderSize(0);
  legKaonTPC->AddEntry(fHistTRDPIDKaonTPC_Clone,"Kaon","l");
  TLegend *legProtonTPC = new TLegend(.5,.75,.7,.85);
  legProtonTPC->SetFillStyle(0);
  legProtonTPC->SetBorderSize(0);
  legProtonTPC->AddEntry(fHistTRDPIDProtonTPC_Clone,"Proton","l");
  TLegend *legDeuteronTPC = new TLegend(.5,.75,.7,.85);
  legDeuteronTPC->SetFillStyle(0);
  legDeuteronTPC->SetBorderSize(0);
  legDeuteronTPC->AddEntry(fHistTRDPIDDeuteronTPC_Clone,"Deuteron","l");
  TLegend *legTritonTPC = new TLegend(.5,.75,.7,.85);
  legTritonTPC->SetFillStyle(0);
  legTritonTPC->SetBorderSize(0);
  legTritonTPC->AddEntry(fHistTRDPIDTritonTPC_Clone,"Triton","l");
  TLegend *legHelium3TPC = new TLegend(.15,.75,.35,.85);
  legHelium3TPC->SetFillStyle(0);
  legHelium3TPC->SetBorderSize(0);
  legHelium3TPC->AddEntry(fHistTRDPIDHelium3TPC_Clone,"Helium3","l");
  TLegend *legAlphaTPC = new TLegend(.15,.75,.32,.85);
  legAlphaTPC->SetFillStyle(0);
  legAlphaTPC->SetBorderSize(0);
  legAlphaTPC->AddEntry(fHistTRDPIDAlphaTPC_Clone,"Alpha","l");

  // TOF legends
  TLegend *legTOF = new TLegend(.6,.45,.8,.85);
  legTOF->SetFillStyle(0);
  legTOF->SetBorderSize(0);
  legTOF->AddEntry(fHistTRDPIDElectronTOF,"Electron","l");
  legTOF->AddEntry(fHistTRDPIDPionTOF,"Pion","l");
  legTOF->AddEntry(fHistTRDPIDKaonTOF,"Kaon","l");
  legTOF->AddEntry(fHistTRDPIDProtonTOF,"Proton","l");
  legTOF->AddEntry(fHistTRDPIDDeuteronTOF,"Deuteron","l");
  legTOF->AddEntry(fHistTRDPIDTritonTOF,"Triton","l");
  legTOF->AddEntry(fHistTRDPIDHelium3TOF,"Helium3","l");
  legTOF->AddEntry(fHistTRDPIDAlphaTOF,"Alpha","l");
  TLegend *legElectronTOF = new TLegend(.5,.75,.7,.85);
  legElectronTOF->SetFillStyle(0);
  legElectronTOF->SetBorderSize(0);
  legElectronTOF->AddEntry(fHistTRDPIDElectronTOF_Clone,"Electron","l");
  TLegend *legPionTOF = new TLegend(.5,.75,.7,.85);
  legPionTOF->SetFillStyle(0);
  legPionTOF->SetBorderSize(0);
  legPionTOF->AddEntry(fHistTRDPIDPionTOF_Clone,"Pion","l");
  TLegend *legKaonTOF = new TLegend(.5,.75,.7,.85);
  legKaonTOF->SetFillStyle(0);
  legKaonTOF->SetBorderSize(0);
  legKaonTOF->AddEntry(fHistTRDPIDKaonTOF_Clone,"Kaon","l");
  TLegend *legProtonTOF = new TLegend(.5,.75,.7,.85);
  legProtonTOF->SetFillStyle(0);
  legProtonTOF->SetBorderSize(0);
  legProtonTOF->AddEntry(fHistTRDPIDProtonTOF_Clone,"Proton","l");
  TLegend *legDeuteronTOF = new TLegend(.5,.75,.7,.85);
  legDeuteronTOF->SetFillStyle(0);
  legDeuteronTOF->SetBorderSize(0);
  legDeuteronTOF->AddEntry(fHistTRDPIDDeuteronTOF_Clone,"Deuteron","l");
  TLegend *legTritonTOF = new TLegend(.5,.75,.7,.85);
  legTritonTOF->SetFillStyle(0);
  legTritonTOF->SetBorderSize(0);
  legTritonTOF->AddEntry(fHistTRDPIDTritonTOF_Clone,"Triton","l");
  TLegend *legHelium3TOF = new TLegend(.5,.75,.7,.85);
  legHelium3TOF->SetFillStyle(0);
  legHelium3TOF->SetBorderSize(0);
  legHelium3TOF->AddEntry(fHistTRDPIDHelium3TOF_Clone,"Helium3","l");
  TLegend *legAlphaTOF = new TLegend(.5,.75,.67,.85);
  legAlphaTOF->SetFillStyle(0);
  legAlphaTOF->SetBorderSize(0);
  legAlphaTOF->AddEntry(fHistTRDPIDAlphaTOF_Clone,"Alpha","l");

  // TPC & TOF legends
  TLegend *legTPCTOF = new TLegend(.6,.45,.8,.85);
  legTPCTOF->SetFillStyle(0);
  legTPCTOF->SetBorderSize(0);
  legTPCTOF->AddEntry(fHistTRDPIDElectronTPCTOF,"Electron","l");
  legTPCTOF->AddEntry(fHistTRDPIDPionTPCTOF,"Pion","l");
  legTPCTOF->AddEntry(fHistTRDPIDKaonTPCTOF,"Kaon","l");
  legTPCTOF->AddEntry(fHistTRDPIDProtonTPCTOF,"Proton","l");
  legTPCTOF->AddEntry(fHistTRDPIDDeuteronTPCTOF,"Deuteron","l");
  legTPCTOF->AddEntry(fHistTRDPIDTritonTPCTOF,"Triton","l");
  legTPCTOF->AddEntry(fHistTRDPIDHelium3TPCTOF,"Helium3","l");
  legTPCTOF->AddEntry(fHistTRDPIDAlphaTPCTOF,"Alpha","l");
  TLegend *legElectronTPCTOF = new TLegend(.5,.75,.7,.85);
  legElectronTPCTOF->SetFillStyle(0);
  legElectronTPCTOF->SetBorderSize(0);
  legElectronTPCTOF->AddEntry(fHistTRDPIDElectronTPCTOF_Clone,"Electron","l");
  TLegend *legPionTPCTOF = new TLegend(.5,.75,.7,.85);
  legPionTPCTOF->SetFillStyle(0);
  legPionTPCTOF->SetBorderSize(0);
  legPionTPCTOF->AddEntry(fHistTRDPIDPionTPCTOF_Clone,"Pion","l");
  TLegend *legKaonTPCTOF = new TLegend(.5,.75,.7,.85);
  legKaonTPCTOF->SetFillStyle(0);
  legKaonTPCTOF->SetBorderSize(0);
  legKaonTPCTOF->AddEntry(fHistTRDPIDKaonTPCTOF_Clone,"Kaon","l");
  TLegend *legProtonTPCTOF = new TLegend(.5,.75,.7,.85);
  legProtonTPCTOF->SetFillStyle(0);
  legProtonTPCTOF->SetBorderSize(0);
  legProtonTPCTOF->AddEntry(fHistTRDPIDProtonTPCTOF_Clone,"Proton","l");
  TLegend *legDeuteronTPCTOF = new TLegend(.5,.75,.7,.85);
  legDeuteronTPCTOF->SetFillStyle(0);
  legDeuteronTPCTOF->SetBorderSize(0);
  legDeuteronTPCTOF->AddEntry(fHistTRDPIDDeuteronTPCTOF_Clone,"Deuteron","l");
  TLegend *legTritonTPCTOF = new TLegend(.5,.75,.7,.85);
  legTritonTPCTOF->SetFillStyle(0);
  legTritonTPCTOF->SetBorderSize(0);
  legTritonTPCTOF->AddEntry(fHistTRDPIDTritonTPCTOF_Clone,"Triton","l");
  TLegend *legHelium3TPCTOF = new TLegend(.15,.75,.35,.85);
  legHelium3TPCTOF->SetFillStyle(0);
  legHelium3TPCTOF->SetBorderSize(0);
  legHelium3TPCTOF->AddEntry(fHistTRDPIDHelium3TPCTOF_Clone,"Helium3","l");
  TLegend *legAlphaTPCTOF = new TLegend(.15,.75,.32,.85);
  legAlphaTPCTOF->SetFillStyle(0);
  legAlphaTPCTOF->SetBorderSize(0);
  legAlphaTPCTOF->AddEntry(fHistTRDPIDAlphaTPCTOF_Clone,"Alpha","l");

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create List, add all histograms and save in a .root file
  ///___________________________________
  TList *fOutputContainer1 = new TList();
  // PID vs. pT
  fOutputContainer1->Add(fHistTRDpTvPID);
  // TPC preselection
  fOutputContainer1->Add(fHistTRDpTvPIDvElectronTPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvPionTPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvKaonTPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvProtonTPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvDeuteronTPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvDeuteronTPC_clear);
  fOutputContainer1->Add(fHistTRDpTvPIDvTritonTPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvTritonTPC_clear);
  fOutputContainer1->Add(fHistTRDpTvPIDvHelium3TPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvAlphaTPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvHelium3AlphaTPC);

  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptDeuteronTPC_clear);
  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptTritonTPC_clear);
  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptHelium3TPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptAlphaTPC);
  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptHelium3AlphaTPC);
  // TOF preselection
  fOutputContainer1->Add(fHistTRDpTvPIDvElectronTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvPionTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvKaonTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvProtonTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvDeuteronTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvTritonTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvHelium3TOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvAlphaTOF);
  // TPC & TOF preseletion
  fOutputContainer1->Add(fHistTRDpTvPIDvElectronTPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvPionTPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvKaonTPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvProtonTPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvDeuteronTPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvDeuteronTPCTOF_clear);
  fOutputContainer1->Add(fHistTRDpTvPIDvTritonTPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvTritonTPCTOF_clear);
  fOutputContainer1->Add(fHistTRDpTvPIDvHelium3TPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvAlphaTPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvHelium3AlphaTPCTOF);

  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear);
  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear);
  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptHelium3TPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptAlphaTPCTOF);
  fOutputContainer1->Add(fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF);

  TList *fOutputContainer2 = new TList();
  // Proj PID vs. pT
  // TPC preselection
  fOutputContainer2->Add(fHistTRDPIDElectronTPC);
  fOutputContainer2->Add(fHistTRDPIDPionTPC);
  fOutputContainer2->Add(fHistTRDPIDKaonTPC);
  fOutputContainer2->Add(fHistTRDPIDProtonTPC);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPC);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPC_clear);
  fOutputContainer2->Add(fHistTRDPIDTritonTPC);
  fOutputContainer2->Add(fHistTRDPIDTritonTPC_clear);
  fOutputContainer2->Add(fHistTRDPIDHelium3TPC);
  fOutputContainer2->Add(fHistTRDPIDAlphaTPC);
  fOutputContainer2->Add(fHistTRDPIDHelium3AlphaTPC);

  fOutputContainer2->Add(fHistTRDPIDAllexceptDeuteronTPC_clear);
  fOutputContainer2->Add(fHistTRDPIDAllexceptTritonTPC_clear);
  fOutputContainer2->Add(fHistTRDPIDAllexceptHelium3TPC);
  fOutputContainer2->Add(fHistTRDPIDAllexceptAlphaTPC);
  fOutputContainer2->Add(fHistTRDPIDAllexceptHelium3AlphaTPC);

  fOutputContainer2->Add(fHistTRDPIDallTPC);
  fOutputContainer2->Add(fHistTRDPIDallDeuteronTPC);
  fOutputContainer2->Add(fHistTRDPIDallTritonTPC);
  fOutputContainer2->Add(fHistTRDPIDallHelium3TPC);
  fOutputContainer2->Add(fHistTRDPIDallAlphaTPC);
  fOutputContainer2->Add(fHistTRDPIDallHelium3AlphaTPC);

  fOutputContainer2->Add(fHistTRDPIDElectronTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDPionTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDKaonTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDProtonTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPC_clear_Clone);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPC_clear_Clone2);
  fOutputContainer2->Add(fHistTRDPIDTritonTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDTritonTPC_clear_Clone);
  fOutputContainer2->Add(fHistTRDPIDTritonTPC_clear_Clone2);
  fOutputContainer2->Add(fHistTRDPIDHelium3TPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDHelium3TPC_Clone2);
  fOutputContainer2->Add(fHistTRDPIDAlphaTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDAlphaTPC_Clone2);
  // TOF preselection
  fOutputContainer2->Add(fHistTRDPIDElectronTOF);
  fOutputContainer2->Add(fHistTRDPIDPionTOF);
  fOutputContainer2->Add(fHistTRDPIDKaonTOF);
  fOutputContainer2->Add(fHistTRDPIDProtonTOF);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTOF);
  fOutputContainer2->Add(fHistTRDPIDTritonTOF);
  fOutputContainer2->Add(fHistTRDPIDHelium3TOF);
  fOutputContainer2->Add(fHistTRDPIDAlphaTOF);

  fOutputContainer2->Add(fHistTRDPIDallTOF);

  fOutputContainer2->Add(fHistTRDPIDElectronTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDPionTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDKaonTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDProtonTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDTritonTPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDHelium3TPC_Clone);
  fOutputContainer2->Add(fHistTRDPIDAlphaTPC_Clone);
  // TPC & TOF preselection
  fOutputContainer2->Add(fHistTRDPIDElectronTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDPionTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDKaonTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDProtonTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPCTOF_clear);
  fOutputContainer2->Add(fHistTRDPIDTritonTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDTritonTPCTOF_clear);
  fOutputContainer2->Add(fHistTRDPIDHelium3TPCTOF);
  fOutputContainer2->Add(fHistTRDPIDAlphaTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDHelium3AlphaTPCTOF);

  fOutputContainer2->Add(fHistTRDPIDAllexceptDeuteronTPCTOF_clear);
  fOutputContainer2->Add(fHistTRDPIDAllexceptTritonTPCTOF_clear);
  fOutputContainer2->Add(fHistTRDPIDAllexceptHelium3TPCTOF);
  fOutputContainer2->Add(fHistTRDPIDAllexceptAlphaTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDAllexceptHelium3AlphaTPCTOF);

  fOutputContainer2->Add(fHistTRDPIDallTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDallDeuteronTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDallTritonTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDallHelium3TPCTOF);
  fOutputContainer2->Add(fHistTRDPIDallAlphaTPCTOF);
  fOutputContainer2->Add(fHistTRDPIDallHelium3AlphaTPCTOF);

  fOutputContainer2->Add(fHistTRDPIDElectronTPCTOF_Clone);
  fOutputContainer2->Add(fHistTRDPIDPionTPCTOF_Clone);
  fOutputContainer2->Add(fHistTRDPIDKaonTPCTOF_Clone);
  fOutputContainer2->Add(fHistTRDPIDProtonTPCTOF_Clone);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPCTOF_Clone);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPCTOF_clear_Clone);
  fOutputContainer2->Add(fHistTRDPIDDeuteronTPCTOF_clear_Clone2);
  fOutputContainer2->Add(fHistTRDPIDTritonTPCTOF_Clone);
  fOutputContainer2->Add(fHistTRDPIDTritonTPCTOF_clear_Clone);
  fOutputContainer2->Add(fHistTRDPIDTritonTPCTOF_clear_Clone2);
  fOutputContainer2->Add(fHistTRDPIDHelium3TPCTOF_Clone);
  fOutputContainer2->Add(fHistTRDPIDHelium3TPCTOF_Clone2);
  fOutputContainer2->Add(fHistTRDPIDAlphaTPCTOF_Clone);
  fOutputContainer2->Add(fHistTRDPIDAlphaTPCTOF_Clone2);

  TList *fOutputContainer3 = new TList();
  // dE/dx vs. Rigidity
  fOutputContainer3->Add(fHistdEdxvRig);

  fOutputContainer3->Add(fHistdEdxvRigvElectronTPC);
  fOutputContainer3->Add(fHistdEdxvRigvPionTPC);
  fOutputContainer3->Add(fHistdEdxvRigvKaonTPC);
  fOutputContainer3->Add(fHistdEdxvRigvProtonTPC);
  fOutputContainer3->Add(fHistdEdxvRigvDeuteronTPC);
  fOutputContainer3->Add(fHistdEdxvRigvTritonTPC);
  fOutputContainer3->Add(fHistdEdxvRigvHelium3TPC);
  fOutputContainer3->Add(fHistdEdxvRigvHelium3TPC_bb);
  fOutputContainer3->Add(fHistdEdxvRigvAlphaTPC);
  fOutputContainer3->Add(fHistdEdxvRigvAlphaTPC_bb);

  fOutputContainer3->Add(fHistdEdxvRigvElectronTOF);
  fOutputContainer3->Add(fHistdEdxvRigvPionTOF);
  fOutputContainer3->Add(fHistdEdxvRigvKaonTOF);
  fOutputContainer3->Add(fHistdEdxvRigvProtonTOF);
  fOutputContainer3->Add(fHistdEdxvRigvDeuteronTOF);
  fOutputContainer3->Add(fHistdEdxvRigvTritonTOF);
  fOutputContainer3->Add(fHistdEdxvRigvHelium3TOF);
  fOutputContainer3->Add(fHistdEdxvRigvAlphaTOF);

  fOutputContainer3->Add(fHistdEdxvRigvElectronTPCTOF);
  fOutputContainer3->Add(fHistdEdxvRigvPionTPCTOF);
  fOutputContainer3->Add(fHistdEdxvRigvKaonTPCTOF);
  fOutputContainer3->Add(fHistdEdxvRigvProtonTPCTOF);
  fOutputContainer3->Add(fHistdEdxvRigvDeuteronTPCTOF);
  fOutputContainer3->Add(fHistdEdxvRigvTritonTPCTOF);
  fOutputContainer3->Add(fHistdEdxvRigvHelium3TPCTOF);
  fOutputContainer3->Add(fHistdEdxvRigvHelium3TPCTOF_bb);
  fOutputContainer3->Add(fHistdEdxvRigvAlphaTPCTOF);
  fOutputContainer3->Add(fHistdEdxvRigvAlphaTPCTOF_bb);


  TFile *coutput = 0;
  if(LessStatisticFactor == 1) {
    coutput = new TFile(Form("/u/brudnyj/Desktop/Bachelor/codes/%s/bbrudnyj_ReadTreeTestTRDTrigger2/bbrudnyj_TestTRDTrigger2%s%iSigma.root",Collisions,Collisions,Sigma),"UPDATE");
    coutput->WriteTObject(fOutputContainer1,Form("bbrudnyj_TestTRDTrigger2%s%iSigmaPID",Collisions,Sigma),"kOverwrite");
    coutput->WriteTObject(fOutputContainer2,Form("bbrudnyj_TestTRDTrigger2%s%iSigmaProj",Collisions,Sigma),"kOverwrite");
    coutput->WriteTObject(fOutputContainer3,Form("bbrudnyj_TestTRDTrigger2%s%iSigmadEdx",Collisions,Sigma),"kOverwrite");
  }
  else {
    coutput = new TFile(Form("/u/brudnyj/Desktop/Bachelor/codes/%s/bbrudnyj_ReadTreeTestTRDTrigger2/bbrudnyj_TestTRDTrigger2%s%iSigma1:%i.root",Collisions,Collisions,Sigma,LessStatisticFactor),"UPDATE");
    coutput->WriteTObject(fOutputContainer1,Form("bbrudnyj_TestTRDTrigger2%s%iSigma1:%iPID",Collisions,Sigma,LessStatisticFactor),"kOverwrite");
    coutput->WriteTObject(fOutputContainer2,Form("bbrudnyj_TestTRDTrigger2%s%iSigma1:%iProj",Collisions,Sigma,LessStatisticFactor),"kOverwrite");
    coutput->WriteTObject(fOutputContainer3,Form("bbrudnyj_TestTRDTrigger2%s%iSigma1:%idEdx",Collisions,Sigma,LessStatisticFactor),"kOverwrite");
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create and fill canvases
  ///_____________________________________________
  printf("Creating and saving 9 canvases...\n\n");
  // PID vs. pT
    TCanvas *c1 = new TCanvas("c1",Form("PIDpT w TPC preselection (Sigma=%i) for separated particles",Sigma));
    c1->Divide(3,3);
    c1->cd(1)->SetLogz();
      fHistTRDpTvPID->Draw("colz");
    c1->cd(2)->SetLogz();
      fHistTRDpTvPIDvElectronTPC->Draw("colz");
    c1->cd(3)->SetLogz();
      fHistTRDpTvPIDvPionTPC->Draw("colz");
    c1->cd(4)->SetLogz();
      fHistTRDpTvPIDvKaonTPC->Draw("colz");
    c1->cd(5)->SetLogz();
      fHistTRDpTvPIDvProtonTPC->Draw("colz");
    c1->cd(6)->SetLogz();
    //fHistTRDpTvPIDvDeuteronTPC->Draw("colz");
      fHistTRDpTvPIDvDeuteronTPC_clear->Draw("colz");
    c1->cd(7)->SetLogz();
    //fHistTRDpTvPIDvTritonTPC->Draw("colz");
      fHistTRDpTvPIDvTritonTPC_clear->Draw("colz");
    c1->cd(8)->SetLogz();
      fHistTRDpTvPIDvHelium3TPC->Draw("colz");
    c1->cd(9)->SetLogz();
      fHistTRDpTvPIDvAlphaTPC->Draw("colz");

    TCanvas *c2 = new TCanvas("c2",Form("PIDpT w TOF preselection (Sigma=%i) for separated particles",Sigma));
    c2->Divide(3,3);
    c2->cd(1)->SetLogz();
      fHistTRDpTvPID->Draw("colz");
    c2->cd(2)->SetLogz();
      fHistTRDpTvPIDvElectronTOF->Draw("colz");
    c2->cd(3)->SetLogz();
      fHistTRDpTvPIDvPionTOF->Draw("colz");
    c2->cd(4)->SetLogz();
      fHistTRDpTvPIDvKaonTOF->Draw("colz");
    c2->cd(5)->SetLogz();
      fHistTRDpTvPIDvProtonTOF->Draw("colz");
    c2->cd(6)->SetLogz();
      fHistTRDpTvPIDvDeuteronTOF->Draw("colz");
    c2->cd(7)->SetLogz();
      fHistTRDpTvPIDvTritonTOF->Draw("colz");
    c2->cd(8)->SetLogz();
      fHistTRDpTvPIDvHelium3TOF->Draw("colz");
    c2->cd(9)->SetLogz();
      fHistTRDpTvPIDvAlphaTOF->Draw("colz");

    TCanvas *c3 = new TCanvas("c3",Form("PIDpT w TPC & TOF preselection (Sigma=%i) for separated particles",Sigma));
    c3->Divide(3,3);
    c3->cd(1)->SetLogz();
      fHistTRDpTvPID->Draw("colz");
    c3->cd(2)->SetLogz();
      fHistTRDpTvPIDvElectronTPCTOF->Draw("colz");
    c3->cd(3)->SetLogz();
      fHistTRDpTvPIDvPionTPCTOF->Draw("colz");
    c3->cd(4)->SetLogz();
      fHistTRDpTvPIDvKaonTPCTOF->Draw("colz");
    c3->cd(5)->SetLogz();
      fHistTRDpTvPIDvProtonTPCTOF->Draw("colz");
    c3->cd(6)->SetLogz();
    //fHistTRDpTvPIDvDeuteronTPCTOF->Draw("colz");
      fHistTRDpTvPIDvDeuteronTPCTOF_clear->Draw("colz");
    c3->cd(7)->SetLogz();
    //fHistTRDpTvPIDvTritonTPCTOF->Draw("colz");
      fHistTRDpTvPIDvTritonTPCTOF_clear->Draw("colz");
    c3->cd(8)->SetLogz();
      fHistTRDpTvPIDvHelium3TPCTOF->Draw("colz");
    c3->cd(9)->SetLogz();
      fHistTRDpTvPIDvAlphaTPCTOF->Draw("colz");

  // PID vs. pT projection histograms
    TCanvas *c4 = new TCanvas("c4",Form("Proj PIDpT w TPC preselection (Sigma=%i) for separated particles",Sigma));
    c4->Divide(3,3);
    c4->cd(1);
      fHistTRDPIDallTPC->Draw(); // clone of fHistTRDPIDPionTPC because of most entries
      fHistTRDPIDAlphaTPC->Draw("same");
      fHistTRDPIDElectronTPC->Draw("same");
      //fHistTRDPIDDeuteronTPC->Draw("same");
      fHistTRDPIDDeuteronTPC_clear->Draw("same");
      //fHistTRDPIDTritonTPC->Draw("same");
      fHistTRDPIDTritonTPC_clear->Draw("same");
      fHistTRDPIDProtonTPC->Draw("same");
      fHistTRDPIDKaonTPC->Draw("same");
      fHistTRDPIDHelium3TPC->Draw("same");
      legTPC->Draw("same");
    c4->cd(2);
      fHistTRDPIDElectronTPC_Clone->Draw();
      legElectronTPC->Draw("same");
    c4->cd(3);
      fHistTRDPIDPionTPC_Clone->Draw();
      legPionTPC->Draw("same");
    c4->cd(4);
      fHistTRDPIDKaonTPC_Clone->Draw();
      legKaonTPC->Draw("same");
    c4->cd(5);
      fHistTRDPIDProtonTPC_Clone->Draw();
      legProtonTPC->Draw("same");
    c4->cd(6);
    //fHistTRDPIDDeuteronTPC_Clone->Draw();
      fHistTRDPIDDeuteronTPC_clear_Clone->Draw();
      legDeuteronTPC->Draw("same");
    c4->cd(7);
    //fHistTRDPIDTritonTPC_Clone->Draw();
      fHistTRDPIDTritonTPC_clear_Clone->Draw();
      legTritonTPC->Draw("same");
    c4->cd(8);
      fHistTRDPIDHelium3TPC_Clone->Draw();
      legHelium3TPC->Draw("same");
    c4->cd(9);
      fHistTRDPIDAlphaTPC_Clone->Draw();
      legAlphaTPC->Draw("same");

    TCanvas *c5 = new TCanvas("c5",Form("Proj PIDpT w TOF preselection (Sigma=%i) for separated particles",Sigma));
    c5->Divide(3,3);
    c5->cd(1);
      fHistTRDPIDallTOF->Draw(); // clone of fHistTRDPIDPionTOF because of most entries
      fHistTRDPIDAlphaTOF->Draw("same");
      fHistTRDPIDElectronTOF->Draw("same");
      fHistTRDPIDDeuteronTOF->Draw("same");
      fHistTRDPIDTritonTOF->Draw("same");
      fHistTRDPIDProtonTOF->Draw("same");
      fHistTRDPIDKaonTOF->Draw("same");
      fHistTRDPIDHelium3TOF->Draw("same");
      legTOF->Draw("same");
    c5->cd(2);
      fHistTRDPIDElectronTOF_Clone->Draw();
      legElectronTOF->Draw("same");
    c5->cd(3);
      fHistTRDPIDPionTOF_Clone->Draw();
      legPionTOF->Draw("same");
    c5->cd(4);
      fHistTRDPIDKaonTOF_Clone->Draw();
      legKaonTOF->Draw("same");
    c5->cd(5);
      fHistTRDPIDProtonTOF_Clone->Draw();
      legProtonTOF->Draw("same");
    c5->cd(6);
      fHistTRDPIDDeuteronTOF_Clone->Draw();
      legDeuteronTOF->Draw("same");
    c5->cd(7);
      fHistTRDPIDTritonTOF_Clone->Draw();
      legTritonTOF->Draw("same");
    c5->cd(8);
      fHistTRDPIDHelium3TOF_Clone->Draw();
      legHelium3TOF->Draw("same");
    c5->cd(9);
      fHistTRDPIDAlphaTOF_Clone->Draw();
      legAlphaTOF->Draw("same");

    TCanvas *c6 = new TCanvas("c6",Form("Proj PIDpT w TPC & TOF preselection (Sigma=%i) for separated particles",Sigma));
    c6->Divide(3,3);
    c6->cd(1);
      fHistTRDPIDallTPCTOF->Draw(); // clone of fHistTRDPIDPionTPCTOF because of most entries
      fHistTRDPIDAlphaTPCTOF->Draw("same");
      fHistTRDPIDElectronTPCTOF->Draw("same");
      //fHistTRDPIDDeuteronTPCTOF->Draw("same");
      fHistTRDPIDDeuteronTPCTOF_clear->Draw("same");
      //fHistTRDPIDTritonTPCTOF->Draw("same");
      fHistTRDPIDTritonTPCTOF_clear->Draw("same");
      fHistTRDPIDProtonTPCTOF->Draw("same");
      fHistTRDPIDKaonTPCTOF->Draw("same");
      fHistTRDPIDHelium3TPCTOF->Draw("same");
      legTPCTOF->Draw("same");
    c6->cd(2);
      fHistTRDPIDElectronTPCTOF_Clone->Draw();
      legElectronTPCTOF->Draw("same");
    c6->cd(3);
      fHistTRDPIDPionTPCTOF_Clone->Draw();
      legPionTPCTOF->Draw("same");
    c6->cd(4);
      fHistTRDPIDKaonTPCTOF_Clone->Draw();
      legKaonTPCTOF->Draw("same");
    c6->cd(5);
      fHistTRDPIDProtonTPCTOF_Clone->Draw();
      legProtonTPCTOF->Draw("same");
    c6->cd(6);
    //fHistTRDPIDDeuteronTPCTOF_Clone->Draw();
      fHistTRDPIDDeuteronTPCTOF_clear_Clone->Draw();
      legDeuteronTPCTOF->Draw("same");
    c6->cd(7);
    //fHistTRDPIDTritonTPCTOF_Clone->Draw();
      fHistTRDPIDTritonTPCTOF_clear_Clone->Draw();
      legTritonTPCTOF->Draw("same");
    c6->cd(8);
      fHistTRDPIDHelium3TPCTOF_Clone->Draw();
      legHelium3TPCTOF->Draw("same");
    c6->cd(9);
      fHistTRDPIDAlphaTPCTOF_Clone->Draw();
      legAlphaTPCTOF->Draw("same");

  // dE/dx vs. Rigidity
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
  ///____________________________________________________________________________________________
      /*  c1 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_PIDTPC%iSigma%s.pdf",Sigma,Sigma,Collisions));
  c2 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_PIDTOF%iSigma%s.pdf",Sigma,Sigma,Collisions));
  c3 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_PIDTPCTOF%iSigma%s.pdf",Sigma,Sigma,Collisions));
  c4 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_ProjPIDTPC%iSigma%s.pdf",Sigma,Sigma,Collisions));
  c5 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_ProjPIDTOF%iSigma%s.pdf",Sigma,Sigma,Collisions));
  c6 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_ProjPIDTPCTOF%iSigma%s.pdf",Sigma,Sigma,Collisions));
  c7 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_dEdxTPC%iSigma%s.pdf",Sigma,Sigma,Collisions));
  c8 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_dEdxTOF%iSigma%s.pdf",Sigma,Sigma,Collisions));
  c9 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma/bbrudnyj_dEdxTPCTOF%iSigma%s.pdf",Sigma,Sigma,Collisions));
*/
  printf("\n\ndone!\n");
}
