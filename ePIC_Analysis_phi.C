// Run this macro, from the Linux / Terminal command line, as
// root -l -q ePIC_Analysis.C
// The only condition is that there is a subdirectory called "data" that contains the input MC file
// (provided by Caroline on Box / copied from SDCC)
// The name of this input MC file (variable "strang") is hardcoded in this macro and must match with the file name
// CR 2024-08-14/20

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <iostream>

void ePIC_Analysis_phi(){

  gSystem->Exec("date");
  
  // define nHCal acceptance (2024-07-16) - make sure this is consistent with your *Plotting.C macros:
  double eta_min_nHCal = -4.05;
  double eta_max_nHCal = -1.2;
  double eta_min_barrel = -1.2;
  double eta_max_barrel = 1.18;
  double eta_min_LFHCal = 1.18;
  double eta_max_LFHCal = 4.2;
  //original: -4.14 to -1.18
  //nHCal: -4.05 to -1.2
  //Barrel: -1.2 to 1.18
  //LFHCal: 1.18 to 4.2
  //
  
  ////////////////////////////////////////////////////
  //// String definitions - modify here as needed ////
  /////////////////////////////////////////////////////
  
  // Define input directory if reading locally:
  const char indir[]="data";
  cout << "Input directory is: " << indir << " \n";
    
  // Define name of local input MC file:
  const char strang[]="sartre_bnonsat_Au_phi_ab_eAu_1.3998.eicrecon.tree.edm4eic"; // "strang = data string"
  
  // define flavor of this macro:  
  const char flavor[]="ePIC";
  
  ////////////////////////////////////
  //// end of string definitions ////
  ///////////////////////////////////

  // If reading from a locally copied file:
  TString infile_ram=indir + TString("/") + strang + TString(".root");  // remains in RAM
  const char *infile=infile_ram.Data(); // points to a valid address in RAM
  cout << "Analyzed MC file will be: " << infile << " \n";
  

  TString outfile_ram= TString("out.") + strang + TString("-") + flavor + TString(".root");
  const char *outfile=outfile_ram.Data();
  TFile *ofile = TFile::Open(outfile,"RECREATE"); // RECREATE overwrites an existing file of the same name
  
  TChain *mychain = new TChain("events");

  // if reading a single file:
  mychain->Add(infile);
  
  ///////////////////////////////////////
  //// end of automated definitions ////
  //////////////////////////////////////  

  // Initialize reader
  TTreeReader tree_reader(mychain);
  
  
  //TTreeReaderArray<data type> name(tree_reader, "file name in input file");
  
  //Get Event-Level Information
  TTreeReaderArray<float> xBjorken(tree_reader, "InclusiveKinematicsTruth.x");
  TTreeReaderArray<float> QSquared(tree_reader, "InclusiveKinematicsTruth.Q2");
  
  
		

  // Get Particle Information
  TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
  TTreeReaderArray<float> partMomX(tree_reader, "MCParticles.momentum.x");
  TTreeReaderArray<float> partMomY(tree_reader, "MCParticles.momentum.y");
  TTreeReaderArray<float> partMomZ(tree_reader, "MCParticles.momentum.z");
  TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");

  // Get Reconstructed Track Information
  TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
  TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
  TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");

  // Get Associations Between MCParticles and ReconstructedChargedParticles
  TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
  TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");

  // Get parent and daugther information
  TTreeReaderArray<int> parents_index(tree_reader, "_MCParticles_parents.index");
  TTreeReaderArray<unsigned int> parents_begin(tree_reader, "MCParticles.parents_begin");
  TTreeReaderArray<unsigned int> parents_end(tree_reader, "MCParticles.parents_end");
  
  TTreeReaderArray<int> daughters_index(tree_reader, "_MCParticles_daughters.index");
  TTreeReaderArray<unsigned int> daughters_begin(tree_reader, "MCParticles.daughters_begin");
  TTreeReaderArray<unsigned int> daughters_end(tree_reader, "MCParticles.daughters_end");  
  
  // Define Histograms

//TH2D *File name = new TH2D("File name" , "plot title ; x-axis ; y_axis", x bins, left, right, y bins, lower, upper);

  //generatorStatus
  TH1D *generatorStatus = new TH1D("generatorStatus","Status of generated particles, all; generatorStatus",101,0,100);
  
  
  //Event-Level Information
  
  
  int xBbins = 300; //Set number of bins for xB
  int Qbins = 300; //Set number of bins for Q Squared
  TH1D *xBjorkenPlot = new TH1D("xBjorkenPlot", "xBjorken ; xBjorken" , 300, 0., .1);
   //xBjorkenPlot->SetLogx(1);
  TH1D *QSquaredPlot = new TH1D("QSquaredPlot", "QSquared ; Q2" , 300, 1., 20.);
   //QSquaredPlot->SetLogx(1);
   
   
  TH1D *Acceptance_VS_xBjorken_nHCal_0 = new TH1D("Acceptance_VS_xBjorken_nHCal_0","Acceptance vs x_Bjorken; x_Bjorken; nHCal Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_xBjorken_nHCal_1 = new TH1D("Acceptance_VS_xBjorken_nHCal_1","Acceptance vs x_Bjorken; x_Bjorken; nHCal Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_xBjorken_nHCal_2 = new TH1D("Acceptance_VS_xBjorken_nHCal_2","Acceptance vs x_Bjorken; x_Bjorken; nHCal Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_QSquared_nHCal_0 = new TH1D("Acceptance_VS_QSquared_nHCal_0","Acceptance vs QSquared; QSquared; nHCal Acceptance (%)", 300, 1, 20);
  TH1D *Acceptance_VS_QSquared_nHCal_1 = new TH1D("Acceptance_VS_QSquared_nHCal_1","Acceptance vs QSquared; QSquared; nHCal Acceptance (%)", 300, 1, 20);
  TH1D *Acceptance_VS_QSquared_nHCal_2 = new TH1D("Acceptance_VS_QSquared_nHCal_2","Acceptance vs QSquared; QSquared; nHCal Acceptance (%)", 300, 1, 20);
  TH2D *Acceptance_VS_xBjorken_nHCal = new TH2D("Acceptance_VS_xBjorken_nHCal","Acceptance vs x_Bjorken; x_Bjorken; nHCal Acceptance (%)", 300, 0, .1, 120, 0, 120);
  TH2D *Acceptance_VS_QSquared_nHCal = new TH2D("Acceptance_VS_QSquared_nHCal","Acceptance vs Q Squared; Q Squared; nHCal Acceptance (%)", 300, 0, 20, 120, 0, 120);
  
  
  
  TH1D *Acceptance_VS_xBjorken_barrel_0 = new TH1D("Acceptance_VS_xBjorken_barrel_0","Acceptance vs x_Bjorken; x_Bjorken; barrel Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_xBjorken_barrel_1 = new TH1D("Acceptance_VS_xBjorken_barrel_1","Acceptance vs x_Bjorken; x_Bjorken; barrel Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_xBjorken_barrel_2 = new TH1D("Acceptance_VS_xBjorken_barrel_2","Acceptance vs x_Bjorken; x_Bjorken; barrel Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_QSquared_barrel_0 = new TH1D("Acceptance_VS_QSquared_barrel_0","Acceptance vs QSquared; QSquared; barrel Acceptance (%)", 300, 1, 20);
  TH1D *Acceptance_VS_QSquared_barrel_1 = new TH1D("Acceptance_VS_QSquared_barrel_1","Acceptance vs QSquared; QSquared; barrel Acceptance (%)", 300, 1, 20);
  TH1D *Acceptance_VS_QSquared_barrel_2 = new TH1D("Acceptance_VS_QSquared_barrel_2","Acceptance vs QSquared; QSquared; barrel Acceptance (%)", 300, 1, 20);
  TH2D *Acceptance_VS_xBjorken_barrel = new TH2D("Acceptance_VS_xBjorken_barrel","Acceptance vs x_Bjorken; x_Bjorken; barrel Acceptance (%)", 300, 0, .1, 120, 0, 120);
  TH2D *Acceptance_VS_QSquared_barrel = new TH2D("Acceptance_VS_QSquared_barrel","Acceptance vs Q Squared; Q Squared; barrel Acceptance (%)", 300, 0, 20, 100, 0, 100);
  
  
  
  
  
  TH1D *Acceptance_VS_xBjorken_LFHCal_0 = new TH1D("Acceptance_VS_xBjorken_LFHCal_0","Acceptance vs x_Bjorken; x_Bjorken; LFHCal Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_xBjorken_LFHCal_1 = new TH1D("Acceptance_VS_xBjorken_LFHCal_1","Acceptance vs x_Bjorken; x_Bjorken; LFHCal Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_xBjorken_LFHCal_2 = new TH1D("Acceptance_VS_xBjorken_LFHCal_2","Acceptance vs x_Bjorken; x_Bjorken; LFHCal Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_QSquared_LFHCal_0 = new TH1D("Acceptance_VS_QSquared_LFHCal_0","Acceptance vs QSquared; QSquared; LFHCal Acceptance (%)", 300, 1, 20);
  TH1D *Acceptance_VS_QSquared_LFHCal_1 = new TH1D("Acceptance_VS_QSquared_LFHCal_1","Acceptance vs QSquared; QSquared; LFHCal Acceptance (%)", 300, 1, 20);
  TH1D *Acceptance_VS_QSquared_LFHCal_2 = new TH1D("Acceptance_VS_QSquared_LFHCal_2","Acceptance vs QSquared; QSquared; LFHCal Acceptance (%)", 300, 1, 20);
  TH2D *Acceptance_VS_xBjorken_LFHCal = new TH2D("Acceptance_VS_xBjorken_LFHCal","Acceptance vs x_Bjorken; x_Bjorken; LFHCal Acceptance (%)", 300, 0, .1, 120, 0, 120);
  TH2D *Acceptance_VS_QSquared_LFHCal = new TH2D("Acceptance_VS_QSquared_LFHCal","Acceptance vs Q Squared; Q Squared; LFHCal Acceptance (%)", 300, 0, 20, 120, 0, 120);
  
  

  
  TH1D *Acceptance_VS_xBjorken_any_0 = new TH1D("Acceptance_VS_xBjorken_any_0","Acceptance vs x_Bjorken; x_Bjorken; Any HCal Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_xBjorken_any_1 = new TH1D("Acceptance_VS_xBjorken_any_1","Acceptance vs x_Bjorken; x_Bjorken; Any HCal Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_xBjorken_any_2 = new TH1D("Acceptance_VS_xBjorken_any_2","Acceptance vs x_Bjorken; x_Bjorken; Any HCal Acceptance (%)", 300, 0, .1);
  TH1D *Acceptance_VS_QSquared_any_0 = new TH1D("Acceptance_VS_QSquared_any_0","Acceptance vs QSquared; QSquared; Any HCal Acceptance (%)", 300, 1, 20);
  TH1D *Acceptance_VS_QSquared_any_1 = new TH1D("Acceptance_VS_QSquared_any_1","Acceptance vs QSquared; QSquared; Any HCal Acceptance (%)", 300, 1, 20);
  TH1D *Acceptance_VS_QSquared_any_2 = new TH1D("Acceptance_VS_QSquared_any_2","Acceptance vs QSquared; QSquared; Any HCal Acceptance (%)", 300, 1, 20);
  TH2D *Acceptance_VS_xBjorken_any = new TH2D("Acceptance_VS_xBjorken_any","Acceptance vs x_Bjorken; x_Bjorken; Any HCal Acceptance (%)", 300, 0, .1, 120, 0, 120);
  TH2D *Acceptance_VS_QSquared_any = new TH2D("Acceptance_VS_QSquared_any","Acceptance vs Q Squared; Q Squared; Any HCal Acceptance (%)", 300, 0, 20, 120, 0, 120);
  
  
  TH2D *QSquared_VS_xBjorken = new TH2D("QSquared_VS_xBjorken","QSquared vs xBjorken; x_Bjorken; Q Squared", 300, 0, .1, 300, 1., 20);
   QSquared_VS_xBjorken->SetStats(0);
   //QSquared_VS_xBjorken->SetLogx(1);
  
  // eta (pseudorapidity)
  TH1D *partEta = new TH1D("partEta","Eta of thrown particles; #eta",120,-6.,6.);
  TH1D *recEta = new TH1D("recEta","Eta of reconstructed tracks that have matching thrown particle; #eta",120,-6.,6.);

  TH1D *electronEta = new TH1D("electronEta", "Eta of thrown e; #eta",120,-6.,6.);
  TH1D *electronRecEta = new TH1D("electronRecEta","Eta of reco e;#eta",120,-6.,6.);
  
  TH1D *protonEta = new TH1D("protonEta","Eta of thrown p;#eta",120,-6.,6.);
  TH1D *protonRecEta = new TH1D("protonRecEta","Eta of reco p;#eta",120,-6.,6.);
  
  TH1D *muonEta = new TH1D("muonEta","Eta of thrown #mu;#eta",120,-6.,6.);
  TH1D *muonRecEta = new TH1D("muonRecEta","Eta of reco #mu;#eta",120,-6.,6.);
  
  TH1D *pionEta = new TH1D("pionEta","Eta of thrown #pi;#eta",120,-6.,6.);
  TH1D *pionRecEta = new TH1D("pionRecEta","Eta of reco #pi;#eta",120,-6.,6.);

  TH1D *pi0Eta = new TH1D("pi0Eta","Eta of thrown #pi;#eta",120,-6.,6.);
  
  TH1D *kaonEta = new TH1D("kaonEta","Eta of thrown K;#eta",120,-6.,6.);
  TH1D *kaonRecEta = new TH1D("kaonRecEta","Eta of reco K;#eta",120,-6.,6.);
  
  TH1D *rho0Eta = new TH1D("rho0Eta","Eta of thrown #rho^{0};#eta",120,-6.,6.);
  TH1D *pipmfromrho0Eta = new TH1D("pipmfromrho0Eta","generated #eta of pi^{#pm} from rho^{0} decay;#eta",120,-6.,6.);
  TH1D *pipmfromrho0RecEta = new TH1D("pipmfromrho0RecEta","reconstructed #eta of pi^{#pm} from rho^{0} decay;#eta",120,-6.,6.);
  
  // phi(1020)
  TH1D *phiEta = new TH1D("phiEta","Eta of thrown #phi(1020);#eta",120,-6.,6.);
  TH1D *kpmfromphiEta = new TH1D("kpmfromphiEta","generated #eta of K^{#pm} from #phi(1020) decay;#eta",120,-6.,6.);
  TH1D *kpmfromphiRecEta = new TH1D("kpmfromphiRecEta","reconstructed #eta of K^{#pm} from #phi(1020) decay;#eta",120,-6.,6.);
  
  // momentum
  TH1D *partMom = new TH1D("partMom","Mom of thrown particles; P [GeV]",150,0.,150.);
  TH1D *recP = new TH1D("recP","Momentum of reconstructed tracks; P [GeV]",150,0.,150.);
  
  // phi
  TH1D *partPhi = new TH1D("partPhi","Phi of thrown charged particles; #phi [rad]",150,-3.2,3.2);
  TH1D *recPhi = new TH1D("recPhi","Phi of reconstructed tracks; #phi [rad]",150,-3.2,3.2);
  
  // count events:
  int ievgen = 0;
  // count generated particles:
  int ngen_electrons = 0; //+-
  int ngen_protons = 0; //+-
  int ngen_muons = 0; //+-
  int ngen_pions = 0; //+-
  int ngen_pi0=0;
  int ngen_kaons = 0; //+-
  int ngen_rho0 = 0;
  int ngen_rhop = 0;
  int ngen_phi = 0;
  int ngen_omega = 0;
  int ngen_jpsi = 0;
  int ngen_upsilon = 0;
  int ngen_d0 = 0;
  int ngen_b0 = 0;
  // count reconstructed particles (+-):
  int nrec_electrons = 0;
  int nrec_protons = 0;
  int nrec_muons = 0;
  int nrec_pions = 0;
  int nrec_kaons = 0;
  // count number of decays:
  int ndecay_pi0_gg = 0;
  int ndecay_rho0_pp = 0;
  int ndecay_rho0_mumu = 0;
  int ndecay_rho0_ee = 0;
  int ndecay_phi_kk = 0;
  // count number of decay particles (reco level) in nHCal acceptance:
  int ndecay_rho0_pionpm_nHCal = 0;
  int ndecay_phi_kaonpm_nHCal = 0;
  int ndecay_phi_one_nHCal = 0;
  int ndecay_phi_zero_nHCal = 0;
  int ndecay_phi_two_nHCal = 0;
  // count number of decay particles (reco level) in barrel acceptance:
  int ndecay_rho0_pionpm_barrel = 0;
  int ndecay_phi_kaonpm_barrel = 0;
  int ndecay_phi_one_barrel = 0;
  int ndecay_phi_zero_barrel = 0;
  int ndecay_phi_two_barrel = 0;
  // count number of decay particles (reco level) in LFHCal acceptance:
  int ndecay_rho0_pionpm_LFHCal = 0;
  int ndecay_phi_kaonpm_LFHCal = 0;
  int ndecay_phi_zero_LFHCal = 0;
  int ndecay_phi_one_LFHCal = 0;
  int ndecay_phi_two_LFHCal = 0;
  // count number of decay particles (reco level) in either hadronic calorimeter acceptance:
  int ndecay_rho0_pionpm_both = 0;
  int ndecay_phi_kaonpm_both = 0;
  int ndecay_phi_zero_both = 0;
  int ndecay_phi_one_both = 0;
  int ndecay_phi_two_both = 0;
  // count number of decay particles (reco level) in any hadronic calorimeter acceptance:
  int ndecay_rho0_pionpm_any = 0;
  int ndecay_phi_kaonpm_any = 0;
  int ndecay_phi_zero_any = 0;
  int ndecay_phi_one_any = 0;
  int ndecay_phi_two_any = 0;
  //array counter for specific k1 and k2
  int k2_k1[3][3] = {
	{0,0,0},
	{0,0,0},
	{0,0,0}
  };
  
  // tag decays on generated particle level:
  int is_rho0decay_pp = 0; // 0 or 1 for a given generated particle 
  int is_rho0decay_mumu = 0;
  int is_rho0decay_ee = 0;
  int is_phidecay_kk = 0;

  cout << "+ Ready to run over events... \n"; 
  
  
  while(tree_reader.Next()) { // Loop over events

    ievgen++;
	
	//uncomment to print every event loop
    cout << "+ Entering event #: " << ievgen << " \n";    
    
    //cout << "Event #: " << ievgen << ", " << partGenStat.GetSize() << " gen particles, " << parents_index.GetSize() << " parent particles, " << daughters_index.GetSize() << " daughter particles \n";   // parent_index and daughter_index must be of the same length since they are in the same tree (is that what pushback does?)


//Fill Event-Level Plots


xBjorkenPlot->Fill(xBjorken[0]);
QSquaredPlot->Fill(QSquared[0]);
QSquared_VS_xBjorken->Fill(xBjorken[0], QSquared[0]);

//reset kaon counters for this event:
int k_count_nHCal = 0; 
int k_count_barrel = 0;
int k_count_LFHCal = 0;
int k_count_both = 0;
int k_count_any = 0;
int k_select_p = 0;
int k_select_m = 0;



    for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over generated particles
      {
	//cout << "++ Entering generated particle #: " << i << " \n";
	
	int pdg = TMath::Abs(partPdg[i]);
	TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);
	float trueEta = trueMom.PseudoRapidity();
	float truePhi = trueMom.Phi();

	generatorStatus->Fill(partGenStat[i]);

	// reset particle decays for each new generated particle:
	is_rho0decay_pp = 0;
	is_rho0decay_mumu = 0;
	is_rho0decay_ee = 0;
	is_phidecay_kk = 0;
	
	int i_parents_begin = parents_begin[i];
        int i_parents_end = parents_end[i];
	int i_parents = parents_end[i] - parents_begin[i];

	int i_daughters_begin = daughters_begin[i]; 
        int i_daughters_end = daughters_end[i] - 1;	
	int i_daughters = daughters_end[i] - daughters_begin[i];

	//Consider only selected generated particles:
	if( (partGenStat[i] == 1) || (partGenStat[i] == 2) ) // Select only stable or decay particles
	{
	  //cout << "Ev#: " << ievgen << ", P-index: " << i <<", PDG: " << partPdg[i] << ", GenStatus:" << partGenStat[i] << ", i_parents: "<< i_parents<<", i_daughters: " << i_daughters << ", pb: " << parents_index[i_parents_begin] << ", pe: " << parents_index[i_parents_end] <<  ", db: " << daughters_index[i_daughters_begin] << ", de: " << daughters_index[i_daughters_end] << " \n";

	  // rho0 decays
	  if( partPdg[i] == 113 )
	    {
	      //cout << "Event " << ievgen << " with gen rho0 #: " << partPdg[i] << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";

	      // count the 2-body rho0 decays:
	      if( i_daughters == 2 )
		{
		  if( (partPdg[daughters_index[i_daughters_begin]] == 211 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -211 ) || (partPdg[daughters_index[i_daughters_begin]] == -211 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 211) )
		    {
		      //cout << "-> Event " << ievgen << " found rho0 decayed into pi+ pi-: "  << partPdg[daughters_index[i_daughters_begin]] << " and " << partPdg[daughters_index[i_daughters_begin]+1] << "  \n";
		      
		      // count it:
		      ndecay_rho0_pp++;
		      // tag it:
		      is_rho0decay_pp = 1;
		      
		      TVector3 trueMom_rho0_pi1(partMomX[daughters_index[i_daughters_begin]],partMomY[daughters_index[i_daughters_begin]],partMomZ[daughters_index[i_daughters_begin]]);
		      TVector3 trueMom_rho0_pi2(partMomX[daughters_index[i_daughters_begin]+1],partMomY[daughters_index[i_daughters_begin]+1],partMomZ[daughters_index[i_daughters_begin]+1]);
		      TVector3 trueMom_rho0_pi12=trueMom_rho0_pi1 + trueMom_rho0_pi2;
		      
		      float trueEta_rho0_pi1 = trueMom_rho0_pi1.PseudoRapidity();
		      float trueEta_rho0_pi2 = trueMom_rho0_pi2.PseudoRapidity();
		      pipmfromrho0Eta->Fill(trueEta_rho0_pi1);
		      pipmfromrho0Eta->Fill(trueEta_rho0_pi2);

		      //	      cout << "--> Event " << ievgen << " rho0 decay to 2pi: generated rho0 eta: " << trueEta << ", pi1: " << trueEta_rho0_pi1 << ", pi2: " << trueEta_rho0_pi2 << "  \n";
		      //  cout << "            trueMomrho0 X: " << trueMom.X() << ", trueMomrho0 Y: " << trueMom.Y() <<", trueMomrho0 Z: " << trueMom.Z() << "  \n";
		      // cout << "            trueMom_rho0_pi12 X: " << trueMom_rho0_pi12.X() << ", trueMom_rho0_pi12 Y: " << trueMompi12.Y() <<", trueMom_rho0_pi12 Z: " << trueMom_rho0_pi12.Z() << "  \n";
		    
		    } // end of rho0 to pi+pi- decays
		  else if( (partPdg[daughters_index[i_daughters_begin]] == 13 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -13 ) || (partPdg[daughters_index[i_daughters_begin]] == -13 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 13) )
		    {
		      ndecay_rho0_mumu++;
		    } // end of rho0 to mumu decays
		  else if( (partPdg[daughters_index[i_daughters_begin]] == 11 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -11 ) || (partPdg[daughters_index[i_daughters_begin]] == -11 ) && ( partPdg[daughters_index[i_daughters_begin]] == 11) )
		    {
		      ndecay_rho0_ee++;
		    } // end of ee decays
		} // end of 2-body decays
	    } // end of rho0 decays
	  // pi0 decays
	  else if( partPdg[i] == 111 )
	    {
	      //cout << "Event with gen pi0 #: " << ievgen << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";
	      if( (partPdg[daughters_index[i_daughters_begin]] == 22) && (partPdg[daughters_index[i_daughters_begin]+1] == 22) )
		{
		  ndecay_pi0_gg++;
		} // end of gg decays
	    }// end of pi0 decay

	  // phi(1020) decays:
      	  if( partPdg[i] == 333 )
	    {
	      //cout << "Event " << ievgen << " with gen phi(1020): " << partPdg[i] << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";

	      // count the 2-body phi(1020) decays:
	      if( i_daughters == 2 )
		{
		  if( (partPdg[daughters_index[i_daughters_begin]] == 321 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -321 ) || (partPdg[daughters_index[i_daughters_begin]] == -321 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 321) )
		    {
		      //cout << "-> Event " << ievgen << " found phi(1020) decayed into K+ K-: "  << partPdg[daughters_index[i_daughters_begin]] << " and " << partPdg[daughters_index[i_daughters_begin]+1] << "  \n";
		      
		      // count it:
		      ndecay_phi_kk++;
		      // tag it:
		      is_phidecay_kk = 1;

		      TVector3 trueMom_phi_k1(partMomX[daughters_index[i_daughters_begin]],partMomY[daughters_index[i_daughters_begin]],partMomZ[daughters_index[i_daughters_begin]]);
		      TVector3 trueMom_phi_k2(partMomX[daughters_index[i_daughters_begin]+1],partMomY[daughters_index[i_daughters_begin]+1],partMomZ[daughters_index[i_daughters_begin]+1]);
		      TVector3 trueMom_phi_k12=trueMom_phi_k1 + trueMom_phi_k2;
		      
		      float trueEta_phi_k1 = trueMom_phi_k1.PseudoRapidity();
		      float trueEta_phi_k2 = trueMom_phi_k2.PseudoRapidity();
		      kpmfromphiEta->Fill(trueEta_phi_k1);
		      kpmfromphiEta->Fill(trueEta_phi_k2);

		      //cout << "--> Event " << ievgen << " phi(1020) decay to 2pi: generated phi eta: " << trueEta << ", pi1: " << trueEta_phi_k1 << ", pi2: " << trueEta_phi_k2 << "  \n";
		      //cout << "            trueMomphi X: " << trueMom.X() << ", trueMomphi Y: " << trueMom.Y() <<", trueMomphi Z: " << trueMom.Z() << "  \n";
		      //cout << "            trueMom_phi_k12 X: " << trueMom_phi_k12.X() << ", trueMom_phi_k12 Y: " << trueMom_phi_k12.Y() <<", trueMom_phi_k12 Z: " << trueMom_phi_k12.Z() << "  \n";
		    
		    } // end of phi to pi+pi- decays
		} // end of 2-body decays
	    } // end of phi meson decays
	  
	  // count any generated particles:
	  if( pdg == 11){
	    ngen_electrons++;
	    electronEta->Fill(trueEta);
	  }// electrons                                                                                                     
	  else if( pdg == 13){
	    ngen_muons++;
	    muonEta->Fill(trueEta);
	  }// muons                                                                                                         
	  else if( pdg == 211){
	    ngen_pions++;
	    pionEta->Fill(trueEta);
	  }//pions_pm
	  else if( pdg == 111){
	    ngen_pi0++;
	    pi0Eta->Fill(trueEta);
	  }//pions_pm  
	  else if( pdg == 321 ){
	    ngen_kaons++;
	    kaonEta->Fill(trueEta);
	  } // kaons_pm                                                                                                     
	  else if( pdg == 113){
	    ngen_rho0++;
	    rho0Eta->Fill(trueEta);
	  } // rho(770)    
	  else if( pdg == 443){
	    ngen_jpsi++;
	  } // J/Psi(1S)                                                                                                      
	  else if( pdg == 2212){
	    ngen_protons++;
	    protonEta->Fill(trueEta);
	  }// protons
	  else if( pdg == 213){
	    ngen_rhop++;
	  }// rhop
	  else if( pdg == 333){
	    ngen_phi++;
	    phiEta->Fill(trueEta);
	  }// phi(1020)
	  else if( pdg == 223){
	    ngen_omega++;
	  }// omega(982)
	  else if( pdg == 553){
	    ngen_upsilon++;
	  }// Upsilon(1S)
	  else if( pdg == 421){
	    ngen_d0++;
	  }// D0
	  else if( pdg == 511){
	    ngen_b0++;
	  }// B0
	    
	    
	  //Fill all true eta: 
	  partEta->Fill(trueEta);
	    
	  // Fill all true momentum:
	  partMom->Fill(trueMom.Mag());
	  
	  // Fill all true phi: 
	  partPhi->Fill(truePhi);
	  
	  
	  // Loop over associations to find matching ReconstructedChargedParticle
	  for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
	    {
	      //cout << "*** Event " << ievgen << ", generated particle " << i << ", simID " << j << " \n";    
	      
	      if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
		{
		  //cout << "***** Event " << ievgen << ", found association index: " << simuAssoc[j] << " \n";
		  
		  TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle                                                       
		  
		  float CPartEta = recMom.PseudoRapidity();
		  float CPartPhi = recMom.Phi();
		  
		  recEta->Fill(CPartEta);
		  recPhi->Fill(CPartPhi);
		  recP->Fill(recMom.Mag());
		  
		  //cout << "Particle is pdg: " << pdg << " .\n";
		  
		  if( pdg == 11){
		    nrec_electrons++;
		    electronRecEta->Fill(CPartEta);
		  }// electrons	       
		  else if( pdg == 13){
		    nrec_muons++;
		    muonRecEta->Fill(CPartEta);
		  }// muons			
		  else if( pdg == 211){
		    nrec_pions++;
		    pionRecEta->Fill(CPartEta);
		  }//pions
		  else if( pdg == 321){
		    nrec_kaons++;
		    kaonRecEta->Fill(CPartEta);
		  }//pions
		  else if( pdg == 2212){
		    nrec_protons++;
		    protonRecEta->Fill(CPartEta);
		  }// protons       
		} // end of matched association gen to rec
	      
	      // Match the decay particles to their recos:
	      // rho0 to pi+ pi-
	      if( is_rho0decay_pp )
		{
		  if( simuAssoc[j] == daughters_index[i_daughters_begin] ) // get the reco decay pi1 of the gen rho0, by accessing the correct MCParticle index
		    {
		      TVector3 recMom_rho0_pi1(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
		      float recEta_rho0_pi1 = recMom_rho0_pi1.PseudoRapidity();
		      float recPhi_rho0_pi1 = recMom_rho0_pi1.Phi();
		      pipmfromrho0RecEta->Fill(recEta_rho0_pi1);

		      // count the decay pions (reco level) that are within the nHCal acceptance, here pion1:
		      if( recEta_rho0_pi1 >= eta_min_nHCal && recEta_rho0_pi1 <= eta_max_nHCal )
			{
			  ndecay_rho0_pionpm_nHCal++;
			}
		      //cout << "---> Event " << ievgen << " rho0 decay, reco index rho0: " << j << " \n";
		      //cout << "          reco daughter-1 eta: " << recEta_rho0_pi1 << ", reco index daughter-1: " << daughters_index[i_daughters_begin] << " \n";
		    }// end of rho0 decay pi1
		  else if( simuAssoc[j] == daughters_index[i_daughters_begin]+1 ) // get the reco decay pi2 of the gen rho0
		    {
		      TVector3 recMom_rho0_pi2(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
		      float recEta_rho0_pi2 = recMom_rho0_pi2.PseudoRapidity();
		      float recPhi_rho0_pi2 = recMom_rho0_pi2.Phi();
		      pipmfromrho0RecEta->Fill(recEta_rho0_pi2);

		      // count the decay pions (reco level) that are within the nHCal acceptance, here pion2:
		      if( recEta_rho0_pi2 >= eta_min_nHCal && recEta_rho0_pi2 <= eta_max_nHCal )
			{
			  ndecay_rho0_pionpm_nHCal++;
			}
		      
		      //cout << "          reco daughter-2 eta: " << recEta_rho0_pi2  << ", reco index daughter-2: " << daughters_index[i_daughters_begin]+1 << " \n\n";
		    }// end of rho0 decay pi2
		} // end of rho0 decay into pp
	      // phi to K+ K-
	      else if( is_phidecay_kk )
		{
		  
		  if( simuAssoc[j] == daughters_index[i_daughters_begin] ) // get the reco decay pi1 of the gen rho0, by accessing the correct MCParticle index
		    {
		      TVector3 recMom_phi_k1(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
		      float recEta_phi_k1 = recMom_phi_k1.PseudoRapidity();
		      float recPhi_phi_k1 = recMom_phi_k1.Phi();
		      kpmfromphiRecEta->Fill(recEta_phi_k1);

		      // count the decay pions (reco level) that are within the nHCal acceptance, here pion1:
		      if( recEta_phi_k1 >= eta_min_nHCal && recEta_phi_k1 <= eta_max_nHCal )
			{
			  ndecay_phi_kaonpm_nHCal++;
			  k_count_nHCal++;
			  k_select_p = 1;
			}
			  if( recEta_phi_k1 >= eta_min_barrel && recEta_phi_k1 <= eta_max_barrel )
			{
			  ndecay_phi_kaonpm_barrel++;
			  k_count_barrel++;
			  k_select_p = 2;
			}
			  if( recEta_phi_k1 >= eta_min_LFHCal && recEta_phi_k1 <= eta_max_LFHCal )
			{
			  ndecay_phi_kaonpm_LFHCal++;
			  k_count_LFHCal++;
			  k_select_p = 3;
			}

			
		      //cout << "---> Event " << ievgen << " phi(1020) decay, reco index phi(1020): " << j << " \n";
		      //cout << "          reco daughter-1 eta: " << recEta_phi_k1 << ", reco index daughter-1: " << daughters_index[i_daughters_begin] << " \n";
		    }// end of phi(1020) decay K1
		  else if( simuAssoc[j] == daughters_index[i_daughters_begin]+1 ) // get the reco decay pi2 of the gen rho0
		    {
		      TVector3 recMom_phi_k2(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
		      float recEta_phi_k2 = recMom_phi_k2.PseudoRapidity();
		      float recPhi_phi_k2 = recMom_phi_k2.Phi();
		      kpmfromphiRecEta->Fill(recEta_phi_k2);

		      // count the decay pions (reco level) that are within the nHCal acceptance, here pion2:
		      if( recEta_phi_k2 >= eta_min_nHCal && recEta_phi_k2 <= eta_max_nHCal )
			{
			  ndecay_phi_kaonpm_nHCal++;
			  k_count_nHCal++;
			  k_select_m = 1;
			  
			}
			  if( recEta_phi_k2 >= eta_min_barrel && recEta_phi_k2 <= eta_max_barrel )
			{
			  ndecay_phi_kaonpm_barrel++;
			  k_count_barrel++;
			  k_select_m = 2;

			}  
			  if( recEta_phi_k2 >= eta_min_LFHCal && recEta_phi_k2 <= eta_max_LFHCal )
			{
			  ndecay_phi_kaonpm_LFHCal++;
			  k_count_LFHCal++;
			  k_select_m = 3;
			  
			}

		      //cout << "          reco daughter-2 eta: " << recEta_phi_k2  << ", reco index daughter-2: " << daughters_index[i_daughters_begin]+1 << " \n\n";
		    }// end of phi(1020) decay K2
			
		} // end of phi(1020) decay into KK
	      
	    }// End loop over associations    
	} // End stable or decay particles condition
      } // End loop over thrown particles, within that event
	  
	  
	  
	  //Indiscriminate kaon numbers per phi
	  k_count_both = k_count_nHCal + k_count_LFHCal;
	  k_count_any = k_count_nHCal + k_count_LFHCal + k_count_barrel;
	  
    if(k_count_nHCal == 0)
		   {
			   ndecay_phi_zero_nHCal++;
			   Acceptance_VS_xBjorken_nHCal_0->Fill(xBjorken[0]);
			   Acceptance_VS_QSquared_nHCal_0->Fill(QSquared[0]);
		   }
		   
		   //is only one kaon in range?
		   else if(k_count_nHCal == 1)
		   {
			   ndecay_phi_one_nHCal++;
			   Acceptance_VS_xBjorken_nHCal_1->Fill(xBjorken[0]);
			   Acceptance_VS_QSquared_nHCal_1->Fill(QSquared[0]);

		   }
		   
		   else if(k_count_nHCal == 2)
		   {
			   ndecay_phi_two_nHCal++;
			   Acceptance_VS_xBjorken_nHCal_2->Fill(xBjorken[0]);
			   Acceptance_VS_QSquared_nHCal_2->Fill(QSquared[0]);
		   }
	if(k_count_LFHCal == 0)
	       {
			   ndecay_phi_zero_LFHCal++;
		   }
		   else if(k_count_LFHCal == 1)
		   {
			   ndecay_phi_one_LFHCal++;
		   }
		   else if(k_count_LFHCal == 2)
		   {
			   ndecay_phi_two_LFHCal++;
		   }
	if(k_count_barrel == 0)
	       {
			   ndecay_phi_zero_barrel++;
		   }
		   else if(k_count_barrel == 1)
		   {
			   ndecay_phi_one_barrel++;
		   }
		   else if(k_count_barrel == 2)
		   {
			   ndecay_phi_two_barrel++;
		   }
	if(k_count_both == 0)
	       {
			   ndecay_phi_zero_both++;
		   }
		   else if(k_count_both == 1)
		   {
			   ndecay_phi_one_both++;
		   }
		   else if(k_count_both == 2)
		   {
			   ndecay_phi_two_both++;
		   }
    if(k_count_any == 0)
	       {
			   ndecay_phi_zero_any++;
		   }
		   else if(k_count_any == 1)
		   {
			   ndecay_phi_one_any++;
		   }
		   else if(k_count_any == 2)
		   {
			   ndecay_phi_two_any++;
		   }
		   
		   
		   
		   
		   //Kaons per phi sensitive to p/m
		   if(k_select_p > 0 && k_select_m > 0)
		   {
			   k2_k1[k_select_m - 1][k_select_p - 1]++;
		   }
		
		   
		
    
    //for(unsigned int k=0; k<trackMomX.GetSize(); k++){ // Loop over all reconstructed tracks, thrown or not                                 

    //  TVector3 recMom(trackMomX[k], trackMomY[k], trackMomZ[k]);
	
    //} // End loop over all reconstructed tracks, within that event 

    // now go to next event
    
  } // End loop over events
  
  //combine the totals of combinations of detector parts:
  ndecay_phi_kaonpm_both = ndecay_phi_kaonpm_nHCal + ndecay_phi_kaonpm_LFHCal;
  ndecay_phi_kaonpm_any = ndecay_phi_kaonpm_nHCal + ndecay_phi_kaonpm_LFHCal + ndecay_phi_kaonpm_barrel;
  
  //Fill the Acceptance Plots for  nHCal
  for(int i=2; i<302; i++) {
	  
	  
	  float binVal = xBjorkenPlot->GetBinLowEdge(i);
	  
	  
	  float bin_data_0 = Acceptance_VS_xBjorken_nHCal_0->GetBinContent(i);
	  float bin_data_1 = Acceptance_VS_xBjorken_nHCal_1->GetBinContent(i);
	  float bin_data_2 = Acceptance_VS_xBjorken_nHCal_2->GetBinContent(i);
	  float bin_total = bin_data_0 + bin_data_1 + bin_data_2;
	  
	  
	  Acceptance_VS_xBjorken_nHCal->Fill(binVal, 100*(bin_data_0/bin_total));
	  
	  Acceptance_VS_xBjorken_nHCal->Fill(binVal, 100*(bin_data_1/bin_total));
	  Acceptance_VS_xBjorken_nHCal->Fill(binVal, 100*(bin_data_1/bin_total));
	  
	  
	  Acceptance_VS_xBjorken_nHCal->Fill(binVal, 100*(bin_data_2/bin_total));
	  Acceptance_VS_xBjorken_nHCal->Fill(binVal, 100*(bin_data_2/bin_total));
	  Acceptance_VS_xBjorken_nHCal->Fill(binVal, 100*(bin_data_2/bin_total));
	  
	  
	  //By repeating the fill for each data point, the different frequencies will automatically create a color gradient
	  //This is a very temporary solution before I transfer this to the plotting macro
  }
  
  
  
  

  // Calculate fractions:
  double fraction_rho0_pionpm_nHCal = 0.;
  fraction_rho0_pionpm_nHCal = ndecay_rho0_pp?(double(ndecay_rho0_pionpm_nHCal)/(2*double(ndecay_rho0_pp))):0;
  double fraction_phi_kaonpm_nHCal = 0.;
  fraction_phi_kaonpm_nHCal = ndecay_phi_kk?(double(ndecay_phi_kaonpm_nHCal)/(2*double(ndecay_phi_kk))):0;
  double fraction_rho0_pionpm_barrel = 0.;
  fraction_rho0_pionpm_barrel = ndecay_rho0_pp?(double(ndecay_rho0_pionpm_barrel)/(2*double(ndecay_rho0_pp))):0;
  double fraction_phi_kaonpm_barrel = 0.;
  fraction_phi_kaonpm_barrel = ndecay_phi_kk?(double(ndecay_phi_kaonpm_barrel)/(2*double(ndecay_phi_kk))):0;
  double fraction_rho0_pionpm_LFHCal = 0.;
  fraction_rho0_pionpm_LFHCal = ndecay_rho0_pp?(double(ndecay_rho0_pionpm_LFHCal)/(2*double(ndecay_rho0_pp))):0;
  double fraction_phi_kaonpm_LFHCal = 0.;
  fraction_phi_kaonpm_LFHCal = ndecay_phi_kk?(double(ndecay_phi_kaonpm_LFHCal)/(2*double(ndecay_phi_kk))):0;
  double fraction_rho0_pionpm_both = 0.;
  fraction_rho0_pionpm_both = ndecay_rho0_pp?(double(ndecay_rho0_pionpm_both)/(2*double(ndecay_rho0_pp))):0;
  double fraction_phi_kaonpm_both = 0.;
  fraction_phi_kaonpm_both = ndecay_phi_kk?(double(ndecay_phi_kaonpm_both)/(2*double(ndecay_phi_kk))):0;
  double fraction_rho0_pionpm_any = 0.;
  fraction_rho0_pionpm_any = ndecay_rho0_pp?(double(ndecay_rho0_pionpm_any)/(2*double(ndecay_rho0_pp))):0;
  double fraction_phi_kaonpm_any = 0.;
  fraction_phi_kaonpm_any = ndecay_phi_kk?(double(ndecay_phi_kaonpm_any)/(2*double(ndecay_phi_kk))):0;
  
  
  
  double fraction_phi_zero_nHCal = 0.;
  fraction_phi_zero_nHCal = ndecay_phi_kk?(double(ndecay_phi_zero_nHCal)/(double(ndecay_phi_kk))):0;
  double fraction_phi_one_nHCal = 0.;
  fraction_phi_one_nHCal = ndecay_phi_kk?(double(ndecay_phi_one_nHCal)/(double(ndecay_phi_kk))):0;
  double fraction_phi_two_nHCal = 0.;
  fraction_phi_two_nHCal = ndecay_phi_kk?(double(ndecay_phi_two_nHCal)/(double(ndecay_phi_kk))):0;
  
  
  double fraction_phi_zero_barrel = 0.;
  fraction_phi_zero_barrel = ndecay_phi_kk?(double(ndecay_phi_zero_barrel)/(double(ndecay_phi_kk))):0;
  double fraction_phi_one_barrel = 0.;
  fraction_phi_one_barrel = ndecay_phi_kk?(double(ndecay_phi_one_barrel)/(double(ndecay_phi_kk))):0;
  double fraction_phi_two_barrel = 0.;
  fraction_phi_two_barrel = ndecay_phi_kk?(double(ndecay_phi_two_barrel)/(double(ndecay_phi_kk))):0;
  
  
  double fraction_phi_zero_LFHCal = 0.;
  fraction_phi_zero_LFHCal = ndecay_phi_kk?(double(ndecay_phi_zero_LFHCal)/(double(ndecay_phi_kk))):0;
  double fraction_phi_one_LFHCal = 0.;
  fraction_phi_one_LFHCal = ndecay_phi_kk?(double(ndecay_phi_one_LFHCal)/(double(ndecay_phi_kk))):0;
  double fraction_phi_two_LFHCal = 0.;
  fraction_phi_two_LFHCal = ndecay_phi_kk?(double(ndecay_phi_two_LFHCal)/(double(ndecay_phi_kk))):0;
  
  
  double fraction_phi_zero_both = 0.;
  fraction_phi_zero_both = ndecay_phi_kk?(double(ndecay_phi_zero_both)/(double(ndecay_phi_kk))):0;
  double fraction_phi_one_both = 0.;
  fraction_phi_one_both = ndecay_phi_kk?(double(ndecay_phi_one_both)/(double(ndecay_phi_kk))):0;
  double fraction_phi_two_both = 0.;
  fraction_phi_two_both = ndecay_phi_kk?(double(ndecay_phi_two_both)/(double(ndecay_phi_kk))):0;
  
  
  double fraction_phi_zero_any = 0.;
  fraction_phi_zero_any = ndecay_phi_kk?(double(ndecay_phi_zero_any)/(double(ndecay_phi_kk))):0;
  double fraction_phi_one_any = 0.;
  fraction_phi_one_any = ndecay_phi_kk?(double(ndecay_phi_one_any)/(double(ndecay_phi_kk))):0;
  double fraction_phi_two_any = 0.;
  fraction_phi_two_any = ndecay_phi_kk?(double(ndecay_phi_two_any)/(double(ndecay_phi_kk))):0;
  
  double fraction_k2_k1[3][3] = {
	  {0,0,0},
	  {0,0,0},
	  {0,0,0}
      };
	
	for(int i=0;i<=2;i++)
	{
		for(int j=0;j<=2;j++)
		{
			fraction_k2_k1[i][j] = double(k2_k1[i][j]) / double(ndecay_phi_kk);
			
		}
	}
  
  
  
  
  //
  
  cout << "Number of generated events: " << ievgen << " \n\n";
  cout << "Number of generated electrons +-: " << ngen_electrons << " \n";
  cout << "Number of generated protons +-: " << ngen_protons << " \n";
  cout << "Number of generated muons +-: " << ngen_muons << " \n";
  cout << "Number of generated pions +-: " << ngen_pions << " \n";
  cout << "Number of generated pi0: " << ngen_pi0 << ", of which decay into 2 gamma: " << ndecay_pi0_gg <<  " \n";
  cout << "Number of generated kaons +-: " << ngen_kaons << " \n";
  cout << "Number of generated rho0: " << ngen_rho0 << ", of which decay into pi+ pi-: " << ndecay_rho0_pp << ", into mu+ mu-: " << ndecay_rho0_mumu << ", into e+ e-: " << ndecay_rho0_ee << " \n";
  cout << "        " << ndecay_rho0_pionpm_nHCal << " pi+ pi- make it into the nHCal acceptance, with corresponds to a fraction " << fraction_rho0_pionpm_nHCal << " \n";
  cout << "Number of generated rho+: " << ngen_rhop << " \n";
  cout << "Number of generated phi: " << ngen_phi <<", of which decay into K+ K-: " << ndecay_phi_kk << " \n";
  cout << "        " << ndecay_phi_kaonpm_nHCal << " K+ K- make it into the nHCal acceptance, which corresponds to a fraction " << fraction_phi_kaonpm_nHCal << " \n";
  cout << "        " << ndecay_phi_kaonpm_barrel << " K+ K- make it into the barrel acceptance, which corresponds to a fraction " << fraction_phi_kaonpm_barrel << " \n";
  cout << "        " << ndecay_phi_kaonpm_LFHCal << " K+ K- make it into the LFHCal acceptance, which corresponds to a fraction " << fraction_phi_kaonpm_LFHCal << " \n";
  cout << "        " << ndecay_phi_kaonpm_both << " K+ K- make it into either end's acceptance, which corresponds to a fraction " << fraction_phi_kaonpm_both << " \n";
  cout << "        " << ndecay_phi_kaonpm_any << " K+ K- make it into any part's acceptance, which corresponds to a fraction " << fraction_phi_kaonpm_any << " \n";
  cout << "Number of phi leaving zero kaons in nHCal: " << ndecay_phi_zero_nHCal << ", as fraction " << fraction_phi_zero_nHCal <<  "\n";
  cout << "Number of phi leaving one kaon in nHCal: " << ndecay_phi_one_nHCal << ", as fraction " << fraction_phi_one_nHCal << " \n";
  cout << "Number of phi leaving two kaons in nHCal: " << ndecay_phi_two_nHCal << ", as fraction " << fraction_phi_two_nHCal << " \n";
  cout << "Number of phi leaving zero kaons in the barrel: " << ndecay_phi_zero_barrel << ", as fraction " << fraction_phi_zero_barrel << " \n";
  cout << "Number of phi leaving one kaon in the barrel: " << ndecay_phi_one_barrel << ", as fraction " << fraction_phi_one_barrel << " \n";
  cout << "Number of phi leaving two kaons in the barrel: " << ndecay_phi_two_barrel << ", as fraction " << fraction_phi_two_barrel << " \n";
  cout << "Number of phi leaving zero kaons in LFHCal: " << ndecay_phi_zero_LFHCal << ", as fraction " << fraction_phi_zero_LFHCal << " \n";
  cout << "Number of phi leaving one kaon in LFHCal: " << ndecay_phi_one_LFHCal << ", as fraction " << fraction_phi_one_LFHCal << " \n";
  cout << "Number of phi leaving two kaons in LFHCal: " << ndecay_phi_two_LFHCal << ", as fraction " << fraction_phi_two_LFHCal << " \n";
  cout << "Number of phi leaving zero kaons in either end: " << ndecay_phi_zero_both << ", as fraction " << fraction_phi_zero_both << " \n";
  cout << "Number of phi leaving one kaon in either end: " << ndecay_phi_one_both << ", as fraction " << fraction_phi_one_both << " \n";
  cout << "Number of phi leaving two kaons in either end: " << ndecay_phi_two_both << ", as fraction " << fraction_phi_two_both << " \n";
  cout << "Number of phi leaving zero kaons in any detector range: " << ndecay_phi_zero_any << ", as fraction " << fraction_phi_zero_any << " \n";
  cout << "Number of phi leaving one kaon in any detector range: " << ndecay_phi_one_any << ", as fraction " << fraction_phi_one_any << " \n";
  cout << "Number of phi leaving two kaons in any detector range: " << ndecay_phi_two_any << ", as fraction " << fraction_phi_two_any << " \n";
  cout << "k2 vs k1 detected per phi (nhCal, Barrel, LFHCal): \n" << k2_k1[0][0] << ",  " << k2_k1[0][1] << ",  " << k2_k1[0][2] << "\n";
  cout << k2_k1[1][0] << ",  " << k2_k1[1][1] << ",  " << k2_k1[1][2] << "\n";
  cout << k2_k1[2][0] << ",  " << k2_k1[2][1] << ",  " << k2_k1[2][2] << "\n";
  cout << "Table as fraction of total phi mesons generated: \n";
  cout << fraction_k2_k1[0][0] << ",  " << fraction_k2_k1[0][1] << ",  " << fraction_k2_k1[0][2] << "\n";
  cout << fraction_k2_k1[1][0] << ",  " << fraction_k2_k1[1][1] << ",  " << fraction_k2_k1[1][2] << "\n";
  cout << fraction_k2_k1[2][0] << ",  " << fraction_k2_k1[2][1] << ",  " << fraction_k2_k1[2][2] << "\n";
  cout << "Number of generated omega: " << ngen_omega << " \n";
  cout << "Number of generated J/Psi: " << ngen_jpsi << " \n";
  cout << "Number of generated Upsilon: " << ngen_upsilon << " \n";
  cout << "Number of generated D0: " << ngen_d0 << " \n";
  cout << "Number of generated B0: " << ngen_b0 << " \n\n";
  cout << "Number of reconstructed electrons +-: " << nrec_electrons << " \n";
  cout << "Number of reconstructed protons +-: " << nrec_protons << " \n";
  cout << "Number of reconstructed muons +-: " << nrec_muons << " \n";
  cout << "Number of reconstructed pions +-: " << nrec_pions << " \n";
  cout << "Number of reconstructed kaons +-: " << nrec_kaons << " \n";
  
  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file

  cout << "Output histograms written in: " << outfile << " \n";
  
  cout << "Run again with root -l -q ePIC_Analysis_phi.C \n";
  
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
}
