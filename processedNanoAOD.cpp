#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <string>
#include <string_view>
#include <vector>
#include <TLatex.h>
#include <map>
#include <cstdlib>

#include "samples.h"





bool ssWWTT = false;
bool ssWWLL = true;


// std::map<std::string, double> good_ele, good_muon, good_jet, good_fatjet;
// std::map<std::string, double> V_to_lep, V_to_had, VBS_presel;

std::map<TString, float> good_ele{};
std::map<std::string, float> good_muon, good_jet, good_fatjet;
std::map<TString, float> V_to_lep, V_to_had, VBS_presel;

good_ele["pt"]  = 15.0;
good_ele["eta"] = 4.7;
good_ele["iso"] = 0.06;

good_muon["pt"]  = 15.0;
good_muon["eta"] = 4.7;
good_muon["iso"] = 0.25;

good_jet["pt"]  = 30.0;
good_jet["eta"] = 4.7;
good_jet["wp_btag"] = 0.0479; // loose working point

good_fatjet["pt"]  = 30.0;
good_fatjet["eta"] = 4.7;
good_fatjet["msd"] = 25.0;

V_to_lep["ele_pt"]  = 30.0;
V_to_lep["muon_pt"] = 30.0;
V_to_lep["MET_pt"]  = 30.0;

V_to_had["m_sd_min"] = 60.0;
V_to_had["m_sd_max"] = 100.0;
V_to_had["tau21"]    = 0.45;

VBS_presel["m_jj"]       = 400.0;
VBS_presel["pt_leading"] = 50.0;
VBS_presel["deltaEta"]   = 2.5;


// --- Map bin -> description ---
std::map<int, std::string> cutDescriptions = {
    {1, "Exactly 1 e or 1 #mu"},
    {2, Form("good ele: p_{T} > %.0f , |#eta|< %.1f, isolation WP80 , iso < %.3f", good_ele["pt"], good_ele["eta"], good_ele["iso"])},
    {3, Form("good muon: p_{T} > %.0f , |#eta|< %.1f , tight ID , iso < %.3f", good_muon["pt"], good_muon["eta"], good_muon["iso"])},
    {4, Form("good AK4: p_{T}>%.0f , #Delta R(lep,jet)>0.4, |#eta|<%.1f", good_jet["pt"], good_jet["eta"])},
    {5, Form("good AK8: p_{T}>%.0f , m_{SD}>%.0f , #Delta R(l,FatJet)>0.8", good_fatjet["pt"], good_fatjet["msd"])},
    {6, ""},
    {7, ""},
    {8, Form("VBS kin : |#Delta #eta| > %.1f , m_{jj} > %.2f , p^{lead.}_{T} > %.0f", VBS_presel["deltaEta"], VBS_presel["m_jj"], VBS_presel["pt_leading"])},
    {9, Form("V to lep : p^e_{T} > %.0f or p^{#mu}_{T} > %.0f , MET > %.0f GeV", V_to_lep["ele_pt"], V_to_lep["muon_pt"], V_to_lep["MET_pt"])},
    {10, ""},
    {11, Form("V to had: %.0f < m_{sd} < %.0f , #tau_{21} < %.2f (if 1 AK8)", V_to_had["m_sd_min"], V_to_had["m_sd_max"], V_to_had["tau21"])},
    {12, Form("b-tagging: WP: %.4f", good_jet["wp_btag"])},
    {13, ""}
};


const double mW = 80.4; // W boson mass

double deltaPhi(double phi1, double phi2) {
    double dphi = phi1 - phi2;
    while (dphi > M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    return dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dEta = eta1 - eta2;
    double dPhi = deltaPhi(phi1, phi2);
    return std::sqrt(dEta*dEta + dPhi*dPhi);
}

double deltaEta(double eta1, double eta2) { return eta1 - eta2; }

double invariantMass(double pt1, double eta1, double phi1,
                     double pt2, double eta2, double phi2,
                     double m1=0., double m2=0.) 
{
    double theta1 = 2.0 * atan(exp(-eta1));
    double theta2 = 2.0 * atan(exp(-eta2));

    double px1 = pt1*cos(phi1), py1 = pt1*sin(phi1), pz1 = pt1/tan(theta1), E1 = sqrt(px1*px1 + py1*py1 + pz1*pz1 + m1*m1);
    double px2 = pt2*cos(phi2), py2 = pt2*sin(phi2), pz2 = pt2/tan(theta2), E2 = sqrt(px2*px2 + py2*py2 + pz2*pz2 + m2*m2);

    double px = px1 + px2, py = py1 + py2, pz = pz1 + pz2, E = E1 + E2;
    return sqrt(E*E - px*px - py*py - pz*pz);
}

void processNanoAOD() {

    TChain chain("Events");
    if (ssWWLL) {
        for (auto &f : samples_ssWWLL) chain.Add(f.c_str());
    } else {
        for (auto &f : samples_ssWWTT) chain.Add(f.c_str());
    }

    // Maximum sizes
    const int maxJets = 20;
    const int maxFatJets = 10;
    const int maxElectrons = 10;
    const int maxMuons = 10;

    // Branch arrays
    Float_t Jet_pt[maxJets], Jet_eta[maxJets], Jet_phi[maxJets], Jet_mass[maxJets], Jet_btagRobustParTAK4B[maxJets];
    Float_t FatJet_pt[maxFatJets], FatJet_eta[maxFatJets], FatJet_phi[maxFatJets], FatJet_mass[maxFatJets], FatJet_msoftdrop[maxFatJets], FatJet_tau1[maxFatJets], FatJet_tau2[maxFatJets];
    Float_t Electron_pt[maxElectrons], Electron_eta[maxElectrons], Electron_phi[maxElectrons], Electron_pfRelIso03_all[maxElectrons];
    Bool_t Electron_mvaIso_WP90[maxElectrons], Electron_mvaIso_WP80[maxElectrons];
    Float_t Muon_pt[maxMuons], Muon_eta[maxMuons], Muon_phi[maxMuons], Muon_pfRelIso04_all[maxMuons];
    Bool_t Muon_tightId[maxMuons];
    Float_t MET_pt, MET_phi;

    Int_t nJet, nFatJet, nElectron, nMuon;

    // Set branch addresses
    chain.SetBranchAddress("Jet_pt", Jet_pt); chain.SetBranchAddress("nJet", &nJet);
    chain.SetBranchAddress("Jet_eta", Jet_eta);
    chain.SetBranchAddress("Jet_phi", Jet_phi);
    chain.SetBranchAddress("Jet_mass", Jet_mass);

    chain.SetBranchAddress("FatJet_pt", FatJet_pt); chain.SetBranchAddress("nFatJet", &nFatJet);
    chain.SetBranchAddress("FatJet_eta", FatJet_eta);
    chain.SetBranchAddress("FatJet_phi", FatJet_phi);
    chain.SetBranchAddress("FatJet_mass", FatJet_mass);
    chain.SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop);
    chain.SetBranchAddress("FatJet_tau1", FatJet_tau1);
    chain.SetBranchAddress("FatJet_tau2", FatJet_tau2);

    chain.SetBranchAddress("Electron_pt", Electron_pt); 
    chain.SetBranchAddress("nElectron", &nElectron);
    chain.SetBranchAddress("Electron_pfRelIso03_all", &Electron_pfRelIso03_all);
    chain.SetBranchAddress("Electron_eta", Electron_eta);
    chain.SetBranchAddress("Electron_phi", Electron_phi);
    chain.SetBranchAddress("Electron_mvaIso_WP90", Electron_mvaIso_WP90);
    chain.SetBranchAddress("Electron_mvaIso_WP80", Electron_mvaIso_WP80);

    chain.SetBranchAddress("Muon_pt", Muon_pt); chain.SetBranchAddress("nMuon", &nMuon);
    chain.SetBranchAddress("Muon_eta", Muon_eta);
    chain.SetBranchAddress("Muon_phi", Muon_phi);
    chain.SetBranchAddress("Muon_tightId", Muon_tightId);
    chain.SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all);
    chain.SetBranchAddress("MET_pt", &MET_pt);
    chain.SetBranchAddress("MET_phi", &MET_phi);
    chain.SetBranchAddress("Jet_btagRobustParTAK4B", Jet_btagRobustParTAK4B);

    // === Histograms ===
    TH1F *h_nJet = new TH1F("h_nJet","Number of clean jets;N_{jet};Events",15,0,15);
    TH1F *h_nFatJet = new TH1F("h_nFatJet","Number of clean fatjets;N_{fatjet};Events",10,0,10);
    TH1F *h_max_mjj = new TH1F("h_max_mjj","Max m_{jj} of VBS jets; m_{jj} [GeV]; Events",50,0,2000);
    TH1F *h_deltaEta_vbs = new TH1F("h_deltaEta_vbs","#Delta#eta of VBS jets;#Delta#eta;Events",50,0,10);
    TH1F *h_MT_lepton = new TH1F("h_MT_lepton","Lepton-MET MT; MT [GeV]; Events",50,0,500);
    TH2F *h_nJet_vs_nFatJet = new TH2F("h_nJet_vs_nFatJet","NJet vs NFatJet;N_{jet};N_{fatjet}",15,0,15,10,0,10);

    TH1I *h_singleLepton = new TH1I("h_singleLepton", "h_singleLepton", 2,0,2);
    TH1I *h_singleCleanElectron_isowp90 = new TH1I("hsingleElectron_iswp90", "h_singleElectron_isowp90", 2, 0, 2);
    TH1I *h_singleCleanElectron_isowp80 = new TH1I("hsingleElectron_iswp80", "h_singleElectron_iswp80", 2, 0, 2);
    TH1I *h_singleCleanMuon = new TH1I("h_singleCleanMuon", "h_singleCleanMuon", 2, 0, 2);
    TH1I *h_singleCleanLepton = new TH1I("h_singleCleanLepton", "h_singleCleanLepton", 2, 0, 2);
    TH1I *h_CleanFatJet = new TH1I("h_CleanFatJet", "h_CleanFatJet", 2, 0, 2);
    TH1I *h_vbs_topology = new TH1I("h_vbs_topology", "h_vbs_topology", 2,0, 2);
    TH1I *h_resolvedWjet = new TH1I("h_resolvedWjet", "h_resolvedWjet", 2, 0, 2);
    TH1I *h_boostedWjet = new TH1I("h_boostedWjet", "h_boostedWjet", 2, 0, 2);

    // --- Cutflow histograms ---
    TH1F *h_cutflow = new TH1F("h_cutflow", "Cutflow;Selection step;Events", 13, 0.5, 10.5);
    h_cutflow->GetXaxis()->SetBinLabel(1, "Single lepton");
    h_cutflow->GetXaxis()->SetBinLabel(2, "Good Electron");
    h_cutflow->GetXaxis()->SetBinLabel(3, "Good Muon");
    h_cutflow->GetXaxis()->SetBinLabel(4, "Good Jet");
    h_cutflow->GetXaxis()->SetBinLabel(5, "Good FatJet");
    h_cutflow->GetXaxis()->SetBinLabel(6, ">= 4 AK4");
    h_cutflow->GetXaxis()->SetBinLabel(7, "1AK8, >= 2AK4");
    h_cutflow->GetXaxis()->SetBinLabel(8, "VBS kinematics");
    h_cutflow->GetXaxis()->SetBinLabel(9, "Lep.Decay");
    h_cutflow->GetXaxis()->SetBinLabel(10, "Had.Decay (2AK4)");
    h_cutflow->GetXaxis()->SetBinLabel(11, "Had.Decay (1AK8)");
    h_cutflow->GetXaxis()->SetBinLabel(12, "No b-tag ");
    h_cutflow->GetXaxis()->SetBinLabel(13, "selected events");

    Long64_t nentries = chain.GetEntries();
    std::cout << "Processing " << nentries << " events..." << std::endl;

    double Wmass = 80.4;
    std::cout<<nentries<<std::endl;
    for(Long64_t i=0; i<nentries; ++i){
        chain.GetEntry(i);
        if(i % 1000 == 0) std::cout << "Processing event nr. " << i << std::endl;

        bool singleLepton = true;
        if(nElectron != 1 && nMuon != 1){
            singleLepton = false;
            h_singleLepton->Fill(singleLepton);
            continue;
        }
        h_cutflow->AddBinContent(1); // Passed single lepton

        // good muons and electrons
        bool found_ele_wp80 = false, found_ele_wp90 = false, found_muon = false;
        if(nElectron == 1 && nMuon == 0){
            if(Electron_mvaIso_WP80[0] && abs(Electron_eta[0]) < good_ele["eta"] && Electron_pt[0] > good_ele["pt"] && Electron_pfRelIso03_all[0] < good_ele["iso"]){  h_cutflow->AddBinContent(2); found_ele_wp80 = true;}
            if(Electron_mvaIso_WP90[0] && abs(Electron_eta[0]) < good_ele["eta"] && Electron_pt[0] > good_ele["pt"] && Electron_pfRelIso03_all[0] < good_ele["iso"]) found_ele_wp90 = true; 
            h_singleCleanElectron_isowp80->Fill(found_ele_wp80);
            h_singleCleanElectron_isowp90->Fill(found_ele_wp90);
        } 
        if(nElectron == 0 && nMuon == 1){
            if(Muon_pt[0] > good_muon["pt"] && Muon_tightId[0] && Muon_pfRelIso04_all[0] < good_muon["iso"] && abs(Muon_eta[0]) < good_muon["eta"]) {h_cutflow->AddBinContent(3); found_muon = true; }
            h_singleCleanMuon->Fill(found_muon);
        }
        bool singleGoodLepton = true;
        if(!found_muon && !found_ele_wp80){
            singleGoodLepton = false;
            h_singleCleanLepton->Fill(singleGoodLepton);
            continue;
        }

        int cleanfatjet = 0;
        std::vector<int> cleanFatJet_idx;
        for(Int_t fj=0; fj<nFatJet; fj++){
            if(FatJet_pt[fj] < good_fatjet["pt"] || FatJet_msoftdrop[fj] < good_fatjet["msd"] || FatJet_eta[fj] > good_fatjet["eta"]){ h_CleanFatJet->Fill(false); continue; }
            bool overlaps = false;
            for(Int_t ee=0; ee<nElectron; ee++){
                if(deltaR(Electron_eta[ee],Electron_phi[ee],FatJet_eta[fj],FatJet_phi[fj]) < 0.8) overlaps = true;
            }
            for(Int_t mu=0; mu<nMuon; mu++){
                if(deltaR(Muon_eta[mu],Muon_phi[mu],FatJet_eta[fj],FatJet_phi[fj]) < 0.8) overlaps = true;
            }
            if(overlaps){ h_CleanFatJet->Fill(false); continue; }
            h_CleanFatJet->Fill(true);
            ++cleanfatjet;
            cleanFatJet_idx.emplace_back(fj);
        }

        // --- Clean jets
        int cleanjet = 0;
        std::vector<int> cleanjet_idx;
        for(Int_t j=0; j<nJet; j++){
            if(Jet_pt[j] < good_jet["pt"] || Jet_eta[j] > good_jet["eta"]) continue;
            bool overlaps=false;
            for(Int_t ee=0; ee<nElectron; ee++){
                if(deltaR(Electron_eta[ee],Electron_phi[ee],Jet_eta[j],Jet_phi[j]) < 0.4) overlaps=true;
            }
            for(Int_t mu=0; mu<nMuon; mu++){
                if(deltaR(Muon_eta[mu],Muon_phi[mu],Jet_eta[j],Jet_phi[j]) < 0.4) overlaps=true;
            }
            if(overlaps) continue;
            ++cleanjet;
            cleanjet_idx.push_back(j);
        }

        h_nJet->Fill(cleanjet);
        h_nFatJet->Fill(cleanfatjet);
        h_nJet_vs_nFatJet->Fill(cleanjet, cleanfatjet);
        if(cleanjet>0) h_cutflow->AddBinContent(4);
        if(cleanfatjet>0) h_cutflow->AddBinContent(5);

        bool resolved_Wjet = false;
        bool boosted_Wjet = false;
        if(cleanfatjet == 1 && cleanjet >= 2){
            boosted_Wjet=true;
            h_cutflow->AddBinContent(7);
        }
        if(cleanfatjet == 0 && cleanjet >= 4){
         h_cutflow->AddBinContent(6);
         resolved_Wjet=true;
        }
        if(!resolved_Wjet && !boosted_Wjet) continue;
        std::pair<double,double> vbs_jets;
        std::pair<int,int> vbs_jets_idx;
        double max_mjj = -1.0;
        for(Int_t j=0;j<nJet;j++){
            for(Int_t jj=j+1;jj<nJet;jj++){
                double m_jj = invariantMass(Jet_pt[j],Jet_eta[j],Jet_phi[j],
                                            Jet_pt[jj],Jet_eta[jj],Jet_phi[jj]);
                if(m_jj>max_mjj){
                    max_mjj=m_jj;
                    vbs_jets={deltaEta(Jet_eta[j],Jet_eta[jj]), m_jj};
                    vbs_jets_idx={j,jj};
                }
            }
        }
        h_max_mjj->Fill(max_mjj);
        h_deltaEta_vbs->Fill(std::abs(vbs_jets.first));

        bool vbs_topology = ((cleanjet>=2 && cleanfatjet==1) || (cleanjet>=4 && cleanfatjet==0));
        bool vbs_kinematics = (std::abs(vbs_jets.first)> VBS_presel["eta"] && vbs_jets.second > VBS_presel["m_jj"]);

        if(!vbs_kinematics) continue;
        h_cutflow->AddBinContent(8); // Passed VBS kinematics

        // --- Leptonic MT
        double MT=0;
        if(found_ele_wp80 || found_ele_wp90){
            MT = sqrt(2*Electron_pt[0]*MET_pt*(1-cos(Electron_phi[0]-MET_phi)));
            if(Electron_pt[0] > V_to_lep["ele_pt"] && MET_pt > V_to_lep["MET_pt"]) h_cutflow->AddBinContent(9);
        } else if(found_muon){
            MT = sqrt(2*Muon_pt[0]*MET_pt*(1-cos(Muon_phi[0]-MET_phi)));
            if(Muon_pt[0] > V_to_lep["muon_pt"] && MET_pt > V_to_lep["MET_pt"]) h_cutflow->AddBinContent(9);
        }
        if(MT > 0) h_MT_lepton->Fill(MT);

        // --- W reconstruction
 
        if(cleanfatjet == 0 && cleanjet>=4){
            std::pair<int,int> wjets_idx;
            double closest_mw_diff = 1e6;
            double w_mjj = -1.0;

            for (int j = 0; j < nJet; j++) {
                if (j == vbs_jets_idx.first || j == vbs_jets_idx.second) continue;
                for (int jj = j + 1; jj < nJet; jj++) {
                    if (jj == vbs_jets_idx.first || jj == vbs_jets_idx.second) continue;
                    double m_jj = invariantMass(Jet_pt[j], Jet_eta[j], Jet_phi[j],
                                                Jet_pt[jj], Jet_eta[jj], Jet_phi[jj]);
                    double diff = std::abs(m_jj - Wmass);
                    if (diff < closest_mw_diff) {
                        closest_mw_diff = diff;
                        w_mjj = m_jj;
                        wjets_idx = {j, jj};
                    }
                }
            }
            if(w_mjj > V_to_had["m_sd_min"] && w_mjj < V_to_had["m_sd_max"]){ h_cutflow->AddBinContent(10); resolved_Wjet = true;}
        } 
        else if(cleanfatjet==1 && cleanjet>=2) {
            if(FatJet_msoftdrop[cleanFatJet_idx[0]] > V_to_had["m_sd_min"] && FatJet_msoftdrop[cleanFatJet_idx[0]] < V_to_had["m_sd_max"] && (FatJet_tau2[cleanFatJet_idx[0]]/FatJet_tau1[cleanFatJet_idx[0]]) < V_to_had["tau21"]){
                boosted_Wjet = true;
                h_cutflow->AddBinContent(11);
            }
        }

        bool b_tagged = false;
        for(int jj= 0; jj<cleanjet_idx.size(); jj++){
            if(Jet_btagRobustParTAK4B[cleanjet_idx[jj]] > good_jet["wp_btag"] && abs(Jet_eta[cleanjet_idx[jj]]) < 2.5){
                b_tagged=true;
                break;
            }
        }
        if(!boosted_Wjet && !resolved_Wjet) continue;
        h_cutflow->AddBinContent(12);
        if(b_tagged) continue; 

        h_cutflow->AddBinContent(13);

        h_resolvedWjet->Fill(resolved_Wjet);
        h_boostedWjet->Fill(boosted_Wjet);
  
    }

    // --- Efficiency histogram ---
    TH1F *h_efficiency = (TH1F*)h_cutflow->Clone("h_efficiency");
    h_efficiency->SetTitle("Cut efficiencies");
    h_efficiency->SetYTitle("Efficiency");

    // --- Save output ---
    TFile* outputFile = new TFile("output.root","RECREATE");
    h_nJet->Write();
    h_nFatJet->Write();
    h_max_mjj->Write();
    h_deltaEta_vbs->Write();
    h_MT_lepton->Write();
    h_nJet_vs_nFatJet->Write();
    h_singleLepton->Write();
    h_singleCleanElectron_isowp90->Write();
    h_singleCleanElectron_isowp80->Write();
    h_singleCleanMuon->Write();
    h_singleCleanLepton->Write();
    h_CleanFatJet->Write();
    h_vbs_topology->Write();
    h_resolvedWjet->Write();
    h_boostedWjet->Write();
    h_cutflow->Write();
    h_efficiency->Write();
    outputFile->Close();
    TCanvas *c_cutflow = new TCanvas("c_cutflow", "Cutflow", 1200, 600);
h_cutflow->SetStats(0);
h_cutflow->SetFillColor(kAzure-9);
h_cutflow->SetBarWidth(0.9);
h_cutflow->SetBarOffset(0.05);
h_cutflow->Draw("hist bar");

std::map<int, std::string> eff = {
    {1,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(1) / nentries)},
    {2,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(2) / h_cutflow->GetBinContent(1))},
    {3,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(3) / h_cutflow->GetBinContent(1))},
    {4,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(4) / (h_cutflow->GetBinContent(2) + h_cutflow->GetBinContent(3)))},
    {5,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(5) / (h_cutflow->GetBinContent(2) + h_cutflow->GetBinContent(3)))},
    {6,  Form("#epsilon_{single lep.} = %.2f", h_cutflow->GetBinContent(6) / h_cutflow->GetBinContent(1))},
    {7,  Form("#epsilon_{single lep.} = %.2f", h_cutflow->GetBinContent(7) / h_cutflow->GetBinContent(1))},
    {8,  Form("#epsilon_{single lep.} = %.2f", h_cutflow->GetBinContent(8) / h_cutflow->GetBinContent(1))},
    {9,  Form("#epsilon_{single lep.} = %.2f", h_cutflow->GetBinContent(9) / h_cutflow->GetBinContent(1))},
    {10, Form("#epsilon_{single lep.} = %.2f", h_cutflow->GetBinContent(10) / h_cutflow->GetBinContent(1))},
    {11, Form("#epsilon_{single lep.} = %.2f", h_cutflow->GetBinContent(11) / h_cutflow->GetBinContent(1))},
    {12, Form("#epsilon_{single lep.} = %.2f", h_cutflow->GetBinContent(9) / (h_cutflow->GetBinContent(11) + h_cutflow->GetBinContent(10)))},
    {13, Form("#epsilon = %.2f", h_cutflow->GetBinContent(13) / nentries)}
};



// --- Add TLatex labels ---
TLatex latex;
latex.SetTextSize(0.03);
latex.SetTextAlign(12);  
latex.SetNDC(false);    
int count_text=0;
for (int bin = 1; bin <= 10; ++bin) {
    double x = h_cutflow->GetBinCenter(bin);
    double y = h_cutflow->GetBinContent(bin) * 1.6; 

 if(cutDescriptions[bin] != ""){ count_text++;
    latex.DrawLatexNDC(0.2, 0.9-0.05*count_text, cutDescriptions[bin].c_str());
 }  
    std::cout<<eff[bin]<<std::endl;
    latex.DrawLatexNDC(x - 0.3, 0.3, eff[bin].c_str());
 
}

c_cutflow->SetTicks(1, 1);
c_cutflow->Update();
c_cutflow->SaveAs("cuts.pdf");

    std::cout << "Histograms written to output.root" << std::endl;
}
