#Generation of .h5 files to give as a input to the DNN
import ROOT
import pandas as pd
import os
import uproot
from concurrent.futures import ProcessPoolExecutor, as_completed
from func import *



#directory_sig = "/mnt/c/users/chiar/Desktop/ProgettoCERN/Dati_CERN/Segnale/"
#directory_bkg = "/mnt/c/users/chiar/Desktop/ProgettoCERN/Dati_CERN/Fondo/"
directory_sig = "/mnt/d/ProgettoCERN/Dati/Segnale"
directory_bkg = "/mnt/d/ProgettoCERN/Dati/Fondo"
columns_for_DNN = ["Dimuon_mass", "Sigma_Dimuon_mass", "Dimuon_pt", "Dimuon_y", "phi_CS", "cos_theta_CS", "Jet_pt", "Jet_eta", "Jet_phi", "jj_m", "Jet_DeltaEta", "min_j_dimuon_DeltaEta", "Zeppenfield", "pt_balance_ratio", "Jet_QGL", "label", "GenWeight"]



ROOT.gInterpreter.Declare("""
    #include <cmath>

    double EvaluateDeltaPhi(double phi1, double phi2) {
        double deltaPhi = std::fabs(phi1 - phi2);
        if (deltaPhi > TMath::Pi()) {
            deltaPhi = 2*TMath::Pi()-deltaPhi;
        }
        return deltaPhi;
    }
""")


# finds the first couple of muons with opposite charge
ROOT.gInterpreter.Declare("""
std::pair<int,int> FindOppositeChargePair(const ROOT::VecOps::RVec<float> &pt,
                                                const ROOT::VecOps::RVec<int> &charge) {
    // Lista di indici ordinati per pT decrescente
    std::vector<int> idx(pt.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int i, int j){ return pt[i] > pt[j]; });

    // Cerca la prima coppia a carica opposta nell'ordine di pT
    for (size_t i = 0; i < idx.size(); ++i) {
        for (size_t j = i+1; j < idx.size(); ++j) {
            if (charge[idx[i]] * charge[idx[j]] == -1) {
                return { idx[i], idx[j] }; // Leading, subleading
            }
        }
    }
    return { -1, -1 }; // Nessuna coppia trovata
}
""")


ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>
#include "TLorentzVector.h"

// Calcola tutte le coppie di TLorentzVector dei jet
std::vector<TLorentzVector> makeJetLV(const ROOT::VecOps::RVec<float> &pt,
                                      const ROOT::VecOps::RVec<float> &eta,
                                      const ROOT::VecOps::RVec<float> &phi,
                                      const ROOT::VecOps::RVec<float> &mass) {
    std::vector<TLorentzVector> jets;
    for (size_t i=0; i<pt.size(); ++i) {
        TLorentzVector v;
        v.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
        jets.push_back(v);
    }
    return jets;
}

// Calcola le variabili per tutte le coppie di jet
std::vector<float> Zeppenfield_per_coppia2(const std::vector<TLorentzVector> &jets,
                                          float Dimuon_y) {
    std::vector<float> Zeppen;
    for (size_t i=0; i<jets.size(); ++i) {
        for (size_t j=i+1; j<jets.size(); ++j) { // coppie senza ripetizioni
            float yL = jets[i].Rapidity();
            float yS = jets[j].Rapidity();
            float val = (Dimuon_y - (yL - yS)/2.0)/std::abs(yL - yS);
            Zeppen.push_back(val);
        }
    }
    return Zeppen;
}
                          
std::vector<float> Zeppenfield_per_coppia(const std::vector<TLorentzVector> &jets,
                                          float Dimuon_y) {
    std::vector<float> Zeppen;
    for (size_t i = 0; i < jets.size(); ++i) {
        for (size_t j = i + 1; j < jets.size(); ++j) { // coppie senza ripetizioni
            float yL = jets[i].Rapidity();
            float yS = jets[j].Rapidity();

            // Evita divisione per zero
            float denom = std::fabs(yL - yS);
            if (denom == 0.0f) {
                return {}; // array vuoto
            }

            float val = (Dimuon_y - (yL - yS) / 2.0f) / denom;

            // Controllo NaN, Inf, valori troppo grandi
            if (std::isnan(val) || std::isinf(val) || std::fabs(val) > 100000.0f) {
                return {}; // restituisce array vuoto
            }

            Zeppen.push_back(val);
        }
    }
    return Zeppen;
}

// Pt balance per coppia
std::vector<float> pt_balance_per_coppia(const std::vector<TLorentzVector> &jets,
                                         const TLorentzVector &Dimuon_LV) {
    std::vector<float> ratios;
    for (size_t i=0; i<jets.size(); ++i) {
        for (size_t j=i+1; j<jets.size(); ++j) {
            TLorentzVector jjLV = jets[i] + jets[j];
            TLorentzVector totalLV = jjLV + Dimuon_LV;
            float ratio = totalLV.Pt() / (Dimuon_LV.Pt() + jets[i].Pt() + jets[j].Pt());
            ratios.push_back(ratio);
        }
    }
    return ratios;
}

// DeltaEta minimo per coppia
std::vector<float> minDeltaEta_per_coppia(const std::vector<TLorentzVector> &jets,
                                          float Dimuon_eta) {
    std::vector<float> minDEta;
    for (size_t i=0; i<jets.size(); ++i) {
        for (size_t j=i+1; j<jets.size(); ++j) {
            float d1 = std::abs(Dimuon_eta - jets[i].Eta());
            float d2 = std::abs(Dimuon_eta - jets[j].Eta());
            minDEta.push_back(std::min(d1,d2));
        }
    }
    return minDEta;
}
                          
// Calcola DeltaEta e DeltaPhi per tutte le coppie di jet
std::vector<float> JetDeltaEta_per_coppia(const std::vector<TLorentzVector> &jets) {
    std::vector<float> deltaEta;
    for (size_t i=0; i<jets.size(); ++i) {
        for (size_t j=i+1; j<jets.size(); ++j) {
            deltaEta.push_back(std::abs(jets[i].Eta() - jets[j].Eta()));
        }
    }
    return deltaEta;
}

std::vector<float> JetDeltaPhi_per_coppia(const std::vector<TLorentzVector> &jets) {
    std::vector<float> deltaPhi;
    for (size_t i=0; i<jets.size(); ++i) {
        for (size_t j=i+1; j<jets.size(); ++j) {
            float dphi = std::fabs(jets[i].Phi() - jets[j].Phi());
            if (dphi > M_PI) dphi = 2*M_PI - dphi;
            deltaPhi.push_back(dphi);
        }
    }
    return deltaPhi;
}
                          
std::vector<float> jj_mass_per_coppia(const std::vector<TLorentzVector> &jets) {
    std::vector<float> masses;
    for (size_t i = 0; i < jets.size(); ++i) {
        for (size_t j = i+1; j < jets.size(); ++j) {
            TLorentzVector jj = jets[i] + jets[j];
            masses.push_back(jj.M());
        }
    }
    return masses;
}

std::vector<int> findAssociatedParticlesWithMotherIndexZero(
    const ROOT::VecOps::RVec<int>& pdgId,
    const ROOT::VecOps::RVec<int>& motherIdx,
    const ROOT::VecOps::RVec<float>& eta,
    const ROOT::VecOps::RVec<float>& pt
) {
    // 1. Cerca l'Higgs con madre all'indice 0.
    bool higgs_found = false;
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] == 25 && motherIdx[i] == 0) {
            higgs_found = true;
            break;
        }
    }

    if (!higgs_found) {
        return {}; // Se l'Higgs non ha madre all'indice 0, restituisci un array vuoto.
    }

    // 2. Cerca i quark candidati che hanno la stessa madre all'indice 0.
    std::vector<int> candidates;
    for (size_t i = 0; i < pdgId.size(); ++i) {
        // La condizione è che sia un quark (PdgId 1-6) e che abbia madre all'indice 0.
        if (std::abs(pdgId[i]) >= 1 && std::abs(pdgId[i]) <= 6 && motherIdx[i] == 0) {
            candidates.push_back(i);
        }
    }

    // Se non ci sono almeno due quark candidati, restituisci un array vuoto.
    if (candidates.size() < 2) {
        return {};
    }

    // 3. Trova la coppia di quark con il massimo DeltaEta.
    int max_idx1 = -1, max_idx2 = -1;
    float max_deta = -1.0;

    for (size_t i = 0; i < candidates.size(); ++i) {
        for (size_t j = i + 1; j < candidates.size(); ++j) {
            float deta = std::abs(eta[candidates[i]] - eta[candidates[j]]);
            if (deta > max_deta) {
                max_deta = deta;
                max_idx1 = candidates[i];
                max_idx2 = candidates[j];
            }
        }
    }

    // Restituisce gli indici dei due quark trovati.
    if (max_idx1 != -1 && max_idx2 != -1) {
        // Aggiunta la logica di ordinamento per pT
        if (pt[max_idx1] > pt[max_idx2]) {
            return {max_idx1, max_idx2};
        } else {
            return {max_idx2, max_idx1};
        }
    }

    return {};
}

#include <tuple>
std::tuple<std::vector<int>, double, double> match_vbf_quarks_to_jets(
    double q1_eta, double q1_phi,
    double q2_eta, double q2_phi,
    const ROOT::VecOps::RVec<float>& Jet_eta,
    const ROOT::VecOps::RVec<float>& Jet_phi) {

    std::vector<int> matched_jet_indices(2, -1);
    
    // Primo quark
    double min_dr1 = 999.0;
    double min_dr12 = 999.0;
    int best_jet_idx1 = -1;
    int best_jet_idx12 = -1;
    for (size_t i = 0; i < Jet_eta.size(); ++i) {
        double d_phi = EvaluateDeltaPhi(q1_phi, Jet_phi[i]);
        double d_eta = std::abs(q1_eta - Jet_eta[i]);
        double dr = std::sqrt(d_eta * d_eta + d_phi * d_phi);
        if (dr < min_dr1) {
            min_dr12 = min_dr1;
            best_jet_idx12 = best_jet_idx1;
            min_dr1 = dr;
            best_jet_idx1 = i;
        } else if (dr < min_dr12 && dr >= min_dr1){
            min_dr12 = dr;
            best_jet_idx12 = i;
        }
    }

    // Secondo quark
    double min_dr2 = 999.0;
    double min_dr22 = 999.0;
    int best_jet_idx2 = -1;
    int best_jet_idx22 = -1;
    for (size_t i = 0; i < Jet_eta.size(); ++i) {
        double d_phi = EvaluateDeltaPhi(q2_phi, Jet_phi[i]);
        double d_eta = std::abs(q2_eta - Jet_eta[i]);
        double dr = std::sqrt(d_eta * d_eta + d_phi * d_phi);
        if (dr < min_dr2) {
            min_dr22 = min_dr2;
            best_jet_idx22 = best_jet_idx2;
            min_dr2 = dr;
            best_jet_idx2 = i;
        } else if (dr < min_dr22 && dr >= min_dr2){
            min_dr22 = dr;
            best_jet_idx22 = i;
        }
    }      
    
    // controlla che i due quark non siano stati associati allo stesso jet e risolve le eventuali ambiguità
    if (best_jet_idx1 == best_jet_idx2){
        if (min_dr1 < min_dr2){
            best_jet_idx2 = best_jet_idx22;
            min_dr2 = min_dr22;
        } else {
            best_jet_idx1 = best_jet_idx12;
            min_dr1 = min_dr12;
        }
    }
                          
    matched_jet_indices[0] = best_jet_idx1;                     
    matched_jet_indices[1] = best_jet_idx2;
    
    return std::make_tuple(matched_jet_indices, min_dr1, min_dr2);
}

                      
ROOT::VecOps::RVec<int> makeJetLabels(const ROOT::VecOps::RVec<float> &Jet_pt, int idx1, int idx2) {
    ROOT::VecOps::RVec<int> labels(Jet_pt.size(), 0);
    if (idx1 >= 0 && idx2 >= 0 && idx1 < (int)Jet_pt.size() && idx2 < (int)Jet_pt.size()) {
        labels[idx1] = 1;
        labels[idx2] = 1;
    }
    return labels;
}

""")

ROOT.gInterpreter.Declare("""
float ComputeSigmaMassDimuon(
    float pt1,  float eta1,  float phi1,
    float s_pt1, float s_eta1, float s_phi1,
    float pt2,  float eta2,  float phi2,
    float s_pt2, float s_eta2, float s_phi2
) {
    // Differenze cinematiche
    float d_eta = eta1 - eta2;
    float d_phi = phi1 - phi2;

    // Invariante: m^2 = 2 * pt1 * pt2 * (cosh(deta) - cos(dphi))
    float f = std::cosh(d_eta) - std::cos(d_phi);
    float m = std::sqrt(2.0f * pt1 * pt2 * f);

    // Evita divisione per zero
    if (m <= 0.0f) {
        return 0.0f;
    }

    // Derivate parziali
    float dm_dpt1  = (pt2 * f) / m;
    float dm_dpt2  = (pt1 * f) / m;
    float dm_deta1 = (pt1 * pt2 * std::sinh(d_eta)) / m;
    float dm_deta2 = -(pt1 * pt2 * std::sinh(d_eta)) / m;
    float dm_dphi1 = (pt1 * pt2 * std::sin(d_phi)) / m;
    float dm_dphi2 = -(pt1 * pt2 * std::sin(d_phi)) / m;

    // Propagazione errori (senza termini di covarianza)
    float sigma_m = std::sqrt(
        std::pow(dm_dpt1  * s_pt1,  2) +
        std::pow(dm_dpt2  * s_pt2,  2) +
        std::pow(dm_deta1 * s_eta1, 2) +
        std::pow(dm_deta2 * s_eta2, 2) +
        std::pow(dm_dphi1 * s_phi1, 2) +
        std::pow(dm_dphi2 * s_phi2, 2)
    );

    return sigma_m;
}
""")


ROOT.gInterpreter.Declare("""
#include <cmath>
#include <algorithm>

const float MUON_MASS = 0.1056583745; // GeV

float wrap_dphi(float dphi) {
    // Wrap Δφ in [-π, π]
    const float pi = 3.14159265358979323846f;
    return std::fmod(dphi + pi, 2*pi) - pi;
}

float ComputeSigmaMassDimuon_ptErrOnly(
    float pt1,  float eta1,  float phi1,  float s_pt1,
    float pt2,  float eta2,  float phi2,  float s_pt2,
    bool massless
) {
    // Differenze cinematiche
    float d_eta = eta1 - eta2;
    float d_phi = wrap_dphi(phi1 - phi2);

    const float eps = 1e-12f;
    float mass = 0.0f;
    float denom = eps;
    float sigma_mass = 0.0f;

    if (massless) {
        // Massa zero
        float f = std::cosh(d_eta) - std::cos(d_phi);
        float m2 = 2.0f * pt1 * pt2 * f;
        mass = std::sqrt(std::max(m2, 0.0f));
        denom = (mass > eps) ? mass : eps;

        // Derivate rispetto a pt
        float dm_dpt1 = (pt2 * f) / denom;
        float dm_dpt2 = (pt1 * f) / denom;

        // Propagazione errore (solo pt)
        sigma_mass = std::sqrt(
            std::pow(dm_dpt1 * s_pt1, 2) +
            std::pow(dm_dpt2 * s_pt2, 2)
        );
    } else {
        // Massa reale muone
        float ET1 = std::sqrt(pt1*pt1 + MUON_MASS*MUON_MASS);
        float ET2 = std::sqrt(pt2*pt2 + MUON_MASS*MUON_MASS);
        float m2 = 2*ET1*ET2*std::cosh(d_eta) - 2*pt1*pt2*std::cos(d_phi) + 2*MUON_MASS*MUON_MASS;
        mass = std::sqrt(std::max(m2, 0.0f));
        denom = (mass > eps) ? mass : eps;

        float dET1_dpt1 = pt1 / ET1;
        float dET2_dpt2 = pt2 / ET2;

        float dm_dpt1 = (ET2 * dET1_dpt1 * std::cosh(d_eta) - pt2 * std::cos(d_phi)) / denom;
        float dm_dpt2 = (ET1 * dET2_dpt2 * std::cosh(d_eta) - pt1 * std::cos(d_phi)) / denom;

        sigma_mass = std::sqrt(
            std::pow(dm_dpt1 * s_pt1, 2) +
            std::pow(dm_dpt2 * s_pt2, 2)
        );
    }

    return sigma_mass;
}
""")


ROOT.gInterpreter.Declare("""
#include <cmath>
#include <array>

std::array<double,4> fourvec_from_ptetaphi(double pt, double eta, double phi, double mass = 0.1056583745)
{
    double px = pt * std::cos(phi);
    double py = pt * std::sin(phi);
    double pz = pt * std::sinh(eta);
    double p2 = px*px + py*py + pz*pz;
    double E  = std::sqrt(p2 + mass*mass);
    return {E, px, py, pz};
}
""")


ROOT.gInterpreter.Declare("""
#include <array>

std::array<double,4> boost_fourvec(const std::array<double,4>& p, const std::array<double,3>& beta)
{
    double bx = beta[0], by = beta[1], bz = beta[2];
    double b2 = bx*bx + by*by + bz*bz;
    if(b2 == 0.0) return p;

    double gamma = 1.0 / std::sqrt(1.0 - b2);
    double bp = bx*p[1] + by*p[2] + bz*p[3];
    double gamma2 = (gamma - 1.0)/b2;

    double E_prime = gamma*(p[0] - bp);
    double px_prime = p[1] + (-gamma*bx*p[0] + gamma2*bp*bx);
    double py_prime = p[2] + (-gamma*by*p[0] + gamma2*bp*by);
    double pz_prime = p[3] + (-gamma*bz*p[0] + gamma2*bp*bz);

    return {E_prime, px_prime, py_prime, pz_prime};
}
""")


ROOT.gInterpreter.Declare("""
#include <array>
#include <cmath>

std::array<double,3> unit(const std::array<double,3>& v, double eps=1e-12)
{
    double n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if(n < eps) return {0.0,0.0,0.0};
    return {v[0]/n, v[1]/n, v[2]/n};
}
""")


ROOT.gInterpreter.Declare("""
#include <cmath>
#include <limits>

double rapidity_from_ptetaphi(double pt, double eta, double phi, double mass = 0.1056583745)
{
    double px = pt * std::cos(phi);
    double py = pt * std::sin(phi);
    double pz = pt * std::sinh(eta);
    double E  = std::sqrt(px*px + py*py + pz*pz + mass*mass);

    double denom = E - pz;
    if(!(denom > 0.0)) return std::numeric_limits<double>::quiet_NaN();

    return 0.5 * std::log((E + pz)/denom);
}
""")


ROOT.gInterpreter.Declare("""
#include <array>
#include <cmath>

std::array<double,2> collins_soper_angles(
    double pt1, double eta1, double phi1,
    double pt2, double eta2, double phi2,
    double m_mu = 0.1056583745,
    bool use_mu_minus = true
)
{
    // quattro-vettori muoni
    auto p1 = fourvec_from_ptetaphi(pt1, eta1, phi1, m_mu);
    auto p2 = fourvec_from_ptetaphi(pt2, eta2, phi2, m_mu);

    // vettore totale dileptone
    std::array<double,4> q = {p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2], p1[3]+p2[3]};
    std::array<double,3> beta_to_rest = {q[1]/q[0], q[2]/q[0], q[3]/q[0]};

    // boost dei muoni nel frame dileptone
    auto p1_rf = boost_fourvec(p1, beta_to_rest);
    auto p2_rf = boost_fourvec(p2, beta_to_rest);

    // beam along +z e -z
    std::array<double,4> pa_lab = {1.0,0.0,0.0,1.0};
    std::array<double,4> pb_lab = {1.0,0.0,0.0,-1.0};

    auto pa_rf = boost_fourvec(pa_lab, beta_to_rest);
    auto pb_rf = boost_fourvec(pb_lab, beta_to_rest);

    auto ka = unit({pa_rf[1], pa_rf[2], pa_rf[3]});
    auto kb = unit({pb_rf[1], pb_rf[2], pb_rf[3]});

    std::array<double,3> z_cs = {ka[0]-kb[0], ka[1]-kb[1], ka[2]-kb[2]};
    z_cs = unit(z_cs);

    std::array<double,3> x_tmp = {ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2]};
    double dot = x_tmp[0]*z_cs[0] + x_tmp[1]*z_cs[1] + x_tmp[2]*z_cs[2];
    std::array<double,3> x_cs = {x_tmp[0]-dot*z_cs[0], x_tmp[1]-dot*z_cs[1], x_tmp[2]-dot*z_cs[2]};
    x_cs = unit(x_cs);

    std::array<double,3> y_cs = {z_cs[1]*x_cs[2] - z_cs[2]*x_cs[1],
                                 z_cs[2]*x_cs[0] - z_cs[0]*x_cs[2],
                                 z_cs[0]*x_cs[1] - z_cs[1]*x_cs[0]};

    auto l_rf = use_mu_minus ? std::array<double,3>{p1_rf[1],p1_rf[2],p1_rf[3]} :
                               std::array<double,3>{p2_rf[1],p2_rf[2],p2_rf[3]};
    double l_norm = std::sqrt(l_rf[0]*l_rf[0] + l_rf[1]*l_rf[1] + l_rf[2]*l_rf[2]);
    std::array<double,3> l_rf_unit = {l_rf[0]/l_norm, l_rf[1]/l_norm, l_rf[2]/l_norm};

    double cos_theta_cs = l_rf_unit[0]*z_cs[0] + l_rf_unit[1]*z_cs[1] + l_rf_unit[2]*z_cs[2];
    double phi_cs = std::atan2(l_rf_unit[0]*y_cs[0] + l_rf_unit[1]*y_cs[1] + l_rf_unit[2]*y_cs[2],
                               l_rf_unit[0]*x_cs[0] + l_rf_unit[1]*x_cs[1] + l_rf_unit[2]*x_cs[2]);

    return {cos_theta_cs, phi_cs};
}
""")

ROOT.gInterpreter.Declare("""

ROOT::VecOps::RVec<int> returnBackground(const ROOT::VecOps::RVec<float> &pt) {
    ROOT::VecOps::RVec<int> result(pt.size(), 0);
    return result;
}

ROOT::VecOps::RVec<int> returnErrorLabel(const ROOT::VecOps::RVec<float> &pt) {
    ROOT::VecOps::RVec<int> result(pt.size(), -1);
    return result;
}
""")



def process_file(input_file, outdir, columns_for_DNN, prefix="Signal_VBF"):
    print(f"Processing {input_file}")
    # creates a ROOT dataframe 
    dataframe = ROOT.RDataFrame("Events", input_file)

    # applies filters to the signal dataset
    df_DNN = (
        dataframe.Filter("HLT_IsoMu24 == true", "HLT_IsoMu24 cut")
        .Filter("All(Muon_tightId == true)", "tightId cut")
        .Filter("All(Muon_pfIsoId >=4)", "pfIsoId tight cut")
        .Filter("All(Muon_eta < 2.4) && All (Muon_eta > -2.4)", "Eta cut 2.4")
        .Filter("All(Jet_eta < 5) && All (Jet_eta > -5)", "Eta cut 5")
        #.Filter("All(Jet_btagPNetB < 0.6)", "btag cut")
        .Filter("nMuon >= 2", "two or more muons")
        .Filter("nJet >= 2", "two or more jets")
        .Filter("Any(Muon_charge == 1) && Any(Muon_charge == -1)", "2 muons with opposite charge")
        .Filter("Muon_pt[0] > 25", "Muon1_pt > 25")

        .Define("opp_pair", "FindOppositeChargePair(Muon_pt, Muon_charge)")
        # Indici separati
        .Define("idx_mu1", "opp_pair.first")
        .Define("idx_mu2", "opp_pair.second")
        # Filtra solo eventi dove coppia esiste
        .Filter("idx_mu1 >= 0 && idx_mu2 >= 0", "Opposite charge pair found")

        # adds to the dataframe columns containing the leading muon variables
        .Define("Muon1_pt", "Muon_pt[idx_mu1]")
        .Define("Muon1_eta", "Muon_eta[idx_mu1]")
        .Define("Muon1_phi", "Muon_phi[idx_mu1]")
        .Define("Muon1_mass", "Muon_mass[idx_mu1]")
        # adds to the dataframe columns containing the subleading muon variables
        .Define("Muon2_pt", "Muon_pt[idx_mu2]")
        .Define("Muon2_eta", "Muon_eta[idx_mu2]")
        .Define("Muon2_phi", "Muon_phi[idx_mu2]")
        .Define("Muon2_mass", "Muon_mass[idx_mu2]")
        # adds to the dataframe the differences in pt, eta, phi and R
        .Define("Muon_DeltaPt", "Muon1_pt - Muon2_pt")
        .Define("Muon_DeltaPhi", "EvaluateDeltaPhi(Muon1_phi, Muon2_phi)")
        .Define("Muon_DeltaEta", "abs(Muon1_eta - Muon2_eta)")
        .Define("Muon_DeltaR", "sqrt(pow(Muon_DeltaEta, 2) + pow(Muon_DeltaPhi, 2))")
        # adds to the dataframe the variables related to the dimuon system
        .Define("Dimuon_LV", "sum_TLorentzVector(Muon1_pt, Muon1_eta, Muon1_phi, Muon1_mass, Muon2_pt, Muon2_eta, Muon2_phi, Muon2_mass)")
        .Define("Dimuon_pt", "Dimuon_LV.Pt()")
        .Define("Dimuon_eta", "Dimuon_LV.Eta()")
        .Define("Dimuon_phi", "Dimuon_LV.Phi()")
        .Define("Dimuon_mass", "Dimuon_LV.M()")

        # adds to the dataframe columns containing the variables of the leading and subleading jets
        .Define("JetL_pt", "Jet_pt[0]")
        .Define("JetL_eta", "Jet_eta[0]")
        .Define("JetL_phi", "Jet_phi[0]")
        .Define("JetL_mass", "Jet_mass[0]")
        .Define("JetS_pt", "Jet_pt[1]")
        .Define("JetS_eta", "Jet_eta[1]")
        .Define("JetS_phi", "Jet_phi[1]")
        .Define("JetS_mass", "Jet_mass[1]")
        # adds to the dataframe the differences in pt, eta, phi and R
        .Define("JetLS_DeltaPt", "JetL_pt - JetS_pt")
        .Define("JetLS_DeltaEta", "abs(JetL_eta - JetS_eta)")
        .Define("JetLS_DeltaPhi", "EvaluateDeltaPhi(JetL_phi, JetS_phi)")
        .Define("JetLS_DeltaR", "sqrt(pow(JetLS_DeltaEta,2) + pow(JetLS_DeltaPhi,2))")
        # adds to the dataframe the variables related to the dijet system
        .Define("jjLS_LV", "sum_TLorentzVector(JetL_pt, JetL_eta, JetL_phi, JetL_mass, JetS_pt, JetS_eta, JetS_phi, JetS_mass)")
        .Define("jjLS_pt", "jjLS_LV.Pt()")
        .Define("jjLS_eta", "jjLS_LV.Eta()")
        .Define("jjLS_phi", "jjLS_LV.Phi()")
        .Define("jjLS_m", "jjLS_LV.M()")

        # adds to the dataframe the uncertainties on pt, eta and phi for the two muons
        .Define("Muon1_ptErr", "Muon_ptErr[idx_mu1]")
        .Define("Muon2_ptErr", "Muon_ptErr[idx_mu2]")
        # adds to the dataframe the uncertainty on the dimuon mass
        .Define("Sigma_Dimuon_mass", "ComputeSigmaMassDimuon_ptErrOnly(Muon1_pt, Muon1_eta, Muon1_phi, Muon1_ptErr, Muon2_pt, Muon2_eta, Muon2_phi, Muon2_ptErr, false)")
        # Filtra solo eventi dove l'incertezza è > 0
        .Filter("Sigma_Dimuon_mass > 0", "sigma dimuon mass found")
        # adds to the dataframe the rapidity of the dimuon system
        .Define("Dimuon_y", "Dimuon_LV.Rapidity()")
        #calcola cos(theta_CS) e phi_CS usando i pt, eta, phi dei due muoni
        .Define("CS_angles", "collins_soper_angles(Muon1_pt, Muon1_eta, Muon1_phi, Muon2_pt, Muon2_eta, Muon2_phi)")
        #separa le due componenti in colonne distinte
        .Define("cos_theta_CS", "CS_angles[0]")
        .Define("phi_CS", "CS_angles[1]")

        .Filter("JetL_pt> 35 && JetS_pt > 25", "jet_pt cut")
        .Filter("abs(JetL_eta) < 4.7 && abs(JetS_eta) < 4.7", "eta_LS < 4.7")
        .Filter("jjLS_m > 400")
        .Filter("JetLS_DeltaEta > 2.5")

        .Define("Jets_LV", "makeJetLV(Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
        .Define("Jet_DeltaEta", "JetDeltaEta_per_coppia(Jets_LV)")
        .Define("Jet_DeltaPhi", "JetDeltaPhi_per_coppia(Jets_LV)")
        .Define("jj_m", "jj_mass_per_coppia(Jets_LV)")
        
        
        .Define("Zeppenfield", "Zeppenfield_per_coppia(Jets_LV, Dimuon_y)")
        .Filter("Zeppenfield.size() > 0", "Scarta eventi con Zeppenfield vuoto")
        .Define("pt_balance_ratio", "pt_balance_per_coppia(Jets_LV, Dimuon_LV)")
        .Define("min_j_dimuon_DeltaEta", "minDeltaEta_per_coppia(Jets_LV, Dimuon_eta)")

        # computes the quark-glon likelihood for the jets
        .Define("Jet_QGL", "Jet_btagPNetQvG")
        .Filter("All(Jet_btagPNetQvG > -20 && Jet_btagPNetQvG < 20)", "QGL valid")

        # adds to the dataframe the weights of the Monte Carlo sampling
        .Define("GenWeight", "Generator_weight")
    )

    if prefix.startswith("S"):
        df_DNN = (
            df_DNN
            # adds the labels to the dataframe
            .Define("vbf_quark_indices", "findAssociatedParticlesWithMotherIndexZero(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_pt)")
            .Filter("vbf_quark_indices.size() == 2", "VBF Higgs event found")
            .Define("vbf_q1_idx", "vbf_quark_indices[0]")
            .Define("vbf_q2_idx", "vbf_quark_indices[1]")
            
            .Define("q1_pt", "GenPart_pt[vbf_q1_idx]")
            .Define("q1_eta", "GenPart_eta[vbf_q1_idx]")
            .Define("q1_phi", "GenPart_phi[vbf_q1_idx]")
            .Define("q1_mass", "GenPart_mass[vbf_q1_idx]")
            .Define("q2_pt", "GenPart_pt[vbf_q2_idx]")
            .Define("q2_eta", "GenPart_eta[vbf_q2_idx]")
            .Define("q2_phi", "GenPart_phi[vbf_q2_idx]")
            .Define("q2_mass", "GenPart_mass[vbf_q2_idx]")

            .Define("vbf_match_tuple", "match_vbf_quarks_to_jets(q1_eta, q1_phi, q2_eta, q2_phi, Jet_eta, Jet_phi)")
            
            .Define("VBFJet_indices", "std::get<0>(vbf_match_tuple)")
            .Define("qj1_DeltaR", "(double) std::get<1>(vbf_match_tuple)")
            .Define("qj2_DeltaR", "(double) std::get<2>(vbf_match_tuple)")
            
            .Define("VBFJet1_idx", "VBFJet_indices[0]")
            .Define("VBFJet2_idx", "VBFJet_indices[1]")

            .Filter("VBFJet1_idx >= 0 && VBFJet2_idx >=0", "Valid matched jets")
            .Filter("qj1_DeltaR <= 0.4 && qj2_DeltaR <= 0.4", "Delta_R_qj <= 0.4")

            .Define("label", "makeJetLabels(Jet_pt, VBFJet1_idx, VBFJet2_idx)")
        )
    elif prefix.startswith("B") or prefix.startswith("F"):
        df_DNN = (
            df_DNN.Define("label", "returnBackground(Jet_pt)")
        )
    else:
        df_DNN = (
            df_DNN.Define("label", "returnErrorLabel(Jet_pt)")
        )


    # Converte le colonne selezionate in un dizionario di array NumPy
    dict_data = df_DNN.AsNumpy(columns_for_DNN)
    # Converte in Pandas DataFrame
    pandas_df = pd.DataFrame(dict_data)

    filename_out = os.path.join(outdir, f"{prefix}_{os.path.basename(input_file).replace('.root', '.h5')}")
    pandas_df.to_hdf(filename_out, key="Events", mode="w")
    print(f"Saved {filename_out}")

    # calcola la somma dei genWeight su tutto il file
    file_opened = uproot.open(input_file)

    # Accedi alla TTree Runs
    runs = file_opened["Runs"]

    # Estrai e somma la variabile genEventSumw
    genEventSumw = runs["genEventSumw"].array(library="np").sum()
    return genEventSumw



outdir_sig = "/mnt/c/users/chiar/Desktop/ProgettoCERN/Dati_DNN_Jet_VBF/Segnale/"
outdir_bkg = "/mnt/c/users/chiar/Desktop/ProgettoCERN/Dati_DNN_Jet_VBF/Fondo/"

def process_directory(directory, outdir, columns_for_DNN, prefix="Signal"):
    files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".root")]
    genWeight_total = 0

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_file, f, outdir, columns_for_DNN, prefix) for f in files]
        for future in as_completed(futures):
            try:
                genWeight = future.result()  # Se process_file fallisce, qui viene sollevata l'eccezione
                genWeight_total += genWeight
            except Exception as e:
                print(f"Errore durante il processing di un file: {e}")
    print(f"Somma totale del pesi per la classe {prefix} = {genWeight_total}")

# Parallelizza i file di segnale
process_directory(directory_sig, outdir_sig, columns_for_DNN, prefix="Signal_VBF")

# Parallelizza i file di fondo
process_directory(directory_bkg, outdir_bkg, columns_for_DNN, prefix="Background_VBF")
