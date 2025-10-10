#####################
# funzioni utilizzate
#####################

import ROOT
import math


# funzione per il calcolo di DeltaPhi
ROOT.gInterpreter.Declare("""
    double DeltaPhiFunc(double phi1, double phi2) {
    double DeltaPhi = abs(phi1 - phi2);
        if (DeltaPhi > TMath::Pi()) {
            DeltaPhi = 2*TMath::Pi() - DeltaPhi;
        }
        return DeltaPhi;
    }
""")

def compute_DeltaPhi(phi1, phi2):
    DeltaPhi = abs(phi1 - phi2)
    if DeltaPhi > math.pi:
        DeltaPhi = 2*math.pi - DeltaPhi
    return DeltaPhi



# somma di due TLorentzVector dati pt, eta, phi, mass
ROOT.gInterpreter.Declare("""
    TLorentzVector sum_TLorentzVector(double pt1, double eta1, double phi1, double m1, 
                          double pt2, double eta2, double phi2, double m2) {
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
    v2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
    return v1 + v2;
    }
""")



# calcola DeltaR

ROOT.gInterpreter.Declare("""
    double compute_DeltaR(double eta1, double phi1, double eta2, double phi2) {
        return sqrt(pow(abs(eta1 - eta2), 2) + pow(DeltaPhiFunc(phi1, phi2), 2));
    }
""")

def compute_DeltaR(eta1, phi1, eta2, phi2):
    return math.sqrt(abs(eta1 - eta2)**2 + compute_DeltaPhi(phi1, phi2)**2)




# crea un nuovo branch a partire due branch facendo un append del secondo

ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    RVecD append_RVecD(const RVecD &v1, const RVecD &v2) {
        RVecD v(v1.size() + v2.size());
        for (int i = 0; i < v.size(); i++) {
            if (i < v1.size()) {
                v[i] = v1[i];
            } else {;
                v[i] = v2[i - v1.size()];
            }
            
        }
        return v;
    }

""")



# ritrona il DeltaR tra i due oggetti del generatore (eg. quark, muoni) e quelli ricostruiti (eg. jet)

ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    RVecD compute_DeltaR_gen_reco(double gen1_eta, double gen1_phi,
                            double gen2_eta, double gen2_phi,
                          const RVecD &reco_eta, const RVecD &reco_phi) {
        RVecD dR(2 * reco_eta.size());
        RVecD gen = {gen1_eta, gen2_eta, gen1_phi, gen2_phi};
        double temp = 0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < reco_eta.size(); j++) {
                dR[temp] = compute_DeltaR(gen[i], gen[i+2], reco_eta[j], reco_phi[j]);
                temp++;
            }
        }
        return dR;
    }
""")



# calcola la massa invariante
#ROOT.gInterpreter.Declare("""
    #double invariant_mass_func(double pt1, double pt2, double delta_eta, double delta_phi) {
        #return sqrt(2*pt1*pt2*(cosh(delta_eta) - cos(delta_phi)));
    #}
#""")




### DA MIGLIORARE ###
# calcola gli indici della prima coppia di muoni con carica opposta, tightId e pfIsoId >= 4
ROOT.gInterpreter.Declare("""
    using ROOT::RVecI;
    using ROOT::RVecB;
    RVecI compute_idxFirstCouple_charge_tightId_pfIsoId(const RVecI &charge, 
                          const RVecB &tightId, const RVecI &pfIsoId) {
        RVecI idx {-1, -1};
        bool exit = false;
        for (int i = 0; i < charge.size(); i++) {
            for (int j = i+1; j < charge.size(); j++) {
                if (charge[i] + charge[j] == 0 && 
                    tightId[i] == true && tightId[j] == true &&
                    pfIsoId[i] >= 4 && pfIsoId[j] >= 4) {
                    idx[0] = i;
                    idx[1] = j;
                    exit = true;
                    break;
                }
            }
            if (exit) break;
        }
        return idx;
    }
""")

# calcola gli indici della prima coppia di muoni con carica opposta
ROOT.gInterpreter.Declare("""
    using ROOT::RVecI;
    using ROOT::RVecB;
    RVecI compute_idxFirstCouple_charge(const RVecI &charge) {
        RVecI idx {-1, -1};
        bool exit = false;
        for (int i = 0; i < charge.size(); i++) {
            for (int j = i+1; j < charge.size(); j++) {
                if (charge[i] + charge[j] == 0) {
                    idx[0] = i;
                    idx[1] = j;
                    exit = true;
                    break;
                }
            }
            if (exit) break;
        }
        return idx;
    }
""")

# calcola la massa invariante della prima coppia di muoni con carica opposta
ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    using ROOT::RVecI;
    double compute_InvMass(const RVecD &pt, const RVecD &eta, const RVecD &phi,
                          const RVecD &mass, const RVecI &charge) {
        for (int i = 0; i < charge.size(); i++) {
            for (int j = i+1; j < charge.size(); j++) {
                if (charge[i] + charge[j] == 0) {
                    TLorentzVector v1, v2;
                    v1.SetPtEtaPhiM(pt[0], eta[0], phi[0], mass[0]); 
                    v2.SetPtEtaPhiM(pt[1], eta[1], phi[1], mass[1]);
                    return (v1 + v2).M();
                }
            }  
        }
        return -1;
    }
""")



# calcola gli indici dei jet VBF da un match geometrico con i quark VBF
ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    using ROOT::RVecI;
    RVecI compute_idxVBFJet(const RVecD &Jet_eta, const RVecD &Jet_phi, 
                          double q1_eta, double q1_phi, double q2_eta, 
                          double q2_phi, double dR_min=0.4) {
        RVecI idx {-1, -1};
        RVecD q = {q1_eta, q2_eta, q1_phi, q2_phi};
        double temp_dR = dR_min;
        for (int i = 0; i < 2; i++) {
            double temp_idx = -1;
            for (int j = 0; j < Jet_eta.size(); j++) {
                double dR = sqrt(pow(abs(q[i] - Jet_eta[j]), 2) + 
                            pow(DeltaPhiFunc(q[i+2], Jet_phi[j]), 2));
                if (dR < temp_dR) {
                    temp_dR = dR;
                    temp_idx = j;
                }
            }
            idx[i] = temp_idx;
            temp_dR = dR_min;
        }
        // if idx[0] == idx[1] idx[0] = -1;
        return idx;
    }
""")



# calcola gli indici dei muoni decaduti dall'Higgs usando le info del generatore
  
ROOT.gInterpreter.Declare("""
    using ROOT::RVecI;
    RVecI compute_idxMuon_gen(const RVecI &GenPart_pdgId, 
                    const RVecI &GenPart_genPartIdxMother,
                    const RVecI &GenPart_status) {

        RVecI idx = {-1, -1};
    
        // cerchiamo i muoni che hanno madre l'Higgs
        for (int i = 0; i < idx.size(); i++) {
            bool found = false;
            int j = 0;
            int idx_mother = -1;
            while (!found) {
                idx_mother = GenPart_genPartIdxMother[j];
                if (GenPart_pdgId[j] == pow(-1, i)*13 && 
                    GenPart_pdgId[idx_mother] == 25) {
                    idx[i] = j;
                    found = true;
                }
                j++;
            }
        }

        if (Any(idx == -1)) return idx;

        // cerchiamo le copie finali dei muoni trovati (status == 1)
        for (int i = 0; i < idx.size(); i++) {
            int j = idx[i] + 1;
            while(GenPart_status[idx[i]] != 1) {
                if (GenPart_pdgId[j] == pow(-1, i)*13 &&
                    GenPart_genPartIdxMother[j] == idx[i]) {
                    idx[i] = j;
                    j++;
                } else {
                    j++;
                    continue;
                }
            }
        }

        return idx;
    }
""")



# ordina gli indici delle gen particle in pt

ROOT.gInterpreter.Declare("""
    using ROOT::RVecI;
    using ROOT::RVecD;
    RVecI ptSort_GenPart(const RVecD &pt, int idx1, int idx2) {
        RVecI idx_sort {-1, -1};
        if (pt[idx1] > pt[idx2]) {
            idx_sort[0] = idx1;
            idx_sort[1] = idx2;
        } else if (pt[idx1] <  pt[idx2]) {
            idx_sort[0] = idx2;
            idx_sort[1] = idx1;
        }
        return idx_sort;
    }
""")




# calcola dR, Deta, pt per coppie di jet esclusa quella VBF

ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    using ROOT::RVecI;
    RVecD compute_DeltaR_allCouplesExVBF(const RVecD &eta, const RVecD &phi,
                          const RVecI &idx_vbf) {
        RVecD dR;
        for (int i = 0; i < eta.size(); i++) {
            for (int j = i+1; j < eta.size(); j++) {
                if (Any(idx_vbf == i) && Any(idx_vbf == j)) {
                    continue;
                }
                dR.push_back(compute_DeltaR(eta[i], phi[i], eta[j], phi[j]));
            }
        }
        return dR;
    }
""")

ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    using ROOT::RVecI;
    RVecD compute_InvMass_allCouplesExVBF(const RVecD &pt, const RVecD &eta, 
                          const RVecD &phi, const RVecD &mass, const RVecI &idx_vbf) {
        RVecD InvMass;
        TLorentzVector v1, v2;
        for (int i = 0; i < eta.size(); i++) {
            v1.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]); 
            for (int j = i+1; j < eta.size(); j++) {
                if (Any(idx_vbf == i) && Any(idx_vbf == j)) {
                    continue;
                }
                v2.SetPtEtaPhiM(pt[j], eta[j], phi[j], mass[j]);
                InvMass.push_back((v1 + v2).M());
            }
        }
        return InvMass;
    }
""")

ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    using ROOT::RVecI;
    RVecD compute_pt_allCouplesExVBF(const RVecD &pt, const RVecD &eta, 
                          const RVecD &phi, const RVecD &mass, const RVecI &idx_vbf) {
        RVecD pt_jj;
        TLorentzVector v1, v2;
        for (int i = 0; i < eta.size(); i++) {
            v1.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]); 
            for (int j = i+1; j < eta.size(); j++) {
                if (Any(idx_vbf == i) && Any(idx_vbf == j)) {
                    continue;
                }
                v2.SetPtEtaPhiM(pt[j], eta[j], phi[j], mass[j]);
                pt_jj.push_back((v1 + v2).Pt());
            }
        }
        return pt_jj;
    }
""")



# calcola la molteplicità di jet con btag

ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    int compute_nJet_btag(double WP, const RVecD &btag) {
        int njet = 0;
        for (int i = 0; i < btag.size(); i++) {
            if (btag[i] >= WP) njet++;
        }
        return njet;
    }
""")
ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    int compute_nJet_btag_pt_eta(double WP, const RVecD &btag, 
                          const RVecD &pt, const RVecD &eta) {
        int njet = 0;
        for (int i = 0; i < btag.size(); i++) {
            if (btag[i] >= WP && pt[i] > 25 && abs(eta[i]) < 2.5) njet++;
        }
        return njet;
    }
""")




# veto su btag come da paper

ROOT.gInterpreter.Declare("""
using ROOT::RVecD;
bool veto_btag(const RVecD &pt, const RVecD &eta, const RVecD &btag) {
    bool temp = true;
    for (int i = 0; i < pt.size(); i++) {
        // un jet tight WP (0.6734)
        if (btag[i] >= 0.6734 && pt[i] > 25 && abs(eta[i]) < 2.5) {
            temp = false;
            break;
        }
        // due jet medium WP (0.245)
        for (int j = 0; j < pt.size(); j++) {
            if (btag[i] >= 0.245 && btag[j] >= 0.245 &&
                pt[i] > 25 && pt[j] > 25 &&
                abs(eta[i]) < 2.5 && abs(eta[j]) < 2.5) {
                temp = false;
                break;
            }
        }
        if (!temp) break;
    }
    return temp;
}
""")




# veto su VH come da paper (circa)

ROOT.gInterpreter.Declare("""
    using ROOT::RVecD;
    using ROOT::RVecI;
    bool veto_VH(const RVecD &Muon_pt, const RVecD &Muon_eta, 
                const RVecD &Muon_Id, const RVecD &Muon_Iso,
                const RVecD &Electron_pt, const RVecD &Electron_eta, 
                const RVecD &Electron_Id, const RVecD &Electron_Iso,
                const RVecI &idx) {
        bool temp = true;
        for (int i = 0; i < Muon_pt.size(); i++) {
            if (Any(idx == i)) continue;
            // Muon_mvaMaID tight (2) and Muon_pfIsoId tight or more
            if (Muon_pt[i] > 20 && Muon_eta[i] < 2.4 &&
                Muon_Id[i] == 2 && Muon_Iso[i] >= 4) {
                temp = false;
                break;
            }
        }
        for (int i = 0; i < Electron_pt.size(); i++) {
            // Electron_mvaIso_WP90 and Electron_cutBased tight (4)
            if (Electron_pt[i] > 20 && Electron_eta[i] < 2.5 &&
                Electron_Id[i] == 4 && Electron_Iso[i] == 4) {
                temp = false;
                break;
            }
        }
        return temp; 
    }
""")






# calcola gli indici dei jet con il minor pt balance

ROOT.gInterpreter.Declare("""
using ROOT::RVecI;
using ROOT::RVecD;
RVecI compute_idxLowPtBalanceJets(const TLorentzVector &Dimuon_lv, 
        const RVecD &Jet_pt, const RVecD &Jet_eta, 
        const RVecD &Jet_phi, const RVecD &Jet_mass,
        const RVecD &Jet_btagPNetQvG, const double pt_bal_min=40, 
        const double pt=25, const double mass=250, const double qvg=0.19,
        const double Deta=0) {

    RVecI idx = {-1, -1};
    double pt_bal = 0;
    TLorentzVector v1, v2;

    for (int i = 0; i < Jet_pt.size(); i++) {

        if (Jet_pt[i] < pt) continue;
        if (Jet_btagPNetQvG[i] < qvg) continue;
        // if (abs(Jet_eta[i]) >= 4.7) continue;

        v1.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);

        for (int j = i+1; j < Jet_pt.size(); j++) {

            if (Jet_pt[j] < pt) continue;
            if (Jet_btagPNetQvG[j] < qvg) continue;
            // if (abs(Jet_eta[j]) >= 4.7) continue;

            v2.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);

            pt_bal = (v1 + v2 + Dimuon_lv).Pt();
            if (abs(Jet_eta[i] - Jet_eta[j]) < Deta) continue;
            if ((v1 + v2).M() < mass) continue;
            if (pt_bal < pt_bal_min) {
                idx[0] = i;
                idx[1] = j;
            }
        }
    }
    return idx;

}
""")



# calcola l'incertezza sulla massa invariante dei due muoni

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




# calcola l'incertezza sulla massa invariante dei due muoni utilizzando
# solo l'incertezza su pt

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





# ---------- COLLINS-SOPER REST FRAME -------------

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









# ------------------- PLOT --------------------------

fontsize_labelCMS = 0.04
fontsize_labelLum = 0.04
fontsize_legend = 0.038
fontsize_title = 0.04
fontsize_label = 0.038

Xtick_length = 0.017
Ytick_length = 0.017

Xtitle_offset = 1.1
Ytitle_offset = 1.6
Ztitle_offset = 1.1

canvas_leftMargin = 0.14
canvas_rightMargin = 0.08 # 0.03
canvas_topMargin = 0.06
canvas_bottomMargin = 0.13

h_linewidth = 2

x_posCMS = 0.17
y_posCMS = 0.87
x_posLum = 0.753
y_posLum = 0.96

# ritorna label di CMS e luminosità
def gen_CMSlabel(x_posCMS=x_posCMS, y_posCMS=y_posCMS, x_posLum=x_posLum, y_posLum=y_posLum, 
                 fontsize_CMS=fontsize_labelCMS,
                 fontsize_Lum=fontsize_labelLum):
    label = ROOT.TLatex() 
    label.SetNDC(True) # coord. normalizzate [0, 1]
    label.SetTextFont(42) # Helvetica
    label.SetTextSize(fontsize_CMS); label.DrawLatex(x_posCMS, y_posCMS, 
                            "#scale[1.3]{#bf{CMS}} #it{Simulation Private Work}")
    label.SetTextSize(fontsize_Lum); label.DrawLatex(x_posLum, y_posLum, "(13.6 TeV)")
    return label

# ritorna legend
def gen_legend(x1=0.65, y1=0.70, x2=0.82, y2=0.88, fontsize=fontsize_legend):
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(fontsize)
    return leg

# ritorna canvas
def gen_canvas(dim=[450, 450], rightMargin=canvas_rightMargin, 
               leftMargin=canvas_leftMargin,
               topMargin=canvas_topMargin,
               bottomMargin=canvas_bottomMargin):
    c = ROOT.TCanvas("", "", dim[0], dim[1])
    c.SetRightMargin(rightMargin)
    c.SetLeftMargin(leftMargin)
    c.SetBottomMargin(bottomMargin)
    c.SetTopMargin(topMargin)
    ROOT.gPad.SetTicks(1,1)
    ROOT.TGaxis.SetMaxDigits(3) 
    return c


def setup_histo(h, x_title, y_title, max):
    h.SetStats(False)
    h.GetXaxis().SetTitle(x_title)
    h.GetXaxis().SetTitleOffset(Xtitle_offset)
    h.GetXaxis().SetTitleSize(fontsize_title)
    h.GetXaxis().SetLabelSize(fontsize_label)
    h.GetXaxis().SetTickLength(Xtick_length)
    h.GetYaxis().SetTitle(y_title)
    h.GetYaxis().SetTitleSize(fontsize_title)
    h.GetYaxis().SetTitleOffset(Ytitle_offset)
    h.GetYaxis().SetLabelSize(fontsize_label)
    h.GetYaxis().SetTickLength(Ytick_length)
    h.SetMaximum(max*1.2)




def plotHisto1D_df(df, col, nbin, x_inf, x_sup, x_label, y_label, 
                   title="", filename="none", line_color=ROOT.kBlue,
                   path="~/Documents/tesi/fig/plot/", fill_color=0,
                   x_vline1='Inf', x_vline2='Inf'):

    c = gen_canvas()

    h = df.Histo1D(("", "", nbin, x_inf, x_sup), col)
    h.SetTitle(title)
    setup_histo(h, x_label, y_label, h.GetMaximum())
    h.SetLineColor(line_color)
    h.SetFillColor(fill_color)
    h.SetLineWidth(h_linewidth)
    h.Draw("HIST")
    
    label = gen_CMSlabel()

    if not math.isinf(float(x_vline1)):
        const1 = ROOT.TLine(x_vline1, h.GetMinimum(), x_vline1, 
                            h.GetMaximum() / 1.15)
        const1.SetLineColor(ROOT.kBlack)
        const1.SetLineWidth(2)
        const1.SetLineStyle(2)
        const1.Draw("L SAME")
        ROOT.SetOwnership(const1, False)
    if not math.isinf(float(x_vline2)):
        const2 = ROOT.TLine(x_vline2, h.GetMinimum(), x_vline2, 
                            h.GetMaximum())
        const2.SetLineColor(ROOT.kBlack)
        const2.SetLineWidth(2)
        const2.SetLineStyle(2)
        const2.Draw("L SAME")
        ROOT.SetOwnership(const2, False)
 
    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(h, False)
    ROOT.SetOwnership(label, False)



# Funzione per plot di due Histo1D dato un RDataFrame 

def plot2Histo1D_df(df1, df2, 
                   col1, col2, 
                   nbin1, x_inf1, x_sup1, 
                   nbin2, x_inf2, x_sup2, 
                   x_label, y_label,
                   name1_leg, name2_leg,
                   x1_leg=0.65, y1_leg=0.70, x2_leg=0.82, y2_leg=0.88,
                   color1=ROOT.kRed, color2=ROOT.kBlue,
                   title="", filename="none", 
                   path="~/Documents/tesi/fig/plot/"):

    c = gen_canvas()

    h1 = df1.Histo1D(("", "", nbin1, x_inf1, x_sup1), col1)
    h2 = df2.Histo1D(("", "", nbin2, x_inf2, x_sup2), col2)
    h1.SetTitle(title)
    setup_histo(h1, x_label, y_label, max(h1.GetMaximum(), h2.GetMaximum()))
    h1.SetLineColor(color1)
    h1.SetLineWidth(h_linewidth)
    h2.SetLineColor(color2)
    h2.SetLineWidth(h_linewidth)
    h1.Draw("HIST")
    h2.Draw("HIST SAME")

    leg = ROOT.TLegend(x1_leg, y1_leg, x2_leg, y2_leg)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(fontsize_legend)
    leg.AddEntry(h1.GetValue(), name1_leg, "l")
    leg.AddEntry(h2.GetValue(), name2_leg, "l")
    leg.Draw()
    
    label = gen_CMSlabel()
 
    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(h1, False)
    ROOT.SetOwnership(h2, False)
    ROOT.SetOwnership(leg, False)
    ROOT.SetOwnership(label, False)

"""
# Funzione per plot di tre Histo1D dato un RDataFrame 

def plot3Histo1D_df(df1, df2, df3,
                   col1, col2, col3,
                   nbin1, x_inf1, x_sup1, 
                   nbin2, x_inf2, x_sup2, 
                   nbin3, x_inf3, x_sup3,
                   x_label, y_label,
                   name1_leg, name2_leg, name3_leg,
                   x1_leg=0.65, y1_leg=0.70, x2_leg=0.82, y2_leg=0.88,
                   x_posCMS=0.12, y_posCMS=0.84, 
                   color1=ROOT.kRed, color2=ROOT.kBlue, color3=ROOT.kOrange,
                   sim=True, title="", filename="none", path=""):

    c = ROOT.TCanvas("", "", 450, 350)
    #c.SetLeftMargin(0.1)
    ROOT.gPad.SetTicks(1,1)

    h1 = df1.Histo1D(("", "", nbin1, x_inf1, x_sup1), col1)
    h2 = df2.Histo1D(("", "", nbin2, x_inf2, x_sup2), col2)
    h3 = df3.Histo1D(("", "", nbin3, x_inf3, x_sup3), col3)
    h1.SetTitle(title)
    h1.GetXaxis().SetTitle(x_label)
    h1.GetXaxis().SetTickLength(Xtick_length)
    h1.GetXaxis().SetTitleOffset(Xtitle_offset)
    h1.GetXaxis().SetTitleSize(fontsize_general)
    h1.GetYaxis().SetTitle(y_label)
    h1.GetYaxis().SetTickLength(Ytick_length)
    h1.GetYaxis().SetTitleSize(fontsize_general)
    h1.SetMaximum(max(h1.GetMaximum(), h2.GetMaximum(), h3.GetMaximum())*1.15)
    h1.SetStats(False)
    h1.SetLineColor(color1)
    h2.SetLineColor(color2)
    h3.SetLineColor(color3)
    h1.Draw("HIST")
    h2.Draw("HIST SAME")
    h3.Draw("HIST SAME")

    leg = ROOT.TLegend(x1_leg, y1_leg, x2_leg, y2_leg)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(fontsize_general)
    leg.AddEntry(h1.GetValue(), name1_leg, "l")
    leg.AddEntry(h2.GetValue(), name2_leg, "l")
    leg.AddEntry(h3.GetValue(), name3_leg, "l")
    leg.Draw()
    
    label = ROOT.TLatex() 
    label.SetNDC(True)
    label.SetTextFont(42) # Helvetica
    label.SetTextSize(fontsize_labelCMS); label.DrawLatex(x_posCMS, y_posCMS, 
                                             "#it{Private work} (#scale[1.15]{#bf{CMS}} #it{Simulation})")
    if sim:
        label.SetTextSize(fontsize_labelLum); label.DrawLatex(0.76, 0.92, "(13.6 TeV)")
    else:
        label.SetTextSize(fontsize_labelLum); label.DrawLatex(0.65, 0.92, "59.8 fb^{-1} (13.6 TeV)")
 
    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(h1, False)
    ROOT.SetOwnership(h2, False)
    ROOT.SetOwnership(h3, False)
    ROOT.SetOwnership(leg, False)
"""



# plot N Histo1D

def plotNHisto1D_df(dfs, cols, bins, x_infs, x_sups, 
                        x_label, y_label, leg_names, colors,
                        x1_leg=0.65, y1_leg=0.70, x2_leg=0.82, y2_leg=0.88,
                        filename="none", 
                        path="~/Documents/tesi/fig/plot/"):
    
    c = gen_canvas()

    leg = gen_legend(x1_leg, y1_leg, x2_leg, y2_leg)

    histos = []

    for i in range(len(dfs)):
        histos.append(dfs[i].Histo1D(("", "", bins[i], x_infs[i], x_sups[i]), cols[i]))
        histos[i].SetLineColor(colors[i])
        histos[i].SetLineWidth(h_linewidth)
        leg.AddEntry(histos[i].GetValue(), leg_names[i], "l")

    setup_histo(histos[0], x_label, y_label, max(histos[i].GetMaximum() \
                             for i in range(len(histos))))
    histos[0].Draw("HIST")
    for i in range(1, len(histos)):
        histos[i].Draw("HIST SAME")

    leg.Draw()

    label = gen_CMSlabel()
    label.Draw()

    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))
    
    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(leg, False)
    ROOT.SetOwnership(label, False)
    for i in range(len(histos)):
        ROOT.SetOwnership(histos[i], False)




# plot N Histo1D normalizzati 

def plotNHisto1DNorm_df(dfs, cols, bins, x_infs, x_sups, 
                        x_label, y_label, leg_names, colors,
                        x1_leg=0.65, y1_leg=0.70, x2_leg=0.82, y2_leg=0.88,
                        filename="none", 
                        path="~/Documents/tesi/fig/plot/"):
    
    c = gen_canvas()

    leg = gen_legend(x1_leg, y1_leg, x2_leg, y2_leg)

    histos = []

    for i in range(len(dfs)):
        histos.append(dfs[i].Histo1D(("", "", bins[i],
                        x_infs[i], x_sups[i]), cols[i]))
        histos[i].Scale(1.0 / histos[i].Integral())
        histos[i].SetLineColor(colors[i])
        histos[i].SetLineWidth(h_linewidth)
        leg.AddEntry(histos[i].GetValue(), leg_names[i], "l")


    setup_histo(histos[0], x_label, y_label, max(histos[i].GetMaximum() \
                             for i in range(len(histos))))
    histos[0].SetMinimum(min(histos[i].GetMinimum() \
                             for i in range(len(histos)))*0.9)
    histos[0].Draw("HIST")
    for i in range(1, len(histos)):
        histos[i].Draw("HIST SAME")

    leg.Draw()

    label = gen_CMSlabel()
    label.Draw()

    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))
    
    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(leg, False)
    ROOT.SetOwnership(label, False)
    for i in range(len(histos)):
        ROOT.SetOwnership(histos[i], False)

# Funzione per il plot Histo2D dato un RDataFrame

def plotHisto2D_df(df, 
                   col1, col2, 
                   nbin1, x_inf1, x_sup1,
                   nbin2, x_inf2, x_sup2,
                   x_label, y_label,
                   colz_title, mode="COLZ",
                   x_posCMS=0.14, y_posCMS=0.94,
                   title="", filename="none", 
                   path="~/Documents/tesi/fig/plot/"):

    c = gen_canvas(leftMargin=0.14, rightMargin=0.16, 
                   bottomMargin=0.11, topMargin=0.08)

    h = df.Histo2D(("", "", nbin1, x_inf1, x_sup1, nbin2, x_inf2, x_sup2), col1, col2)
    h.SetTitle(title)
    h.GetXaxis().SetTitle(x_label)
    h.GetXaxis().SetTitleSize(fontsize_title)
    h.GetXaxis().SetLabelSize(fontsize_label)
    h.GetXaxis().SetTitleOffset(Xtitle_offset)
    h.GetYaxis().SetTitle(y_label)
    h.GetYaxis().SetTitleSize(fontsize_title)
    h.GetYaxis().SetTitleOffset(Ytitle_offset)
    h.GetYaxis().SetLabelSize(fontsize_label)
    if mode == "COLZ":
        h.GetZaxis().SetTitle(colz_title)
        h.GetZaxis().SetTitleSize(fontsize_title)
        h.GetZaxis().SetTitleOffset(Ztitle_offset)
    ROOT.gStyle.SetPalette(109) # kViridis: 112, kInvertedDarkBodyRadiator: 56, kLightTemperature: 87
    h.SetStats(0)
    h.Draw(mode)

    label = gen_CMSlabel(x_posCMS=x_posCMS, y_posCMS=y_posCMS, 
                         x_posLum=0.67, y_posLum=0.94)
    #label = ROOT.TLatex() 
    #label.SetNDC(True)
    #label.SetTextFont(42) # Helvetica
    #label.SetTextSize(fontsize_labelCMS); label.DrawLatex(x_posCMS, y_posCMS, 
                                             #"#it{Private work} (#scale[1.15]{#bf{CMS}} #it{Simulation})")
    #if h.GetMaximum() > 1e3:
        #if sim:
            #label.SetTextSize(fontsize_labelLum); label.DrawLatex(0.68, 0.92, "(13.6 TeV)")
        #else:
            #label.SetTextSize(fontsize_labelLum); label.DrawLatex(0.59, 0.92, "59.8 fb^{-1} (13.6 TeV)")
    #else:
        #if sim:
            #label.SetTextSize(fontsize_labelLum); label.DrawLatex(0.72, 0.92, "(13.6 TeV)")
        #else:
            #label.SetTextSize(fontsize_labelLum); label.DrawLatex(0.59, 0.92, "59.8 fb^{-1} (13.6 TeV)")

    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(h, False)
    ROOT.SetOwnership(label, False)







# ritorna hist e canvas per il fit

def returnHisto1D_df(df, col, nbin, x_inf, x_sup, x_label, y_label, 
                   x_posCMS=0.16, y_posCMS=0.87):

    c = gen_canvas()

    h = df.Histo1D(("", "", nbin, x_inf, x_sup), col)
    setup_histo(h, x_label, y_label, h.GetMaximum())
    h.SetLineColor(ROOT.kBlue)
    h.Draw()
    
    label = gen_CMSlabel(x_posCMS=x_posCMS, y_posCMS=y_posCMS)

    return c, h, label

def fitHisto1D_df(h, fit_func, fit_inf, fit_sup, fit_opt="SNQR", color=ROOT.kGreen, 
                  line_width=2, init_par=[]):
    f = ROOT.TF1("f", fit_func, fit_inf, fit_sup)
    #f = ROOT.TF1("f", fit_func, fit_inf - (fit_inf*0.1), fit_sup + (fit_sup*0.1))
    if len(init_par) != 0:
        for i, par in enumerate(init_par):
            f.SetParameter(i, par)
    f.SetLineColor(color)
    f.SetLineWidth(line_width)
    f.SetLineStyle(2)
    #h.Fit("f", fit_opt, "", fit_inf, fit_sup)
    fitRes = h.Fit("f", fit_opt, "", fit_inf, fit_sup)
    
    return f, fitRes







# ratio plot di un histo1D da rdataframe
def plot2Histo1D_Ratio(df1, df2, col1, col2, nbins, x_inf, x_sup, x_label, 
                        y_label_upper, y_label_lower, leg_name1, leg_name2, 
                        x1_leg=0.63, y1_leg=0.69, x2_leg=0.80, y2_leg=0.88,
                        x_labelCMS=0.16, y_labelCMS=0.87, x_labelLum=0.78, 
                        y_labelLum=0.96, sc1=1.3, sc2=3,
                        filename="none", path="~/Documents/tesi/fig/plot/"):
    
    c = gen_canvas()

    h_all = df1.Histo1D(("", "", nbins, x_inf, x_sup), col1)
    h_HLT = df2.Histo1D(("", "", nbins, x_inf, x_sup), col2)

    h_all.SetMaximum(max(h_all.GetMaximum(), h_HLT.GetMaximum())*1.2)
    h_all.SetStats(False)

    h_all.GetXaxis().SetLabelSize(0)
    h_all.GetXaxis().SetTitleSize(fontsize_title * sc1)
    h_all.GetXaxis().SetTickLength(Xtick_length * sc1)

    h_all.GetYaxis().SetLabelSize(fontsize_label * sc1)
    h_all.GetYaxis().SetTickLength(Ytick_length)
    h_all.GetYaxis().SetTitle(y_label_upper)
    h_all.GetYaxis().SetTitleSize(fontsize_title * sc1)

    h_HLT.SetLineColor(ROOT.kRed)
    h_HLT.SetLineWidth(h_linewidth)
    h_all.SetLineWidth(h_linewidth)

    pad1 = ROOT.TPad("pad1", "Upper pad", 0, 0.3, 1, 1)
    pad1.Draw()
    pad1.SetBottomMargin(0.02)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.cd()

    h_all.Draw("HIST")
    h_HLT.Draw("HIST SAME")

    label = gen_CMSlabel(x_labelCMS, y_labelCMS, x_labelLum, y_labelLum, 
                         fontsize_labelCMS*sc1, fontsize_label*sc1)
    label.Draw()

    leg1 = ROOT.TLegend(x1_leg, y1_leg, x2_leg, y2_leg)
    leg1.SetFillColor(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(fontsize_legend * sc1)
    leg1.AddEntry(h_all.GetValue(), leg_name1, "l")
    leg1.AddEntry(h_HLT.GetValue(), leg_name2, "l")
    leg1.Draw()
 
    c.cd()

    pad2 = ROOT.TPad("pad2", "Lower pad", 0, 0, 1, 0.3)
    pad2.Draw()
    pad2.SetTopMargin(0.015)
    pad2.SetBottomMargin(0.3)
    pad2.cd()
    pad2.SetTickx(1)
    pad2.SetTicky(1)

    # più semplice utilizzando .Divide()
    h_ratio = ROOT.TH1D("", "", nbins, x_inf, x_sup)
    for i in range(1, nbins+1):
        ratio = 0
        if h_all.GetBinContent(i) > 0:
            ratio = h_HLT.GetBinContent(i) / h_all.GetBinContent(i)
        elif h_all.GetBinContent(i) <= 0 and h_HLT.GetBinContent(i) > 0:
            ratio = 1
        h_ratio.SetBinContent(i, ratio)

    h_ratio.SetLineColor(ROOT.kGray+3)
    h_ratio.SetLineWidth(h_linewidth)
    h_ratio.GetXaxis().SetTitle(x_label)
    h_ratio.GetXaxis().SetTitleSize(fontsize_title * sc2)
    h_ratio.GetXaxis().SetTitleOffset(1.1)
    h_ratio.GetXaxis().SetLabelOffset(0.02)
    h_ratio.GetXaxis().SetTickLength(Xtick_length * sc2)
    h_ratio.GetXaxis().SetLabelSize(fontsize_label * sc2)

    h_ratio.GetYaxis().SetTitle(y_label_lower)
    h_ratio.GetYaxis().SetTitleOffset(0.35)
    h_ratio.GetYaxis().SetTitleSize(fontsize_title * sc2)
    h_ratio.GetYaxis().SetTickLength(Ytick_length)
    h_ratio.GetYaxis().SetLabelSize(fontsize_label * sc2)
    #h_ratio.GetYaxis().SetRangeUser(0, 1)
    h_ratio.GetYaxis().SetNdivisions(510)

    h_ratio.SetStats(False)
    h_ratio.Draw("HIST")

    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(h_all, False)
    ROOT.SetOwnership(h_HLT, False)
    ROOT.SetOwnership(h_ratio, False)
    ROOT.SetOwnership(label, False)
    ROOT.SetOwnership(leg1, False)



# ratio plot di N Histo1D da rdataframe
def plotNHisto1D_Ratio(dfs, cols, nbins, x_inf, x_sup, x_label, 
                        y_label_upper, y_label_lower, leg_names, colors,
                        x1_leg=0.70, y1_leg=0.67, x2_leg=0.90, y2_leg=0.88,
                        x_labelCMS=0.17, y_labelCMS=0.82, x_labelLum=0.753, 
                        y_labelLum=0.95, sc1=1.43, sc2=3.3, 
                        filename="none", 
                        path="~/Documents/tesi/fig/plot/"):

    c = gen_canvas()

    leg = gen_legend(x1_leg, y1_leg, x2_leg, y2_leg, fontsize_legend * sc1)

    histos = []
    for i in range(len(dfs)):
        histos.append(dfs[i].Histo1D(("", "", nbins, x_inf, x_sup), cols[i]))
        histos[i].SetLineColor(colors[i])
        histos[i].SetLineWidth(h_linewidth)
        leg.AddEntry(histos[i].GetValue(), leg_names[i], "l")

    histos[0].SetMaximum(max(histos[i].GetMaximum()\
                              for i in range(len(histos)))*1.2)
    histos[0].SetMinimum(0)
    histos[0].SetStats(False)
    histos[0].GetXaxis().SetLabelSize(0)
    histos[0].GetXaxis().SetTitleSize(fontsize_title * sc1)
    histos[0].GetXaxis().SetTickLength(Xtick_length * sc1)
    histos[0].GetYaxis().SetLabelSize(fontsize_label * sc1)
    histos[0].GetYaxis().SetTickLength(Ytick_length)
    histos[0].GetYaxis().SetTitle(y_label_upper)
    histos[0].GetYaxis().SetTitleSize(fontsize_title * sc1)
    #histos[0].GetYaxis().SetTitleOffset(Ytitle_offset * sc1)


    pad1 = ROOT.TPad("pad1", "Upper pad", 0, 0.3, 1, 1)
    pad1.Draw()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(canvas_leftMargin)
    pad1.SetRightMargin(canvas_rightMargin)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.cd()

    histos[0].Draw("HIST")
    for i in range(1, len(histos)):
        histos[i].Draw("HIST SAME")

    label = gen_CMSlabel(x_labelCMS, y_labelCMS, x_labelLum, y_labelLum,
                         fontsize_labelCMS * sc1, fontsize_labelLum * sc1)
    label.Draw()

    leg.Draw()
 
    c.cd()

    pad2 = ROOT.TPad("pad2", "Lower pad", 0, 0, 1, 0.3)
    pad2.Draw()
    pad2.SetTopMargin(0.015)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(canvas_leftMargin)
    pad2.SetRightMargin(canvas_rightMargin)
    pad2.cd()
    pad2.SetTickx(1)
    pad2.SetTicky(1)

    # più semplice utilizzando .Divide()
    ratios = []
    for i in range(1, len(histos)):
        ratios.append(ROOT.TH1D("", "", nbins, x_inf, x_sup))
        ratios[i-1].SetLineColor(colors[i])
        ratios[i-1].SetLineWidth(h_linewidth)
        for j in range(1, nbins+1):
            ratio = 0
            if histos[0].GetBinContent(j) > 0:
                ratio = histos[i].GetBinContent(j) / histos[0].GetBinContent(j)
            elif histos[0].GetBinContent(j) <= 0 \
                and histos[i].GetBinContent(j) > 0:
                ratio = 1
            ratios[i-1].SetBinContent(j, ratio)

    ratios[0].SetMaximum(max(ratios[i].GetMaximum() \
                             for i in range(len(ratios)))*1.01)
    ratios[0].SetMinimum(min(ratios[i].GetMinimum() \
                             for i in range(len(ratios)))*0.99)
    ratios[0].SetStats(False)
    
    ratios[0].GetXaxis().SetTitle(x_label)
    ratios[0].GetXaxis().SetTitleSize(fontsize_title * sc2)
    ratios[0].GetXaxis().SetTitleOffset(1.1)
    ratios[0].GetXaxis().SetLabelOffset(0.02)
    ratios[0].GetXaxis().SetTickLength(Xtick_length * sc2)
    ratios[0].GetXaxis().SetLabelSize(fontsize_label * sc2)

    ratios[0].GetYaxis().SetTitle(y_label_lower)
    ratios[0].GetYaxis().SetTitleOffset(0.45)
    ratios[0].GetYaxis().SetTitleSize(fontsize_title * sc2)
    ratios[0].GetYaxis().SetTickLength(Ytick_length)
    ratios[0].GetYaxis().SetLabelSize(fontsize_label * sc2)
    ratios[0].GetYaxis().SetNdivisions(505)

    ratios[0].Draw("HIST")
    for i in range(1, len(ratios)):
        ratios[i].Draw("HIST SAME")

    const = ROOT.TLine(x_inf, 1, x_sup, 1)
    const.SetLineColor(colors[0])
    const.SetLineWidth(1)
    const.SetLineStyle(2)
    const.Draw("L SAME")

    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(leg, False)
    ROOT.SetOwnership(label, False)
    ROOT.SetOwnership(const, False)
    for i in range(len(histos)):
        ROOT.SetOwnership(histos[i], False)
    for i in range(len(ratios)):
        ROOT.SetOwnership(ratios[i], False)




# ratio plot di N Histo1D normalizzati da rdataframe
def plotNHisto1DNorm_Ratio(dfs, cols, nbins, x_inf, x_sup, x_label, 
                        y_label_upper, y_label_lower, leg_names, colors, 
                        colors_ratio, x1_leg=0.70, y1_leg=0.67, x2_leg=0.90, 
                        y2_leg=0.88, x_labelCMS=0.17, y_labelCMS=0.82, 
                        x_labelLum=0.753, y_labelLum=0.95, sc1=1.43, sc2=3.3, 
                        filename="none", fill_style=0,
                        path="~/Documents/tesi/fig/plot/"):

    c = gen_canvas()

    leg = gen_legend(x1_leg, y1_leg, x2_leg, y2_leg, fontsize_legend * sc1)

    histos = []
    for i in range(len(dfs)):
        histos.append(dfs[i].Histo1D(("", "", nbins, x_inf, x_sup), cols[i]))
        if fill_style == 0:
            histos[i].SetLineColor(colors[i])
            histos[i].SetLineWidth(h_linewidth)
            histos[i].Scale(1. / histos[i].Integral())
            leg.AddEntry(histos[i].GetValue(), leg_names[i], "l")
        elif fill_style == 1:
            histos[i].SetLineColor(ROOT.kBlack)
            histos[i].SetLineWidth(1)
            histos[i].Scale(1. / histos[i].Integral())
            histos[i].SetFillColor(colors[i])
            histos[i].SetFillStyle(3002+i*3)
            leg.AddEntry(histos[i].GetValue(), leg_names[i], "f")

    histos[0].SetMaximum(max(histos[i].GetMaximum()\
                              for i in range(len(histos)))*1.2)
    histos[0].SetMinimum(0)
    histos[0].SetStats(False)
    histos[0].GetXaxis().SetLabelSize(0)
    histos[0].GetXaxis().SetTitleSize(fontsize_title * sc1)
    histos[0].GetXaxis().SetTickLength(Xtick_length * sc1)
    histos[0].GetYaxis().SetLabelSize(fontsize_label * sc1)
    histos[0].GetYaxis().SetTickLength(Ytick_length)
    histos[0].GetYaxis().SetTitle(y_label_upper)
    histos[0].GetYaxis().SetTitleSize(fontsize_title * sc1)


    pad1 = ROOT.TPad("pad1", "Upper pad", 0, 0.3, 1, 1)
    pad1.Draw()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(canvas_leftMargin)
    pad1.SetRightMargin(canvas_rightMargin)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.cd()

    histos[0].Draw("HIST")
    for i in range(1, len(histos)):
        histos[i].Draw("HIST SAME")

    label = gen_CMSlabel(x_labelCMS, y_labelCMS, x_labelLum, y_labelLum,
                         fontsize_labelCMS * sc1, fontsize_labelLum * sc1)
    label.Draw()

    leg.Draw()
 
    c.cd()

    pad2 = ROOT.TPad("pad2", "Lower pad", 0, 0, 1, 0.3)
    pad2.Draw()
    pad2.SetTopMargin(0.015)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(canvas_leftMargin)
    pad2.SetRightMargin(canvas_rightMargin)
    pad2.cd()
    pad2.SetTickx(1)
    pad2.SetTicky(1)

    # più semplice utilizzando .Divide()
    ratios = []
    for i in range(1, len(histos)):
        ratios.append(ROOT.TH1D("", "", nbins, x_inf, x_sup))
        ratios[i-1].SetLineColor(colors_ratio[i])
        ratios[i-1].SetLineWidth(h_linewidth)
        for j in range(1, nbins+1):
            ratio = 0
            if histos[0].GetBinContent(j) > 0:
                ratio = histos[i].GetBinContent(j) / histos[0].GetBinContent(j)
            elif histos[0].GetBinContent(j) <= 0 \
                and histos[i].GetBinContent(j) > 0:
                ratio = 1
            ratios[i-1].SetBinContent(j, ratio)

    
    ratios[0].SetMaximum(max(ratios[i].GetMaximum() \
                             for i in range(len(ratios)))*1.01)
    ratios[0].SetMinimum(0)
    ratios[0].SetMinimum(min(ratios[i].GetMinimum() \
                             for i in range(len(ratios)))*0.99)
    ratios[0].SetStats(False)
    
    ratios[0].GetXaxis().SetTitle(x_label)
    ratios[0].GetXaxis().SetTitleSize(fontsize_title * sc2)
    ratios[0].GetXaxis().SetTitleOffset(1.1)
    ratios[0].GetXaxis().SetLabelOffset(0.02)
    ratios[0].GetXaxis().SetTickLength(Xtick_length * sc2)
    ratios[0].GetXaxis().SetLabelSize(fontsize_label * sc2)

    ratios[0].GetYaxis().SetTitle(y_label_lower)
    ratios[0].GetYaxis().SetTitleOffset(0.45)
    ratios[0].GetYaxis().SetTitleSize(fontsize_title * sc2)
    ratios[0].GetYaxis().SetTickLength(Ytick_length)
    ratios[0].GetYaxis().SetLabelSize(fontsize_label * sc2)
    ratios[0].GetYaxis().SetNdivisions(505)

    ratios[0].Draw("HIST")
    for i in range(1, len(ratios)):
        ratios[i].Draw("HIST SAME")

    const = ROOT.TLine(x_inf, 1, x_sup, 1)
    const.SetLineColor(colors_ratio[0])
    const.SetLineWidth(2)
    const.SetLineStyle(2)
    const.Draw("L SAME")

    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(leg, False)
    ROOT.SetOwnership(label, False)
    ROOT.SetOwnership(const, False)
    for i in range(len(histos)):
        ROOT.SetOwnership(histos[i], False)
    for i in range(len(ratios)):
        ROOT.SetOwnership(ratios[i], False)


# fit di un histo1D da rdataframe e plot del pull
def plotHisto1D_fitAndPull(df, col, nbins, x_inf, x_sup, fit_func, fit_inf, 
                        fit_sup, x_label, y_label_upper, y_label_lower, 
                        leg_name_hist, leg_name_fit, init_par=[], yie=False, 
                        fit_opt="SNQR", fit_color=ROOT.kRed, h_color=ROOT.kP6Blue,
                        x1_leg=0.70, y1_leg=0.67, x2_leg=0.90, y2_leg=0.88, 
                        x_labelCMS=0.17, y_labelCMS=0.82, x_labelLum=0.753, 
                        y_labelLum=0.95, sc1=1.43, sc2=3.3,
                        filename="none", path="~/Documents/tesi/fig/plot/"):
    c = gen_canvas()

    pad1 = ROOT.TPad("pad1", "Upper pad", 0.0, 0.3, 1.0, 1.0)
    pad1.Draw()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(canvas_leftMargin)
    pad1.SetRightMargin(canvas_rightMargin)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.cd()

    h = df.Histo1D(("", "", nbins, x_inf, x_sup), col)
    h.SetMaximum(h.GetMaximum()*1.2)
    h.SetStats(False)
    h.GetXaxis().SetLabelSize(0)
    h.GetXaxis().SetTitleSize(fontsize_label * sc1)
    h.GetXaxis().SetTickLength(Xtick_length * sc1)
    h.GetYaxis().SetLabelSize(fontsize_label * sc1)
    h.GetYaxis().SetTickLength(Ytick_length)
    h.GetYaxis().SetTitle(y_label_upper)
    h.GetYaxis().SetTitleSize(fontsize_title * sc1)
    #h.GetYaxis().SetTitleOffset(Ytitle_offset * sc1)
    h.SetLineWidth(h_linewidth)
    h.SetLineColor(h_color)
    h.Draw("HIST")

    f, fitRes = fitHisto1D_df(h, fit_func, fit_inf, fit_sup, fit_opt, fit_color,
                              3, init_par)
    f.Draw("SAME")
    if yie:
        integral = f.Integral(fit_inf, fit_sup)
        integral_err = f.IntegralError(fit_inf, fit_sup, f.GetParameters(), 
                                    fitRes.GetCovarianceMatrix().GetMatrixArray())
    else:
        integral = 0
        integral_err = 0


    label = gen_CMSlabel(x_labelCMS, y_labelCMS, x_labelLum, y_labelLum, 
                         fontsize_labelCMS*sc1, fontsize_labelLum*sc1)
    label.Draw()

    leg1 = ROOT.TLegend(x1_leg, y1_leg, x2_leg, y2_leg)
    leg1.SetFillColor(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(fontsize_legend * sc1)
    leg1.AddEntry(h.GetValue(), leg_name_hist, "l")
    leg1.AddEntry(f, leg_name_fit, "l")
    leg1.Draw()


    c.cd()


    pad2 = ROOT.TPad("pad2", "Upper pad", 0.0, 0.0, 1.0, 0.3)
    pad2.Draw()
    pad2.SetTopMargin(0.015)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(canvas_leftMargin)
    pad2.SetRightMargin(canvas_rightMargin)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.cd()

    pull = ROOT.TH1F("", "", h.GetNbinsX(), 
                    h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    pull.SetStats(False)
    pull.SetLineColor(ROOT.kGray+2)
    pull.SetLineWidth(2)
    pull.GetXaxis().SetTitle(x_label)
    pull.GetXaxis().SetTitleSize(fontsize_title * sc2)
    pull.GetXaxis().SetTitleOffset(1.1)
    pull.GetXaxis().SetLabelOffset(0.02)
    pull.GetXaxis().SetTickLength(Xtick_length * sc2)
    pull.GetXaxis().SetLabelSize(fontsize_label * sc2)
    pull.GetYaxis().SetTitle(y_label_lower)
    pull.GetYaxis().SetTitleOffset(0.35)
    pull.GetYaxis().SetTitleSize(fontsize_title * sc2)
    pull.GetYaxis().SetTickLength(Ytick_length)
    pull.GetYaxis().SetLabelSize(fontsize_label * sc2)
    pull.GetYaxis().SetNdivisions(505)

    for i in range(1, h.GetNbinsX()+1):
        obs = h.GetBinContent(i)
        err = h.GetBinError(i)
        x = h.GetBinCenter(i)
        if x >= fit_inf and x <= fit_sup:
            fit_val = f.Eval(x)
            if err > 0:
                p = (obs - fit_val) / err
                pull.SetBinContent(i, p)
        else:
            pull.SetBinContent(i, 0)

    pull.Draw("HIST")

    const = ROOT.TLine(x_inf, 0, x_sup, 0)
    const.SetLineColor(ROOT.kBlack)
    const.SetLineWidth(2)
    const.SetLineStyle(2)
    const.Draw("L SAME")


    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(h, False)
    ROOT.SetOwnership(f, False)
    ROOT.SetOwnership(pull, False)
    ROOT.SetOwnership(label, False)
    ROOT.SetOwnership(leg1, False)
    ROOT.SetOwnership(const, False)

    return fitRes, integral, integral_err



# ratio plot a coppie di N Histo1D da rdataframe
def plotNHisto1D_RatioCouples(dfs, cols, nbins, x_inf, x_sup, x_label, 
                        y_label_upper, y_label_lower, leg_names, colors, 
                        colors_ratio, x1_leg=0.70, y1_leg=0.67, 
                        x2_leg=0.90, y2_leg=0.88, x_labelCMS=0.17, 
                        y_labelCMS=0.82, x_labelLum=0.753, 
                        y_labelLum=0.95, sc1=1.43, sc2=3.3, 
                        filename="none", 
                        path="~/Documents/tesi/fig/plot/"):

    c = gen_canvas()

    leg = gen_legend(x1_leg, y1_leg, x2_leg, y2_leg, fontsize_legend * sc1)

    histos = []
    for i in range(len(dfs)):
        histos.append(dfs[i].Histo1D(("", "", nbins, x_inf, x_sup), cols[i]))
        histos[i].SetLineColor(colors[i])
        histos[i].SetLineWidth(h_linewidth)
        leg.AddEntry(histos[i].GetValue(), leg_names[i], "l")

    histos[0].SetMaximum(max(histos[i].GetMaximum()\
                              for i in range(len(histos)))*1.2)
    histos[0].SetMinimum(0)
    histos[0].SetStats(False)
    histos[0].GetXaxis().SetLabelSize(0)
    histos[0].GetXaxis().SetTitleSize(fontsize_title * sc1)
    histos[0].GetXaxis().SetTickLength(Xtick_length * sc1)
    histos[0].GetYaxis().SetLabelSize(fontsize_label * sc1)
    histos[0].GetYaxis().SetTickLength(Ytick_length)
    histos[0].GetYaxis().SetTitle(y_label_upper)
    histos[0].GetYaxis().SetTitleSize(fontsize_title * sc1)
    #histos[0].GetYaxis().SetTitleOffset(Ytitle_offset * sc1)


    pad1 = ROOT.TPad("pad1", "Upper pad", 0, 0.3, 1, 1)
    pad1.Draw()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(canvas_leftMargin)
    pad1.SetRightMargin(canvas_rightMargin)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.cd()

    histos[0].Draw("HIST")
    for i in range(1, len(histos)):
        histos[i].Draw("HIST SAME")

    label = gen_CMSlabel(x_labelCMS, y_labelCMS, x_labelLum, y_labelLum,
                         fontsize_labelCMS * sc1, fontsize_labelLum * sc1)
    label.Draw()

    leg.Draw()
 
    c.cd()

    pad2 = ROOT.TPad("pad2", "Lower pad", 0, 0, 1, 0.3)
    pad2.Draw()
    pad2.SetTopMargin(0.015)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(canvas_leftMargin)
    pad2.SetRightMargin(canvas_rightMargin)
    pad2.cd()
    pad2.SetTickx(1)
    pad2.SetTicky(1)

    # più semplice utilizzando .Divide()
    ratios = []
    for i in range(0, int(len(histos)/2)):
        ratios.append(ROOT.TH1D("", "", nbins, x_inf, x_sup))
        ratios[i].SetLineColor(colors_ratio[i+1])
        ratios[i].SetLineWidth(h_linewidth)
        for j in range(1, nbins+1):
            ratio = 0
            if histos[i*2].GetBinContent(j) > 0:
                ratio = histos[i*2+1].GetBinContent(j) / histos[i*2].GetBinContent(j)
            elif histos[i*2].GetBinContent(j) <= 0 \
                and histos[i*2+1].GetBinContent(j) > 0:
                ratio = 1
            ratios[i].SetBinContent(j, ratio)

    ratios[0].SetMaximum(max(ratios[i].GetMaximum() \
                             for i in range(len(ratios)))*1.01)
    ratios[0].SetMinimum(min(ratios[i].GetMinimum() \
                             for i in range(len(ratios)))*0.99)
    ratios[0].SetStats(False)
    
    ratios[0].GetXaxis().SetTitle(x_label)
    ratios[0].GetXaxis().SetTitleSize(fontsize_title * sc2)
    ratios[0].GetXaxis().SetTitleOffset(1.1)
    ratios[0].GetXaxis().SetLabelOffset(0.02)
    ratios[0].GetXaxis().SetTickLength(Xtick_length * sc2)
    ratios[0].GetXaxis().SetLabelSize(fontsize_label * sc2)

    ratios[0].GetYaxis().SetTitle(y_label_lower)
    ratios[0].GetYaxis().SetTitleOffset(0.45)
    ratios[0].GetYaxis().SetTitleSize(fontsize_title * sc2)
    ratios[0].GetYaxis().SetTickLength(Ytick_length)
    ratios[0].GetYaxis().SetLabelSize(fontsize_label * sc2)
    ratios[0].GetYaxis().SetNdivisions(505)

    ratios[0].Draw("HIST")
    for i in range(1, len(ratios)):
        ratios[i].Draw("HIST SAME")

    const = ROOT.TLine(x_inf, 1, x_sup, 1)
    const.SetLineColor(colors_ratio[0])
    const.SetLineWidth(1)
    const.SetLineStyle(2)
    const.Draw("L SAME")

    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(leg, False)
    ROOT.SetOwnership(label, False)
    ROOT.SetOwnership(const, False)
    for i in range(len(histos)):
        ROOT.SetOwnership(histos[i], False)
    for i in range(len(ratios)):
        ROOT.SetOwnership(ratios[i], False)




# plot func confronto jet VBF

def plotNHisto1D_VBF(dfs, cols, nbins, x_inf, x_sup, x_label, 
                    y_label_upper, y_label_lower, leg_names, leg_names_fill, 
                    vbf_filter, colors, leg_ncol=1,
                    x1_leg=0.70, y1_leg=0.67, x2_leg=0.90, y2_leg=0.88,
                    x_labelCMS=0.17, y_labelCMS=0.82, x_labelLum=0.753, 
                    y_labelLum=0.95, sc1=1.43, sc2=3.3, 
                    filename="none", 
                    path="~/Documents/tesi/fig/plot/"):

    c = gen_canvas()

    leg = gen_legend(x1_leg, y1_leg, x2_leg, y2_leg, fontsize_legend * sc1)
    leg.SetNColumns(leg_ncol)

    histos = []
    for i in range(len(dfs)):
        histos.append(dfs[i].Histo1D(("", "", nbins, x_inf, x_sup), cols[i]))
        histos[i].SetLineColor(colors[i])
        histos[i].SetLineWidth(h_linewidth)
        if i == 0:
            leg.AddEntry(histos[i].GetValue(), leg_names[i], "f")
        else: 
            leg.AddEntry(histos[i].GetValue(), leg_names[i], "l")

    histos[0].SetMaximum(max(histos[i].GetMaximum()\
                              for i in range(len(histos)))*1.2)
    histos[0].SetMinimum(0)
    histos[0].SetStats(False)
    histos[0].GetXaxis().SetLabelSize(0)
    histos[0].GetXaxis().SetTitleSize(fontsize_title * sc1)
    histos[0].GetXaxis().SetTickLength(Xtick_length * sc1)
    histos[0].GetYaxis().SetLabelSize(fontsize_label * sc1)
    histos[0].GetYaxis().SetTickLength(Ytick_length)
    histos[0].GetYaxis().SetTitle(y_label_upper)
    histos[0].GetYaxis().SetTitleSize(fontsize_title * sc1)
    #histos[0].GetYaxis().SetTitleOffset(Ytitle_offset * sc1)
    histos[0].SetFillColor(colors[0])
    histos[0].SetFillStyle(3013) #3001

    filled_histos = []
    for i in range(1, len(dfs)):
        filled_histos.append(dfs[i].Filter(vbf_filter[i-1])\
                             .Histo1D(("", "", nbins, x_inf, x_sup), cols[i]))
        filled_histos[i-1].SetFillColor(colors[i])
        filled_histos[i-1].SetFillStyle(3013)
        filled_histos[i-1].SetLineWidth(0)
        #filled_histos[i-1].SetLineColor(ROOT.kBlack)
        leg.AddEntry(filled_histos[i-1].GetValue(), leg_names_fill[i-1], "f")

    pad1 = ROOT.TPad("pad1", "Upper pad", 0, 0.3, 1, 1)
    pad1.Draw()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(canvas_leftMargin)
    pad1.SetRightMargin(canvas_rightMargin)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.cd()

    histos[0].Draw("HIST")
    for i in range(1, len(histos)):
        histos[i].Draw("HIST SAME")
    for i in range(len(filled_histos)):
        filled_histos[i].Draw("HIST SAME")

    label = gen_CMSlabel(x_labelCMS, y_labelCMS, x_labelLum, y_labelLum,
                         fontsize_labelCMS * sc1, fontsize_labelLum * sc1)
    label.Draw()

    leg.Draw()
 
    c.cd()

    pad2 = ROOT.TPad("pad2", "Lower pad", 0, 0, 1, 0.3)
    pad2.Draw()
    pad2.SetTopMargin(0.015)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(canvas_leftMargin)
    pad2.SetRightMargin(canvas_rightMargin)
    pad2.cd()
    pad2.SetTickx(1)
    pad2.SetTicky(1)

    # più semplice utilizzando .Divide()
    ratios = []
    for i in range(1, len(histos)):
        ratios.append(ROOT.TH1D("", "", nbins, x_inf, x_sup))
        ratios[i-1].SetLineColor(colors[i])
        ratios[i-1].SetLineWidth(h_linewidth)
        for j in range(1, nbins+1):
            ratio = 0
            if histos[0].GetBinContent(j) > 0:
                ratio = histos[i].GetBinContent(j) / histos[0].GetBinContent(j)
            elif histos[0].GetBinContent(j) <= 0 \
                and histos[i].GetBinContent(j) > 0:
                ratio = 1
            ratios[i-1].SetBinContent(j, ratio)

    ratios[0].SetMaximum(max(ratios[i].GetMaximum() \
                             for i in range(len(ratios)))*1.01)
    ratios[0].SetMinimum(min(ratios[i].GetMinimum() \
                             for i in range(len(ratios)))*0.99)
    ratios[0].SetStats(False)
    
    ratios[0].GetXaxis().SetTitle(x_label)
    ratios[0].GetXaxis().SetTitleSize(fontsize_title * sc2)
    ratios[0].GetXaxis().SetTitleOffset(1.1)
    ratios[0].GetXaxis().SetLabelOffset(0.02)
    ratios[0].GetXaxis().SetTickLength(Xtick_length * sc2)
    ratios[0].GetXaxis().SetLabelSize(fontsize_label * sc2)

    ratios[0].GetYaxis().SetTitle(y_label_lower)
    ratios[0].GetYaxis().SetTitleOffset(0.45)
    ratios[0].GetYaxis().SetTitleSize(fontsize_title * sc2)
    ratios[0].GetYaxis().SetTickLength(Ytick_length)
    ratios[0].GetYaxis().SetLabelSize(fontsize_label * sc2)
    ratios[0].GetYaxis().SetNdivisions(505)

    ratios[0].Draw("HIST")
    for i in range(1, len(ratios)):
        ratios[i].Draw("HIST SAME")

    const = ROOT.TLine(x_inf, 1, x_sup, 1)
    const.SetLineColor(colors[0])
    const.SetLineWidth(2)
    const.SetLineStyle(2)
    const.Draw("L SAME")

    c.Draw()

    if filename != "none":
        c.SaveAs("".join([path, filename]))

    ROOT.SetOwnership(c, False)
    ROOT.SetOwnership(leg, False)
    ROOT.SetOwnership(label, False)
    ROOT.SetOwnership(const, False)
    for i in range(len(histos)):
        ROOT.SetOwnership(histos[i], False)
    for i in range(len(filled_histos)):
        ROOT.SetOwnership(filled_histos[i], False)
    for i in range(len(ratios)):
        ROOT.SetOwnership(ratios[i], False)
