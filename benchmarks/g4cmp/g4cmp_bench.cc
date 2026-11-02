// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 TC Yu and the CMPMC contributors
//
// G4CMP side of the CMPMC-vs-G4CMP speed comparison.
//
// Geometry, surfaces and sensor film are those of G4CMP's examples/phonon:
// a Ge cylinder (R = 3.81 cm, half-height 1.27 cm, [100] along z), Al Kaplan films
// on both faces with the example's film parameters, and the example's polynomial
// reflection model on all surfaces. One deliberate change: the sensor faces' simple
// (non-film) phonon absorption probability is 0 instead of 0.3, so that a sensor
// hit is processed identically in both codes (film entry with probability 0.2,
// otherwise reflection).
// Nothing is written to disk; the timed region is BeamOn() only.
//
// Usage: g4cmp_bench <mode> <events> [seed] [voltage_V]
//   mode = phonon   one 7.5 meV phonon per event, mode drawn from the Ge DOS
//          charge   one e-/h+ pair per event at the centre, uniform field,
//                   no Luke phonons produced
//          chargefull  as "charge", but every Luke phonon is produced and tracked
//   Charge modes produce no recombination phonons when a carrier is absorbed.

#include "G4Box.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPFieldManager.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPPhysicsList.hh"
#include "G4CMPStackingAction.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPUtils.hh"
#include "G4Event.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4RandomDirection.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Tubs.hh"
#include "G4UniformElectricField.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VProcess.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace {

constexpr G4double GHz = 1e9 * hertz;  // not defined by G4SystemOfUnits

std::string g_mode = "phonon";
double g_voltage = 4.0 * volt;

struct Counters {
    unsigned long long steps = 0, phonon_steps = 0, charge_steps = 0, tracks = 0, luke_phonons = 0;
    double edep_eV = 0.0;
} g_count;

class Detector : public G4VUserDetectorConstruction {
public:
    G4VPhysicalVolume* Construct() override {
        auto* nist = G4NistManager::Instance();
        auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
        auto* ge = nist->FindOrBuildMaterial("G4_Ge");
        auto* al = nist->FindOrBuildMaterial("G4_Al");

        auto* worldLV = new G4LogicalVolume(new G4Box("World", 16. * cm, 16. * cm, 16. * cm), vacuum, "World");
        auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0);

        auto* geLV = new G4LogicalVolume(new G4Tubs("Ge", 0., 3.81 * cm, 1.27 * cm, 0., 360. * deg), ge, "Ge");
        auto* gePV = new G4PVPlacement(nullptr, {}, geLV, "Ge", worldLV, false, 0);

        auto* LM = G4LatticeManager::GetLatticeManager();
        auto* lattice = new G4LatticePhysical(LM->LoadLattice(ge, "Ge"));
        lattice->SetMillerOrientation(1, 0, 0);
        LM->RegisterLattice(gePV, lattice);

        auto* alLV = new G4LogicalVolume(new G4Tubs("Al", 0., 3.81 * cm, 0.01 * cm, 0., 360. * deg), al, "Al");
        auto* topPV = new G4PVPlacement(nullptr, {0., 0., 1.28 * cm}, alLV, "AlTop", worldLV, false, 0);
        auto* botPV = new G4PVPlacement(nullptr, {0., 0., -1.28 * cm}, alLV, "AlBot", worldLV, false, 1);

        const std::vector<G4double> anh = {0, 0, 0, 0, 0, 1.51e-14};
        const std::vector<G4double> diff = {5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11};
        const std::vector<G4double> spec = {0.928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13};
        const G4double anhCutoff = 520., reflCutoff = 350.;

        auto sensor = [&](const char* name) {
            // pAbsProb = 0 (the stock example uses 0.3): a phonon reaching a sensor face
            // then enters the Kaplan film with probability filmAbsorption and otherwise
            // reflects, which is exactly CMPMC's sensor model with prob_abs = 1.
            auto* prop = new G4CMPSurfaceProperty(name, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
            prop->AddScatteringProperties(anhCutoff, reflCutoff, anh, diff, spec, GHz, GHz, GHz);
            auto* mpt = prop->GetPhononMaterialPropertiesTablePointer();
            G4CMP::UpdateMPT(mpt, "filmAbsorption", 0.20);
            G4CMP::UpdateMPT(mpt, "filmThickness", 600. * nm);
            G4CMP::UpdateMPT(mpt, "gapEnergy", 173.715e-6 * eV);
            G4CMP::UpdateMPT(mpt, "lowQPLimit", 3.);
            G4CMP::UpdateMPT(mpt, "phononLifetime", 242. * ps);
            G4CMP::UpdateMPT(mpt, "phononLifetimeSlope", 0.29);
            G4CMP::UpdateMPT(mpt, "vSound", 3.26 * km / s);
            G4CMP::UpdateMPT(mpt, "subgapAbsorption", 0.1);
            prop->SetPhononElectrode(new G4CMPPhononElectrode);
            return prop;
        };
        auto* wall = new G4CMPSurfaceProperty("Wall", 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        wall->AddScatteringProperties(anhCutoff, reflCutoff, anh, diff, spec, GHz, GHz, GHz);

        new G4CMPLogicalBorderSurface("detTop", gePV, topPV, sensor("Top"));
        new G4CMPLogicalBorderSurface("detBot", gePV, botPV, sensor("Bot"));
        new G4CMPLogicalBorderSurface("detWall", gePV, worldPV, wall);

        if (g_mode != "phonon") {
            auto* field = new G4UniformElectricField(G4ThreeVector(0., 0., -g_voltage / (2.54 * cm)));
            geLV->SetFieldManager(new G4CMPFieldManager(field), true);
        }
        return worldPV;
    }
};

class Primary : public G4VUserPrimaryGeneratorAction {
public:
    Primary() : gun_(1) { gun_.SetParticlePosition({0., 0., 0.}); }
    void GeneratePrimaries(G4Event* ev) override {
        if (g_mode == "phonon") {
            const G4double u = G4UniformRand();
            G4ParticleDefinition* pd = G4PhononLong::Definition();
            if (u < 0.53539) pd = G4PhononTransSlow::Definition();
            else if (u < 0.90217) pd = G4PhononTransFast::Definition();
            gun_.SetParticleDefinition(pd);
            gun_.SetParticleEnergy(7.5e-3 * eV);
            gun_.SetParticleMomentumDirection(G4RandomDirection());
            gun_.GeneratePrimaryVertex(ev);
        } else {
            gun_.SetParticleEnergy(1e-6 * eV);
            gun_.SetParticleMomentumDirection(G4RandomDirection());
            gun_.SetParticleDefinition(G4CMPDriftElectron::Definition());
            gun_.GeneratePrimaryVertex(ev);
            gun_.SetParticleMomentumDirection(G4RandomDirection());
            gun_.SetParticleDefinition(G4CMPDriftHole::Definition());
            gun_.GeneratePrimaryVertex(ev);
        }
    }
private:
    G4ParticleGun gun_;
};

class Stepping : public G4UserSteppingAction {
public:
    void UserSteppingAction(const G4Step* step) override {
        ++g_count.steps;
        const auto* pd = step->GetTrack()->GetParticleDefinition();
        if (G4CMP::IsPhonon(pd)) ++g_count.phonon_steps;
        else if (G4CMP::IsChargeCarrier(pd)) ++g_count.charge_steps;
        g_count.edep_eV += step->GetNonIonizingEnergyDeposit() / eV;  // as G4CMP::FillHit
    }
};

class Tracking : public G4UserTrackingAction {
public:
    void PreUserTrackingAction(const G4Track* track) override {
        ++g_count.tracks;
        const auto* creator = track->GetCreatorProcess();
        if (creator && creator->GetProcessName().find("Luke") != std::string::npos) ++g_count.luke_phonons;
    }
};

class Actions : public G4VUserActionInitialization {
public:
    void Build() const override {
        SetUserAction(new Primary);
        SetUserAction(new G4CMPStackingAction);  // attaches G4CMP track info (as in the examples)
        SetUserAction(new Stepping);
        SetUserAction(new Tracking);
    }
};

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s phonon|charge|chargefull <events> [seed] [voltage_V]\n", argv[0]);
        return 2;
    }
    g_mode = argv[1];
    const int events = std::atoi(argv[2]);
    const long seed = argc > 3 ? std::atol(argv[3]) : 12345;
    if (argc > 4) g_voltage = std::atof(argv[4]) * volt;

    CLHEP::HepRandom::setTheSeed(seed);
    const auto t0 = std::chrono::steady_clock::now();
    auto* rm = new G4RunManager;
    // Charge modes: no recombination phonons at charge absorption (CMPMC's default
    // ChargeRelease::RecordOnly emits none); "charge" also produces no Luke phonons.
    if (g_mode != "phonon") G4CMPConfigManager::SetGenPhonons(0.);
    if (g_mode == "charge") G4CMPConfigManager::SetLukeSampling(0.);
    rm->SetUserInitialization(new Detector);
    auto* physics = new G4CMPPhysicsList;
    physics->SetCuts();
    rm->SetUserInitialization(physics);
    rm->SetUserInitialization(new Actions);
    rm->SetVerboseLevel(0);
    rm->Initialize();
    const auto t1 = std::chrono::steady_clock::now();
    rm->BeamOn(events);
    const auto t2 = std::chrono::steady_clock::now();

    const double init_s = std::chrono::duration<double>(t1 - t0).count();
    const double run_s = std::chrono::duration<double>(t2 - t1).count();
    std::printf("code=G4CMP mode=%s events=%d init_s=%.3f run_s=%.3f us_per_event=%.1f "
                "steps=%llu phonon_steps=%llu charge_steps=%llu tracks=%llu edep_eV=%.6f luke_phonons=%llu\n",
                g_mode.c_str(), events, init_s, run_s, 1e6 * run_s / events, g_count.steps,
                g_count.phonon_steps, g_count.charge_steps, g_count.tracks, g_count.edep_eV,
                g_count.luke_phonons);
    delete rm;
    return 0;
}
