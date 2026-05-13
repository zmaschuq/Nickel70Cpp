#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>

#include "TFile.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TCutG.h"
#include "TApplication.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TBranch.h"

using namespace std;

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);

    string RawRootFile = "run_20_Raw.root";
    TFile *inputFile = TFile::Open(RawRootFile.c_str(), "READ");
    if (inputFile->IsZombie()) {
        cerr << "This Root File Doesn't Exist!" << endl;
        exit(1);
    } 

    TH2D *Beam = new TH2D("Beam", "Beam Profile; Obj Timing;XFP - Obj Timing", 60, -100, -25, 100, 200, 300);

    vector<float> *ObjVector = nullptr; vector<float> *XFPVector = nullptr;
    float ObjValue; float XFPValue; float dE;

    TTree *tree = (TTree*)inputFile->Get("RawData");
    tree->SetBranchAddress("Object", &ObjVector);
    tree->SetBranchAddress("XFP", &XFPVector);
    tree->SetBranchAddress("EnergyLoss", &dE);

    int NumberEntries = tree->GetEntries();
    for (int i = 0; i < NumberEntries; i++) {
        tree->GetEntry(i);
        ObjValue = 0; XFPValue = 0;

        for (int j = 0; j < ObjVector->size(); j++) {
            if (ObjVector->at(j) > -100 && ObjVector->at(j) < -30) {ObjValue = ObjVector->at(j);}
        }
        for (int j = 0; j < XFPVector->size(); j++) {
            if (XFPVector->at(j) > 100 && XFPVector->at(j) < 300) {XFPValue = XFPVector->at(j);}
        }

        float XFP_Obj = XFPValue - ObjValue;
        Beam->Fill(ObjValue, XFP_Obj);

    }

    TCanvas *c1 = new TCanvas();
    h2->Draw("colz");

    app.Run();
    return 0;
}