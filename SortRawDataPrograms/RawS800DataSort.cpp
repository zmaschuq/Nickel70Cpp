#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>
#include "H5Cpp.h"

#include "TFile.h"
#include "TH2D.h"
#include "TH2F.h"
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

#include "S800Calc.hh"
#include "S800Calibration.hh"
#include "S800.hh"

struct S800Data {

    float ICSum;
    float CRDC1_yraw;
    float CRDC1_xraw;
    float CRDC2_yraw;
    float CRDC2_xraw;
    float E1UpScint;
    float E1DownScint;
    
    vector<float> ObjScint;
    vector<float> xfpScint;
};

int main(int argc, char** argv) {

    std::string RawRun = argv[1];
    std::string S800FileName = argv[2];
    std::string h5filename = argv[3];

    // ======================== Reading in data to convert pad to coordinates ===============================
    const int num_of_pads = 10240;
    std::vector<float> x_pad_mapping(num_of_pads, 0.0f);
    std::vector<float> y_pad_mapping(num_of_pads, 0.0f);
    std::vector<float> r_pad_mapping(num_of_pads, 0.0f);
    int linesRead = 0;

    // Read Mapping File
    std::ifstream MappingFile("center_of_mass.txt");
    std::string line;

    while (getline(MappingFile, line)) {
        std::istringstream iss(line);
        int pad;
        float x, y;

        if (iss >> pad >> x >> y) {
            x_pad_mapping[pad] = x;
            y_pad_mapping[pad] = y; 
            r_pad_mapping[pad] = sqrt(x * x + y * y);
            linesRead++;
        }
        if (linesRead == num_of_pads) {break;}
    }

    //std::string RawRun = "139";
    std::string RawRunRootFile = "run_" + RawRun + "_Raw.root";
    TFile *RawOutput = TFile::Open(RawRunRootFile.c_str(), "RECREATE");
    TTree *rawTree = new TTree("RawData", "Raw Data");

    //std::string S800FileName = "run-2139-00.root";
    std::string root_path = "/groups/tahn1/data/70Ni_NSCL/rootS800/cal/" + S800FileName;
    TFile *inputFile = TFile::Open(root_path.c_str(), "READ");
    if (inputFile->IsZombie()) {
        cerr << "Root File Does Not Exist!" << endl;
        exit(1);
    }

    TTree *tree = (TTree*)inputFile->Get("caltree");

    Long64_t entries = tree->GetEntries();
    S800Calc *s800cal = new S800Calc();
    TBranch *bS800cal = tree->GetBranch("s800calc");
    bS800cal->SetAddress(&s800cal);

    std::map<long long, S800Data> s800MapData;

    for (Int_t i = 0; i < entries; i++) {
        s800cal->Clear();
        bS800cal->GetEntry(i);

        Long64_t S800TS = s800cal->GetTS();

        Float_t FirstObj = s800cal->GetMultiHitTOF()->GetFirstObjHit();
        Float_t FirstXf = s800cal->GetMultiHitTOF()->GetFirstXfHit();
        Float_t FirstE1Up = s800cal->GetMultiHitTOF()->GetFirstE1UpHit();
        Float_t FirstE1Down = s800cal->GetMultiHitTOF()->GetFirstE1DownHit();
        Float_t FirstRf = s800cal->GetMultiHitTOF()->GetFirstRfHit();


        std::vector<Float_t> fObj = s800cal->GetMultiHitTOF()->GetMTDCObj();
        std::vector<Float_t> fXf = s800cal->GetMultiHitTOF()->GetMTDCXf();

        // Obtaining raw CRDC values and Calibrating
        Float_t x0_fit = s800cal->GetCRDC(0)->GetXfit();
        Float_t x1_fit = s800cal->GetCRDC(1)->GetXfit();
        Float_t y0_raw = s800cal->GetCRDC(0)->GetTAC();
        Float_t y1_raw = s800cal->GetCRDC(1)->GetTAC();

        // Extracting Energy Loss
        Float_t EnergyLoss = s800cal->GetIC()->GetSum();

        S800Data data;

        data.ICSum = EnergyLoss;
        data.CRDC1_xraw = x0_fit;
        data.CRDC1_yraw = y0_raw;
        data.CRDC2_xraw = x1_fit;
        data.CRDC2_yraw = y1_raw;
        data.E1UpScint = FirstE1Up;
        data.E1DownScint = FirstE1Down;

        data.ObjScint = fObj;
        data.xfpScint = fXf;

        s800MapData[S800TS] = data;
    }

    std::vector<float> raw_data; std::vector<float> traceValues;
    std::vector<std::vector<double>> data;
    std::vector<float> firstTwentyTB; std::vector<float> traces;
    firstTwentyTB.reserve(20);

    std::vector<Float_t> Q;
    std::vector<Int_t> tb;
    std::vector<Float_t> x_list;
    std::vector<Float_t> y_list;
    std::vector<Float_t> r_list;
    Int_t Event;
    Long64_t ATTPC_TimeStamp;

     // Branches to hold ATTPC Data
    rawTree->Branch("Event", &Event, "Event/I");
    rawTree->Branch("ATTPCTS", &ATTPC_TimeStamp, "ATTPC_TimeStamp/L");
    rawTree->Branch("X", &x_list);
    rawTree->Branch("Y", &y_list);
    rawTree->Branch("TimeBucket", &tb);
    rawTree->Branch("Q", &Q);
    rawTree->Branch("R", &r_list);

    // Branches to Hold S800 Data
    Long64_t S800TimeStamp;
    Float_t IC;
    Float_t ftac1, ftac2;
    Float_t xfit1, xfit2;
    Float_t E1UpSignal;
    Float_t E1DownSignal;
    std::vector<Float_t> *ObjSignal = nullptr;
    std::vector<Float_t> *xfpSignal = nullptr;

    rawTree->Branch("S800TimeStamp", &S800TimeStamp, "S800TimeStamp/L");
    rawTree->Branch("EnergyLoss", &IC, "IC/F");
    rawTree->Branch("CRDC1_tac", &ftac1, "ftac1/F");
    rawTree->Branch("CRDC1_xfit", &xfit1, "xfit1/F");
    rawTree->Branch("CRDC2_tac", &ftac2, "ftac2/F");
    rawTree->Branch("CRDC2_xfit", &xfit2, "xfit2/F");
    rawTree->Branch("E1Up", &E1UpSignal, "E1UpSignal/F");
    rawTree->Branch("E1Down", &E1DownSignal, "E1DownSignal/F");
    rawTree->Branch("Object", &ObjSignal);
    rawTree->Branch("XFP", &xfpSignal);

    //std::string h5filename = "run_0139.h5";
    std::string H5FilePath = "/groups/tahn1/data/70Ni_NSCL/h5/" + h5filename;
    std::ifstream H5FileCheck(H5FilePath);
    if (!H5FileCheck) {
        cerr << "The HDF5 File Does Not Exist!" << endl;
        return 1;
    }
    H5FileCheck.close();

    H5::H5File file(H5FilePath.c_str(), H5F_ACC_RDONLY);
    H5::Group EventGroup = file.openGroup("/get");
    
    string FirstObjName = EventGroup.getObjnameByIdx(0);
    int evtPos = FirstObjName.find("evt");
    int UnderScorePos = FirstObjName.find("_");

    std::string EventNumberString = FirstObjName.substr(evtPos + 3, UnderScorePos - (evtPos + 3));
    int CurrentEvent = std::stoi(EventNumberString);

    int num_objs = EventGroup.getNumObjs();
    int num_events = num_objs / 2;
    int processed_events = 0;
    int noMatch = 0;
    int success = 0;

    while (processed_events < num_events) {
        
        Q.clear(); tb.clear(); x_list.clear(); y_list.clear(); r_list.clear();
        firstTwentyTB.clear(); traces.clear();
        Event = CurrentEvent;

        try {

            std::string DataObjName = "evt" + std::to_string(Event) + "_data";
            std::string HeaderObjName = "evt" + std::to_string(Event) + "_header";

            if (!EventGroup.nameExists(DataObjName) || !EventGroup.nameExists(HeaderObjName)) {
                CurrentEvent++;
                continue;
            }

            H5::DataSet HeaderData = EventGroup.openDataSet(HeaderObjName);
            std::vector<long long> hb(3);
            HeaderData.read(hb.data(), H5::PredType::NATIVE_LLONG);
            ATTPC_TimeStamp = hb[2];

            long long calcS800TS = 0;
            long long offsetted = ATTPC_TimeStamp - 1492;
            long long offsetted1 = ATTPC_TimeStamp - 1493;
            long long offsetted2 = ATTPC_TimeStamp + 281345310;
            long long offsetted3 = ATTPC_TimeStamp + 281345309;
            long long offsetted4 = ATTPC_TimeStamp + 77236264;
            long long offsetted5 = ATTPC_TimeStamp + 77236263;
            long long offsetted6 = ATTPC_TimeStamp + 105566589;
            long long offsetted7 = ATTPC_TimeStamp + 105566588;
            
            if (s800MapData.find(offsetted) != s800MapData.end()) {calcS800TS = offsetted;}
            else if (s800MapData.find(offsetted1) != s800MapData.end()) {calcS800TS = offsetted1;}
            else if (s800MapData.find(offsetted2) != s800MapData.end()) {calcS800TS = offsetted2;}
            else if (s800MapData.find(offsetted3) != s800MapData.end()) {calcS800TS = offsetted3;}
            else if (s800MapData.find(offsetted4) != s800MapData.end()) {calcS800TS = offsetted4;}
            else if (s800MapData.find(offsetted5) != s800MapData.end()) {calcS800TS = offsetted5;}
            else if (s800MapData.find(offsetted6) != s800MapData.end()) {calcS800TS = offsetted6;}
            else if (s800MapData.find(offsetted7) != s800MapData.end()) {calcS800TS = offsetted7;}

            if (calcS800TS == 0) {
                CurrentEvent++;
                processed_events++;
                noMatch++;
                continue;
            }

            H5::DataSet dataset = EventGroup.openDataSet(DataObjName);
            H5::DataSpace dataspace = dataset.getSpace();

            const int RANK = dataspace.getSimpleExtentNdims();
            hsize_t dims_out[RANK];
            dataspace.getSimpleExtentDims(dims_out, NULL);
            hsize_t n_rows = dims_out[0];
            hsize_t n_cols = dims_out[1];

            raw_data.resize(n_rows * n_cols);
            dataset.read(raw_data.data(), H5::PredType::NATIVE_UINT16);

            for (int i = 0; i < n_rows; i++) {
                traceValues.clear();
                
                auto row_start = raw_data.begin() + (i * n_cols);
                auto row_end = row_start + n_cols;
                traceValues.assign(row_start + 5, row_end);

                // Obtaining Mean Trace Values from First Twenty TBs
                firstTwentyTB.assign(traceValues.begin(), traceValues.begin() + 20);
                float mean_baseline = std::accumulate(firstTwentyTB.begin(), firstTwentyTB.end(), 0.0) / firstTwentyTB.size();

                traces.assign(traceValues.begin() + 5, traceValues.begin() + 495);
                if (traces.empty()) {continue;}

                // Obtaining adc max and corresponding timebucket
                float adc_max = *std::max_element(traces.begin(), traces.end());
                float adc = adc_max - mean_baseline;
                int tb_max = std::distance(traces.begin(), std::max_element(traces.begin(), traces.end())) + 5;

                // const double drift_vel = 6.391e+6; // Units in  mm/s
                // const double frequency = 3.125e+6; // Units in Hz
                // double z_pos = drift_vel * tb_max / frequency;
                        
                // int cobo = static_cast<int>(row[0]);
                // int asad = static_cast<int>(row[1]);
                // int aget = static_cast<int>(row[2]);
                // int channel = static_cast<int>(row[3]);
                int pad = static_cast<int>(*(row_start + 4));

                x_list.push_back(x_pad_mapping[pad]);
                y_list.push_back(y_pad_mapping[pad]);
                r_list.push_back(r_pad_mapping[pad]);
                Q.push_back(adc);
                tb.push_back(tb_max);
            
            }
            S800TimeStamp = calcS800TS;
            IC = s800MapData[calcS800TS].ICSum;
            ftac1 = s800MapData[calcS800TS].CRDC1_yraw;
            xfit1 = s800MapData[calcS800TS].CRDC1_xraw;
            ftac2 = s800MapData[calcS800TS].CRDC2_yraw;
            xfit2 = s800MapData[calcS800TS].CRDC2_xraw;
            E1UpSignal = s800MapData[calcS800TS].E1UpScint;
            E1DownSignal = s800MapData[calcS800TS].E1DownScint;
            ObjSignal = &(s800MapData[calcS800TS].ObjScint);
            xfpSignal = &(s800MapData[calcS800TS].xfpScint);
            

        } catch (...) {
            CurrentEvent++;
            processed_events++;
            continue;
        }

        rawTree->Fill();
        CurrentEvent++;
        processed_events++;
        success++;
    }

    cout << "Total # of Events in Run: " << num_events << endl;
    cout << "Successes: " << success << " No Matches: " << noMatch << endl;
    RawOutput->Write();
    RawOutput->Close();

    return 0;
}