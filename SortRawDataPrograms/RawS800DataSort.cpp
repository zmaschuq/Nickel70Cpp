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

    double ICSum;
    double CRDC1_yraw;
    double CRDC1_xraw;
    double CRDC2_yraw;
    double CRDC2_xraw;
    double E1UpScint;
    double E1DownScint;
    
    vector<float> ObjScint;
    vector<float> xfpScint;
};

int main(int argc, char** argv) {

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

    std::string RawRun = "140";
    std::string RawRunRootFile = "run_" + RawRun + "_Raw.root";
    TFile *RawOutput = TFile::Open(RawRunRootFile.c_str(), "RECREATE");
    TTree *rawTree = new TTree("RawData", "Raw Data");

    std::string S800FileName = "run-2140-00.root";
    std::string root_path = "/groups/tahn1/data/70Ni_NSCL/rootS800/cal/" + S800FileName;
    TFile *inputFile = TFile::Open(root_path.c_str(), "READ");
    TTree *tree = (TTree*)inputFile->Get("caltree");

    Long64_t entries = tree->GetEntries();
    S800Calc *s800cal = new S800Calc();
    TBranch *bS800cal = tree->GetBranch("s800calc");
    bS800cal->SetAddress(&s800cal);

    std::map<long long, S800Data> s800MapData;
    std::vector<long long> storeS800TS;
    storeS800TS.reserve(entries);

    for (Int_t i = 0; i < entries; i++) {
        s800cal->Clear();
        bS800cal->GetEntry(i);

        long long int S800TS = s800cal->GetTS();
        storeS800TS.push_back(S800TS);

        Double_t FirstObj = s800cal->GetMultiHitTOF()->GetFirstObjHit();
        Double_t FirstXf = s800cal->GetMultiHitTOF()->GetFirstXfHit();
        Double_t FirstE1Up = s800cal->GetMultiHitTOF()->GetFirstE1UpHit();
        Double_t FirstE1Down = s800cal->GetMultiHitTOF()->GetFirstE1DownHit();
        Double_t FirstRf = s800cal->GetMultiHitTOF()->GetFirstRfHit();


        std::vector<Float_t> fObj = s800cal->GetMultiHitTOF()->GetMTDCObj();
        std::vector<Float_t> fXf = s800cal->GetMultiHitTOF()->GetMTDCXf();

        // Obtaining raw CRDC values and Calibrating
        Double_t x0_fit = s800cal->GetCRDC(0)->GetXfit();
        Double_t x1_fit = s800cal->GetCRDC(1)->GetXfit();
        Double_t y0_raw = s800cal->GetCRDC(0)->GetTAC();
        Double_t y1_raw = s800cal->GetCRDC(1)->GetTAC();

        // Extracting Energy Loss
        Double_t EnergyLoss = s800cal->GetIC()->GetSum();

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

    std::vector<Double_t> raw_data; std::vector<Double_t> traceValues;
    std::vector<Double_t> Q;
    std::vector<Int_t> tb;
    std::vector<Double_t> x_list;
    std::vector<Double_t> y_list;
    std::vector<Double_t> r_list;
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
    Double_t IC;
    Double_t ftac1, ftac2;
    Double_t xfit1, xfit2;
    Double_t E1UpSignal;
    Double_t E1DownSignal;
    std::vector<Float_t> *ObjSignal = nullptr;
    std::vector<Float_t> *xfpSignal = nullptr;

    rawTree->Branch("S800TimeStamp", &S800TimeStamp, "S800TimeStamp/L");
    rawTree->Branch("EnergyLoss", &IC, "IC/D");
    rawTree->Branch("CRDC1_tac", &ftac1, "ftac1/D");
    rawTree->Branch("CRDC1_xfit", &xfit1, "xfit1/D");
    rawTree->Branch("CRDC2_tac", &ftac2, "ftac2/D");
    rawTree->Branch("CRDC2_xfit", &xfit2, "xfit2/D");
    rawTree->Branch("E1Up", &E1UpSignal, "E1UpSignal/D");
    rawTree->Branch("E1Down", &E1DownSignal, "E1DownSignal/D");
    rawTree->Branch("Object", &ObjSignal);
    rawTree->Branch("XFP", &xfpSignal);

    std::string h5filename = "run_0140.h5";
    std::string H5FilePath = "/groups/tahn1/data/70Ni_NSCL/h5/" + h5filename;
    std::sort(storeS800TS.begin(), storeS800TS.end());

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

    while (processed_events < num_events) {
        
        Q.clear(); tb.clear(); x_list.clear(); y_list.clear(); r_list.clear();
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

            auto it = std::lower_bound(storeS800TS.begin(), storeS800TS.end(), ATTPC_TimeStamp);

            long long minDelta = -1;

            if (it == storeS800TS.begin()) {
                minDelta = std::abs(ATTPC_TimeStamp - *it);
            } else if (it == storeS800TS.end()) {
                minDelta = std::abs(ATTPC_TimeStamp - *(it - 1));
            } else {
                // It's between *it and *(it-1), find which is closer
                minDelta = std::min(std::abs(ATTPC_TimeStamp - *it), std::abs(ATTPC_TimeStamp - *(it - 1)));
            }
                
            long long calcS800TS = ATTPC_TimeStamp - minDelta;
            if (s800MapData.find(calcS800TS) != s800MapData.end()) {

                H5::DataSet dataset = EventGroup.openDataSet(DataObjName);
                H5::DataSpace dataspace = dataset.getSpace();

                const int RANK = dataspace.getSimpleExtentNdims();
                hsize_t dims_out[RANK];
                dataspace.getSimpleExtentDims(dims_out, NULL);
                hsize_t n_rows = dims_out[0];
                hsize_t n_cols = dims_out[1];

                raw_data.resize(n_rows * n_cols);
                dataset.read(raw_data.data(), H5::PredType::NATIVE_DOUBLE);

                // Reshape into 2D vector
                std::vector<vector<double>> data(n_rows, vector<double>(n_cols));
                for (size_t i = 0; i < n_rows; ++i) {
                    for (size_t j = 0; j < n_cols; ++j) {
                        data[i][j] = raw_data[i * n_cols + j];
                    }
                }

                for (const auto& row : data) {
                    traceValues.clear();
                    
                    for (auto rowElement = row.begin() + 5; rowElement != row.end(); rowElement++) {
                        traceValues.push_back(*rowElement);
                    }

                    // Obtaining Mean Trace Values from First Twenty TBs
                    std::vector<double> firstTwentyTB(traceValues.begin(), traceValues.begin() + 20);
                    double mean_baseline = accumulate(firstTwentyTB.begin(), firstTwentyTB.end(), 0.0) / firstTwentyTB.size();

                    std::vector<double> traces(traceValues.begin() + 5, traceValues.begin() + 495);
                    if (traces.empty()) {continue;}

                    // Obtaining adc max and corresponding timebucket
                    double adc_max = *max_element(traces.begin(), traces.end());
                    double adc = adc_max - mean_baseline;
                    int tb_max = distance(traces.begin(), max_element(traces.begin(), traces.end())) + 5;

                    // const double drift_vel = 6.391e+6; // Units in  mm/s
                    // const double frequency = 3.125e+6; // Units in Hz
                    // double z_pos = drift_vel * tb_max / frequency;

                    if (adc > 110) {
                        
                        int cobo = static_cast<int>(row[0]);
                        int asad = static_cast<int>(row[1]);
                        int aget = static_cast<int>(row[2]);
                        int channel = static_cast<int>(row[3]);
                        int pad = static_cast<int>(row[4]);

                        x_list.push_back(x_pad_mapping[pad]);
                        y_list.push_back(y_pad_mapping[pad]);
                        r_list.push_back(r_pad_mapping[pad]);
                        Q.push_back(adc);
                        tb.push_back(tb_max);
                    }
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
            }
        } catch (...) {
            CurrentEvent++;
            continue;
        }

        rawTree->Fill();
        CurrentEvent++;
        processed_events++;
    }

    RawOutput->Write();
    RawOutput->Close();

    return 0;
}