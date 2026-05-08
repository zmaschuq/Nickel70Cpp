#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <sstream>
#include <unordered_set>
#include "H5Cpp.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TApplication.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"

using namespace std;
#define PI 3.14159

double *extrapolate(vector<double> &, double, const double);
vector<double> g(vector<float> &, vector<double> &, float, float &, float &);
double BasisSplines(int, int, double, const vector<double> &);
double EvaluateSpline(double, const vector<double> &, const vector<double> &, int);
vector<double> ComputingCoeff(vector<double> &, vector<double> &, vector<double> &, int);
void ConvolvedBraggModel(vector<double> &, vector<double> &, vector<double> &, 
const vector<double> &, int, int, double);

int main() {
    
    int argc = 0;
    char *argv[] = {};
    TApplication theApp("App", &argc, argv);
    
    // ======================== Reading in data to convert pad to coordinates ===============================
    const int num_of_pads = 10240;
    vector<float> x_pad_mapping(num_of_pads, 0.0f);
    vector<float> y_pad_mapping(num_of_pads, 0.0f);
    vector<float> r_pad_mapping(num_of_pads, 0.0f);
    int linesRead = 0;

    // Read Mapping File
    ifstream MappingFile("center_of_mass.txt");
    string line;

    while (getline(MappingFile, line)) {
        istringstream iss(line);
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

   // ======================================== Inputing HDF5 File =============================================

    string RunFileName;
    cout << "Input run file name:" << endl;
    cin >> RunFileName; //HDF5 File Name: run_00XX.h5

    string H5FilePath;
    H5FilePath = "/groups/tahn1/data/70Ni_NSCL/h5/" + RunFileName;

    ifstream myfile(H5FilePath);
    if (!myfile) {
        cerr << "This File Does Not Exist" << endl;
        return 1;
    }
    myfile.close();
    
    // ===================== Defining Output File and Reading & Reading In S800 PID Events =========================
    
    ofstream Ni70_Results("70Ni_Results.txt", ios::trunc);
    Ni70_Results << "Event" << "  " << "Angle" << "  " << "Er" << "  " << "Vertex" << "  " << "  " << "Ex" <<
    "  " << "ChiBragg" << "  " << "Intx" << "  " << "Scale" << "  " << "LSQChi" << endl;

    ifstream S800EventResults("TrackEventsPID.txt");
    int PIDEvents;
    vector<int> EventsPID;

    while (S800EventResults >> PIDEvents) {
        EventsPID.push_back(PIDEvents);
    }

    // ================================= Reading In Energy Loss File (SRIM) ======================================

    ifstream SRIMEnergyLossFile("energy_loss_He_11MeV.txt");
    double energyLoss, distanceTravel;

    vector<double> dEdx;
    vector<double> alphaDistanceTravel;
    vector<double> reverseEnergyLoss;
    vector<double> reverseAlphaDistance;

    while (SRIMEnergyLossFile >> distanceTravel >> energyLoss) {

        dEdx.push_back(energyLoss);
        alphaDistanceTravel.push_back(distanceTravel);
    }

    vector<double>::reverse_iterator rit;
    vector<double>::reverse_iterator rit1;

    // SRIM calculates energy loss FROM end of track, so want to reverse
    for (rit = dEdx.rbegin(); rit != dEdx.rend(); ++rit) {reverseEnergyLoss.push_back(*rit);}
    for (rit1 = alphaDistanceTravel.rbegin(); rit1 != alphaDistanceTravel.rend(); ++rit1) {
        reverseAlphaDistance.push_back(*rit1);
    }

    // Extrapolating the Bragg Curve with Linear Model
    double start = reverseAlphaDistance.back() + 1.0;
    double end   = reverseAlphaDistance.back() + 200.0;
    int nPoints  = 150;

    vector<double> extendAlphaDistance;
    extendAlphaDistance.reserve(nPoints);

    for (int i = 0; i < nPoints; i++) {
        double val = start + (end - start) * i / (nPoints - 1);
        val = std::round(val * 100.0) / 100.0; // round to 2 decimals
        extendAlphaDistance.push_back(val);
    }

    // Concatenate to original vector
    reverseAlphaDistance.insert(reverseAlphaDistance.end(),extendAlphaDistance.begin(),extendAlphaDistance.end());

    // Calculating new Energy Loss from Extrapolation
    double EPSlope = -6.6713232e-5; double EPIntrcpt = 0.248613273; // Calculated using Python
    double *dEdx_Extension = extrapolate(extendAlphaDistance, EPSlope, EPIntrcpt);
    
    for (int i = 0; i < 150; i++) {
        reverseEnergyLoss.push_back(dEdx_Extension[i]);
    }
    delete [] dEdx_Extension;
    dEdx_Extension = nullptr;
    
    // ================================== Spling the Original, Extrpolated Bragg Curve ==================================
    
    // Splining the Original Bragg Curve
    int splineDeg = 3;
    vector<double> AlphaX(reverseAlphaDistance);
    vector<double> s_knots = {AlphaX[0], AlphaX[0], AlphaX[0], AlphaX[0], AlphaX[2], AlphaX[5], AlphaX[7], 
        AlphaX[10], AlphaX[13], AlphaX[16], AlphaX[23], AlphaX[28], AlphaX[33], AlphaX[36], AlphaX[50], 
        AlphaX[85], AlphaX[96], AlphaX[104], AlphaX[200], AlphaX.back(), AlphaX.back(), AlphaX.back(), AlphaX.back()};

    vector<double> splineCoeff = {1.72986300e-03,  1.79761991e-03,  2.02189255e-03,  2.35494371e-03, 
    2.97410001e-03,  3.47673751e-03,  4.27879205e-03,  5.15550036e-03, 
    6.40442249e-03,  8.30772500e-03,  9.06274315e-03,  7.81611792e-03, 
    4.52578335e-03,  1.53405384e-03, -2.41494076e-05,  2.15095570e-05, 
    -1.58384820e-05,  5.72638308e-06, -2.60945180e-06};

    auto splineFunc = [&](double *x, double *) {
        return EvaluateSpline(x[0], s_knots, splineCoeff, splineDeg);
    };

    double xmin = AlphaX.front(); double xmax = AlphaX.back();
    TF1 *fSpline = new TF1("Spline", splineFunc, xmin, xmax, 0);

    // ================================== Convolving and Splining the Bragg Curve ========================================

    // Discretizing Bragg Curve for Convolution
    vector<double> xDiscretization;
    for (double x = 0.0; x <= *std::max_element(AlphaX.begin(), AlphaX.end()); x += 5.0) {xDiscretization.push_back(x);}
        
    vector<double> energyLossDiscret;
    for (double x : xDiscretization) {energyLossDiscret.push_back(fSpline->Eval(x));}
    delete fSpline;

    // Convolving the Bragg Curve
    vector<double> BraggCenters;
    vector<double> ConvBraggSum;
    vector<float> window_range;
    vector<float> Bragg_tau_range;
    float idxConv = 0.0; float idxConv1 = 0.0; float step_size = 5.0;

    while (idxConv < *std::max_element(xDiscretization.begin(), xDiscretization.end()) + 100.0) {
        window_range.push_back(idxConv);
        idxConv += 0.1;
    }

    while (idxConv1 < *std::max_element(window_range.begin(), window_range.end()) - 20.0) {
        Bragg_tau_range.push_back(idxConv1);
        idxConv1 += step_size;
    }

    for (float tau : Bragg_tau_range) {

        float BraggStart, BraggEnd;
        vector<double> Bragg_g_vector = g(window_range, energyLossDiscret, tau, BraggStart, BraggEnd);
        BraggCenters.push_back((BraggStart + BraggEnd) / 2.0);
        
        float BraggSum = 0;
        for (int i = 0; i < xDiscretization.size(); i++) {
            if (xDiscretization.at(i) >= BraggStart && xDiscretization.at(i) <= BraggEnd) {
                BraggSum += energyLossDiscret.at(i);
            }
        }
        ConvBraggSum.push_back(BraggSum);
    }

    // Splining Convolved Bragg Curve
    vector<double> s1_knots = {BraggCenters[0], BraggCenters[0], BraggCenters[0], BraggCenters[0], BraggCenters[175], 
    BraggCenters[450], BraggCenters[580], BraggCenters[660], BraggCenters[711], BraggCenters[717], BraggCenters[723], 
    BraggCenters[728], BraggCenters[735], BraggCenters[742], BraggCenters[748], BraggCenters[793], BraggCenters.back(), 
    BraggCenters.back(), BraggCenters.back(), BraggCenters.back()};

    vector<double> CoeffConvolvedBragg = ComputingCoeff(BraggCenters, ConvBraggSum, s1_knots, splineDeg);

    // ================================ Reading in Recoil Energy Data (SRIM) =====================================

    // Obtaining the Recoil Energy of Alphas
    ifstream EvXFile("EvsX_11MeV.txt");
    double EvX_Energy, EvX_Position;
    vector<double> x_R; vector<double> energy_R;
    vector<double> xR_rev; vector<double> energyR_rev;

    while (EvXFile >> EvX_Position >> EvX_Energy) {
        x_R.push_back(EvX_Position);
        energy_R.push_back(EvX_Energy);
    }

    vector<double>::reverse_iterator xRit;
    vector<double>::reverse_iterator eRit;

    for (xRit = x_R.rbegin(); xRit != x_R.rend(); ++xRit) {xR_rev.push_back(*xRit);}
    for (eRit = energy_R.rbegin(); eRit != energy_R.rend(); ++eRit) {energyR_rev.push_back(*eRit);}

    vector<double> EvsXKnots = {xR_rev[0], xR_rev[0], xR_rev[0], xR_rev[0], xR_rev[2], xR_rev[5], xR_rev[7], 
    xR_rev[10], xR_rev[13], xR_rev[16], xR_rev[23], xR_rev[28], xR_rev[33], xR_rev[36], xR_rev[50], xR_rev.back(), 
    xR_rev.back(), xR_rev.back(), xR_rev.back()};

    vector<double> EvXcoeff = ComputingCoeff(xR_rev, energyR_rev, EvsXKnots, splineDeg);
    
    auto EvXsplineFunc = [&](double *x, double *) {
        return EvaluateSpline(x[0], EvsXKnots, EvXcoeff, splineDeg);
    };
    
    // =============================== Declaring Histogram for Result Visualization ===================================

    //TH2D *hEhT = new TH2D("hEhT", "Recoil Energy vs Angle", 181, 0, 180, 181, 0, 3);
    TH1D *h1 = new TH1D("h1", "Excitation Spectrum", 50, 0, 50);
    // TH1D *hChi = new TH1D("hChi", "Chi-Squared (Normalized)", 2101, -100, 4000);

    // ============================ Declaring Constants, Vectors to Store Data =========================================

    // Vectors to store raw data
    vector<double> raw_data; vector<double> traceValues;
    vector<double> Q;
    vector<double> z;
    vector<double> x_list;
    vector<double> y_list;
    vector<double> r_list;

    // Hough Transform Constants
    const float theta_high = 90.0; const float theta_low = -90.0; const float theta_increment = 0.5;
    const float r_increment = 2.0;
    const float theta_diff = theta_high - theta_low;
    const int thetabins = round(theta_diff/theta_increment) + 1.0;

    const float AT_TPC_Radius = 275.0;
    const float AT_TPC_Length = 1400.0;
    const float rhoBounds = sqrt(AT_TPC_Radius*AT_TPC_Radius + AT_TPC_Length*AT_TPC_Length);
    const float rhoMin = -rhoBounds, rhoMax = rhoBounds, rhoRange = rhoMax - rhoMin;
    const int rbins = round(rhoRange/r_increment) + 1.0;
        
    const int threshold = 10;
    vector<vector<int>> accumulator(rbins, vector<int> (thetabins));
    vector<double> slopesHT; vector<double> interceptHT; vector<int> intersectionsHT;

    // Vectors to store data from Hough Transform
    vector<double> isolated_r;
    vector<double> isolated_Q;
    vector<double> isolated_z;
    vector<double> isolated_x;
    vector<double> isolated_y;

    // Vectors storing convolved ADC data
    vector<double> xDataConv;
    vector<double> yDataConv;
    vector<float> step_func_range;
    vector<float> tau_range;

    vector<int> shiftRange; vector<int> scaleYRange; vector<double> QError; vector<double> yModel;
    for (int i = 3400; i < 3707; i++) {shiftRange.push_back(i);}
    for (int i = 400000; i < 600000; i += 10000) {scaleYRange.push_back(i);}
    vector<vector<double>> chiSquaredMatrix(shiftRange.size(), vector<double>(scaleYRange.size()));

    // ======================================== Event Loop =================================================================

    // Opening HDF5 File
    H5::H5File file(H5FilePath, H5F_ACC_RDONLY);

    for (int i = 0; i < EventsPID.size(); i++) {

        Q.clear(); z.clear(); x_list.clear(); y_list.clear(); r_list.clear();
        slopesHT.clear(); interceptHT.clear(); intersectionsHT.clear();
        isolated_r.clear(); isolated_Q.clear(); isolated_z.clear(); isolated_x.clear(); isolated_y.clear();
        xDataConv.clear(); yDataConv.clear(); step_func_range.clear(); tau_range.clear();
        
        int event = EventsPID[i]; 
        string datasetPath = "get/evt" + to_string(event) + "_data";

        // Check to see if event number exists
        if (!file.nameExists(datasetPath)) {
            cerr << "WARNING: Dataset not found for event " << event << ". Skipping." << endl;
            continue; 
        }

        H5::DataSet dataset = file.openDataSet(datasetPath);
        H5::DataSpace dataspace = dataset.getSpace();

        const int RANK = dataspace.getSimpleExtentNdims();
        hsize_t dims_out[RANK];
        dataspace.getSimpleExtentDims(dims_out, NULL);
        hsize_t n_rows = dims_out[0];
        hsize_t n_cols = dims_out[1];

        raw_data.resize(n_rows * n_cols);
        dataset.read(raw_data.data(), H5::PredType::NATIVE_DOUBLE);

        // Reshape into 2D vector
        vector<vector<double>> data(n_rows, vector<double>(n_cols));
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
            vector<double> firstTwentyTB(traceValues.begin(), traceValues.begin() + 20);
            double mean_baseline = accumulate(firstTwentyTB.begin(), firstTwentyTB.end(), 0.0) / firstTwentyTB.size();

            vector<double> traces(traceValues.begin() + 5, traceValues.begin() + 495);
            if (traces.empty()) {continue;}

            // Obtaining adc max and corresponding timebucket
            double adc_max = *max_element(traces.begin(), traces.end());
            double adc = adc_max - mean_baseline;
            int tb_max = distance(traces.begin(), max_element(traces.begin(), traces.end())) + 5;

            const double drift_vel = 6.391e+6; // Units in  mm/s
            const double frequency = 3.125e+6; // Units in Hz
            double z_pos = drift_vel * tb_max / frequency;

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
                z.push_back(z_pos);
            }
        }

        // Some events are empty, skip event
        if (r_list.empty()) {
            Ni70_Results << event << " " << "No Data!" << endl;
            continue;
        }

        // Here must use Hough Transform to isolate particle tracks and subsequently analyze them
        for (int i = 0; i < rbins; i++) {
            for (int j = 0; j < thetabins; j++) {
                accumulator[i][j] = 0;
            }
        }

        for (int i = 0; i < r_list.size(); i++) {
            for (float j = theta_low; j <= theta_high; j+=theta_increment) {

                float radians = j * PI/180.0;
                float rho = z[i]*cos(radians) + r_list[i]*sin(radians);
                
                int rIndex = round((rho - rhoMin)/r_increment);
                int thetaIndex = round((j - theta_low)/theta_increment);
                accumulator[rIndex][thetaIndex]++;
            }
        }

        // Extracting the theta and r values that correspond to the highest voted cell
        for (int i = 0; i < rbins; i++) {
            for (int j = 0; j < thetabins; j++) {

                if (accumulator[i][j] >= threshold) {

                    float r_peak = rhoMin + i*r_increment;
                    float theta_peak = (theta_low + j*theta_increment) * PI / 180.0;
                    float slope = -cos(theta_peak) / sin(theta_peak);
                    float intercept = r_peak / sin(theta_peak);
                    float roundSlope = std::round(slope * 100.0) / 100.0;
                    float roundIntercept = std::round(intercept * 100.0) / 100.0;

                    //if (roundSlope > 0.0) {
                    slopesHT.push_back(roundSlope);
                    interceptHT.push_back(roundIntercept);
                    intersectionsHT.push_back(accumulator[i][j]);
                    //} 
                } 
            }
        }

        // Extracting identified track from Hough Transform slope and intercept
        int max_intersections = -1;
        double best_slope = 0.0; double best_intercept = 0.0;

        for (int i = 0; i < slopesHT.size(); i++) {
            if (intersectionsHT[i] > max_intersections) {

                max_intersections = intersectionsHT[i];
                best_slope = slopesHT[i];
                best_intercept = interceptHT[i];
            }
        }
        
        // Defining region perpendicular to identified HT line
        double perpDist = 7.5;
        double slope = best_slope;
        double intercept = best_intercept;
        double distDenom = sqrt(1 + slope*slope);
        bool InfiniteSlope = false;

        for (int i = 0; i < z.size(); i++) {
            double r_i = r_list[i];
            double z_i = z[i];

            double dist;
            if (isinf(slope)) {
                dist = fabs(z_i - intercept);  // vertical line
                InfiniteSlope = true;
            } else {
                dist = fabs(r_i - slope * z_i - intercept) / distDenom;
            }

            if (dist < perpDist) {
                isolated_r.push_back(r_i);
                isolated_Q.push_back(Q[i]);
                isolated_z.push_back(z_i);
                isolated_x.push_back(x_list[i]);
                isolated_y.push_back(y_list[i]);
            }
        }

        // Condition that does not meet definition of line (at least 10 points)
        if (isolated_r.empty()) {
            Ni70_Results << event << " " << "No Track!" << endl;
            continue;
        }

        // Performin LSQ Fit to Extract Scattering Angle
        double weight;
        double weightLow = 2.5; double weightHigh = 5.0;
        double A, B, C, D, E, F;
        A = 0.0; B = 0.0; C = 0.0; D = 0.0; E = 0.0; F = 0.0;

        for (int i = 0; i < isolated_r.size(); i++) {
            weight = (isolated_r[i] < 150) ? weightLow : weightHigh;

            A += isolated_z[i] / (weight*weight);
            B += 1 / (weight*weight);
            C += isolated_r[i] / (weight*weight);
            D += (isolated_z[i]*isolated_z[i]) / (weight*weight);
            E += (isolated_z[i]*isolated_r[i]) / (weight*weight);
            F += (isolated_r[i]*isolated_r[i]) / (weight*weight);
        }

        double track_slope = (E*B - C*A) / (D*B - A*A);
        double track_intercept = (D*C - E*A) / (D*B - A*A);
        double slope_err = B / (B*D - A*A);
        double intercept_err = D / (B*D - A*A);
        double covar = -A / (B*D - A*A);
        
        double interceptZ = -track_intercept / track_slope;
        double LabAngle;
        double uncertaintyAngle = 1 / (1.0 + track_slope*track_slope) * sqrt(slope_err) * 180.0/PI;

        LabAngle = (InfiniteSlope == true) ? 90.0 : atan2(track_slope, 1.0) * 180.0/PI;

        // Conditions to exclude events where the vertex is not within active volume and angles too high
        if (interceptZ < 0.0 || interceptZ > 1400.0) {
            Ni70_Results << event << " " << "No Track!" << endl;
            continue;
        }

        if (LabAngle > 90.0 || LabAngle < 0) {
            Ni70_Results << event << " " << "No Track!" << endl;
            continue;
        }

        double LSQChiVal = 0.0; double weightLSQ;
        for (int i = 0; i < isolated_r.size(); i++) {
            weightLSQ = (isolated_r[i] < 150) ? weightLow : weightHigh;
            
            double modelVal = track_slope * isolated_z[i] + track_intercept;
            double LSQDiff = isolated_r[i] - modelVal;
            LSQChiVal += (LSQDiff*LSQDiff) / (weightLSQ*weightLSQ);
        }
        double reducedLSQChi = LSQChiVal / (isolated_r.size() - 2);
    
        // Convolving the ADC Data
        float idx = 0.0; float idx1 = 21.75;
        while (idx < *std::max_element(isolated_r.begin(), isolated_r.end()) + 100.0) {
            step_func_range.push_back(idx);
            idx += 0.1;
        }

        while (idx1 < *std::max_element(step_func_range.begin(), step_func_range.end()) - 20.0) {
            tau_range.push_back(idx1);
            idx1 += step_size;
        }

        for (float tau : tau_range) {

            float g_start, g_end;
            vector<double> g_vector = g(step_func_range, isolated_Q, tau, g_start, g_end);
            xDataConv.push_back(round((g_start + g_end) * 1000.0 / 2.0) / 1000.0);
            
            float sum_charge = 0;
            for (int i = 0; i < isolated_r.size(); i++) {
                if (isolated_r[i] >= g_start && isolated_r[i] <= g_end) {
                    sum_charge += isolated_Q[i];
                }
            }
            yDataConv.push_back(sum_charge);
        }

        float scaleFactor = sin(LabAngle*PI / 180.0);
        // Performing Bragg Curve Fit on Convolved Data
        QError.assign(yDataConv.size(), 200.0);
        
        // Must reset all elements in ChiSquared Matrix per event
        for (auto &row : chiSquaredMatrix) {
            std::fill(row.begin(), row.end(), NAN);
        }

        double QThreshold = 500.0;
        for (int idxShift = 0; idxShift < shiftRange.size(); idxShift++) {
            
            ConvolvedBraggModel(yModel, xDataConv, s1_knots, CoeffConvolvedBragg, 
                splineDeg, shiftRange[idxShift], scaleFactor);

            for (int idxScale = 0; idxScale < scaleYRange.size(); idxScale++) {

                double chi2 = 0.0;
                for (int k = 0; k < yDataConv.size(); k++) {
                    if (yDataConv[k] > QThreshold) {
                        double diff = yDataConv[k] - yModel[k]*scaleYRange[idxScale];
                        chi2 += diff*diff / (QError[k]*QError[k]);
                    }
                }
                chiSquaredMatrix[idxShift][idxScale] = chi2;
            }
        }

        // Obtaining best shift and scale value from minimum chi-squared
        double maxChi = std::numeric_limits<double>::infinity();
        int shiftIdx = -1; int scaleIdx = -1;

        for (int i = 0; i < chiSquaredMatrix.size(); i++) {
            for (int j = 0; j < chiSquaredMatrix[i].size(); j++) {

                double chiValue = chiSquaredMatrix[i][j];
                if (!std::isnan(chiValue) && chiValue < maxChi) {

                    maxChi = chiValue;
                    shiftIdx = i;
                    scaleIdx = j;
                }
            }
        }

        int bestShift = shiftRange[shiftIdx]; int bestScale = scaleYRange[scaleIdx];

        double maxChi_norm = maxChi / yDataConv.size();
        if (maxChi_norm > 5.0 || maxChi_norm < 0.30) {
            Ni70_Results << event << " " << "Bad Chi Squared!" << endl;
            continue;
        }

        // hChi->Fill(maxChi_norm);

        // cout << "Chi Value: " << maxChi << endl;
        // cout << "Shift: " << bestShift << " Scale: " << bestScale << endl;

        // Obtain Recoil Energy and Store in File
        double xRmin = xR_rev.front(); double xRmax = xR_rev.back();
        TF1 *EvXSpline = new TF1("Spline", EvXsplineFunc, xRmin, xRmax, 0);
        float EnergyofAlpha = EvXSpline->Eval(bestShift);
        
        // Computing Excitation Energy of Nickel
        const float constant = 931.5;
        const int T2 = 4760;
        double T3 = EnergyofAlpha;
        double m1, m2, m3, m4;
        double ExcitationEnergy;

        m1 = 4.0022602 * constant; m3 = m1;
        m2 = 69.93614 * constant;

        double FirstTerm = m1 + m2 + T2 - m3 - T3;
        double SecondTerm = sqrt(T2*T2 + 2*m2*T2) * sqrt(T3*T3 + 2*m3*T3);

        m4 = (FirstTerm*FirstTerm - (T3*T3 + T2*T2) - 2*(m2*T2 + m3*T3) + 2*SecondTerm*cos(LabAngle * PI/180.0));
        ExcitationEnergy = sqrt(m4) - m2;

        Ni70_Results << event << "  " << LabAngle << "  " << EnergyofAlpha << "  " << interceptZ << "  " << ExcitationEnergy <<
        "  " << "  " << maxChi_norm << "  " << max_intersections << "  " << bestScale << 
        "  " << reducedLSQChi << endl;

        delete EvXSpline;

        h1->Fill(ExcitationEnergy);
        //hEhT->Fill(LabAngle, EnergyofAlpha);
    }

    h1->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    h1->GetYaxis()->SetTitle("Counts");
    TCanvas *c11 = new TCanvas();
    h1->Draw();

    // hChi->GetXaxis()->SetTitle("Chi Squared");
    // hChi->GetYaxis()->SetTitle("Counts");
    // TCanvas *c12 = new TCanvas();
    // hChi->Draw();

    // hEhT->GetXaxis()->SetTitle("Lab Angle");
    // hEhT->GetYaxis()->SetTitle("Recoil Energy");
    // TCanvas *c10 = new TCanvas();
    // hEhT->Draw();

    theApp.Run();
    return 0;
}

double *extrapolate(vector<double> &x, double a, const double b) {
    
    double result;
    double *reverseAlphaExtension = new double [x.size()];

    for (int i = 0; i < x.size(); i++) {
        result = x.at(i)*a + b;
        result = std::max(result, 0.0); // y-values that are negative, keep at 0
        reverseAlphaExtension[i] = result;
    }

    return reverseAlphaExtension;
}

vector<double> g(vector<float> &t, vector<double> &Q, float tau, float &start_step, float &end_step) {

    int size = t.size();
    vector<double> convolutionVector (size, 0.0);

    float step_func_width = 20.0;
    start_step = tau;
    end_step = start_step + step_func_width;

    double maxQ = *max_element(Q.begin(), Q.end());
    for (int i = 0; i < size; i++) {
        if (t[i] >= start_step && t[i] <= end_step) {
            convolutionVector.at(i) = maxQ + 50.0;
        }
    }
    return convolutionVector;
}

double BasisSplines(int i, int deg, double BraggPos, const vector<double> &knots) {

    // Implementing Cox de Boor recursion
    if (deg == 0) {
        if (i == knots.size() - 2 && BraggPos == knots.back()) return 1.0;
        return (BraggPos >= knots[i] && BraggPos < knots[i+1]) ? 1.0 : 0.0;
    }

    double denom1 = knots[i + deg] - knots[i];
    double denom2 = knots[i + deg + 1] - knots[i + 1];
    double inv1 = (denom1 != 0.0) ? 1.0 / denom1 : 0.0;
    double inv2 = (denom2 != 0.0) ? 1.0 / denom2 : 0.0;
    
    double numerator1 = BraggPos - knots[i];
    double numerator2 = knots[i + deg + 1] - BraggPos;

    return numerator1 * inv1 * BasisSplines(i, deg - 1, BraggPos, knots) + 
    numerator2 * inv2 * BasisSplines(i + 1, deg - 1, BraggPos, knots);
}

double EvaluateSpline(double BraggPos, const vector<double> &knots, const vector<double> &coeff, int deg) {

    double sum = 0.0;
    BraggPos = std::max(knots.front(), std::min(BraggPos, knots.back()));

    int nBasis = knots.size() - deg - 1;
    for (int i = 0; i < nBasis; i++) {sum += coeff.at(i) * BasisSplines(i, deg, BraggPos, knots);}

    return sum;
}

vector<double> ComputingCoeff(vector<double> &xDisc, vector<double> &ELossDisc, vector<double> &knots, int deg) {
    
    // Constructing design matrix to compute coefficients (c): Ac = y
    int NumBasis = knots.size() - deg - 1;
    int NumData = xDisc.size();

    TMatrixD DesignMtx(NumData, NumBasis);
    for (int i = 0; i < NumData; i++) {
        for (int j = 0; j < NumBasis; j++) {

            DesignMtx(i, j) = BasisSplines(j, deg, xDisc[i], knots);
        }
    }

    TVectorD y(NumData);
    for (int i = 0; i < NumData; i++) {y(i) = ELossDisc[i];}

    TDecompSVD svd(DesignMtx);
    TVectorD Coeff(y);
    Bool_t ok = svd.Solve(Coeff);

    vector<double> Coefficients;
    for (int i = 0; i < NumBasis; i++) {Coefficients.push_back(Coeff(i));}

    return Coefficients;
}

void ConvolvedBraggModel(vector<double> &y, vector<double> &x,
vector<double> &knots, const vector<double> &SplineCoeff, int deg, int ParamShift, double angScale) {

    y.resize(x.size());
    for (int i = 0; i < x.size(); i++) {

        double xOriginal = x.at(i) / angScale + ParamShift;
        double evaluation = EvaluateSpline(xOriginal, knots, SplineCoeff, deg);
        y[i] = evaluation;
    }
}