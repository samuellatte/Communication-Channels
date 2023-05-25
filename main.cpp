#include "Constants.h"
#include "Objects.h"
#include <cfloat>
#include <cmath>
#include <chrono>
#include <fstream>
#include <random>

complex<double> E0 = 0., E1 = 0., E2 = 0., EGnd = 0., EDiff = 0., Voc = 0.;
const complex<double> j(0.0,1.0);
array<double,2> TX = {40,300};
double MatrixPrxdBm[300][80];
double MatrixSNRdB[300][80];
double MatrixDelaySpread[300][80];
double MatrixRiceFactor[300][80];
double PrxWatt = 0., PrxdBm = 0., SNRdB = 0., delaySpread = 0., K = 0., a0 = 0.;
double ClosestDist = DBL_MAX;
double deltaTau = 1/(2*BW);
double He = -Lambda/Pi;
double n = 0.;
double p0 = 0.;
double PathLoss = 0.;
double sigmaL = 0.;
double meanL = 0.;
double Lm = 142.;
double sensitivity = SnrUE + RXnoiseFig + 10*log10(kBoltz*Temp*BW) + 30; // dBm
int NbrReflections = 2;
int lMax = maxTau/deltaTau; // Number of taps
vector<complex<double>> TDL(lMax);
vector<int> tap = {};
vector<double> tdl = {};
vector<double> delays = {};
vector<double> RiceFactor = {};
vector<double> CIR = {};
vector<double> VectorPrxdBm = {};
vector<double> VectorSNRdB = {};
vector<double> VectorDelaySpread = {};
vector<double> VectorRiceFactor = {};
vector<double> VectorPathLoss = {};
vector<double> LogDistance = {};
vector<double> Distance = {};



// Initialize the walls
Wall HorizD({0,0},{80,0}),
VerticalLeft({20,10},{20,320}),
VerticalRight1({60,10},{60,70}),
VerticalRight2({60,80},{60,180}),
VerticalRight3({60,190},{60,270}),
VerticalRight4({60,280},{60,300}),

StreetLeftU({0,10},{20,10}),
StreetRight1U({60,10},{80,10}),
StreetRight2U({60,80},{80,80}),
StreetRight2D({60,70},{80,70}),
StreetRight3U({60,190},{80,190}),
StreetRight3D({60,180},{80,180}),
StreetRight4U({60,280},{80,280}),
StreetRight4D({60,270},{80,270});

// Initialize corners where there is a diffraction

Corner LeftCorner({20,10}),
RightCorner1({60,10}),
RightCorner2({60,80}),
RightCorner3({60,190}),
RightCorner4({60,280});

// Initialize outside zones

Rectangle rectL({0,10},{20,320}),
rectR1({60,10},{80,70}),
rectR2({60,80},{80,180}),
rectR3({60,190},{80,270}),
rectR4({60,280},{80,300});

// Initialize points at the intersection with each of the crossing streets

MiddleStreet mid1({20,5}),
mid2({60,5}),
mid3({60,75}),
mid4({60,185}),
mid5({60,275});

// Add all the parameters to the list Walls
array<Wall,14> Walls {{ {HorizD} , {VerticalLeft}, {VerticalRight1}, {VerticalRight2}, {VerticalRight3}, {VerticalRight4}, {StreetLeftU}, {StreetRight1U}, {StreetRight2U}, {StreetRight2D}, {StreetRight3U}, {StreetRight3D}, {StreetRight4U}, {StreetRight4D} }};

// Add all corners to a list
array<Corner,5> Corners {{ {LeftCorner}, {RightCorner1}, {RightCorner2}, {RightCorner3}, {RightCorner4} }};

// Add all outside zones to the list
array<Rectangle,5> Zones {{ {rectL}, {rectR1}, {rectR2}, {rectR3}, {rectR4} }};

// Add all intersection points to the list
array<MiddleStreet,5> Middle {{ {mid1} , {mid2} , {mid3} , {mid4} , {mid5} }};


void LineOfSight(double NormD0){
    double tauLOS = 0.;
    int l = 0;
    
    tauLOS = NormD0/C;
    E0 = sqrt(60 * EIRPmax) * (exp(-Beta * NormD0 * j) / NormD0);
    Voc += He * E0;
    a0 = abs(He * E0);
    delays.push_back(tauLOS);
    CIR.push_back(a0);
    l = ceil((tauLOS/maxTau) * lMax);
    TDL[l] += He * E0 * sinc(2*BW*(tauLOS - l*deltaTau));
}

void GroundReflection(double NormD0){
    array<double,2> TXvert, RXvert, sGnd, imageGnd, dGnd, dGndUnit;
    double tauGnd = 0., CosThetaI = 0., SinThetaI = 0., RC = 0., NormDGnd = 0., HeGND = 0., ThetaI = 0., ThetaGnd = 0., Gtx = 0., Ptx = 0., a = 0.;
    int l = 0;
    
    TXvert = {0,2};RXvert = {NormD0,2};
    Wall ground({0,0},{NormD0,0});
    sGnd = ArraySubstract(ground.Origin, TXvert);
    imageGnd = ArraySubstract(ConstArrayProduct(2*InnerProduct(sGnd, ground.n), ground.n), TXvert);
    dGnd = ArraySubstract(imageGnd, RXvert);
    NormDGnd = NormArrays(dGnd);
    tauGnd = NormDGnd/C;
    dGndUnit = UnitaryVector(dGnd, NormDGnd);
    CosThetaI = abs(InnerProduct(dGndUnit, ground.n));
    SinThetaI = sqrt(1 - pow(CosThetaI,2));
    ThetaI = asin(SinThetaI);
    ThetaGnd = 180 - ThetaI;
    RC = (CosThetaI - 1/sqrt(EpsilonR) * sqrt(1 - 1/EpsilonR * pow(SinThetaI,2))) / (CosThetaI + 1/sqrt(EpsilonR) * sqrt(1 - 1/EpsilonR * pow(SinThetaI,2)));
    Gtx = (16/3*Pi) * pow(sin(ThetaGnd),3);
    Ptx = EIRPmax/Gtx;
    EGnd = RC * sqrt(60 * Ptx * Gtx) * (exp(-Beta * NormDGnd * j) / NormDGnd);
    HeGND = (-Lambda/Pi) * cos((Pi/2) * cos(ThetaI))/pow(SinThetaI,2);
    Voc += HeGND * EGnd;
    a = abs(HeGND * EGnd);
    delays.push_back(tauGnd);
    RiceFactor.push_back(pow(a,2));
    CIR.push_back(a);
    l = ceil((tauGnd/maxTau) * lMax);
    TDL[l] += HeGND * EGnd * sinc(2*BW*(tauGnd - l*deltaTau));;
}

double ComputeReflectionCoeff(array<double,2> d, double NormD,  int i){
    double CosThetaI = 0., SinThetaI = 0., RC = 0.;
    array<double,2> dUnit;
    
    dUnit = UnitaryVector(d, NormD);
    CosThetaI = abs(InnerProduct(dUnit, Walls[i].n));
    SinThetaI = sqrt(1 - pow(CosThetaI,2));
    RC = (CosThetaI - sqrt(EpsilonR) * sqrt(1 - 1/EpsilonR * pow(SinThetaI,2))) / (CosThetaI + sqrt(EpsilonR) * sqrt(1 - 1/EpsilonR * pow(SinThetaI,2)));
    return RC;
}

complex<double> ComputeDiffractionCoeff(double deltaR){
    double nu = 0., ModulusFnudB = 0., ModulusFnu = 0., AngleFnu = 0.;
    complex<double> Fnu = 0.;
    
    nu = sqrt((2/Pi) * Beta * deltaR);
    ModulusFnudB = -6.9 - 20 * log10(sqrt(pow(nu-0.1, 2)+1)+nu-0.1); // in dB
    ModulusFnu = sqrt(pow(10, ModulusFnudB/20)); // in linear scale
    AngleFnu = -(Pi/4) - (Pi * pow(nu,2) / 2);
    Fnu = polar(ModulusFnu,AngleFnu);
     
    return Fnu;
}

bool CheckTransmission(array<double,2> O, array<double,2> E){
    array<double,2> s,d;
    double t;
    
    d = ArraySubstract(O, E);
    for (int i = 0; i < 14; i++){
        if(Walls[i].u[0] * d[1] != Walls[i].u[1] * d[0]){
            t = (d[1] * (O[0] - Walls[i].Origin[0]) - d[0] * (O[1] - Walls[i].Origin[1])) / (Walls[i].u[0] * d[1] - Walls[i].u[1] * d[0]);
            s = ArraySubstract(Walls[i].Origin, O);
            if(t > 0 && t < NormArrays((ArraySubstract(Walls[i].Origin, Walls[i].End)))){
                if(Sign(InnerProduct(s, Walls[i].n)) != Sign(InnerProduct(ArraySubstract(Walls[i].Origin, E), Walls[i].n)) && Sign(InnerProduct(ArraySubstract(Walls[i].Origin, E), Walls[i].n)) != 0 && Sign(InnerProduct(s, Walls[i].n)) != 0){
                    return true;
                }
            }
        }
    }
    return false;
}

bool notInZone(array<Rectangle,5> Zones, array<double,2> RX){
    for(int i = 0; i < 5; i++){
        if (RX[0] >= Zones[i].bottomLeft[0] and RX[0] <= Zones[i].topRight[0] and RX[1] >= Zones[i].bottomLeft[1] and RX[1] <= Zones[i].topRight[1])
            return false;
    }
    return true;
}

void ComputeWaves(array<double,2> RX){
    array<double,2> d0;
    double NormD0 = 0., NormDiff = 0., a = 0., s1 = 0., s2 = 0., tauDiff = 0., deltaR = 0.;
    complex<double> DiffractionCoeff = 0.;
    bool TransmissionTXRX;
    int ind = 0, l = 0;
    
    if(NbrReflections >= 0){ //0 Reflection
        
        d0 = ArraySubstract(TX, RX);
        NormD0 = NormArrays(d0);
        
        TransmissionTXRX = CheckTransmission(TX, RX);
        if(not TransmissionTXRX){
            LineOfSight(NormD0);
            GroundReflection(NormD0);
        }
        else{
            for(int i = 0; i < 5; i++){
                NormDiff = NormArrays(ArraySubstract(RX, Corners[i].Point));
                if(NormDiff < ClosestDist){
                    ClosestDist = NormDiff;
                    ind = i;
                }
            }
            s1 = NormArrays(ArraySubstract(TX, Corners[ind].Point));
            s2 = NormArrays(ArraySubstract(Corners[ind].Point, RX));
            tauDiff = (s1 + s2)/C;
            deltaR = s1 + s2 - NormD0;
            DiffractionCoeff = ComputeDiffractionCoeff(deltaR);
            EDiff = DiffractionCoeff * sqrt(60.0 * EIRPmax) * (exp(-Beta * NormD0 * j) / (NormD0));
            Voc += He * EDiff;
//            a = abs(He * EDiff);
//            CIR.push_back(a); // No phase shift for the moment)
            delays.push_back(DBL_MAX);
            delays.push_back(0);
            
            a0 = 0;
            ClosestDist = DBL_MAX;
            ind = 0;
        }
        
        if(1 <= NbrReflections){ //1 Reflection
            
            for(int i = 0; i < 14; i++){
                array<double,2> s1, image1, d1, P1;
                double NormD1 = 0., t1 = 0., tau1 = 0., ReflectionCoeff1 = 0.;
                bool TransmissionTXP1, TransmissionP1RX;
                
                s1 = ArraySubstract(Walls[i].Origin, TX);
                image1 = ArraySubstract(ConstArrayProduct(2*InnerProduct(s1, Walls[i].n), Walls[i].n), TX);
                d1 = ArraySubstract(image1, RX);
                NormD1 = NormArrays(d1);
                tau1 = NormD1/C;
                if(Walls[i].u[0] * d1[1] != Walls[i].u[1] * d1[0]){
                    t1 = (d1[1] * (image1[0] - Walls[i].Origin[0]) - d1[0] * (image1[1] - Walls[i].Origin[1])) / (Walls[i].u[0] * d1[1] - Walls[i].u[1] * d1[0]);
                    if((0 < t1) && (t1 < NormArrays(ArraySubstract(Walls[i].Origin, Walls[i].End)))
                       && (Sign(InnerProduct(s1, Walls[i].n)) == Sign(InnerProduct(ArraySubstract(Walls[i].Origin, RX),Walls[i].n)))){
                        P1[0] = Walls[i].Origin[0] + t1 * Walls[i].u[0];
                        P1[1] = Walls[i].Origin[1] + t1 * Walls[i].u[1];
                        TransmissionTXP1 = CheckTransmission(TX, P1);
                        TransmissionP1RX = CheckTransmission(P1, RX);
                        if(not TransmissionTXP1 and not TransmissionP1RX){
                            ReflectionCoeff1 = ComputeReflectionCoeff(d1, NormD1, i);
                            E1 = ReflectionCoeff1 * sqrt(60.0 * EIRPmax) * (exp(-Beta * NormD1 * j) / (NormD1));
                            Voc += He * E1;
                            a = abs(He * E1);
                            CIR.push_back(a); // No phase shift for the moment)
                            RiceFactor.push_back(pow(a,2));
                            delays.push_back(tau1);
                            l = ceil((tau1/maxTau) * lMax);
                            TDL[l] += He * E1 * sinc(2*BW*(tau1 - l*deltaTau));
                        }
                    }
                            
                    if(2 <= NbrReflections){ //2 Reflection
                        
                        for(int k = 0; k < 14; k++){
                            if(k != i){
                                array<double,2> s2, image2, d2, P2;
                                double NormD2 = 0., t2 = 0., tau2 = 0.;
                                complex<double> ReflectionCoeff1 = 0., ReflectionCoeff2 = 0.;
                                bool TransmissionTXP1, TransmissionP1P2, TransmissionP2RX;
                                
                                s2 = ArraySubstract(Walls[k].Origin, image1);
                                image2 = ArraySubstract(ConstArrayProduct(2*InnerProduct(s2, Walls[k].n), Walls[k].n), image1);
                                d2 = ArraySubstract(image2, RX);
                                NormD2 = NormArrays(d2);
                                tau2 = NormD2/C;
                                if(Walls[k].u[0] * d2[1] != Walls[k].u[1] * d2[0]){
                                    t2 = (d2[1] * (image2[0] - Walls[k].Origin[0]) - d2[0] * (image2[1] - Walls[k].Origin[1])) / (Walls[k].u[0] * d2[1] - Walls[k].u[1] * d2[0]);
                                    if((0 < t2) && (t2 < NormArrays(ArraySubstract(Walls[k].Origin, Walls[k].End)))
                                       && (Sign(InnerProduct(s2, Walls[k].n)) == Sign(InnerProduct(ArraySubstract( Walls[k].Origin, RX),Walls[k].n)))){
                                        P2[0] = Walls[k].Origin[0] + t2 * Walls[k].u[0];
                                        P2[1] = Walls[k].Origin[1] + t2 * Walls[k].u[1];
                                        d1 = ArraySubstract(image1, P2);
                                        t1 = (d1[1] * (image1[0] - Walls[i].Origin[0]) - d1[0] * (image1[1] - Walls[i].Origin[1])) / (Walls[i].u[0] * d1[1] - Walls[i].u[1] * d1[0]);
                                        NormD1 = NormArrays(d1);
                                        if((0 < t1) && (t1 < NormArrays(ArraySubstract(Walls[i].Origin, Walls[i].End))) && (Sign(InnerProduct(s1, Walls[i].n)) == Sign(InnerProduct(ArraySubstract(Walls[i].Origin, P2),Walls[i].n)))){
                                            P1[0] = Walls[i].Origin[0] + t1 * Walls[i].u[0];
                                            P1[1] = Walls[i].Origin[1] + t1 * Walls[i].u[1];
                                            TransmissionTXP1 = CheckTransmission(TX, P1);
                                            TransmissionP1P2 = CheckTransmission(P1, P2);
                                            TransmissionP2RX = CheckTransmission(P2, RX);
                                            if(not TransmissionTXP1 and not TransmissionP1P2 and not TransmissionP2RX){
                                                ReflectionCoeff1 = ComputeReflectionCoeff(d1, NormD1, i);
                                                ReflectionCoeff2 = ComputeReflectionCoeff(d2, NormD2, k);
                                                E2 = ReflectionCoeff1 * ReflectionCoeff2 * sqrt(60.0 * EIRPmax) * (exp(-Beta * NormD2 * j) / (NormD2));
                                                Voc += He * E2;
                                                a = abs(He * E2);
                                                CIR.push_back(a); // No phase shift for the moment)
                                                RiceFactor.push_back(pow(a,2));
                                                delays.push_back(tau2);
                                                l = ceil((tau2/maxTau) * lMax);
                                                TDL[l] += He * E2 * sinc(2*BW*(tau2 - l*deltaTau));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    PrxWatt = pow(abs(Voc),2)/(8*Resistance);
}

void ComputeDelaySpread(){
    delaySpread = *max_element(delays.begin(), delays.end()) - *min_element(delays.begin(), delays.end());
}

void ComputeRiceFactor(){
    double s = 0.;
    
    s = ArraySum(RiceFactor);
    if(s == 0){
        s = 1;
    }
    K = 10*log10(pow(a0,2)/s);
}

void ComputeTDL(){
    double s = TDL.size();
    
    for(int i = 0; i < s; i++){
        if(abs(TDL[i]) > 0){
            tap.push_back(i);
            tdl.push_back(abs(TDL[i]));
        }
    }
}

void ComputePathLoss(){
    double s = LogDistance.size(), meanDist = 0., meanPrxdBm = 0., cov = 0., var = 0., beta0 = 0., beta1 = 0.;
    ofstream FilePathLoss ("FilePathLoss.txt"), FileParameters ("FileParameters.txt"), FileFading ("FileFading.txt");
    vector<double> fading = {};
    
    meanDist = ArraySum(LogDistance)/s; meanPrxdBm = ArraySum(VectorPrxdBm)/s;
    
    if (FilePathLoss.is_open() && FileParameters.is_open() && FileFading){
        for(int i = 0; i < s; i++){
            cov += (LogDistance[i]-meanDist)*(VectorPrxdBm[i]-meanPrxdBm);
            var += pow((LogDistance[i]-meanDist),2);
        }
        beta1 = cov/var;
        beta0 = meanPrxdBm - beta1 * meanDist;
        n = -beta1/10;
        p0 = beta0;
        FileParameters << n << endl;
        FileParameters << p0 << endl;
        for(int i = 0; i < s; i++){
            PathLoss = -10*n * LogDistance[i] + p0;
            VectorPathLoss.push_back(PathLoss);
            FilePathLoss << PathLoss << endl;
        }
        
        for(int i = 0; i < s; i++){
            sigmaL += pow((VectorPrxdBm[i]-VectorPathLoss[i]),2);
            
        }
        sigmaL = sqrt(sigmaL/s);
        FileParameters << sigmaL << endl;
        
        for(int i = 0; i < s; i++){
            fading.push_back(VectorPrxdBm[i]-VectorPathLoss[i]);
            FileFading << fading[i] << endl;
        }
        meanL = ArraySum(fading)/s;
        FileParameters << meanL << endl;
        FilePathLoss.close(); FileParameters.close(); FileFading.close();
    }
    else cout << "Unable to open file";
}

void ComputeCellRange(){
    double s = VectorPrxdBm.size();
    vector<double> gamma = {}, probability = {};
    ofstream FileCellRange ("FileCellRange.txt");
            
    if (FileCellRange.is_open()){
        for(int i = 0; i < s; i++){
            gamma.push_back(VectorPathLoss[i]-sensitivity);
            probability.push_back(1. - 0.5*erfc(gamma[i] / (sigmaL*sqrt(2))));
            FileCellRange << probability[i] << endl;
        }
        FileCellRange.close();
    }
    else cout << "Unable to open file";
}

void ComputeOkumuraHataModel(){
    double s = Distance.size(), hTX = 30., n = 0., d0 = 1000, f = F/1000000, a = 0., sigmaLOH = 0., meanL0 = 0.;
    vector<double> L0 = {}, gamma = {}, probability = {}, fading = {};
    ofstream FilePathLossOH ("FilePathLossOH.txt"), FileParametersOH ("FileParametersOH.txt"), FileFadingOH ("FileFadingOH.txt"), FileCellRangeOH ("FileCellRangeOH.txt");
    
    a = 3.2*pow(log10(11.75*H), 2) - 4.97;
    n = 4.49 - 0.655*log10(hTX);
    
    if (FilePathLossOH.is_open() && FileParametersOH.is_open() && FileFadingOH && FileCellRangeOH){
        for(int i = 0; i < s; i++){
            L0.push_back(-abs(69.55 + 10*n*log10(Distance[i]/d0) + 26.16*log10(f) - 13.82*log10(hTX) - a + 30.));
            FilePathLossOH << L0[i] << endl;
        }
        
        sigmaLOH = 0.65*pow(log10(f), 2) - 1.3*log10(f) + 5.2;
        FileParametersOH << n << endl;
        FileParametersOH << L0[0] << endl;
        FileParametersOH << sigmaLOH << endl;
        
        for(int i = 0; i < s; i++){
            gamma.push_back(VectorPathLoss[i]-sensitivity);
            probability.push_back(1. - (1./2.) * erfc(gamma[i] / (sigmaLOH*sqrt(2))));
            FileCellRangeOH << probability[i] << endl;
        }
        
        for(int i = 0; i < s; i++){
            fading.push_back(VectorPrxdBm[i]-L0[i]);
            FileFadingOH << fading[i] << endl;
        }
        meanL0 = ArraySum(fading)/s;
        FileParametersOH << meanL0 << endl;
        
        FilePathLossOH.close(); FileParametersOH.close(); FileFadingOH.close(); FileCellRangeOH.close();
    }
    else cout << "Unable to open file";
    
}

void ComputeWholeCellCoverage(){
    double s = VectorPathLoss.size(), a = 0., b = 0., Lr = 0.;
    vector<double> Fu = {}, gamma = {};
    ofstream FileGamma ("FileGamma.txt"), FileWholeCellCoverage ("FileWholeCellCoverage.txt");
    
    if (FileGamma.is_open() && FileWholeCellCoverage.is_open()){
        for(int i = 0; i < s; i++){
            gamma.push_back(VectorPathLoss[i]-sensitivity);
            FileGamma << gamma[i] << endl;
            Lr = Lm - gamma[i];
            a = (Lm - Lr)/(sqrt(2)*sigmaL);
            b = (1/(sqrt(2)*sigmaL)) * 10*n*log10(exp(1));
            Fu.push_back(1 - 0.5*erfc(a) + 0.5*exp( (2*a/b) + (1/(pow(b, 2))) )*erfc(a + (1/b)));
            FileWholeCellCoverage << Fu[i] << endl;
        }
        FileGamma.close(); FileWholeCellCoverage.close();
    }
    else cout << "Unable to open file";
}

void ComputeParameters1D(){
    int i;
    array<double,2> RX;
    ofstream FileVectorPrxdBm ("FileVectorPrxdBm.txt"), FileVectorSNR ("FileVectorSNR.txt"), FileVectorDelaySpread ("FileVectorDelaySpread.txt"), FileVectorRiceFactor ("FileVectorRiceFactor.txt");

    if (FileVectorPrxdBm.is_open() && FileVectorSNR.is_open() && FileVectorDelaySpread.is_open() && FileVectorRiceFactor.is_open()){
        for(double y = 289.5; y > 0; y--){ // Iterate over columns
            RX = {TX[0],y};
//            cout << RX[0] << " ";
//            cout << RX[1] << endl;
            i = y-0.5;
            if(NormArrays(ArraySubstract(TX, RX)) > 10){
                ComputeWaves(RX);
                PrxdBm = 10 * log10(PrxWatt * 1e3);
                VectorPrxdBm.push_back(PrxdBm);
                VectorSNRdB.push_back(SNRdB);
                Distance.push_back(NormArrays(ArraySubstract(TX, RX)));
                LogDistance.push_back(log10(NormArrays(ArraySubstract(TX, RX))));
                SNRdB = PrxdBm - RXnoiseFig - 10*log10(kBoltz*Temp*BW) - 30;
                ComputeDelaySpread();
                ComputeRiceFactor();
                VectorDelaySpread.push_back(delaySpread);
                VectorRiceFactor.push_back(K);
                FileVectorPrxdBm << PrxdBm << endl;
                FileVectorSNR << SNRdB << endl;
                FileVectorDelaySpread << delaySpread << endl;
                FileVectorRiceFactor << K << endl;
            }
            PrxWatt = 0.; PrxdBm = 0.; Voc = 0.; SNRdB = 0.; K = 0.; delaySpread = 0.; a0 = 0.; RiceFactor = {}; delays = {};
        }
        ComputePathLoss();
        ComputeCellRange();
//        ComputeOkumuraHataModel();
        ComputeWholeCellCoverage();

        FileVectorPrxdBm.close(); FileVectorSNR.close(); FileVectorDelaySpread.close(); FileVectorRiceFactor.close();
    }
    else cout << "Unable to open file";
}

void ComputeParameters2D(){
    int i,j;
    array<double,2> RX;
    ofstream FilePrxdBm ("FilePrxdBm.txt"), FileSNRdB ("FileSNRdB.txt"), FileDelaySpread ("FileDelaySpread.txt"), FileRiceFactor ("FileRiceFactor.txt");

    if (FilePrxdBm.is_open() && FileSNRdB.is_open() && FileDelaySpread.is_open() && FileRiceFactor.is_open()){
        for(double x = 0.5; x < 80; x++){ // Iterate over columns
            for(double y = 0.5; y < 300; y++){ // Iterate over rows
                RX = {x,y};
                i = x - 0.5;
                j = y - 0.5;
                
                if(notInZone(Zones, RX) && NormArrays(ArraySubstract(TX, RX)) > 10){
                    ComputeWaves(RX);
                    PrxdBm = 10 * log10(PrxWatt * 1e3);
                    SNRdB = PrxdBm - RXnoiseFig - 10*log10(kBoltz*Temp*BW) - 30;
                    ComputeDelaySpread();
                    ComputeRiceFactor();
                    MatrixPrxdBm[j][i] = PrxdBm;
                    MatrixSNRdB[j][i] = SNRdB;
                    MatrixDelaySpread[j][i] = delaySpread;
                    MatrixRiceFactor[j][i] = K;
                    FilePrxdBm << PrxdBm << endl;
                    FileSNRdB << SNRdB << endl;
                    FileDelaySpread << delaySpread << endl;
                    FileRiceFactor << K << endl;
                }
                else{
                    PrxWatt = 0.; PrxdBm = log10(0);
                    SNRdB = PrxdBm - RXnoiseFig - 10*log10(kBoltz*Temp*BW);
                    delaySpread = DBL_MAX;
                    K = log10(0);
                    MatrixPrxdBm[j][i] = PrxdBm;
                    MatrixSNRdB[j][i] = SNRdB;
                    MatrixDelaySpread[j][i] = delaySpread;
                    MatrixRiceFactor[j][i] = K;
                    FilePrxdBm << PrxdBm << endl;
                    FileSNRdB << SNRdB << endl;
                    FileDelaySpread << delaySpread << endl;
                    FileRiceFactor << K << endl;
                }
                PrxWatt = 0.; PrxdBm = 0.; Voc = 0.; SNRdB = 0.; K = 0.; delaySpread = 0.; a0 = 0.; RiceFactor = {}; delays = {};
            }
            i = 0;
            j = 0;
        }
        FilePrxdBm.close(); FileSNRdB.close(); FileDelaySpread.close();FileRiceFactor.close();
    }
    else cout << "Unable to open file";
}

void ComputeParametersMiddleStreet(){
    
    ComputeWaves(Middle[4].Point);
    ComputeTDL();
    ofstream FileTau ("FileTau.txt"), FileAlpha ("FileAlpha.txt"), FileTDL ("FileTDL.txt"), FileTap ("FileTap.txt");
    double s1 = delays.size(), s2 = tap.size();
    
    // Physical channel model
    if (FileTau.is_open() && FileAlpha.is_open()){
        for(int i = 0; i < s1; i++){
            FileTau << delays[i] << endl;
            FileAlpha << CIR[i] << endl;
        }
        FileTau.close(); FileAlpha.close();
    }
    else cout << "Unable to open file";
    
    // Tapped Delay Line
    if (FileTDL.is_open() && FileTap.is_open()){
        for(int i = 0; i < s2; i++){
            FileTap << tap[i] << endl;
            FileTDL << tdl[i] << endl;
        }
        FileTDL.close(); FileTap.close();
    }
    else cout << "Unable to open file";
}

int main(){
    
    int plot = 3; // Select what to plot
    
    switch(plot){
            
        case 1:
            // 1D plot
            ComputeParameters1D();
            break;
            
        case 2:
            // 2D plot
            ComputeParameters2D();
            break;
            
        case 3:
            // Channel impulse response
            ComputeParametersMiddleStreet();
            break;
    }
    
    return 0;
}

