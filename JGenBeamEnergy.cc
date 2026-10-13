#include "JGenBeamEnergy.h"
#include <TSystem.h>

const double kaonmass=0.49761;
const double neutronmass=0.93957;
const double photonmass=0.0;
//const string KaonFile="BeamProfile_kaons.root";
const string KaonFile=Form("%s/BeamProfile_kaons_2024.root",gSystem->Getenv("KLGEN"));
const string NeutronFile=Form("%s/BeamProfile_neutrons.root",gSystem->Getenv("KLGEN"));
const string PhotonFile=Form("%s/BeamProfile_photons.root",gSystem->Getenv("KLGEN"));
const string KaonHistname="h24E_2";
const string NeutronHistname="BeamProfile";
const string PhotonHistname="BeamProfile";

    // constructor parses -E option
JGenBeamEnergy::JGenBeamEnergy (char * coption_) {
    int   narg = 0;
    char  coption [256];
    char* copt [20];
    char* s = coption;
    char* keyword [] = { "mono", "plain", "histo", "kaon", "neutron", "photon", NULL }; //distribution keyword
                                                                                        //char* keywordT [] = { "kaon", "neutron", "photon", NULL }; //type keyword
    int   jKey=0;
    int   jKeyy=0;
    int   jKeyD = 0; //distribution (mono plain histo)
    int   jKeyT = 0; //Type (kaon, neutron, photon)
    
    strncpy (coption, coption_, 255);   // copy line before manipulating
    int   len = strlen(coption);
    
        //places option into copt
    while ((s=strtok(s, ":")) != NULL) {  // Split string into tokens based on the delimiter ':'
        copt [narg++] = s;
        if ((len -= (strlen (s) + 1)) <= 0) break;
        s += strlen (s) + 1;
    }
    
    if (narg == 1){
        while ( keyword [jKey] && strcasecmp (keyword[jKey], copt [0])) {
            jKey++;
        }
        if (jKey<3){
            jKeyD=jKey;
            jKeyT=6;
        }
        else if (jKey<6){
            jKeyT=jKey;
            jKeyD=6;
        }
        else {
            jKeyD=6;
            jKeyT=6;
        }
    }
    else if (narg > 1){
        while ( keyword [jKey] && strcasecmp (keyword[jKey], copt [0])) {
            jKey++;
        }
        while ( keyword [jKeyy] && strcasecmp (keyword[jKeyy], copt [1])) {
            jKeyy++;
        }
        if (jKey==6){
            jKeyD=6;
            jKeyT=6;
        }
        else if (jKey<3&&jKeyy==6){
            jKeyD=jKey;
            jKeyT=6;
        }
        else if (jKey>2&&jKey<6&&jKeyy==6){
            jKeyT=jKey;
            jKeyD=6;
        }
        else if (jKey<3&&jKeyy>2&&jKeyy<6){
            jKeyD=jKey;
            jKeyT=jKeyy;
        }
        else if (jKey>2&&jKey<6&&jKeyy<3){
            jKeyD=jKeyy;
            jKeyT=jKey;
        }
        else {
            jKeyD=6;
            jKeyT=6;
        }
    }
        //cout<<"Type: "<<jKeyT<<" Dist: "<<jKeyD<<endl;
    
    if (!keyword [jKeyD] && !keyword [jKeyT]) {   // no keyword found, assume float argument
        switch (narg) {
            case 1:
                constructor (kaon, mono, atof(copt[0]), atof(copt[0]));
                break;
            case 2:
                constructor (kaon, plain, atof(copt[0]), atof(copt[1]));
                break;
            default:
                throw XnumberArgument;
        }
    }
    else if (!keyword [jKeyD] && keyword [jKeyT]) {   // no keyword distribution found, assume float argument
        switch ((EdistribBeamType_t) jKeyT) {
            case kaon:
                if (narg == 1)
                    constructor (kaon, KaonFile, KaonHistname, 0, 12);
                else if (narg == 2)
                    constructor (kaon, mono, atof(copt[1]), atof(copt[1]));
                else if (narg == 3)
                    constructor (kaon, plain, atof(copt[1]), atof(copt[2]));
                else
                    throw XnumberArgument;
                break;
            case neutron:
                if (narg == 1)
                    constructor (neutron, NeutronFile, NeutronHistname, 0,12);
                else if (narg == 2)
                    constructor (neutron, mono, atof(copt[1]), atof(copt[1]));
                else if (narg == 3)
                    constructor (neutron, plain, atof(copt[1]), atof(copt[2]));
                else
                    throw XnumberArgument;
                break;
            case photon:
                if (narg == 1)
                    constructor (photon, PhotonFile, PhotonHistname, 0,12);
                else if (narg == 2)
                    constructor (photon, mono, atof(copt[1]), atof(copt[1]));
                else if (narg == 3)
                    constructor (photon, plain, atof(copt[1]), atof(copt[2]));
                else
                    throw XnumberArgument;
                break;
            default:
                throw XnumberArgument;
        }
    }
    else if (keyword [jKeyD] && !keyword [jKeyT]) {   // no keyword distribution found, assume float argument
        switch ((EdistribBeamType_t) jKeyD) {
            case mono:
                if (narg != 2) throw XnumberArgument;
                constructor (kaon, mono, atof(copt[1]), atof(copt[1]));
                break;
            case plain:
                if (narg != 3) throw XnumberArgument;
                constructor (kaon, plain, atof(copt[1]), atof(copt[2]));
                break;
            case histo:
                if (narg == 1 || narg == 2 )
                    constructor (kaon, KaonFile, KaonHistname, 0, 12);
                else if (narg == 3)
                    constructor (kaon, KaonFile, KaonHistname, atof(copt[1]), atof(copt[2]));
                else throw XnumberArgument;
                break;
            default:
                throw XnumberArgument;
        }
    }
    else  {   // keywords found
        switch ((EdistribBeamType_t) jKeyD) {
            case mono:
                if (narg != 3) throw XnumberArgument;
                constructor ((EdistribBeamType_t) jKeyT, (EdistribBeamType_t) jKeyD, atof(copt[2]), atof(copt[2]));
                break;
            case plain:
                if (narg != 4) throw XnumberArgument;
                constructor ((EdistribBeamType_t) jKeyT, (EdistribBeamType_t) jKeyD, atof(copt[2]), atof(copt[3]));
                break;
            case histo:
                if (narg == 2 ||narg == 3){// throw XnumberArgument;
                    switch ((EdistribBeamType_t) jKeyT) {
                        case kaon:
                            constructor (kaon, KaonFile, KaonHistname, 0, 12);
                            break;
                        case neutron:
                            constructor (neutron, NeutronFile, NeutronHistname, 0, 12);
                            break;
                        case photon:
                            constructor (photon, PhotonFile, PhotonHistname, 0, 12);
                            break;
                    }
                }
                else if (narg == 4){// throw XnumberArgument;
                    switch ((EdistribBeamType_t) jKeyT) {
                        case kaon:
                            constructor (kaon, KaonFile, KaonHistname, atof(copt[2]), atof(copt[3]));
                            break;
                        case neutron:
                            constructor (neutron, NeutronFile, NeutronHistname, atof(copt[2]), atof(copt[3]));
                            break;
                        case photon:
                            constructor (photon, PhotonFile, PhotonHistname,  atof(copt[2]), atof(copt[3]));
                            break;
                    }
                }
		//specify histogram name
		else if (narg == 5){// throw XnumberArgument;
                    switch ((EdistribBeamType_t) jKeyT) {
                        case kaon:
                            constructor (kaon, KaonFile, copt[4], atof(copt[2]), atof(copt[3]));
                            break;
                        case neutron:
                            constructor (neutron, NeutronFile, copt[4], atof(copt[2]), atof(copt[3]));
                            break;
                        case photon:
                            constructor (photon, PhotonFile, copt[4],  atof(copt[2]), atof(copt[3]));
                            break;
                    }
                }
		//specify file and histogram name
		else if (narg == 6){// throw XnumberArgument;
                    switch ((EdistribBeamType_t) jKeyT) {
                        case kaon:
			  //cout << "Sampling " << copt[5] << " from " << copt[4] << endl;
                            constructor (kaon, copt[4], copt[5], atof(copt[2]), atof(copt[3]));
                            break;
                        case neutron:
                            constructor (neutron, copt[4], copt[5], atof(copt[2]), atof(copt[3]));
                            break;
                        case photon:
                            constructor (photon, copt[4], copt[5],  atof(copt[2]), atof(copt[3]));
                            break;
                    }
                }		

                break;
                
        }
    }
}

void JGenBeamEnergy::constructor (EdistribBeamType_t typeB_, EdistribBeamType_t typeE_,
                                  Double_t Emin_, Double_t Emax_) {
    typeE = typeE_;
    typeB = typeB_;
    Emin = Emin_;
    Emax = Emax_;
    r = new TRandom3(0);
    if(gRandom) delete gRandom;
    gRandom = r;
}


void JGenBeamEnergy::constructor (EdistribBeamType_t typeB_, string filename, string hstname, Double_t Emin_, Double_t Emax_) {
    typeE = histo;
    typeB = typeB_;
    Emin = Emin_;
    Emax = Emax_;
    genDir = new TDirectory ("gendir", "gendir");
    TFile f(filename.c_str());                    // open a file switches directory
    if (!f.IsOpen ()) throw XfileNotOpen;
    TH1F* dummy = (TH1F*) f.Get(hstname.c_str()); // access to histogram
    if (dummy == NULL) throw XhistNotFound;
    
    genDir -> cd ();                     // change dir back to root-memory
    hgen = new TH1F();                   // allocate memory
    *hgen = * ((TH1F*) dummy->Clone ()); // copy histogram
    r = new TRandom3(0);
    if(gRandom) delete gRandom;
    gRandom = r;
    
}

JGenBeamEnergy::JGenBeamEnergy(EdistribBeamType_t typeB_, EdistribBeamType_t typeD_, Double_t Emin_, Double_t Emax_){
    switch ((EdistribBeamType_t) typeD_) {
        case mono:
        case plain:
            constructor (typeB_, typeD_, Emin_, Emax_);
            break;
        case histo:
            switch ((EdistribBeamType_t) typeB_) {
                case kaon:
                    constructor (kaon, KaonFile, KaonHistname, Emin_, Emax_);
                    break;
                case neutron:
                    constructor (neutron, NeutronFile, NeutronHistname, Emin_, Emax_);
                    break;
                case photon:
                    constructor (photon, PhotonFile, PhotonHistname, Emin_, Emax_);
                    break;
            }
            break;
    }
}

Double_t JGenBeamEnergy::Generate () {
    
    switch (typeE) {
        case mono:
            Evalue = Emin;
            break;
        case plain:
            Evalue = Emin + (Emax - Emin) * r->Rndm();
            break;
        case histo:
            Evalue=-1;
            while (Evalue<Emin || Evalue>Emax)
                Evalue = hgen->GetRandom();
            break;
    }
    
    switch (typeB) {
        case kaon:
            P4beam.SetXYZM(0,0, TMath::Sqrt(Evalue*Evalue+2.0*Evalue*kaonmass), kaonmass);
            break;
        case neutron:
            P4beam.SetXYZM(0,0, TMath::Sqrt(Evalue*Evalue+2.0*Evalue*neutronmass), neutronmass);
            break;
        case photon:
            P4beam.SetXYZM(0,0, TMath::Sqrt(Evalue*Evalue+2.0*Evalue*photonmass), photonmass);
            break;
    }
    
    return Evalue;
}

