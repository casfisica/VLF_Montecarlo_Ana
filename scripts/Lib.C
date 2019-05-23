#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include <string>
#include <iostream>
#include "TTreeReader.h"
#include <tuple>
#include <vector>


class Error_Propagation{
    protected:

    
    private:
    
    public:
    
    static Double_t SumErr(Double_t _x, Double_t _y,Double_t  _dx, Double_t _dy, Double_t & Error2){
        /*Return de x+y and the error in the reference Error */
        Error2= TMath::Sqrt(_dx*_dx+_dy*_dy);
        return _x+_y;
    }
    
    static Double_t MultErr(Double_t _x, Double_t _y,Double_t  _dx, Double_t _dy, Double_t & Error2){
        /*Return de x*y and the error in the reference Error */
        Error2= (_x*_y)*TMath::Sqrt((_dx/_x)*(_dx/_x)+(_dy/_y)*(_dy/_y));
        return _x*_y;
    }
    
    static Double_t DivErr(Double_t _x, Double_t _y,Double_t  _dx, Double_t _dy, Double_t & Error2){
        /*Return de x/y and the error in the reference Error */
        Error2= (_x/_y)*TMath::Sqrt((_dx/_x)*(_dx/_x)+(_dy/_y)*(_dy/_y));
        return _x/_y;
    }
    
    static Double_t PowErr(Double_t _x, Double_t _y,Double_t  _dx, Double_t & Error2){
        /*Return de x^y and the error in the reference Error */
        Error2= TMath::Power(_x,_y-1.0)*_y*_dx;
        return TMath::Power(_x,_y);
    }
    
};


class Significance {
   private:
    TH1D *    SigTH1;
    TH1D *    BackTH1;
    std::vector<TH1D *> ArBackTH1;
   
    public: 
    /*Constructors*/
    Significance(void):SigTH1(NULL),BackTH1(NULL){}
    Significance(TH1D * _SigTH1, TH1D * _BackTH1):SigTH1(_SigTH1),BackTH1(_BackTH1){}
    Significance(TH1D * _SigTH1, std::vector<TH1D *> _ArBackTH1):SigTH1(_SigTH1),ArBackTH1(_ArBackTH1){
        TH1D *    TH1Copy=(TH1D *) _SigTH1->Clone();
        BackTH1=(TH1D *) _SigTH1->Clone();
        BackTH1->Reset();//To get rid of the signal events
        BackTH1->Sumw2();//To get good error
        for (auto it:_ArBackTH1){
            //Get a sum off all background
            BackTH1->Add(it);
         }
    }
    /*Destructor*/
    ~Significance(){}
    
    /*Set Methods*/
    void SetSignal(TH1D * _SigTH1){
        SigTH1= (TH1D *) _SigTH1->Clone();
    }//end SetBackground
    
    void SetBackground(TH1D * _BackTH1){
        BackTH1= (TH1D *)_BackTH1->Clone();
    }//end SetBackground
    
    /*
    void AddBackground(TH1D * _BackTH1){
        BackTH1->Add(_BackTH1);
        }//end SetBackground
    */
    
    /*Get Methods*/
    std::vector<TH1D *> GetVector(void){
        
        return ArBackTH1;
    }
    
    TH1D * GetBackground(void){
        return BackTH1;
    }
    
    TH1D * GetSigMoreThan(void){
        int Nbines=BackTH1->GetNbinsX();       
        TH1D *    TH1CopyInt=(TH1D *) SigTH1->Clone();
        TH1CopyInt->SetNameTitle("Significance More Than", "Significance");
        gStyle->SetOptStat(0);
        for (int l=0; l<=Nbines; l++){
            Double_t intZig=0;
            Double_t SigErr,BacErr;
            Double_t IntSignal=SigTH1->IntegralAndError(l,Nbines,SigErr);
            Double_t IntBackground=BackTH1->IntegralAndError(l,Nbines,BacErr);
            Double_t sumSB=IntSignal+IntBackground;
            if (sumSB>0){
                intZig=IntSignal/(TMath::Sqrt(sumSB));
            }else{
                intZig=0;
            }
            //std::cout<<"l: "<<l<<std::endl;
            //std::cout<<"intZig: "<<intZig<<std::endl;
            TH1CopyInt->SetBinContent(l,intZig);
            //TH1Copy->SetMarkerStyle(kOpenCircle);           
            //Using Cuadrature to calculate the error
            Double_t IntZigErr=TMath::Sqrt(((((TMath::Power(IntBackground,2))*(TMath::Power(SigErr,2)))/4)+((TMath::Power(IntSignal,2))*(TMath::Power(BacErr,2))))/(TMath::Power(sumSB,3)));
            TH1CopyInt->SetBinContent(l,intZig);
            if (!std::isnan(intZig)){
                TH1CopyInt->SetBinError(l,IntZigErr);
            }else{
                TH1CopyInt->SetBinError(l,0);
                //TH1Copy->SetMarkerStyle(kOpenCircle);
            }
              
        }
        return TH1CopyInt;
    }//End GetSigMoreThan
    
    /*Significance LessThanType*/
    
    TH1D * GetSigLessThan(void){
        int Nbines=BackTH1->GetNbinsX();       
        TH1D *    TH1CopyInt=(TH1D *) SigTH1->Clone();
        TH1CopyInt->SetNameTitle("Significance Less Than", "Significance");
        gStyle->SetOptStat(0);
        for (int l=Nbines; l>=0; l--){
            Double_t intZig=0;
            Double_t SigErr,BacErr;
            Double_t IntSignal=SigTH1->IntegralAndError(0,l,SigErr);
            Double_t IntBackground=BackTH1->IntegralAndError(0,l,BacErr);
            Double_t sumSB=IntSignal+IntBackground;
            if (sumSB>0){
                intZig=IntSignal/(TMath::Sqrt(sumSB));
            }else{
                intZig=0;
            }
            //std::cout<<"l: "<<l<<std::endl;
            //std::cout<<"intZig: "<<intZig<<std::endl;
            TH1CopyInt->SetBinContent(l,intZig);
            //TH1Copy->SetMarkerStyle(kOpenCircle);           
            //Using Cuadrature to calculate the error
            Double_t IntZigErr=TMath::Sqrt(((((TMath::Power(IntBackground,2))*(TMath::Power(SigErr,2)))/4)+((TMath::Power(IntSignal,2))*(TMath::Power(BacErr,2))))/(TMath::Power(sumSB,3)));
            TH1CopyInt->SetBinContent(l,intZig);
            if (!std::isnan(intZig)){
                TH1CopyInt->SetBinError(l,IntZigErr);
            }else{
                TH1CopyInt->SetBinError(l,0);
                //TH1Copy->SetMarkerStyle(kOpenCircle);
            }
              
        }
        return TH1CopyInt; 
    }//End GetSigLessThan
    
    /*Get effMoreThan*/
    std::vector<TH1D *> GetEffMoreThan(void){
        Int_t Nbines=SigTH1->GetNbinsX();       
        std::vector<TH1D *> arrTH1Copy;
        std::vector<TH1D *> arrTH1Imput;
        arrTH1Imput =ArBackTH1; //Put the back
        //arrTH1Imput.push_back(SigTH1);
        arrTH1Imput.insert(arrTH1Imput.begin(), SigTH1);//Put the signal in the front
        //std::cout<<ArBackTH1.size()<<std::endl; 
        //std::cout<<arrTH1Imput.size()<<std::endl; 
        //Int_t cont=0;
        
        for (auto it:arrTH1Imput){
            TH1D * TH1CopyInt=(TH1D *) it->Clone();
            TH1CopyInt->SetNameTitle("Efficiency Less Than", "Efficiency");
            gStyle->SetOptStat(0);
            //std::cout<<cont<<std::endl;
            //cont++;
            for (int l=Nbines; l>=0; l--){
                Double_t intEff=0;
                Double_t EffErr;
                Double_t IntSignal=it->IntegralAndError(l,Nbines,EffErr);
                Double_t IntSignalTot=it->IntegralAndError(0,Nbines,EffErr);
                //std::cout<<IntSignal<<std::endl;
                //std::cout<<IntSignalTot<<std::endl;
                
                if (IntSignalTot>0){
                    intEff=IntSignal/IntSignalTot;
                    //std::cout<<intEff<<std::endl;
                }else{
                    intEff=0;
                }
                //std::cout<<"Before: "<<TH1CopyInt->GetBinContent(l)<<std::endl;
                TH1CopyInt->SetBinContent(l,intEff);
                //std::cout<<"After: "<<TH1CopyInt->GetBinContent(l)<<std::endl;
                Double_t IntEffErr=0;
                if (!std::isnan(IntEffErr)){
                    TH1CopyInt->SetBinError(l,IntEffErr);
                }else{
                    TH1CopyInt->SetBinError(l,0);
                }  
            }
            //std::cout<<TH1CopyInt->GetBinContent(TH1CopyInt->GetMaximumBin())<<std::endl;
            arrTH1Copy.push_back(TH1CopyInt);
           
        }        
        //std::cout<<arrTH1Copy[0]->GetBinContent(arrTH1Copy[0]->GetMaximumBin())<<std::endl;
        return arrTH1Copy;
    }//End GetEffMoreThan
    
    
    std::vector<TH1D *> GetEffLessThan(void){
        //The last element of the array is the signal
        Int_t Nbines=SigTH1->GetNbinsX();       
        std::vector<TH1D *> arrTH1Copy;
        std::vector<TH1D *> arrTH1Imput;
        arrTH1Imput =ArBackTH1; //Put the back
        arrTH1Imput.insert(arrTH1Imput.begin(), SigTH1);
        //arrTH1Imput.push_back(SigTH1);
        //std::cout<<ArBackTH1.size()<<std::endl; 
        //std::cout<<arrTH1Imput.size()<<std::endl; 
        //Int_t cont=0;
        
        for (auto it:arrTH1Imput){
            TH1D * TH1CopyInt=(TH1D *) it->Clone();
            TH1CopyInt->SetNameTitle("Efficiency Less Than", "Efficiency");
            gStyle->SetOptStat(0);
            //std::cout<<cont<<std::endl;
            //cont++;
            for (int l=0; l<=Nbines; l++){
                Double_t intEff=0;
                Double_t EffErr;
                Double_t IntSignal=it->IntegralAndError(0,l,EffErr);
                Double_t IntSignalTot=it->IntegralAndError(0,Nbines,EffErr);
                //std::cout<<IntSignal<<std::endl;
                //std::cout<<IntSignalTot<<std::endl;
                
                if (IntSignalTot>0){
                    intEff=IntSignal/IntSignalTot;
                    //std::cout<<intEff<<std::endl;
                }else{
                    intEff=0;
                }
                //std::cout<<"Before: "<<TH1CopyInt->GetBinContent(l)<<std::endl;
                TH1CopyInt->SetBinContent(l,intEff);
                //std::cout<<"After: "<<TH1CopyInt->GetBinContent(l)<<std::endl;
                Double_t IntEffErr=0;
                if (!std::isnan(IntEffErr)){
                    TH1CopyInt->SetBinError(l,IntEffErr);
                }else{
                    TH1CopyInt->SetBinError(l,0);
                }  
            }
            //std::cout<<TH1CopyInt->GetBinContent(TH1CopyInt->GetMaximumBin())<<std::endl;
            arrTH1Copy.push_back(TH1CopyInt);
           
        }        
        //std::cout<<arrTH1Copy[0]->GetBinContent(arrTH1Copy[0]->GetMaximumBin())<<std::endl;
        return arrTH1Copy;
    }//End GetEffLessThan
    
    
    
    TH2F *GetSignificance(void){
        
        
        Int_t nbinsx = BackTH1->GetNbinsX()+1, nbinsy =BackTH1->GetNbinsX()+1;//One more bin
        Double_t xlow = (Double_t) BackTH1->GetBinLowEdge(0), xup = (Double_t) (BackTH1->GetBinLowEdge(nbinsx)+BackTH1->GetBinWidth(nbinsx)) ; 
        Double_t ylow = (Double_t) BackTH1->GetBinLowEdge(0), yup = (Double_t) (BackTH1->GetBinLowEdge(nbinsx)+BackTH1->GetBinWidth(nbinsx)) ; 
        TH2F *htext3 = new TH2F("htext3","Significance",nbinsx,xlow,xup,nbinsy,ylow,yup);

        gStyle->SetOptStat(0);
        for (int l=nbinsx; l>=0; l--){
            
            for (int m=0; m<=nbinsx; m++){
                Double_t intZig = 0;
                Double_t sumSB = 0;
                Double_t IntSignal = 0;
                if (m > l){
                    Double_t IntSignalLess=SigTH1->Integral(0,l);//Number of events in the bins less than l 
                    Double_t IntBackgroundLess=BackTH1->Integral(0,l);                               
                    Double_t IntSignalMore=SigTH1->Integral(m,nbinsx);//Number of events in the bins more than l 
                    Double_t IntBackgroundMore=BackTH1->Integral(m,nbinsx);
                    sumSB=IntSignalLess+IntBackgroundLess+IntSignalMore+IntBackgroundMore;
                    IntSignal = IntSignalLess+IntSignalMore;
                            
                }else{ 
                    IntSignal=SigTH1->Integral(m,l);//Number of events between the bins m and l 
                    Double_t IntBackground=BackTH1->Integral(m,l);
                    sumSB=IntSignal+IntBackground;  
                }               
                
                if (sumSB>0){
                    intZig=IntSignal/(TMath::Sqrt(sumSB));
                }else{
                    intZig=0;
                }
                htext3->SetBinContent(l,m,intZig);

            }//End For m (more than)
              
        }//End for l(lessthan)

        //htext3->Draw("COLZ");
        return htext3;
    }//END GetSignificance   
       
};


class CompareTH1D {
    private:
        TH1D *    FirstTH1;
        std::vector<TH1D *> ArrSecondTH1;
   
    public: 
        /*Constructors*/
        CompareTH1D(void):FirstTH1(NULL),ArrSecondTH1(){}
        CompareTH1D(TH1D * _FirstTH1, TH1D * _SecondTH1):FirstTH1(_FirstTH1){
            ArrSecondTH1.push_back(_SecondTH1);
        }
        CompareTH1D(TH1D * _FirstTH1, std::vector<TH1D *> _ArrSecondTH1):FirstTH1(_FirstTH1),ArrSecondTH1(_ArrSecondTH1){}
        CompareTH1D(std::vector<TH1D *> _ArrSecondTH1):FirstTH1(NULL),ArrSecondTH1(_ArrSecondTH1){}

    /*Destructor*/
    ~CompareTH1D(){}
    
    
    /*Set Methods*/
    void SetFirst(TH1D * _FirstTH1){
        FirstTH1= _FirstTH1;
    }//end SetFirst
    
    void AddTH1(TH1D * _SecondTH1){
        ArrSecondTH1.push_back(_SecondTH1);
    }//end SetSecond
     
    /*Get Methods*/
    std::vector<TH1D *> CompareChape(Double_t norm =1.0){
        if (!FirstTH1){
            std::cout<<"Compare Only including First TH1, use SetFirst(TH1D * _FirstTH1)"<<std::endl;
            exit (EXIT_FAILURE);
        }
        std::vector<TH1D *> ArrSecondTH1Normal;
        for (auto it:ArrSecondTH1){
            ArrSecondTH1Normal.push_back((TH1D *)it->Clone());
        }
         TH1D * FirstTH1Normal= (TH1D *)FirstTH1->Clone();
        /*if option "width" is specified, the integral is the sum of
        the bin contents multiplied by the bin width in x.*/
        //FirstTH1Normal->Scale(norm/FirstTH1Normal->Integral("width"));
        FirstTH1Normal->Scale(norm/FirstTH1Normal->Integral());
        for (auto it:ArrSecondTH1Normal){
            it->Scale(norm/it->Integral());
            auto nbines=it->GetNbinsX();
            for (auto i=0; i<nbines; i++){
                Double_t temp =FirstTH1Normal->GetBinContent(i)-it->GetBinContent(i);
                it->SetBinContent(i,temp);
            }
        }
        return  ArrSecondTH1Normal;
    }//END CompareChape
    
    
    std::vector<TH1D *> GetNormalizedTH1(Double_t norm =1.0){
        std::vector<TH1D *> ArrSecondTH1Normal;
        ArrSecondTH1Normal.push_back((TH1D *)FirstTH1->Clone());
        for (auto it:ArrSecondTH1){
            ArrSecondTH1Normal.push_back((TH1D *)it->Clone());
        }
        for (auto it:ArrSecondTH1Normal){
            //it->Scale(norm/it->Integral(), "width");
            it->Scale(norm/it->Integral());
        }
        return  ArrSecondTH1Normal;
    }//END GetNormalizedTH1
        
    //std::vector<TH1D *> CompareContent(void){}
    //   return ArBackTH1;
    
   static TH1D * NormalizedTH1(TH1D * _InTH1,Double_t norm =1.0){
        TH1D * _OutTH1 = (TH1D *)_InTH1->Clone();
        _OutTH1->Scale(norm/_OutTH1->Integral());
        return  _OutTH1;
    }
    
};


    /**/
std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}


/**/
std::string NameMaker(std::string str){
    std::string strout =  str;
    strout = ReplaceAll(strout, std::string(" "), std::string("_"));
    strout = ReplaceAll(strout, std::string("<"), std::string("LT"));
    strout = ReplaceAll(strout, std::string(">"), std::string("MT"));
    strout = ReplaceAll(strout, std::string("#"), std::string(""));
    strout = ReplaceAll(strout, std::string("@"), std::string(""));
    strout = ReplaceAll(strout, std::string("("), std::string(""));
    strout = ReplaceAll(strout, std::string(")"), std::string(""));
    strout = ReplaceAll(strout, std::string("||"), std::string("OR"));
    strout = ReplaceAll(strout, std::string("&&"), std::string("AND"));
    strout = ReplaceAll(strout, std::string("."), std::string(""));
    strout = ReplaceAll(strout, std::string("=="), std::string("EQ"));
    strout = ReplaceAll(strout, std::string("!="), std::string("NEQ"));
    strout = ReplaceAll(strout, std::string("*"), std::string(""));
    strout = ReplaceAll(strout, std::string("%"), std::string(""));
    strout = ReplaceAll(strout, std::string(":"), std::string(""));
    strout = ReplaceAll(strout, std::string("+"), std::string("Plus"));
    strout = ReplaceAll(strout, std::string("-"), std::string("Minus"));
    //strout = ReplaceAll(strout, std::string(""), std::string(""));
    
    return strout;
};

class list_files{
    /********************************************************
    Use
    *********************************************
    list_files lista;
    a=lista.GetFileList(dirname.c_str(),ext.c_str());
    a.size()
    a.at(290)
    *********************************************************/
    private:
    // Vector to hold file names
    std::vector<string> filelist = std::vector<string>();
    //std::vector<string> filelist;
    //const char *_dirname;
   
    
    public:
    /*Constructor*/
    //list_files(const char *_dirname):_dirname(dirname){}
    list_files(void){}
    /*Destructor*/
    ~list_files(){}
    
    void collectAllFiles(const char *_dirname="/", const char *ext=".root")
    {
        TSystemDirectory dir(_dirname, _dirname); 
        TList *files = dir.GetListOfFiles(); 
        std::string dirnamestring= (std::string) _dirname;

        if (files) { 
            TSystemFile *file; 
            TString fname; 
            TIter next(files); 
            while ((file=(TSystemFile*)next())) { 
                fname = file->GetName(); 
                if(!file->IsDirectory() && fname.EndsWith(ext)) { 
                    //cout << fname.Data() << endl; 
                    filelist.push_back(dirnamestring+"/"+fname.Data());
                }
            
                // If this which we found now, is a directory, recursively 
                // call the function again
                if (file->IsDirectory() && strcmp(fname,".") != 0 && strcmp(fname,"..")) { 
                    std::string tempdirnamestring=dirnamestring+"/"+fname.Data();
                    collectAllFiles(tempdirnamestring.c_str(),ext);
                    //collectAllFiles(tempdirnamestring.c_str(),ext,filelist);
                    //std::cout<<fname.Data()<<std::endl;
                }      
            } 
        }
    }

    /*Get Method*/
    std::vector<string> GetFileList(const char *_dirname="/", const char *ext=".root"){
        collectAllFiles(_dirname,ext);
        return filelist;
    }
};




int Lib(void){
return 0;
}
