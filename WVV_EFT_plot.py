#!usr/bin/env python
import os
import re
import glob
import math
import datetime
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse   import OptionParser
from time       import gmtime, strftime
from array import array
import copy
import numpy as np
from ROOT import gROOT, TPaveLabel, TPie, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack, TGraph, TGraphErrors,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, TVectorD, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite,TArrow
import ctypes

ROOT.TH1.SetDefaultSumw2()

parser = OptionParser()
parser.add_option('--channel',    action="store",type="string",dest="channel"    ,default=None)
parser.add_option('--MODE',       action="store",type="string",dest="MODE"       ,default="MC" )
parser.add_option('--REGION',     action="store",type="string",dest="REGION"     ,default="PS" )
parser.add_option('--SFs',        action="store",type="int"   ,dest="SFs"        ,default=0    )
parser.add_option('--piechart',   action="store",type="int"   ,dest="piechart"   ,default=0    )
parser.add_option('--tau',        action="store",type="float" ,dest="tau"        ,default=0.4  )
parser.add_option('--y',          action="store",type="string",dest="y"          ,default="16,17,18")
parser.add_option('--FBT',        action="store",type="int"   ,dest="FBT"        ,default=0    )
parser.add_option('--standalone',        action="store",type="string"   ,dest="standalone"        ,default=None    )
parser.add_option('--Temp',        action="store",type="string"   ,dest="Temp"        ,default=None    )
(options, args) = parser.parse_args()

# ====================== for Limit =====================
# ====================== for Limit =====================
# ====================== for Limit =====================

def UnderOverFlow1D(h):
    Bins=h.GetNbinsX();
    h.SetBinContent( 1,  h.GetBinContent(1)+h.GetBinContent(0) );  h.SetBinError(   1,  math.sqrt( h.GetBinError(1)*h.GetBinError(1) + h.GetBinError(0)*h.GetBinError(0)) );
    h.SetBinContent( Bins,  h.GetBinContent(Bins)+h.GetBinContent(Bins+1) );  h.SetBinError(   Bins,  math.sqrt( h.GetBinError(Bins)*h.GetBinError(Bins) + h.GetBinError(Bins+1)*h.GetBinError(Bins+1)) );
    return h;

def Integerization(h):
    Bins=h.GetNbinsX();
    for i in range(1,Bins+1):
        if (h.GetBinContent(i)-int(h.GetBinContent(i)))>0.5 : value=int(h.GetBinContent(i))+1;
        else : value=int(h.GetBinContent(i));
        h.SetBinContent( i, value          );
        h.SetBinError(   i, math.sqrt(value) );
    return h;

def OptimalCut(B,S):
    Bins = B.GetNbinsX(); B0=B.Integral(); S0=S.Integral(); SigMax=0; CutLocation_L=0; 
    LeftCutBin = 0; RightCutBin = 0; 
    CutLocation_R=0; BKGRejection=0; Sig_Eff=0;
    for RightEnd in range(1,Bins+1):
        for LeftEnd in range(1,RightEnd+1):            #print LeftEnd,RightEnd;
            sig = S.Integral( LeftEnd , RightEnd )/((B.Integral( LeftEnd , RightEnd )+1)**0.5);
            if sig>SigMax: 
                SigMax=sig; LeftCutBin=LeftEnd; RightCutBin=RightEnd; 
                BKGrjc  = float(round(100*(1-B.Integral(LeftEnd,RightEnd)/(B0+0.0001)),1) );
                Sig_Eff = float(round(100*S.Integral(LeftEnd,RightEnd)/(S0+0.0001)    ,1) );
                SigMax  = float(round(SigMax,1));
    ResultList = [ LeftCutBin , RightCutBin , S.GetBinLowEdge(LeftCutBin) , S.GetBinLowEdge(RightCutBin+1), str(BKGrjc), str(Sig_Eff), str(SigMax) ];   #print ResultList;
    return ResultList;  

def RationUnc(h_data,h_TotalMC,h_Ratio,MaxY):
    for i in range(1,h_Ratio.GetNbinsX()+1,1):
        D  = h_data.GetBinContent(i);    eD = h_data.GetBinError(i);
        if D==0: eD=0.92;
        B  = h_TotalMC.GetBinContent(i); eB = h_TotalMC.GetBinError(i);
        if B<0.1 and eB>=B : eB=0.92; Err= 0.;
        if B!=0.        :Err=TMath.Sqrt( (eD*eD)/(B*B)  +(D*D*eB*eB)/(B*B*B*B)     ); h_Ratio.SetBinContent(i, D/B   );  h_Ratio.SetBinError(i, Err); #print i,")",h_Ratio.GetNbinsX()+1,")   data:",D," pm ",eD,"     Bkg:",B," pm ",eB,"   R:",D/B," pm ", Err
        if B==0.        :Err=TMath.Sqrt( (eD*eD)/(eB*eB)+(D*D*eB*eB)/(eB*eB*eB*eB) ); h_Ratio.SetBinContent(i, D/0.92);  h_Ratio.SetBinError(i, Err);
        if D==0 and B==0:                                                             h_Ratio.SetBinContent(i, -1);      h_Ratio.SetBinError(i, 0  );
        if h_Ratio.GetBinContent(i)>MaxY:h_Ratio.SetBinContent(i,MaxY); ### To visualise the points above axis... #h_Ratio.Fit("pol1");
    return h_Ratio;

def OptimalCut(B,S):
    Bins = B.GetNbinsX(); B0=B.Integral(); S0=S.Integral(); SigMax=0; InitialSig=S0/((B0+1)**0.5);
    for RightEnd in range(Bins,0,-1):
        for LeftEnd in range(1,RightEnd+1):            #print LeftEnd,RightEnd;
            sig = S.Integral( LeftEnd , RightEnd )/((B.Integral( LeftEnd , RightEnd )+1)**0.5);
            if sig>SigMax: 
                SigMax=sig; LeftCutBin=LeftEnd; RightCutBin=RightEnd; 
                BKGrjc  = int(round(100*(1-B.Integral(LeftCutBin,RightCutBin)/(B0+0.000001))));
                Sig_Eff = int(round(100*   S.Integral(LeftCutBin,RightCutBin)/(S0+0.000001)) );
                SigMax_Print=int(round(100*(SigMax-InitialSig)/InitialSig));
    ResultList = [ LeftCutBin , RightCutBin , S.GetBinLowEdge(LeftCutBin) , S.GetBinLowEdge(RightCutBin+1), str(BKGrjc), str(Sig_Eff), str(SigMax_Print) ];   #print ResultList;
    return ResultList;  

class Histogram_Remap():
    def __init__(self,outFile):
        self.outFile = outFile

        self.ROOTFiles = []
        self.Histogram     = {}
        self.Histogram_out = {}

    def ReadHistogram(self,InFile):
        fn = ROOT.TFile.Open(InFile)
        for e in fn.GetListOfKeys() : 
            if e.GetClassName() == "TDirectoryFile" :
                for ih in fn.Get(e.GetName()).GetListOfKeys() : 
                    h = fn.Get(e.GetName()).Get(ih.GetName())
                    h.SetDirectory(0)
                    self.Histogram[ih.GetName()] = h
            else :
                h = fn.Get(e.GetName())
                h.SetDirectory(0)
                self.Histogram[e.GetName()] = h
        fn.Close()

    def ReadFiles(self,):
        for irootfile in self.ROOTFiles :
            self.ReadHistogram(irootfile)

    def Remap(self,keep,maps):
        for name in keep :
            self.Histogram_out[name] = self.Histogram[name]
        for name in maps :
            self.Histogram_out[name] = self.Histogram[maps[name][0]].Clone(name)
            for ihname in maps[name][1:] :
                self.Histogram_out.Add(self.Histogram[ihname])

    def Write(self,):
        fn = ROOT.TFile.Open(self.outFile,"recreate")
        for name in self.Histogram_out :
            Histogram = self.Histogram_out[name]
            Histogram.SetLineStyle(3); Histogram.SetMarkerStyle(0); Histogram.SetLineWidth(10); Histogram.SetLineColor(1); # view in cernbox
            Histogram.BufferEmpty(-1) 
            Histogram.GetYaxis().SetRangeUser( Histogram.GetBinContent(Histogram.GetMinimumBin())*0.5 , Histogram.GetBinContent(Histogram.GetMaximumBin())*1.3 )
            Histogram.Write()
        fn.Close()

    def Create(self,keep = [],maps = {}):
        self.ReadFiles()
        self.Remap(keep,maps)
        self.Write()

class LIMIT():
    def __init__(self,WorkPath,ViewPath,WorKFolder, LimitName= "gKK", ):
        self.WorkPath   = WorkPath
        self.ViewPath   = ViewPath
        self.WorKFolder = WorKFolder # this should be unique for every limit calculation
        # self.LimitName  = LimitName

        self.ID         = self.WorKFolder.replace("/","_") # this should be unique for every limit calculation

        # create during running
        self.processes = {}
        self.observation       = {}
        self.Rate              = {}
        self.SystematicSources = []
        self.Systematics       = {}

        pass 

    def Init_ROOTFileName(self,ROOTFile):
        self.ROOTFile = ROOTFile

    def Init_CardName(self,CardName = None):
        if CardName :
            self.Card = "%s/%s"%(self.WorkPath,CardName)
        else :
            self.Card = "%s/datacard.txt"%( self.WorkPath )

    def CreateDir(self):
        self.WorkPath  = "%s/limit_%s/"%(self.WorkPath,self.WorKFolder)
        self.ViewPath  = "%s/limit_%s/"%(self.ViewPath,self.WorKFolder) 
        for ipath in [self.WorkPath,self.ViewPath]:
            if os.path.isdir(ipath) : os.system("rm -r %s"%ipath)
            os.makedirs( ipath )

    def Init(self,ROOTFile,CardName=None,Run = False,**kwargs):
        self.Limit_RunCommands = Run 
        self.ParameterName     = kwargs.get("ParameterName","FT0")
        self.Limit_npoints     = kwargs.get("Limit_npoints","10000")
        self.PlotPath          = kwargs.get("PlotPath","")
        self.parameter_range   = kwargs.get("parameter_range","-30000,30000")

        self.CreateDir()
        self.Init_ROOTFileName(ROOTFile)
        self.Init_CardName(CardName)

    def SetHistName(self, ProcessHistName, SYSTEMATICHistName):
        self.ProcessHistName    = ProcessHistName
        self.SYSTEMATICHistName = SYSTEMATICHistName

    def ConvertRootfile(self,namemap,is2D = False) : 
        # namemap should be like : {"rest":["STOP","Rest",WJets], }
        inf  = ROOT.TFile(self.ROOTFile)
        h_list = [] 
        for iout in namemap:
            # check do we have all the input histogram
            IsInputs = True
            for i in namemap[iout] :
                if i not in [e.GetName() for e in inf.GetListOfKeys()] :
                    IsInputs = False
            if not IsInputs : continue
            hout = inf.Get(namemap[iout][0])
            for iin in namemap[iout][1:] :
                hout.Add(inf.Get(iin))
            if is2D : hout = self.h2D_To_h1D(hout)
            h_list.append( hout.Clone(iout) )
        OutFile = "%s/limit.root"%(self.WorkPath)
        outf = ROOT.TFile( OutFile , "recreate")
        for ih in h_list : 
            ih.SetLineStyle(3); ih.SetMarkerStyle(0); ih.SetLineWidth(10); ih.SetLineColor(1); # view in cernbox
            ih.BufferEmpty(-1) 
            ih.GetYaxis().SetRangeUser( ih.GetBinContent(ih.GetMinimumBin())*0.5 , ih.GetBinContent(ih.GetMaximumBin())*1.3 )
            ih.Write()
        outf.Close()
        inf.Close()
        self.ROOTFile = OutFile

    def updateIndex(self):
        # update the index, index is alphabetic
        for iBin in self.processes:
            processes = list(self.processes[iBin].keys())
            processes.sort()
            processes_nosignal = [ i for i in processes if not self.processes[iBin][i]["IsSignal"] ]
            for iprocess in processes:
                if self.processes[iBin][iprocess]["IsSignal"] : self.processes[iBin][iprocess]["Index"] = 0
            for index,iprocess in enumerate(processes)        : self.processes[iBin][iprocess]["Index"] = index

    def updateObservation(self):
        # update observation
        for iBin in self.processes : self.observation[iBin] = self.observation.get(iBin,-1)

    def updateRate(self):
        # update rate 
        for iBin in self.processes : 
            for iprocess in self.processes[iBin]:
                self.Rate[(iBin,iprocess)] = self.Rate.get((iBin,iprocess),-1)

    def updateUN(self):
        # update UN 
        for iUN in self.Systematics :
            for iBin in self.processes : 
                for iprocess in self.processes[iBin]:
                    self.Systematics[iUN][(iBin,iprocess)] = self.Systematics[iUN].get((iBin,iprocess),None)

    def update(self):
        self.updateIndex()
        self.updateObservation()
        self.updateRate()
        self.updateUN()
        
    def AddProcess(self, Bin, process, observation = None, rate = None, UN = [], IsSignal = False):
        # UN should be like: [[ "lumi_pu", "lnN" , 1.027 ],[ "norm_qcd",  "lnN" , 1.3],] 
        if Bin not in self.processes : 
            if observation : self.observation[ibin] = observation
            self.processes[Bin] = {}
        self.processes[Bin][process] = {}
        self.processes[Bin][process]["IsSignal"] = IsSignal
        for iUN in UN:
            self.Systematics[(iUN[0],iUN[1])]                = self.Systematics.get((iUN[0],iUN[1]),{})
            self.Systematics[(iUN[0],iUN[1])][(Bin,process)] = iUN[2]
        if rate : self.Rate[Bin][process]          = rate
        self.update()
    
    def Card_Part1(self,):
        self.TempleteCard_Part1 = '''
imax    {nbins} number of bins
jmax    {nprocesses} number of processes minus 1
kmax    * number of nuisance parameters
--------------------------------------------------------------------------------
        '''
        nbins      = str(len(self.processes.keys()))
        nprocesses = str(len(self.processes[list(self.processes.keys())[0]].keys())-1)
        return self.TempleteCard_Part1.format( nbins = nbins , nprocesses = nprocesses )

    def Card_Part2(self,):
        # histogram : R{REGION}/signal_region/had_SR{REGION}_{variable}_$PROCESS
        # histogram_with_systematics : R{REGION}/signal_region/had_SR{REGION}_{variable}_$PROCESS_$SYSTEMATIC
        self.TempleteCard_Part2 =  '''
# shapes process channel file histogram [histogram_with_systematics]
shapes * * {ROOTFile} {ProcessHistName} {SYSTEMATICHistName}
--------------------------------------------------------------------------------
        '''
        return self.TempleteCard_Part2.format( ROOTFile = self.ROOTFile , ProcessHistName = self.ProcessHistName , SYSTEMATICHistName = self.SYSTEMATICHistName )

    def Card_Part3(self,):
        self.TempleteCard_Part3 =  '''
bin          {bins}
observation      {observations}
--------------------------------------------------------------------------------
        '''
        bins_ = [] ; observations_ = [] ;
        for iBin in self.processes :
            bins_.append(iBin)
            observations_.append(self.observation.get(iBin,-1))
        bins_ = [ i.ljust(15) for i in bins_ ]
        observations_ = [ str(i).ljust(15) for i in observations_ ]
        bins         = "".join(bins_)
        observations = "".join(observations_)
        return self.TempleteCard_Part3.format( bins = bins, observations = observations )

    def Card_Part4(self,):
        self.TempleteCard_Part4 =  '''
{bin_}                          
{process}                               
{process_Index}                                
{rate}                               
----------------------------------------------------------------------------------------------
# Systematics
'''
        bin_Name_ = [] ; process_Name_ = [] ; process_Index_ = [] ; rate_ = []
        Systematics_ = dict( [ (iUN,[]) for iUN in self.Systematics ] )
        for iBin in self.processes :
            for iprocess in self.processes[iBin] :
                bin_Name_.append(iBin) 
                process_Name_.append(iprocess)
                process_Index_.append( str(self.processes[iBin][iprocess]["Index"]) )
                rate_.append( str(self.Rate[(iBin,iprocess)]) )
                for iUN in self.Systematics:
                    if self.Systematics[iUN][(iBin,iprocess)]:
                        Systematics_[iUN].append(str(self.Systematics[iUN][(iBin,iprocess)]))
                    else:
                        Systematics_[iUN].append("-")
        
        bin_          = "bin".ljust(30)     ; 
        for i in bin_Name_     : 
            bin_          += i.ljust(15)
        
        process       = "process".ljust(30) ; 
        for i in process_Name_ : 
            process       += i.ljust(15)
        
        process_Index = "process".ljust(30) ; 
        for i in process_Index_: 
            process_Index += i.ljust(15)
        
        rate          = "rate".ljust(30)    ; 
        for i in rate_         : 
            rate          += i.ljust(15)
        
        Card = self.TempleteCard_Part4.format( bin_ = bin_ , process = process, process_Index = process_Index, rate = rate )
        for iUN in Systematics_:
            iUN_ = iUN[0].ljust(20) + iUN[1].ljust(10) ; 
            for i in Systematics_[iUN] : iUN_ += i.ljust(15)
            Card += "%s \n"%(iUN_)

        return Card

    def Prepare_DataCard(self):
        # currently only support 1 bin
        print self.Card
        with open(self.Card,"w") as f:
            f.write(self.Card_Part1())
            f.write(self.Card_Part2())
            f.write(self.Card_Part3())
            f.write(self.Card_Part4())

    def READLimitLOG(self, write = True):
        self.limit_SimpleLog = "%s/limitSimple.Log"%(self.WorkPath)
        upper_Limit = {}
        with open(self.limit_Log,"r") as f:
            Lines = f.readlines()
            for i,Line in enumerate(Lines):
                if "-- AsymptoticLimits ( CLs ) --" in Line:
                    upper_Limit_ = [j.replace("Expected ","").replace("\n","") for j in Lines[i+2:i+7]]
                    for iline in upper_Limit_:
                        upper_Limit[ iline.split(": r < ")[0].rstrip() ] = float(iline.split(": r < ")[-1].strip())
        if write:
            with open(self.limit_SimpleLog,"w") as f:
                f.write(str(upper_Limit))
        return upper_Limit

    def CopyToViewPath(self):
        tocopy = [self.Card,self.ROOTFile]
        for icopy in tocopy :
            print "cp -r %s %s "%(icopy,self.ViewPath)
            os.system("cp -r %s %s "%(icopy,self.ViewPath))

    def Commands(self):
        self.Commands_Scripts = "%s/limit.sh"%(self.WorkPath)
        self.limit_Log        = "%s/limit.Log"%(self.WorkPath)
        self.limit_debug      = "%s/limit.debug"%(self.WorkPath)
        cmd  = '''
cd {WorkPath} ; eval `scram runtime -sh`
/eos/user/q/qiguo/limit/aQGC/v1/21_5_27/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/scripts/text2workspace.py {Card} -P HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingEFTNegative:analiticAnomalousCouplingEFTNegative --X-allow-no-signal --PO eftOperators=cG -o {workspaceroot} > spacew.log
combine -M MultiDimFit {workspaceroot} --algo=grid --points {npoins} -m 125 -t -1 --redefineSignalPOIs k_cG --freezeParameters r --setParameters r=1 --setParameterRanges k_cG={parameter_range} --verbose -1 >{LOG} 2>{debug}
root -l higgsCombineTest.MultiDimFit.mH125.root  higgsCombineTest.MultiDimFit.mH125.root   /eos/user/q/qiguo/limit/aQGC/v1/21_5_27/CMSSW_10_2_13/src/HiggsAnalysis/AnalyticAnomalousCoupling/test/draw.cxx\(\\"k_cG\\",\\"expect\\",\\"{PlotPath}\\",\\"{ParameterName}.png\\",\\"{ParameterName}.txt\\",\\"{ParameterName}\\",\\"{ParameterName}\\"\)
        '''.format(
            WorkPath        = self.WorkPath,
            Card            = self.Card,
            workspaceroot   = "%s/workspace.root"%(self.WorkPath),
            LOG             = self.limit_Log,
            debug           = self.limit_debug,
            ParameterName   = self.ParameterName,
            npoins          = self.Limit_npoints,
            PlotPath        = self.PlotPath,
            parameter_range = self.parameter_range,
        )
        with open(self.Commands_Scripts,"w") as f:
            f.write(cmd)
        if self.Limit_RunCommands :
            os.system("sh %s"%(self.Commands_Scripts))
            return_list = self.READLimitLOG()
            self.CopyToViewPath()
        else : 
            return_list = []
        return return_list

class ANALYSIS:
    def __init__(self, channel , fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1", input_workspace=None):
        
        self.setTDRStyle();
        
        self.channel    =channel;

        self.color_palet={'data':1, 'QCD':2,  'Rest':62,  'VV':62, 'STop':8, 'TTbar':80, 'ZJets':6, 'WJets':90, 'Signal':1, 'Uncertainty':1, }

        self.Debug_keep_Nentries = None

        self.MODE = options.MODE

        self.AddTxt = None
        self.REGION_text = None

        
        # ================== signal scale ================== 
        # ================== signal scale ================== 
        self.Scale_Signal_Auto = False
        self.Signal_Scale1 = 10000
        self.Signal_Scale2 = 10000
        if options.MODE == "DECO":
            self.Signal_Scale1 = 1
            self.Signal_Scale2 = 1
        # ================== signal scale ================== 
        # ================== signal scale ================== 
        
        self.Optimal = True

        # ================ for decomposite ==================
        # ================ for decomposite ==================
        self.DECO_OnlySignal = False
        self.Signal_To_Draw = None
        self.DECO_py = "config/DECO.py"
        self.signal1_File = "_out_Res1ToRes2GluToGluVV_M1-2500_R-250.root" ; self.signal1_DECO_Label = "(2500, 250) GeV"
        self.signal2_File = "_out_Res1ToRes2GluToGluVV_M1-3000_R-1500.root"
        self.BKG_to_Draw = None
        # self.DECO_Matching = "deep-Wa_flip"
        self.DECO_Matching = "deep-W_a"
        self.DECO_Matching = "deep-W_c"

        self.Unit_Norm = False
        # ================ for decomposite ==================
        # ================ for decomposite ==================
        
        # PlotPath
        self.PlotPath = "/eos/user/q/qiguo/www/gKK/plots/MET_Study/"

        # =============== for Intime Compile =================
        # =============== for Intime Compile =================
        self.Intime_Add_selection_Variables = None
        self.KeepColumn = {"Pt_tag3","Pt_tag1", "MET_et", "Nj8", "weight"}
        self.Plot2DFaster = False
        
        self.MutiThreads = False
        # =============== for Intime Compile =================
        # =============== for Intime Compile =================


        # =============== add Text for significance ================
        # =============== add Text for significance ================
        self.AddText_ToEachBin = False
        self.Textsize = 1.0
        if self.AddText_ToEachBin:
            self.Pad_Split = 0.5
        # =============== add Text for significance ================
        # =============== add Text for significance ================

        # special pad setting
        self.Pad_Split = None

        self.plotsize = None

        # ============== setting for Limit ================
        # ============== setting for Limit ================
        self.Run_limit = False
        self.Signal_limit = "gKK_2500_250"
        self.LIMIT_CODEPATH = "/eos/user/q/qiguo/gKK/limit/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/limit/"
        self.TextColor_limit = 4
        # ============== setting for Limit ================
        # ============== setting for Limit ================

        self.Draw_Data_2D = False

        self.Plotname = None

        if options.MODE in ["MC","MCvsDATA","DECO"] :
            self.DrawData = True
        if options.MODE in ["MC", "DECO"] :
            self.DrawData = False

        # for VVV EFT
        self.Intime_py_list = [
            "weight.py",
            "Events_Level.py",
            "AK8.py",
            "leptonic_W.py",
            "bJet.py",
            "Genparticle_matching.py",
        ]
        self.Intime_py_list = ["config/Intime/"+i for i in self.Intime_py_list]

        self.TreeName = "PKUTree"
        self.MCweight="weight"

        self.BKG_Correct = False
        self.TSFFile = "config/TSF/TF_v1.root"

        # ============= DECO ==============
        self.DECO = ["unmatching","qg","W"]
        self.DECO_Cut = {
            "W" : "1",
            "qg" : "1",
            "unmatching" : "1",
        }
        self.Fill_Style = {
            "W" : 1,
            # "W" : 3001,
            "qg" : 3023,
            "unmatching" : 3010,
        }
        self.Only_Show_Component = True
        self.color_DECO = {
            "W" : 80,
            "qg" : 2,
            "unmatching" : 19,
        }
        self.DECO_Label = {
            "W" : "W",
            "qg" : "qg",
            "unmatching" : "un",
        }
        # ============= DECO ==============

        # ================ for store histogram =================
        # ================ for store histogram =================
        self.Histogram_ToStore = ["h_data","h_QCD","h_WJets","h_TTbar","h_STop","h_WWW","h_WWZ","h_WZZ","h_TotalMC","h_VV", "h_WWW_aQGC_FT1_0p8","h_WWW_aQGC_FT1_0p4","h_WWW_aQGC_FT1_0p2"]
        # ================ for store histogram =================
        # ================ for store histogram =================
        self.init()

    def init(self,):
        self.init_channel()
        self.init_year()
        self.init_region()
        self.init_cut()
        self.init_weight()
        self.init_Intime()
        self.init_padsplit()
        self.init_deco()
        self.init_storeROOTFile()
        self.init_table()
        self.init_InputFile()
        self.init_limit()
        self.init_Compare()

    def init_year(self,):
        self.year = options.y

    def init_channel(self,):
        if not options.channel :
            sys.exit("channel is not given")
        self.channel = options.channel

    def Signal_Name_WeightName(self,signal,iop,igrid):
        Name       = "Limit_%s_%s_%s"%( signal,iop,str(igrid).replace("-","m").replace(".","p") )
        WeightName = "weight_%s_%s_%s"%( signal,iop,str(igrid).replace("-","m").replace(".","p") )
        return Name,WeightName

    def init_Compare(self,):
        self.readHistogram = False
        self.SaveHistograms = []
        self.pad1range = None

        self.TFpath        = "TF/2016/QCD/"
        if not os.path.isdir(self.TFpath) : os.makedirs(self.TFpath)

        self.resolvedCRSR = "center"

        self.HistogramTF = [] ; self.weight_BKGTF = None ; self.weight_DataTF = None ; 

        self.pad2Range = [0,2]

        self.Pad1LogY = False

        self.MarkerSize = {}
        self.LineWidth  = {}

        self.MRBin_ = [50,100,150,200,250,300,350,400]

        self.COM_py = "config/Cuts/com.py"
        self.ComName  = "Center"

    def init_limit(self,):
        self.run_limit       = False
        self.limit_parameter = ["FT0"]
        
        self.signal_File = {
            "WWW" : "WWW_aQGC.root",
            "WWZ" : "WWZ_aQGC.root",
        }
        
        self.Limit_SignalSamples = [
            "WWW","WWZ"
            # "WWW",
            # "WWZ",
        ]
        self.aQGC_Term       = {} 
        self.Limit_Histogram_3Term = []

        self.LIMIT_CODEPATH = "/eos/user/q/qiguo/limit/aQGC/v1/21_5_27/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/Limit/"
        self.LimitFolder_PostFix = "Test"

        self.Limit_namemap_init = {
            "h_BKG"        : ["h_TotalMC"],
            "h_data_obs"   : ["h_data_obs"],
            "h_BKG_StatisticUp"   : ["h_TotalMC_Statistic_up"],
            "h_BKG_StatisticDown" : ["h_TotalMC_Statistic_down"],
        }
        self.Limit_Process = {
            "quad_cG"        : {},
            "sm"             : {},
            "sm_lin_quad_cG" : {},
            "BKG"            : {
                "un" : [
                    ["Statistic","shape",1.],
                ],
            }
        }

        self.Limit_namemap = self.Limit_namemap_init
        self.Limit_ChannelName = "Muon"
        self.Limit_RunCommands = False
        self.Limit_npoints = 10000

        self.Limit_Range = {
            "FS0" : "-3000,3000",
            "FS1" : "-3000,3000",
            "FS2" : "-3000,3000",
            "FM0" : "-30,30",
            "FM1" : "-30,30",
            "FM2" : "-30,30",
            "FM3" : "-30,30",
            "FM4" : "-30,30",
            "FM5" : "-30,30",
            "FM7" : "-30,30",
            "FT0" : "-3,3",
            "FT1" : "-3,3",
            "FT2" : "-3,3",
            "FT3" : "-3,3",
            "FT4" : "-3,3",
            "FT7" : "-3,3",
        }

        self.limit_results = {}

        # store mc's statistic uncertainty as systmatic uncertainty
        self.Limit_BKGStatisticUncertainty = []

        self.histogram_Yields_ratio = None

    def init_InputFile(self,):
        self.InputFileConfig = "InputFile_V3.py"

    def init_region(self,):
        self.REGION = options.REGION

    def init_weight(self,):
        # weight, weight_dnncorr_0
        self.MCweight="weight"

    def init_cut(self,):
        self.init_cut_mu()


    def init_cut_mu(self):
        
        self.PS_Ele = " abs(Lep1fatJet2_LeptonPDGID) == 11 && Lep1fatJet2_LeptonPt > 30 && Lep1fatJet2_FatJet_pt > 200 && Lep1fatJet2_FatJet_pt_2 > 200 && Lep1fatJet2_FatJet_msoftdrop > 40 && Lep1fatJet2_FatJet_msoftdrop_2 > 40 "


    def init_Intime(self,):
        self.Intime_Add_Variables_BKG = None
        self.Intime_Add_Variables_signal = None
        self.Intime_Add_Variables_data = None

    def init_padsplit(self,):
        self.Pad_Split = 0.29
        if self.MODE in ["BKG_DECO_vsDATA","BKG_DECO_compare"]:
            self.plotsize  = (700,700)
        if self.MODE in ["BKG_DECO"]:
            self.plotsize  = (700,565)

    def init_deco(self,):
        self.FillStyle = {}
        self.LineWidth = {}

    def init_storeROOTFile(self,):
        self.Histogram_ToStore        = [] 
        self.Histogram_ToStore_Signal = [] 

    def init_table(self,):
        self.TableRaw = "PS"
        self.text     = "/eos/user/q/qiguo/www/VVV/plots/plots/PS/Test/table/temp/PS_table.txt"

    def Limit_GetBinContent(self,):
        Content          = {}
        Operators = self.limit_parameter + ["SM"]
        iops = [] ; igrids = [] ; ibins = [] ; signals = []
        for iop in Operators:
            for signal in self.Limit_SignalSamples :
                if not iop in self.Parameter_Ranges[signal] : continue
                for igrid in self.Parameter_Ranges[signal][iop] :
                    signal_name,_ = self.Signal_Name_WeightName( signal,iop,igrid )
                    histogram     = self.LimitHistogram[signal_name]
                    for ibin in range(1,histogram.GetNbinsX()+1) :
                        if "Limit_WWZ_SM_0" == signal_name :
                            print "ibin",ibin,histogram.GetBinContent(ibin)
                        iops.append(iop)
                        igrids.append(igrid)
                        ibins.append(ibin)
                        signals.append(signal)
                        Content[(iop,igrid,ibin,signal)] = histogram.GetBinContent(ibin)
        iops    = list(set(iops)) ; 
        igrids  = list(set(igrids)) ; 
        ibins   = list(set(ibins)) ; 
        signals = list(set(signals)) ; 
        self.Limit_BinContent = {}
        for iop in iops :
            for igrid in igrids :
                for ibin in ibins :
                    self.Limit_BinContent[ (iop,igrid,ibin) ] = 0
                    for signal in signals :
                        self.Limit_BinContent[ (iop,igrid,ibin) ] += Content.get((iop,igrid,ibin,signal),0)
        print Content
        print self.Limit_BinContent

    def Get_iopRange(self,iop):
        for signal in self.Parameter_Ranges :
            if iop in self.Parameter_Ranges[signal] :
                return self.Parameter_Ranges[signal][iop]

    def fit_results(self, **kwargs):
        # def fit_results(par_var,ratio,path,filename="fit_parameter.txt"):
        iop   = kwargs.get("iop","")
        ibin  = kwargs.get("ibin","")

        plotpath = "%s/quadratic/"%(self.limit_FilePath)
        if not os.path.isdir( plotpath ):
            os.makedirs( plotpath ) 
        plotname = "%s/%s_%s.png"%(plotpath,iop,ibin)

        ratio = [] ; par_var = []
        for igrid in self.Get_iopRange(iop) :
            value = self.Limit_BinContent[ (iop,igrid,ibin) ]/self.Limit_BinContent[ ("SM",0,ibin) ]
            ratio.append(value)
            par_var.append(igrid)
        x = np.array(par_var).astype(np.double)
        y = np.array(ratio).astype(np.double)
        gr = ROOT.TGraph(len(x), x, y)
        low = min(par_var)
        high = max(par_var)
        fitFunc = ROOT.TF1("fit_result","[0]*(x**2) + [1]*x + 1",low,high) 
        fitFunc.SetLineColor(ROOT.kBlue) 
        r = gr.Fit("fit_result","ESR") 
        gr.SetLineWidth(2) 
        gr.SetLineColor(ROOT.kBlue) ;
        gr.SetMarkerStyle(4) 
        gr.SetTitle("Ratio: [aQGC]/[aQGC->SM]")
        chi2   = r.Chi2() ;
        par0   = fitFunc.GetParameter(0);
        par1   = fitFunc.GetParameter(1);
        err0   = fitFunc.GetParError(0) ;
        err1   = fitFunc.GetParError(1) ;
        self.aQGC_Term[ (iop,ibin) ] = {
            "Quad" : par0,
            "Lin"  : par1,
        }
        # r.Print("V") 
        c1= ROOT.TCanvas("c1","fitFunc",500,500) 
        c1.SetGridx(1) 
        c1.SetGridy(1) 
        gr.Draw("AP")
        c1.SaveAs(plotname)

    def Limit_FitBin(self,):
        histogram = self.LimitHistogram[list(self.LimitHistogram.keys())[0]]
        bins      = range(1,histogram.GetNbinsX()+1)
        for iop in self.limit_parameter:
            if iop == "SM" : continue
            for ibin in bins :
                self.fit_results(iop=iop,ibin=ibin)

    def Limit_Histogram_Signal(self,):
        for index,signal in enumerate(self.Limit_SignalSamples):
            histogram_name,_ = self.Signal_Name_WeightName(signal,"SM",0)
            if index == 0 :
                histogram = self.LimitHistogram[histogram_name].Clone("SM_temp")
            else :
                histogram.Add( self.LimitHistogram[histogram_name] )
        bins      = range(1,histogram.GetNbinsX()+1)
        self.Limit_Histogram_3Term = {}
        for iop in self.limit_parameter:
            if iop == "SM" : continue
            Sm_name        = "h_Sm_%s"%(iop)
            Quad_name      = "h_Quad_%s"%(iop)
            Lin_name       = "h_Lin_%s"%(iop)
            SmLinQuad_name = "h_SmLinQuad_%s"%(iop)
            h_Sm        = histogram.Clone(Sm_name)
            h_Quad      = histogram.Clone(Quad_name)
            h_Lin       = histogram.Clone(Lin_name)
            h_SmLinQuad = histogram.Clone(SmLinQuad_name)
            for ibin in bins :
                Sm        = histogram.GetBinContent(ibin)
                Quad      = Sm*self.aQGC_Term[(iop,ibin)]["Quad"]
                Lin       = Sm*self.aQGC_Term[(iop,ibin)]["Lin"]
                SmLinQuad = Sm*(self.aQGC_Term[(iop,ibin)]["Quad"]+self.aQGC_Term[(iop,ibin)]["Lin"]+1)
                h_Quad.SetBinContent(ibin,Quad) 
                h_Lin.SetBinContent(ibin,Lin) 
                h_SmLinQuad.SetBinContent(ibin,SmLinQuad) 
                self.Limit_Histogram_3Term[Sm_name]        = h_Sm
                self.Limit_Histogram_3Term[Lin_name]       = h_Lin
                self.Limit_Histogram_3Term[Quad_name]      = h_Quad
                self.Limit_Histogram_3Term[SmLinQuad_name] = h_SmLinQuad

    def Limit_Histogram(self,):
        self.Limit_GetBinContent()
        self.Limit_FitBin()
        self.Limit_Histogram_Signal()
        self.Limit_FakeData_Hisgoram()

    def Limit_FakeData_Hisgoram(self,):
        histogram = self.LimitHistogram[self.LimitHistogram.keys()[0]].Clone("h_data_obs")
        FakeBinContent = 1000
        for ibin in range(1,histogram.GetNbinsX()+1) :
            histogram.SetBinContent(ibin,FakeBinContent)
        self.Limit_Histogram_3Term["data_obs"] = histogram

    def Limit_CommandsFolder(self,args):
        templete = "/".join(args) + "/"
        return templete

    def Limit_Run(self,):
        # loop over All the aQGC parameters
        for iop in self.limit_parameter :
            if iop == "SM" : continue
            # create dir for each operator
            self.Limit_Path_iop = "%s/Limit_%s/%s/"%(self.PlotPath,self.LimitFolder_PostFix,iop)
            self.Limit_CreateDir(iop)
            # create histogram needed for parameters
            rootfile = self.Limit_CreateRootfile(iop)
            # create data cards and run scripts
            self.Limit_Datacard_Commands(iop,rootfile)

    def Limit_NameMap(self,iop):
        self.Limit_namemap = self.Limit_namemap_init
        signals = ["Sm","Quad","SmLinQuad"]
        for iname in self.Limit_Histogram_3Term:
            if not iop in iname : continue
            issignal = False 
            for isignal in signals:
                if isignal in iname : issignal = True
            if issignal : 
                if "Sm"        in iname : hname = "h_sm"
                if "Quad"      in iname : hname = "h_quad_cG"
                if "SmLinQuad" in iname : hname = "h_sm_lin_quad_cG"
                self.Limit_namemap[hname] = [iname]

    def Limit_StoreUncertainties_Histogram(self,iop):
        output_ROOTFile = "%s/Uncertainties.root"%(self.Limit_Path_iop)
        outf = ROOT.TFile( output_ROOTFile, "recreate")

        self.histogramROOTfileUntemp = output_ROOTFile

        NSigma = 1

        templete_name_up   = "%s_Statistic_up"
        templete_name_down = "%s_Statistic_down"
        for process in self.Limit_BKGStatisticUncertainty :
            hname               = self.LimitHistogram[process].GetName()
            histogram_name_up   = templete_name_up%( hname )
            histogram_name_down = templete_name_down%( hname )
            histogram_up   = self.LimitHistogram[process].Clone(histogram_name_up)
            histogram_down = self.LimitHistogram[process].Clone(histogram_name_down)
            for i in range(1,1+histogram_up.GetNbinsX()) :
                sigma  = histogram_up.GetBinError(i)
                center = histogram_up.GetBinContent(i)
                up = center + NSigma*sigma ; down = center - NSigma*sigma 
                histogram_up.SetBinContent(   i, up )
                histogram_down.SetBinContent( i, down )
            histogram_up.Write()
            histogram_down.Write()

        outf.Close()

    def Limit_PrepareUncertainty(self,iop):
        self.Limit_StoreUncertainties_Histogram(iop)

    def Limit_CreateDir(self,iop):
        plotpath = self.Limit_Path_iop
        if not os.path.isdir(plotpath) : os.makedirs(plotpath)

    def Limit_CreateRootfile(self,iop,) : 
        self.Limit_PrepareUncertainty(iop)
        self.Limit_NameMap(iop)
        rootfile = "%s/Histogram.root"%(self.Limit_Path_iop)
        Histogram_Remap_ = Histogram_Remap( rootfile )
        Histogram_Remap_.ROOTFiles += [
            self.histogramROOTfile ,
            self.histogramROOTfileUntemp ,
        ]
        Histogram_Remap_.Create(
            maps = self.Limit_namemap
        )
        return rootfile
        
    def Limit_Datacard_Commands(self,iop,rootfile):
        Limit_Folder = self.Limit_CommandsFolder([self.LimitFolder_PostFix,])
        LIMIT_ = LIMIT(self.LIMIT_CODEPATH,self.PlotPath,Limit_Folder)
        plotpath = self.Limit_Path_iop
        
        LIMIT_.Init( 
            rootfile ,
            Run = self.Limit_RunCommands,
            ParameterName = iop,
            Limit_npoints = self.Limit_npoints,
            PlotPath = plotpath,
            parameter_range = self.Limit_Range[iop],
        )
        for process in self.Limit_Process :
            if "data" in process : continue 
            un = self.Limit_Process[process].get("un",[])
            LIMIT_.AddProcess(self.Limit_ChannelName,process,UN = un)
        LIMIT_.SetHistName( "h_$PROCESS", "h_$PROCESS_$SYSTEMATIC")
        LIMIT_.Prepare_DataCard()
        LIMIT_.Commands()
        LIMIT_.CopyToViewPath()
        if self.Limit_RunCommands : 
            Limit_Log = os.path.normpath("%s/%s.txt"%(plotpath,iop))
            with open(Limit_Log,"r") as f : 
                self.limit_results.update(eval(f.read()))

    def View_Selection_With_Plots(self,Text_Name1,variable,selection,comment,View_Selection_Only = True):
        Swith_Char = {};Swith_Char["<"] = "&lt;"; Swith_Char[">"] = "&gt;"; Variable_Style = "<aaa>"; selection_Style = "<aa>"; comment_Style = "<aaaa>"; Switch_Line_Number_variable = 50; Switch_Line_Number_selection = 100; Switch_Line_Number_comment = 80; Text_Content = "";

        for i in Swith_Char:
            variable.replace(i,Swith_Char[i])
            selection.replace(i,Swith_Char[i])

        variable_list = list(variable)
        selection_list = list(selection)
        comment_list = list(comment)
        for i in range(1,int(float(len(variable_list)/Switch_Line_Number_variable))+1):
            variable_list.insert(i*Switch_Line_Number_variable,"<br>")
        for i in range(1,int(float(len(selection_list)/Switch_Line_Number_selection))+1):
            selection_list.insert(i*Switch_Line_Number_selection,"<br>")
        for i in range(1,int(float(len(comment)/Switch_Line_Number_comment))+1):
            selection_list.insert(i*Switch_Line_Number_comment,"<br>")
        variable = ''.join(variable_list)
        selection = ''.join(selection_list)
        comment = ''.join(comment_list)

        Text_Content += Variable_Style + variable + Variable_Style + "<br>"
        Text_Content += comment_Style + comment + comment_Style + "<br>"
        Text_Content += selection_Style + selection + selection_Style

        if View_Selection_Only:
            Text_Content = selection_Style + selection + selection_Style

        with open(Text_Name1,"w") as f:
            f.write(Text_Content)


    #================ SETTINGS FOR Canvas/pads/histos and more ==================
    def setTDRStyle(self):
        self.tdrStyle =TStyle("tdrStyle","Style for P-TDR");  self.tdrStyle.SetCanvasBorderMode(0);        self.tdrStyle.SetCanvasColor(kWhite);        self.tdrStyle.SetCanvasDefH(700);        self.tdrStyle.SetCanvasDefW(700);        self.tdrStyle.SetCanvasDefX(0);          self.tdrStyle.SetCanvasDefY(0);
        self.tdrStyle.SetPadBorderMode(0);             self.tdrStyle.SetPadColor(kWhite);        self.tdrStyle.SetPadGridX(False);        self.tdrStyle.SetPadGridY(False);        self.tdrStyle.SetGridColor(0);        self.tdrStyle.SetGridStyle(3);        self.tdrStyle.SetGridWidth(1);      
        self.tdrStyle.SetFrameBorderMode(0);        self.tdrStyle.SetFrameBorderSize(1);        self.tdrStyle.SetFrameFillColor(0);        self.tdrStyle.SetFrameFillStyle(0);        self.tdrStyle.SetFrameLineColor(1);        self.tdrStyle.SetFrameLineStyle(1);        self.tdrStyle.SetFrameLineWidth(1);
        self.tdrStyle.SetHistLineColor(1);        self.tdrStyle.SetHistLineStyle(0);        self.tdrStyle.SetHistLineWidth(1);        self.tdrStyle.SetEndErrorSize(2);              self.tdrStyle.SetMarkerStyle(20);      self.tdrStyle.SetErrorX(0.);
        self.tdrStyle.SetOptFit(1);        self.tdrStyle.SetFitFormat("5.4g");        self.tdrStyle.SetFuncColor(2);        self.tdrStyle.SetFuncStyle(1);        self.tdrStyle.SetFuncWidth(1);      self.tdrStyle.SetOptDate(0);      
        self.tdrStyle.SetOptFile(0); self.tdrStyle.SetOptStat(0); self.tdrStyle.SetStatColor(kWhite); self.tdrStyle.SetStatFont(42); self.tdrStyle.SetStatFontSize(0.025); self.tdrStyle.SetStatTextColor(1); self.tdrStyle.SetStatFormat("6.4g"); self.tdrStyle.SetStatBorderSize(1); self.tdrStyle.SetStatH(0.1); self.tdrStyle.SetStatW(0.15);
        self.tdrStyle.SetPadTopMargin(0.05);        self.tdrStyle.SetPadBottomMargin(0.13);        self.tdrStyle.SetPadLeftMargin(0.18);        self.tdrStyle.SetPadRightMargin(0.06);      
        self.tdrStyle.SetOptTitle(0);        self.tdrStyle.SetTitleFont(42);        self.tdrStyle.SetTitleColor(1);        self.tdrStyle.SetTitleTextColor(1);        self.tdrStyle.SetTitleFillColor(10);        self.tdrStyle.SetTitleFontSize(0.05);
        self.tdrStyle.SetTitleColor(1, "XYZ");        self.tdrStyle.SetTitleFont(42, "XYZ");        self.tdrStyle.SetTitleSize(0.06, "XYZ");  
        self.tdrStyle.SetTitleXOffset(0.8);        self.tdrStyle.SetTitleYOffset(0.8);      
        self.tdrStyle.SetLabelColor(1, "XYZ");        self.tdrStyle.SetLabelFont(42, "XYZ");        self.tdrStyle.SetLabelOffset(0.007, "XYZ");        self.tdrStyle.SetLabelSize(0.04, "XYZ");
        self.tdrStyle.SetAxisColor(1, "XYZ");        self.tdrStyle.SetStripDecimals(kTRUE);        self.tdrStyle.SetTickLength(0.03, "XYZ");        self.tdrStyle.SetNdivisions(510, "XYZ");        self.tdrStyle.SetPadTickX(1);       self.tdrStyle.SetPadTickY(1);      
        self.tdrStyle.SetOptLogx(0); self.tdrStyle.SetOptLogy(0); self.tdrStyle.SetOptLogz(0);
        self.tdrStyle.SetPaperSize(20.,20.); self.tdrStyle.cd();


    def DefineSelection_0lep(self):#==========[ 0lep REGIONS & SELECTION ]===========================================
        REGION=options.REGION; MODE=options.MODE;
        PS = "1"
        PS3 = "1"
        
        if self.MODE in ["MC","MCvsDATA","DECO","BKG_DECO_compare","BKG_DECO"]: self.Make_Controlplots_for_0lep(eval(REGION),"","");

        # add by Qilong start
        if self.MODE == "2DPlot": 
            if options.standalone:
                with open(options.standalone,"r") as f:
                    exec(f.read())


    def Make_Controlplots_for_0lep(self,selection,selection2,tag,CR=0):
        REGION=options.REGION; Nj=234; MODE=options.MODE; logy=0; tag="";
        if self.MODE in ["MC","MCvsDATA","DECO"]: logy=0;
        if options.standalone:
            with open(options.standalone,"r") as f:
                exec(f.read())
   
    def AddIntimeVar(self,TreeName):
        Tree_RDF = self.decodf_tree[TreeName]
        Type = self.DECO_FILEType.get(TreeName,"BKG")
        
        for Intime_py in self.Intime_py_list:
            with open(Intime_py,"r") as f:
                exec(f.read())

        if Type == "BKG" :
            if self.Intime_Add_Variables_BKG:
                for iselection in self.Intime_Add_Variables_BKG:
                    if iselection not in [str(i) for i in Tree_RDF.GetColumnNames()]:
                        Tree_RDF = Tree_RDF.Define( iselection ,eval(iselection) )

        if Type == "data" :
            if self.Intime_Add_Variables_data:
                for iselection in self.Intime_Add_Variables_data:
                    if iselection not in [str(i) for i in Tree_RDF.GetColumnNames()]:
                        Tree_RDF = Tree_RDF.Define( iselection ,eval(iselection) )
        self.decodf_tree[TreeName] = Tree_RDF
        
    

    def histogramf(self,variable,variable_bin,name,xtitle,ytitle,nbin,min,max,Intime,Tree_RDF,cut,weight) :
        print name,weight
        histogram_prefix = "h_"
        histogram_name = histogram_prefix+name
        if variable_bin :
            histo_bin = array('d',variable_bin)
            h    = TH1D( histogram_name ,histogram_name  +";%s;%s"%(xtitle,ytitle),(len(histo_bin)-1),histo_bin); 
            Bin = ("Plot",";%s;%s"%(xtitle,ytitle),(len(histo_bin)-1),histo_bin)
        else:
            h    = TH1D( histogram_name ,histogram_name +";%s;%s"%(xtitle,ytitle),nbin,min,max);
            Bin = ("Plot",";%s;%s"%(xtitle,ytitle),nbin,min,max)
        h.Sumw2();

        if getattr(self,"corr_weights",None) :
            if weight in self.corr_weights :
                Tree_RDF = Tree_RDF.Define( weight ,self.corr_weights[weight] )

        # weight = "weight"
        if Intime:
            if weight : 
                h = Tree_RDF.Filter(cut).Histo1D(Bin, variable, weight)
            else :
                h = Tree_RDF.Filter(cut).Histo1D(Bin, variable)
            outname = "o_%s"%(name)
            self.histogram_booked[outname] = {
                "h" : h,
                "histogram_name" : histogram_name,
            }
            
        
        h = UnderOverFlow1D(h)
        if "data" in name:
            h.SetBinErrorOption(TH1D.kPoisson);
        
        histogram_on1D = [i["name"] for i in self.histogram_to_draw_pad2]
        histogram_on2D = [i["name"] for i in self.histogram_to_draw_pad1]
        self.color_palet = self.color_palet_Com 
        print self.color_palet
        if name in (histogram_on1D+histogram_on2D):
            h.SetLineColor(self.color_palet[name]); 
            h.SetFillColor(self.color_palet[name]);
            h.SetMarkerColor(self.color_palet[name]);
        h.SetFillStyle(self.FillStyle.get(name,0));
        h.SetLineWidth(self.LineWidth.get(name,6));

        if not Intime : self.histogram[name] = h

        if not self.h_TotalMC:
            if variable_bin:
                self.h_TotalMC = TH1D( "h_TotalMC" ,"h_TotalMC"  +";%s;%s"%(xtitle,ytitle),(len(histo_bin)-1),histo_bin); 
                Bin = ("Plot",";%s;%s"%(xtitle,ytitle),(len(histo_bin)-1),histo_bin)
            else:
                self.h_TotalMC = TH1D( "h_TotalMC" ,"h_TotalMC" +";%s;%s"%(xtitle,ytitle),nbin,min,max);
            self.h_TotalMC.Sumw2();
            
    def Ratio_Histogram(self,hlist,name,**kwargs):
        minY = kwargs.get("minY",0)
        maxY = kwargs.get("maxY",2)
        color = kwargs.get("color",1)
        markerstyle = kwargs.get("markerstyle",8)
        Title = kwargs.get("Title",None)
        h1 = hlist[0] ; h2 = hlist[1] ; 

        h_Ratio = h1.Clone("h_"+name);
        h_Ratio.Divide( h2 );
        if Title:
            h_Ratio.GetYaxis().SetTitle(Title);  
        h_Ratio.SetLineColor(color);
        h_Ratio.SetMarkerStyle(markerstyle);

        h_Ratio.SetLineWidth(2);
        h_Ratio.SetMarkerSize(0.7);

        h_Ratio.GetYaxis().SetNdivisions(504,0);
        h_Ratio.GetYaxis().SetTitleOffset(0.35);  
        h_Ratio.GetYaxis().SetTitleSize(0.13);  
        h_Ratio.GetYaxis().SetTitleSize(0.13);  
        h_Ratio.GetYaxis().SetLabelSize(0.11); 
        h_Ratio.GetXaxis().SetLabelSize(0.1); 
        h_Ratio.GetXaxis().SetTitleOffset(0.7); 
        h_Ratio.GetXaxis().SetTitleSize(0.14);

        h_Ratio = RationUnc(h1,h2,h_Ratio,maxY);
        h_Ratio.GetYaxis().SetRangeUser(minY,maxY);  

        self.histogram[name] = h_Ratio

    def f_Four_arithmetic_histogram(self,**kwargs):
        hlist   = kwargs.get("hlist"  ,None)
        fomula_ = kwargs.get("fomula_",None)
        name    = kwargs.get("name"   ,None)
        print hlist
        print self.histogram.keys()
        print self.DECO_hs.keys()
        for i,ih in enumerate(hlist):
            ihname = ih if ih in self.histogram else ("o_"+ih)
            hi = self.histogram.get(ihname,None)
            if not hi : 
                print "no",ih
            hi.Integral()
            exec("h%s =  hi"%(str(i+1)))
        self.histogram[name] = eval(fomula_)

    def plotdeco_norm(self,**kwargs):
        name      = kwargs.get("name",None)
        hIntegral = kwargs.get("hIntegral",None)
        Integral = 1.
        for j in [name,"o_"+name]:
            if j in self.histogram:
                iname = j
                break
        iname = kwargs.get("iname", iname)
        self.histogram[name] = self.histogram[iname].Clone("h_"+name)
        if hIntegral :
            hIntegral = hIntegral if hIntegral in self.histogram else "o_"+hIntegral
            scale = 1/(self.histogram[hIntegral].Integral())
            self.histogram[name].Scale(scale)
        else :
            self.histogram[name].Scale(Integral/self.histogram[name].Integral())

    def plotdeco_norm_by_A_divide_B(self,**kwargs):
        name      = kwargs.get("name",None)
        inname      = kwargs.get("inname",None)
        hA      = kwargs.get("hA",None)
        hB      = kwargs.get("hB",None)
        for j in [inname,"o_"+inname]:
            if j in self.histogram:
                iname = j
                break
        self.histogram[name] = self.histogram[iname].Clone("h_"+name)
        HistA = self.histogram.get(hA,self.histogram["o_"+hA])
        HistB = self.histogram.get(hB,self.histogram["o_"+hB])
        scale = HistA.Integral()/HistA.Integral()
        self.histogram[name].Scale(scale)
        

    def f_derivedhistogram(self,min,max):
        for i in self.derivedhistogram :
            objtype = i.get("type","histogram")
            if objtype == "histogram":
                iname           = i["name"] 
                functionname    = i.get("function","f_Four_arithmetic_histogram")
                functionnkwargs = i.get("functionargs",{})
                functionnkwargs["name"] = iname
                getattr(self,functionname)(**functionnkwargs)

                histogram_on1D = [i["name"] for i in self.histogram_to_draw_pad2]
                histogram_on2D = [i["name"] for i in self.histogram_to_draw_pad1]
                if iname in (histogram_on1D+histogram_on2D):
                    self.histogram[iname].SetLineColor(self.color_palet[iname]); 
                    self.histogram[iname].SetFillColor(self.color_palet[iname]);
                    self.histogram[iname].SetMarkerColor(self.color_palet[iname]);
                self.histogram[iname].SetFillStyle(self.FillStyle.get(iname,0));
                self.histogram[iname].SetLineWidth(self.LineWidth.get(iname,2));

            if objtype == "fitfunction":
                function = i["functionargs"]["fitfunction"]
                hname    = i["functionargs"]["hist"]
                iname    = i["name"]
                self.FitHistogram(iname,hname,function,min,max)

            if objtype == "corr_histogram":
                name     = i["name"]
                hist     = i["functionargs"]["hist"]
                self.f_corrhistogram(hist=hist,name=name)

    def plotdeco_drawbookedhistogram(self,):
        for idecoh in self.histogram_booked :
            self.histogram[idecoh] = self.histogram_booked[idecoh]["h"].GetValue().Clone( self.histogram_booked[idecoh]["histogram_name"] )
            self.histogram[idecoh].SetLineWidth(self.LineWidth.get(idecoh,6));

    def plotdeco_totalhistogram(self,):
        self.h_fortoal = getattr(self,"h_fortoal",[j for j in self.DECO_hs])
        for iname in self.h_fortoal:
            print iname
            if iname in self.histogram:
                self.h_TotalMC.Add(self.histogram[iname])
        self.histogram["TotalMC"] = self.h_TotalMC
    
    def plotdeco_setRangeUserpad1(self,):
        if self.pad1range :
            for i in self.histogram_to_draw_pad1:
                iname = i["name"]
                if iname in self.histogram :
                    self.histogram[iname].GetYaxis().SetRangeUser(self.pad1range[0], self.pad1range[1])
        else :
            maxY_pad1 = 0
            minY_pad1 = 100000000
            for i in self.histogram_to_draw_pad1:
                iname = i["name"]
                if iname in self.histogram :
                    maxY_pad1 = TMath.Max(self.histogram[iname].GetMaximum(),maxY_pad1)
                    minY_pad1 = TMath.Min(self.histogram[iname].GetMinimum(),minY_pad1)
            for i in self.histogram_to_draw_pad1:
                iname = i["name"]
                if iname in self.histogram :
                    if self.Pad1LogY : 
                        self.histogram[iname].GetYaxis().SetRangeUser(minY_pad1*0.33,maxY_pad1*3)
                    else :
                        self.histogram[iname].GetYaxis().SetRangeUser(minY_pad1*0.8,maxY_pad1*1.3)

    def plotdeco_setlinepad1(self,):
        for i in self.histogram_to_draw_pad1:
            iname = i["name"]
            if iname in self.histogram :
                self.histogram[iname].SetLineWidth(self.LineWidth.get(iname,6));

    def plotdeco_setlinepad2(self,):
        for i in self.histogram_to_draw_pad2:
            iname = i["name"]
            if iname in self.histogram :
                self.histogram[iname].SetMarkerSize(self.MarkerSize.get(iname,2));
        

    def plotdeco_labelnameInorder(self,):
        # return the name of the label in order
        labels = []
        self.label_order = getattr(self,"label_order",[])
        histogram_to_draw1D = [i["name"] for i in self.histogram_to_draw_pad1]
        for i in self.label_order:
            if i in histogram_to_draw1D: labels.append(i)
        for i in histogram_to_draw1D:
            if i not in labels+["un"]: labels.append(i)
        for i in histogram_to_draw1D:
            if i == "un" : labels.append("un")
        return labels

    def plotdeco_AddText(self,):
        if self.histogram_Yields_ratio :
            Text_ = self.histogram_Yields_ratio["Text"]
            name1 = self.histogram_Yields_ratio["Hist"][0]
            name2 = self.histogram_Yields_ratio["Hist"][1]
            Text = Text_%( (self.histogram[name1].Integral()/self.histogram[name2].Integral()) )
            return Text

    def plotdeco_labelspad1(self,iname,percentage,debug):
        for i in self.histogram_to_draw_pad1:
            if iname == i["name"] : 
                Objname    = i.get("Objname","histogram")
                label_     = i["label"]
                
                labelstyle = i.get("label_style","L")
        label = label_
        if Objname == "histogram":
            if percentage :
                label = "%s : %s%%"%(label_,str( round((100*self.histogram[iname]).Integral()/(self.histogram["TotalMC"]).Integral()) ) )
        return label,labelstyle

    def plotdeco_labelspad2(self,iname):
        for i in self.histogram_to_draw_pad2:
            if iname == i["name"] : 
                label_     = i["label"]
                
                labelstyle = i.get("label_style","P")
                return label_,labelstyle

    def plotdeco_defineweight(self,name):
        weight = self.DECO_hs[name].get( "weight", None )
        IsData = False
        RootFiles = self.DECO_FIFLES[self.DECO_hs[name]["df_tree_name"]]
        print name,RootFiles
        for i in RootFiles :
            if "data" in i.lower() : 
                IsData = True
            if "muon" in i.lower() : 
                IsData = True
        if not weight :
            if not IsData : 
                weight = self.weight
        return weight


    def plotdeco_Generatehists(self,variable,variable_bin,xtitle,ytitle,nbin,min,max,Intime,cut,mutivariables,debug):
        with open("config/Input/%s"%(self.InputFileConfig),"r") as f:
            exec(f.read())

        # prepare Tree or RDF
        path = self.FilePath
        self.decodf_tree = {}
        if Intime:
            for iBKG in self.DECO_FIFLES:
                if type(path) == type([1,2]):
                    iBKGs = ROOT.std.vector('string')()
                    for ipath in path :
                        for n in self.DECO_FIFLES[iBKG]: 
                            iFile = self.InputFile[n]
                            iBKGs.push_back(ipath+iFile)
                else :
                    iBKGs = ROOT.std.vector('string')()
                    for n in self.DECO_FIFLES[iBKG]: 
                        iFile = self.InputFile[n]
                        iBKGs.push_back(path+iFile)
                self.decodf_tree[iBKG] = ROOT.RDataFrame("PKUTree", iBKGs) 
                self.AddIntimeVar(iBKG)

        with open(self.COM_py,"r") as f: exec(f.read())

        self.plotdeco_HistogramTF_Corrector()
        
        if not mutivariables:
            for idecoh in self.DECO_hs:
                name = idecoh
                Tree_RDF = self.decodf_tree[self.DECO_hs[idecoh]["df_tree_name"]]
                icut = self.DECO_cut[self.DECO_hs[idecoh].get("cut",idecoh)]
                Intime = Intime
                weight = self.plotdeco_defineweight(idecoh)
                self.histogramf(variable,variable_bin,name,xtitle,ytitle,nbin,min,max,Intime,Tree_RDF,icut,weight)
                print self.histogram.keys()
            if debug:
                name = "ROOTFile"
                Tree_RDF = self.decodf_tree[self.DECO_hs[idecoh]["df_tree_name"]]
                icut = cut 
                Intime = Intime
                weight = self.plotdeco_defineweight(idecoh)
                self.histogramf(variable,variable_bin,name,xtitle,ytitle,nbin,min,max,Intime,Tree_RDF,icut,weight)
            self.plotdeco_drawbookedhistogram()
            print self.histogram.keys()

        if mutivariables:
            for idecoh in self.DECO_hs:
                for index,ivariable in enumerate(mutivariables):
                    name = idecoh+str(index)
                    Tree_RDF = self.decodf_tree[self.DECO_hs[idecoh]["df_tree_name"]]
                    icut = self.DECO_cut[self.DECO_hs[idecoh].get("cut",idecoh)]
                    Intime = Intime
                    weight = self.plotdeco_defineweight(idecoh)
                    self.histogramf(ivariable,variable_bin,name,xtitle,ytitle,nbin,min,max,Intime,Tree_RDF,icut,weight)
            self.plotdeco_drawbookedhistogram()
            for idecoh in self.DECO_hs:
                for index,ivariable in enumerate(mutivariables):
                    if index == 0 : self.histogram[idecoh] = self.histogram[idecoh+str(index)]
                    else : self.histogram[idecoh] += self.histogram[idecoh+str(index)]

    def FitHistogram(self,name,hname,function,min,max):
        histogram = self.histogram[hname]
        polfunc = TF1("polfunc", function, min, max)
        r = ROOT.TFitResultPtr() ; 
        r = histogram.Fit(polfunc,"S,N,0")	
        values = r.GetConfidenceIntervals(0.95, True)
        print "Chi2/Ndf = %0.3f"%(r.Chi2()/r.Ndf())
        polfunc.SetLineColor(1)
        polfunc.SetLineColor(self.histogram[hname].GetLineColor());
        polfunc.SetLineWidth(4)
        polfunc.SetLineStyle(1)

        ptstats = ROOT.TPaveStats(0.61,0.75,0.98,0.95,"brNDC");
        ptstats.SetBorderSize(1);
        ptstats.SetFillColor(0);
        #ptstats.SetTextAlign(12);
        ptstats.SetTextFont(42);
        ptstats.SetTextSize(0.035);
        ptstats.AddText("#chi^{2}/ndf=%0.3f/%s(%0.3f)"%(r.Chi2(),r.Ndf(),r.Chi2()/r.Ndf()));
        
        ptstats.AddText("p0       = %0.1f #pm %0.1f "%(r.Parameter(0),r.ParError(0)));
        ptstats.AddText("p1       = %0.1e #pm %0.1e "%(r.Parameter(1),r.ParError(1)));
        if function in ["pol2","pol3","pol4","pol5","pol6"]: ptstats.AddText("p2       = %0.1e #pm %0.1e "%(r.Parameter(2),r.ParError(2)));
        if function in [       "pol3","pol4","pol5","pol6"]: ptstats.AddText("p3       = %0.1e #pm %0.1e "%(r.Parameter(3),r.ParError(3)));
        if function in [              "pol4","pol5","pol6"]: ptstats.AddText("p4       = %0.1e #pm %0.1e "%(r.Parameter(4),r.ParError(4)));
        if function in [                     "pol5","pol6"]: ptstats.AddText("p5       = %0.1e #pm %0.1e "%(r.Parameter(5),r.ParError(5)));
        if function in [                            "pol6"]: ptstats.AddText("p6       = %0.1e #pm %0.1e "%(r.Parameter(6),r.ParError(6)));
        ptstats.SetOptStat(0);
        ptstats.SetOptFit(111);

        values = r.GetConfidenceIntervals(0.95, True)
        print "value:",values

        hint = histogram.Clone("hint")
        hint.SetMarkerSize(0)
        hint.SetLineWidth(0)

        for i in range(hint.GetNbinsX()):
            hint.SetBinContent(i+1,polfunc.Eval(hint.GetBinCenter(i+1)))
            hint.SetBinError(i+1,values[i])

        hint.SetStats(kFALSE);
        hint.SetFillColor(polfunc.GetLineColor());
        hint.SetFillStyle(3003)

        if not os.path.isdir(self.TFpath) : os.makedirs(self.TFpath)
        TFFile = "%s/%s"%(self.TFpath,self.TFfile)
        outf = ROOT.TFile(TFFile, "recreate")
        polfunc.Write()
        hint.Write()
        outf.Close()

        self.histogram["hint"] = hint

        self.fitfunction[name] = {
            "function"  : polfunc,
            "PaveStats" : ptstats,
        }

    def f_AddHint(self,**kwargs):
        hname = kwargs.get("hist")
        name  = kwargs.get("name")
        TFfile  = kwargs.get("TFfile")
        # name  = kwargs.get("TFFile")

        correction =  "hint"
        fin = ROOT.TFile(TFfile) ; TF_correction  = fin.Get(correction) ; TF_correction.SetDirectory(0)
        fin.Close()

        h_origin = self.histogram[hname]
        h_corr   = h_origin.Clone(name)

        Bins = h_corr.GetNbinsX()
        for i in range(1, Bins+1) :
            center = (float(h_corr.GetBinLowEdge(i))+float(h_corr.GetBinLowEdge(i+1)))/2.
            SF     = TF_correction.GetBinError(TF_correction.FindBin(center))*h_corr.GetBinContent(i)
            print "SF =============== ",SF
            h_corr.SetBinError(i,SF)
        print "h_corr color",h_corr.GetLineColor()
        h_corr.SetStats(kFALSE);
        h_corr.SetLineWidth(0)
        h_corr.SetFillColor(h_corr.GetLineColor());
        h_corr.SetFillStyle(3003)

        self.histogram[name] = h_corr

    def decideLegendPosition(self,):
        legendy_min = 0.9-0.08*6-0.005
        legendy_max = 0.92
        HasFitfunction = False
        for i in self.histogram_to_draw_pad1:
            if i.get("Objname",None) == "fitfunction" : HasFitfunction = True
        if HasFitfunction:
            legendy_max = 0.72
        return legendy_min,legendy_max

    def f_readHistogram(self,):
        self.inf  = ROOT.TFile(self.ROOTFile)
        for iname in [e.GetName() for e in self.inf.GetListOfKeys()] :
            name = "o_"+iname.replace("h_","")
            print name, type(self.inf.Get(iname))
            self.histogram[name] = self.inf.Get(iname)
        return [self.inf]

    def f_corrhistogram(self,**kwargs):
        hname = kwargs.get("hist")
        name  = kwargs.get("name")
        # name  = kwargs.get("TFFile")

        TFfile     = "%s/%s"%(self.TFpath,self.TFfile)
        correction =  "polfunc"
        fin = ROOT.TFile(TFfile) ; TF_correction  = fin.Get(correction)
        fin.Close()

        h_origin = self.histogram[hname]
        h_corr   = h_origin.Clone(name)

        Bins = h_corr.GetNbinsX()
        for i in range(1, Bins+1) :
            center = (float(h_corr.GetBinLowEdge(i))+float(h_corr.GetBinLowEdge(i+1)))/2.
            SF     = TF_correction.Eval(center)
            h_corr.SetBinContent(i,SF*h_origin.GetBinContent(i))

        self.histogram[name] = h_corr

    def plotdeco_HistogramTF(self,):
        TFPath = "/".join(os.path.normpath(self.BKGes_SFFile).split("/")[:-1])
        if not os.path.isdir(TFPath) : os.makedirs(TFPath)
        outf = ROOT.TFile( self.BKGes_SFFile, "recreate")
        for i in self.BKGes_HistogramTF : 
            TFname = i["TFname"]
            name   = i["name"]
            self.histogram[name].Clone(TFname).Write()
        outf.Clone()

    def plotdeco_PlotHistogramTF(self,**kwargs):
        name    = kwargs.get("name")
        File_TF = kwargs.get("File_TF")
        IBin    = kwargs.get("IBin")
        TFY = 1
        for i in File_TF :
            TFfile = i["file"]
            hname  = i["hname"]
            ix     = i.get("x",None)
            fin = ROOT.TFile(TFfile)
            if ix :
                ih   = fin.Get(hname)
                ibin = ih.GetXaxis().FindBin(ix)
                TFY  = TFY * ih.GetBinContent(ibin)
            else :
                histogram = fin.Get(hname).Clone(name)
                histogram.SetDirectory(0)
                print type(histogram)
            fin.Close()
        print type(histogram)
        Bins = histogram.GetNbinsX()
        for i in range(1, Bins + 1):
            histogram.SetBinContent( i, TFY*histogram.GetBinContent(i) );
        self.histogram[name] = histogram

    

    def plotdeco_HistogramTF_Corrector(self):
        if len(self.HistogramTF) == 0 : return 0
        for iHistogramTF in self.HistogramTF :
            File_TF = iHistogramTF
            weightname     = File_TF["weightname"]
            originalweight = File_TF["originalweight"]
            weight_BKGTF   = '''
    float {weightname} = {CenterWeight} ; 
            '''.format( CenterWeight = originalweight, weightname = weightname )
            for i in File_TF["corrs"] :
                TFfile   = i["file"]
                hname    = i["hname"]
                corrname = i["corrname"]
                corrvar  = i["corrvar"]
                fin = ROOT.TFile(TFfile)
                H_Temp = fin.Get(hname)
                # auto {corrname} = {hname}  ; {corrname}->SetDirectory(0) ; 
                ROOT.gInterpreter.ProcessLine(
                    '''
    auto {corrname} = {hname}  ; 
                    '''.format(
                        hname    = hname ,
                        corrname = corrname ,
                    )
                )
                fin.Close()
                weight_BKGTF += '''
    {weightname} = {weightname} * ({corrname}->Eval({corrvar})) ; 
    // std::cout << {corrname}->Eval(1500); << endl;
                '''.format( 
                    corrname = corrname ,
                    corrvar  = corrvar ,
                    weightname = weightname,
                )
            weight_BKGTF  += 'return {weightname} ; '.format(weightname = weightname)
            self.corr_weights = getattr(self,"corr_weights",{})
            self.corr_weights[weightname] = weight_BKGTF 

    def plotratio_fixratioRange(self,histogram,min_,max_):
        OutName = "v_"+histogram.GetName()
        histogram_out = histogram.Clone(OutName)
        for i in range( 1, histogram.GetNbinsX()+1 ) :
            BinC = histogram.GetBinContent(i)
            if BinC < min_ : BinC = min_
            if BinC > max_ : BinC = max_
            histogram_out.SetBinContent( i, BinC )
        return histogram_out

    def plot_deco_Pad2_ratioRange(self,minYpad2,maxYpad2):
        for index,i in enumerate(self.histogram_to_draw_pad2):
            Objname = i.get("Objname","histogram")
            Obj     = getattr(self,Objname)
            if Objname == "histogram" :
                iname = i["name"]
                oname = "v_%s"%(iname)
                hist  = self.histogram[iname]
                self.histogram_to_draw_pad2[index]["name"] = oname
                self.histogram[oname] = self.plotratio_fixratioRange(hist,minYpad2,maxYpad2)
            
    def construct_plot_deco(self,Nj,variable,cut,cut1,tag,nbin,min,max,xtitle="",ytitle="",logy=1,CR=0, Intime = False, variable_bin = None, percentage = False, debug = False , compare  = False , mutivariables = None):

        with open(self.COM_py,"r") as f : exec(f.read())

        self.h_TotalMC = None
        self.TextInfo = []

        self.FitObjDraw = {}

        if self.PlotPath:
            if not os.path.isdir( self.PlotPath ):
                print "mkdir ",self.PlotPath
                os.makedirs( self.PlotPath ) 

        print " -->  MODE:",self.MODE," variable:",variable,"\n       { "+cut+" }\n";

        if options.y=="16" : lumi= 36.3;
        if options.y=="17" : lumi= 41.5;
        if options.y=="18" : lumi= 59.8;
        if options.y=="161718" : lumi= 138;

        if self.MODE in ["BKG_DECO_vsDATA","BKG_DECO_compare"]:
            canvas_controlplot = TCanvas(self.REGION+"_"+variable, self.REGION+"_"+variable, self.plotsize[0],self.plotsize[1]);
            fPads1 = TPad("pad1", "", 0.0, self.Pad_Split, 1.00, 1.00);
            fPads2 = TPad("pad2", ""    , 0.0, 0.00, 1.00, self.Pad_Split);
            fPads1.SetBottomMargin(0.007);fPads1.SetLeftMargin( 0.10);fPads1.SetRightMargin( 0.03);
            fPads2.SetLeftMargin(  0.10 );fPads2.SetRightMargin(0.03);fPads2.SetBottomMargin(0.25);
            fPads1.Draw(); fPads2.Draw(); fPads1.cd();
        if self.MODE in ["BKG_DECO"]:
            canvas_controlplot = TCanvas(self.REGION+"_"+variable, self.REGION+"_"+variable, self.plotsize[0],self.plotsize[1]);
            canvas_controlplot.SetLeftMargin(0.1); canvas_controlplot.SetRightMargin(0.03); 
        
        self.histogram = {} ; self.histogram_booked = {} ;
        infile = []
        if self.readHistogram :
            infile = self.f_readHistogram()
        else : 
            self.plotdeco_Generatehists(variable,variable_bin,xtitle,ytitle,nbin,min,max,Intime,cut,mutivariables,debug)

        self.fitfunction = {}
        self.f_derivedhistogram(min,max)

        self.plotdeco_totalhistogram()
        
        self.plotdeco_setRangeUserpad1()

        self.plotdeco_setlinepad1()

        if self.MODE in ["BKG_DECO_vsDATA","BKG_DECO_compare"]:
            fPads1.cd() 
            if self.Pad1LogY : fPads1.SetLogy()

        for i in self.histogram_to_draw_pad1:
            iname   = i["name"]
            Objname = i.get("Objname","histogram")
            Obj     = getattr(self,Objname)
            if Objname == "histogram" :
                op = i.get("Drawoption","same HIST e")
                # Obj[iname] = self.Histogram_SetMinNone0(Obj[iname])
                Obj[iname].Draw(op)
            if Objname == "fitfunction":
                op_f = i.get("Drawoption_function","same")
                op_t = i.get("Drawoption_text"    ,"same")
                func = Obj[iname].get("function")
                text = Obj[iname].get("PaveStats",None)
                if op_f : func.Draw(op_f)
                if op_t : text.Draw(op_t)

        banner_text = "%s : %s fb^{-1} (13 TeV)"%(self.channel,str(lumi) )
        banner          = TLatex(0.96,0.96,banner_text);   banner.SetNDC();   banner.SetTextSize(0.034);     banner.SetTextFont(42);    banner.SetTextAlign(31);    banner.SetLineWidth(2);    banner.Draw();
        CMS             = TLatex(0.22,0.96,"CMS"                        );      CMS.SetNDC();      CMS.SetTextSize(0.042);        CMS.SetTextFont(42);       CMS.SetTextAlign(31);       CMS.SetLineWidth(2);       CMS.Draw();

        if self.MODE in ["BKG_DECO","BKG_DECO_compare","BKG_DECO_vsDATA"]:
            Extratext   = TLatex(0.24,0.96,"Preliminary" );Extratext.SetNDC();Extratext.SetTextSize(0.034);  Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        Text = self.REGION_text
        if self.plotdeco_AddText() :
            Text += " ,"+self.plotdeco_AddText()
        RegionTxt       = TLatex(0.15,0.88, Text );RegionTxt.SetNDC();RegionTxt.SetTextSize(0.042);  RegionTxt.SetTextFont(42);    RegionTxt.SetLineWidth(2); RegionTxt.Draw();

        legendy_min,legendy_max = self.decideLegendPosition()
        theLeg = TLegend(0.68,legendy_min,0.9,legendy_max, "", "NDC");theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42) ;theLeg.SetTextSize(.04);
        theLeg.SetFillColor(0); theLeg.SetBorderSize(0);theLeg.SetLineColor(0) ;theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42);#theLeg.SetNColumns(2);

        labels_name = self.plotdeco_labelnameInorder()

        for iname in labels_name :
            label,labelstyle = self.plotdeco_labelspad1(iname,percentage,debug)
            if iname in self.histogram :
                theLeg.AddEntry( self.histogram[iname],label,labelstyle);
            if iname in self.fitfunction:
                obj = self.fitfunction[iname]["function"]
                theLeg.AddEntry( obj,label,labelstyle);

        theLeg.SetY1NDC(legendy_min);
        theLeg.SetY1(theLeg.GetY1NDC()); theLeg.Draw();
            
        if self.MODE in ["BKG_DECO_vsDATA","BKG_DECO_compare"] : 
            fPads2.cd(); 
            minYpad2 = self.pad2Range[0]
            maxYpad2 = self.pad2Range[1]

            if self.fixRatio : self.plot_deco_Pad2_ratioRange(minYpad2,maxYpad2)

            self.plotdeco_setlinepad2()
            
            for index,i in enumerate(self.histogram_to_draw_pad2):
                iname   = i["name"]
                Objname = i.get("Objname","histogram")
                Obj     = getattr(self,Objname)
                if Objname == "histogram" :
                    op = i.get("Drawoption","same e")
                    # self.histogram[iname] = self.Histogram_SetMinNone0(self.histogram[iname])
                    self.histogram[iname].GetYaxis().SetRangeUser(minYpad2,maxYpad2);  
                    self.histogram[iname].Draw(op)
                    if index == 0 :
                        Yaxis = self.histogram[iname].GetYaxis();  
                        Yaxis.SetTitle("Ratio");  
                        Yaxis.SetTitleOffset(0.23);  
                        Yaxis.SetTitleSize(0.13);  
                        Yaxis.SetLabelSize(0.11); 
                        Yaxis.SetLabelSize(0.1); 
                        Yaxis.SetTitleOffset(0.7); 
                        Yaxis.SetTitleSize(0.14);
                        Xaxis = self.histogram[iname].GetXaxis()
                        Xaxis.SetLabelSize(0.1); 
                        Xaxis.SetTitleOffset(0.7); 
                        Xaxis.SetTitleSize(0.14);
                if Objname == "fitfunction":
                    op_f = i.get("Drawoption_function","same")
                    op_t = i.get("Drawoption_text"    ,"same")
                    func = Obj[iname].get("function")
                    text = Obj[iname].get("PaveStats",None)
                    if op_f : func.Draw(op_f)
                    if op_t : text.Draw(op_t)

            if variable_bin : 
                min_ = variable_bin[0] ; max_ = variable_bin[-1]
            else :
                min_ = min ; max_ = max
            axis1=TGaxis( min_,1,max_,1, 0,0,0, "L"); axis1.SetLineColor(1); axis1.SetLineWidth(1); axis1.Draw();

            theLeg1 = TLegend(0.75, 0.65, 0.9, 0.92, "", "NDC");theLeg1.SetName("theLeg1end"); theLeg1.SetBorderSize(0); theLeg1.SetLineColor(0); theLeg1.SetFillColor(0);theLeg1.SetFillStyle(0); theLeg1.SetLineWidth(0); theLeg1.SetLineStyle(0); theLeg1.SetTextFont(42);theLeg1.SetTextSize(.07);
            theLeg1.SetFillColor(0);theLeg1.SetBorderSize(0);theLeg1.SetLineColor(0);theLeg1.SetLineWidth(0);theLeg1.SetLineStyle(0);theLeg1.SetTextFont(42);#theLeg1.SetNColumns(2);

            for i in self.histogram_to_draw_pad2:
                iname = i["name"]
                label,labelstyle = self.plotdeco_labelspad2(iname)
                if iname in self.histogram :
                    theLeg1.AddEntry( self.histogram[iname],label,labelstyle);
                if iname in self.fitfunction :
                    obj = self.fitfunction[iname]["function"]
                    theLeg1.AddEntry( obj,label,labelstyle);
            theLeg1.Draw();

        if getattr(self,"BKGes_SFFile",None):
            self.plotdeco_HistogramTF()

        extension   = "";
        if tag    !=  "": extension = extension + "_"+tag;
        if logy         : extension = extension + "_log";
        #----------- Rename variables to a shorter name -----------------
        for c in [".","/","(",")","[","]","*","+",">","<"," ","=",",","deep","dnn","Decorr","jetAK8puppi","ass_tag","t_tag","_tag","|","&"]:variable=variable.replace(c,"_");
        for c in ["__","___","____","_____","______","_"]:variable=variable.replace(c,"");

        Name=self.REGION+"_"+variable+"_"+self.MODE+"_"+options.y+extension+".png"
        if self.Plotname:
            Name = self.Plotname+extension+".png"
        Name = self.PlotPath+Name
        file=TString(Name); 

        output_ROOTFile = Name.replace(".png",".root") ; outf = ROOT.TFile( output_ROOTFile, "recreate")
        for iname in self.histogram:
            if iname in self.SaveHistograms :
                self.histogram[iname].Write()
        outf.Close()

        canvas_controlplot.SaveAs( file.Data() )
    
    def Calculate_Significance(self,S,B):
        s = S.Integral()
        b = B.Integral()
        return (2*((s+b)*math.log((1+(s/b)),math.e)-s))**0.5

    def GenPlotPath(self,*args):
        postpath = "/".join(list(args))
        self.PlotPath = os.path.normpath(postpath)+"/"

    def Yields_Table(self,raw):
        Nsignificant_digits1 = 4
        Nsignificant_digits2 = 2
        content = "%.{N1}g +- %.{N2}g".format(N1=Nsignificant_digits1,N2=Nsignificant_digits2)
        Table = {}
        for iprocess in self.histogram_dict : 
            histogram = self.histogram_dict[iprocess]
            e = ctypes.c_double(0.)
            v = histogram.IntegralAndError(1, histogram.GetNbinsX(), e)
            Table[(raw,iprocess)] = content%(v,float(e.value))
        # for signal in self.signal_List_ToDraw :
        #     iprocess = "significance %s"%(signal)
        #     Table[(raw,iprocess)] = "%.2g"%(self.Calculate_Significance(self.histogram_dict[signal],self.histogram_dict["TotalMC"]))
        textP = "%s/YieldsTable/"%(self.PlotPath)
        if not os.path.isdir(textP) : os.makedirs(textP)
        text  = "%s/table.txt"%(textP)
        self.SaveTxt(text,Table)

    def BinYields_Table(self,):
        Nsignificant_digits1 = 4
        Nsignificant_digits2 = 2
        content = "%.{N1}g +- %.{N2}g".format(N1=Nsignificant_digits1,N2=Nsignificant_digits2)
        Table = {}
        for iprocess in self.histogram_dict :
            histogram = self.histogram_dict[iprocess] 
            for i in range(1,1+histogram.GetNbinsX()):
                raw   = "Bin%s"%(str(i))
                v = histogram.GetBinContent(i)
                e = histogram.GetBinError(i)
                Table[(raw,iprocess)] = content%(v,e)
        textP = "%s/Bin_YieldsTable/"%(self.PlotPath)
        if not os.path.isdir(textP) : os.makedirs(textP)
        text  = "%s/table.txt"%(textP)
        self.SaveTxt(text,Table)

    def Limit_Table(self,raw):
        Table = {}
        Nsignificant_digits = 2
        content = "%.{N}g,%.{N}g".format(N=Nsignificant_digits)
        for iop in self.limit_results :
            lower = self.limit_results[iop][0]
            uper  = self.limit_results[iop][1]
            Table[(raw,iop)] = content%(lower,uper)
        textP = "%s/LimitTable/"%(self.PlotPath)
        if not os.path.isdir(textP) : os.makedirs(textP)
        text  = "%s/table.txt"%(textP)
        self.SaveTxt(text,Table)

    def SaveTxt(self,text,Value):
        Path = "/".join(text.split("/")[:-1])
        if not os.path.isdir(Path) : os.makedirs(Path)
        if type(Value) != type("hh") : Value = str(Value)        
        with open(text,"w") as f : 
            f.write(Value)

    def construct_plot_getcut(self,Region,CutName = "Center"):
        if not CutName == "Center" :
            print "warning : using %s instead of center cut"%(CutName)
        self.CutName = CutName
        cutFile = "./config/Cuts/cut.py"
        with open(cutFile,"r") as f: 
            exec(f.read())
        return self.Cut[Region]

    def construct_plot(self,Nj,variable,cut,cut1,tag,nbin,min,max,xtitle="",ytitle="",logy=1,CR=0, Intime = False, variable_bin = None,**kwargs):

        with open("config/Input/%s"%(self.InputFileConfig),"r") as f:
            exec(f.read())
        
        if self.PlotPath:
            if not os.path.isdir( self.PlotPath ):
                print "mkdir ",self.PlotPath
                os.makedirs( self.PlotPath ) 
            if self.run_limit : 
                self.limit_FilePath = "%s/Limit/"%(self.PlotPath)
                if not os.path.isdir( self.limit_FilePath ):
                    os.makedirs( self.limit_FilePath ) 

        SFs=options.SFs; channel=options.channel; MODE=options.MODE;  REGION=options.REGION; 
        print " -->  MODE:",MODE," variable:",variable,"\n       { "+cut+" }\n"
        
        #----------------- paths to root files -------------------
        if options.y == "16" : path = self.path16 ;  lumi = 36.3;
        if options.y == "17" : path = self.path17 ;  lumi = 41.5;
        if options.y == "18" : path = self.path18 ;  lumi = 59.8;
        if options.y == "161718" : path = self.path161718 ;  lumi = 138 ;

        #====== DEFINE CANVAS ==========================
        if self.MODE in ["MC","MCvsDATA","COMP"]:
            if self.plotsize:
                canvas_controlplot = TCanvas(REGION+"_"+variable, REGION+"_"+variable, self.plotsize[0],self.plotsize[1]);
            else:
                canvas_controlplot = TCanvas(REGION+"_"+variable, REGION+"_"+variable, 700,700);
            fPads1 = TPad("pad1", "", 0.0, 0.29, 1.00, 1.00);
            fPads2 = TPad("pad2", ""    , 0.0, 0.00, 1.00, 0.29);
            if self.Pad_Split:
                fPads1 = TPad("pad1", "", 0.0, self.Pad_Split, 1.00, 1.00);
                fPads2 = TPad("pad2", ""    , 0.0, 0.00, 1.00, self.Pad_Split);
            fPads1.SetBottomMargin(0.007);fPads1.SetLeftMargin( 0.10);fPads1.SetRightMargin( 0.03);
            fPads2.SetLeftMargin(  0.10 );fPads2.SetRightMargin(0.03);fPads2.SetBottomMargin(0.25);
            fPads1.Draw(); fPads2.Draw(); fPads1.cd()

        if self.MODE in ["DECO"]:
            if self.plotsize:
                canvas_controlplot = TCanvas(REGION+"_"+variable, REGION+"_"+variable, self.plotsize[0],self.plotsize[1]);
            else:
                canvas_controlplot = TCanvas(REGION+"_"+variable, REGION+"_"+variable, 700,565);
            canvas_controlplot.SetLeftMargin(0.1); canvas_controlplot.SetRightMargin(0.03);  

        #====================== DEFINE TREES AND HISTOS ======================================
        # ROOT.EnableImplicitMT() # allow to use mutiple core

        for Intime_py in self.Intime_py_list:
            with open(Intime_py,"r") as f:
                exec(f.read())
        
        print "plot ===== > ",
        try:
            print eval(variable)
        except NameError:
            print variable

        if type(path) == type([1,2]) :
            # data
            datas = ROOT.std.vector('string')()
            for ipath in path :
                for idata in self.Data_List: 
                    datas.push_back(ipath+self.InputFile[idata])
            df_data  = ROOT.RDataFrame(self.TreeName, datas )
            # BKG
            for BKG in self.BKG_List:
                BKGs = ROOT.std.vector('string')()
                for ipath in path :
                    BKGs.push_back(ipath+self.InputFile[BKG])
                exec('df_%s  = ROOT.RDataFrame(self.TreeName, BKGs) '%( BKG))
            # signal
            for signal in self.signal_List:
                signals = ROOT.std.vector('string')()
                for ipath in path :
                    signals.push_back(ipath+self.InputFile[signal])
                exec('df_%s  = ROOT.RDataFrame(self.TreeName, signals) '%( signal))
        else :
            datas = ROOT.std.vector('string')()
            for idata in self.Data_List: datas.push_back(path+self.InputFile[idata])
            df_data  = ROOT.RDataFrame(self.TreeName, datas )
            for BKG in self.BKG_List:
                exec('df_%s  = ROOT.RDataFrame(self.TreeName, path + self.InputFile["%s"]) '%( BKG,BKG))
            for signal in self.signal_List:
                exec('df_%s  = ROOT.RDataFrame(self.TreeName, path + self.InputFile["%s"]) '%( signal,signal))
        
        if self.Intime_Add_Variables_BKG:
            for iselection in self.Intime_Add_Variables_BKG:
                for BKG in self.BKG_List:
                    if iselection not in [str(i) for i in eval('df_%s'%(BKG)).GetColumnNames()]:
                        print BKG,"Add ",iselection
                        exec('df_%s      = df_%s.Define(   iselection ,eval(iselection) )'%(BKG,BKG))
        
        if self.Intime_Add_Variables_data:
            for iselection in self.Intime_Add_Variables_data:
                if self.DrawData:
                    if iselection not in [str(i) for i in eval('df_%s'%("data")).GetColumnNames()]:
                        print "data Add ",iselection
                        df_data      = df_data.Define(   iselection ,eval(iselection) )
        
        Intime_Add_Variables_signal = []
        signal_weights = []
        for signal in self.signal_List :
            if signal in self.Signal_Weight : 
                Intime_Add_Variables_signal.append( self.Signal_Weight[signal] )
                signal_weights.append( self.Signal_Weight[signal] )
        if self.Intime_Add_Variables_signal:
            Intime_Add_Variables_signal += self.Intime_Add_Variables_signal
        if Intime_Add_Variables_signal:
            for iselection in Intime_Add_Variables_signal:
                for signal in self.signal_List:
                    signalweight = False
                    define       = True
                    if iselection in signal_weights : 
                        signalweight = True
                        if iselection != self.Signal_Weight[signal] : 
                            define = define & False
                    define = define & ( iselection not in [str(i) for i in eval('df_%s'%(signal)).GetColumnNames()] )
                    if define:
                        print signal,"Add ",iselection
                        if signalweight : 
                            fomula = self.signalweight_defined[iselection]
                            print fomula
                        else :
                            fomula = eval(iselection)
                        exec('df_%s      = df_%s.Define(   iselection , fomula )'%(signal,signal))

        hstack_TotalMC= THStack("hstack_TotalMC","hstack_TotalMC"+";%s;%s"%(xtitle,ytitle));                
        hstack_TotalMC_DECO= THStack("hstack_TotalMC_DECO","hstack_TotalMC_DECO"+";%s;%s"%(xtitle,ytitle));                

        #=================== SET WEIGHTS, SCALE TREES, DEFINE TOTAL AND STACK  =================================================
        weight=self.MCweight

        if not variable_bin:
            Bin = ("Plot",";%s;%s"%(xtitle,ytitle),nbin,min,max)
        if variable_bin:
            histo_bin = array('d',variable_bin)
            Bin = ("Plot",";%s;%s"%(xtitle,ytitle),(len(histo_bin)-1),histo_bin)

        # ======================================= debug =======================================
        if self.Debug_keep_Nentries:
            for signal in self.signal_List:
                exec('df_%s     = df_%s.Filter(    cut).Range(0,self.Debug_keep_Nentries)'%(signal,signal))
            for BKG in self.BKG_List:
                exec('df_%s     = df_%s.Filter(    cut).Range(0,self.Debug_keep_Nentries)'%(BKG,BKG))

        # ======================================= Draw =======================================

        for signal in self.signal_List:
            df = eval("df_%s"%(signal))
            exec('h_%s     = df_%s.Filter(    cut).Histo1D(Bin, variable, self.Signal_Weight[signal]).GetValue().Clone("h_%s"    )'%(signal,signal,signal))
        if self.MODE in ["MC","MCvsDATA"]:
            for BKG in self.BKG_List:
                df = eval("df_%s"%(BKG))
                if "weight" not in df.GetColumnNames() : 
                    print BKG
            for BKG in self.BKG_List:
                df = eval("df_%s"%(BKG))
                exec('h_%s     = df_%s.Filter(    cut).Histo1D(Bin, variable, weight).GetValue().Clone("h_%s"    )'%(BKG,BKG,BKG))
        if self.MODE in ["MCvsDATA"]:
            h_data     = df_data.Filter(    cut).Histo1D(Bin, variable).GetValue().Clone("h_data"    )

        if self.MODE == "DECO":
            for BKG in self.BKG_List:
                df = eval("df_%s"%(BKG))
                for DECO in self.DECO:
                    exec('h_{BKG}_{DECO}     = df_{BKG}.Filter( self.DECO_Cut[DECO] ).Histo1D(Bin, variable, weight).GetValue().Clone("h_{BKG}_{DECO}"    )'.format( BKG = BKG, DECO = DECO ) )
            for BKG in self.BKG_List:
                for DECO in self.DECO:
                    if DECO == self.DECO[0]:
                        exec('h_%s = h_%s_%s.Clone("h_%s")'%(BKG,BKG,DECO,BKG))
                    else:
                        exec('h_%s.Add(h_%s_%s)'%(BKG,BKG,DECO))
        
        # ============= Merge =============
        for mBKG in self.Merge_BKG_Dic:
            for BKG in self.Merge_BKG_Dic[mBKG]:
                if BKG == self.Merge_BKG_Dic[mBKG][0]:
                    exec('h_%s = h_%s.Clone("h_%s")'%(mBKG,BKG,mBKG))
                else:
                    exec('h_%s.Add(h_%s)'%(mBKG,BKG))

        if self.MODE == "DECO":
            for mBKG in self.Merge_BKG_Dic:
                for BKG in self.Merge_BKG_Dic[mBKG]:
                    for DECO in self.DECO:
                        if BKG == self.Merge_BKG_Dic[mBKG][0]:
                            exec('h_%s_%s = h_%s_%s.Clone("h_%s_%s")'%(mBKG,DECO,BKG,DECO,mBKG,DECO))
                        else:
                            exec('h_%s_%s.Add(h_%s_%s)'%(mBKG,DECO,BKG,DECO))


        # ============= UnderOverflow =============
        for h in self.Merge_BKG_Dic:
            exec('h_%s = UnderOverFlow1D(h_%s)'%(mBKG,mBKG))
        for signal in self.signal_List:
            exec('h_%s = UnderOverFlow1D(h_%s)'%(signal,signal))
        if self.MODE == "DECO":
            for h in self.Merge_BKG_Dic:
                for DECO in self.DECO:
                    exec('h_%s_%s = UnderOverFlow1D(h_%s_%s)'%(mBKG,DECO,mBKG,DECO))
        if self.DrawData:
            for h in [h_data] : h = UnderOverFlow1D(h)
        if self.DrawData:
            h_data.SetBinErrorOption(TH1D.kPoisson);
        
        # ============= correct BKG =============
        if self.BKG_Correct:
            h_TTbar_corr = h_TTbar.Clone("h_TTbar_c")
            Norm1 = h_TTbar_corr.Integral()
            h_TTbar = self.Correct_histogram(h_TTbar_corr, postfix = "TTbar",correction = "TTbarpol2_tf", TSFFile = self.TSFFile)
            Norm2 = h_TTbar.Integral()
            h_TTbar.Scale(Norm1/Norm2)
            print type(h_TTbar)

            h_WJets_corr = h_WJets.Clone("h_WJets_c")
            Norm3 = h_WJets_corr.Integral()
            h_WJets = self.Correct_histogram(h_WJets_corr, postfix = "WJets",correction = "WJetspol2_tf", TSFFile = self.TSFFile)
            Norm4 = h_WJets.Integral()
            h_WJets.Scale(Norm3/Norm4)
            print Norm1,"Norm3"
            print Norm2,"Norm4"
            print type(h_WJets)

        # ============= Stack =============
        First = True
        for mBKG in self.Merge_BKG_Dic_order:
            if mBKG not in self.Merge_BKG_Dic: continue
            hstack_TotalMC.Add(eval('h_%s'%(mBKG)))
            if First:
                h_TotalMC = eval('h_%s'%(mBKG)).Clone('h_TotalMC')
                First = False
            else:
                h_TotalMC.Add(eval('h_%s'%(mBKG)))

        # ============= Attach histogram to instance ============
        self.histogram_dict = {}
        if self.DrawData:
            self.histogram_dict["data"]  = h_data
        self.histogram_dict["TotalMC"]   = h_TotalMC
        for mBKG in self.Merge_BKG_Dic:
            self.histogram_dict[mBKG] = eval('h_%s'%(mBKG))
        for signal in self.signal_List_ToDraw:
            self.histogram_dict[signal] = eval('h_%s'%(signal))

        cutname = kwargs.get("cutname","PS")
        self.Yields_Table(cutname)
        self.BinYields_Table()

        if self.MODE == "DECO":
            if self.Only_Show_Component:
                for DECO in self.DECO:
                    First = True
                    for mBKG in self.Merge_BKG_Dic_order:
                        if mBKG not in self.Merge_BKG_Dic: continue
                        if First:
                            exec('h_%s = h_%s_%s.Clone("h_%s")'%(DECO,mBKG,DECO,DECO))
                            First = False
                        else:
                            exec('h_%s.Add(h_%s_%s)'%(DECO,mBKG,DECO))
                    hstack_TotalMC_DECO.Add(eval('h_%s'%(DECO)))
            else:
                for DECO in self.DECO:
                    for mBKG in self.Merge_BKG_Dic_order:
                        if mBKG not in self.Merge_BKG_Dic: continue
                        hstack_TotalMC_DECO.Add(eval('h_%s_%s'%(mBKG,DECO)))
        
        # ============= Norm =============
        if self.DrawData:
            print "data:",h_data.Integral()
            print "TotalMC:",h_TotalMC.Integral()

        for signal in self.signal_List:
            if self.Signal_Scale[signal] == -1:
                self.Signal_Scale[signal] = float("%.1g"%(h_TotalMC.GetMaximum()/eval('h_%s'%(signal)).GetMaximum()))
            exec('h_%s.Scale(self.Signal_Scale["%s"])'%(signal,signal))

        if self.DrawData:
            norm=h_data.Integral()/(h_TotalMC.Integral()+0.00001); 
            print "  norm=",norm;

        if self.MODE in ["MC","MCvsDATA"]: #---------- histogram cosmetics ---------------------
            if self.DrawData:
                h_data.SetLineColor(self.color_palet["data"]); h_data.SetFillColor(self.color_palet["data"]);
            for signal in self.signal_List:
                exec('h_%s.SetLineColor(self.color_palet["%s"])'%(signal,signal))
                exec('h_%s.SetLineWidth(4)'%(signal))
                # exec('h_%s.SetFillColor(self.color_palet["%s"])'%(signal,signal))
            for mBKG in self.Merge_BKG_Dic:
                exec('h_%s.SetLineColor(0)'%(mBKG))
                exec('h_%s.SetLineWidth(0)'%(mBKG))
                exec('h_%s.SetFillColor(self.color_palet["%s"])'%(mBKG,mBKG))
            h_TotalMC.SetLineStyle(3); h_TotalMC.SetMarkerStyle(0); h_TotalMC.SetLineWidth(5); h_TotalMC.SetLineColor(15);

        if self.MODE == "DECO":
            h_TotalMC.SetLineStyle(3); h_TotalMC.SetMarkerStyle(0); h_TotalMC.SetLineWidth(5); h_TotalMC.SetLineColor(15);
            for signal in self.signal_List:
                    exec('h_%s.SetLineColor(self.color_palet["%s"])'%(signal,signal))
                    exec('h_%s.SetLineWidth(4)'%(signal))
            if self.Only_Show_Component:
                for DECO in self.DECO:
                    exec('h_%s.SetLineColor(0)'%(DECO))
                    exec('h_%s.SetLineWidth(0)'%(DECO))
                    exec('h_%s.SetFillColor(self.color_DECO["%s"])'%(DECO,DECO))
            else:
                for mBKG in self.Merge_BKG_Dic:
                    for DECO in self.DECO:
                        exec('h_%s_%s.SetLineColor(0)'%(mBKG,DECO))
                        exec('h_%s_%s.SetLineWidth(0)'%(mBKG,DECO))
                        exec('h_%s_%s.SetFillColor(self.color_palet["%s"])'%(mBKG,DECO,mBKG))
                        exec('h_%s_%s.SetFillStyle(self.Fill_Style["%s"])'%(mBKG,DECO,DECO))
        
        #============ DRAW TOP PAD =====================
        if self.MODE in ["MC"]:
            print "5"; 
            h_TotalMC.Draw("e"); h_TotalMC.GetXaxis().SetNdivisions(509);
            hstack_TotalMC.Draw("same HIST"); # For unc-bars
            h_TotalMC.Draw("same e"   ); 
            for signal in self.signal_List_ToDraw:
                exec('h_%s.Draw("same HIST")'%(signal))

        if self.MODE in ["MCvsDATA"]:
            h_TotalMC.Draw("e"); h_TotalMC.GetXaxis().SetNdivisions(509);   
            h_data.Draw("e same"); h_data.GetXaxis().SetNdivisions(509);
            hstack_TotalMC.Draw("same HIST"); # For unc-bars
            h_data.Draw("same e");  #!needed 2nd time to draw data
            h_TotalMC.Draw("same e"   );
            for signal in self.signal_List:
                exec('h_%s.Draw("same HIST")'%(signal))
            canvas_controlplot.Update(); 

        if self.MODE == "DECO":
            h_TotalMC.Draw("e"); h_TotalMC.GetXaxis().SetNdivisions(509);
            hstack_TotalMC_DECO.Draw("same HIST"); # For unc-bars
            h_TotalMC.Draw("same e"   ); 
            for signal in self.signal_List:
                exec('h_%s.Draw("same HIST")'%(signal))


        #---------------- Add text in top pad -----------------------
        banner_Text = "%s : %s fb^{-1} (13 TeV)"%(self.channel,str(lumi))
        banner          = TLatex(0.96,0.96,banner_Text);   banner.SetNDC();   banner.SetTextSize(0.034);     banner.SetTextFont(42);    banner.SetTextAlign(31);    banner.SetLineWidth(2);    banner.Draw();
        CMS             = TLatex(0.22,0.96,"CMS" );      CMS.SetNDC();      CMS.SetTextSize(0.042);        CMS.SetTextFont(42);       CMS.SetTextAlign(31);       CMS.SetLineWidth(2);       CMS.Draw();
        if self.MODE=="MCvsData":
            Extratext   = TLatex(0.24,0.96,"Preliminary"                );Extratext.SetNDC();Extratext.SetTextSize(0.034);  Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        if self.MODE=="DECO" or MODE=="MC":
            Extratext   = TLatex(0.24,0.96,"Simulation"                 );Extratext.SetNDC();Extratext.SetTextSize(0.034);  Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        if self.REGION_text:
            RegionTxt       = TLatex(0.15,0.88,"%s"%(self.REGION_text)                );RegionTxt.SetNDC();RegionTxt.SetTextSize(0.042);  RegionTxt.SetTextFont(42);    RegionTxt.SetLineWidth(2); RegionTxt.Draw();
        else:
            RegionTxt       = TLatex(0.15,0.88,"%s"%(REGION)                );RegionTxt.SetNDC();RegionTxt.SetTextSize(0.042);  RegionTxt.SetTextFont(42);    RegionTxt.SetLineWidth(2); RegionTxt.Draw();
        if self.AddTxt:
            AddTxt       = TLatex(0.15,0.78,self.AddTxt                );AddTxt.SetNDC();AddTxt.SetTextSize(0.042);  AddTxt.SetTextFont(42);    AddTxt.SetLineWidth(2); AddTxt.Draw();
        if self.MODE in ["MCvsDATA"]:
            D_o_MC_txt  = TLatex(0.15,0.83,"Data / MC = %.2f"%(norm)    );D_o_MC_txt.SetNDC();D_o_MC_txt.SetTextSize(0.042);D_o_MC_txt.SetTextFont(42);   D_o_MC_txt.SetLineWidth(2);D_o_MC_txt.Draw();
        if SFs :
            SFsCorr     = TLatex(0.55, 0.96, "SFs corrected"            );   SFsCorr.SetNDC();   SFsCorr.SetTextSize(0.034);   SFsCorr.SetTextFont(52);   SFsCorr.SetTextAlign(11);   SFsCorr.SetLineWidth(2);   SFsCorr.Draw();
        #canvas_controlplot.Update(); 


        #========== DRAW BOTTOM PAD ============================================
        if self.MODE in ["MCvsDATA"]: #--------- Data / MC on 2nd pad ---------------------
            fPads2.cd(); 
            h_Ratio = h_data.Clone("h_Ratio"); h_Ratio.Divide( h_TotalMC ); MaxY=2; #TMath.Max( 2,  TMath.Min(3,h_Ratio.GetMaximum()*1.1) );
            MaxY =  getattr(self,"RatioMaxY",MaxY) 
            h_Ratio.SetLineColor(1); h_Ratio.SetLineWidth(2); h_Ratio.SetMarkerStyle(8); h_Ratio.SetMarkerSize(0.7); h_Ratio.GetYaxis().SetRangeUser( 0 , MaxY );  h_Ratio.GetYaxis().SetNdivisions(504,0);
            h_Ratio.GetYaxis().SetTitle("Data / MC  ");  h_Ratio.GetYaxis().SetTitleOffset(0.35);  h_Ratio.GetYaxis().SetTitleSize(0.13);  h_Ratio.GetYaxis().SetTitleSize(0.13);  h_Ratio.GetYaxis().SetLabelSize(0.11); h_Ratio.GetXaxis().SetLabelSize(0.1); h_Ratio.GetXaxis().SetTitleOffset(0.7); h_Ratio.GetXaxis().SetTitleSize(0.14); 
            axis1=TGaxis( min,1,max,1, 0,0,0, "L"); axis1.SetLineColor(1); axis1.SetLineWidth(1);  #axis1->SetLabelColor(16); #fPads2.SetGridx(); #fPads2.SetGridy();
            h_Ratio=RationUnc(h_data,h_TotalMC,h_Ratio,MaxY);
            h_Ratio.Draw("e0"); axis1.Draw();
            fPads2.RedrawAxis(); fPads2.Update();
            fPads1.RedrawAxis(); fPads1.Update();

        if self.MODE in ["MC"]:  #--------- Significances on 2nd pad ---------------------
            # fPads2.cd();  fPads2.SetLogy(); MaxY=7;
            fPads2.cd();  MaxY=7;
            axis2=TGaxis( min,2,max,2, 0,0,0, "L"); axis2.SetLineColor(2); axis2.SetLineWidth(1);
            axis3=TGaxis( min,5,max,5, 0,0,0, "L"); axis3.SetLineColor(3); axis3.SetLineWidth(1);
            #-------- build denominator of significance ------------
            MaxY = 0
            for signal in self.signal_List_ToDraw:
                exec('h_Signif_{signal} = self.Significance_Histogram(h_{signal}, h_TotalMC,"h_Signif_{signal}")'.format(signal = signal))
                eval('h_Signif_{signal}'.format(signal = signal)).GetYaxis().SetTitle("significance")
                eval('h_Signif_{signal}'.format(signal = signal)).GetYaxis().SetTitleOffset(0.38)
                eval('h_Signif_{signal}'.format(signal = signal)).GetYaxis().SetTitleSize(0.13)
                eval('h_Signif_{signal}'.format(signal = signal)).GetYaxis().SetLabelSize(0.11)
                eval('h_Signif_{signal}'.format(signal = signal)).GetYaxis().SetLabelSize(0.1)
                eval('h_Signif_{signal}'.format(signal = signal)).GetXaxis().SetTitleOffset(0.7)
                eval('h_Signif_{signal}'.format(signal = signal)).GetXaxis().SetTitleSize(0.14)
                MaxY = TMath.Max(MaxY, eval('h_Signif_{signal}'.format(signal = signal)).GetMaximum())
            MaxY=MaxY*1.5; MinY=0.01;
            for signal in self.signal_List_ToDraw:
                eval('h_Signif_{signal}'.format(signal = signal)).GetYaxis().SetRangeUser(MinY,MaxY);
            for signal in self.signal_List_ToDraw:
                if signal == self.signal_List_ToDraw[0] : 
                    eval('h_Signif_{signal}'.format(signal = signal)).Draw("hist")
                    eval('h_Signif_{signal}'.format(signal = signal)).GetXaxis().SetLabelSize(0.1); eval('h_Signif_{signal}'.format(signal = signal)).GetXaxis().SetTitleOffset(0.7); eval('h_Signif_{signal}'.format(signal = signal)).GetXaxis().SetTitleSize(0.14); 
                else:
                    eval('h_Signif_{signal}'.format(signal = signal)).Draw("hist same")

            if self.Optimal:
                length=0.06*(h_TotalMC.GetBinLowEdge(h_TotalMC.GetNbinsX()+1)-h_TotalMC.GetBinLowEdge(1));
                Latexposition = 0.8 
                arrposition   = 0.7
                for isignal in self.signal_List_ToDraw :
                    Latexposition += -0.2
                    arrposition   += -0.2
                    cmd = '''V{isignal} = OptimalCut(h_TotalMC, h_{isignal});
sig{isignal}_Leftcut =TGaxis( V{isignal}[2], MinY, V{isignal}[2],MaxY, 0,0,0, "L"); sig{isignal}_Leftcut.SetLineColor(self.color_palet[isignal]);   sig{isignal}_Leftcut.SetLineWidth(1);  sig{isignal}_Leftcut.Draw();
sig{isignal}_Rightcut=TGaxis( V{isignal}[3], MinY, V{isignal}[3],MaxY, 0,0,0, "L"); sig{isignal}_Rightcut.SetLineColor(self.color_palet[isignal]);  sig{isignal}_Rightcut.SetLineWidth(1); sig{isignal}_Rightcut.Draw();
PrintValus{isignal} = TLatex( h_TotalMC.GetBinLowEdge(2), MaxY*Latexposition,"    r "+V{isignal}[4]+"%   e "+V{isignal}[5]+"%   s "+V{isignal}[6]+"%"); PrintValus{isignal}.SetTextColor(self.color_palet[isignal]); PrintValus{isignal}.SetTextSize(0.1); PrintValus{isignal}.Draw();
ar_L{isignal} = TArrow( V{isignal}[2] ,MaxY*arrposition, V{isignal}[2]+length, MaxY*arrposition, 0.03 , "|>"); ar_L{isignal}.SetLineWidth(3); ar_L{isignal}.SetFillColor(self.color_palet[isignal]); ar_L{isignal}.SetLineColor(self.color_palet[isignal]); ar_L{isignal}.SetAngle(30); ar_L{isignal}.Draw();
ar_R{isignal} = TArrow( V{isignal}[3] ,MaxY*arrposition, V{isignal}[3]-length, MaxY*arrposition, 0.03 , "|>"); ar_R{isignal}.SetLineWidth(3); ar_R{isignal}.SetFillColor(self.color_palet[isignal]); ar_R{isignal}.SetLineColor(self.color_palet[isignal]); ar_R{isignal}.SetAngle(30); ar_R{isignal}.Draw();
                    '''
                    exec(cmd.format(isignal = isignal))

            # axis2.Draw();axis3.Draw(); 
                
            fPads2.RedrawAxis(); fPads2.Update();
            fPads1.RedrawAxis(); fPads1.Update();


        #============= THE LEGEND SESSION =======================
        if self.MODE in ["MC","MCvsDATA"]:
            theLeg = TLegend(0.48, 0.55, 0.9, 0.9, "", "NDC");theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(.05);
            theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42);#theLeg.SetNColumns(2);
            if self.MODE=="MCvsDATA":theLeg.AddEntry(h_data, "Data "+options.y,"ep");
            for mBKG in self.Merge_BKG_Dic:
                theLeg.AddEntry(eval('h_'+mBKG)    , self.Label[mBKG]             ,"F");
            for signal in self.signal_List_ToDraw :
                theLeg.AddEntry(eval('h_'+signal),  self.Label[signal] + "#times %s"%(self.Signal_Scale[signal]),"L");
            theLeg.SetY1NDC(0.9-0.08*6-0.005);
            theLeg.SetY1(theLeg.GetY1NDC()); fPads1.cd(); theLeg.Draw(); #theLeg.AddEntry(gr_MCStat, "Sys.","F");
            #============ SET MAX Y-AXIS FOR PLOTS ==================
            histsigmax = 0
            histsigmin = 999999999
            for signal in self.signal_List_ToDraw:
                histsigmax = TMath.Max( histsigmax, eval('h_'+signal).GetMaximum() )
                histsigmin = TMath.Max( histsigmin, eval('h_'+signal).GetMinimum() )
            if self.MODE in ["MCvsDATA"]:
                histsigmax = TMath.Max( histsigmax, h_data.GetMaximum() )  
                print histsigmax
            histsigmax = TMath.Max( histsigmax, h_TotalMC.GetMaximum() )         
            histsigmin = TMath.Min( histsigmin, h_TotalMC.GetMinimum() )
            for signal in self.signal_List:
                eval('h_'+signal).GetYaxis().SetRangeUser(0, histsigmax*1.3 )
                h_TotalMC.GetYaxis().SetRangeUser(0, histsigmax*1.3 )
            if self.DrawData:
                h_data.GetYaxis().SetRangeUser(   0, histsigmax*1.3 )
            if logy == 1:
                fPads1.cd();  fPads1.SetLogy();
                if histsigmin<=0: histsigmin=1;
                for signal in self.signal_List:
                    eval('h_'+signal).GetYaxis().SetRangeUser(0.1*histsigmin, histsigmax*
                    3. )
                h_TotalMC.GetYaxis().SetRangeUser(0.1*histsigmin, histsigmax*3. )
                if self.DrawData:
                    h_data.GetYaxis().SetRangeUser(   0.1*histsigmin, histsigmax*3. )

        if self.MODE in ["DECO"]:
            theLeg = TLegend(0.48, 0.55, 0.9, 0.9, "", "NDC");theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(.05);
            theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42);#theLeg.SetNColumns(2);
            if self.Only_Show_Component:
                for DECO in self.DECO:
                    fraction = str(round(100*(eval('h_%s'%(DECO)).Integral()/h_TotalMC.Integral()),1))+"%"
                    theLeg.AddEntry(eval('h_%s'%(DECO))    , self.DECO_Label[DECO]+" "+  fraction   ,"F");
            else:
                theLeg.SetNColumns(len(self.DECO));
                theLeg.SetTextSize(.02);
                for mBKG in self.Merge_BKG_Dic:
                    for DECO in self.DECO:
                        fraction = str(round(100*(eval('h_%s_%s'%(mBKG,DECO)).Integral()/eval('h_%s'%(mBKG)).Integral()),1))+"%"
                        theLeg.AddEntry(eval('h_%s_%s'%(mBKG,DECO))    , "%s,%s %s"%(self.Label[mBKG],self.DECO_Label[DECO],fraction)   ,"F");
            for signal in self.signal_List:
                theLeg.AddEntry(eval('h_'+signal),  self.Label[signal] + "#times %s"%(self.Signal_Scale[signal]),"L");
            theLeg.SetY1NDC(0.9-0.08*6-0.005);
            # theLeg.SetY1(theLeg.GetY1NDC()); fPads1.cd(); theLeg.Draw(); #theLeg.AddEntry(gr_MCStat, "Sys.","F");
            theLeg.Draw()
            #============ SET MAX Y-AXIS FOR PLOTS ==================
            histsigmax = 0
            histsigmin = 999999999
            for signal in self.signal_List:
                histsigmax = TMath.Max( histsigmax, eval('h_'+signal).GetMaximum() )
                histsigmin = TMath.Max( histsigmin, eval('h_'+signal).GetMinimum() )
            histsigmax = TMath.Max( histsigmax, hstack_TotalMC_DECO.GetMaximum() )         
            histsigmin = TMath.Min( histsigmin, hstack_TotalMC_DECO.GetMinimum() )
            for signal in self.signal_List:
                eval('h_'+signal).GetYaxis().SetRangeUser(0, histsigmax*1.3 )
            hstack_TotalMC_DECO.GetYaxis().SetRangeUser(0, histsigmax*1.3 )
            
        #============ SAVE PLOTS IN A DIRECTORY ============================
        extension   = "";
        if tag    !=  "": extension = extension + "_"+tag;
        if logy         : extension = extension + "_log";
        if options.FBT  : extension = extension + "_FBT";
        if SFs          : extension = extension + "_SFsCorr";
        #----------- Rename variables to a shorter name -----------------
        for c in [".","/","(",")","[","]","*","+",">","<"," ","=",",","deep","dnn","Decorr","jetAK8puppi","ass_tag","t_tag","_tag","|","&"]:variable=variable.replace(c,"_");
        for c in ["__","___","____","_____","______","_"]:variable=variable.replace(c,"");
        #----------------- Save and open the plot -----------------------
        Name=REGION+"_"+variable+"_"+self.MODE+"_"+options.y+extension+".png"
        if self.Plotname:
            Name = self.Plotname+extension+".png"
        Name = self.PlotPath+"/"+Name
        file=TString(Name); 

        if self.run_limit : 
            # book all histogram to self
            self.LimitHistogram = {}
            for signal in self.signal_List : self.LimitHistogram[signal] = eval("h_"+signal)
            self.LimitHistogram["TotalMC"] = h_TotalMC 
            # generate histogram for limit caculation
            self.Limit_Histogram()


        output_ROOTFile = Name.replace(".png",".root") ; outf = ROOT.TFile( output_ROOTFile, "recreate")
        self.histogramROOTfile = output_ROOTFile
        for signal in self.signal_List:
            eval("h_%s"%(signal)).Write()
        for mBKG in self.Merge_BKG_Dic:
            eval("h_%s"%(mBKG)).Write()
        for hname in self.Limit_Histogram_3Term:
            self.Limit_Histogram_3Term[hname].Write()
        h_TotalMC.Write()
        outf.Close()

        if self.run_limit : 
            # create datacard and commmands
            self.Limit_Run()
            cutname = kwargs.get("cutname","PS")
            self.Limit_Table(cutname)

        self.comment = "No comment" ; Text_Name = Name.replace(".png",".txt") ; self.View_Selection_With_Plots(Text_Name,variable,cut,self.comment)

        canvas_controlplot.SaveAs( file.Data() )
        os.system("display %s &"%(Name) ); print "\n --> display %s &"%(Name);

        self.Plotname = None

        return_list = {
            "limit_results" : self.limit_results ,
        }

        return return_list

def Significance_Histogram(self, h_signal, h_TotalMC,name):
    h = h_signal.Clone(name)
    for i in range(1,h_TotalMC.GetNbinsX()+1,1):
        s = h_signal.GetBinContent(i)
        b = h_TotalMC.GetBinContent(i)
        if b == 0:
            b = 1
        significance = (2*((s+b)*math.log((1+(s/b)),math.e)-s))**0.5
        h.SetBinContent(i,significance)
    return h

ANALYSIS.Significance_Histogram = Significance_Histogram

################# Main Code ################################
def Draw_Control_Plot( channel ) :
    Instance_ANALYSIS = ANALYSIS( channel ); Instance_ANALYSIS.DefineSelection_0lep();
    #if channel in ["mu","el","lep"] : Instance_ANALYSIS = ANALYSIS( channel ); Instance_ANALYSIS.DefineSelection_1lep();

if __name__ == '__main__' : 
    Beginning = strftime("%H:%M:%S",gmtime())
    print '\n----RUN--------------channel:[',options.channel,']----------Region:[',options.REGION,']----------------[',Beginning,']--------'
    Draw_Control_Plot( options.channel );
    Finishing = strftime("%H:%M:%S",gmtime());
    #========== CALCULATE DURATION OF THE RUN ===========
    MIN=int(Finishing[3:5])-int(Beginning[3:5]); SEC=int(Finishing[6:8])-int(Beginning[6:8]); 
    if SEC<0 and MIN>0 : SEC=60+SEC; MIN=MIN-1;
    if SEC>0 and MIN<0 : MIN=60+MIN;
    if SEC<0 and MIN<0 : SEC=60+SEC; MIN=60+MIN-1;
    print '----END-----------------------------------------------[time:',Finishing,', duration:',MIN,'MIN:',SEC,'SEC]---------\n'
