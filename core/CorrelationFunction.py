import sys

from ROOT import TH1F


def ChangeUnits(hist, multiplier):
    '''Only for histogram with constant binwidth!'''
    nbins = hist.GetNbinsX()
    low_edge = hist.GetBinLowEdge(1)
    up_edge = hist.GetBinLowEdge(nbins+1)
    hSame_new = TH1F('', '', nbins, low_edge*multiplier, up_edge*multiplier)
    for i in range(0, nbins+2):
        hSame_new.SetBinContent(i, hist.GetBinContent(i))
        hSame_new.SetBinError(i, hist.GetBinError(i))
    return hSame_new


class CorrelationFunction:
    '''Class to handle correlation functions and related manipulations'''

    def __init__(self, se, me, norm=None, um='MeV'):
        '''
        Lazy constructor
        '''
        self.se = se.Clone()
        self.me = me.Clone()

        self.se.Sumw2()

        self.cf = None
        self.norm = norm
        self.um = um

    def Rebin(self, rebin):
        if rebin == 1:
            return

        self.se.Rebin(rebin)
        self.me.Rebin(rebin)

    def GetCF(self):
        # change units
        if self.um == "MeV":
            title = ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)'
        elif self.um == "GeV":
            title = ';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)'
        elif self.um == "GeV2MeV":
            self.se = ChangeUnits(self.se, 1000)
            self.me = ChangeUnits(self.me, 1000)
            self.um = 'MeV'
            title = ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)'
        elif self.um == "MeV2GeV":
            self.se = ChangeUnits(self.se, 0.001)
            self.me = ChangeUnits(self.me, 0.001)
            self.um = 'GeV'
            title = ';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)'
        else:
            print("\033[31mError\033[0m]: Units are not valid. Exit!")
            sys.exit()

        # check that the normalization is valid
        if self.norm == None:
            if self.se.Integral() == 0:
                print(
                    "\033[33mWarning\033[0m: the same-event distribution has 0 entries. CF is not computed.")
                return None
            elif self.me.Integral() == 0:
                print(
                    "\033[33mWarning\033[0m: the mixed-event distribution has 0 entries. CF is not computed.")
                return None

            print(
                "\033[33mWarning\033[0m: the normalization range is not specified. CF normalized to yields.")
            self.se.Scale(1./self.se.Integral())
            self.me.Scale(1./self.me.Integral())
        elif isinstance(self.norm, list):
            normBinFirst = self.se.FindBin(self.norm[0])
            normBinLast = self.se.FindBin(self.norm[1]-1.e-6)
            if normBinFirst is not self.me.FindBin(self.norm[0]):
                print("Error: wrong binning")
                sys.exit()
            if normBinLast is not self.me.FindBin(self.norm[1]-1.e-6):
                print("Error: wrong binning")
                sys.exit()

            if self.se.Integral(normBinFirst, normBinLast) == 0:
                print("\033[33mWarning\033[0m: the same-event distribution has 0 entries \
                      in the specified normalization range. CF is not computed.")
                return None
            elif self.me.Integral(normBinFirst, normBinLast) == 0:
                print("\033[33mWarning\033[0m: the mixed-event distribution has 0 entries \
                      in the specified normalization range. CF is not computed.")
                return None
            else:
                self.se.Scale(1./self.se.Integral(normBinFirst, normBinLast))
                self.me.Scale(1./self.me.Integral(normBinFirst, normBinLast))
        else:
            print(
                f"\033[31mError:\033[0m normalization {self.norm} not supported. Exit!")
            sys.exit()

        self.cf = self.se.Clone()
        self.cf.Divide(self.me)
        self.cf.SetTitle(title)

        return self.cf

    def GetSE(self):
        return self.se

    def GetME(self):
        return self.me
