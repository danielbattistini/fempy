import sys

from uproot.models.TH import Model_TH1D_v3

from ROOT import TH1F, TH1


def change_units(hist, multiplier):
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
        if (isinstance(se, TH1) and isinstance(me, TH1)):
            self.se = se.Clone()
            self.me = me.Clone()
        elif (isinstance(se, Model_TH1D_v3) and isinstance(me, Model_TH1D_v3)):
            self.se = se.to_pyroot()
            self.me = me.to_pyroot()

        self.se.Sumw2()

        self.cf = None
        self.norm = norm
        self.um = um

    def rebin(self, rebin):
        if rebin == 1:
            return

        self.se.Rebin(rebin)
        self.me.Rebin(rebin)

    def get_cf(self):
        if self.se == None and self.me == None and self.cf != None:
            return self.cf

        # change units
        if self.um == "MeV":
            title = ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)'
        elif self.um == "GeV":
            title = ';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)'
        elif self.um == "GeV2MeV":
            self.se = change_units(self.se, 1000)
            self.me = change_units(self.me, 1000)
            self.um = 'MeV'
            title = ';#it{k}* (MeV/#it{c});#it{C}(#it{k}*)'
        elif self.um == "MeV2GeV":
            self.se = change_units(self.se, 0.001)
            self.me = change_units(self.me, 0.001)
            self.um = 'GeV'
            title = ';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)'
        else:
            print("\033[31mError\033[0m]: Units are not valid. Exit!")
            sys.exit()

        # check that the normalization is valid
        if self.norm == None:
            if self.se.Integral() == 0:
                print("\033[33mWarning\033[0m: the same-event distribution has 0 entries. CF will be empty.")
                return TH1F('', '', self.se.GetNbinsX(), 0, self.se.GetXaxis().GetXmax())
            elif self.me.Integral() == 0:
                print("\033[33mWarning\033[0m: the mixed-event distribution has 0 entries. CF will be empty.")
                return TH1F('', '', self.se.GetNbinsX(), 0, self.se.GetXaxis().GetXmax())

            print("\033[33mWarning\033[0m: the normalization range is not specified. CF normalized to yields.")
            self.se.Scale(1./self.se.Integral())
            self.me.Scale(1./self.me.Integral())
        elif isinstance(self.norm, list):
            normBinFirst = self.se.FindBin(self.norm[0])
            normBinLast = self.se.FindBin(self.norm[1]-1.e-6)
            # if normBinFirst is not self.me.FindBin(self.norm[0]):
            #     print(f"Error: wrong binning {normBinFirst} {self.se.FindBin(self.norm[0])}")
            #     sys.exit()
            # if normBinLast is not self.me.FindBin(self.norm[1]-1.e-6):
            #     print(f"Error: wrong binning {normBinLast} {self.me.FindBin(self.norm[1]-1.e-6)}")
            #     sys.exit()

            if self.se.Integral(normBinFirst, normBinLast) == 0:
                print("\033[33mWarning\033[0m: the same-event distribution has 0 entries "
                      "in the specified normalization range. CF will be empty.")
                return TH1F('', '', self.se.GetNbinsX(), 0, self.se.GetXaxis().GetXmax())
            if self.me.Integral(normBinFirst, normBinLast) == 0:
                print("\033[33mWarning\033[0m: the mixed-event distribution has 0 entries "
                      "in the specified normalization range. CF will be empty.")
                return TH1F('', '', self.se.GetNbinsX(), 0, self.se.GetXaxis().GetXmax())
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

    def get_se(self):
        return self.se

    def get_me(self):
        return self.me

    # def __add__(self, other):
    #     cf = self.cf.Clone()
    #     if isinstance(other, TH1):
    #         print('isth1')
    #         cf.Add(self.cf, other)
    #     elif isinstance(other, (int, float)):
    #         for iBin in range(cf.GetNbinsX()+2):
    #             cf.SetBinContent(iBin, self.cf.GetBinContent(iBin) + other)
    #     elif isinstance(other, CorrelationFunction):
    #         print('add cf')
    #         cf.Add(other.cf)
    #     else:
    #         print(f'Error: the type {type(other)} is not implemented. Exit!')
    #         sys.exit()
    #     return CorrelationFunction(cf)

    # def __sub__(self, other):
    #     cf = self.cf.Clone()
    #     if isinstance(other, TH1):
    #         cf.Add(self.cf, other, -1)
    #     elif isinstance(other, (int, float)):
    #         for iBin in range(cf.GetNbinsX()+2):
    #             cf.SetBinContent(iBin, self.cf.GetBinContent(iBin) - other)
    #     elif isinstance(other, CorrelationFunction):
    #         print('add cf')
    #         cf.Add(other.cf, -1)
    #     else:
    #         print(f'Error: the type {type(other)} is not implemented. Exit!')
    #         sys.exit()
    #     return CorrelationFunction(cf)

    # def __mul__(self, other):
    #     cf = self.cf.Clone()
    #     if isinstance(other, (int, float)):
    #         for iBin in range(cf.GetNbinsX()+2):
    #             cf.SetBinContent(iBin, self.cf.GetBinContent(iBin) * other)
    #             cf.SetBinError(iBin, self.cf.GetBinError(iBin) * other)
    #     else:
    #         print(f'Error: the type {type(other)} is not implemented. Exit!')
    #         sys.exit()
    #     return CorrelationFunction(cf)

    # def __truediv__(self, other):
    #     cf = self.cf.Clone()
    #     if isinstance(other, (int, float)):
    #         for iBin in range(cf.GetNbinsX()+2):
    #             cf.SetBinContent(iBin, self.cf.GetBinContent(iBin) / other)
    #             cf.SetBinError(iBin, self.cf.GetBinError(iBin) / other)
    #     else:
    #         print(f'Error: the type {type(other)} is not implemented. Exit!')
    #         sys.exit()
    #     return CorrelationFunction(cf)

    # def __radd__(self, other):
    #     return self.__add__(other)

    # def __rsub__(self, other):
    #     return self.__sub__(other)

    # def __rmul__(self, other):
    #     return self.__mul__(other)
