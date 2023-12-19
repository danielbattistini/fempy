import fempy

from ROOT import gStyle, gROOT, EColor

colors = {
    'kWhite': EColor.kWhite,
    'kBlack': EColor.kBlack,
    'kGray': EColor.kGray,
    'kRed': EColor.kRed,
    'kGreen': EColor.kGreen,
    'kBlue': EColor.kBlue,
    'kYellow': EColor.kYellow,
    'kMagenta': EColor.kMagenta,
    'kCyan': EColor.kCyan,
    'kOrange': EColor.kOrange,
    'kSpring': EColor.kSpring,
    'kTeal': EColor.kTeal,
    'kAzure': EColor.kAzure,
    'kViolet': EColor.kViolet,
    'kPink': EColor.kPink,
}


def GetColor(color:str):
    '''
    Converts a color from string to ROOT.

    Parameters
    ----------
    color : str
        ROOT color in str format

    Returns
    -------
    int
        The corresponding ROOT color.

    Examples
    --------
    >>> GetColor('kRed-3')
    629
    '''
    for name, value in colors.items():
        if name == color:
            return value
        if name in color:
            break

    for shade in range(-10, 11):
        abshade = abs(shade)
        if f'+{abshade}' in color or f' + {abshade}' in color or f' +{abshade}' in color or f'+ {abshade}' in color or \
           f'-{abshade}' in color or f' - {abshade}' in color or f' -{abshade}' in color or f'- {abshade}' in color:
            return value + shade

    fempy.logger.warn('Color %s could not be loaded. Using kBlack instead', color)
    return EColor.kBlack


def SetStyle(**kwargs):
    '''Sets the general style'''
    gStyle.SetAxisColor(1, 'xy')
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetCanvasBorderSize(0)
    gStyle.SetCanvasColor(0)
    gStyle.SetCanvasDefH(600)
    gStyle.SetCanvasDefW(600)
    gStyle.SetCanvasDefX(topx=10)
    gStyle.SetCanvasDefY(topy=10)
    gStyle.SetCapLinePS(0)
    gStyle.SetColorModelPS(0)
    gStyle.SetDateX(0.01)
    gStyle.SetDateY(0.01)
    gStyle.SetDrawBorder(0)
    gStyle.SetEndErrorSize(0)
    gStyle.SetErrorX(0.5)
    gStyle.SetFrameBorderMode(0)
    gStyle.SetFrameBorderSize(0)
    gStyle.SetFrameFillColor(1)
    gStyle.SetFrameFillStyle(0)
    gStyle.SetFrameLineColor(1)
    gStyle.SetFrameLineStyle(0)
    gStyle.SetFrameLineWidth(1)
    gStyle.SetFuncColor(1)
    gStyle.SetFuncStyle(1)
    gStyle.SetFuncWidth(4)
    gStyle.SetGridColor(0)
    gStyle.SetGridStyle(3)
    gStyle.SetGridWidth(1)
    gStyle.SetHistFillColor(1)
    gStyle.SetHistFillStyle(0)
    gStyle.SetHistLineColor(1)
    gStyle.SetHistLineStyle(0)
    gStyle.SetHistLineWidth(1)
    gStyle.SetHistMinimumZero(True)
    gStyle.SetHistTopMargin(0.05)
    gStyle.SetIsReading(True)
    gStyle.SetJoinLinePS(0)
    gStyle.SetLabelColor(1, 'xy')
    gStyle.SetLabelFont(42, 'xy')
    gStyle.SetLabelOffset(0.005, 'xy')
    gStyle.SetLabelSize(0.045, 'xyz')
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendFillColor(0)
    gStyle.SetLegendFont(42)
    gStyle.SetLegendTextSize(0.035)
    gStyle.SetLegoInnerR(0.5)
    gStyle.SetLineScalePS(3)
    gStyle.SetNdivisions(510, 'xy')
    gStyle.SetNumberContours(20)
    gStyle.SetOptDate(0)
    gStyle.SetOptFile(0)
    gStyle.SetOptFit(1)
    gStyle.SetOptLogx(0)
    gStyle.SetOptLogy(0)
    gStyle.SetOptLogz(0)
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    gStyle.SetPadBorderMode(1)
    gStyle.SetPadBorderSize(1)
    gStyle.SetPadColor(19)
    gStyle.SetPadGridX(False)
    gStyle.SetPadGridY(False)
    gStyle.SetPadTopMargin(kwargs.get('topMargin', 0.02))
    gStyle.SetPadBottomMargin(kwargs.get('bottomMargin', 0.11))
    gStyle.SetPadRightMargin(kwargs.get('rightMargin', 0.02))
    gStyle.SetPadLeftMargin(kwargs.get('leftMargin', 0.15))
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetPaintTextFormat('g')
    gStyle.SetPalette(55)
    gStyle.SetScreenFactor(1)
    gStyle.SetStatBorderSize(2)
    gStyle.SetStatColor(19)
    gStyle.SetStatFont(42)
    gStyle.SetStatFontSize(0)
    gStyle.SetStatFormat('6.4g')
    gStyle.SetStatH(0.1)
    gStyle.SetStatStyle(1001)
    gStyle.SetStatTextColor(1)
    gStyle.SetStatW(0.19)
    gStyle.SetStatX(0)
    gStyle.SetStatY(0)
    gStyle.SetStripDecimals(True)
    gStyle.SetTickLength(0.03, 'xy')
    gStyle.SetTitleAlign(13)
    gStyle.SetTitleBorderSize(2)
    gStyle.SetTitleColor(1, 'xy')
    gStyle.SetTitleFillColor(1)
    gStyle.SetTitleFont(42, 'xy')
    gStyle.SetTitleFontSize(0)
    gStyle.SetTitleH(0)
    gStyle.SetTitleOffset(1.4, 'yz')
    gStyle.SetTitleSize(0.045, 'xyz')
    gStyle.SetTitleStyle(1001)
    gStyle.SetTitleTextColor(1)
    gStyle.SetTitleW(0)
    gStyle.SetTitleX(0)
    gStyle.SetTitleXOffset(1)
    gStyle.SetTitleXSize(0.05)
    gStyle.SetTitleY(0.985)
    gStyle.SetTitleYOffset(1.4)
    gStyle.SetTitleYSize(0.05)

    # TAttLine
    gStyle.SetLineColor(1)
    gStyle.SetLineColorAlpha(1, 1)
    gStyle.SetLineWidth(1)

    # TAttMarker
    gStyle.SetMarkerColor(1)
    gStyle.SetMarkerColorAlpha(1, 1)
    gStyle.SetMarkerSize(1.5)
    gStyle.SetMarkerStyle(1)

    # TAttText
    gStyle.SetTextColor(1)
    gStyle.SetTextColorAlpha(1, 1)
    gStyle.SetTextFont(42)
    gStyle.SetTextSize(1)

    gROOT.ForceStyle()
