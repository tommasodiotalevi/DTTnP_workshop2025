import numpy as np
import os
from ROOT import TFile, TEfficiency

# Histogram booking functions

def book_1Dhistos(df, histograms, histo_dict, histo_settings):
    for key, value in histo_dict.items():
        print("Booked histo: " + key)
        if key.startswith('eff'):
            if 'ptBins' in histo_settings[key[:-3]]:
                histograms[key] = [df.Histo1D((key, histo_settings[key[:-3]]['title'], histo_settings[key[:-3]]['nbins'], np.array(histo_settings[key[:-3]]['ptBins'])), value + "passd"), 
                                df.Histo1D((key+"total", histo_settings[key[:-3]]['title'], histo_settings[key[:-3]]['nbins'], np.array(histo_settings[key[:-3]]['ptBins'])), value + "total")]
            else:
                histograms[key] = [df.Histo1D((key, histo_settings[key[:-3]]['title'], histo_settings[key[:-3]]['nbins'], histo_settings[key[:-3]]['xlow'], histo_settings[key[:-3]]['xhigh']), value + "passd"), 
                               df.Histo1D((key+"total", histo_settings[key[:-3]]['title'], histo_settings[key[:-3]]['nbins'], histo_settings[key[:-3]]['xlow'], histo_settings[key[:-3]]['xhigh']), value + "total")]
        else:
            if "MB" in key:
                histograms[key] = df.Histo1D((key, histo_settings[key[:-3]]['title'], histo_settings[key[:-3]]['nbins'], histo_settings[key[:-3]]['xlow'], histo_settings[key[:-3]]['xhigh']), value)
            else:
                histograms[key] = df.Histo1D((key, histo_settings[key]['title'], histo_settings[key]['nbins'], histo_settings[key]['xlow'], histo_settings[key]['xhigh']), value)


def book_2Dhistos(df, histograms, histo_dict, histo_settings):
    for key, value in histo_dict.items():
        print("Booked histo: " + key)
        if key.startswith('eff'):
            if "MB" in key:
                histograms[key] = [df.Histo2D((key, histo_settings[key[:-3]]['title'], histo_settings[key[:-3]]['nbinsx'], histo_settings[key[:-3]]['xlow'], histo_settings[key[:-3]]['xhigh'], 
                                           histo_settings[key[:-3]]['nbinsy'], histo_settings[key[:-3]]['ylow'], histo_settings[key[:-3]]['yhigh']), value[0] + "passd", value[1] + "passd"), 
                               df.Histo2D((key+"total", histo_settings[key[:-3]]['title'], histo_settings[key[:-3]]['nbinsx'], histo_settings[key[:-3]]['xlow'], histo_settings[key[:-3]]['xhigh'],
                                           histo_settings[key[:-3]]['nbinsy'], histo_settings[key[:-3]]['ylow'], histo_settings[key[:-3]]['yhigh']), value[0] + "total", value[1] + "total")]
            else:
                histograms[key] = [df.Histo2D((key, histo_settings[key]['title'], histo_settings[key]['nbinsx'], histo_settings[key]['xlow'], histo_settings[key]['xhigh'], 
                                           histo_settings[key]['nbinsy'], histo_settings[key]['ylow'], histo_settings[key]['yhigh']), value[0] + "passd", value[1] + "passd"), 
                               df.Histo2D((key+"total", histo_settings[key]['title'], histo_settings[key]['nbinsx'], histo_settings[key]['xlow'], histo_settings[key]['xhigh'],
                                           histo_settings[key]['nbinsy'], histo_settings[key]['ylow'], histo_settings[key]['yhigh']), value[0] + "total", value[1] + "total")]
        else:
            if "MB" in key:
                histograms[key] = df.Histo2D((key, histo_settings[key[:-3]]['title'], histo_settings[key[:-3]]['nbinsx'], histo_settings[key[:-3]]['xlow'], histo_settings[key[:-3]]['xhigh'], 
                                              histo_settings[key[:-3]]['nbinsy'], histo_settings[key[:-3]]['ylow'], histo_settings[key[:-3]]['yhigh']), value[0], value[1])
            else:
                histograms[key] = df.Histo2D((key, histo_settings[key]['title'], histo_settings[key]['nbinsx'], histo_settings[key]['xlow'], histo_settings[key]['xhigh'], 
                                              histo_settings[key]['nbinsy'], histo_settings[key]['ylow'], histo_settings[key]['yhigh']), value[0], value[1])

def book_profile(df, histograms, histo_dict, histo_settings):
    for key, value in histo_dict.items():
        print("Booked histo: " + key)
        if "MB" in key:
            histograms[key] = df.Profile1D((key, histo_settings[key[:-3]]['title'], histo_settings[key[:-3]]['nbinsx'], histo_settings[key[:-3]]['xlow'], histo_settings[key[:-3]]['xhigh'], 
                                        histo_settings[key[:-3]]['ylow'], histo_settings[key[:-3]]['yhigh']), value[0], value[1])
        else:
            histograms[key] = df.Profile1D((key, histo_settings[key]['title'], histo_settings[key]['nbinsx'], histo_settings[key]['xlow'], histo_settings[key]['xhigh'], histo_settings[key]['nbinsy'], 
                                        histo_settings[key]['ylow'], histo_settings[key]['yhigh']), value[0], value[1])

# Write histograms to file

def write_histos(histograms, config):
    if not os.path.exists(config['Data']['outputDirNameRDF']):
        os.makedirs(config['Data']['outputDirNameRDF'])
    file = TFile(os.path.join(config['Data']['outputDirNameRDF'], config['Data']['outputFileName']), 'recreate')
    for name, histo in histograms.items():
        if name.startswith('eff'):
            eff = TEfficiency(histo[0].GetValue(), histo[1].GetValue())
            eff.SetName(histo[0].GetName())
            eff.SetTitle(histo[0].GetTitle())
            eff.Write()
            #histo[0].Divide(histo[1].GetPtr())
            #histo[0].Write()
        else:
            histo.Write()
        print("Created Histo: " + name)
    file.Close()