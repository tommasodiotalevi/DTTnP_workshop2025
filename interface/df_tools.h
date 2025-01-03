#ifndef DF_TOOLS
#define DF_TOOLS
#include <iostream>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "boost/lexical_cast.hpp"
#include "boost/property_tree/ini_parser.hpp"
#include "boost/property_tree/ptree.hpp"

using namespace ROOT;
using namespace ROOT::VecOps;

using rvec_f = const RVec<Float_t>&;
using rvec_i = const RVec<Int_t>&;
using rvec_b = const RVec<bool>&;
using rvec_s = const RVec<Short_t>&;
using rvec_ui = const RVec<UInt_t>&;
using rvec_uc = const RVec<UChar_t>&;

const Float_t MU_MASS = 0.106;

bool configBoolCaster(std::string value) { return (value == "true" || value == "True") ? true : false; }

std::string readIniValue(const std::string& filename, const std::string& section, const std::string& key) {
  boost::property_tree::ptree pt;

  try {
    boost::property_tree::ini_parser::read_ini(filename, pt);
  } catch (const boost::property_tree::ini_parser_error& ex) {
    std::cerr << "Error parsing INI file: " << ex.what() << std::endl;
    return "";
  }

  try {
    return pt.get<std::string>(section + "." + key);
  } catch (const boost::property_tree::ptree_bad_path& ex) {
    std::cerr << "Error reading value from INI file: " << ex.what() << std::endl;
    return "";
  }
}

template <typename T> std::vector<T> toArray(const std::string& entries) {
  std::vector<T> results;
  std::stringstream sentries(entries);
  std::string item;

  while (std::getline(sentries, item, ',')) {
    T result;
    if constexpr (std::is_floating_point<T>::value) {
      result = std::atof(item.c_str());
    } else if constexpr (std::is_integral<T>::value) {
      result = std::atoi(item.c_str());
    } else {
      result = T{(item.c_str())};
    }

    results.push_back(result);
  }

  return results;
}

// Load external variables from configuration file - configFile declared in notebook (ignore error)
const bool conf_tag_useIsoHltPath = configBoolCaster(readIniValue(configFile, "TagAndProbe", "tag_useIsoHltPath"));
const Float_t conf_tag_isoCut = std::stof(readIniValue(configFile, "TagAndProbe", "tag_isoCut"));
const Float_t conf_tag_minPt = std::stof(readIniValue(configFile, "TagAndProbe", "tag_minPt"));
const Int_t conf_probe_minPixelHits = std::stoi(readIniValue(configFile, "TagAndProbe", "probe_minPixelHits"));
const Int_t conf_probe_minTrkLayers = std::stoi(readIniValue(configFile, "TagAndProbe", "probe_minTrkLayers"));
const Int_t conf_probe_minNMatchedSeg = std::stoi(readIniValue(configFile, "TagAndProbe", "probe_minNMatchedSeg"));
const Int_t conf_probe_minNRPCLayers = std::stoi(readIniValue(configFile, "TagAndProbe", "probe_minNRPCLayers"));
const Float_t conf_probe_isoCut = std::stof(readIniValue(configFile, "TagAndProbe", "probe_isoCut"));
const Float_t conf_probe_minPt = std::stof(readIniValue(configFile, "TagAndProbe", "probe_minPt"));
const Float_t conf_pair_maxAbsDz = std::stof(readIniValue(configFile, "TagAndProbe", "pair_maxAbsDz"));
const Float_t conf_pair_minInvMass = std::stof(readIniValue(configFile, "TagAndProbe", "pair_minInvMass"));
const Float_t conf_pair_maxInvMass = std::stof(readIniValue(configFile, "TagAndProbe", "pair_maxInvMass"));
const Float_t conf_pair_minDr = std::stof(readIniValue(configFile, "TagAndProbe", "pair_minDr"));
const Float_t conf_passing_probe_maxTkSegDr =
    std::stof(readIniValue(configFile, "TagAndProbe", "passing_probe_maxTkSegDr"));
const Float_t conf_passing_probe_maxTkSegDx =
    std::stof(readIniValue(configFile, "TagAndProbe", "passing_probe_maxTkSegDx"));
const Float_t conf_passing_probe_maxTkSegDy =
    std::stof(readIniValue(configFile, "TagAndProbe", "passing_probe_maxTkSegDy"));
const std::vector<Float_t> conf_probe_maxAbsEta =
    toArray<Float_t>(readIniValue(configFile, "TagAndProbe", "probe_maxAbsEta"));
const Float_t conf_probe_maxBorderDx = std::stof(readIniValue(configFile, "TagAndProbe", "probe_maxBorderDx"));
const Float_t conf_probe_maxBorderDy = std::stof(readIniValue(configFile, "TagAndProbe", "probe_maxBorderDy"));

static constexpr Int_t N_WHEELS = 5;
static constexpr Int_t N_SECTORS = 14;
static constexpr Int_t N_STATIONS = 4;
static constexpr Int_t N_PHI_LAYERS = 8;
static constexpr std::size_t NULL_IDX = 999;

template <typename T>
RVec<RVec<T>> setProbeVars(Int_t nmuon, rvec_i muon_charge, rvec_b muon_isTight, rvec_b muon_isRPC, rvec_b muon_firesTrig,
                           rvec_b muon_firesIsoTrig, rvec_b muon_isTracker, rvec_s muon_trk_numberOfValidPixelHits,
                           rvec_s muon_trk_numberOfValidTrackerLayers, rvec_s muon_rpcMu_numberOfMatchedRPCLayers,
                           rvec_s muon_trk_origAlgo, rvec_f muon_pfIso04, rvec_f muon_pt, rvec_f muon_eta, rvec_f muon_phi,
                           rvec_f muon_trk_dz) {
  RVec<T> probeNPixelHits, probeNTrkLayers, probeNRPCLayers, probeOrigAlgo, probeReliso, pairMass, probePt, pairDr,
      pair_tag_id, pair_probe_id, pairDz;
  RVec<RVec<T>> probeIntVars, probeFloatVars;

  for (std::size_t iTag{0}; iTag != nmuon; ++iTag) {
    TLorentzVector tagVec;
    tagVec.SetPtEtaPhiM(muon_pt[iTag], muon_eta[iTag], muon_phi[iTag], MU_MASS);
    bool hasTrigger = conf_tag_useIsoHltPath ? muon_firesIsoTrig[iTag] : muon_firesTrig[iTag];
    bool tQ = muon_isTight[iTag] == true && muon_pfIso04[iTag] < conf_tag_isoCut && tagVec.Pt() > conf_tag_minPt;
    if (tQ && hasTrigger) {
      for (std::size_t iProbe{0}; iProbe != nmuon; ++iProbe) {
        if (iTag == iProbe)
          continue;
        TLorentzVector probeVec;
        probeVec.SetPtEtaPhiM(muon_pt[iProbe], muon_eta[iProbe], muon_phi[iProbe], MU_MASS);
        bool pQ = (muon_isTracker[iProbe] == true || muon_isRPC[iProbe] == 1) &&
                  muon_trk_origAlgo[iProbe] != 14 &&  // the track is not created out of a STA mu based seeding
                  muon_trk_numberOfValidPixelHits[iProbe] >= conf_probe_minPixelHits &&
                  muon_trk_numberOfValidTrackerLayers[iProbe] >= conf_probe_minTrkLayers &&
                  muon_pfIso04[iProbe] < conf_probe_isoCut && probeVec.Pt() > conf_probe_minPt;
        if constexpr (std::is_same<T, Int_t>::value) {
          probeNPixelHits.emplace_back(muon_trk_numberOfValidPixelHits[iProbe]);
          probeNTrkLayers.emplace_back(muon_trk_numberOfValidTrackerLayers[iProbe]);
          probeNRPCLayers.emplace_back(muon_isRPC[iProbe] ? muon_rpcMu_numberOfMatchedRPCLayers[iProbe] : 0);
          probeOrigAlgo.emplace_back(muon_trk_origAlgo[iProbe]);
        } else if constexpr (std::is_same<T, Float_t>::value) {
          probeReliso.emplace_back(muon_pfIso04[iProbe]);
        }

        if (pQ) {
          Float_t mass = (tagVec + probeVec).M();
          Float_t tnpDr = tagVec.DeltaR(probeVec);
          Float_t tnpDz = muon_trk_dz[iTag] - muon_trk_dz[iProbe];
          if constexpr (std::is_same<T, Float_t>::value) {
            pairMass.emplace_back(mass);
            probePt.emplace_back(probeVec.Pt());
            pairDr.emplace_back(tnpDr);
            pairDz.emplace_back(tnpDz);
          }
          if (std::abs(tnpDz) < conf_pair_maxAbsDz && muon_charge[iTag] * muon_charge[iProbe] == -1 &&
              mass > conf_pair_minInvMass && mass < conf_pair_maxInvMass && tnpDr > conf_pair_minDr) {
            if constexpr (std::is_same<T, Int_t>::value) {
              pair_tag_id.emplace_back(iTag);
              pair_probe_id.emplace_back(iProbe);
            }
            break;  // just one probe per tag
          }
        }
      }
    }
  }
  if constexpr (std::is_same<T, Int_t>::value) {
    probeIntVars.emplace_back(probeNPixelHits);
    probeIntVars.emplace_back(probeNTrkLayers);
    probeIntVars.emplace_back(probeNRPCLayers);
    probeIntVars.emplace_back(probeOrigAlgo);
    probeIntVars.emplace_back(pair_tag_id);
    probeIntVars.emplace_back(pair_probe_id);
    return probeIntVars;
  } else if constexpr (std::is_same<T, Float_t>::value) {
    probeFloatVars.emplace_back(probeReliso);
    probeFloatVars.emplace_back(pairMass);
    probeFloatVars.emplace_back(probePt);
    probeFloatVars.emplace_back(pairDr);
    probeFloatVars.emplace_back(pairDz);
    return probeFloatVars;
  }
}

Int_t nMatchedCh(Int_t iMu, Int_t iCh, rvec_uc muon_trkMu_stationMask) {
  Int_t nMatchedCh{};
  UInt_t chMask{muon_trkMu_stationMask[iMu]};

  for (int index{}; index != 8; ++index) {
    if ((chMask & 1 << index) && index != iCh - 1) {
      ++nMatchedCh;
    }
  }

  return nMatchedCh;
}

template <typename T> T getXY(T& arr, int x, int y) {
  return static_cast<T>((*((TVectorT<float>*)((arr)[x])))[y]);
};

std::size_t getPassingProbeInCh(Int_t stMu, Int_t secMu, Int_t whMu, Int_t xMu, Int_t yMu, Int_t ndtSegment,
                                rvec_s dtSegment_wheel, rvec_s dtSegment_sector, rvec_s dtSegment_station, rvec_f dtSegment_seg4D_posLoc_x,
                                rvec_f dtSegment_seg4D_posLoc_y, rvec_b dtSegment_seg4D_hasPhi, rvec_b dtSegment_seg4D_hasZed, rvec_s dtSegment_seg2D_phi_nHits,
                                rvec_s dtSegment_seg2D_z_nHits) {
  Int_t iBestSeg{NULL_IDX};
  Int_t bestSegNHits{};

  for (std::size_t iSeg{}; iSeg < ndtSegment; ++iSeg) {
    Int_t whSeg = dtSegment_wheel[iSeg];
    Int_t secSeg = dtSegment_sector[iSeg];
    Int_t stSeg = dtSegment_station[iSeg];
    Float_t xSeg = dtSegment_seg4D_posLoc_x[iSeg];
    Float_t ySeg = dtSegment_seg4D_posLoc_y[iSeg];
    bool hasPhi = dtSegment_seg4D_hasPhi[iSeg];
    bool hasTheta = dtSegment_seg4D_hasZed[iSeg];

    Float_t nHitsSeg{};
    if (hasPhi) {
      nHitsSeg += dtSegment_seg2D_phi_nHits[iSeg];
    }
    if (hasTheta) {
      nHitsSeg += dtSegment_seg2D_z_nHits[iSeg];
    }

    Float_t dX = std::abs(xSeg - xMu);
    Float_t dY = std::abs(ySeg - yMu);
    Float_t dR = sqrt(dX * dX + dY * dY);

    if (whMu == whSeg && secMu == secSeg && stMu == stSeg && nHitsSeg > bestSegNHits &&
        dR < conf_passing_probe_maxTkSegDr &&
        ((hasPhi && (dX < conf_passing_probe_maxTkSegDx)) || (hasTheta && (dY < conf_passing_probe_maxTkSegDy)))) {
      iBestSeg = iSeg;
      bestSegNHits = nHitsSeg;
    }
  }

  return iBestSeg;
}

std::pair<std::size_t, std::size_t> getPassingProbe(std::size_t iMu, Int_t iCh, Int_t ndtSegment,
                                                    rvec_s muon_matches_station, rvec_s muon_matches_wheel,
                                                    rvec_s muon_matches_sector, rvec_f muon_matches_x,
                                                    rvec_f muon_matches_y, rvec_ui muon_matches_begin, rvec_ui muon_matches_end, 
                                                    rvec_s dtSegment_wheel, rvec_s dtSegment_sector, rvec_s dtSegment_station, rvec_f dtSegment_seg4D_posLoc_x,
                                                    rvec_f dtSegment_seg4D_posLoc_y, rvec_b dtSegment_seg4D_hasPhi, rvec_b dtSegment_seg4D_hasZed,
                                                    rvec_s dtSegment_seg2D_phi_nHits, rvec_s dtSegment_seg2D_z_nHits) {
  Int_t iBestMatch{NULL_IDX};
  Int_t iBestSeg{NULL_IDX};
  Float_t bestSegDr{999.};

  auto nMatches{muon_matches_end[iMu] - muon_matches_begin[iMu]};

  for (std::size_t iMatch{0}; iMatch < nMatches; ++iMatch) {
    auto i{muon_matches_begin[iMu] + iMatch};
    Int_t stMu = muon_matches_station[i];

    if (stMu == iCh) {
      auto whMu = muon_matches_wheel[i];
      auto secMu = muon_matches_sector[i];

      auto xMu = muon_matches_x[i];
      auto yMu = muon_matches_y[i];

      auto iSeg = getPassingProbeInCh(stMu, secMu, whMu, xMu, yMu, ndtSegment, dtSegment_wheel, dtSegment_sector, dtSegment_station,
                                      dtSegment_seg4D_posLoc_x, dtSegment_seg4D_posLoc_y, dtSegment_seg4D_hasPhi, dtSegment_seg4D_hasZed, dtSegment_seg2D_phi_nHits, dtSegment_seg2D_z_nHits);

      if (iSeg == NULL_IDX)
        continue;

      Float_t xSeg = dtSegment_seg4D_posLoc_x[iSeg];
      Float_t ySeg = dtSegment_seg4D_posLoc_y[iSeg];

      Float_t dX = std::abs(xSeg - xMu);
      Float_t dY = std::abs(ySeg - yMu);
      Float_t dR = sqrt(dX * dX + dY * dY);

      if (dR < bestSegDr) {
        iBestMatch = iMatch;
        iBestSeg = iSeg;
        bestSegDr = dR;
      }
    }
  }
  return std::pair<std::size_t, std::size_t>(iBestSeg, iBestMatch);
}

template <typename T>
std::map<std::string, RVec<T>> segmentEfficiency(
    rvec_i pair_probe_id, rvec_s muon_rpcMu_numberOfMatchedRPCLayers, UInt_t run,
    rvec_f muon_pt, rvec_f muon_eta, rvec_f muon_phi, rvec_i muon_charge, rvec_uc muon_trkMu_stationMask, Int_t ndtSegment,
    rvec_s muon_matches_station, rvec_s muon_matches_wheel, rvec_s muon_matches_sector,
    rvec_f muon_matches_x, rvec_f muon_matches_y, rvec_f muon_matches_edgeX, rvec_f muon_matches_edgeY,
    rvec_ui muon_matches_begin, rvec_ui muon_matches_end, rvec_s dtSegment_wheel, rvec_s dtSegment_sector, rvec_s dtSegment_station, rvec_f dtSegment_seg4D_posLoc_x,
    rvec_f dtSegment_seg4D_posLoc_y, rvec_f dtSegment_seg2D_phi_t0, rvec_b dtSegment_seg4D_hasPhi, rvec_b dtSegment_seg4D_hasZed, rvec_s dtSegment_seg2D_phi_nHits,
    rvec_s dtSegment_seg2D_z_nHits /*, Int_t environment_instLumi, Short_t environment_nPV */) {
  std::map<std::string, RVec<T>> segmentEffVars;

  TLorentzVector probeVec;
  for (const auto& pair : pair_probe_id) {
    probeVec.SetPtEtaPhiM(muon_pt[pair], muon_eta[pair], muon_phi[pair], MU_MASS);
    for (Int_t iCh{1}; iCh <= N_STATIONS; ++iCh) {
      std::stringstream iChTag;
      iChTag << "MB" << iCh;

      const auto nMatchInOtherCh = nMatchedCh(pair, iCh, muon_trkMu_stationMask);
      if constexpr (std::is_same<T, Int_t>::value) {
        segmentEffVars["nMatchInOtherCh"].emplace_back(nMatchInOtherCh);
      } else if constexpr (std::is_same<T, Float_t>::value) {
        segmentEffVars["probeMuonEta"].emplace_back(muon_eta[pair]);
      }
      if (nMatchInOtherCh >= conf_probe_minNMatchedSeg ||
          muon_rpcMu_numberOfMatchedRPCLayers[pair] >= conf_probe_minNRPCLayers) {
        auto iPassingProbe{getPassingProbe(pair, iCh, ndtSegment, muon_matches_station, muon_matches_wheel,
                                           muon_matches_sector, muon_matches_x, muon_matches_y, muon_matches_begin, muon_matches_end, 
                                           dtSegment_wheel, dtSegment_sector, dtSegment_station, dtSegment_seg4D_posLoc_x, dtSegment_seg4D_posLoc_y, dtSegment_seg4D_hasPhi, 
                                           dtSegment_seg4D_hasZed, dtSegment_seg2D_phi_nHits, dtSegment_seg2D_z_nHits)};
        auto iPassingSeg{iPassingProbe.first};
        auto iPassingMatch{iPassingProbe.second};

        bool hasPhi{(iPassingSeg != NULL_IDX) && dtSegment_seg4D_hasPhi[iPassingSeg]};

        std::string hName = "effAccVsEta" + iChTag.str();
        if constexpr (std::is_same<T, Float_t>::value) {
          if (hasPhi) {
            segmentEffVars[hName + "passd"].emplace_back(muon_eta[pair]);
          }
          segmentEffVars[hName + "total"].emplace_back(muon_eta[pair]);

          if (muon_charge[pair] == 1) {
            hName = "effAccPhiVsEtaPlus" + iChTag.str();
            if constexpr (std::is_same<T, Float_t>::value) {
              if (hasPhi) {
                segmentEffVars[hName + "etapassd"].emplace_back(muon_eta[pair]);
                segmentEffVars[hName + "phipassd"].emplace_back(muon_phi[pair]);
              }
              segmentEffVars[hName + "etatotal"].emplace_back(muon_eta[pair]);
              segmentEffVars[hName + "phitotal"].emplace_back(muon_phi[pair]);
            }
          } else if (muon_charge[pair] == -1) {
            hName = "effAccPhiVsEtaMinus" + iChTag.str();
            if constexpr (std::is_same<T, Float_t>::value) {
              if (hasPhi) {
                segmentEffVars[hName + "etapassd"].emplace_back(muon_eta[pair]);
                segmentEffVars[hName + "phipassd"].emplace_back(muon_phi[pair]);
              }
              segmentEffVars[hName + "etatotal"].emplace_back(muon_eta[pair]);
              segmentEffVars[hName + "phitotal"].emplace_back(muon_phi[pair]);
            }
          }
        }
        if (std::abs(probeVec.Eta()) < conf_probe_maxAbsEta[iCh - 1]) {
          hName = "probePt" + iChTag.str();
          if constexpr (std::is_same<T, Float_t>::value) {
            segmentEffVars[hName].emplace_back(probeVec.Pt());
          }
          hName = "probeEta" + iChTag.str();
          if constexpr (std::is_same<T, Float_t>::value) {
            segmentEffVars[hName].emplace_back(probeVec.Eta());
          }
          hName = "probePhi" + iChTag.str();
          if constexpr (std::is_same<T, Float_t>::value) {
            segmentEffVars[hName].emplace_back(probeVec.Phi());
          }
          if (hasPhi) {
            auto i{muon_matches_begin[pair] + iPassingMatch};
            auto xMu = muon_matches_x[i];
            auto yMu = muon_matches_y[i];

            Float_t xSeg = dtSegment_seg4D_posLoc_x[iPassingSeg];
            Float_t ySeg = dtSegment_seg4D_posLoc_y[iPassingSeg];

            Float_t dX = std::abs(xSeg - xMu);
            Float_t dY = std::abs(ySeg - yMu);
            Float_t dR = sqrt(dX * dX + dY * dY);
            hName = "probePt_hasPhi" + iChTag.str();
            segmentEffVars[hName].emplace_back(probeVec.Pt());
            hName = "probeDrVsPt" + iChTag.str();
            if constexpr (std::is_same<T, Float_t>::value) {
              segmentEffVars[hName].emplace_back(dR);
            }
            hName = "probeDxVsPt" + iChTag.str();
            if constexpr (std::is_same<T, Float_t>::value) {
              segmentEffVars[hName].emplace_back(dX);
            }
            hName = "probeDyVsPt" + iChTag.str();
            if constexpr (std::is_same<T, Float_t>::value) {
              segmentEffVars[hName].emplace_back(dY);
            }
          }
          if (muon_charge[pair] == 1) {
            hName = "effAccVsPhiPlus" + iChTag.str();
            if constexpr (std::is_same<T, Float_t>::value) {
              if (hasPhi) {
                segmentEffVars[hName + "passd"].emplace_back(muon_phi[pair]);
              }
              segmentEffVars[hName + "total"].emplace_back(muon_phi[pair]);
            }
          } else if (muon_charge[pair] == -1) {
            hName = "effAccVsPhiMinus" + iChTag.str();
            if constexpr (std::is_same<T, Float_t>::value) {
              if (hasPhi) {
                segmentEffVars[hName + "passd"].emplace_back(muon_phi[pair]);
              }
              segmentEffVars[hName + "total"].emplace_back(muon_phi[pair]);
            }
          }
          hName = "effAccVsPt" + iChTag.str();
          if constexpr (std::is_same<T, Float_t>::value) {
            if (hasPhi) {
              segmentEffVars[hName + "passd"].emplace_back(probeVec.Pt());
            }
            segmentEffVars[hName + "total"].emplace_back(probeVec.Pt());
          }

          auto nMatches{muon_matches_end[pair] - muon_matches_begin[pair]};
          for (std::size_t iMatch = 0; iMatch < nMatches; ++iMatch) {
            auto i{muon_matches_begin[pair] + iMatch};
            auto stMu = muon_matches_station[i];
            auto xBorderMu = muon_matches_edgeX[i];
            auto yBorderMu = muon_matches_edgeY[i];

            if (stMu == iCh && xBorderMu < conf_probe_maxBorderDx && yBorderMu < conf_probe_maxBorderDy) {
              auto whMu = muon_matches_wheel[i];
              auto secMu = muon_matches_sector[i];

              auto xMu = muon_matches_x[i];
              auto yMu = muon_matches_y[i];

              iPassingSeg =
                  getPassingProbeInCh(stMu, secMu, whMu, xMu, yMu, ndtSegment, dtSegment_wheel, dtSegment_sector, dtSegment_station,
                                      dtSegment_seg4D_posLoc_x, dtSegment_seg4D_posLoc_y, dtSegment_seg4D_hasPhi, dtSegment_seg4D_hasZed, dtSegment_seg2D_phi_nHits, dtSegment_seg2D_z_nHits);
              hasPhi = (iPassingSeg != NULL_IDX) && dtSegment_seg4D_hasPhi[iPassingSeg];

              hName = "effVsPt" + iChTag.str();
              if constexpr (std::is_same<T, Float_t>::value) {
                if (hasPhi) {
                  segmentEffVars[hName + "passd"].emplace_back(probeVec.Pt());
                }
                segmentEffVars[hName + "total"].emplace_back(probeVec.Pt());
              }
              hName = "effSecVsWh" + iChTag.str();
              if constexpr (std::is_same<T, Float_t>::value) {
                if (hasPhi) {
                  segmentEffVars[hName + "secmupassd"].emplace_back(secMu + 0.5);
                  segmentEffVars[hName + "whmupassd"].emplace_back(whMu + 0.5);
                }
                segmentEffVars[hName + "secmutotal"].emplace_back(secMu + 0.5);
                segmentEffVars[hName + "whmutotal"].emplace_back(whMu + 0.5);
              }
              Float_t secBinOne = secMu == 13 ? 4 : secMu == 14 ? 10 : secMu;
              secBinOne += (stMu % 2 == 1) ? -0.1 : 0.1;

              Float_t whBinOne = whMu + ((stMu - 1) / 2 == 0 ? -0.1 : 0.1);

              hName = "effSecVsWhAllOne";
              if constexpr (std::is_same<T, Float_t>::value) {
                if (hasPhi) {
                  segmentEffVars[hName + "secpassd"].emplace_back(secBinOne);
                  segmentEffVars[hName + "whpassd"].emplace_back(whBinOne);
                }
                segmentEffVars[hName + "sectotal"].emplace_back(secBinOne);
                segmentEffVars[hName + "whtotal"].emplace_back(whBinOne);
              }
              Float_t secBinTwo = secMu == 13 ? 4 : secMu == 14 ? 10 : secMu;
              secBinTwo += 0.1;
              Float_t whBinTwo = whMu + 3.1 + ((stMu - 1) * 5);

              hName = "effSecVsWhAllTwo";
              if constexpr (std::is_same<T, Float_t>::value) {
                if (hasPhi) {
                  segmentEffVars[hName + "secpassd"].emplace_back(secBinTwo);
                  segmentEffVars[hName + "whpassd"].emplace_back(whBinTwo);
                }
                segmentEffVars[hName + "sectotal"].emplace_back(secBinTwo);
                segmentEffVars[hName + "whtotal"].emplace_back(whBinTwo);
              }

              hName = "effVsRun" + iChTag.str();
              if constexpr (std::is_same<T, Int_t>::value) {
                if (hasPhi) {
                  segmentEffVars[hName + "passd"].emplace_back(run);
                }
                segmentEffVars[hName + "total"].emplace_back(run);
              }

              hName = "effVsNHitsPhi" + iChTag.str();
              Int_t nPhiHits{hasPhi ? dtSegment_seg2D_phi_nHits[iPassingSeg] : 0};
              for (Int_t phiHits = 1; phiHits <= N_PHI_LAYERS; ++phiHits) {
                if constexpr (std::is_same<T, Int_t>::value) {
                  if (nPhiHits >= phiHits) {
                    segmentEffVars[hName + "passd"].emplace_back(phiHits);
                  }
                  segmentEffVars[hName + "total"].emplace_back(phiHits);
                }
              }
              if (hasPhi) {
                auto t0{dtSegment_seg2D_phi_t0[iPassingSeg]};
                auto secBin = secMu + ((whMu + 2) * 14);
                hName = "t0" + iChTag.str();
                if constexpr (std::is_same<T, Float_t>::value) {
                  segmentEffVars[hName + "secBin"].emplace_back(secBin);
                  segmentEffVars[hName + "t0"].emplace_back(t0);
                }
              }
            }
          }
        }
      }
    }
  }
  return segmentEffVars;
}

#endif