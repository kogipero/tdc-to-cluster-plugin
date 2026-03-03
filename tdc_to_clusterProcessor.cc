#include "tdc_to_clusterProcessor.h"

#include <JANA/JApplication.h>
#include <JANA/JEvent.h>

#include "algorithms/geo.h"

#include <Evaluator/DD4hepUnits.h>
#include <DD4hep/Objects.h>

#include <edm4eic/RawTrackerHit.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);
  app->Add(new tdc_to_clusterProcessor());
}
}

void tdc_to_clusterProcessor::clearBranches() {
  b_cellID.clear();

  b_pos_x.clear(); b_pos_y.clear(); b_pos_z.clear();
  b_posErr_xx.clear(); b_posErr_yy.clear(); b_posErr_zz.clear();

  b_time.clear();
  b_timeErr.clear();

  b_edep.clear();
  b_edepErr.clear();
}

void tdc_to_clusterProcessor::Init() {
  auto nn = algorithms::GeoSvc::instance().cellIDPositionConverter();
  m_converter = nn.get();

  // ---- output file ----
  m_outfile = std::make_unique<TFile>("calib_out.root", "RECREATE");
  m_tree = new TTree("events", "events");

  // ---- branches ----
  m_tree->Branch("TOFBarrelCalibratedHits.cellID", &b_cellID);

  m_tree->Branch("TOFBarrelCalibratedHits.position.x", &b_pos_x);
  m_tree->Branch("TOFBarrelCalibratedHits.position.y", &b_pos_y);
  m_tree->Branch("TOFBarrelCalibratedHits.position.z", &b_pos_z);

  m_tree->Branch("TOFBarrelCalibratedHits.positionError.xx", &b_posErr_xx);
  m_tree->Branch("TOFBarrelCalibratedHits.positionError.yy", &b_posErr_yy);
  m_tree->Branch("TOFBarrelCalibratedHits.positionError.zz", &b_posErr_zz);

  m_tree->Branch("TOFBarrelCalibratedHits.time", &b_time);
  m_tree->Branch("TOFBarrelCalibratedHits.timeError", &b_timeErr);

  m_tree->Branch("TOFBarrelCalibratedHits.edep", &b_edep);
  m_tree->Branch("TOFBarrelCalibratedHits.edepError", &b_edepErr);
}

void tdc_to_clusterProcessor::Process(const std::shared_ptr<const JEvent>& event) {
  using dd4hep::mm;

  clearBranches();

  auto raw_ptrs = event->Get<edm4eic::RawTrackerHit>("TOFBarrelADCTDC");
  if (raw_ptrs.empty()) {
    m_tree->Fill(); 
    return;
  }

  for (const auto* raw_ptr : raw_ptrs) {
    const auto& TDCADC_hit = *raw_ptr;

    auto id = TDCADC_hit.getCellID();

    dd4hep::Position pos;
    try {
      pos = m_converter->position(id);
    } catch (const std::exception&) {
      continue;
    }

    double ADC = TDCADC_hit.getCharge();
    double TDC = TDCADC_hit.getTimeStamp();

    // adc to charge
    double charge = ADC * m_cfg.c_slope + m_cfg.c_intercept;
    // TDC to time
    float time = static_cast<float>(TDC * m_cfg.t_slope + m_cfg.t_intercept);

    auto cellSize = m_converter->cellDimensions(id);

    // sqrt(12) factor convertes ranges of uniform distribution to it's standard deviation
    double varX = cellSize[0] / mm / std::sqrt(12.);
    varX *= varX;
    double varY = cellSize[1] / mm / std::sqrt(12.);
    varY *= varY;
    double varZ = cellSize.size() > 2 ? cellSize[2] / mm / std::sqrt(12.) : 0;
    varZ *= varZ;

    // time/charge error
    float timeErr  = static_cast<float>(m_cfg.t_slope / std::sqrt(12.));
    float edep     = static_cast<float>(std::max(0.0, charge));      
    float edepErr  = static_cast<float>(m_cfg.c_slope / std::sqrt(12.));

    // ---- fill vectors ----
    b_cellID.push_back(static_cast<std::uint64_t>(id));

    b_pos_x.push_back(static_cast<float>(pos.x()));
    b_pos_y.push_back(static_cast<float>(pos.y()));
    b_pos_z.push_back(static_cast<float>(pos.z()));

    b_posErr_xx.push_back(static_cast<float>(varX));
    b_posErr_yy.push_back(static_cast<float>(varY));
    b_posErr_zz.push_back(static_cast<float>(varZ));

    b_time.push_back(time);
    b_timeErr.push_back(timeErr);

    b_edep.push_back(edep);
    b_edepErr.push_back(edepErr);
  }

  m_tree->Fill();
}

void tdc_to_clusterProcessor::Finish() {
  if (!m_outfile) return;

  m_outfile->cd();
  if (m_tree) m_tree->Write();
  m_outfile->Close();
}