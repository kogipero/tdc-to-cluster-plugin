#pragma once

#include <JANA/JEventProcessor.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>

#include <DD4hep/VolumeManager.h>   
#include <DDRec/CellIDPositionConverter.h> 

class tdc_to_clusterProcessor : public JEventProcessor {
public:
  tdc_to_clusterProcessor() = default;

  void Init() override;
  void Process(const std::shared_ptr<const JEvent>& event) override;
  void Finish() override;

private:
  struct LGADHitCalibrationCfg {
    double c_slope     = 1.0;
    double c_intercept = 0.0;
    double t_slope     = 1.0;
    double t_intercept = 0.0;
  };

  LGADHitCalibrationCfg m_cfg;

  const dd4hep::rec::CellIDPositionConverter* m_converter = nullptr;

  // ---- output ROOT ----
  std::unique_ptr<TFile> m_outfile;
  TTree* m_tree = nullptr;

  // ---- branches  ----
  std::vector<std::uint64_t> b_cellID;

  std::vector<float> b_pos_x, b_pos_y, b_pos_z;
  std::vector<float> b_posErr_xx, b_posErr_yy, b_posErr_zz;

  std::vector<float> b_time;
  std::vector<float> b_timeErr;

  std::vector<float> b_edep;
  std::vector<float> b_edepErr;

  void clearBranches();
};