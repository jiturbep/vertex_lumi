
#ifndef _IVTXTMATHALGTOOL_
#define _IVTXTMATHALGTOOL_

#include "GaudiKernel/IAlgTool.h"

#include "InDetTruthVertexValidation/InDetVertexTruthMatch.h"
#include <vector>

static const InterfaceID IID_IVertexTruthMatchAthAlgTool( "IVertexTruthMatchAthAlgTool", 1, 0 );

class IVertexTruthMatchAthAlgTool : virtual public IAlgTool {

public:
  
  static const InterfaceID& interfaceID() { return IID_IVertexTruthMatchAthAlgTool; };

  virtual StatusCode process() = 0;

  virtual std::vector<VertexMatchInfo> getMatches() = 0;

  virtual InDetVertexTruthMatch * getMatcher() = 0;

  //accessors for trimmed vectors used to run the matching

  virtual int get_ngenevt()  = 0;
  virtual std::vector<short> get_genpileuptype()  = 0;
  virtual std::vector<int> get_ngenparticle()  = 0;

  virtual int get_ngenvtx()  = 0;
  virtual std::vector<int> get_vtxevtindex()  = 0;
  virtual std::vector<float> get_genvtxx()  = 0;
  virtual std::vector<float> get_genvtxy()  = 0;
  virtual std::vector<float> get_genvtxz()  = 0;

  virtual int get_ngenpart()  = 0;
  virtual std::vector<int> get_partevtindex()  = 0;
  virtual std::vector<int> get_partprodvtxindex()  = 0;
  virtual std::vector<int> get_partmctype()  = 0;
  virtual std::vector<float> get_partpt()  = 0;
  virtual std::vector<float> get_parteta()  = 0;
  virtual std::vector<int> get_partbarcode()  = 0;
  virtual std::vector<int> get_partstatus()  = 0;

  virtual int get_ntrks()  = 0;
  virtual std::vector<float> get_trk_pt()  = 0;
  virtual std::vector<float> get_trk_eta()  = 0;
  virtual std::vector<float> get_trk_d0_wrtPV()  = 0;
  virtual std::vector<float> get_trk_z0_wrtPV()  = 0;
  virtual std::vector<float> get_trk_mc_probability()  = 0;
  virtual std::vector<int> get_trk_mc_index()  = 0;

  virtual int get_nvtx()  = 0;
  virtual std::vector<int> get_vx_trk_n()  = 0;
  virtual std::vector<std::vector<float> > get_vx_trk_weight()  = 0;
  virtual std::vector<std::vector<int> > get_vx_trk_index()  = 0;
  virtual std::vector<float> get_vx_x()  = 0;
  virtual std::vector<float> get_vx_y()  = 0;
  virtual std::vector<float> get_vx_z()  = 0;
  virtual std::vector<float> get_vx_cov_x()  = 0;
  virtual std::vector<float> get_vx_cov_y()  = 0;
  virtual std::vector<float> get_vx_cov_z()  = 0;

};

#endif
