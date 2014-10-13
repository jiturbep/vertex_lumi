#ifndef _VTXTMATHALGTOOL_
#define _VTXTMATHALGTOOL_

#include "InDetTruthVertexValidation/IVertexTruthMatchAthAlgTool.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IPartPropSvc.h"

#include "GeneratorObjects/McEventCollection.h"
#include "TruthD3PDAnalysis/IGenObjectsFilterTool.h"

#include <string>
#include <vector>




class VertexTruthMatchAthAlgTool : virtual public AthAlgTool,
			virtual public IVertexTruthMatchAthAlgTool,
			virtual public IIncidentListener {

public:

  VertexTruthMatchAthAlgTool(const std::string& type, const std::string& name,
		      const IInterface* pParent);

  virtual ~VertexTruthMatchAthAlgTool();

  StatusCode queryInterface( const InterfaceID& riid, void** ppvIf );

  virtual StatusCode initialize();
  virtual StatusCode process();
  virtual StatusCode finalize();
  
  virtual void handle(const Incident& inc);
  
  virtual std::vector<VertexMatchInfo> getMatches() {return m_vtxtm.matchedVtx;}

  virtual InDetVertexTruthMatch * getMatcher() {return &m_vtxtm;}

  //accessors for trimmed vectors used to run the matching

  virtual int get_ngenevt() {return ngenevt;}
  virtual std::vector<short> get_genpileuptype() {return genpileuptype;}
  virtual std::vector<int> get_ngenparticle() {return ngenparticle;}

  virtual int get_ngenvtx() {return ngenvtx;}
  virtual std::vector<int> get_vtxevtindex() {return vtxevtindex;}
  virtual std::vector<float> get_genvtxx() {return genvtxx;}
  virtual std::vector<float> get_genvtxy() {return genvtxy;}
  virtual std::vector<float> get_genvtxz() {return genvtxz;}

  virtual int get_ngenpart() {return ngenpart;}
  virtual std::vector<int> get_partevtindex() {return partevtindex;}
  virtual std::vector<int> get_partprodvtxindex() {return partprodvtxindex;}
  virtual std::vector<int> get_partmctype() {return partmctype;}
  virtual std::vector<float> get_partpt() {return partpt;}
  virtual std::vector<float> get_parteta() {return parteta;}
  virtual std::vector<int> get_partbarcode() {return partbarcode;}
  virtual std::vector<int> get_partstatus() {return partstatus;}

  virtual int get_ntrks() {return ntrks;}
  virtual std::vector<float> get_trk_pt() {return trk_pt;}
  virtual std::vector<float> get_trk_eta() {return trk_eta;}
  virtual std::vector<float> get_trk_d0_wrtPV() {return trk_d0_wrtPV;}
  virtual std::vector<float> get_trk_z0_wrtPV() {return trk_z0_wrtPV;}
  virtual std::vector<float> get_trk_mc_probability() {return trk_mc_probability;}
  virtual std::vector<int> get_trk_mc_index() {return trk_mc_index;}

  virtual int get_nvtx() {return nvtx;}
  virtual std::vector<int> get_vx_trk_n() {return vx_trk_n;}
  virtual std::vector<std::vector<float> > get_vx_trk_weight() {return vx_trk_weight;}
  virtual std::vector<std::vector<int> > get_vx_trk_index() {return vx_trk_index;}
  virtual std::vector<float> get_vx_x() {return vx_x;}
  virtual std::vector<float> get_vx_y() {return vx_y;}
  virtual std::vector<float> get_vx_z() {return vx_z;}
  virtual std::vector<float> get_vx_cov_x() {return vx_cov_x;}
  virtual std::vector<float> get_vx_cov_y() {return vx_cov_y;}
  virtual std::vector<float> get_vx_cov_z() {return vx_cov_z;}

private:

  InDetVertexTruthMatch m_vtxtm;

  //Implement directly pass methods from GenObjectsFilterTool until better sol'n is found
  /// Function selecting GenEvent objects
  bool pass( const HepMC::GenEvent* evt,
	     const McEventCollection* coll = 0 ) const;
  /// Function selecting GenParticle objects
  bool pass( const HepMC::GenParticle* part,
		     const McEventCollection* coll = 0 ) const;
  /// Function selecting GenVertex objects
  bool pass( const HepMC::GenVertex* vtx,
	     const McEventCollection* coll = 0 ) const;

  //names
  std::string m_vxContainerName;
  std::string m_trackContainerName;
  std::string m_mcContainerName;
  std::string m_TruthMap;

  ServiceHandle< IPartPropSvc > m_partPropSvc;

  //config vars
  double m_ptMinCut; ///< Minimum track pT to be selected (in MeV)
  double m_etaCut; ///< Maximum track eta. Set to negative not to apply the cut
  
  bool m_removeEmptyEvents;
  bool m_removeDummyEvents;
  bool m_removeInTimePileUp;
  bool m_remove2BCPileUp;
  bool m_remove800nsPileUp;
  bool m_removeCavernBkg;
  bool m_selectTruthTracks; ///< Only select stable charged particles
  bool m_addOnlyFirstVertex;  ///< Only first vertex per gen event

  //vars to run truth matching with
  int ngenevt;
  std::vector<short> genpileuptype;
  std::vector<int> ngenparticle;

  int ngenvtx;
  std::vector<int> vtxevtindex;
  std::vector<float> genvtxx;
  std::vector<float> genvtxy;
  std::vector<float> genvtxz;

  int ngenpart;
  std::vector<int> partevtindex;
  std::vector<int> partprodvtxindex;
  std::vector<int> partmctype;
  std::vector<float> partpt;
  std::vector<float> parteta;
  std::vector<int> partbarcode;
  std::vector<int> partstatus;

  int ntrks;
  std::vector<float> trk_pt;
  std::vector<float> trk_eta;
  std::vector<float> trk_d0_wrtPV;
  std::vector<float> trk_z0_wrtPV;
  std::vector<float> trk_mc_probability;
  std::vector<int> trk_mc_index;

  int nvtx;
  std::vector<int> vx_trk_n;
  std::vector<std::vector<float> > vx_trk_weight;
  std::vector<std::vector<int> > vx_trk_index;
  std::vector<float> vx_x;
  std::vector<float> vx_y;
  std::vector<float> vx_z;
  std::vector<float> vx_cov_x;
  std::vector<float> vx_cov_y;
  std::vector<float> vx_cov_z;

};

#endif
