#include "./VertexTruthMatchAthAlgTool.h"

#include "HepPDT/ParticleDataTable.hh"

// Helper functors:
#include "TruthHelper/IsGenStable.h"
#include "TruthHelper/IsGenInteracting.h"

// EDM include(s):
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"

#include "EventInfo/EventID.h"
#include "EventInfo/EventType.h"
#include "EventInfo/EventInfo.h"
#include "VxVertex/VxContainer.h"
#include "VxVertex/VxTrackAtVertex.h"
#include "Particle/TrackParticleContainer.h"
#include "ParticleTruth/TrackParticleTruthCollection.h"
#include "TrkParticleBase/LinkToTrackParticleBase.h"

VertexTruthMatchAthAlgTool::VertexTruthMatchAthAlgTool(const std::string& type, const std::string& name,
				     const IInterface* pParent) :
  AthAlgTool(type, name, pParent),
  m_partPropSvc( "PartPropSvc", name )
{

  declareInterface< IVertexTruthMatchAthAlgTool >(this);

  declareProperty("vxContainerName",m_vxContainerName="VxPrimaryCandidate");
  declareProperty("trackContainerName",m_trackContainerName="TrackParticleCandidate");
  declareProperty("mcContainerName",m_mcContainerName="GEN_AOD");
  declareProperty ("truthMapName", m_TruthMap = "TrackParticleTruthCollection");


  declareProperty( "PartPropSvc", m_partPropSvc,
		   "Handle to the particle property service" );

  // Declare the properties of the tool:
  declareProperty( "PtMin" , m_ptMinCut = 100.,
		   "Minimum track pT to be selected (in MeV)" );
  declareProperty( "EtaMax" , m_etaCut = 5.,
		   "Maximum track eta. Set to negative not to apply the cut" );
  declareProperty( "RemoveEmptyEvents" , m_removeEmptyEvents = true );
  declareProperty( "RemoveDummyEvents" , m_removeDummyEvents = false );
  declareProperty( "RemoveInTimePileUp" , m_removeInTimePileUp = false );
  declareProperty( "Remove2BCPileUp" , m_remove2BCPileUp = false );
  declareProperty( "Remove800nsPileUp" , m_remove800nsPileUp = false );
  declareProperty( "RemoveCavernBkg" , m_removeCavernBkg = false );
  declareProperty( "SelectTruthTracks", m_selectTruthTracks = false,
		   "Only select stable charged particles" );
  declareProperty( "AddOnlyFirstVertex", m_addOnlyFirstVertex=false,
		   "add only first vertex per gen event");


}

VertexTruthMatchAthAlgTool::~VertexTruthMatchAthAlgTool()
{

}

StatusCode VertexTruthMatchAthAlgTool::queryInterface( const InterfaceID& riid, void** ppvIf )
{
   if ( riid == IVertexTruthMatchAthAlgTool::interfaceID() )  {
      *ppvIf = (IVertexTruthMatchAthAlgTool*)this;
      addRef();
      return StatusCode::SUCCESS;
   }

   return AthAlgTool::queryInterface( riid, ppvIf );
}

StatusCode VertexTruthMatchAthAlgTool::initialize() {
  msg(MSG::INFO) << "Initialising VertexTruthMatchAthAlgTool" << endreq;

  IIncidentSvc* incsvc;
  StatusCode sc = service("IncidentSvc", incsvc);
  int priority = 100;
  if( sc.isSuccess() ) {
     incsvc->addListener( this, "BeginEvent", priority);
  }

  if( m_partPropSvc.retrieve().isFailure()){
    msg(MSG::WARNING) << "Failed to get particle properties service." << endreq;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode VertexTruthMatchAthAlgTool::process() {

  // This gets the EventInfo object from StoreGate
  const EventInfo* myEventInfo = 0;
  if(evtStore()->retrieve(myEventInfo/*,"MyEvent"*/).isFailure()) {
    // Key "MyEvent" is optional, usually not specified for EventInfo because
    // there'll be only one. When not specified, just takes the first container.
    msg(MSG::ERROR) << "Failed to retrieve event information" << endreq;
    return StatusCode::FAILURE;
  }

  unsigned int ei_RunNumber = myEventInfo->event_ID()->run_number();
  unsigned int ei_EventNumber =myEventInfo->event_ID()->event_number();

  const McEventCollection* mcCollptr;
  if ( evtStore()->retrieve( mcCollptr, m_mcContainerName ).isFailure() ) {
    ATH_MSG_WARNING( "Could not retrieve McEventCollection" );
    return StatusCode::SUCCESS;
  }    
  ATH_MSG_DEBUG( "HepMC info loaded" );

  const VxContainer* vxContainer;
  if (evtStore()->retrieve(vxContainer,m_vxContainerName).isFailure() ) {
    ATH_MSG_WARNING ("Could not retrieve primary vertex container with key "+m_vxContainerName);
    return StatusCode::SUCCESS;
  }
  ATH_MSG_DEBUG("Vx container loaded");

  const Rec::TrackParticleContainer* trackParticleContainer;
  if (evtStore()->retrieve(trackParticleContainer,m_trackContainerName).isFailure() ) {
    ATH_MSG_WARNING ("Could not retrieve TrackParticleBaseCollection container with key "+m_trackContainerName);
    return StatusCode::SUCCESS;
  }
  ATH_MSG_DEBUG("Trk container loaded");

  const TrackParticleTruthCollection* tm;
  if(evtStore()->retrieve(tm, m_TruthMap).isFailure()) {
    ATH_MSG_WARNING ("Could not retrieve TrackParticleTruthCollection container with key "+m_TruthMap);
    return StatusCode::SUCCESS;
  }

  ATH_MSG_DEBUG("Trk truth loaded");

  //SETUP GEN INFO
  //Set up variables to store generator information to pass to matching tool
  int countevt=0;
  ngenevt=0;
  short pileuptype;
  genpileuptype;
  ngenparticle;

  int countvtx=0;
  ngenvtx=0;
  vtxevtindex.clear();
  genvtxx.clear();
  genvtxy.clear();
  genvtxz.clear();

  ngenpart=0;
  partevtindex.clear();
  partprodvtxindex.clear();
  partmctype.clear();
  partpt.clear();
  parteta.clear();
  partbarcode.clear();
  partstatus.clear();
  std::vector<HepMC::GenParticle*> tmpgencol;

  // loop over all events in McEventCollection
  for ( McEventCollection::const_iterator itr = mcCollptr->begin(); itr != mcCollptr->end(); ++itr ) {
    HepMC::GenEvent *evt = (*itr);
    if(!pass( evt, mcCollptr)) continue;

    //determine pileuptype
    pileuptype=-1;
    if( ( evt->event_number() == -1 ) && ( evt->signal_process_id() == 0 ) ) {
      pileuptype = 5;
    } else if (itr == mcCollptr->begin()) {
      pileuptype=0;//HS
    } else {
      pileuptype=1;
      for(McEventCollection::const_iterator itr2 = mcCollptr->begin(); itr2 != itr; ++itr2 ) {
	if( ( ( ( *itr2 )->event_number() == -1 ) &&
	      ( ( *itr2 )->signal_process_id() == 0 ) ) ) {
	  ++pileuptype;
	}
      }
    }

    genpileuptype.push_back(pileuptype);
    countevt++;


    HepMC::GenEvent::vertex_const_iterator vitr;
    for (vitr = evt->vertices_begin(); vitr != evt->vertices_end(); ++vitr ) {
      HepMC::GenVertex *vtx = (*vitr);
      if(!pass(vtx,mcCollptr)) continue;

      countvtx++;
      vtxevtindex.push_back(countevt-1);
      genvtxx.push_back( vtx->point3d().x() );
      genvtxy.push_back( vtx->point3d().y() );
      genvtxz.push_back( vtx->point3d().z() );


    }

    HepMC::GenEvent::particle_const_iterator pitr;

    for (pitr = evt->particles_begin(); pitr != evt->particles_end(); ++pitr ) {
      HepMC::GenParticle *part = (*pitr);

      if(!pass(part, mcCollptr)) continue;
      

	partevtindex.push_back(countevt-1);
	//	partprodvtxindex.push_back(countvtx-1);
	partmctype.push_back(part->pdg_id());
	partpt.push_back(part->momentum().perp());
	parteta.push_back(part->momentum().eta());
	partbarcode.push_back(part->barcode());
	partstatus.push_back(part->status());
	tmpgencol.push_back(part);

    }
  
    ngenparticle.push_back( evt->particles_size() );

  }
   
  ngenevt=countevt;
  ngenvtx=countvtx;
  ngenpart=partpt.size();


  //SETUP TRACKS
  //variables to use
  ntrks=0;
  trk_pt.clear();
  trk_eta.clear();
  trk_d0_wrtPV.clear();
  trk_z0_wrtPV.clear();
  trk_mc_probability.clear();
  trk_mc_index.clear();

  //loop over all tracks
  for(Rec::TrackParticleContainer::const_iterator tkit=trackParticleContainer->begin();
      tkit!=trackParticleContainer->end(); tkit++) {
    Rec::TrackParticle * tk = *tkit;
    //count tracks individually - could add cuts here if needed later
    ntrks++;
    trk_pt.push_back( tk->measuredPerigee()->pT() );
    trk_eta.push_back( tk->measuredPerigee()->eta() );
    trk_d0_wrtPV.push_back( tk->measuredPerigee()->parameters()[Trk::d0] );
    trk_z0_wrtPV.push_back( tk->measuredPerigee()->parameters()[Trk::z0] );

    bool genfound=false;
    ElementLink<Rec::TrackParticleContainer> tlink;
    tlink.setElement(tk);
    tlink.setStorableObject(*trackParticleContainer);
    TrackParticleTruthCollection::const_iterator found = tm->find(tlink);
    if(found != tm->end() ) { //found a matching particle

      for(uint i=0; i<tmpgencol.size(); i++) {
	if(tmpgencol.at(i) == found->second.particleLink().cptr()) {
	  //	if(partbarcode.at(i)== found->second.particleLink().barcode() ) {
	  if(genfound)
	    ATH_MSG_WARNING("Double match of gen particles to track");
	  else{
	    //ATH_MSG_DEBUG("Found gen particle");
	    genfound = true;
	    trk_mc_index.push_back(i);
	    trk_mc_probability.push_back(found->second.probability());
	    break;
	  }
	}
      }
    }
    if(!genfound) { //add dummy entry
      trk_mc_index.push_back(-1);
      trk_mc_probability.push_back(0);
    }


  }


  //SETUP VERTICES
  //variables to use
  nvtx=0;
  vx_trk_n.clear();
  vx_trk_weight.clear();
  vx_trk_index.clear();
  vx_x.clear();
  vx_y.clear();
  vx_z.clear();
  vx_cov_x.clear();
  vx_cov_y.clear();
  vx_cov_z.clear();


  //loop over all vertices
  //ATH_MSG_DEBUG("Start loop over vertices");
  for(VxContainer::const_iterator vit=vxContainer->begin(); vit!=vxContainer->end(); vit++) {
    Trk::VxCandidate *v = *vit;

    nvtx++;
    vx_trk_n.push_back( v->vxTrackAtVertex()->size() );
    
    vx_x.push_back( v->recVertex().position().x() );
    vx_y.push_back( v->recVertex().position().y() );
    vx_z.push_back( v->recVertex().position().z() );

    vx_cov_x.push_back( v->recVertex().errorPosition().covValue( Trk::x ) );
    vx_cov_y.push_back( v->recVertex().errorPosition().covValue( Trk::y ) );
    vx_cov_z.push_back( v->recVertex().errorPosition().covValue( Trk::z ) );

    //loop over tracks in vertex
    //ATH_MSG_DEBUG("Loop over linked tracks");
    std::vector<float> tmp;
    std::vector<int> tmp2;
    for(std::vector<Trk::VxTrackAtVertex*>::const_iterator vtit = v->vxTrackAtVertex()->begin();
	vtit!=v->vxTrackAtVertex()->end(); vtit++) {

      Trk::VxTrackAtVertex* vt = *vtit;
      tmp.push_back( vt->weight() );

      //find the corresponding regular track
      const Trk::LinkToTrackParticleBase * linkToTrackParticle = dynamic_cast<const Trk::LinkToTrackParticleBase *>( vt->trackOrParticleLink() ); 
      bool tkfound=false;
      if(linkToTrackParticle) {
	//ATH_MSG_DEBUG("Look for matching track");
	for(int i=0; i<ntrks; i++) {
	  const Rec::TrackParticle * trk = trackParticleContainer->at(i);
	  const Trk::TrackParticleBase *particle_base = **linkToTrackParticle;
	  if(  trk == particle_base ) {
	    tkfound=true;
	    tmp2.push_back( i );
	    break;
	    //ATH_MSG_DEBUG("track pointer check worked");
	  }
	}
      }
      if(!tkfound)
	tmp2.push_back(-1);

    }
    vx_trk_weight.push_back(tmp);
    vx_trk_index.push_back(tmp2);
  }

  //DEBUG stuff about what is being passed to truth matching algorithm
  msg(MSG::DEBUG) << "Start truth matching with this information:" << endreq;
  msg(MSG::DEBUG) << "Num GenEvents = " << ngenevt << endreq;
  for(uint i=0; i<genpileuptype.size();i++) {
    msg(MSG::DEBUG) << "    Evt "<< i << "- type " << genpileuptype.at(i) << ", " << ngenparticle.at(i) << " particles" << endreq;

    //temporary just to compare 1 event in detail
    if(i==1 || i==4 || i==10 || i==19 ) {
      for(int j=0; j<ntrks; j++) {
	if(trk_mc_index.at(j)!= -1 && partevtindex.at(trk_mc_index.at(j)) == i) 
	  msg(MSG::DEBUG) << "        Track " << j << " pt = " << trk_pt.at(j) << ", prob = " << trk_mc_probability.at(j) << endreq;
      }
    }


    if(i>=80)
      break;
  }
  msg(MSG::DEBUG) << "Num GenVertex = " << ngenvtx << endreq;
  for(uint i=0; i<ngenvtx; i++) {

    msg(MSG::DEBUG) << "    Vtx " << i << " - gen evt " << vtxevtindex.at(i) << " = (" << genvtxx.at(i) << "," << genvtxy.at(i) << "," << genvtxz.at(i) << ")" << endreq;

  }

  msg(MSG::DEBUG) << "Num GenPart   = " << ngenpart << endreq;

  msg(MSG::DEBUG) << "Num RecoVertex= " << nvtx << " " << endreq;
  for(uint i=0; i<nvtx; i++) {

    msg(MSG::DEBUG) << "    Vtx " << i << " Ntrk = " << vx_trk_n.at(i) << ", x = (" << vx_x.at(i) << "," << vx_y.at(i) << "," << vx_z.at(i) << ")" << endreq;


    for(uint j=0; j< vx_trk_index.at(i).size(); j++)
      msg(MSG::DEBUG) << "        Track index " << vx_trk_index.at(i).at(j) << endreq;


  }


  msg(MSG::DEBUG) << "Num Tracks    = " << ntrks << endreq;


  //Setup truth matching algo
  m_vtxtm.SetRecoVtxInfo(nvtx, 
			 &vx_trk_n, &vx_trk_weight, &vx_trk_index);//,
// 		      &vx_x, &vx_y, &vx_z,
// 		      &vx_cov_x, &vx_cov_y, &vx_cov_z);
  m_vtxtm.SetRecoTrkInfo(ntrks,
		      &trk_mc_probability, &trk_mc_index,
			 &trk_pt);//, &trk_eta, &trk_d0_wrtPV, &trk_z0_wrtPV);
  m_vtxtm.SetGenPartInfo(ngenpart,
			 &partevtindex, &partmctype,
			 &partpt, &parteta, &partbarcode, &partstatus);
  m_vtxtm.SetGenVtxInfo(ngenvtx,&vtxevtindex);//,
// 			&genvtxx, &genvtxy, &genvtxz);
  m_vtxtm.SetGenEventsInfo(ngenevt, &genpileuptype, &ngenparticle);

  m_vtxtm.MatchVertices(ei_RunNumber, ei_EventNumber);

  std::vector<VertexMatchInfo> matchedVtx = m_vtxtm.matchedVtx;

//   msg(MSG::DEBUG) << matchedVtx.size() << " matches" << endreq;
//   for(uint i=0; i<matchedVtx.size(); i++)
//     msg(MSG::DEBUG) << "vtx " << i << " type = " << matchedVtx[i].GetType() << endreq;

  ATH_MSG_DEBUG("End of process");

  return StatusCode::SUCCESS;
}

StatusCode VertexTruthMatchAthAlgTool::finalize(){
  msg(MSG::INFO) << "Finalising VertexTruthMatchAthAlgTool." << endreq;
  return StatusCode::SUCCESS;
}

void VertexTruthMatchAthAlgTool::handle(const Incident& inc) {
  // Get the messaging service, print where you are

  msg(MSG::DEBUG) << "entering handle(), incidence type " << inc.type()
		  << " from " << inc.source() << endreq;

  // Only call fillIOV for EndEvent incident
  if (inc.type() != "BeginEvent") return;
  if(1) //for now just run every time w/o setup
    if(process().isFailure()) msg(MSG::FATAL) << "Process of VtxTm failed" << endreq;

  msg(MSG::DEBUG) << "end event handle" << endreq;

}

bool VertexTruthMatchAthAlgTool::pass( const HepMC::GenEvent* evt,
                                 const McEventCollection* coll ) const {

   bool isEmpty = ( evt->particles_size() == 0 );
   bool isDummy = ( ( evt->event_number() == -1 ) &&
                    ( evt->signal_process_id() == 0 ) );
   if( isDummy ) isEmpty = false;

   if( m_removeEmptyEvents && isEmpty ) return false;
   if( m_removeDummyEvents && isDummy ) return false;

   // We can't do the further selection without the full truth record:
   if( ! coll ) return true;

   McEventCollection::const_iterator iter = coll->begin();
   McEventCollection::const_iterator end = coll->end();

   if( *iter == evt ) return true; /// always keep signal

   int gotzero = 1;
   for( ; iter != end; ++iter ) {
      if( ( ( ( *iter )->event_number() == -1 ) &&
            ( ( *iter )->signal_process_id() == 0 ) ) ) {
         ++gotzero;
      }
      if( evt == *iter ) break;
   }

   if( ( gotzero == 1 ) && m_removeInTimePileUp ) return false;
   if( ( gotzero == 2 ) && m_remove2BCPileUp ) return false;
   if( ( gotzero == 3 ) && m_remove800nsPileUp ) return false;
   if( ( gotzero == 4 ) && m_removeCavernBkg ) return false;

   return true;
}

bool VertexTruthMatchAthAlgTool::pass( const HepMC::GenParticle* part,
                                 const McEventCollection* coll ) const {

   // Check if the particle is coming from a "good" GenEvent:
   if( ! pass( part->parent_event(), coll ) ) return false;

   // Execute the 4-momenum cuts:
   const HepMC::FourVector& p4 = part->momentum();
   double pt = p4.perp();
   double eta = p4.eta();
   if( pt < m_ptMinCut ) return false;
   if( ( std::abs( eta ) > m_etaCut ) && ( m_etaCut >= 0. ) ) return false;

   // If we don't want to specifically select charged truth tracks, then this
   // is already good enough:
   if( ! m_selectTruthTracks ) return true;

   if (part->barcode() < 200000) {
     if( ! TruthHelper::IsGenStable()( part ) ) return false;
     if( ! TruthHelper::IsGenInteracting()( part ) ) return false;
   }


   int pdg = part->pdg_id();

   /// remove gluons and quarks of status 2 that pass IsGenStable!!!
   if( abs(pdg) < 7 || abs(pdg) == 21 ) return false;

   const HepPDT::ParticleData* pd = m_partPropSvc->PDT()->particle( abs( pdg ) );
   if( ! pd ) {
     ATH_MSG_DEBUG( "Could not get particle data for pdg = " << pdg 
                      << " status " << part->status() << " barcode " <<part->barcode()
                      << " process id " <<part->parent_event()->signal_process_id());
      return false;
   }
   float charge = pd->charge();

   if( std::abs( charge ) < 1E-5 ) return false;

   return true;
}

bool VertexTruthMatchAthAlgTool::pass( const HepMC::GenVertex* vtx,
                                 const McEventCollection* coll ) const {

  HepMC::GenEvent* event = vtx->parent_event();

   // Check if the vertex is coming from a "good" GenEvent:
   if( ! pass(event , coll ) ) return false;

   /// always add first vertex in event
   HepMC::GenEvent::vertex_const_iterator vtxfirst =  event->vertices_begin(); 
   ///can't be invalid - already have vertices in event
   if(*vtxfirst == vtx){
     return true;
   }
   
   if(m_addOnlyFirstVertex) return false;

   // Check if any of the incoming particles pass our selection cuts:
   std::vector< HepMC::GenParticle* >::const_iterator iter =
      vtx->particles_in_const_begin();
   std::vector< HepMC::GenParticle* >::const_iterator end =
      vtx->particles_in_const_end();
   for( ; iter != end; ++iter ) {
      if( pass( *iter ) ) return true;
   }

   // Check if any of the outgoing particles pass our selection cuts:
   iter = vtx->particles_out_const_begin();
   end = vtx->particles_out_const_end();
   for( ; iter != end; ++iter ) {
      if( pass( *iter ) ) return true;
   }

   // If there aren't any good particles associated with the vertex, then
   // it shouldn't be selected:
   return false;
}
