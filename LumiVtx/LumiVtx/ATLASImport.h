//Small helper namespace with useful conversion units

namespace Units {
  const Double_t GeV = 1./1000; //MeV -> GeV

}

namespace Trk {
  enum VertexType {
    NoVtx   = 0,      //!< Dummy vertex, TrackParticle was not used in vertex fit
    PriVtx  = 1,      //!< Primary Vertex
    SecVtx  = 2,      //!< Secondary Vertex
    PileUp  = 3,      //!< Pile Up Vertex
    ConvVtx = 4,      //!< Converstion Vertex
    V0Vtx   = 5,      //!< Vertex from V0 Decay
    KinkVtx = 6,      //!< Kink Vertex
    V0Lambda = 7,     //!< Temporary addition for V0 Lambda
    V0LambdaBar = 8,  //!< Temporary addition for V0 LambdaBar
    V0KShort = 9,     //!< Temporary addition for KShort
    NotSpecified = -99 //!< this is the default
  };
}
