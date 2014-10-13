#include "InDetTruthVertexValidation/VertexTree.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

  //USAGE : TestTmD3PD [INFILE]


  //will need better argument handling later

  if(argc < 2 ) {
    cerr << "Provide input D3PD file" << endl;
    return 1;
  }

  TFile fin(argv[1]);

  const char * treename = (argc > 2) ? argv[2] : "InDetTrackTree";


  TTree * tr = (TTree*) fin.Get(treename);

  VertexTree * vloop = new VertexTree(tr);
  vloop->Loop();
  cout << "Done loop" << endl;
}
