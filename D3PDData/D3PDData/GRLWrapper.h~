#ifndef GRL_WRAPPER_H
#define GRL_WRAPPER_H

#include "GoodRunsLists/TGoodRunsListReader.h"
#include "GoodRunsLists/TGoodRunsList.h"

#include <string>

class GRLWrapper {

public:

  GRLWrapper( const std::string &GRLxmlPath );
  ~GRLWrapper();

  bool passedGRL( const unsigned long &runNumber, const unsigned long &lbn ) const;

private:  
  Root::TGoodRunsListReader *m_reader;
  Root::TGoodRunsList *m_list;

};

#endif
