#include "D3PDData/GRLWrapper.h"
#include <iostream>

GRLWrapper::GRLWrapper( const std::string &GRLxmlPath ) : m_reader(0), m_list(0) {
  m_reader = new Root::TGoodRunsListReader();
  std::cout << "GoodRunsLists: using file " << GRLxmlPath << std::endl;
  m_reader->SetXMLFile(GRLxmlPath.c_str());
  m_reader->Interpret();
  std::cout << "Line 9" << std::endl;
  m_list = new Root::TGoodRunsList((m_reader->GetMergedGoodRunsList()));
  m_list->Summary(false);
}

GRLWrapper::~GRLWrapper() {
  if( m_reader ) { delete m_reader; }
  if( m_list ) { delete m_list; }
}

bool GRLWrapper::passedGRL( const unsigned long &runNumber, const unsigned long &lbn ) const {
  return m_list->HasRunLumiBlock(runNumber, lbn);
}
