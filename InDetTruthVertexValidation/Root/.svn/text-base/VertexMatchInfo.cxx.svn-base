#include "InDetTruthVertexValidation/VertexMatchInfo.h"

using namespace std;

VertexMatchInfo::VertexMatchInfo() 
{ 
  m_type = VtxTM_NMatch; m_recoVtx = -1;
}

VertexMatchInfo::VertexMatchInfo(int p_index, float p_weight, int ntrk) 
{
  m_type = VtxTM_NMatch;
  Add(p_index, p_weight, ntrk);
}

VertexMatchInfo::~VertexMatchInfo()
{
}

VertexMatchInfo::VtxTMatch VertexMatchInfo::GetType() 
{
  return m_type;
}

void VertexMatchInfo::Add(int p_index, float p_weight, int ntrk) 
{
  m_matchList.push_back(make_pair<int, float>(p_index, p_weight));
  m_ntrk.push_back( ntrk );
}

int VertexMatchInfo::GetNMatched(double threshold) 
{
  if (m_matchList.size() == 0)
    return -1; //should really never happen..
  if (m_matchList[0].first == -1 && m_matchList[0].second == -1)
    return 0;
  int numMatchedThreshold=0;
  float totW=0;
  while (totW < threshold) {
    // if, for sum reason, we do not have enough vtx (should never happen), return -1
    if (numMatchedThreshold >= m_matchList.size())
      return -1;
    totW += m_matchList[numMatchedThreshold].second;
    numMatchedThreshold++;
  }
  return numMatchedThreshold;
}

