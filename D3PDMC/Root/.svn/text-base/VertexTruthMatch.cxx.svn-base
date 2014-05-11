#include "D3PDMC/VertexTruthMatch.h"

using namespace std;

VertexTruthMatch::VertexTruthMatch() {
  m_type = VtxTM_NMatch;
  m_recoVtx = -1;
}

VertexTruthMatch::VertexTruthMatch(int p_index, float p_weight) {
  m_type = VtxTM_NMatch;
  Add(p_index, p_weight);
}

VertexTruthMatch::~VertexTruthMatch() {
}

VertexTruthMatch::VtxTMatch VertexTruthMatch::GetType() {
  return m_type;
}

void VertexTruthMatch::Add(int p_index, float p_weight) {
  m_matchList.push_back(make_pair<int, float>(p_index, p_weight));
}

int VertexTruthMatch::GetNMatched(double threshold) {
  if (m_matchList.size() == 0) {
    return -1;  //should really never happen..
  }
  if (m_matchList[0].first == -1 && m_matchList[0].second == -1) {
    return 0;
  }
  unsigned int numMatchedThreshold=0;
  float totW=0;
  while (totW < threshold) {
    // if, for sum reason, we do not have enough vtx (should never happen), return -1
    if (numMatchedThreshold >= m_matchList.size()) {
      return -1;
    }
    totW += m_matchList[numMatchedThreshold].second;
    numMatchedThreshold++;
  }
  return numMatchedThreshold;
}

