/** @file Vertex truth matched EDM */

#ifndef VertexTruthMatch_h
#define VertexTruthMatch_h

#include <vector>

/// Describes one reconstructed vertex with the truth-matching to generated (or fake) vertices
class VertexTruthMatch {
  public:
    enum VtxTMatch {
      VtxTM_Match,
      VtxTM_Merge,
      VtxTM_Split,
      VtxTM_Fake,
      VtxTM_Others,
      VtxTM_NMatch //no valid matching
    };

    /** Index of matched generated vertex.
     *  Fraction of weight of these tracks w.r.t total weight of the vertex
     */
    typedef std::pair<int, float> GenVtxMatch;

    /** Store matched genVertices.
     * sorted by the weight they carry on the reconstructed vertex
     */
    std::vector<GenVtxMatch> m_matchList;
    /** keep track of reconstructed vertex index to which we refer to (safety) */

    int m_recoVtx;
    VtxTMatch m_type;

  public:
    VertexTruthMatch();
    VertexTruthMatch(int p_index, float p_weight);
    ~VertexTruthMatch();

    VtxTMatch GetType();

    void Add(int p_index, float p_weight);

    /// Return number of matched gen vertices until their sum of weights reaches the threshold
    int GetNMatched(double threshold);
};


#endif
