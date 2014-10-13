// Store global settings to be shared by all components
// S. Pagan Griso <spagangriso@lbl.gov>
// TODO: Move to run-time loading of configuration from file

#ifndef GlobalSettings_h
#define GlobalSettings_h

#include <TString.h>
#include <string>

/** class to store global settings.
 * Will allow public members for settings, for easy usage/modification.
 *
 * Available sections:
 *
 * - Overall path/file names settings.
 * Also specifies the expected directory structure. Fields in {} refer to settings defined just below,
 * while fields in [] are set automatically inside this program by loop indexes or command-line arguments,
 * with usually obvious meaning (if not, specified). Fields after a ':' mean objects stored in the ROOT file
 * specified before the ':' itself.
 *
 *
 */
class GlobalSettings {
  public:
    /** Init global settings.
    * @param inputConfigFile input configuration file, if not set use defaults hard-coded
    */
    GlobalSettings(std::string inputConfigFile="");

    /// Destructor, nothing special to be done here.
    ~GlobalSettings();

    //--- Global constants, enums, etc..
    enum AnalysisMethods {
      kVtxC, ///< Vertex-counting (default)
      kEvtC, ///< Event-counting
      kUnfC, ///< Unfolding NVtx
      kNAnalysisMethods
    };

    // --- Overall path/file names settings
    static const std::string getPrefix( const TString &runNumber );

    /** Raw counts from D3PDData. Assume the following structure:
     *  {path_inputRawCount}/[run]/[settings]/{path_inputRawCount_v}{fname_inputRawCount_Tree/Histo}
     */
    static const std::string path_inputRawCount; ///< path for output tree and histograms
    static const std::string path_inputRawCount_v; ///<optional version (include "/" if directory)
    static const std::string fname_inputRawCount_Tree; ///< filename
    static const std::string fname_inputRawCount_Histo; ///< filename

    /** Storage location of pileup corrections
      */
    static const std::string path_D3PDMCResults;
    static const std::string path_D3PDMCResults_v;
    static const std::string path_fakeCorrection;
    static const std::string path_maskingCorrection;


    /** Unfolding settings
     *  {path_unf_ResponseMatrix}/{fname_unf_ResponseMatrix}:{other_unf_ResponseMatrixHNameBase}_NTrk[nTrkCut]
     */
    static const std::string path_unf_ResponseMatrix;
    static const std::string fname_unf_ResponseMatrix;
    static const std::string other_unf_ResponseMatrixHNameBase; ///<_NTrkX will be appended

    /** Output paths
     *  {path_outputVdM}/[run]/[settings]/[systematics]/N[method]/{fname_outputVdM}
     */
    static const std::string path_outputVdM;
    static const std::string fname_outputVdM;
    static const std::string fname_outputDebugVdM;

    /** Luminosity measurement results
      * {path_lumiVtx}/[run]/[settings]/lumi_r[run].root
      */
    static const std::string path_lumiVtx;

    /** VdM external information: timestamps
     * {path_timestamps}/[run]/scan/[x1,x2,...,y1,y2,...]_timestamps.dat
     */
    static const std::string path_timestamps;

    /** External information: deadtime corrections
     * {path_deadtime}/dt_[run].root
     */
    static const std::string path_deadtime;

    /** External information: Official lumi ntuples from COOL info
    * {path_lumiNtuples}/[lumiNtupleName]
    * where [lumiNtupleName] is the standard one for the corresponding scan:
    *  May '11: 2011MayScan1Raw-v8.root
    *  Apr '12: 2012AprScan1Raw-v10.root and 2012AprScan2Raw-v10.root
    */
    static const std::string path_lumiNtuples;

    /** External information: Official lumi ntuples for physics runs
     * {path_PhysRun_lumiNtuples}/r[Run].root
     */
    static const std::string path_PhysRun_lumiNtuples;

    /** Prefixes
     */
    static const std::string path_VdM_prefix;
    static const std::string path_mu_prefix;
	static const std::string path_data_prefix;

    /** External information: Luminosity file comparisons.
     * {path_lumiTxt}/{fname_lumiBcmH}   //vdM calibrations
     * {path_lumiTxt}/{fname_lumiBcmV}   //vdM calibrations
     * {path_lumiTxt}/{fname_muScanAll}  //muScan 2011 results
     * {path_lumiTxt}/{fname_r191373}  //Results for run 191373
     */
    static const std::string path_lumiTxt;
    static const std::string fname_lumiBcmH;
    static const std::string fname_lumiBcmV;
    static const std::string fname_muScanAll;
    static const std::string fname_r191373;

    /** Skip-list of pLB for various special runs. */
    static const std::string path_plb_skiplist;

};

#endif
