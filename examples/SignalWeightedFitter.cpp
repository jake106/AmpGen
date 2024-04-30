#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "AmpGen/Integrator.h"
#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/BackgroundPdf.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Generator.h"
#include "AmpGen/ErrorPropagator.h"
#ifdef _OPENMP
  #include <omp.h>
  #include <thread>
#endif

#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>

using namespace AmpGen;

void randomizeStartingPoint( MinuitParameterSet& mps, TRandom3& rand, double multiplier)
{
  for (auto& param : mps) {
    if ( ! param->isFree() || param->name() == "Px" || param->name() == "Py" || param->name() == "Pz" ) continue;
    double min = param->minInit();
    double max = param->maxInit();
    //double multiplier = 100;
    double new_value = rand.Uniform(param->mean()-multiplier*param->stepInit(),param->mean()+multiplier*param->stepInit());
    //double new_value = rand.Uniform(param->mean()-18*param->stepInit(),param->mean()+18*param->stepInit());
    if( min != 0 && max != 0 )
      new_value = rand.Uniform(min,max);
    param->setInit( new_value );
    param->setCurrentFitVal( new_value );
    INFO( param->name() << "  = " << param->mean() << " " << param->stepInit() );
  }
}

double calcDervivative(Minimiser mini, MinuitParameter* param, EventList current_dataset, SumPDF<EventList_type, CoherentSum&> pdf, double h, size_t j)
{
  double m = param->mean();
  param->setCurrentFitVal( m + h);
  double newLL = mini.FCN();
  double pdf_plus_variation = pdf.operator()(current_dataset[j]);
  param->setCurrentFitVal( m - h);
  newLL = mini.FCN();
  double pdf_minus_variation = pdf.operator()(current_dataset[j]);
  // Reset param to best fit value
  param->setCurrentFitVal(m);
  return (pdf_plus_variation - pdf_minus_variation) / (2*h);
}

void calcAsymptoticCorrectedCovariance(std::vector<EventList_type> data, std::vector<SumPDF<EventList_type, CoherentSum&>> pdfs, MinuitParameterSet& MPS, FitResult* fr, Minimiser mini)
{
  // Step 1 - get number of floated params and covariance matrix (reduced)
  TMatrixD covReduced = fr->getReducedCovariance();
  int nFloated = covReduced.GetNcols();
  //TMatrixD covFull = fr->cov();

  // Step 2 - initialise new cov matrix
  //TMatrixTSym<Double_t> num(nFloated);
  TMatrixD num(nFloated, nFloated);
  for (int k = 0; k < nFloated; k++) {
    for (int l = 0; l < nFloated; l++) {
      num(k, l) = 0.0;
      }
   }

  // Step 3 - get list of floated params and make sure they match with cov matrix index
  int nTot = MPS.size();
  std::vector<int> floatingIndex(nFloated,0);
  int idx = 0;
  for (int p = 0; p < nTot; p++) {
    auto param = MPS.at(p);
    if (param->isFree() == 1) {
      floatingIndex[idx] = p;
      idx += 1;
    }
  }

  // Perform calculation of D matrix (num) here
  for (size_t i = 0; i < data.size(); i++) {
    auto& current_dataset = data[i];
    for (size_t j = 0; j < current_dataset.size(); j++) {
    //for (size_t j = 0; j < 3; j++) {

      double weightSquared = current_dataset.weight(j)*current_dataset.weight(j);
      double pdf_val = pdfs[i].operator()(current_dataset[j]);
      std::vector<double> diffs(nFloated, 0.0);

      for (int k = 0; k < nFloated; k++) {
        // Here we need to populate diffs with correct derivatives
        auto parameter = MPS.at(floatingIndex[k]);

        // Use a Richardson extrapolation with two terms for derivative estimate
        // What is the best step size (h1) to use here??
        double h1 = std::sqrt( covReduced(k, k) );
        //double h1 = parameter->stepInit();
        //double h1 = 0.001;
        double h2 = h1/2;
        double nh1 = calcDervivative(mini, parameter, current_dataset, pdfs[i], h1, j);
        double nh2 = calcDervivative(mini, parameter, current_dataset, pdfs[i], h2, j);
        diffs[k] = diffs[k] + nh2 + (nh2-nh1)/3;
      }
    for (int k = 0; k < nFloated; k++) {
      for (int l = 0; l < nFloated; l++) {
        num(k, l) += weightSquared * diffs[k] * diffs[l] / (pdf_val*pdf_val);
        }
      }
    }
  }
  //num.Similarity(covFull);
  TMatrixD covNew(nFloated, nFloated);
  // Perform matrix multiplication CDC here
  for (int i = 0; i < nFloated; i++) {
    for (int j = 0; j < nFloated; j++) {
      for (int k = 0; k < nFloated; k++) {
        for (int l = 0; l < nFloated; l++) {
          covNew(i,j) += covReduced(i,k) * num(k,l) * covReduced(l,j);
        }
      }
    }
  }

  // Covariance matrix at given index
  std::cout << "Covariance matrix:" << std::endl;
  for (size_t i = 0; i < (size_t)covReduced.GetNrows(); ++i){
    for (size_t j = 0; j < (size_t)covReduced.GetNrows(); ++j) std::cout << covReduced[i][j] << " ";
   std::cout << "\n";
  }
  std::cout << "New covariance matrix" << std::endl;
  for (size_t i = 0; i < (size_t)covNew.GetNrows(); ++i){
    for (size_t j = 0; j < (size_t)covNew.GetNrows(); ++j) std::cout << covNew[i][j] << " ";
   std::cout << "\n";
  }
  //for (size_t i = 0; i < (size_t)num.GetNrows(); ++i){
  //  for (size_t j = 0; j < (size_t)num.GetNrows(); ++j) std::cout << num[i][j] << " ";
  //std::cout << "\n";
  //}

  // Errors are just sqrt of the covariance matrix diagonals!
}

//template <typename PDF> FitResult* doFit( PDF&& pdf, EventList_type& data, EventList_type& mc, MinuitParameterSet& MPS );
FitResult* doFit( std::vector<EventList_type> data, std::vector<EventList_type> mc, MinuitParameterSet& MPS, size_t NBins );

int main( int argc, char* argv[] )
{
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );

  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  const auto datasets        = NamedParameter<std::string>("Datasets","",
      "List of data/simulated samples to fit, in the format \
      \033[3m data[0] sim[0] data[1] sim[1] ... \033[0m. \nIf a simulated sample is specified FLAT, uniformly generated phase-space events are used for integrals ").getVector();
  //std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  //std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "Fitter.log",     "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root",     "Name of the output plot file");
  const size_t NBins   = NamedParameter<size_t>     ("nBins"     , 100         , "Number of bins used for plotting.");
  const std::string weight_branch         = NamedParameter<std::string>("WeightBranch","","Name of branch containing event weights.");
  const std::string mc_weight_branch      = NamedParameter<std::string>("MCWeightBranch","","Name of branch containing event weights.");

  std::string outOptFile = NamedParameter<std::string>("OutputOptionFile", ""  , "Name of output option file updated with the best-fit parameters");
  std::string inOptFile = NamedParameter<std::string>("InputOptionFile", argv[1] , "Name of input option file to use as template for OutputOptionFile");

  // Add named parameter to allow fixing of a given parameter
  const auto scan_string = NamedParameter<std::string>("ScanParameter", "", "Name and values for parameter to overwrite in options file");  
  
  auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()
              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
  auto MCbNames = NamedParameter<std::string>("MCBranches", std::vector<std::string>(),
                                              "List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
  auto pNames = NamedParameter<std::string>("EventType" , ""    
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
  double multiplier = NamedParameter<double>("Multiplier", 100, "Multiplier to apply to step size to randomize starting point");
 
  // [[maybe_unused]]
  size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );
   
  //if( datasets == "" ) FATAL("Must specify input with option " << italic_on << "DataSample" << italic_off );
  if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);

  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

  INFO("Outputs: LogFile: " << logFile << "; Plots: " << plotFile << "; Options: " << outOptFile);

#if ENABLE_AVX
//  if(!idbranch.empty() || !weight_branch.empty() || !mcidbranch.empty() || !mc_weight_branch.empty()){
 //   ERROR("Vectorized version currently not supported when adding extra branches");
 //   return 1;
//  }
#endif
  
#ifdef _OPENMP
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif

  /* A MinuitParameterSet is (unsurprisingly) a set of fit parameters, and can be loaded from 
     the parsed options. For historical reasons, this is referred to as loading it from a "Stream" */
  MinuitParameterSet MPS;
  MPS.loadFromStream();

  //typedef std::vector<std::string> ScanVector;
  //ScanVector scan_vector;
  //boost::split( scan_vector, scan_string, boost::is_any_of(" ") );

  if ( scan_string.size() > 1 ){
  double newval = ::atof(scan_string.getVal(1).c_str());
  MPS[scan_string.getVal(0)]->setVal(newval);
  MPS[scan_string.getVal(0)]->fix();
  }

  if ( NamedParameter<bool>("RandomizeStartingPoint",false) ){
    randomizeStartingPoint(MPS,rndm,multiplier);
    std::cout << "Starting Point Randomized" << std::endl;
  }

  /* An EventType specifies the initial and final state particles as a vector that will be described by the fit. 
     It is typically loaded from the interface parameter EventType. */
  EventType evtType(pNames);
  /* A CoherentSum is the typical amplitude to be used, that is some sum over quasi two-body contributions 
     weighted by an appropriate complex amplitude. The CoherentSum is generated from the couplings described 
     by a set of parameters (in a MinuitParameterSet), and an EventType, which matches these parameters 
     to a given final state and a set of data. A common set of rules can be matched to multiple final states, 
     i.e. to facilitate the analysis of coupled channels. 
     The CoherentSum is only appropriate for decays involving only (pseudo)scalars in the inital / final state, 
     otherwise the sum must also be over initial / final spin states. In this case, as PolarisedSum should be used. 
     See FitterWithPolarisation for an example of this use case.    
  */
  //CoherentSum sig(evtType, MPS);
  
  /* Events are read in from ROOT files. If only the filename and the event type are specified, 
     the file is assumed to be in the specific format that is defined by the event type, 
     unless the branches to load are specified in the user options */
  std::vector<EventList_type> events;
  /* Generate events to normalise the PDF with. This can also be loaded from a file, 
     which will be the case when efficiency variations are included. Default number of normalisation events 
     is 5 million. */
  std::vector<EventList_type> eventsMC;
  
  for(size_t i=0;i < datasets.size() ; i+=2 ){
    events.emplace_back( datasets[i], evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weight_branch)  );
    if( datasets[i+1] == "FLAT" ) eventsMC.emplace_back( Generator<>(evtType, &rndm).generate(2.5e6) );
    else eventsMC.emplace_back( datasets[i+1], evtType, Branches(MCbNames), GetGenPdf(true), WeightBranch(mc_weight_branch) );
  }

  //EventList_type events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weight_branch) );
  //EventList_type eventsMC = intFile == "" ? Generator<>(evtType, &rndm).generate(2.5e6) : EventList_type(intFile, evtType, Branches(MCbNames), GetGenPdf(true), WeightBranch(mc_weight_branch));
  
  //sig.setMC( eventsMC );

  TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
  
  /* Do the fit and return the fit results, which can be written to the log and contains the 
     covariance matrix, fit parameters, and other observables such as fit fractions */
  FitResult* fr = doFit( events, eventsMC, MPS, NBins );
  /* Calculate the `fit fractions` using the signal model and the error propagator (i.e. 
     fit results + covariance matrix) of the fit result, and write them to a file. 
   */

  fr->writeToFile( logFile.c_str() );
  if ( outOptFile != "" ) fr->writeOptions( outOptFile, inOptFile );
  
  output->Close();
}

FitResult* doFit(std::vector<EventList_type> data, std::vector<EventList_type> mc, MinuitParameterSet& MPS, size_t NBins )
{
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();
  /* Minimiser is a general interface to Minuit1/Minuit2, 
     that is constructed from an object that defines an operator() that returns a double 
     (i.e. the likielihood, and a set of MinuitParameters. */
 
  SimFit likelihood;

  bool sig_only; 
  std::vector<CoherentSum> fcs(data.size());
  std::vector<BackgroundPdf> fcb(data.size());

  std::vector<SumPDF<EventList_type, CoherentSum&, BackgroundPdf&>> pdfs;
  std::vector<SumPDF<EventList_type, CoherentSum&>> pdfs_sig;

  if ( NamedParameter<bool>("BkgPDF",false) ) {
    sig_only = false;
    pdfs.reserve(data.size());
    for(size_t i = 0; i < data.size(); ++i){
      fcs[i] = CoherentSum(data[i].eventType(), MPS);
      fcb[i] = BackgroundPdf(data[i].eventType(), MPS);
      fcs[i].setWeight(MPS["fPDF"]);
      fcb[i].setWeight(MPS["fBkg"]);
      fcs[i].setMC(mc[i]);
      fcb[i].setMC(mc[i]);
      pdfs.emplace_back( make_pdf(fcs[i], fcb[i]) );
      pdfs[i].setEvents(data[i]);
      likelihood.add( pdfs[i] );
    }
  } else {
    sig_only = true;
    pdfs_sig.reserve(data.size());
    for(size_t i = 0; i < data.size(); ++i){
      fcs[i] = CoherentSum(data[i].eventType(), MPS);
      pdfs_sig.emplace_back( make_pdf<EventList_type>(fcs[i]) );
      pdfs_sig[i].setEvents(data[i]);
      auto& mci = mc[i];
      for_each( pdfs_sig[i].pdfs(), [&mci](auto& pdf){pdf.setMC(mci);});
      likelihood.add( pdfs_sig[i] );
    }
  }
  
  Minimiser mini( likelihood, &MPS );
  //mini->maxCalls = 200000
  if ( NamedParameter<bool>("Verbose",false) ) { 
    mini.setPrintLevel( PrintLevel::VeryVerbose );
  }
  mini.doFit();
  FitResult* fr = new FitResult(mini);

  //std::ofstream myfile;
  //myfile.open("fit_pdfs.csv");
  //for(size_t i = 0; i < data.size(); ++i){
    //std::cout << "######### PRINTING PDFS GETVAL ###########" << std::endl;

//    auto& di = data[i];
//    for(size_t j = 0; j < di.size(); ++j){
//    auto& dj = di[j];
//    myfile << pdfs[i].operator()(dj);
//    myfile << ", ";
    //std::cout << pdfs[i].operator()(dj) << std::endl;
   // }
  //}
  //myfile.close();

  auto fitFractions = fcs[0].fitFractions( fr->getErrorPropagator() );
  fr->addFractions( fitFractions );

  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );

  //calcAsymptoticCorrectedCovariance(data, pdfs, MPS, fr, mini);
 
  //Parameter error = sqrt(cov[i,i])
  //Unsure how to recalculate, so just make sure cov is correct and get fit errors from there
  

  //auto lhood = pdfs[0].make 


  /* Estimate the chi2 using an adaptive / decision tree based binning, 
     down to a minimum bin population of 15, and add it to the output. */
  
  //Chi2Estimator* chi2;
  //chi2 = new Chi2Estimator( data[0], mc[0], fcs[0] );
  //chi2->writeBinningToFile("chi2_binning.txt");
  //fr->addChi2( chi2->chi2(), chi2->nBins() );
  //fr->print();

  //FitResult(mini).writeToFile(logFile.c_str());
  //TFile* output_plots = TFile::Open( plotFile.c_str(), "RECREATE");
  for( size_t i = 0 ; i < data.size(); ++i )
  {
    INFO("Making figures for sample: " << i << " ...");
    auto evaluator_per_component = ( fcs[0] ).componentEvaluator(&mc[i]);
    for( auto proj : data[i].eventType().defaultProjections(NBins) )
    {
      proj(data[i], PlotOptions::Prefix("Data"+std::to_string(i)), PlotOptions::AutoWrite() );
    }
    for( auto proj : data[i].eventType().defaultProjections(NBins) )
    {

    if( NamedParameter<bool>("AllComponents",false ) ) {
      proj(mc[i], evaluator_per_component, PlotOptions::Prefix("amp"+std::to_string(i)), PlotOptions::Norm(data[i].size()), PlotOptions::AutoWrite() );
      }
    if( sig_only ) {
      proj(mc[i]  , pdfs_sig[i].componentEvaluator(&mc[i]), PlotOptions::Prefix("pdf"+std::to_string(i)), PlotOptions::Norm(data[i].size()), PlotOptions::AutoWrite() );
    } else {
      proj(mc[i]  , pdfs[i].componentEvaluator(&mc[i]), PlotOptions::Prefix("pdf"+std::to_string(i)), PlotOptions::Norm(data[i].size()), PlotOptions::AutoWrite() );
    }
    }    
  }
  return fr;
}
