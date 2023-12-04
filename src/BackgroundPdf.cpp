#include "AmpGen/BackgroundPdf.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/EventType.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/ErrorPropagator.h"

/// STL
#include <iomanip>
#include <chrono>
#include <ctime>
#include <fstream>

#ifdef __USE_OPENMP__
  #include <omp.h>
#endif

using namespace AmpGen; 

void BackgroundPdf::addMatrixElement( std::pair<Particle,Coupling>& particleWithCoupling ){

  //auto& protoParticle = particleWithCoupling.first ; 
  //auto& coupling = particleWithCoupling.second; 

  //if( !protoParticle.isStateGood() ){ 
    //ERROR("Decay tree not configured correctly for " << protoParticle.uniqueString() );
    //m_stateIsGood = false;
    //return; 
  //}
  ////if( CPConjugate ) protoParticle.CPConjugateThis();
  ////if( FlavConjugate ) protoParticle.setConj(true);
  //const std::string name = protoParticle.uniqueString();
  //std::vector<DBSYMBOL> dbExpressions;
  //INFO( "Matrix Element = " << name ); 
  //const Expression expression = 
    //protoParticle.getExpression(m_dbThis?&dbExpressions:NULL);    
  //DEBUG("Got expression for this tree");
  //m_matrixElements.emplace_back(
      //std::make_shared<Particle>(protoParticle),
      //coupling,
      //CompiledExpression<std::complex<double>>( 
        //expression , 
        //name , 
        //m_evtType.getEventFormat(), 
        //m_dbThis?&dbExpressions:NULL ) );
}

BackgroundPdf::BackgroundPdf( const EventType& type , 
    AmpGen::MinuitParameterSet& mps,
    const std::string& prefix ) : 
  m_protoAmplitudes(mps),
  m_events(0),
  m_sim(0),
  m_evtType(type),
  m_weight(1),
  m_weightParam(nullptr),
  m_prepareCalls(0), 
  m_lastPrint(0),
  m_printFreq(0), 
  m_prefix(prefix),
  m_stateIsGood(true)
{
  //bool useCartesian = AmpGen::NamedParameter<unsigned int>("BackgroundPdf::UseCartesian",true);
  //m_printFreq       = AmpGen::NamedParameter<unsigned int>("BackgroundPdf::PrintFrequency",100);
  //m_dbThis          = AmpGen::NamedParameter<unsigned int>("BackgroundPdf::Debug",false);
  //m_verbosity       = AmpGen::NamedParameter<unsigned int>("BackgroundPdf::Verbosity",0);
  //auto rules = m_protoAmplitudes.rulesForDecay( type.mother() );
  //for( auto& p : rules ){
    //if( p.prefix() != m_prefix  ) continue ; 
    //std::vector<std::pair<AmpGen::Particle,Coupling > > tmpParticles;
    //auto fs = type.finalStates();
    //tmpParticles.emplace_back( AmpGen::Particle( p.name(), fs), p.makeCoupling(useCartesian ) );
    //do { 
      //std::vector<std::pair<AmpGen::Particle,Coupling >> newTmpParticles; 
      //for( auto& particleWithCoupling : tmpParticles ){
        //auto protoParticle = particleWithCoupling.first;
        //auto coupling = particleWithCoupling.second; 
        //auto protoFinalStates = protoParticle.getFinalStateParticles();
        //if( protoFinalStates.size() == type.size() ){
          //addMatrixElement( particleWithCoupling );
          //continue ; /// this particle is fully expanded 
        //}
        //std::string nameToExpand = protoParticle.uniqueString(); 
        //for( auto& ifs : protoFinalStates ){
          //auto expandedRules = 
            //m_protoAmplitudes.rulesForDecay( ifs->name() ); /// get rules for decaying particle 
          //if( expandedRules.size() == 0 ) continue ; 
          //for( auto& subTree : expandedRules ){
            //auto expanded_amplitude = replaceAll( nameToExpand , ifs->name(), subTree.name() );
            //auto fs2 = type.finalStates();
            //newTmpParticles.emplace_back( 
                //AmpGen::Particle( expanded_amplitude, fs2 ), Coupling(coupling, subTree, useCartesian) ) ; 
          //}
          //break;  // we should only break if there are rules to be expanded ... 
        //}
      //}
      //tmpParticles = newTmpParticles; 
    //} while ( tmpParticles.size() != 0 ) ; 
  //}
  //for( auto& p : m_matrixElements ){
    //p.pdf.resolveParameters( mps );
    //m_lib.add( &p.pdf );
  //}
  //m_isConstant = isFixedPDF();
  //m_coefficients.resize( m_matrixElements.size() );
  //m_normalisations.resize( m_matrixElements.size() , m_matrixElements.size() );
}

void BackgroundPdf::prepare(){
  if( m_weightParam != 0 ) m_weight = m_weightParam->mean();
  //if( m_isConstant && m_prepareCalls != 0 ) return ;
  //if( m_prepareCalls == 0 ) resync(); // ensure pointers are syncronised //  
  //preprepare();
  //std::vector<unsigned int> changedPdfIndices;
  //auto tStartEval = std::chrono::high_resolution_clock::now();
  //bool printed=false;

  //for( unsigned int i = 0 ; i < m_matrixElements.size(); ++i){
    //auto& pdf = m_matrixElements[i].pdf;
    //if( m_prepareCalls != 0 && ! pdf.hasExternalsChanged() ) continue ;
    //auto t_start = std::chrono::high_resolution_clock::now();
    //if( m_events != 0 ){
      //if( i >= m_cacheAddresses.size() ) 
        //m_cacheAddresses.push_back( m_events->registerExpression(pdf) );
      //m_events->updateCache( pdf, m_cacheAddresses[i] );
    //}
    //else if( i == 0 && m_verbosity ){
      //WARNING("No data events specified for " << this);
    //}
    //auto t_end = std::chrono::high_resolution_clock::now();
    //auto time  = std::chrono::duration<double, std::milli>(t_end-t_start).count() ;
    //if( m_verbosity && ( m_prepareCalls > m_lastPrint + m_printFreq || m_prepareCalls == 0 ) ){
      //INFO(pdf.name() << " (t = " << time << " ms, nCalls = " << m_prepareCalls << ")" );    
      //printed=true;
    //}
    //changedPdfIndices.push_back( i ); 
    //pdf.resetExternals();
  //}
  //auto tStartIntegral = std::chrono::high_resolution_clock::now();
  
  //if( m_sim != 0 ){
    //updateNorms(changedPdfIndices); 
  //}
  //else if( m_verbosity ) { 
    //WARNING( "No simulated sample specified for " << this );
  //}

  //if( m_verbosity && printed ){
    //auto tNow   = std::chrono::high_resolution_clock::now();
    //double timeEval  = std::chrono::duration<double, std::milli>(tStartIntegral-tStartEval).count() ;
    //double timeIntg  = std::chrono::duration<double, std::milli>(tNow-tStartIntegral).count() ;
    //double timeTotal = std::chrono::duration<double, std::milli>(tNow-tStartEval).count() ;
    //INFO( "Time Performance: " 
        //<< "Eval = "     << timeEval   << " ms"
        //<< ", Integral = " << timeIntg   << " ms"
        //<< ", Total = "    << timeTotal  << " ms" );
    //m_lastPrint=m_prepareCalls; 
  //}
  //m_norm = norm(); /// update normalisation 
  //m_prepareCalls++; 
//}

//void BackgroundPdf::updateNorms( const std::vector<unsigned int>& changedPdfIndices ){
  //std::vector<bool> integralHasChanged( size()*size() ,0);
 //// INFO( "PDFs to update : " << vectorToString( changedPdfIndices ) ); 
  //for( auto& i : changedPdfIndices ){
    //auto& pdf = m_matrixElements[i].pdf;
    //m_integralDispatch.prepareExpression( pdf );
  //}

  //for( auto& i : changedPdfIndices ){
    //for( unsigned int j = 0 ; j < size(); ++j){
      //if( integralHasChanged[i*size()+j] ) continue ; 
      //integralHasChanged[i*size()+j] = true;
      //integralHasChanged[j*size()+i] = true;  
      //m_integralDispatch.addIntegral( m_matrixElements[i].pdf , m_matrixElements[j].pdf ,
          //[i,j,this]( const std::complex<double>& val ){
          ////m_integralDispatch.addIntegral(i,j);
          //DEBUG( i << " , " << j << " " << " = " << val );
          //this->m_normalisations.set( i , j , val);
          //this->m_normalisations.set( j , i , std::conj( val ) );
          //} );
    //}
  //}
  //m_integralDispatch.flush(); /// compute remaining integrals in the buffer /// 
}

void BackgroundPdf::debug( const unsigned int& N, const std::string& nameMustContain){ 
  //for( auto& pdf : m_matrixElements ) pdf.pdf.resetExternals();
  //if( nameMustContain == "" ) 
    //for( auto& pdf : m_matrixElements ){
      //pdf.pdf.debug( m_events->getEvent(N) );
    //}
  //else 
    //for( auto& pdf: m_matrixElements ) 
      //if( pdf.pdf.name().find( nameMustContain ) != std::string::npos ) 
        //pdf.pdf.debug( m_events->getEvent(N) );

  //prepare();
  //INFO( "Pdf = " << prob( m_events->at(N) )); 
}

std::vector<FitFraction> BackgroundPdf::fitFractions(
    const LinearErrorPropagator&  linProp )
{
  //struct processCalculator {
    //std::vector<FFCalculator> calculators; 
    //std::vector<FitFraction>    fractions;
    //FitFraction                       sumV;
    //std::string                      name; 
    //double sum()  {
      //double F =0;
      //for( auto& calc : calculators ) F+=calc();
      //return F;
    //};
    
    //size_t size() const { return calculators.size() + 1; }
  //};

  std::vector<FitFraction> outputFractions ; 
  //std::vector<processCalculator> AllCalculators( m_protoAmplitudes.rules().size() );

  //size_t counter=0;
  //size_t     pos=0;
  //for( auto& processes : m_protoAmplitudes.rules() ){
    //auto& pCalc = AllCalculators[counter++];
    //pCalc.name  = processes.first;
    //for( auto& process : processes.second ){  
      //if( process.head() == m_evtType.mother() && process.prefix() != m_prefix ) continue ; 
      //std::string parentProcessName = getParentProcess(process.name());
      //if( parentProcessName == "" ) continue; 
      //std::string pName, ppName; 
      //auto pIndex = processIndex(process.name());
      //auto parentIndex = processIndex( parentProcessName );
      //pCalc.calculators.emplace_back( process.name(), this, pIndex, parentIndex );
    //}
    //pos += pCalc.size();
  //}

  //bool hardcore = NamedParameter<bool>("Hardcore",false);
  //bool interference = NamedParameter<bool>("Interference",false);
  //auto FitFractions = [this,&AllCalculators,&hardcore](){
    //if( hardcore ) this->prepare();
    //else this->transferParameters();
    //std::vector<double> rv; 
    //for( auto& pCalc : AllCalculators ){
      //for( auto& calc : pCalc.calculators ){
        //rv.push_back( calc() );
      //}
      //rv.push_back( pCalc.sum() );
    //}
    //return rv;
  //};

  //auto values = FitFractions();
  //auto errors = linProp.getVectorError( FitFractions, values.size() );
  //counter=0;
  //for( auto& pCalc : AllCalculators ){
    //for( auto& calc : pCalc.calculators ){
      //pCalc.fractions.emplace_back( calc.name, values[counter], errors[counter] );
      //counter++; 
    //}
    //std::sort ( pCalc.fractions.begin(), pCalc.fractions.end(), []( 
        //const FitFraction& f1, const FitFraction& f2){ return fabs(f1.val()) > fabs(f2.val()) ; } );
    //for( auto& f : pCalc.fractions ) 
      //outputFractions.push_back(f);
    //pCalc.sumV = FitFraction( "Sum_"+pCalc.name, values[counter], errors[counter] );
    //outputFractions.push_back( pCalc.sumV );  
    //counter++; 
  //}
  //INFO("Calculating interference fractions");
  //if( hardcore && interference ){ 
  //std::vector<FitFraction> interferenceFractions ; 
  //auto ffForHead = m_protoAmplitudes.rulesForDecay( m_evtType.mother()  );
  
  //for( unsigned int i = 0 ; i < ffForHead.size(); ++i ){
    //auto process_i = ffForHead[i];
    //if( process_i.prefix() != m_prefix ) continue ; 
    //for( unsigned int j=i+1; j < ffForHead.size(); ++j ){
      //auto process_j = ffForHead[j];
      //if( process_j.prefix() != m_prefix ) continue;
      //std::string parent_process_name = getParentProcess(process_i.name());
      //FFCalculator iCalc( process_i.name() + "x"+process_j.name(),
          //this,
          //processIndex( process_i.name() ),
          //processIndex( process_j.name() ),
          //processIndex( parent_process_name )); 
      //interferenceFractions.emplace_back( iCalc.name,iCalc() , linProp.getError( [this,&iCalc](){ this->prepare() ; return iCalc(); } ) );
    //}
  //}
  //std::sort ( interferenceFractions.begin(), interferenceFractions.end(), []( 
        //const FitFraction& f1, const FitFraction& f2){ return fabs(f1.val()) > fabs(f2.val()) ; } );
  //for( auto& f : interferenceFractions ) outputFractions.push_back( f );
  //}
  //for( auto& p : outputFractions ){
    //INFO( std::setw(100) << p.name() << " " << std::setw(5) << round( p.val()*100,3) << " ± " << round( p.err()*100,3) << " %" );
  //}
  return outputFractions; 
}

void BackgroundPdf::makeBinary( const std::string& fname, const double& normalisation ){
  //std::ofstream stream( fname );
  //stream << "#include <complex>" << std::endl;
  //stream << "#include <vector>" << std::endl; 
  //for( auto& p : m_matrixElements ) p.pdf.compile( stream );
  //transferParameters();
  //stream << std::setprecision(10) ; 
  //for( auto& p : m_matrixElements ) p.pdf.compileWithParameters( stream );

  //stream << "extern \"C\" double FCN( double* E , const int& parity){" << std::endl;
  //stream << " std::complex<double> amplitude = " << std::endl;
  //for( unsigned int i = 0 ; i < size(); ++i ){
    //auto& p = m_matrixElements[i];
    //int parity = p.decayTree->finalStateParity();
    //if( parity == -1 ) stream << " double(parity) * ";
    //stream << "std::complex<double>"<< p.coupling() <<" * ";
    //stream << "r" << p.pdf.hash() << "( E )";
    //stream << ( i==size()-1 ? ";" : "+" ) << std::endl;  
  //}
  //stream << " return std::norm(amplitude) / "<< normalisation << " ; }" << std::endl; 
  //stream.close();
}


std::vector<unsigned int> BackgroundPdf::processIndex(const std::string& label) const { 
  std::vector<unsigned int> indices ;
  //for( unsigned int i = 0 ; i < m_matrixElements.size();++i ){
    //bool couplingIncludesThis = false ; 
    //for( auto& reAndIm : m_matrixElements[i].coupling.couplings )
      //if( reAndIm.first->name().find( label) != std::string::npos ) 
        //couplingIncludesThis = true ; 
    //if( couplingIncludesThis ) indices.push_back(i);
  //}
  return indices; 
}

std::string BackgroundPdf::getParentProcess( const std::string& label ) const {
  //auto pI = processIndex(label);
  //if( pI.size() == 0 ) return "";
  //auto coupling = m_matrixElements[pI[0]].coupling.couplings;
  //for( unsigned int i = 0 ; i < coupling.size(); ++i ){
    //if( coupling[i].first->name().find(label) != std::string::npos ){
      //return i==0 ? m_evtType.mother() : coupling[i-1].first->name();
    //}
  //}
  return "";
}


unsigned int BackgroundPdf::getPdfIndex( const std::string& name ) const {
  //for( unsigned int i = 0 ; i < size() ; ++i ){
    //if( m_matrixElements[i].decayTree->uniqueString() == name ) return i;
  //}
  //ERROR("Component " << name << " not found");
  return 999; 
}

bool BackgroundPdf::isFixedPDF() const {
  //for( auto& p : m_matrixElements ){
    //auto params = p.getDependencies();
    //for( auto& param : params )
      //if( param->iFixInit() == 0 ) return false; 
  //}
  return true; 
}

void BackgroundPdf::PConjugate(){
  //for( auto& amp : m_matrixElements ){
    //if( amp.decayTree->finalStateParity() == - 1 ){
      //auto& top_coupling = *amp.coupling.couplings.begin();
      //top_coupling.second->setCurrentFitVal( top_coupling.second->mean() + M_PI );
    //} 
  //}  
}

std::complex<double> BackgroundPdf::getValNoCache( const Event& evt  ) {
  std::complex<double> value(0,0);
  //for( unsigned int i = 0 ; i < m_matrixElements.size(); ++i ){
    //value += m_coefficients[i] * m_matrixElements[i].pdf(evt);
  //}
  return value ; 
}

void BackgroundPdf::preprepare(){
  //if( ! m_lib.isReady() ){
    
    //std::string libname = NamedParameter<std::string>("Lib::" + m_prefix+
        //m_matrixElements[0].decayTree->name() , "");
    
    //INFO( libname );
    //m_lib.link( FCNLibrary::RECOMPILE   | ( m_dbThis & FCNLibrary::DEBUG )  , libname ); 
  //}
  //transferParameters();

  //for( auto& mE : m_matrixElements ){
    //mE.pdf.prepare();
  //}
}


void BackgroundPdf::reset(bool resetEvents){
  m_prepareCalls=0;
  m_lastPrint=0;
  m_cacheAddresses.clear();
  if( resetEvents ){ 
    m_events = 0;
    m_sim = 0 ; 
  };
}
void BackgroundPdf::setEvents( EventList& list ){
  if( m_verbosity ) INFO("Setting events to size = " << list.size() << " for " << this ); 
  reset();
  m_events=&(list); 
  //m_cache.allocate(m_events, m_matrixElements);
}

void BackgroundPdf::setMC(EventList& sim ){
  if( m_verbosity ) INFO("Setting MC = " << &sim << " for " << this ); 
  reset();
  m_integrator = Integrator( &sim, m_matrixElements );
}

double BackgroundPdf::norm() const {
  std::complex<double> acc(0,0);
  //for( unsigned int i=0;i<m_coefficients.size();++i){
    //for( unsigned int j=0;j<m_coefficients.size();++j){
      //auto val = m_normalisations.get(i,j)*m_coefficients[i]*std::conj(m_coefficients[j]);
      //acc += val ;
    //}
  //}
  return acc.real();
}

double BackgroundPdf::norm(const Bilinears& norms ) const {
  std::complex<double> acc(0,0);
  //for( unsigned int i=0;i<m_coefficients.size();++i){
    //for( unsigned int j=0;j<m_coefficients.size();++j){
      //auto val = norms.get(i,j)*m_coefficients[i]*std::conj(m_coefficients[j]);
      //acc += val ;
    //}
  //}
  return acc.real();
}

std::complex<double> BackgroundPdf::norm ( const unsigned int& x, const unsigned int & y ) const { 
  return m_normalisations.get(x,y) ; 
}

void BackgroundPdf::transferParameters(){
  //for( unsigned int i=0;i<m_matrixElements.size();++i ){
    //m_coefficients[i]= m_matrixElements[i].coupling() ;
  //}
}
void BackgroundPdf::printVal( const Event& evt,bool isSim ){
  //for( unsigned int i = 0 ; i < m_coefficients.size(); ++i){
    //unsigned int address = m_cacheAddresses[i];

    //std::cout <<  m_matrixElements[i].decayTree->uniqueString() 
      //<< " = " << m_coefficients[i] << " x " << 
      //evt.getCache(address) << " address = " << address << " " << m_matrixElements[i].pdf(evt)  << std::endl; 
    //std::cout << "Couplings: " << std::endl; 
    //for( auto& ci : m_matrixElements[i].coupling.couplings ){
      //std::cout << ci.first->name() << " " << ci.first->mean() << " " << ci.second->mean() << std::endl; 
    //}

       //std::cout << "================================" << std::endl;  
  //}
}

void BackgroundPdf::resync() {
  //for( auto& m : m_matrixElements ){
    //m_lib.updatePtr( &m.pdf );
  //}
}

