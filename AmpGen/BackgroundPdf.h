#ifndef AMPGEN_BACKGROUNDPDF_H
#define AMPGEN_BACKGROUNDPDF_H 1 

#include <memory.h>
#include <stddef.h>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/MatrixElement.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/Types.h"
#include "AmpGen/Event.h"
#include "AmpGen/Projection.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/Store.h"
#include "AmpGen/KeyedFunctors.h"

namespace AmpGen { 

  class CompiledExpressionBase;
  class MinuitParameter;
  class MinuitParameterSet;
  class LinearErrorPropagator;
  class FitFraction;
  class Particle;

  class BackgroundPdf
    {
    protected:
      std::vector<MatrixElement> m_matrixElements;
      std::vector<std::complex<double>> m_coefficients;
      Bilinears m_normalisations; //// bilinear normalisation terms ////
      AmplitudeRules m_protoAmplitudes;
      double m_norm; 
      EventList* m_events;
      EventList* m_sim; 
      EventType  m_evtType;
      double m_weight; /// global weight, i.e. the yield ///
      MinuitParameter* m_weightParam; 
      std::vector<unsigned int> m_cacheAddresses; /// the addresses in the event cache for each PDF /// 
      std::vector<unsigned int> m_cacheAddrInt;
      int m_prepareCalls;
      int m_lastPrint;
      int m_printFreq;
      Integrator m_integrator;
      std::string m_prefix;
      bool       m_stateIsGood;
      bool       m_isConstant; 
      bool       m_dbThis; 
      bool       m_verbosity; 
      void addMatrixElement( std::pair<Particle,Coupling>& particleWithCoupling ) ; 
      bool isFixedPDF() const ;
    public:
      using EventList_type  = EventList;
      BackgroundPdf( const EventType& type , 
          AmpGen::MinuitParameterSet& mps ,
          const std::string& prefix="" ) ; 
         BackgroundPdf();
      virtual ~BackgroundPdf(); 
      double operator()( const Event& evt) const {
        return prob(evt);
      }

      AmplitudeRules protoAmplitudes(){ return m_protoAmplitudes ; }
      std::vector<MatrixElement> matrixElements() { return m_matrixElements ; } 

      std::vector<unsigned int> processIndex(const std::string& label) const ;
      std::string getParentProcess( const std::string& label ) const ;

      unsigned int getPdfIndex( const std::string& name ) const ;
      void PConjugate(); 

      MatrixElement operator[]( const unsigned int& index ){
        return m_matrixElements[index]; 
      }

      std::string prefix() const { return m_prefix ; } 

      inline std::complex<double> amplitude( unsigned int i ) const {
        return m_coefficients[i];
      }
      double getWeight() const { return m_weight; }
      void setWeight( const double& weight ){ m_weight = weight ; }
      void setWeight( MinuitParameter* param ){
        m_weightParam = param ; 
      }
      unsigned int size() const { return m_matrixElements.size(); }

      void reset(bool resetEvents=false); 
      void setEvents( EventList_type& list );
      void setMC(EventList_type& sim ); 

      double norm(const Bilinears& norms) const ;
      double norm() const ;
      std::complex<double> norm ( const unsigned int& x, const unsigned int & y ) const  ;
      void transferParameters(); 
      void preprepare();
      void prepare() ;
      void printVal( const Event& evt,bool isSim=false );
      void updateNorms( const std::vector<unsigned int>& changedPdfIndices );
      std::vector<unsigned int> cacheAddresses( const EventList& evts ) const {
        std::vector<unsigned int> addresses;
        //for( auto& mE : m_matrixElements ){
          //addresses.push_back(  evts.getCacheIndex(mE.pdf) );
        //}
        return addresses;
      }

      std::complex<double> getVal( const Event& evt ) const {
        std::complex<double> value(0.,0.);
        //for( unsigned int i=0;i<m_coefficients.size();++i){
          //value += m_coefficients[i]*evt.getCache(m_cacheAddresses[i]);
        //}
        return value;
      }

      std::complex<double> getval( const Event& evt, 
          const std::vector<unsigned int>& cacheaddresses ) const {
        std::complex<double> value(0.,0.);
        //for( unsigned int i=0;i<m_coefficients.size();++i){
        ////  info( m_coefficients[i] << " " << cacheaddresses.size() << " " << evt.cachesize() );
          //value += m_coefficients[i]*evt.getCache(cacheAddresses[i]);
        //}
        return value;
      }
      std::complex<double> getValNoCache( const Event& evt ) ;

      bool isStateGood(){ return m_stateIsGood ; } 
      double prob( const Event& evt ) const {
        //return m_weight * std::norm( getVal(evt) ) / m_norm ; 
        return m_weight * evt[18] ; 
      }
      double prob_unnormalised( const Event& evt ) const {
        //return std::norm( getVal( evt ) ) ;
        return evt[18] ;
      }
      void debug( const unsigned int& N=0, const std::string& nameMustContain="") ; 

      std::vector<FitFraction> fitFractions(const LinearErrorPropagator& linProp); 

      void makeBinary( const std::string& fname , const double& normalisation=1) ; 

      double getNorm( const Bilinears& normalisations); 
      std::function<real_t(const Event&)> evaluator(const EventList_type* = nullptr) const; 
      std::function<complex_t(const Event&)> amplitudeEvaluator(const EventList_type* = nullptr) const; 
      KeyedFunctors<double(Event)> componentEvaluator(const EventList_type* = nullptr) const; 
      void resync() ;  
  }; 
} 

#endif

