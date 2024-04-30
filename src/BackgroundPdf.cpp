#include "AmpGen/BackgroundPdf.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ratio>
#include <thread>

#include "AmpGen/CompiledExpression.h"
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Particle.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/CompiledExpressionBase.h"
#include "AmpGen/CompilerWrapper.h"
#include "AmpGen/ThreadPool.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/simd/utils.h"
#include "AmpGen/Array.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AmpGen;
BackgroundPdf::BackgroundPdf() = default; 

BackgroundPdf::BackgroundPdf( const EventType& type, const MinuitParameterSet& mps, const std::string& prefix )
  :   m_evtType  (type)
      , m_printFreq(NamedParameter<size_t>(     "BackgroundPdf::PrintFrequency", 100)  )
      , m_dbThis   (NamedParameter<bool>(       "BackgroundPdf::Debug"         , false))
      , m_verbosity(NamedParameter<bool>(       "BackgroundPdf::Verbosity"     , 0)    )
  , m_objCache (NamedParameter<std::string>("BackgroundPdf::ObjectCache"   ,"")    )
  , m_prefix   (prefix)
      , m_mps(&mps) 
{
}

void BackgroundPdf::prepare()
{
}

void BackgroundPdf::updateNorms()
{
}

void BackgroundPdf::debug( const Event& evt, const std::string& nameMustContain )
{
}

std::vector<FitFraction> BackgroundPdf::fitFractions(const LinearErrorPropagator& linProp)
{
  std::vector<FitFraction> outputFractions;
  return outputFractions;
}

complex_t BackgroundPdf::getValNoCache( const Event& evt ) const
{
  std::complex<double> value(0,0);
  return value; 
}

void BackgroundPdf::reset( bool resetEvents )
{
  m_prepareCalls                                     = 0;
  m_lastPrint                                        = 0;
  if ( resetEvents ){ 
    m_events = nullptr;
    m_integrator = Integrator();
  }
}

void BackgroundPdf::setEvents( const EventList_type& list )
{
  DEBUG( "Setting event list with:" << list.size() << " events for " << this );
  reset();
  for( auto& me : m_matrixElements ){ DEBUG("Registering: " << me.name() ) ; }
  if( m_ownEvents && m_events != nullptr ) delete m_events; 
  m_events = &list;
  m_cache.allocate(m_events, m_matrixElements); 
}


void BackgroundPdf::setMC( const EventList_type& sim )
{
  if ( m_verbosity ) INFO( "Setting norm. event list with:" << sim.size() << " events for " << this );
  reset();
  m_integrator = Integrator( &sim, m_matrixElements );
}

real_t BackgroundPdf::norm() const
{
  std::complex<double> acc(0,0);

  //return norm(m_normalisations);
  return acc.real();
}

real_t BackgroundPdf::norm(const Bilinears& norms) const
{
  complex_t acc(0, 0);
  //for ( size_t i = 0; i < size(); ++i ) {
  // for ( size_t j = 0; j < size(); ++j ) {
       //INFO( i << " " << j << " " << m_matrixElements[i].coefficient * std::conj(m_matrixElements[j].coefficient) << " " <<  ( i > j ? std::conj(norm(j,i)) : norm(i,j) ) );
  //    acc += m_matrixElements[i].coefficient * std::conj(m_matrixElements[j].coefficient)* ( i > j ? std::conj(norm(j,i)) : norm(i,j) );
    //}
  //}
  return acc.real();
}

complex_t BackgroundPdf::norm(const size_t& x, const size_t& y) const
{
  return m_normalisations.get(x, y);
}

void BackgroundPdf::transferParameters()
{
  //for ( auto& mE : m_matrixElements ) mE.coefficient = mE.coupling();
  //m_weight.update();
}

void BackgroundPdf::printVal(const Event& evt)
{
  /*
     for ( auto& mE : m_matrixElements ) {
     unsigned int address = std::distance( &mE , &m_matrixElements[0] );
     std::cout << mE.decayTree.decayDescriptor() << " = " << mE.coefficient << " x " << m_cache( evt.index() / utils::size<real_v>::value, address )
     << " address = " << address << " " << mE( evt ) << std::endl;
     if( mE.coupling.size() != 1 ){
     std::cout << "CouplingConstants: " << std::endl;
     mE.coupling.print();
     std::cout << "================================" << std::endl;
     }
     }
     */
}

complex_t BackgroundPdf::getVal( const Event& evt ) const
{
  complex_v value( 0., 0. );
  //for (unsigned int i = 0 ; i != m_matrixElements.size(); ++i ) {
  //  value += complex_v( m_matrixElements[i].coefficient ) * m_cache(evt.index() / utils::size<real_v>::value, i );
  //}
#if ENABLE_AVX
  //return utils::at(value, evt.index() % utils::size<real_v>::value);
#else 
  return value;
#endif
}

real_v BackgroundPdf::operator()( const real_v* /*evt*/, const unsigned block ) const 
{
  complex_v value( 0., 0. );
  for ( const auto& mE : m_matrixElements ) 
  {
    unsigned address = &mE - &m_matrixElements[0];
    value += complex_v(mE.coefficient) * m_cache(block, address); 
  }
  //return (m_weight/m_norm ) * utils::norm(value);
  return m_weight;
}

#if ENABLE_AVX

//double BackgroundPdf::operator()( const double* /*evt*/, const unsigned block ) const 
//{
//  return operator()((const real_v*)nullptr, block / utils::size<real_v>::value ).at( block % utils::size<real_v>::value );
//}
#endif


std::function<real_t(const Event&)> BackgroundPdf::evaluator(const EventList_type* ievents) const 
{
  auto events = ievents == nullptr ? m_integrator.events<EventList_type>() : ievents;  
  FunctionCache<EventList_type, complex_v, Alignment::AoS> store(events, m_matrixElements);
  for( auto& me : m_matrixElements ) store.update(me);
  std::vector<double> values( events->aligned_size() );
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for( unsigned int block = 0 ; block < events->nBlocks(); ++block )
  {
    utils::store( values.data() + block * utils::size<real_v>::value,  (m_weight * events->bkgPDF(block))  );
  }
  return arrayToFunctor<real_t, typename EventList_type::value_type>(values);
}

std::function<complex_t(const Event&)> BackgroundPdf::amplitudeEvaluator(const EventList_type* ievents) const 
{
  auto events = ievents == nullptr ? m_integrator.events<EventList_type>() : ievents;  
  FunctionCache<EventList_type, complex_v, Alignment::AoS> store(events, m_matrixElements);
  for( auto& me : m_matrixElements ) store.update( me );
  std::vector<complex_t> values( events->aligned_size() );
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for( unsigned int block = 0 ; block < events->nBlocks(); ++block )
  {
    complex_v amp(0.,0.);
    for( unsigned j = 0 ; j != m_matrixElements.size(); ++j ) 
      amp = amp + complex_v(m_matrixElements[j].coefficient) * store(block, j);
    for( unsigned k = 0; k != utils::size<complex_v>::value; ++k )
    {
      values[ block * utils::size<complex_v>::value + k] = utils::at( amp, k ); 
    }
  }
  return arrayToFunctor<complex_t, typename EventList_type::value_type>(values);
}


BackgroundPdf::~BackgroundPdf()

{
  if( m_ownEvents && m_events !=nullptr ) delete m_events; 
}
