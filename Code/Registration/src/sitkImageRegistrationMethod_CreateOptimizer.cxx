#include "sitkImageRegistrationMethod.h"

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkScaledRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentOptimizer.h"
#include "itkVersorTransformOptimizer.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkConjugateGradientOptimizer.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkNormalVariateGenerator.h"
#include "itkExhaustiveOptimizer.h"
#include "itkAmoebaOptimizer.h"
#include "itkLBFGSOptimizer.h"

#include "itkVersorRigid3DTransform.h"

namespace {

template< typename TValue, typename TType>
itk::Array<TValue> sitkSTLVectorToITKArray( const std::vector< TType > & in )
{
  itk::Array<TValue> out(in.size());
  for( unsigned int i = 0; i < in.size(); ++i )
    {
    out[i] = in[i];
    }
  return out;
}


class OptimizerIterationUpdate
  : public itk::Command
{
public:
  typedef OptimizerIterationUpdate   Self;
  typedef itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;


  itkNewMacro( Self );

protected:
  OptimizerIterationUpdate()
    : m_Iteration(0)
    {}

  unsigned int m_Iteration;

public:

  typedef itk::SingleValuedNonLinearOptimizer      OptimizerType;
  typedef const OptimizerType*                     OptimizerPointer;

  unsigned int GetIteration( ) const
    {
      return this->m_Iteration;
    }

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object *, const itk::EventObject & event)
  {
    if ( itk::StartEvent().CheckEvent( &event ) )
      {
      m_Iteration=0;
      return;
      }

    if( !itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    ++m_Iteration;

  }
};


std::vector<double> _GetOptimizerPosition_Vnl( const itk::SingleValuedNonLinearVnlOptimizer *opt )
{
    typedef itk::SingleValuedNonLinearOptimizer::ParametersType ParametersType;

    const ParametersType &p = opt->GetCachedCurrentPosition();
    return std::vector<double>(p.begin(),p.end());
}

std::vector<double> _GetOptimizerPosition( const itk::SingleValuedNonLinearOptimizer *opt )
{
    typedef itk::SingleValuedNonLinearOptimizer::ParametersType ParametersType;

    const ParametersType &p = opt->GetCurrentPosition();
    return std::vector<double>(p.begin(),p.end());
}


}

namespace itk
{
namespace simple
{



  itk::SingleValuedNonLinearOptimizer*
  ImageRegistrationMethod::CreateOptimizer( )
  {
    OptimizerIterationUpdate::Pointer observer = OptimizerIterationUpdate::New();

    itk::SingleValuedNonLinearOptimizer::ScalesType scales(m_OptimizerScales.size());
    std::copy( m_OptimizerScales.begin(), m_OptimizerScales.end(), scales.begin() );

    if ( m_OptimizerType == GradientDescent )
      {
      typedef itk::GradientDescentOptimizer _OptimizerType;
      _OptimizerType::Pointer      optimizer     = _OptimizerType::New();
      optimizer->SetLearningRate( this->m_OptimizerLearningRate );
      optimizer->SetNumberOfIterations( this->m_OptimizerNumberOfIterations  );
      optimizer->SetMinimize( this->m_OptimizerMinimize );
      if (scales.GetSize()) optimizer->SetScales(scales);

      this->m_pfGetMetricValue = std::tr1::bind(&_OptimizerType::GetValue,optimizer);
      this->m_pfGetOptimizerIteration = std::tr1::bind(&_OptimizerType::GetCurrentIteration,optimizer);
      this->m_pfGetOptimizerPosition = std::tr1::bind(&_GetOptimizerPosition,optimizer);

      optimizer->Register();
      return optimizer.GetPointer();
      }
    else if ( m_OptimizerType == RegularStepGradientDescent )
      {
      typedef itk::RegularStepGradientDescentBaseOptimizer _OptimizerType;
      _OptimizerType::Pointer      optimizer;
      if ( dynamic_cast<const itk::VersorRigid3DTransform<double> *>(this->m_Transform.GetITKBase()) )
        {
        optimizer = VersorRigid3DTransformOptimizer::New();
        }
      else if( dynamic_cast<const itk::VersorTransform<double> *>(this->m_Transform.GetITKBase()) )
        {
        optimizer = VersorTransformOptimizer::New();
        }
      else
        {
        optimizer = itk::ScaledRegularStepGradientDescentOptimizer::New();
        }
      optimizer->SetMaximumStepLength( this->m_OptimizerMaximumStepLength );
      optimizer->SetMinimumStepLength( this->m_OptimizerMinimumStepLength );
      optimizer->SetNumberOfIterations( this->m_OptimizerNumberOfIterations  );
      optimizer->SetMinimize( this->m_OptimizerMinimize );
      optimizer->SetRelaxationFactor( this->m_OptimizerRelaxationFactor );
      if (scales.GetSize()) optimizer->SetScales(scales);
      optimizer->Register();

      this->m_pfGetMetricValue = std::tr1::bind(&_OptimizerType::GetValue,optimizer);
      this->m_pfGetOptimizerIteration = std::tr1::bind(&_OptimizerType::GetCurrentIteration,optimizer);
      this->m_pfGetOptimizerPosition = std::tr1::bind(&_GetOptimizerPosition,optimizer);

      return optimizer.GetPointer();
      }
    if ( m_OptimizerType == ConjugateGradient )
      {
      typedef itk::ConjugateGradientOptimizer _OptimizerType;
      _OptimizerType::Pointer      optimizer     = _OptimizerType::New();
      if (scales.GetSize()) optimizer->SetScales(scales);

      this->m_pfGetMetricValue = std::tr1::bind(&_OptimizerType::GetCachedValue,optimizer);
      optimizer->AddObserver( itk::IterationEvent(), observer );
      optimizer->AddObserver( itk::StartEvent(), observer );
      this->m_pfGetOptimizerIteration = std::tr1::bind(&OptimizerIterationUpdate::GetIteration,observer);
      this->m_pfGetOptimizerPosition = std::tr1::bind(&_GetOptimizerPosition_Vnl,optimizer);

      optimizer->Register();
      return optimizer.GetPointer();
      }

    else if ( m_OptimizerType == OnePlusOneEvolutionary )
      {
      typedef itk::OnePlusOneEvolutionaryOptimizer _OptimizerType;
      _OptimizerType::Pointer      optimizer     = _OptimizerType::New();
      optimizer->SetMaximumIteration( this->m_OptimizerNumberOfIterations  );
      optimizer->SetMinimize( this->m_OptimizerMinimize );
      if (scales.GetSize()) optimizer->SetScales(scales);

      typedef itk::Statistics::NormalVariateGenerator  GeneratorType;
      GeneratorType::Pointer generator = GeneratorType::New();
      generator->Initialize(12345);
      optimizer->SetNormalVariateGenerator( generator );

      optimizer->Initialize( this->m_OptimizerInitialRadius );
      optimizer->SetEpsilon( this->m_OptimizerEpsilon );
      optimizer->SetGrowthFactor( this->m_OptimizerGrowthFactor );
      optimizer->SetShrinkFactor( this->m_OptimizerShrinkFactor );

      this->m_pfGetMetricValue = std::tr1::bind(&_OptimizerType::GetValue,optimizer);
      this->m_pfGetOptimizerIteration = std::tr1::bind(&_OptimizerType::GetCurrentIteration,optimizer);
      this->m_pfGetOptimizerPosition = std::tr1::bind(&_GetOptimizerPosition,optimizer);

      optimizer->Register();
      return optimizer.GetPointer();
      }
    else if ( m_OptimizerType == Exhaustive )
      {
      typedef itk::ExhaustiveOptimizer _OptimizerType;
      _OptimizerType::Pointer      optimizer     = _OptimizerType::New();
      optimizer->SetStepLength( this->m_OptimizerStepLength );
      optimizer->SetNumberOfSteps( sitkSTLVectorToITKArray<_OptimizerType::StepsType::ValueType>(this->m_OptimizerNumberOfSteps));
      if (scales.GetSize()) optimizer->SetScales(scales);

      this->m_pfGetMetricValue = std::tr1::bind(&_OptimizerType::GetCurrentValue,optimizer);
      optimizer->AddObserver( itk::IterationEvent(), observer );
      optimizer->AddObserver( itk::StartEvent(), observer );
      this->m_pfGetOptimizerIteration = std::tr1::bind(&OptimizerIterationUpdate::GetIteration,observer);
      this->m_pfGetOptimizerPosition = std::tr1::bind(&_GetOptimizerPosition,optimizer);

      optimizer->Register();
      return optimizer.GetPointer();
      }
    else if ( m_OptimizerType == Amoeba )
      {
      typedef itk::AmoebaOptimizer _OptimizerType;
      _OptimizerType::Pointer      optimizer     = _OptimizerType::New();
      optimizer->SetMaximumNumberOfIterations( this->m_OptimizerNumberOfIterations  );
      optimizer->SetParametersConvergenceTolerance(this->m_OptimizerParametersConvergenceTolerance);
      optimizer->SetFunctionConvergenceTolerance(this->m_OptimizerFunctionConvergenceTolerance);
      if (scales.GetSize()) optimizer->SetScales(scales);

      _OptimizerType::ParametersType simplexDelta(this->m_Transform.GetParameters().size() );
      simplexDelta.Fill( this->m_OptimizerSimplexDelta );
      optimizer->SetInitialSimplexDelta( simplexDelta );

      this->m_pfGetMetricValue = std::tr1::bind(&_OptimizerType::GetCachedValue,optimizer);
      optimizer->AddObserver( itk::IterationEvent(), observer );
      optimizer->AddObserver( itk::StartEvent(), observer );
      this->m_pfGetOptimizerIteration = std::tr1::bind(&OptimizerIterationUpdate::GetIteration,observer);
      this->m_pfGetOptimizerPosition = std::tr1::bind(&_GetOptimizerPosition_Vnl,optimizer);

      optimizer->Register();
      return optimizer.GetPointer();
      }
    else if ( m_OptimizerType == LBFGS )
      {
      typedef itk::LBFGSOptimizer _OptimizerType;
      _OptimizerType::Pointer      optimizer     = _OptimizerType::New();

       optimizer->SetGradientConvergenceTolerance( this->m_OptimizerGradientConvergenceTolerance );
       optimizer->SetLineSearchAccuracy( this->m_OptimizerLineSearchAccuracy );
       optimizer->SetDefaultStepLength( this->m_OptimizerDefaultStepLength );
       //optimizer->TraceOn();
       optimizer->SetMaximumNumberOfFunctionEvaluations( this->m_OptimizerMaximumNumberOfFunctionEvaluations );

      if (scales.GetSize()) optimizer->SetScales(scales);

      this->m_pfGetMetricValue = std::tr1::bind(&_OptimizerType::GetCachedValue,optimizer);
      optimizer->AddObserver( itk::IterationEvent(), observer );
      optimizer->AddObserver( itk::StartEvent(), observer );
      this->m_pfGetOptimizerIteration = std::tr1::bind(&OptimizerIterationUpdate::GetIteration,observer);
      this->m_pfGetOptimizerPosition = std::tr1::bind(&_GetOptimizerPosition_Vnl,optimizer);

      optimizer->Register();
      return optimizer.GetPointer();
      }
    else
      {
      sitkExceptionMacro("LogicError: Unexpected case!");
      }
  }

}
}
