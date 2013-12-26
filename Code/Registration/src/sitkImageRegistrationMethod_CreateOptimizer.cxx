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

}

namespace itk
{
namespace simple
{

  itk::SingleValuedNonLinearOptimizer*
  ImageRegistrationMethod::CreateOptimizer( )
  {
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
      return optimizer.GetPointer();
      }
    if ( m_OptimizerType == ConjugateGradient )
      {
      typedef itk::ConjugateGradientOptimizer _OptimizerType;
      _OptimizerType::Pointer      optimizer     = _OptimizerType::New();
      if (scales.GetSize()) optimizer->SetScales(scales);
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
