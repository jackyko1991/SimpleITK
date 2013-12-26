#include "sitkImageRegistrationMethod.h"

#include "sitkCreateInterpolator.hxx"
#include "itkImage.h"
#include "itkImageRegistrationMethod.h"


#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkMeanReciprocalSquareDifferenceImageToImageMetric.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMatchCardinalityImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkKullbackLeiblerCompareHistogramImageToImageMetric.h"
#include "itkMeanSquaresHistogramImageToImageMetric.h"


#include "itkRegularStepGradientDescentBaseOptimizer.h"
#include "itkGradientDescentOptimizer.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkExhaustiveOptimizer.h"
#include "itkSingleValuedNonLinearVnlOptimizer.h"

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

namespace itk
{
namespace simple
{



class CommandIterationUpdate
  : public itk::Command
{
public:
  typedef CommandIterationUpdate   Self;
  typedef itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;

  static Pointer New(itk::simple::ImageRegistrationMethod &sitkIRM)
    {
      Pointer smartPtr = new Self(sitkIRM);
      smartPtr->UnRegister();
      return smartPtr;
    }

protected:
  CommandIterationUpdate(itk::simple::ImageRegistrationMethod &sitkIRM)
    : m_IterationNumber(0), m_sitkIRM(sitkIRM)
    {}

  unsigned int m_IterationNumber;
  itk::simple::ImageRegistrationMethod &m_sitkIRM;

public:

  typedef itk::SingleValuedNonLinearOptimizer      OptimizerType;
  typedef const OptimizerType*                     OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer =
                         dynamic_cast< OptimizerPointer >( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    ++m_IterationNumber;


    const itk::GradientDescentOptimizer *gdOpt
      = dynamic_cast<const itk::GradientDescentOptimizer *>(optimizer);
    const itk::RegularStepGradientDescentBaseOptimizer *rsgdOpt
      = dynamic_cast<const itk::RegularStepGradientDescentBaseOptimizer *>(optimizer);
    const itk::OnePlusOneEvolutionaryOptimizer *opoeOpt
      = dynamic_cast<const itk::OnePlusOneEvolutionaryOptimizer *>(optimizer);
    const itk::ExhaustiveOptimizer *eOpt
      = dynamic_cast<const itk::ExhaustiveOptimizer *>(optimizer);
    const itk::SingleValuedNonLinearVnlOptimizer *vnlOpt
      = dynamic_cast<const itk::SingleValuedNonLinearVnlOptimizer *>(optimizer);

    if ( gdOpt )
      {
      std::cout << gdOpt->GetCurrentIteration() << " = ";
      std::cout << gdOpt->GetValue() << " : ";
      m_sitkIRM.m_MetricValue = gdOpt->GetValue();
      m_sitkIRM.m_Iteration = gdOpt->GetCurrentIteration();
      }
    else if ( rsgdOpt )
      {
      std::cout << rsgdOpt->GetCurrentIteration() << " = ";
      std::cout << rsgdOpt->GetValue() << " : ";
      std::cout << rsgdOpt->GetCurrentStepLength() << "   ";
      m_sitkIRM.m_MetricValue =  rsgdOpt->GetValue();
      m_sitkIRM.m_Iteration = rsgdOpt->GetCurrentIteration();
      }
    else if( opoeOpt )
      {
      std::cout << opoeOpt->GetCurrentIteration() << "   ";
      std::cout << opoeOpt->GetValue() << " : ";
      std::cout << opoeOpt->GetFrobeniusNorm() << "   ";
      m_sitkIRM.m_MetricValue =  opoeOpt->GetValue();
      m_sitkIRM.m_Iteration = opoeOpt->GetCurrentIteration();
      }
    else if ( eOpt )
      {
      std::cout << eOpt->GetCurrentIndex() << "   ";
      std::cout << eOpt->GetCurrentValue() << " : ";
      m_sitkIRM.m_MetricValue =  eOpt->GetCurrentValue();
      m_sitkIRM.m_Iteration = m_IterationNumber;
      }
    else if ( vnlOpt )
      {
      std::cout << m_IterationNumber << "   ";
      std::cout << vnlOpt->GetCachedValue() << " :  ";
      std::cout << vnlOpt->GetCachedCurrentPosition() << std::endl;
      m_sitkIRM.m_MetricValue =  vnlOpt->GetCachedValue();
      m_sitkIRM.m_Iteration = m_IterationNumber;
      return;
      }
    else
      {
      std::cout << m_IterationNumber << " = ";
      m_sitkIRM.m_Iteration = m_IterationNumber;
      }

    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }

};


  ImageRegistrationMethod::ImageRegistrationMethod()
    : m_ActiveOptimizer(NULL)
  {
    m_MemberFactory.reset( new  detail::MemberFunctionFactory<MemberFunctionType>( this ) );

    m_MemberFactory->RegisterMemberFunctions< BasicPixelIDTypeList, 3 > ();
    m_MemberFactory->RegisterMemberFunctions< BasicPixelIDTypeList, 2 > ();

    //m_MemberFactory->RegisterMemberFunctions< RealPixelIDTypeList, 3 > ();
    //m_MemberFactory->RegisterMemberFunctions< RealPixelIDTypeList, 2 > ();


    m_Interpolator = sitkLinear;
    m_OptimizerMinimize = true;
  }


  ImageRegistrationMethod::~ImageRegistrationMethod()
  {
  }

  std::string  ImageRegistrationMethod::ToString() const
  {
    std::ostringstream out;
    out << "itk::simple" << this->GetName() << std::endl;

    out << "  Interpolator: ";
    this->ToStringHelper(out, this->m_Interpolator);
    out << std::endl;

    out << "  Transform: ";
    this->ToStringHelper(out, this->m_Transform.ToString());
    out << std::endl;

    out << "  FixedImageRegionSize: ";
    this->ToStringHelper(out, this->m_FixedImageRegionSize);
    out << std::endl;

     out << "  FixedImageRegionIndex: ";
    this->ToStringHelper(out, this->m_FixedImageRegionIndex);
    out << std::endl;

    return out.str();
  }


  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetMetricAsMeanSquares( uint64_t numberOfSpatialSamples )
  {
    m_MetricType = MeanSquares;
    m_MetricNumberOfSpatialSamples = numberOfSpatialSamples;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod:: SetMetricAsNormalizedCorrelation( bool subtractMean )
  {
    m_MetricType = NormalizedCorrelation;
    m_MetricSubtractMean = subtractMean;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetMetricAsMeanReciprocalSquareDifference( double _lambda, double delta )
  {
    m_MetricType = MeanReciprocalSquareDifference;
    m_MetricLambda = _lambda;
    m_MetricDelta = delta;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetMetricAsMutualInformation( double fixedImageStandardDeviation,
                                                          double movingImageStandardDeviation,
                                                          uint64_t numberOfSpatialSamples)
  {
   m_MetricType = MutualInformation;
   m_MetricFixedImageStandardDeviation = fixedImageStandardDeviation;
   m_MetricMovingImageStandardDeviation = movingImageStandardDeviation;
   m_MetricNumberOfSpatialSamples = numberOfSpatialSamples;
   return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetMetricAsMattesMutualInformation( unsigned int numberOfHistogramBins,
                                                               bool useExplicitPDFDerivatives,
                                                               uint64_t numberOfSpatialSamples  )
  {
    m_MetricType = MattesMutualInformation;
    m_MetricNumberOfHistogramBins = numberOfHistogramBins;
    m_MetricUseExplicitPDFDerivatives = useExplicitPDFDerivatives;
    m_MetricNumberOfSpatialSamples = numberOfSpatialSamples;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetMetricAsMatchCardinality(  bool measureMatches,
                                                        uint64_t numberOfSpatialSamples  )
  {
    m_MetricType = MatchCardinality;
    m_MetricMeasureMatches = measureMatches;
    m_MetricNumberOfSpatialSamples = numberOfSpatialSamples;
    return *this;
  }



  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetMetricAsKullbackLeiblerCompareHistogram( double epsilon,
                                                                       const std::vector<unsigned int> &histogramSize )
  {
    m_MetricType = KullbackLeiblerCompareHistogram;
    m_MetricEpsilon = epsilon;
    m_MetricHistogramSize = histogramSize;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetMetricAsNormalizedMutualInformationHistogram( const std::vector<unsigned int> &histogramSize )
  {
    m_MetricType = NormalizedMutualInformationHistogram;
    m_MetricHistogramSize = histogramSize;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetMetricAsMeanSquaresHistogram( const std::vector<unsigned int> &histogramSize )
  {
    m_MetricType = MeanSquaresHistogram;
    m_MetricHistogramSize = histogramSize;
    return *this;
  }


  template <class TImageType>
  itk::ImageToImageMetric<TImageType,TImageType>*
  ImageRegistrationMethod::CreateMetric( )
  {
    typedef TImageType     FixedImageType;
    typedef TImageType     MovingImageType;


    switch (m_MetricType)
      {
      case MeanSquares:
      {
        typedef itk::MeanSquaresImageToImageMetric< FixedImageType, MovingImageType >    _MetricType;
        typename _MetricType::Pointer         metric        = _MetricType::New();
        metric->UseAllPixelsOn();
        if ( this->m_MetricNumberOfSpatialSamples )
          {
          metric->SetNumberOfSpatialSamples( this->m_MetricNumberOfSpatialSamples );
          }
        this->m_OptimizerMinimize = true;

        metric->Register();
        return metric.GetPointer();
      }
      case NormalizedCorrelation:
      {
        typedef itk::NormalizedCorrelationImageToImageMetric< FixedImageType, MovingImageType >    _MetricType;
        typename _MetricType::Pointer         metric        = _MetricType::New();
        metric->SetSubtractMean(this->m_MetricSubtractMean);

        this->m_OptimizerMinimize = true;

        metric->Register();
        return metric.GetPointer();
        break;
      }
      case MeanReciprocalSquareDifference:
      {
        typedef itk::MeanReciprocalSquareDifferenceImageToImageMetric< FixedImageType, MovingImageType >    _MetricType;
        typename _MetricType::Pointer         metric        = _MetricType::New();
        metric->SetLambda(this->m_MetricLambda);
        metric->SetDelta(this->m_MetricDelta);

        this->m_OptimizerMinimize = true;

        metric->Register();
        return metric.GetPointer();
        break;
      }

      case MutualInformation:
      {
        typedef itk::MutualInformationImageToImageMetric< FixedImageType, MovingImageType >    _MetricType;
        typename _MetricType::Pointer         metric        = _MetricType::New();
        metric->SetFixedImageStandardDeviation( this->m_MetricFixedImageStandardDeviation );
        metric->SetMovingImageStandardDeviation( this->m_MetricMovingImageStandardDeviation );
        if ( this->m_MetricNumberOfSpatialSamples )
          {
          metric->SetNumberOfSpatialSamples( this->m_MetricNumberOfSpatialSamples );
          }

        this->m_OptimizerMinimize = false;

        metric->Register();
        return metric.GetPointer();
        break;
      }
      case MattesMutualInformation:
      {
        typedef itk::MattesMutualInformationImageToImageMetric< FixedImageType, MovingImageType >    _MetricType;
        typename _MetricType::Pointer         metric        = _MetricType::New();
        metric->SetNumberOfHistogramBins( this->m_MetricNumberOfHistogramBins );
        metric->SetUseExplicitPDFDerivatives( this->m_MetricUseExplicitPDFDerivatives );
        metric->UseAllPixelsOn();
        if ( this->m_MetricNumberOfSpatialSamples )
          {
          metric->SetNumberOfSpatialSamples( this->m_MetricNumberOfSpatialSamples );
          }

        this->m_OptimizerMinimize = true;

        metric->Register();
        return metric.GetPointer();
        break;
      }
      case MatchCardinality:
      {
        typedef itk::MatchCardinalityImageToImageMetric< FixedImageType, MovingImageType >    _MetricType;
        typename _MetricType::Pointer         metric        = _MetricType::New();
        metric->SetMeasureMatches(m_MetricMeasureMatches);
        metric->UseAllPixelsOn();
        if ( this->m_MetricNumberOfSpatialSamples )
          {
          metric->SetNumberOfSpatialSamples( this->m_MetricNumberOfSpatialSamples );
          }

        this->m_OptimizerMinimize = !m_MetricMeasureMatches;

        metric->Register();
        return metric.GetPointer();
        break;
      }
      case KullbackLeiblerCompareHistogram:
      {
        typedef itk::KullbackLeiblerCompareHistogramImageToImageMetric< FixedImageType, MovingImageType >    _MetricType;
        typename _MetricType::Pointer         metric        = _MetricType::New();
        metric->SetEpsilon(this->m_MetricEpsilon);
        metric->SetHistogramSize(sitkSTLVectorToITKArray<typename _MetricType::HistogramSizeType::ValueType>(this->m_MetricHistogramSize));
        this->m_OptimizerMinimize = true;

        metric->Register();
        return metric.GetPointer();
        break;
      }
      case NormalizedMutualInformationHistogram:
      {
        typedef itk::NormalizedMutualInformationHistogramImageToImageMetric< FixedImageType, MovingImageType >    _MetricType;
        typename _MetricType::Pointer         metric        = _MetricType::New();
        metric->SetHistogramSize(sitkSTLVectorToITKArray<typename _MetricType::HistogramSizeType::ValueType>(this->m_MetricHistogramSize));
        this->m_OptimizerMinimize = false;

        metric->Register();
        return metric.GetPointer();
        break;
      }
      case MeanSquaresHistogram:
      {
        typedef itk::MeanSquaresHistogramImageToImageMetric< FixedImageType, MovingImageType >    _MetricType;
        typename _MetricType::Pointer         metric        = _MetricType::New();
        metric->SetHistogramSize(sitkSTLVectorToITKArray<typename _MetricType::HistogramSizeType::ValueType>(this->m_MetricHistogramSize));
        this->m_OptimizerMinimize = true;

        metric->Register();
        return metric.GetPointer();
        break;
      }
      default:
        sitkExceptionMacro("LogicError: Unexpected case!");
      }

 }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetOptimizerAsRegularStepGradientDescent( double maxStep,
                                                                     double minStep,
                                                                     unsigned int numberOfIteratons,
                                                                     double relaxationFactor )
  {
    m_OptimizerType = RegularStepGradientDescent;
    m_OptimizerMaximumStepLength = maxStep;
    m_OptimizerMinimumStepLength = minStep;
    m_OptimizerNumberOfIterations = numberOfIteratons;
    m_OptimizerRelaxationFactor = relaxationFactor;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetOptimizerAsGradientDescent( double learningRate, unsigned int numberOfIteratons )
  {
    m_OptimizerType = GradientDescent;
    m_OptimizerLearningRate = learningRate;
    m_OptimizerNumberOfIterations = numberOfIteratons;

    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetOptimizerAsConjugateGradient( )
  {
    m_OptimizerType = ConjugateGradient;

    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetOptimizerAsOnePlusOneEvolutionary( double initialRadius,
                                                                 double epsilon,
                                                                 unsigned int numberOfIteratons,
                                                                 double growthFactor,
                                                                 double shrinkFactor )
  {
    m_OptimizerType = OnePlusOneEvolutionary;
    m_OptimizerInitialRadius = initialRadius;
    m_OptimizerEpsilon = epsilon;
    m_OptimizerNumberOfIterations = numberOfIteratons;
    m_OptimizerGrowthFactor = growthFactor;
    m_OptimizerShrinkFactor = shrinkFactor;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetOptimizerAsExhaustive( double stepLength,
                                                     const std::vector<unsigned int> &numberOfSteps )
  {
    m_OptimizerType = Exhaustive;
    m_OptimizerStepLength = stepLength;
    m_OptimizerNumberOfSteps = numberOfSteps;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetOptimizerAsAmoeba( double simplexDelta,
                                                 double parametersConvergenceTolerance,
                                                 double functionConvergenceTolerance,
                                                 unsigned int numberOfIterations)
  {
    m_OptimizerType = Amoeba;
    m_OptimizerSimplexDelta = simplexDelta;
    m_OptimizerParametersConvergenceTolerance = parametersConvergenceTolerance;
    m_OptimizerFunctionConvergenceTolerance = functionConvergenceTolerance;
    m_OptimizerNumberOfIterations = numberOfIterations;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetOptimizerAsLBFGS( double defaultStepLength,
                                                double gradientConvergenceTolerance,
                                                double lineSearchAccuracy,
                                                unsigned int maximumNumberOfFunctionEvaluations )
  {
    m_OptimizerType = LBFGS;
    m_OptimizerDefaultStepLength = defaultStepLength;
    m_OptimizerGradientConvergenceTolerance = gradientConvergenceTolerance;
    m_OptimizerLineSearchAccuracy = lineSearchAccuracy;
    m_OptimizerMaximumNumberOfFunctionEvaluations = maximumNumberOfFunctionEvaluations;

    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetOptimizerScales( const std::vector<double> &scales)
  {
    this->m_OptimizerScales = scales;
    return *this;
  }

  ImageRegistrationMethod::Self&
  ImageRegistrationMethod::SetFixedImageRegion( const std::vector<unsigned int> &size,
                                                    const std::vector<unsigned int> &index)
  {
    this->m_FixedImageRegionSize = size;
    this->m_FixedImageRegionIndex = index;
    return *this;
  }


  std::string ImageRegistrationMethod::GetOptimizerStopConditionDescription() const
  {
    return this->m_StopConditionDescription;
  }

  unsigned int ImageRegistrationMethod::GetOptimizerIteration() const
  {
    return this->m_Iteration;
  }


  std::vector<double> ImageRegistrationMethod::GetOptimizerPosition() const
  {
    if(this->m_ActiveOptimizer==NULL)
      {
      return std::vector<double>();
      }

    typedef itk::SingleValuedNonLinearOptimizer::ParametersType ParametersType;

    const ParametersType &p = this->m_ActiveOptimizer->GetCurrentPosition();
    return std::vector<double>(p.begin(),p.end());
  }

  double ImageRegistrationMethod::GetMetricValue() const
  {
    return this->m_MetricValue;
  }


  Transform ImageRegistrationMethod::Execute ( const Image &fixed, const Image & moving )
  {
    const PixelIDValueType fixedType = fixed.GetPixelIDValue();
    const unsigned int fixedDim = fixed.GetDimension();
    if ( fixed.GetPixelIDValue() != moving.GetPixelIDValue() )
      {
      sitkExceptionMacro ( << "Fixed and moving images must be the same datatype! Got "
                           << fixed.GetPixelIDValue() << " and " << moving.GetPixelIDValue() );
      }

    if ( fixed.GetDimension() != moving.GetDimension() )
      {
      sitkExceptionMacro ( << "Fixed and moving images must be the same dimensionality! Got "
                           << fixed.GetDimension() << " and " << moving.GetDimension() );
      }

    if (this->m_MemberFactory->HasMemberFunction( fixedType, fixedDim ) )
      {
      return this->m_MemberFactory->GetMemberFunction( fixedType, fixedDim )( fixed, moving );
      }

    sitkExceptionMacro( << "Filter does not support fixed image type: " << itk::simple::GetPixelIDValueAsString (fixedType) );

  }

    template<class TImageType>
    Transform ImageRegistrationMethod::ExecuteInternal ( const Image &inFixed, const Image &inMoving )
    {
      typedef TImageType     FixedImageType;
      typedef TImageType     MovingImageType;


      typedef itk::ImageRegistrationMethod<FixedImageType, MovingImageType>  RegistrationType;
      typename RegistrationType::Pointer   registration  = RegistrationType::New();

      // Get the pointer to the ITK image contained in image1
      typename FixedImageType::ConstPointer fixed = this->CastImageToITK<FixedImageType>( inFixed );
      typename MovingImageType::ConstPointer moving = this->CastImageToITK<MovingImageType>( inMoving );

      registration->SetFixedImage( fixed );
      registration->SetMovingImage( moving );

      typedef itk::InterpolateImageFunction< MovingImageType, double > InterpolatorType;
      typename InterpolatorType::Pointer   interpolator  = CreateInterpolator(moving.GetPointer(), m_Interpolator);
      registration->SetInterpolator( interpolator );

      typename itk::ImageToImageMetric<FixedImageType, MovingImageType>::Pointer metric = this->CreateMetric<FixedImageType>();
      registration->SetMetric( metric );
      metric->UnRegister();

      typename itk::SingleValuedNonLinearOptimizer::Pointer optimizer = this->CreateOptimizer();
      optimizer->UnRegister();


      CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New(*this);
      optimizer->AddObserver( itk::IterationEvent(), observer );

      registration->SetOptimizer( optimizer );

      typename RegistrationType::TransformType *itkTx;
      if ( !(itkTx = dynamic_cast<typename RegistrationType::TransformType *>(this->m_Transform.GetITKBase())) )
        {
        sitkExceptionMacro( "Unexpected error converting transform! Possible miss matching dimensions!" );
        }
      registration->SetTransform( itkTx );

      registration->SetInitialTransformParameters(itkTx->GetParameters());

      if (m_FixedImageRegionSize.size() && m_FixedImageRegionIndex.size())
        {
        typedef typename FixedImageType::RegionType RegionType;
        RegionType r( sitkSTLVectorToITK<typename RegionType::IndexType>(m_FixedImageRegionIndex),
                      sitkSTLVectorToITK<typename RegionType::SizeType>(m_FixedImageRegionSize) );
        registration->SetFixedImageRegion( r );
        }


      this->PreUpdate( registration.GetPointer() );
      if (this->GetDebug())
        {
        std::cout << optimizer;
        std::cout << metric;
        }

      registration->Update();

      // update measurements
      this->m_StopConditionDescription = registration->GetOptimizer()->GetStopConditionDescription();

      return this->m_Transform;
    }



}
}
