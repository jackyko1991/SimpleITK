#ifndef __sitkImageRegistrationMethod_h
#define __sitkImageRegistrationMethod_h

#include "sitkRegistration.h"

#include "sitkDetail.h"
#include "sitkImage.h"
#include "sitkPixelIDTokens.h"
#include "sitkMemberFunctionFactory.h"
#include "sitkProcessObject.h"

#include "sitkInterpolator.h"
#include "sitkTransform.h"

#ifndef SWIGPYTHON
#define LAMBDA lambda
#else
#define LAMBDA _lambda
#endif



namespace itk
{

class SingleValuedNonLinearOptimizer;
template<class T, class U> class ImageToImageMetric;

namespace simple
{

  class CommandIterationUpdate;

  class SITKRegistration_EXPORT ImageRegistrationMethod
    : public ProcessObject
  {
  public:

    typedef ImageRegistrationMethod Self;

    ImageRegistrationMethod();
    virtual ~ImageRegistrationMethod();

    std::string GetName() const { return std::string("ImageRegistrationMethod"); }
    std::string ToString() const;


    InterpolatorEnum GetInterpolator()
      { return this->m_Interpolator; }
    Self& SetInterpolator ( InterpolatorEnum Interpolator )
      { this->m_Interpolator = Interpolator; return *this; }

    Self& SetTransform ( const Transform &Transform )
      { this->m_Transform = Transform; return *this; }
    Transform GetTransform()
      { return this->m_Transform; }

    Self& SetMetricAsMeanSquares( uint64_t numberOfSpatialSamples = 0 );
    Self& SetMetricAsNormalizedCorrelation( bool subtractMean = false );
    Self& SetMetricAsMeanReciprocalSquareDifference( double LAMBDA=1.0,
                                                     double delta=0.00011 );
    Self& SetMetricAsMutualInformation( double fixedImageStandardDeviation=0.4,
                                        double movingImageStandardDeviation=0.4,
                                        uint64_t numberOfSpatialSamples = 50 );
    Self& SetMetricAsMattesMutualInformation( unsigned int numberOfHistogramBins = 50,
                                              bool useExplicitPDFDerivatives=true,
                                              uint64_t numberOfSpatialSamples = 0 );
    Self& SetMetricAsMatchCardinality( bool measureMatches = true,
                                       uint64_t numberOfSpatialSamples = 0 );
    Self& SetMetricAsKullbackLeiblerCompareHistogram( double epsilon = 1e-12,
                                                      const std::vector<unsigned int> &histogramSize =  std::vector<unsigned int>(2,256u) );
    Self& SetMetricAsNormalizedMutualInformationHistogram( const std::vector<unsigned int> &histogramSize =  std::vector<unsigned int>(2,256u) );
    Self& SetMetricAsMeanSquaresHistogram( const std::vector<unsigned int> &histogramSize =  std::vector<unsigned int>(2,256u) );


    Self& SetOptimizerAsRegularStepGradientDescent( double maxStep,
                                                    double minStep,
                                                    unsigned int numberOfIterations,
                                                    double relaxationFactor=0.5);
    Self& SetOptimizerAsGradientDescent( double learningRate,
                                         unsigned int numberOfIterations );
    Self& SetOptimizerAsConjugateGradient( );
    Self& SetOptimizerAsOnePlusOneEvolutionary( double initialRadius,
                                                double epsilon=1.5e-4,
                                                unsigned int numberOfIterations=100,
                                                double growthFactor=1.05,
                                                double shrinkFactor=0.9878);
    Self& SetOptimizerAsExhaustive( double stepLength,
                                    const std::vector<unsigned int> &numberOfSteps );
    Self& SetOptimizerAsAmoeba( double simplexDelta,
                                double parametersConvergenceTolerance=1e-8,
                                double functionConvergenceTolerance=1e-4,
                                unsigned int numberOfIterations=100);
    Self& SetOptimizerAsLBFGS( double defaultStepLength,
                               double gradientConvergenceTolerance,
                               double lineSearchAccuracy,
                               unsigned int maximumNumberOfFunctionEvaluations );
    Self& SetOptimizerScales( const std::vector<double> &scales);


    Self& SetFixedImageRegion( const std::vector<unsigned int> &size, const std::vector<unsigned int> &index);
    std::vector<unsigned int> GetFixedImageRegionSize() const { return this->m_FixedImageRegionSize;}
    std::vector<unsigned int> GetFixedImageRegionIndex() const { return this->m_FixedImageRegionIndex;}

    Transform Execute ( const Image &fixed, const Image & moving );

     /**
      * Active measurements which can be obtained during call backs.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
    unsigned int GetOptimizerIteration() const;
    std::vector<double> GetOptimizerPosition() const;
    double GetMetricValue() const;

    /** Measurement updated at the end of execution.
      */
    std::string GetOptimizerStopConditionDescription() const;

  protected:

    template<class TImage>
    Transform ExecuteInternal ( const Image &fixed, const Image &moving );

    itk::SingleValuedNonLinearOptimizer* CreateOptimizer( );

    template <class TImageType>
      itk::ImageToImageMetric<TImageType,TImageType>* CreateMetric( );


  private:

    typedef Transform (ImageRegistrationMethod::*MemberFunctionType)( const Image &fixed, const Image &moving );
    friend struct detail::MemberFunctionAddressor<MemberFunctionType>;
    std::auto_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;

    InterpolatorEnum  m_Interpolator;
    Transform  m_Transform;

    // optimizer
    enum OptimizerType { RegularStepGradientDescent,
                         GradientDescent,
                         ConjugateGradient,
                         OnePlusOneEvolutionary,
                         Exhaustive,
                         Amoeba,
                         LBFGS
    };
    OptimizerType m_OptimizerType;
    double m_OptimizerLearningRate;
    double m_OptimizerMaximumStepLength;
    double m_OptimizerMinimumStepLength;
    unsigned int m_OptimizerNumberOfIterations;
    double m_OptimizerRelaxationFactor;
    bool m_OptimizerMinimize;
    std::vector<double> m_OptimizerScales;
    double m_OptimizerInitialRadius;
    double m_OptimizerEpsilon;
    double m_OptimizerGrowthFactor;
    double m_OptimizerShrinkFactor;
    std::vector<unsigned int> m_OptimizerNumberOfSteps;
    double m_OptimizerStepLength;
    double m_OptimizerSimplexDelta;
    double m_OptimizerParametersConvergenceTolerance;
    double m_OptimizerFunctionConvergenceTolerance;
    double m_OptimizerDefaultStepLength;
    double m_OptimizerGradientConvergenceTolerance;
    double m_OptimizerLineSearchAccuracy;
    unsigned int m_OptimizerMaximumNumberOfFunctionEvaluations;

    // metric
    enum MetricType { MeanSquares,
                      NormalizedCorrelation,
                      MeanReciprocalSquareDifference,
                      MutualInformation,
                      MattesMutualInformation,
                      MatchCardinality,
                      KullbackLeiblerCompareHistogram,
                      NormalizedMutualInformationHistogram,
                      MeanSquaresHistogram
    };
    MetricType m_MetricType;
    uint64_t  m_MetricNumberOfSpatialSamples;
    double m_MetricFixedImageStandardDeviation;
    double m_MetricMovingImageStandardDeviation;
    unsigned int m_MetricNumberOfHistogramBins;
    std::vector<unsigned int> m_MetricHistogramSize;
    bool m_MetricUseExplicitPDFDerivatives;
    bool m_MetricMeasureMatches;
    bool m_MetricSubtractMean;
    double m_MetricLambda;
    double m_MetricDelta;
    double m_MetricEpsilon;
    // metric seed

    std::vector<unsigned int> m_FixedImageRegionSize;
    std::vector<unsigned int> m_FixedImageRegionIndex;


    std::string m_StopConditionDescription;
    double m_MetricValue;
    unsigned int m_Iteration;

    friend class CommandIterationUpdate;

    itk::SingleValuedNonLinearOptimizer *m_ActiveOptimizer;
  };

}
}

#endif // __sitkImageRegistrationMethod_h
