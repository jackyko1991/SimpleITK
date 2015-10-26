## Needs to be of a form appropriate for ctest

library(SimpleITK)

# arguments after --args
args <- commandArgs( TRUE )
a1 <- '/gpfs/M2Home/rbeare/Monash002_scratch/Devel/SITK/SimpleITK-build/ExternalData/Testing/Data/Input/BrainProtonDensitySliceBorder20.png'
a2 <- '/gpfs/M2Home/rbeare/Monash002_scratch/Devel/SITK/SimpleITK-build/ExternalData/Testing/Data/Input/BrainProtonDensitySliceBSplined10.png'
#gctorture(on=TRUE)
fixed = ReadImage(a1, 'sitkFloat32')

moving = ReadImage(a2, 'sitkFloat32')

commandIteration <- function(method)
{
    res <- function() {
     msg <- paste("Optimizer iteration number ", method$GetOptimizerIteration(), "\n")
     cat(msg)
    }
    return(res)
}

AT=AffineTransform(fixed$GetDimension())
initialTx = CenteredTransformInitializer(fixed, moving, AT)



Reg = ImageRegistrationMethod()

Reg$SetShrinkFactorsPerLevel(c(3,2,1))
Reg$SetSmoothingSigmasPerLevel(c(2,1,1))

Reg$SetMetricAsJointHistogramMutualInformation(20)
Reg$MetricUseFixedImageGradientFilterOff()
Reg$MetricUseMovingImageGradientFilterOff()

Reg$SetOptimizerAsGradientDescent(learningRate=1.0,
                                numberOfIterations=100,
                                convergenceMinimumValue=1e-6,
                                convergenceWindowSize=10,
                                estimateLearningRate = 'EachIteration')

Reg$SetOptimizerScalesFromPhysicalShift()

Reg$SetInitialTransform(initialTx,inPlace=TRUE)

Reg$SetInterpolator('sitkLinear')

Reg$AddRCommand('sitkIterationEvent', commandIteration(Reg))


outTx = Reg$Execute(fixed, moving)
