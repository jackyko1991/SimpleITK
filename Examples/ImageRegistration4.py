#!/usr/bin/env python
#=========================================================================
#
#  Copyright Insight Software Consortium
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0.txt
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#=========================================================================

from __future__ import print_function

import SimpleITK as sitk
import sys
import os


if len ( sys.argv ) < 4:
    print( "Usage: {0} <fixedImageFilter> <movingImageFile> <outputImageFile> <numberOfBins> <numberOfSamples> <useExplicitPDFderivatives>".format(sys.argv[0]))
    sys.exit ( 1 )


fixedInput = sitk.ReadImage(sys.argv[1])
fixed = sitk.VectorIndexSelectionCast(fixedInput,0,sitk.sitkFloat32)


movingInput = sitk.ReadImage(sys.argv[2])
moving = sitk.VectorIndexSelectionCast(movingInput,0,sitk.sitkFloat32)

numberOfBins = 24
numberOfSamples = 10000
useExplicitPDFDerivatives = True

if len ( sys.argv ) > 4:
    numberOfBins = int(sys.argv[4])
if len ( sys.argv ) > 5:
    numberOfSamples = int(sys.argv[5])
if len ( sys.argv ) > 6:
    useExplicitPDFDerivatives = bool(sys.argv[6])

print( numberOfBins, numberOfSamples, useExplicitPDFDerivatives)

R = sitk.ImageRegistrationMethod()
R.SetMetricAsMattesMutualInformation(numberOfBins,useExplicitPDFDerivatives,numberOfSpatialSamples=numberOfSamples,)
R.SetOptimizerAsRegularStepGradientDescent(1.0,.001,200)
R.SetTransform(sitk.Transform(fixed.GetDimension(), sitk.sitkTranslation))
R.SetInterpolator(sitk.sitkLinear)
R.DebugOn()
outTx = R.Execute(fixed, moving)


resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixed);
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(100)
resampler.SetTransform(outTx)

outImg = resampler.Execute(movingInput)

#if ( not "SITK_NOSHOW" in os.environ ):

sitk.Show(sitk.Compose(sitk.VectorIndexSelectionCast(outImg, 0, sitk.sitkUInt8),
                       sitk.Cast(fixed, sitk.sitkUInt8),
                       sitk.VectorIndexSelectionCast(outImg, 0, sitk.sitkUInt8)
                       ), "ImageRegistration4" )


sitk.WriteImage(outImg, sys.argv[3])
