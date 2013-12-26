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
    print( "Usage: {0} <fixedImageFilter> <movingImageFile> <outputImageFile>".format(sys.argv[0]))
    sys.exit ( 1 )

pixelType = sitk.sitkUInt8

fixedInput = sitk.ReadImage(sys.argv[1])
if fixedInput.GetNumberOfComponentsPerPixel() > 1:
    fixed = sitk.VectorIndexSelectionCast(fixedInput,0,pixelType)
else:
    fixed = sitk.Cast(fixedInput,pixelType)

movingInput = sitk.ReadImage(sys.argv[2])
if movingInput.GetNumberOfComponentsPerPixel() > 1:
    moving = sitk.VectorIndexSelectionCast(movingInput,0,pixelType)
else:
    moving = sitk.Cast(movingInput,pixelType)

tx = sitk.Transform(fixed.GetDimension(), sitk.sitkTranslation)
tx.SetParameters( [-11,-13] )

R = sitk.ImageRegistrationMethod()
R.SetMetricAsMatchCardinality(False, 5000)
R.SetOptimizerAsAmoeba(5.0,0.25,.001,500)
#R.SetOptimizerScales([.2,.2])
R.SetTransform(tx)
R.SetInterpolator(sitk.sitkNearestNeighbor)
R.DebugOn()
outTx = R.Execute(fixed, moving)

print("-------")
print(outTx)
print("Optimizer stop condition: {0}".format(R.GetOptimizerStopConditionDescription()))
print(" Iteration: {0}".format(R.GetOptimizerIteration()))
print(" Metric value: {0}".format(R.GetMetricValue()))

resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixed);
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(100)
resampler.SetTransform(outTx)

outImg = resampler.Execute(moving)

#if ( not "SITK_NOSHOW" in os.environ ):
#sitk.Show( sitk.CheckerBoard( outImg, fixed, [10,10] ) , "CheckerBoard" )
#sitk.WriteImage(outImg, sys.argv[3])
