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


fixedInput = sitk.ReadImage(sys.argv[1])
if fixedInput.GetNumberOfComponentsPerPixel() > 1:
    fixed = sitk.VectorIndexSelectionCast(fixedInput,0,sitk.sitkFloat32)
else:
    fixed = sitk.Cast(fixedInput,sitk.sitkFloat32)

movingInput = sitk.ReadImage(sys.argv[2])
if movingInput.GetNumberOfComponentsPerPixel() > 1:
    moving = sitk.VectorIndexSelectionCast(movingInput,0,sitk.sitkFloat32)
else:
    moving = sitk.Cast(movingInput,sitk.sitkFloat32)


R = sitk.ImageRegistrationMethod()
R.SetMetricAsMeanSquares()
R.SetOptimizerAsRegularStepGradientDescent(1, .0001, 400, relaxationFactor=0.6 )
R.SetOptimizerScales([5.0, 5.0, 5.0, .01, .01, .01])
R.SetTransform(sitk.Transform(fixed.GetDimension(), sitk.sitkVersorRigid))
R.SetInterpolator(sitk.sitkLinear)

outTx = R.Execute(fixed, moving)

print("-------")
print(outTx)

resampler = sitk.ResampleImageFilter()
resampler.SetReferenceImage(fixed);
resampler.SetInterpolator(sitk.sitkLinear)
resampler.SetDefaultPixelValue(100)
resampler.SetTransform(outTx)

outImg = resampler.Execute(movingInput)

#if ( not "SITK_NOSHOW" in os.environ ):
sitk.Show( sitk.Compose(sitk.Cast(outImg, sitk.sitkUInt8),
                        sitk.Cast(fixed, sitk.sitkUInt8),
                        sitk.Cast(outImg, sitk.sitkUInt8) ) , "CheckerBoard" )
sitk.WriteImage(outImg, sys.argv[3])
