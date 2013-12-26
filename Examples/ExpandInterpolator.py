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

if len (sys.argv) < 3:
    print("Usage: %s <input> <output>" %sys.argv[0])
    sys.exit (1)


image = sitk.ReadImage(sys.argv[1])

size = image.GetSize();
image = image[(size[0]//2-25):(size[0]//2+25),(size[1]//2-25):(size[1]//2+25)]


iterps = [sitk.sitkNearestNeighbor,
          sitk.sitkLinear,
          sitk.sitkBSpline,
          sitk.sitkGaussian,
          sitk.sitkHammingWindowedSinc,
          sitk.sitkCosineWindowedSinc,
          sitk.sitkWelchWindowedSinc,
          sitk.sitkLanczosWindowedSinc,
          sitk.sitkBlackmanWindowedSinc]

eFactor=5

image_list = []

for i in iterps:
    image_list.append( sitk.Expand( image, [eFactor]*3, i ))

tiles = sitk.Tile( image_list, [3,0] )


sitk.WriteImage(tiles, sys.argv[2] )

if ( not "SITK_NOSHOW" in os.environ ):
    sitk.Show( tiles, "Expand" )


image_list = []

resample = sitk.ResampleImageFilter()
resample.SetReferenceImage(image)
resample.SetSize([image.GetSize()[i]*eFactor for i in range(len(size))])
resample.SetOutputSpacing([image.GetSpacing()[i]/eFactor for i in range(len(size))])
resample.SetOutputOrigin([ image.GetOrigin()[i]-image.GetSpacing()[i]*(0.5-1.0/(2.0*eFactor)) for i in range(len(size))])
resample.SetOutputDirection(image.GetDirection())

for i in iterps:
    resample.SetInterpolator( i )
    image_list.append( resample.Execute( image ))

tiles = sitk.Tile( image_list, [3,0] )

sitk.Show( tiles, "Resample")
