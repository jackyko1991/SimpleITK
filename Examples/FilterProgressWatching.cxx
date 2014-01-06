/*=========================================================================
*
*  Copyright Insight Software Consortium
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*         http://www.apache.org/licenses/LICENSE-2.0.txt
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*=========================================================================*/

// This one header will include all SimpleITK filters and external
// objects.
#include <SimpleITK.h>

#include <iostream>
#include <cstdlib>
#include <csignal>

// create convenient namespace alias
namespace sitk = itk::simple;

sig_atomic_t g_Abort = false;
void (*g_OldHandler)(int) = NULL;

void sitkabort_handler(int sig)
{
  if (g_Abort++)
    {
    if (g_OldHandler)
      {
      g_OldHandler(sig);
      }
    else
      {
      std::abort();
      }
    }
}

class MyFilterWatcher
  : public sitk::BasicFilterWatcher
{
public:

  typedef BasicFilterWatcher  Superclass;
  typedef MyFilterWatcher      Self;


  MyFilterWatcher(sitk::ProcessObject& o, const std::string &comment="")
    :BasicFilterWatcher(o,comment)
    {
    }

protected:

  virtual void OnAbortEvent( )
    {
      Superclass::OnAbortEvent();
    }

  virtual void OnProgressEvent( )
    {
      Superclass::OnProgressEvent();
      if (g_Abort)
        {
        this->GetProcessObject().Abort();
        }
      else
        {
        usleep(50000);
        }
    }

  virtual void OnIterationEvent( )
    {
      Superclass::OnIterationEvent();
    }
private:

};



int main ( int argc, char* argv[] ) {

  if ( argc < 4 ) {
    std::cerr << "Usage: " << argv[0] << " <input> <sigma> <output>\n";
    return 1;
  }

  g_OldHandler = std::signal(SIGINT, sitkabort_handler);

  // Read the image file
  sitk::ImageFileReader reader;
  reader.SetFileName ( std::string ( argv[1] ) );
  sitk::Image image = reader.Execute();

  sitk::Image blurredImage;
  try
    {
    // This filters perform a gaussian bluring with sigma in physical
    // space. The output image will be of real type.
    sitk::SmoothingRecursiveGaussianImageFilter gaussian;
    gaussian.SetSigma( atof ( argv[2] ) );
    gaussian.SetNumberOfThreads(1);
    gaussian.DebugOn();

    MyFilterWatcher watcher(gaussian);

    blurredImage = gaussian.Execute ( image );
    }
  catch (std::exception &e)
    {
    std::cerr << std::endl;
    std::cerr << "Exception during execution!" << std::endl;
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
    }
  catch (...)
    {
    std::cerr << std::endl;
    std::cerr << "Unknown Exception during execution!" << std::endl;
    throw;
    }

  // Covert the real output image back to the original pixel type, to
  // make writing easier, as many file formats don't support real
  // pixels.
  sitk::CastImageFilter caster;
  caster.SetOutputPixelType( image.GetPixelIDValue() );
  sitk::Image outputImage = caster.Execute( blurredImage );

  // write the image
  sitk::ImageFileWriter writer;
  writer.SetFileName ( std::string ( argv[3] ) );
  writer.Execute ( outputImage );

  return 0;
}
