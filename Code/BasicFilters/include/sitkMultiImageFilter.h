#ifndef __sitkMultiImageFilter_h
#define __sitkMultiImageFilter_h

#include "sitkImageFilter.h"

namespace itk {
  namespace simple {

  /** \class MultiImageFilter
   * \brief The base interface for SimpleITK filters that take two input images
   *
   * All SimpleITK filters which take two input images should inherit from this
   * class
   */
  class MultiImageFilter :
      protected NonCopyable
  {
    public:

      /**
       * Default Constructor that takes no arguments and initializes
       * default parameters
       */
      MultiImageFilter();


      // Print ourselves out
      virtual std::string ToString() const = 0;

      virtual Image Execute ( const std::vector<Image> &) = 0;

     virtual ~MultiImageFilter();

    private:

    };

  } // end namespace simple
}// end namespace itk
#endif //__sitkMultiImageFilter_h
