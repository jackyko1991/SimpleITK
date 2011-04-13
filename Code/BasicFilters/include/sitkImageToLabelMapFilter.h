#ifndef __sitkImageToLabelMapFilter_h
#define __sitkImageToLabelMapFilter_h

#include "sitkMacro.h"
#include "sitkMemberFunctionFactory.h"
#include "sitkImage.h"

namespace itk {
  namespace simple {

    /** \class ImageToLabelMapFilter
     * \brief Convert an binary image into multiple labels
     */
    class ImageToLabelMapFilter {
    public:
      typedef ImageToLabelMapFilter Self;

      // function pointer type
      typedef Image (Self::*MemberFunctionType)( const Image& );

      // this filter works with all itk::Image and itk::VectorImage
      // types.
      typedef typelist::MakeTypeList<BasicPixelID<uint8_t>,
                                     BasicPixelID<uint16_t>,
                                     BasicPixelID<uint32_t> >::Type
            PixelIDTypeList;

      ImageToLabelMapFilter();

      // Print ourselves out
      std::string ToString() const;

      Image Execute ( const Image& );


    private:

      template <class TImageType> Image ExecuteInternal ( const Image&  image );

      // friend to get access to executeInternal member
      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::auto_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;
    };

    Image ImageToLabelMap( const Image& image );
  }
}
#endif
