#ifndef __sitkImageTypeTraits_h
#define __sitkImageTypeTraits_h

namespace itk
{
namespace simple
{

template <typename TImageType>
struct ImageTypeTraits;

/**
* example usage:
* typedef typename ImageTypeTraits<TImageType>::template RebindPixelType< float >::Type OutputImageType;
*
*/
template <typename TPixelType, unsigned int TImageDimension>
struct ImageTypeTraits< itk::Image<TPixelType, TImageDimension> >
{
  typedef itk::Image<TPixelType, TImageDimension> Type;
  typedef TPixelType                              PixelType;
  static const unsigned int ImageDimension = TImageDimension;

  template <typename UPixelType, unsigned int UImageDimension = ImageDimension>
  struct Rebind
    {
      typedef itk::Image<UPixelType, UImageDimension>  Type;
    };

  typedef itk::Image<typename itk::NumericTraits<PixelType>::AccumulateType, ImageDimension > AccumulateType;
  typedef itk::Image<typename itk::NumericTraits<PixelType>::AbsType, ImageDimension >        AbsType;
  typedef itk::Image<typename itk::NumericTraits<PixelType>::FloatType, ImageDimension >      FloatType;
  typedef itk::Image<typename itk::NumericTraits<PixelType>::RealType, ImageDimension >       RealType;

};

template <typename TPixelType, unsigned int TImageDimension>
struct ImageTypeTraits< itk::VectorImage<TPixelType, TImageDimension> >
{
  typedef itk::VectorImage<TPixelType, TImageDimension> Type;
  typedef TPixelType                                    PixelType;
  static const unsigned int ImageDimension = TImageDimension;

  template <typename UPixelType, unsigned int UImageDimension = TImageDimension>
  struct Rebind
  {
    typedef itk::VectorImage<UPixelType, UImageDimension>  Type;
  };
};


template <typename TPixelType, unsigned int TImageDimension,
          template <typename SPixelType, unsigned int SImageDimension>  class TLabelObject >
struct ImageTypeTraits< itk::LabelMap< TLabelObject<TPixelType, TImageDimension> > >
{
  typedef itk::LabelMap< TLabelObject<TPixelType, TImageDimension> > Type;
  typedef TPixelType                                                 PixelType;
  static const unsigned int ImageDimension = TImageDimension;

  template <typename UPixelType, unsigned int UImageDimension = TImageDimension>
  struct Rebind
  {
    typedef itk::LabelMap< TLabelObject<UPixelType, UImageDimension> > Type;
  };
};


} // end namespace simple
} // end namespace itk

#endif // __sitkImageTypeTraits_h
