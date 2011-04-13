
#include "sitkImageToLabelMapFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"


namespace itk {
  namespace simple {
    ImageToLabelMapFilter::ImageToLabelMapFilter () {
      this->m_MemberFactory.reset( new detail::MemberFunctionFactory<MemberFunctionType>( this ) );

      this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
      this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();
    }

    std::string ImageToLabelMapFilter::ToString() const {
      std::ostringstream out;
      out << "itk::simple::ImageToLabelMapFilter\n";
      return out.str();
    }

    Image ImageToLabelMapFilter::Execute ( const Image& image ) {

      PixelIDValueType type = image.GetPixelIDValue();
      unsigned int dimension = image.GetDimension();

      return this->m_MemberFactory->GetMemberFunction( type, dimension )( image );
    }

    template <class TImageType>
    Image ImageToLabelMapFilter::ExecuteInternal ( const Image& inImage )
    {
      typedef TImageType InputImageType;

      typedef itk::LabelObject< uint32_t, InputImageType::ImageDimension  > LabelObjectType;

      typedef itk::LabelMap< LabelObjectType > OutputImageType;

      typename InputImageType::ConstPointer image =
        dynamic_cast <const InputImageType*> ( inImage.GetImageBase() );

      typedef itk::BinaryImageToLabelMapFilter<InputImageType, OutputImageType> FilterType;
      typename FilterType::Pointer filter = FilterType::New();
      filter->SetInput( image );

      filter->Update();

      return Image( filter->GetOutput() );
    }

    Image ImageToLabelMap ( const Image& image ) {
      return ImageToLabelMapFilter().Execute ( image );
    }
  }
}
