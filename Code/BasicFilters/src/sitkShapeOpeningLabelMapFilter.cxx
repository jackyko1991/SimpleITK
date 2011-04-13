#include "sitkShapeOpeningLabelMapFilter.h"

#include "itkShapeLabelMapFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkConvertLabelMapFilter.h"

#include "itkLabelMap.h"
#include "itkLabelObject.h"
#include "itkShapeLabelObject.h"

namespace itk {
namespace simple {

ShapeOpeningLabelMapFilter::ShapeOpeningLabelMapFilter()
{

  typedef LabelPixelIDTypeList PixelIDTypeList;
  this->m_MemberFactory.reset( new detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();
}

// Print ourselves out
std::string ShapeOpeningLabelMapFilter::ToString() const
{
  std::ostringstream out;
  out << "itk::simple::ShapeOpeningLabelMapFilter\n";
  return out.str();
}

Image ShapeOpeningLabelMapFilter::Execute ( const Image& image )
{
  PixelIDValueType type = image.GetPixelIDValue();
  unsigned int dimension = image.GetDimension();

  return this->m_MemberFactory->GetMemberFunction( type, dimension )( image );
}


template <class TImageType> Image ShapeOpeningLabelMapFilter::ExecuteInternal ( const Image& inImage )
{
  typedef TImageType InputImageType;
  typedef itk::ShapeLabelObject<typename InputImageType::LabelType, InputImageType::ImageDimension > ShapeLabelObjectType;
  typedef itk::LabelMap< ShapeLabelObjectType  > ShapeLabelMapType;

  typename InputImageType::ConstPointer image =
    dynamic_cast <const InputImageType*> ( inImage.GetImageBase() );

  if ( image.IsNull() )
    {
    sitkExceptionMacro( << "Could not cast input image to proper type" );
    }

  typedef itk::ConvertLabelMapFilter<  InputImageType, ShapeLabelMapType  > ConvertFilterType;
  typename ConvertFilterType::Pointer convert = ConvertFilterType::New();
  convert->SetInput( image );
  convert->SetNumberOfThreads(1);

  typedef itk::ShapeLabelMapFilter<ShapeLabelMapType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(convert->GetOutput() );
  filter->InPlaceOn();
  filter->Update();

  typename ShapeLabelMapType::Pointer labelMap = filter->GetOutput();
  for ( size_t label = 1; label < labelMap->GetNumberOfLabelObjects(); ++label )
    {
    const ShapeLabelObjectType *lo = labelMap->GetLabelObject( label );
    std::cout << "label: " << label << "\t" << lo->GetPhysicalSize() << "\t" << lo->GetRoundness() << std::endl;
    }

  typedef itk::ShapeOpeningLabelMapFilter< ShapeLabelMapType> OpeningFilterType;
  typename OpeningFilterType::Pointer opening = OpeningFilterType::New();
  opening->SetInput( filter->GetOutput() );
  opening->InPlaceOn();

  typedef itk::ConvertLabelMapFilter< ShapeLabelMapType, InputImageType > ConvertBackFilterType;
  typename ConvertBackFilterType::Pointer convertBack = ConvertBackFilterType::New();
  convertBack->SetInput( opening->GetOutput() );

  return Image( convertBack->GetOutput() );
}


Image ShapeOpeningLabelMap( const Image& image )
{
  return ShapeOpeningLabelMapFilter().Execute( image );
}
}
}

