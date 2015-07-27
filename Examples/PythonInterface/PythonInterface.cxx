#include "Python.h"
#include "SimpleITK.h"
#include "itkImage.h"


namespace
{

template<unsigned int Dimension>
void DoIt(itk::simple::Image & img)
{
  typedef itk::Image<float, Dimension> ImageType;
  typename ImageType::Pointer inputImage = dynamic_cast<ImageType*>( img.GetITKBase() );

  std::cout << "ITK Image: " << (void *) inputImage.GetPointer() << std::endl;

  if ( inputImage.IsNull() )
    {
    std::cerr << "Failure to convert to ITK Image!";
    }

  std::cout << inputImage;

}
}

extern "C" {

extern itk::simple::Image * convertSwigSimpleITKImage( PyObject * pyObj );



static PyObject * method( PyObject *self, PyObject *obj)
{
  namespace sitk = itk::simple;

  const char thisStr[] = "this";

  //first we need to get the this attribute from the Python Object
  if (!PyObject_HasAttrString(obj, thisStr))
    return NULL;

  PyObject* thisAttr = PyObject_GetAttrString(obj, thisStr);
  if (thisAttr == NULL)
    return NULL;

  const sitk::Image *pyImage = convertSwigSimpleITKImage( thisAttr );

  if (pyImage == NULL)
    {
    Py_DECREF(thisAttr);

    std::cerr << "Conversion Failure!" << std::endl;
    // failure to convert
    Py_INCREF(Py_None);
    return Py_None;
    }

  // perform shallow copy, to do sitk memory management and not python
  sitk::Image img(*pyImage);


  if (img.GetPixelID() != sitk::sitkFloat32)
    {
    img = sitk::Cast(img, sitk::sitkFloat32);
    }

  switch (img.GetDimension())
    {
    case 2:
      DoIt<2>(img);
      break;

    case 3:
      DoIt<3>(img);
      break;

    case 4:
      DoIt<4>(img);
      break;

    default:
      std::cerr << "Unexpected image dimension!" << std::endl;
    }

  Py_DECREF(thisAttr);

  Py_INCREF(Py_None);
  return Py_None;
}




static PyMethodDef module_methods[] = {
   { "method", (PyCFunction)method, METH_O, NULL },
   { NULL }
};


void initSimpleITKPythonInterface(void)
{
    Py_InitModule("SimpleITKPythonInterface", module_methods );
}

}
