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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"

#include "itkBSplineDeformableTransform.h"

#include <fstream>

//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
class CommandProgressUpdate : public itk::Command
{
public:
  typedef  CommandProgressUpdate      Self;
  typedef  itk::Command               Superclass;
  typedef itk::SmartPointer<Self>     Pointer;
  itkNewMacro( Self );
protected:
  CommandProgressUpdate() {};
public:
  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    const itk::ProcessObject * filter = dynamic_cast< const itk::ProcessObject * >( object );
    if( ! itk::ProgressEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout << filter->GetProgress() << std::endl;
    }
};


template <unsigned int VSplineOrder>
class BSplineDeformableTransformTest2Helper
{
public:
static int RunTest(int argc, char * argv [] )
{
  const     unsigned int   ImageDimension = 2;

  typedef   unsigned char                             PixelType;
  typedef   itk::Image< PixelType, ImageDimension >   FixedImageType;
  typedef   itk::Image< PixelType, ImageDimension >   MovingImageType;

  typedef   itk::ImageFileReader< FixedImageType  >   FixedReaderType;
  typedef   itk::ImageFileReader< MovingImageType >   MovingReaderType;

  typedef   itk::ImageFileWriter< MovingImageType >   MovingWriterType;


  typename FixedReaderType::Pointer fixedReader = FixedReaderType::New();
  fixedReader->SetFileName( argv[2] );

  try
    {
    fixedReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }


  typename MovingReaderType::Pointer movingReader = MovingReaderType::New();
  typename MovingWriterType::Pointer movingWriter = MovingWriterType::New();

  movingReader->SetFileName( argv[3] );
  movingWriter->SetFileName( argv[4] );


  typename FixedImageType::ConstPointer fixedImage = fixedReader->GetOutput();


  typedef itk::ResampleImageFilter< MovingImageType,
                                    FixedImageType  >  FilterType;

  typename FilterType::Pointer resampler = FilterType::New();

  typedef itk::LinearInterpolateImageFunction<
                       MovingImageType, double >  InterpolatorType;

  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

  resampler->SetInterpolator( interpolator );

  typename FixedImageType::SpacingType   fixedSpacing    = fixedImage->GetSpacing();
  typename FixedImageType::PointType     fixedOrigin     = fixedImage->GetOrigin();
  typename FixedImageType::DirectionType fixedDirection  = fixedImage->GetDirection();

  resampler->SetOutputSpacing( fixedSpacing );
  resampler->SetOutputOrigin(  fixedOrigin  );
  resampler->SetOutputDirection(  fixedDirection  );


  typename FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
  typename FixedImageType::SizeType   fixedSize =  fixedRegion.GetSize();
  resampler->SetSize( fixedSize );
  resampler->SetOutputStartIndex(  fixedRegion.GetIndex() );


  resampler->SetInput( movingReader->GetOutput() );

  movingWriter->SetInput( resampler->GetOutput() );


  const unsigned int SpaceDimension = ImageDimension;
  typedef double CoordinateRepType;

  typedef itk::BSplineDeformableTransform<
                            CoordinateRepType,
                            SpaceDimension,
                            VSplineOrder >     TransformType;

  typename TransformType::Pointer bsplineTransform = TransformType::New();


  typedef typename TransformType::RegionType RegionType;
  RegionType bsplineRegion;
  typename RegionType::SizeType   size;

  const unsigned int numberOfGridNodesOutsideTheImageSupport = VSplineOrder;

  const unsigned int numberOfGridNodesInsideTheImageSupport = 5;

  const unsigned int numberOfGridNodes =
                        numberOfGridNodesInsideTheImageSupport +
                        numberOfGridNodesOutsideTheImageSupport;

  const unsigned int numberOfGridCells =
                        numberOfGridNodesInsideTheImageSupport - 1;

  size.Fill( numberOfGridNodes );
  bsplineRegion.SetSize( size );

  typedef typename TransformType::SpacingType SpacingType;
  SpacingType spacing;

  typedef typename TransformType::OriginType OriginType;
  OriginType origin;

  spacing[0] = fixedSpacing[0] * fixedSize[0]  / numberOfGridCells;
  spacing[1] = fixedSpacing[1] * fixedSize[1]  / numberOfGridCells;

  const unsigned int orderShift = VSplineOrder / 2;

  origin[0] = fixedOrigin[0] - orderShift * spacing[0] - fixedSpacing[0] / 2.0;
  origin[1] = fixedOrigin[1] - orderShift * spacing[1] - fixedSpacing[1] / 2.0;

  bsplineTransform->SetGridSpacing( spacing );
  bsplineTransform->SetGridOrigin( origin );
  bsplineTransform->SetGridRegion( bsplineRegion );
  bsplineTransform->SetGridDirection( fixedImage->GetDirection() );


  typedef typename TransformType::ParametersType     ParametersType;

  const unsigned int numberOfParameters = bsplineTransform->GetNumberOfParameters();


  ParametersType parameters( numberOfParameters );

  std::ifstream infile;

  infile.open( argv[1] );

  const unsigned int numberOfNodes = numberOfParameters / SpaceDimension;
  for( unsigned int n=0; n < numberOfNodes; n++ )
    {
    infile >>  parameters[n];
    infile >>  parameters[n+numberOfNodes];
    }

  infile.close();

  bsplineTransform->SetParameters( parameters );

  typename CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();

  resampler->AddObserver( itk::ProgressEvent(), observer );

  resampler->SetTransform( bsplineTransform );

  try
    {
    movingWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }


  typedef itk::Point<  float, ImageDimension >        PointType;
  typedef itk::Vector< float, ImageDimension >        VectorType;
  typedef itk::Image< VectorType, ImageDimension >    DeformationFieldType;

  typename DeformationFieldType::Pointer field = DeformationFieldType::New();
  field->SetRegions( fixedRegion );
  field->SetOrigin( fixedOrigin );
  field->SetSpacing( fixedSpacing );
  field->Allocate();

  typedef itk::ImageRegionIterator< DeformationFieldType > FieldIterator;
  FieldIterator fi( field, fixedRegion );

  fi.GoToBegin();

  typename TransformType::InputPointType  fixedPoint;
  typename TransformType::OutputPointType movingPoint;
  typename DeformationFieldType::IndexType index;

  VectorType displacement;

  while( ! fi.IsAtEnd() )
    {
    index = fi.GetIndex();
    field->TransformIndexToPhysicalPoint( index, fixedPoint );
    movingPoint = bsplineTransform->TransformPoint( fixedPoint );
    displacement[0] = movingPoint[0] - fixedPoint[0];
    displacement[1] = movingPoint[1] - fixedPoint[1];
    fi.Set( displacement );
    ++fi;
    }


  typedef itk::ImageFileWriter< DeformationFieldType >  FieldWriterType;
  typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

  fieldWriter->SetInput( field );

  if( argc >= 6 )
    {
    fieldWriter->SetFileName( argv[5] );
    try
      {
      fieldWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }
    }
  return EXIT_SUCCESS;
}
};


int itkBSplineDeformableTransformTest2( int argc, char * argv[] )
{

  if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " coefficientsFile fixedImage ";
    std::cerr << "movingImage deformedMovingImage" << std::endl;
    std::cerr << "[deformationField][spline order 2,3]" << std::endl;
    return EXIT_FAILURE;
    }

   unsigned int splineOrder = 3;

   if( argc > 6 )
     {
     splineOrder = atoi( argv[6] );
     }

   int status = 0;

   switch( splineOrder )
    {
    case 1:
      {
      status |= BSplineDeformableTransformTest2Helper< 1 >::RunTest( argc, argv );
      break;
      }
    case 2:
      {
      status |= BSplineDeformableTransformTest2Helper< 2 >::RunTest( argc, argv );
      break;
      }
    case 3:
      {
      status |= BSplineDeformableTransformTest2Helper< 3 >::RunTest( argc, argv );
      break;
      }
    }

   if( status )
    {
    return EXIT_FAILURE;
    }

   return EXIT_SUCCESS;
}
