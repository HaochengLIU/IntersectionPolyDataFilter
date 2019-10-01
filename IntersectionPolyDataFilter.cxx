#include <vtkActor.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkCellCenters.h>
#include <vtkGlyph3D.h>
#include <vtkGlyph3DMapper.h>
#include <vtkType.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkLabeledDataMapper.h>
#include <vtkNamedColors.h>
#include <vtkTextProperty.h>
#include <vtkActor2D.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include "vtkRegressionTestImage.h"

void tripoint(double P1[3], double P2[3], double P3[3], double u, double v, double P[3])
{
	double w = 1 - u - v;
	for (int i = 0; i<3; i++)
		P[i] = P1[i] * w + P2[i] * u + P3[i] * v;
}

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////////////////////
  // Build two spheres and make the intersection
  /////////////////////////////////////////////////////////////

  vtkSmartPointer<vtkSphereSource> sphereSource1 = vtkSmartPointer<vtkSphereSource>::New();
  sphereSource1->SetCenter(0.0, 0.0, 0.0);
  sphereSource1->SetRadius(2.0f);
  sphereSource1->SetPhiResolution(6);
  sphereSource1->SetThetaResolution(8);
  sphereSource1->Update();
  vtkSmartPointer<vtkPolyDataMapper> sphere1Mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  //sphere1Mapper->SetInputConnection( sphereSource1->GetOutputPort() );
  //sphere1Mapper->ScalarVisibilityOff();
  //vtkSmartPointer<vtkActor> sphere1Actor = vtkSmartPointer<vtkActor>::New();
  //sphere1Actor->SetMapper( sphere1Mapper );
  ////sphere1Actor->GetProperty()->SetOpacity(.3);
  //sphere1Actor->GetProperty()->SetColor(1,0,0);
  //sphere1Actor->GetProperty()->SetEdgeVisibility(1);

  vtkSmartPointer<vtkSphereSource> sphereSource2 = vtkSmartPointer<vtkSphereSource>::New();
  sphereSource2->SetCenter(1.0, 0.0, 0.0);
  sphereSource2->SetRadius(2.0f);
  sphereSource2->SetPhiResolution(6);
  sphereSource2->SetThetaResolution(8);
  vtkSmartPointer<vtkPolyDataMapper> sphere2Mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  //sphere2Mapper->SetInputConnection( sphereSource2->GetOutputPort() );
  //sphere2Mapper->ScalarVisibilityOff();
  //vtkSmartPointer<vtkActor> sphere2Actor = vtkSmartPointer<vtkActor>::New();
  //sphere2Actor->SetMapper( sphere2Mapper );
  ////sphere2Actor->GetProperty()->SetOpacity(.3);
  //sphere2Actor->GetProperty()->SetColor(0,1,0);
  //sphere2Actor->GetProperty()->SetEdgeVisibility(1);

  vtkSmartPointer<vtkIntersectionPolyDataFilter> intersectionPolyDataFilter =
    vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
  intersectionPolyDataFilter->SetInputConnection( 0, sphereSource1->GetOutputPort() );
  intersectionPolyDataFilter->SetInputConnection( 1, sphereSource2->GetOutputPort() );
  
  ///////////////////////////////////////////////////////
  // Get output from first port
  ///////////////////////////////////////////////////////

  /*******************************************/
  // You must call update first before trying to fetch the data from vtk filter. It's because VTK uses
  // lazy evaluation and downstream users needs to ask data from the upstream providers.
  intersectionPolyDataFilter->Update();

  vtkSmartPointer<vtkPolyDataMapper> intersectionMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  intersectionMapper->SetInputConnection( intersectionPolyDataFilter->GetOutputPort() );
  intersectionMapper->ScalarVisibilityOff();

  vtkSmartPointer<vtkActor> intersectionActor = vtkSmartPointer<vtkActor>::New();
  intersectionActor->SetMapper( intersectionMapper );
  intersectionActor->GetProperty()->SetColor(1, 1, 1);
  intersectionActor->GetProperty()->SetLineWidth(2);

  // Output data of the first port
  vtkSmartPointer<vtkPolyData> output = intersectionPolyDataFilter->GetOutput(0);
  std::cout << "********** Output data of first port ******************" << std::endl;
  std::cout << "Number of outputPoints=" << output->GetNumberOfPoints() << std::endl;
  std::cout << "Number of outputLines=" << output->GetNumberOfLines() << std::endl;
  std::cout << "Number of outputCells=" << output->GetNumberOfCells() << std::endl;
  std::cout << "Number of arrays in point data=" << output->GetPointData()->GetNumberOfArrays() << std::endl;
  auto pd = output->GetPointData();
  for (int i = 0; i < pd->GetNumberOfArrays(); i++)
  {
	  auto array = pd->GetArray(i);
	  std::cout << "  Print data in point array " << array->GetName() << std::endl;
	  // Let's say we have an array stores data as [[1,2],[2,3],[3,4],[1,2]]. It has 4 tuples and numberOfComponet
	  // is 2.
	  for (int j = 0; j < pd->GetNumberOfTuples(); j++)
	  {
		  for (int k = 0; k < pd->GetNumberOfComponents(); k++)
		  {
			  std::cout << array->GetComponent(j, k) << " ";
		  }
	  }
	  std::cout << std::endl;
  }
  std::cout << "Number of arrays in cell data=" << output->GetCellData()->GetNumberOfArrays() << std::endl;
  auto cd = output->GetCellData();
  for (int i = 0; i < cd->GetNumberOfArrays(); i++)
  {
	  auto array = cd->GetArray(i);
	  std::cout << "  Print data in cell array " << array->GetName() << std::endl;
	  std::cout << "  Number Of Components " << array->GetNumberOfComponents() << std::endl;
	  // You can get the array or print out the underlying data using the same logic above.
	  for (int j = 0; j < array->GetNumberOfTuples(); j++)
	  {
		  for (int k = 0; k < array->GetNumberOfComponents(); k++)
		  {
			  std::cout << array->GetComponent(j, k) << " ";
		  }
	  }
	  std::cout << std::endl;
  }

  // Same logic to get output data from second and third port
  /*******************************************/

  ///////////////////////////////////////////////////////
  // Get output data from second port
  ///////////////////////////////////////////////////////

  vtkSmartPointer<vtkPolyData> output1 = intersectionPolyDataFilter->GetOutput(1); // For the third port change it to 2
  output1->BuildCells();
  vtkSmartPointer<vtkPoints> points1 = output1->GetPoints();
  std::cout << "********** Output data of second port ******************" << std::endl;
  std::cout << "Number of outputPoints=" << output1->GetNumberOfPoints() << std::endl;
  std::cout << "Number of outputLines=" << output1->GetNumberOfLines() << std::endl;
  std::cout << "Number of outputCells=" << output1->GetNumberOfCells() << std::endl;
  std::cout << "Number of points=" << 
	  points1->GetNumberOfPoints() << std::endl;
  double pt1[3];
  for (int i = 0; i < points1->GetNumberOfPoints(); i++)
  {
	  points1->GetPoint(i, pt1);
	  for (int j = 0; j < 3; j++)
		  std::cout << pt1[j] << " ";
	  std::cout << std::endl;
  }

  // The triangles and their interpolated centers
  std::cout << "There are " << output1->GetNumberOfPolys() << " polys." << std::endl;
  vtkIdType *pts1;
  for (int i = 0; i<output1->GetNumberOfPolys(); i++)
  {
	  output1->GetCell(i, pts1);
	  std::cout << "Poly " << i << " is:  ";
	  for (int j = 0; j < pts1[0]; j++)
		  std::cout << pts1[j + 1] << " ";
	  std::cout << std::endl;
  }

  // The intersected triangles
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
  auto labelScalars = vtkSmartPointer<vtkDoubleArray>::New();
  labelScalars->SetNumberOfComponents(1);
  vtkIdType pid[1];
  double p1[3], p2[3], p3[3], p[3], u = 1.0 / 3.0, v = u;
  std::cout << std::endl << "******* The intersected triangles ***********" << std::endl;
  auto array = cd->GetArray(2); // For the third port change it to 3
  std::cout << "  Print data in cell array " << array->GetName() << std::endl;
  std::cout << "  Number Of Components " << array->GetNumberOfComponents() << std::endl;
  points->SetNumberOfPoints((array->GetNumberOfTuples())*(array->GetNumberOfComponents()));
  int i = 0;
  for (int j = 0; j < array->GetNumberOfTuples(); j++)
  {
	  for (int k = 0; k < array->GetNumberOfComponents(); k++)
	  {
		  pid[0] = i;
		  vertices->InsertNextCell(1, pid);
		  labelScalars->InsertNextTuple1(array->GetComponent(j, k));
		  std::cout << array->GetComponent(j, k) << " ";
		  output1->GetCell(array->GetComponent(j, k), pts1);
		  points1->GetPoint(pts1[1], p1);
		  points1->GetPoint(pts1[2], p2);
		  points1->GetPoint(pts1[3], p3);
		  tripoint(p1, p2, p3, u, v, p);
		  points->SetPoint(i, p);
		  i++;
	  }
  }
  std::cout << std::endl;
  std::cout << std::endl << "******* points ***********" << std::endl;
  points->Print(std::cout);
  std::cout << std::endl << "******* vertices ***********" << std::endl;
  vertices->Print(std::cout);

  //////////////////////////////////////////////
  // The first split surface
  //////////////////////////////////////////////
  sphere1Mapper->SetInputConnection(intersectionPolyDataFilter->GetOutputPort(1));
  sphere1Mapper->ScalarVisibilityOff();
  vtkSmartPointer<vtkActor> sphere1Actor =
	  vtkSmartPointer<vtkActor>::New();
  sphere1Actor->SetMapper(sphere1Mapper);
  //sphere1Actor->GetProperty()->SetOpacity(.3);
  sphere1Actor->GetProperty()->SetColor(0.0, 0.0, 1);
  sphere1Actor->GetProperty()->SetEdgeVisibility(1);
  sphere1Actor->GetProperty()->SetInterpolationToFlat();

  //////////////////////////////////////////////
  // The second split surface
  //////////////////////////////////////////////

  sphere2Mapper->SetInputConnection(intersectionPolyDataFilter->GetOutputPort(2));
  sphere2Mapper->ScalarVisibilityOff();
  vtkSmartPointer<vtkActor> sphere2Actor =
	  vtkSmartPointer<vtkActor>::New();
  sphere2Actor->SetMapper(sphere2Mapper);
  //sphere2Actor->GetProperty()->SetOpacity(.3);
  sphere2Actor->GetProperty()->SetColor(1, 0.0, 0.0);
  sphere2Actor->GetProperty()->SetEdgeVisibility(1);
  sphere2Actor->GetProperty()->SetInterpolationToFlat();

  /////////////////////////////////////////////////////
  //  Create a point set
  /////////////////////////////////////////////////////

  vtkSmartPointer<vtkPolyData> intersect =
	  intersectionPolyDataFilter->GetOutput(0);
  vtkSmartPointer<vtkPolyData> outSphere1 =
	  intersectionPolyDataFilter->GetOutput(1);
  vtkSmartPointer<vtkPolyData> outSphere2 =
	  intersectionPolyDataFilter->GetOutput(2);
  vtkSmartPointer<vtkPolyData> pointsData =
	  vtkSmartPointer<vtkPolyData>::New();
  //pointsData->SetPoints(outSphere1->GetPoints());
  //pointsData->SetVerts(outSphere1->GetVerts());
  pointsData->SetPoints(points);
  pointsData->SetVerts(vertices);
  pointsData->GetPointData()->SetScalars(labelScalars);

  vtkSmartPointer<vtkSphereSource> sphere =
	  vtkSmartPointer<vtkSphereSource>::New();
  sphere->SetPhiResolution(11);
  sphere->SetThetaResolution(11);
  sphere->SetRadius(0.03);

  // Create a mapper and actor
  vtkSmartPointer<vtkCellCenters> centers =
	  vtkSmartPointer<vtkCellCenters>::New();
  centers->SetInputData(pointsData);
  vtkSmartPointer<vtkGlyph3DMapper> pointMapper =
	  vtkSmartPointer<vtkGlyph3DMapper>::New();
  pointMapper->SetInputConnection(centers->GetOutputPort());
  pointMapper->SetSourceConnection(sphere->GetOutputPort());

  vtkSmartPointer<vtkLabeledDataMapper> labelMapper =
	  vtkSmartPointer<vtkLabeledDataMapper>::New();
  labelMapper->SetInputData(pointsData);
  labelMapper->SetLabelModeToLabelScalars();
  labelMapper->SetLabelFormat("%6.0f");

  vtkSmartPointer<vtkActor2D> labelActor =
	  vtkSmartPointer<vtkActor2D>::New();
  labelActor->SetMapper(labelMapper);

  vtkSmartPointer<vtkActor> pointActor =
	  vtkSmartPointer<vtkActor>::New();
  pointActor->SetMapper(pointMapper);
  pointActor->GetProperty()->SetPointSize(8);
  pointActor->GetProperty()->SetColor(0, 1, 0);

  /////////////////////////////////////////////////////
  //  Rendering
  /////////////////////////////////////////////////////
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddViewProp(sphere1Actor);
  //renderer->AddViewProp(sphere2Actor);
  renderer->AddViewProp(intersectionActor);
  renderer->AddViewProp(labelActor);
  renderer->AddViewProp(pointActor);
  renderer->SetBackground(0, 0, 0);

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer( renderer );

  vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renWinInteractor->SetRenderWindow( renderWindow );

  std::cout << "******* intersectionPolyDataFilter ***********" << std::endl;
  intersectionPolyDataFilter->Print(std::cout);
  //std::cout << std::endl << "******* intersect ***********" << std::endl;
  //intersect->Print(std::cout);

  renderWindow->Render();
  //renWinInteractor->Start();

  int retVal = vtkRegressionTestImage(renderWindow);
  retVal = 3;
  if (retVal == vtkRegressionTester::DO_INTERACTOR)
  {
	  renWinInteractor->Start();
  }

  return !retVal;
}
