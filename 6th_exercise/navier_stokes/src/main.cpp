#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkPlaneSource.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkImageCanvasSource2D.h>
#include <vtkImageData.h>
#include <vtkOggTheoraWriter.h>
#include <vtkAVIWriter.h>
#include <vtkPNGWriter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include "..\include\viridis.h"

extern "C" {
#include "..\include\navier_stokes.h"
}


#define VIDEO   1
#define PNG     0


const size_t Nx = 1024;
const size_t Ny = 1024;
const size_t Ntot = Nx * Ny;
const size_t frames = 3600;
const unsigned int steps = 10;
const double dt = .05;
const double nu = .0;


double initial_func(double x, double y);

template<typename vtkSmartPointerT>
void set_colors(double *data, size_t npoints,
        vtkSmartPointer<vtkLookupTable> lut,
        vtkSmartPointer<vtkUnsignedCharArray> colors,
        vtkSmartPointerT source);


int main(int argc, char **argv)
{
    size_t i, j;
    double xmin, xmax, ymin, ymax;
    double zmin, zmax;
    double *x, *y, *z;

    x = (double *)malloc(sizeof(*x) * Nx);
    y = (double *)malloc(sizeof(*y) * Ny);
    z = (double *)malloc(sizeof(*z) * Ntot);


    //***************************************************************
    // set up values
    //***************************************************************
    xmin = -2*M_PI; xmax = 2*M_PI;
    ymin = -2*M_PI; ymax = 2*M_PI;
	linspace(xmin, xmax, Nx, x);
	linspace(ymin, ymax, Ny, y);
    for (i = 0; i < Ny; ++i)
        for (j = 0; j < Nx; ++j)
            z[j+i*Nx] = initial_func(x[j], y[i]);

    zmin = *std::min_element(z, z + Ntot);
    zmax = *std::max_element(z, z + Ntot);


    //***************************************************************
    // set up PDE workspace
    //***************************************************************
    Params p = {Nx, Ny, dt, nu};
    PDE *pde = init(p, z);

    // those are not needed anymore
    free(x);
    free(y);
    free(z);


    //***************************************************************
    // set up vtk
    //***************************************************************
    // create color lookup table
    vtkSmartPointer<vtkLookupTable> lut =
        vtkSmartPointer<vtkLookupTable>::New();
    lut->SetTableRange(zmin, zmax);
    lut->SetNumberOfTableValues(256);
    for (i = 0; i < 256; ++i)
        lut->SetTableValue(i, viridis[i]);
    lut->Build();

    // vtkCharArray to store colors in [0..255] format
    vtkSmartPointer<vtkUnsignedCharArray> colors =
        vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetNumberOfTuples(Ntot);


#if VIDEO || PNG

    // image source
    vtkSmartPointer<vtkImageCanvasSource2D> source =
        vtkSmartPointer<vtkImageCanvasSource2D>::New();
    source->SetScalarTypeToUnsignedChar();
    source->SetNumberOfScalarComponents(3);
    source->SetExtent(1, Nx, 1, Ny, 0, 0);
	source->Update();
	// std::cerr << "Ntot = " << Ntot << "\n"
	//     << "npoints = " << source->GetOutput()->GetNumberOfPoints() << std::endl;
    assert(Ntot == source->GetOutput()->GetNumberOfPoints());


#if VIDEO

    // set up video writer
    vtkSmartPointer<vtkOggTheoraWriter> writer =
        vtkSmartPointer<vtkOggTheoraWriter>::New();
    // vtkSmartPointer<vtkAVIWriter> writer =
    //     vtkSmartPointer<vtkAVIWriter>::New();
    writer->SetInputConnection(source->GetOutputPort());
    writer->SetQuality(2);
    writer->SetRate(60);
    writer->SetFileName("turbolence.avi");
    writer->Start();


    //***************************************************************
    // mainloop
    //***************************************************************
    set_colors<vtkSmartPointer<vtkImageCanvasSource2D>>(
            pde->o, Ntot, lut, colors, source);
    source->Update();
	writer->Write();

    for (i = 0; i < frames; ++i) {
        time_step(steps, pde);
        set_colors<vtkSmartPointer<vtkImageCanvasSource2D>>(
                pde->o, Ntot, lut, colors, source);
        source->Update();
        writer->Write();
        if (i % 100 == 0)
            printf("frame %zu\n", i);
    }

    //***************************************************************
    // cleanup
    //***************************************************************
    writer->End();


#elif PNG

    vtkSmartPointer<vtkPNGWriter> writer =
        vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName("test.png");
    writer->SetInputConnection(source->GetOutputPort());

    time_step(steps, pde);
    set_colors<vtkSmartPointer<vtkImageCanvasSource2D>>(
            pde->o, Ntot, lut, colors, source);
    source->Update();
    writer->Write();


#endif  /* VIDEO */

#else   /* VIDEO || PNG */

    vtkSmartPointer<vtkPlaneSource> plane=
        vtkSmartPointer<vtkPlaneSource>::New();
    plane->SetResolution(Nx-1, Ny-1);
    plane->SetOrigin(xmin, ymin, 0.);
    plane->SetPoint1(xmax, ymin, 0.);
    plane->SetPoint2(xmin, ymax, 0.);
    plane->Update();
    assert(Ntot == plane->GetOutput()->GetNumberOfPoints());

    time_step(steps, pde);

    set_colors<vtkSmartPointer<vtkPlaneSource>>(
            pde->o, Ntot, lut, colors, plane);

    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(plane->GetOutput());

    vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> renderer =
        vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow =
        vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> interactor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();

    renderer->AddActor(actor);
    renderer->SetBackground(.5, .5, .5);
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("Demo");
    interactor->SetRenderWindow(renderWindow);

    renderWindow->Render();
    interactor->Start();


#endif  /* VIDEO || PNG */


    cleanup(pde);


    return EXIT_SUCCESS;
}


double initial_func(double x, double y)
{
    return exp(-4 * (SQUARE(x-1.5) + SQUARE(y))) \
           +exp(-4 * (SQUARE(x+1.5) + SQUARE(y)));
           // -exp(-4 * (SQUARE(x-1.5) + SQUARE(y+.75)))
           // +exp(-4 * (SQUARE(x+1.5) + SQUARE(y-.75)));
}


template<typename vtkSmartPointerT>
void set_colors(double *data, size_t npoints,
        vtkSmartPointer<vtkLookupTable> lut,
        vtkSmartPointer<vtkUnsignedCharArray> colors,
        vtkSmartPointerT source)
{
    size_t i, j;
    double dcolor[3];
    unsigned char ucolor[3];

    for (i = 0; i < npoints; ++i) {
		lut->GetColor(data[i], dcolor);
        for (j = 0; j < 3; ++j)
            ucolor[j] = static_cast<unsigned char>(255. * dcolor[j]);
		colors->SetTypedTuple(i, ucolor);
    }

    source->GetOutput()->GetPointData()->SetScalars(colors);
}
