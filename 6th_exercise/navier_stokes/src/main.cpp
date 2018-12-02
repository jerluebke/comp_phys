#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <cstdio>
#include <windows.h>

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

/* TODO: write array to file, read initial condition from file
 *     use meta files (including target file name, resolution) */


#define VIDEO   1
#define PNG     0
#define STAT    1
#define PRINT_PROGRESS      1
#define PROGRESS_FREQUENCY  100
#define FRAMERATE   60


const size_t Nx = 2048;
const size_t Ny = 2048;
const size_t Ntot = Nx * Ny;
const size_t frames = 5400;
const unsigned int steps = 10;
const double dt = .05;
const double nu = .0;
const double xmin = -M_PI, xmax = M_PI,
             ymin = -M_PI, ymax = M_PI;


double initial_func(double x, double y);

void mean_std_ul(uint64_t *arr, size_t size, uint64_t *mean, double *std);

template<typename vtkSmartPointerT>
void set_colors(double *data, size_t npoints,
        vtkSmartPointer<vtkLookupTable> lut,
        vtkSmartPointer<vtkUnsignedCharArray> colors,
        vtkSmartPointerT source);


int main(int argc, char **argv)
{
    size_t i, j;
    double zmin, zmax;
    double *x, *y, *z;
    char filename[256];
#if STAT
    uint64_t t0_all, t1_all, t0_ts, t_m, t1_write,
             ts_max, ts_min, write_max, write_min,
             ts_mean, write_mean;
    double ts_std, write_std;
    uint64_t *ts_ts, *ts_write;
    ts_ts       = (uint64_t *)malloc(sizeof(*ts_ts) * frames);
    ts_write    = (uint64_t *)malloc(sizeof(*ts_write) * frames);
#endif

    x = (double *)malloc(sizeof(*x) * Nx);
    y = (double *)malloc(sizeof(*y) * Ny);
    z = (double *)malloc(sizeof(*z) * Ntot);


#if VIDEO || PNG
    printf("Enter filename (251 chars, without file ending): ");
    scanf_s("%s", filename, 251);
    printf("\n");
#if VIDEO
    sprintf(filename, "%s.avi", filename);
#elif PNG
    sprintf(filename, "%s.png", filename);
#endif
#endif


    //***************************************************************
    // set up values
    //***************************************************************
	linspace(xmin-1., xmax, Nx, x);
	linspace(ymin-1., ymax, Ny, y);
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
    writer->SetRate(FRAMERATE);
    writer->SetFileName(filename);
    writer->Start();


    //***************************************************************
    // mainloop
    //***************************************************************
#if STAT
    t0_all = GetTickCount64();
#endif
    set_colors<vtkSmartPointer<vtkImageCanvasSource2D>>(
            pde->o, Ntot, lut, colors, source);
    source->Update();
	writer->Write();

    for (i = 0; i < frames; ++i) {
#if STAT
        t0_ts = GetTickCount64();
#endif
        time_step(steps, pde);
#if STAT
        t_m = GetTickCount64();
#endif
        set_colors<vtkSmartPointer<vtkImageCanvasSource2D>>(
                pde->o, Ntot, lut, colors, source);
        source->Update();
        writer->Write();
#if STAT
        t1_write = GetTickCount64();
        ts_ts[i] = t_m - t0_ts;
        ts_write[i] = t1_write - t_m;
#endif
#if PRINT_PROGRESS
        if (i % PROGRESS_FREQUENCY == 0)
            printf("frame %zu\n", i);
#endif
    }

#if STAT
    t1_all = GetTickCount64();
    ts_max = *std::max_element(ts_ts, ts_ts + frames);
    ts_min = *std::min_element(ts_ts, ts_ts + frames);
    write_max = *std::max_element(ts_write, ts_write + frames);
    write_min = *std::min_element(ts_write, ts_write + frames);
    mean_std_ul(ts_ts, frames, &ts_mean, &ts_std);
    mean_std_ul(ts_write, frames, &write_mean, &write_std);

    printf("\n\nSUMMARY:\n"\
           "time step of PDE    : %llu +- %lf ns\n"     \
           "    max / min       : %llu / %llu\n"        \
           "writing per frame   : %llu +- %lf ns\n"     \
           "    max / min       : %llu / %llu\n\n"      \
           "total time          : %llu s\n\n",
            ts_mean, ts_std, ts_max, ts_min,
            write_mean, write_std, write_max, write_min,
            (t1_all - t0_all) / 1000lu);
#endif


#elif PNG

    vtkSmartPointer<vtkPNGWriter> writer =
        vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(filename);
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


    //***************************************************************
    // cleanup
    //***************************************************************
#if VIDEO
    writer->End();
#endif

    cleanup(pde);


    return EXIT_SUCCESS;
}


double initial_func(double x, double y)
{
    return  exp(-4 * (SQUARE(x-1.) + SQUARE(y-.5))) \
           +exp(-4 * (SQUARE(x+1.) + SQUARE(y+.5))) \
           -exp(-4 * (SQUARE(x-1.) + SQUARE(y+.5))) \
           -exp(-4 * (SQUARE(x+1.) + SQUARE(y-.5)));
}


void mean_std_ul(uint64_t *arr, size_t size, uint64_t *mean, double *std)
{
    uint64_t sum, sum_var, *ptr, *end;
    sum = sum_var = 0;

    ptr = arr, end = arr + size;
    while (ptr != end)
        sum += *ptr++;
    *mean = sum / size;

    ptr = arr;
    while (ptr != end) {
        sum_var += SQUARE(*ptr - *mean);
        ++ptr;
    }
    *std = std::sqrt((double) sum_var / ((double) size-1));
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
