#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

void writeVtiCells(float* u, std::string destination, int lx, int hx, int ly, int hy)
{
    int ox = 0, oy = 0;
    int sx = 1, sy = 1;
    int size = (hx-lx)*(hy-ly);
    int offset = 0;
    int binSize = size*4 + 4;

    // std::cout << "size = " << size << std::endl;

    // header
    {
	std::ofstream os(destination.c_str());

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	os << "<ImageData WholeExtent=\"" << lx << " " << hx << " " << ly<< " "
	    << hy << " 0 0\" Origin=\"" << ox << " " << oy
	    << " 0\" Spacing=\"" << sx << " " << sy << " 1\">\n";
	os << "<Piece Extent=\"" << lx << " " << hx << " " << ly << " " << hy << " 0 0\">\n";

	os << "<CellData Scalars=\""<< "output" << "\">\n";

	os << "<DataArray type=\"Float32\" Name=\"u\" format=\"appended\" offset=\"" << offset << "\"/>\n";
	//offset += binSize;

	os << "</CellData>\n";
	os << "</Piece>\n";
	os << "</ImageData>\n";
	os << "<AppendedData encoding=\"raw\">\n";
	os << "_" ;

	os.close();
    }

    // binary data
    {
	int N = 4*size;
	std::ofstream os(destination.c_str(), std::ios::binary | std::ios::app);
	// write the number of bytes
	os.write(reinterpret_cast<const char*>(&N),sizeof(int));
	for(int i = 0; i < size; ++i) os.write(reinterpret_cast<const char*>(&u[i]), sizeof(float));

	os.close();
    }

    // footer
    {
	std::ofstream os(destination.c_str(), std::ios::app);

	os << "\n</AppendedData>\n";
	os << "</VTKFile>\n";

	os.close();
    }
}

void writePVtiCells(std::string destination, int time,
	int HX, int HY, int nX, int nY)
{
    std::stringstream sFile;
    sFile << destination << "/output" << "_" << time << ".pvti";
    std::ofstream os(sFile.str().c_str());
    os << "<?xml version=\"1.0\" ?>\n";
    os << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "<PImageData WholeExtent=\"" << 0 << " " << HX << " " << 0 << " "
	<< HY << " 0 0\" GhostLevel=\"2\" Origin=\"" << 0 << " " <<  0 << " " <<  0
	<< "\" Spacing=\"" << 1 << " "<<  1 << " "<<  1 << "\">\n";
    os << "<PCellData Scalars=\"" << "output" << "\">\n";

    os << "<PDataArray type=\"Float32\" Name=\"u\"/>\n";

    os << "</PCellData>\n";
    for(int hx = nX; hx < HX+1; hx += nX)
    {
	for(int hy = nY; hy < HY+1; hy += nY)
	{
	    os << "<Piece Extent=\"" << hx - nX << " " << hx << " " << hy - nY
		<< " " << hy <<" 0 0\" Source=\""
		<< "vti/" << "output" << time << "_"  << hx/nX -1 << "_"
		<< hy/nY -1 << ".vti\"/>\n";

	}
    }
    os << "</PImageData>\n";
    os << "</VTKFile>\n";

    os.close();
}

extern "C" void write_vti_cells_(float* u, int* pPosX, int * pPosY, int* time, char* prefix,
	int* lx, int* hx, int* ly, int* hy)
{
    std::stringstream sDest;
    sDest << prefix << "/" << "output" << *time << "_" << *pPosX << "_" << *pPosY << ".vti";
    std::string destination = sDest.str();
    writeVtiCells(u, destination, *lx, *hx, *ly, *hy);
}

extern "C" void write_pvti_cells_( char* destination,
	int* time, int* HX, int* HY, int* nX, int* nY)
{
    writePVtiCells(destination, *time, *HX, *HY, *nX, *nY);
}





/*void writeVtiCells(float* u, std::string destination, int lx, int hx, int ly, int hy)
{
    int ox = 0, oy = 0;
    int sx = 1, sy = 1;
    int size = (hx-lx)*(hy-ly);
    int offset = 0;
    int binSize = size*4 + 4;

    // std::cout << "size = " << size << std::endl;

    // header
    {
	std::ofstream os(destination.c_str());

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	os << "<ImageData WholeExtent=\"" << lx << " " << hx << " " << ly<< " "
	    << hy << " 0 0\" Origin=\"" << ox << " " << oy
	    << " 0\" Spacing=\"" << sx << " " << sy << " 1\">\n";
	os << "<Piece Extent=\"" << lx << " " << hx << " " << ly << " " << hy << " 0 0\">\n";

	os << "<CellData Scalars=\""<< "output" << "\">\n";

	os << "<DataArray type=\"Float32\" Name=\"u\" format=\"appended\" offset=\"" << offset << "\"/>\n";

	os << "</CellData>\n";
	os << "</Piece>\n";
	os << "</ImageData>\n";
	os << "<AppendedData encoding=\"raw\">\n";
	os << "_" ;

	os.close();
    }

    // binary data
    {
	int N = 4*size;
	std::ofstream os(destination.c_str(), std::ios::binary | std::ios::app);
	// write the number of bytes
	os.write(reinterpret_cast<const char*>(&N),sizeof(int));
	for(int i = 0; i < size; ++i) os.write(reinterpret_cast<const char*>(&u[i]), sizeof(float));

	os.close();
    }

    // footer
    {
	std::ofstream os(destination.c_str(), std::ios::app);

	os << "\n</AppendedData>\n";
	os << "</VTKFile>\n";

	os.close();
    }
}

void writePVtiCells(std::string destination, int time,
	int HX, int HY, int nX, int nY)
{
    std::stringstream sFile;
    sFile << destination << "/output" << "_" << time << ".pvti";
    std::ofstream os(sFile.str().c_str());
    os << "<?xml version=\"1.0\" ?>\n";
    os << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "<PImageData WholeExtent=\"0 " << HX << " 0 "
	<< HY << " 0 0" << "\" GhostLevel=\"2\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    os << "<PCellData Scalars=\"" << "output" << "\">\n";

    os << "<PDataArray type=\"Float32\" Name=\"u\"/>\n";

    os << "</PCellData>\n";
    for(int hx = nX; hx <= HX; hx += nX)
    {
	os << "<Piece Extent=\"" << hx - nX << " " << hx << " 0 " << HY <<" 0 0\" Source=\""
	    << "vti/" << "output" << time << "_"  << hx/nX -1 << ".vti\"/>\n";

    }
    os << "</PImageData>\n";
    os << "</VTKFile>\n";

    os.close();
}

extern "C" void write_vti_cells_(float* u, int* pPosX, int* pPosY, int* time, char* prefix,
	int* lx, int* hx, int* ly, int* hy)
{
    std::stringstream sDest;
    sDest << prefix << "/" << "output" << *time << "_" << *pPosX "_" << *pPosY << ".vti";
    std::string destination = sDest.str();
    writeVtiCells(u, destination, *lx, *hx, *ly, *hy);
}

extern "C" void write_pvti_cells_(char* destination, int* time, int* HX, int* HY, int* nX, int* nY)
{
    writePVtiCells(destination, *time, *HX, *HY, *nX, nY);
}*/
