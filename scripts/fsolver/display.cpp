
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "fsolver.h"

FluidSolver *fsolver = NULL;
int displayVelocities = 0;
bool running = false;
bool button1Down = false;
bool button2Down = false;
int lastX = 0, lastY = 0;
double angleX = 0.0f, angleY = 0.0f;
double distZ = -140;
int displayMode = 0;
int width = 768;
int height = 768;
int frameCounter = 0;



void exportDensity(int i) {
	char fname[128], cmd[128];
	snprintf(fname, sizeof(fname), "output/density-%04i.vol", i);
	snprintf(cmd, sizeof(fname), "gzip -f output/density-%04i.vol&", i);
	cout << "    + Saving \"" << fname << "\"" << endl;
	std::ofstream os(fname);

	int xres = fsolver->getGridSizeX(),
	    yres = fsolver->getGridSizeY(),
		zres = fsolver->getGridSizeZ();
	const vector<Density> &density = fsolver->getDensity();

	Float scale = 1.0f / std::max(std::max(xres, yres), zres);

	os.write("VOL", 3);
	char version = 3;
	os.write((char *) &version, sizeof(char));
	int value = 1;
	os.write((char *) &value, sizeof(int));
	os.write((char *) &xres, sizeof(int));
	os.write((char *) &yres, sizeof(int));
	os.write((char *) &zres, sizeof(int));
	value = 1;
	os.write((char *) &value, sizeof(int));
	
	float minX = -xres/2.0f*scale;
	float minY = -yres/2.0f*scale;
	float minZ = -zres/2.0f*scale;
	float maxX = xres/2.0f*scale;
	float maxY = yres/2.0f*scale;
	float maxZ = zres/2.0f*scale;

	os.write((char *) &minX, sizeof(float));
	os.write((char *) &minY, sizeof(float));
	os.write((char *) &minZ, sizeof(float));
	os.write((char *) &maxX, sizeof(float));
	os.write((char *) &maxY, sizeof(float));
	os.write((char *) &maxZ, sizeof(float));
	for (size_t i=0; i<density.size(); ++i) {
		float value = density[i].density;
		os.write((char *) &value, sizeof(float));
	}
	os.close();
	
}


int main(int argc, char **argv) {
	int resX = 80, resY = 80, resZ = 40;
	Float size = 0.5f/60 * 5;
	fsolver = new FluidSolver(resX, resY, resZ, size);
	int frameCounter = 0;
	for(int i= 0; i <= 10; i++) {
		fsolver->step(0.1f);
		exportDensity(frameCounter);
		frameCounter++;
	}
	delete fsolver;
	return 0;
}
