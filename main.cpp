#include "Object.h"
#include "Ray.h"
#include "Camera.h"
#include "Sphere.h"
#include "Plane.h"
#include "read_json.h"
#include "write_ppm.h"
#include "viewing_ray.h"
#include "first_hit.h"
#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>
#include <limits>
#include <functional>
#include <math.h>


double E = 540;
double mu = .35;

// helper functions:
double separation(double z)
{
	return static_cast<int>(E*((1.-z*mu)/(2.-z*mu))+.5);
}


int main(int argc, char * argv[])
{

	Camera camera;
	std::vector< std::shared_ptr<Object> > objects;
	// Read a camera and scene description from given .json file
	read_json(argc<=1?"../shared/data/bunny.json":argv[1],camera,objects);

	int width = 640;
	int height = 360;
	std::vector<unsigned char> depth_image(width*height);
	std::vector<double> zbuffer(width*height);

	double D = camera.d; // Focal length

	// For each pixel (i,j)
	for(unsigned i=0; i<height; ++i)
	{
		for(unsigned j=0; j<width; ++j)
		{
			// Set background color
			depth_image[j+width*i] = 0;
			// Compute viewing ray
			Ray ray;
			viewing_ray(camera,i,j,width,height,ray);
			// Find first visible object hit by ray and its surface normal n
			double t;
			Eigen::Vector3d n;
			int hit_id;

			if(first_hit(ray,1.0,objects,hit_id,t,n))
			{
				// depth image
				const double zNear = camera.d;
				double linearized_depth = zNear/(t);
				linearized_depth = linearized_depth<1?linearized_depth:1;
				depth_image[j+width*i] = 255.0*(linearized_depth);

				// autostereogram
				zbuffer[j+width*i] = linearized_depth;
			}
		}
	}

	// Generate a noise image
	std::vector<unsigned char> framebuffer(width*height);
	for (int j=0; j<height; j++) { // generate a random-ish image
		for (int i=0; i<width; i++) {
			framebuffer[(i+j*width)] = (rand()%256)/2.0;
		}
	}

	for(int i = 0; i < height; i++)
	{
		// std::cout << "Entered row " << i << "\n";
		int same[width];
		int s;
		int left, right;
		// Initialize;
		for(int j = 0; j < width; j++)
		{
			same[j] = j;
		}

		for(int j = 0; j < width; j++)
		{

			s = separation(zbuffer[j+width*i]);
			left = j - s/2;
			right = left + s;
			if(0<=left && right<width){
				int visible;
				int t=1;
				float zt;
				do{
					zt = zbuffer[i*width+j]+2*(2-mu*zbuffer[i*width+j])*t/(mu*E);
					visible = zbuffer[width*(i-t) + j]<zt && zbuffer[width*(i+t) + j]<zt;
					t++;
				} while(visible && zt<1);
				if(visible){
					int l = same[left];
					while(l != left && l != right){
						if( l < right ){
							left = l;
							l = same[left];
						}
						else{
							same[left] = right;
							left = right;
							l = same[left];
							right = l;
						}
					}
				same[left] = right;
				}
			}
		}

		for(int j = width-1; j>=0; j--)
		{
			if(same[j] != j){
				// std::cout << "got here";
				framebuffer[i*width + j] = framebuffer[i*width + same[j]];
			}
		}
	}
	write_ppm("depth.ppm",depth_image,width,height,1);
	write_ppm("autostereogram.ppm",framebuffer,width,height,1);
}
