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
#include <string>
#include <stdbool.h>


double E = 540;
double mu = .35;
int width = 640;
int height = 360;
int spaces = 80; // For Perlin Noise granularity

// helper functions:
double separation(double z)
{
	return static_cast<int>(E*((1.-z*mu)/(2.-z*mu))+.5);
}

float interpolate(float x, float y)
{
    float fadeX1 = ( 3.0 - 2.0 * x ) * x * x;
    float fadeY1 = ( 3.0 - 2.0 * y ) * y * y;
    float fadeX0 = 1.0 - fadeX1;
    float fadeY0 = 1.0 - fadeY1;

    float x0 = x;
    float x1 = x - 1.0;
    float y0 = y;
    float y1 = y - 1.0;

    return ( x0 * gradX00 + y0 * gradY00 ) * fadeX0 * fadeY0 +
           ( x1 * gradX10 + y0 * gradY10 ) * fadeX1 * fadeY0 +
           ( x0 * gradX01 + y1 * gradY01 ) * fadeX0 * fadeY1 +
           ( x1 * gradX11 + y1 * gradY11 ) * fadeX1 * fadeY1; 
}

float computePerlinValue( 
	std::vector<unsigned char> &rgbBuffer,
	int width, int height,
	int pixMinX, int pixMinY, 				// Image coordinates of start of current cell
	int pixMaxX, int pixMaxY,				// Image coordinates of end of current cell
{

	float domMinX = 0.0; float domMinY = 0.0;
	float domMaxX = 1.0; float domMaxY = 1.0;

	int pixSpanX = pixMaxX - pixMinX; // spaces
    int pixSpanY = pixMaxY - pixMinY; // spaces

    float pixToDomX = ( domMaxX - domMinX ) * 1.0 / pixSpanX;
    float pixToDomY = ( domMaxY - domMinY ) * 1.0 / pixSpanY;
    
    float pixels[pixSpanX][pixSpanY];

    for ( int pixY=0; pixY<pixSpanY; ++pixY) 
    {
		float domY = pixY * pixToDomY + domMinY;
		for ( int pixX = 0; pixX < pixSpanX; ++pixX ) 
		{
			float domX = pixX * pixToDomX + domMinX;

			float value = std::max( std::min( interpolate( domX, domY ), 1.0 ), -1.0 );
			float grey = 128 + 126 * value | 0;

			rgbBuffer[ 3*((pixMinY+pixY)*width + pixMinX + pixX ) ] 	= grey;
			rgbBuffer[ 3*((pixMinY+pixY)*width + pixMinX + pixX ) + 1 ] = grey;
			rgbBuffer[ 3*((pixMinY+pixY)*width + pixMinX + pixX ) + 2 ] = grey;
		}
    }
}

std::vector<unsigned char> drawPerlinNoise(
	std::vector<unsigned char> &rgbBuffer,	// RGB image buffer to populate
	int width, int height, 	  				// Width and height of image
	int spaces,				  				// Granularity of grid for perlin, ideally divides evenly into width and height of img
	float gradsX[][],		  				// X-Coord for each gradoemt in grid
	float gradsY[][])		  				// Y-Coord for each gradient in grid
{
	int cellMinX = 0;
	int cellMinY = 0;
    int cellMaxX = ceil ( width * 1.0 /  spaces);
    int cellMaxY = ceil ( height * 1.0 /  spaces);

    for(int cellY=cellMinY; cellY<cellMaxY; cellY++)
    {
    	for(int cellX=cellMinX; cellX<cellMaxX; cellX++)
    	{

    		// Get all the gradients for the four current grid vertices
    		int pixX = cellX * space;
	        int pixY = cellY * space;
	        float gradX00 = gradsX[ cellY + 0 ][ cellX + 0 ];
	        float gradY00 = gradsY[ cellY + 0 ][ cellX + 0 ];
	        float gradX10 = gradsX[ cellY + 0 ][ cellX + 1 ];
	        float gradY10 = gradsY[ cellY + 0 ][ cellX + 1 ];
	        float gradX01 = gradsX[ cellY + 1 ][ cellX + 0 ];
	        float gradY01 = gradsY[ cellY + 1 ][ cellX + 0 ];
	        float gradX11 = gradsX[ cellY + 1 ][ cellX + 1 ];
	        float gradY11 = gradsY[ cellY + 1 ][ cellX + 1 ];

    		computePerlinValue(
    			rgbBuffer,
    			width, height,
    			pixX, pixY, 
    			pixX + space, 
    			pixY + space);
    	}
    }
}


int main(int argc, char * argv[])
{

	Camera camera;
	std::vector< std::shared_ptr<Object> > objects;
	// Read a camera and scene description from given .json file
	read_json(argc<=1?"../shared/data/bunny.json":argv[1],camera,objects);

	std::vector<unsigned char> depth_image(width*height);
	std::vector<double> zbuffer(width*height);
	double D = camera.d; // Focal length
	// For each pixel (i,j)
	for(int i=0; i<height; ++i)
	{
		for(int j=0; j<width; ++j)
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

				// used for generating autostereogram
				zbuffer[j+width*i] = linearized_depth;
			}
		}
	}

	// Generate GS noise image
	std::vector<unsigned char> framebuffer(width*height);
	for (int j=0; j<height; j++) 
	{
		for (int i=0; i<width; i++) 
		{
			framebuffer[(i+j*width)] = (rand()%256)/2.0;
		}
	}
	// Generate RGB noise image
	std::vector<unsigned char> framebufferRGB(width*height*3);
	for (int j=0; j<height; j++) 
	{
		for (int i=0; i<width; i++) 
		{
			framebufferRGB[3*(i+j*width)] 	= (rand()%256)/1.1;
			framebufferRGB[3*(i+j*width) + 1] = (rand()%256)/1.1;
			framebufferRGB[3*(i+j*width) + 2] = (rand()%256)/1.1;
		}
	}
	// Generate RGB buffer for Perlin Noise
	int cellsX = width / space;
    int cellsY = height / space;
    // Adding 1 to array size because we are defining gradients for each corner of the grid, inclusive of the endpoints
    float gradsX[cellsY+1][cellsX+1];
    float gradsY[cellsY+1][cellsX+1];

    for(int cellY=0; cellY<=cellsY; cellY++)
    {
    	for(int cellX=0; cellX<=cellsX; cellX++)
    	{
    		var dir = rand()*2.0*M_PI;
	        gradsX[ cellY ][ cellX ] = cos( dir );
	        gradsY[ cellY ][ cellX ] = sin( dir );
    	}
    }
	std::vector<unsigned char> framebufferPerlinRGB(width*height*3);
	drawNoise(
		framebufferPerlinRGB,
		width, height,
		spaces, 
		gradsX, gradsY);


	for(int i = 0; i < height; i++)
	{
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
			if(0<=left && right<width)
			{
				int visible;
				int t=1;
				double zt;
				do
				{
					zt = zbuffer[i*width+j]+2*(2-mu*zbuffer[i*width+j])*t/(mu*E);
					visible = (zbuffer[width*(i-t) + j] < zt) && (zbuffer[width*(i+t) + j] < zt);
					t++;
				} while(visible && zt<1);
				if(visible)
				{
					int l = same[left];
					while(l != left && l != right)
					{
						if( l < right )
						{
							left = l;
							l = same[left];
						}
						else
						{
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
			if(same[j] != j)
			{
				framebuffer[i*width + j] = framebuffer[i*width + same[j]];
				framebufferRGB[3*(i*width + j)] = framebufferRGB[3*(i*width + same[j])];
				framebufferRGB[3*(i*width + j) + 1] = framebufferRGB[3*(i*width + same[j]) + 1];
				framebufferRGB[3*(i*width + j) + 2] = framebufferRGB[3*(i*width + same[j]) + 2];
			}
		}
	}

	
	write_ppm("../images/depth_map_bunny.ppm",depth_image,width,height,1);
	write_ppm("../images/autostereogram_bunny_GS.ppm",framebuffer,width,height,1);
	write_ppm("../images/autostereogram_bunny_RGB.ppm",framebufferRGB,width,height,3);
}
