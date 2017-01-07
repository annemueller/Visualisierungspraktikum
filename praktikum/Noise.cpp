//Perlin Noise Textur

//
// Created by yves on 05.01.17.
//
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <fantom/algorithm.hpp>
#include <fantom/register.hpp>
#include <fantom/graphics.hpp>
#include <fantom/fields.hpp>
#include <random>

#define FASTFLOOR(x) ( ((x)>0) ? ((int)x) : (((int)x)-1) )

using namespace fantom;
namespace
{
    class LIC : public VisAlgorithm
    {
        std::unique_ptr< Primitive > noiseTexture;
        std::unique_ptr< Primitive > licTexture;

    public:
        struct Options : public VisAlgorithm::Options
        {
            Options( fantom::Options::Control& control )
                    : VisAlgorithm::Options( control )
            {
                add< TensorFieldContinuous< 3, Vector3 > >( "Field", "A 3D vector field" );

            }
        };
        struct VisOutputs : public VisAlgorithm::VisOutputs
        {
            VisOutputs( fantom::VisOutputs::Control& control )
                    : VisAlgorithm::VisOutputs( control )
            {
                addGraphics( "noiseTexture" );
                addGraphics( "licTexture" );
            }
        };

        LIC( InitData& data )
                : VisAlgorithm( data )
        {
        }
/*--------------------------------------------------------------------------------------------
--------------------Simplex Noise-------------------------------------------------------------
--------------------------------------------------------------------------------------------*/
		std::vector<int>  PermutationVektor(){

			std::vector<int> p;
			p.resize(256);
			std::iota(p.begin(), p.end(), 0);
			std::default_random_engine engine(200);
			std::shuffle(p.begin(), p.end(), engine);
			p.insert(p.end(), p.begin(), p.end());
			return p;
		}

		float  grad( int hash, float x, float y ) {
		    int h = hash & 7;      // Convert low 3 bits of hash code
		    float u = h<4 ? x : y;  // into 8 simple gradient directions,
		    float v = h<4 ? y : x;  // and compute the dot product with (x,y).
		    return ((h&1)? -u : u) + ((h&2)? -2.0f*v : 2.0f*v);
		}

		float noise(float x, float y){
			std::vector<int> p = PermutationVektor(); 

			#define F2 0.366025403f // F2 = 0.5*(sqrt(3.0)-1.0)
			#define G2 0.211324865f // G2 = (3.0-Math.sqrt(3.0))/6.0

			float n0, n1, n2; //Beitrag der 3 Ecken des Simplex

			//Gitter verzehren
			float s = (x+y)*F2;
			float xs = x + s;
			float ys = y + s;
			int i = FASTFLOOR(xs);
			int j = FASTFLOOR(ys);


			float t = (float)(i+j)*G2;
			float X0 = i-t;
			float Y0 = j-t;
			float x0 = x - X0;
			float y0 = y - Y0;

			//Welcher Simplex?
			int i1, j1;
			if (x0>y0)// lower triangle, XY order: (0,0)->(1,0)->(1,1)
			{
				i1 = 1;
				j1 = 0;
			}
			else// upper triangle, YX order: (0,0)->(0,1)->(1,1)
			{
				i1 = 0;
				j1 = 1;
			}

			// Wrap the integer indices at 256, to avoid indexing p[] out of bounds
		    int ii = i & 0xff;
		    int jj = j & 0xff;

		    float x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
		    float y1 = y0 - j1 + G2;
		    float x2 = x0 - 1.0f + 2.0f * G2; // Offsets for last corner in (x,y) unskewed coords
		    float y2 = y0 - 1.0f + 2.0f * G2;

		    //Farbbeitrag von Ecke 0
		    float t0 = 0.5f - x0*x0 - y0*y0;
		    if (t0 < 0.0f)
		    {
		    	n0 = 0.0f;		
		    }
		    else
		    {
		    	t0 *= t0;
		    	n0 = t0 * t0 * grad(p[ii+p[jj]], x0, y0);
		    }

		    //Farbbeitrag von Ecke 1
		    float t1 = 0.5f - x1*x1 - y1*y1;
		    if (t1 < 0.0f)
		    {
		    	n1 = 0.0f;		
		    }
		    else
		    {
		    	t1 *= t1;
		    	n1 = t1 * t1 * grad(p[ii+i1+p[jj+j1]], x1, y1);
		    }

		    //Farbbeitrag von Ecke 2
		    float t2 = 0.5f - x2*x2 - y2*y2;
		    if (t2 < 0.0f)
		    {
		    	n2 = 0.0f;		
		    }
		    else
		    {
		    	t2 *= t2;
		    	n2 = t2 * t2 * grad(p[ii+1+p[jj+1]], x2, y2);
		    }

			return 40.0f * (n0 + n1 + n2); 
		}

        std::unique_ptr< Texture > createNoiseTexture(size_t w, size_t h, bool noise=false){
            std::unique_ptr< Texture > texture = makeTexture(false, w, h, 1);

           	
            for( size_t j = 0; j < texture->height(); ++j )
            {
                for( size_t i = 0; i < texture->width(); ++i )
                {
                    if (noise) {
                        double color = noise(i,j);
                        texture->set({color, color, color}, i, j);
                    } else{
                        texture->set({1, 0, 0}, i, j);
                    }
                }
            }

            
            return texture;
        }


        std::vector<Point3> compute_streamline(std::shared_ptr< const TensorFieldContinuous< 3, Vector3 > > &field, double x, double y){
            std::vector<Point3> points = runge_kutta(field, 50, 0.1, {x,y,0});
            return points;
        }

        double compute_convolution(std::shared_ptr< Texture > &noise, std::vector<Point3> points){
            double sum = 0;
            for(size_t i=0; i<points.size(); i++){
                Point3 p = points[i];
                sum += noise->get(((p[0]+10)*10.)/2.,((p[1]+10)*10.)/2., 0).r();
            }
            sum /= points.size();
            return sum;
        }

        std::vector<Point3> runge_kutta(std::shared_ptr< const TensorFieldContinuous< 3, Vector3 > > &field, size_t steps, double h, Point3 start){
            // foreward
            auto evaluator = field->makeEvaluator();

            std::vector<Point3> points;
            //return points;
            Point3 point = start;
            points.push_back(point);

            for (size_t i = 0; i < steps; i++) {
                if (evaluator->reset(point)) {
                    // Berechne q1-q4
                    auto v = evaluator->value();
                    Point3 q1 = h * v;
                    if (evaluator->reset(point + 0.5 * q1)) {
                        v = evaluator->value();
                        Point3 q2 = h * v;
                        if (evaluator->reset(point + 0.5 * q2)) {
                            v = evaluator->value();
                            Point3 q3 = h * v;
                            if (evaluator->reset(point + 0.5 * q3)) {
                                v = evaluator->value();
                                Point3 q4 = h * v;
                                point = point + (1.0 / 6.0) * (q1 + 2 * q2 + 2 * q3 + q4);
                                //infoLog()<< int(point[0]) << ";" << int(point[1]) << ";" << point[2] << std::endl;
                                points.push_back(point);
                            } else {
                                break; // nicht mehr in Domain
                            }

                        } else {
                            break; // nicht mehr in Domain
                        }
                    } else {
                        break; // nicht mehr in Domain
                    }
                } else {
                    break; // nicht mehr in Domain
                }
            }
            // backward?
            return points;
        }

        virtual void execute( const Algorithm::Options& options, const volatile bool& abortFlag ) override
        {
            std::shared_ptr< const TensorFieldContinuous< 3, Vector3 > > field= options.get< TensorFieldContinuous< 3, Vector3 > >( "Field" );

            if( !field )
            {
                debugLog() << "Input Field not set!" << std::endl;
                return;
            }

            noiseTexture = getGraphics( "noiseTexture" ).makePrimitive();
            licTexture = getGraphics( "licTexture" ).makePrimitive();

            std::shared_ptr< Texture > noise = createNoiseTexture(100, 100, true);
            std::shared_ptr< Texture > lic = createNoiseTexture(100, 100, false);

            // LIC
            for( size_t y = 0; y < lic->height(); ++y ) {
                for (size_t x = 0; x < lic->width(); ++x){
                    std::vector<Point3> points = compute_streamline(field, ((x/10.)*2.)-10, ((y/10.)*2.)-10);
                    double sum = compute_convolution(noise, points);
                    lic->set( {sum, sum, sum}, x, y);
                }
            }

            std::vector< Point3 > cube( 4 );
            cube[0] = Point3( -10, -10, 15 );
            cube[1] = Point3( 10, -10, 15 );
            cube[2] = Point3( 10, 10, 15 );
            cube[3] = Point3( -10, 10, 15 );

            std::vector< Point3 > cube2( 4 );
            cube2[0] = Point3( -10, -10, -0.5 );
            cube2[1] = Point3( 10, -10, -0.5 );
            cube2[2] = Point3( 10, 10, -0.5 );
            cube2[3] = Point3( -10, 10, -0.5 );

            std::vector< Point3 > texCoords( 4 );
            texCoords[0] = Point3( 0.0, 0.0, 1.0 );
            texCoords[1] = Point3( 1.0, 0.0, 1.0 );
            texCoords[2] = Point3( 1.0, 1.0, 1.0 );
            texCoords[3] = Point3( 0.0, 1.0, 1.0 );

            std::vector< unsigned int > sides( 4 );
            sides[0] = 0;
            sides[1] = 1;
            sides[2] = 2;
            sides[3] = 3;

            noiseTexture->setTexture( 0, *noise );
            noiseTexture->add( Primitive::QUADS ).setTexCoords( 0, texCoords ).setVertices( cube, sides );

            licTexture->setTexture( 0, *lic );
            licTexture->add( Primitive::QUADS ).setTexCoords( 0, texCoords ).setVertices( cube2, sides );
        }
    };
    AlgorithmRegister< LIC > reg( "VisPraktikum/LIC", "LIC");
}

