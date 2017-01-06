//
// Created by yves on 05.01.17.
//

#include <stdexcept>
#include <vector>
#include <fantom/algorithm.hpp>
#include <fantom/register.hpp>
#include <fantom/graphics.hpp>
#include <fantom/fields.hpp>
#include <random>

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

        std::unique_ptr< Texture > createNoiseTexture(size_t w, size_t h, bool noise=false){
            std::unique_ptr< Texture > texture = makeTexture(false, w, h, 1);
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0, 1);
            for( size_t j = 0; j < texture->height(); ++j )
            {
                for( size_t i = 0; i < texture->width(); ++i )
                {
                    if (noise) {
                        double color = dis(gen);
                        texture->set({color, color, color}, i, j);
                    } else{
                        texture->set({1, 0, 0}, i, j);
                    }
                }
            }
            return texture;
        }

        std::vector<Point3> compute_streamline(std::shared_ptr< const TensorFieldContinuous< 3, Vector3 > > &field, double x, double y){
            //infoLog() << x << ";" << y << std::endl;
            std::vector<Point3> points = runge_kutta(field, 50, 0.1, {x,y,0});
            return points;
        }
int bla = 0;
        double compute_convolution(std::shared_ptr< Texture > &noise, std::vector<Point3> points){
            double sum = 0;
            for(size_t i=0; i<points.size(); i++){
                Point3 p = points[i];
                if (bla == 0) {
                    infoLog() << p[0] << ";" << p[1] << std::endl;
                    infoLog() <<((p[0]+10)*10)/2 << ";" << ((p[1]+10)*10)/2 << std::endl;
                    bla++;
                }
                sum += noise->get(((p[0]+10)*10.)/2.,((p[1]+10)*10.)/2., 0).r();
               // break;
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
            if (bla == 0) {
                infoLog() << "k:" << point[0] << ";" << point[1] << std::endl;
            }
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
                    if (bla == 0) {
                        infoLog()<< "lic: " << ((x/10.)*2.)-10 << ";" << ((y/10.)*2.)-10 << std::endl;
                    }
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

