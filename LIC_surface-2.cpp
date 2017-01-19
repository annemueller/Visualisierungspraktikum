//
// Created by yves on 05.01.17.
//
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fantom/algorithm.hpp>
#include <fantom/register.hpp>
#include <fantom/graphics.hpp>
#include <fantom/fields.hpp>
#include "perlin.h"

using namespace fantom;
namespace {
    class LIC : public VisAlgorithm {
        std::unique_ptr<Primitive> noiseTexture;
        std::unique_ptr<Primitive> licTexture;
      
       

        //const size_t minHitsStreamLines = 1;

        //const size_t L = 40; // wie viele Pixel maximal in Berechnung nutzen?
        //const int B = 20; // wie viele Nachbarn zu {sx,sy} duerfen einen Wert bekommen?
        //const size_t max_hits = 20;

    public:
        struct Options : public VisAlgorithm::Options {
            Options(fantom::Options::Control &control)
                    : VisAlgorithm::Options(control) {
	        add<TensorFieldContinuous<3, Vector3> >("Field", "A 3D vector field");
                add< Grid <3 > >( "Grid", "3D Grid of the ICE" );
                add< int >("Number of pixel", "How many pixel per triangle?", 10, &acceptNumber);

            }

            static int acceptNumber(const int& i){
                return std::max(i , 5);
            }
        };

        struct VisOutputs : public VisAlgorithm::VisOutputs {
            VisOutputs(fantom::VisOutputs::Control &control)
                    : VisAlgorithm::VisOutputs(control) {
                addGraphics("noiseTexture");
		addGraphics("licTexture");
            }
        };

        LIC(InitData &data)
                : VisAlgorithm(data) {
        }

        std::unique_ptr<Texture> createTexture(size_t w, size_t h, bool isNoise = false) {
            std::unique_ptr<Texture> texture = makeTexture(false, w, h, 1);
            //std::shared_ptr< Texture > texture2D = makeTexture( resourcePath() + "/baal.jpg" );
            //std::shared_ptr< Texture > texture2D = makeTexture( resourcePath() + "/aaa.jpg" );
            for (size_t j = 0; j < texture->height(); ++j) {
                for (size_t i = 0; i < texture->width(); ++i) {
                    if (isNoise) {
                        float color = (noise(i, j) + 1.0) * 0.5; // siehe perlin.h/cpp
                        texture->set({color, color, color}, i, j);
                        //Color color = texture2D.get()->get(i, j);
                        //texture->set(color, i, j);
                    } else {
                        texture->set({1.0, 0.0, 0.0}, i, j);
                    }
                }
            }
            return texture;
        }

 
        virtual void execute(const Algorithm::Options &options, const volatile bool &abortFlag) override {
            std::shared_ptr<const TensorFieldContinuous<3, Vector3> > field = options.get<TensorFieldContinuous<3, Vector3> >(
                    "Field");

            if (!field) {
                debugLog() << "Input Field not set!" << std::endl;
                return;
            }

            std::shared_ptr< const Grid < 3 > > grid = options.get< Grid < 3 > >( "Grid" );

            if( !grid )
            {
                throw std::logic_error( "Wrong type of grid!" );
            }
        

            noiseTexture = getGraphics("noiseTexture").makePrimitive();
	    licTexture = getGraphics("licTexture").makePrimitive();

            size_t numCells = grid->numCells();
            int numPixels = options.get<int>("Number of pixel");
           
            std::vector<std::vector<double>> values;
            std::vector<std::vector<double>> hits;

             auto& points = grid->points();
             int numRows = numCells/500 + 1;
             double partition = 1.00/numRows;

            std::shared_ptr< Texture > noise = createTexture(numPixels*500+1, numRows*numPixels+1 , true);
            noiseTexture->setTexture(0, *noise);
            std::shared_ptr<Texture> lic = createTexture(numPixels*500+1, numRows*numPixels+1 , false);
	    licTexture->setTexture(0, *lic);

            // Init: summierte Farbwerte und Anzahl Treffer pro Pixel
            /*for (size_t y = 0; y < numRows*numPixels+1 ; y++) {
                std::vector<double> xval;
                std::vector<double> xhits;
                for (size_t x = 0; x < numPixels*500+1 ; x++) {
                    xval.push_back(0);
                    xhits.push_back(0);
                }
                values.push_back(xval);
                hits.push_back(xhits);
            }

            infoLog() << "Starte LIC" << std::endl;

            // LIC
            for (size_t y = 0; y < numRows*numPixels+1; ++y) {
                for (size_t x = 0; x < numPixels*500+1; ++x) {
                    Point2 fcoord = coord_to_field(x, y);
                    if (hits[y][x] >= minHitsStreamLines) {
                        continue;
                    }
                    std::vector<std::vector<Point3>> points = compute_streamline(field, fcoord[0], fcoord[1]);
                    compute_convolution(noise, points, x, y, values, hits);
                }
            }

            infoLog() << "Berechne Farbwerte" << std::endl;
            for (size_t y = 0; y < numRows*numPixels+1; y++) {
                for (size_t x = 0; x < numPixels*500+1; x++) {
                    // Color = Summierter Farbwert durch Anzahl Treffer
                    float sum = values[y][x] / float(hits[y][x]);
                    /*
                    if (hits[y][x] <= 1 and values[y][x] < 0.1){
                        sum = values[x][y] / double(hits[x][y]);
                    }
                    //double sum = values[y][x] / double(hits[y][x]);
                     */
                    //lic->set({sum, sum, sum}, x, y);
                //}
           // }

        int start = 0;
        int end = 500;
       
        std::vector< Point3 > texCoords;
        std::vector< Point3 > vertices;

        for(int j = 0; j < numRows; j++ ){
          // std::cout << end << std::endl;
          int cellRow = 0;
          
            for(int i = start; i < end ; i++){
          
              Cell cell = grid->cell(i);
              
              vertices.push_back( points[cell.index(0)] );
              vertices.push_back( points[cell.index(2)] );
              vertices.push_back( points[cell.index(1)] );
        
             
              texCoords.push_back( Point3(cellRow/500.00,j*partition,0) );//1
              //std::cout << texCoords.back();
              texCoords.push_back( Point3(cellRow/500.00,j*partition+partition,0) );//2
              //std::cout << texCoords.back();
              texCoords.push_back( Point3((cellRow+1)/500.00, j*partition+partition,0) );//3
              //std::cout << texCoords.back() << std::endl;
              cellRow += 1;
              // std::cout << cellRow << std::endl;
            }
        
            start = end;
            end = start + 500;

            if(end > numCells){
                end = numCells;
            }
        //std::cout << j << std::endl;
        }

            
            noiseTexture->add( Primitive::TRIANGLES ).setTexCoords( 0, texCoords ).setVertices(vertices);            
	    licTexture->add( Primitive::TRIANGLES ).setTexCoords( 0, texCoords ).setVertices(vertices);
        }
    };

  AlgorithmRegister<LIC> reg("VisPraktikum/LIC_surface", "");
}
