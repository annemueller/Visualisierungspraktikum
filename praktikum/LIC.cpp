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
#include "runge_kutta.h"

using namespace fantom;
namespace {
    class LIC : public VisAlgorithm {
        std::unique_ptr<Primitive> noiseTexture;
        //std::unique_ptr<Primitive> noiseTexture2D;
        std::unique_ptr<Primitive> licTexture;

        size_t anzPixelTexturex = 0;
        size_t anzPixelTexturey = 0;
        const size_t minHitsStreamLines = 1;

        const size_t L = 40; // wie viele Pixel maximal in Berechnung nutzen?
        const int B = 20; // wie viele Nachbarn zu {sx,sy} duerfen einen Wert bekommen?
        const size_t max_hits = 20;

    public:
        struct Options : public VisAlgorithm::Options {
            Options(fantom::Options::Control &control)
                    : VisAlgorithm::Options(control) {
                add<TensorFieldContinuous<3, Vector3> >("Field", "A 3D vector field");
                add<Grid<3> >("Grid", "A 3D vector field");
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
                //addGraphics("noiseTexture2D");
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
                        float color = (mynoise(i, j) + 1.0) * 0.5; // siehe perlin.h/cpp
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

        inline Point3 coord_to_field(double v, double k, Cell c, const ValueArray<Point3> &points){
            Point3 p0 = points[c.index(0)];
            Point3 p1 = points[c.index(1)];
            Point3 p2 = points[c.index(2)];
            Point3 a = p2-p0;
            Point3 b = p1-p0;
            Point3 x = v * a + k * b + p0;
            return x;
        }

        inline Point3 coord_to_texture(double v, double k, size_t index, double part){
            /*Point3 p0 = {(index%500)/500.0      , (index/500)*part, 0};
            Point3 p1 = {(index%500)/500.0      , (index/500)*part+part, 0};
            Point3 p2 = {((index%500)+1)/500.0  , (index/500)*part+part, 0};
*/
            Point3 p0 = {(index)%500*part      , size_t(index/500)*part, 0};
            Point3 p1 = {(index%500)*part + part     , size_t(index/500)*part, 0};
            Point3 p2 = {((index)%500)*part + 0.5*part  , size_t(index/500)*part+part, 0};
            Point3 b = p2-p0;
            Point3 a = p1-p0;
            Point3 x = v * a + k * b + p0;
            return x;
        }

        Point3 getNormale(Cell c, std::shared_ptr<const Grid<3 >> grid){
            Point3 p0 = grid->points()[c.index(0)];
            Point3 p1 = grid->points()[c.index(1)];
            Point3 p2 = grid->points()[c.index(2)];

            Point3 a = p1 - p0;
            Point3 b = p2 - p0;

            double a1 = a[0];
            double a2 = a[1];
            double a3 = a[2];

            double b1 = b[0];
            double b2 = b[1];
            double b3 = b[2];

            Point3 normale;
            normale[0] = a2*b3 - a3*b2;
            normale[1] = a3*b1 - a1*b3;
            normale[2] = a1*b2 - a2*b1;

            return normale;
        }

        std::vector<std::vector<Point3>> compute_streamline(std::shared_ptr<const TensorFieldContinuous<3, Vector3> > &field, Point3 fcoord, Cell c, std::shared_ptr<const Grid<3 >> grid) {
            Point3 normale = getNormale(c, grid);
            std::vector<Point3> points = runge_kutta(field, 20, 0.04, fcoord, true, normale, grid, c); // foreward
            std::vector<Point3> points2 = runge_kutta(field, 20, 0.04, fcoord, false, normale, grid, c); // bachward
            return {points, points2};
        }

        void compute_convolution(std::shared_ptr<Texture> &noise, std::vector<std::vector<Point3>> &pointsStreamline, size_t sx, size_t sy, std::vector<std::vector<double>> &values, std::vector<std::vector<double>> &hits, size_t index_cell, Cell c, double part, const ValueArray<Point3> &points) {
            /*
            const size_t L = 40; // wie viele Pixel maximal in Berechnung nutzen?
            const size_t B = 15; // wie viele Nachbarn zu {sx,sy} duerfen einen Wert bekommen?
            const size_t max_hits = 10;
             */

            // Berechne von Stromlinie getroffene Pixel, nach Bresenham

            // Liste mit 2 Eintraegen (fore- and backward), die jeweils eine Liste mit Vektoren enthalten, den Punkten
            std::vector<std::vector<std::vector<size_t>>> all_pixels;
            for(size_t fb = 0; fb < 2; fb++) {
                // fb = 0 oder 1 = foreward 0, backward 1
                size_t anz_pixel = 0;
                std::vector<std::vector<size_t>> pixels; // Liste mit Punkten

                for (size_t i = 0; i < pointsStreamline[fb].size() - 1; i++) {
                    Point3 pp0 = points[c.index(0)];
                    Point3 pp1 = points[c.index(1)];
                    Point3 pp2 = points[c.index(2)];
                    Point3 aa = pp2-pp0;
                    Point3 bb = pp1-pp0;

                    double pos_x = pointsStreamline[fb][i][0];
                    double pos_y = pointsStreamline[fb][i][1];
                    double vv = (-bb[1]*pos_x + bb[0]*pos_y + bb[1]*pp0[0] -bb[0]*pp0[1]) / (bb[0]*aa[1]-bb[1]*aa[0]);
                    double kk = (pos_x-pp0[0] - vv*aa[0]) / (bb[0]);

                    Point3 ppp0 = coord_to_texture(vv, kk, index_cell, part);
                    //infoLog() << "pos_x:" << pos_x << " pos_y:" << pos_y << " v:" << vv << " k:" << kk << " p:" << ppp0 << std::endl;

                    pos_x = pointsStreamline[fb][i][0];
                    pos_y = pointsStreamline[fb][i][1];
                    vv = (-bb[1]*pos_x + bb[0]*pos_y + bb[1]*pp0[0] -bb[0]*pp0[1]) / (bb[0]*aa[1]-bb[1]*aa[0]);
                    kk = (pos_x-pp0[0] - vv*aa[0]) / (bb[0]);
                    Point3 ppp1 = coord_to_texture(vv, kk, index_cell, part);

                    //Point2 p0 = coord_to_texture(pointsStreamline[fb][i]); // Punkt1 aus Stromlinie
                    //Point2 p1 = coord_to_texture(pointsStreamline[fb][i+1]); // Punkt2 aus Stromlinie
                    // Coord sind in Texturcoord und erzeuge Liste mit getroffenen Pixeln der Geraden aus p1 und p2

                    //std::vector<std::vector<size_t>> pixel = bresenham(p1[0], p1[1], p2[0], p2[1]);
                    //Point3 col = bresenham(p1[0], p1[1], p2[0], p2[1]);

                    size_t x0 = ppp0[0]; size_t y0 = ppp0[1];
                    size_t x1 = ppp1[0]; size_t y1 = ppp1[1];

                    int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
                    int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
                    int err = dx+dy, e2; // error value e_xy

                    pixels.push_back({x0, y0}); // Startpixel
                    anz_pixel++;

                    while(anz_pixel < L/2){
                        size_t x_old = x0;
                        size_t y_old = y0;
                        e2 = 2*err;
                        if (e2 > dy) { err += dy; x0 += sx; }
                        if (e2 < dx) { err += dx; y0 += sy; }
                        // Ende erreicht? Abbrechen und nicht speichern,
                        // da das der gleiche Punkt vom naechsten Segment am Anfang ist |a...b| - |b...c| usw.
                        if (x0==x1 and y0==y1) break;
                        // push naechsten, verhindert Doppelte
                        if (x0 == x_old and y0 == y_old)
                            continue;
                        pixels.push_back({x0, y0});
                        anz_pixel++;
                    }

                    if(anz_pixel >= L/2){
                        break;
                    }
                }

                all_pixels.push_back(pixels);

            }

            std::vector<std::vector<size_t>> pixel_list = all_pixels[1]; // erst die backward
            int pos_mid = pixel_list.size() -1; // merke 'Mitte'
            // Reihenfolge tauschen, vorher: -1,-2,-3; jetzt: -3,-2,-1
            std::reverse(pixel_list.begin(), pixel_list.end());
            // jetzt foreward hinzufuegen: -3,-2,-1, 2,3,4
            // begin()+1, da dieser Pixel in beiden Listen ist
            if (all_pixels[0].size() > 0){
                pixel_list.insert(pixel_list.end(), all_pixels[0].begin()+1, all_pixels[0].end());
            }

            size_t start = 0;//(pos_mid - B) >= 0 ? pos_mid-B : 0; // linker Rand
            size_t ende = pixel_list.size();//(pos_mid + B) > pixel_list.size() ? pixel_list.size() : pos_mid+B; // rechter Rand

            for(size_t ne=start; ne < ende; ne++) {
                size_t xpos = pixel_list[ne][0]>=anzPixelTexturex ? anzPixelTexturex-1 : pixel_list[ne][0];
                size_t ypos = pixel_list[ne][1]>=anzPixelTexturey ? anzPixelTexturey-1 : pixel_list[ne][1];

                if (hits[ypos][xpos] >= max_hits){
                    //infoLog() << "continue;" << std::endl;
                    continue;
                }

                double val = 0;
                int anzahl = 0;
                int pos = ne;//pos_mid;

                for (int i = -B; i <= B; i++) {
                    if (pos + i < 0 or pos + i > pixel_list.size() - 1)
                        continue;
                    //infoLog() << "size:" << pixel_list.size() << " i:" << i << " pos:" << pos << std::endl;
                    //infoLog() << pixel_list[pos + 1][0] << "; " << pixel_list[pos + i][1] << std::endl;
                    if (pixel_list[pos + i][0] < anzPixelTexturex and pixel_list[pos + i][1] < anzPixelTexturey) {
                        val += noise->get(pixel_list[pos + i][0], pixel_list[pos + i][1], 0).r();
                        anzahl++;
                       // infoLog() << "ja" << std::endl;
                    }
                }
                double sumr = val / anzahl;
                // nun Farbe an Stelle von pos in output speichern, gehe weiter nach links und rechts und mache das gleiche
                if (ne == pos_mid) {
                    values[sy][sx] += sumr;
                    hits[sy][sx] += 1;
                } else {
                    values[ypos][xpos] += sumr;
                    hits[ypos][xpos] += 1;
                }
            }
        }

        virtual void execute(const Algorithm::Options &options, const volatile bool &abortFlag) override {
            std::shared_ptr<const TensorFieldContinuous<3, Vector3> > field = options.get<TensorFieldContinuous<3, Vector3> >(
                    "Field");
            std::shared_ptr<const Grid<3 >> grid = options.get<Grid<3> >("Grid");

            if (!field or !grid) {
                debugLog() << "Input not set!" << std::endl;
                return;
            }

            // Variablen anlegen
            const ValueArray<Point3> &points = grid->points();
            size_t numCells = grid->numCells();
            //numCells = 501;
            int numPixels = options.get<int>("Number of pixel");

            int numRows = numCells/500.0 + 1.0;
            double partition = 1.00/numRows;
            anzPixelTexturex = numPixels*500+1;
            anzPixelTexturey = numRows*numPixels+1;

            infoLog() << "Schritt 1: Texturen erzeugen." << std::endl;

            noiseTexture = getGraphics("noiseTexture").makePrimitive();
            //noiseTexture2D = getGraphics("noiseTexture2D").makePrimitive();
            licTexture = getGraphics("licTexture").makePrimitive();

            std::shared_ptr< Texture > noise = createTexture(anzPixelTexturex, anzPixelTexturey, true);
            noiseTexture->setTexture(0, *noise);
            std::shared_ptr< Texture > noise2D = createTexture(anzPixelTexturex, anzPixelTexturey, true);
            //noiseTexture2D->setTexture(0, *noise2D);
            std::shared_ptr<Texture> lic = createTexture(anzPixelTexturex, anzPixelTexturey, false);
            licTexture->setTexture(0, *lic);

            std::vector<std::vector<double>> values;
            std::vector<std::vector<double>> hits;

            // Init: summierte Farbwerte und Anzahl Treffer pro Pixel = 0
            for (size_t y = 0; y < anzPixelTexturey; y++) {
                std::vector<double> xval;
                std::vector<double> xhits;
                for (size_t x = 0; x < anzPixelTexturex; x++) {
                    xval.push_back(0);
                    xhits.push_back(0);
                }
                values.push_back(xval);
                hits.push_back(xhits);
            }

            infoLog() << "Schritt 2: Starte LIC." << std::endl;

            size_t d2_in = 0;
            size_t d2_out = 0;
            size_t d3_in = 0;
            size_t d3_out = 0;

            size_t akt_cell = 0;
            for (size_t yy=0; yy < 1/*numRows*/; yy++){
                for (size_t xx=0; xx < 500/*500*/; xx++){
                    // numPixels * numPixels fuer das eine Dreieck
                    size_t index_cell = akt_cell++;
                    if(index_cell >= numCells){
                        continue;
                    }

                    // Welches Dreieck ist das?
                    Cell c = grid->cell(index_cell);
                    //infoLog() << "Cell: " << index_cell << std::endl;
                    //infoLog() << "Eckpunkte: " << points[c.index(0)] << ";" << points[c.index(1)] << ";" << points[c.index(2)] << std::endl;
                    Point3 p0 = {xx*numPixels, yy*numPixels, 0};
                    Point3 p1 = {xx*numPixels+numPixels, yy*numPixels, 0};
                    Point3 p2 = {xx*numPixels+size_t(0.5*numPixels), yy*numPixels + numPixels, 0};
                    Point3 b = p2-p0;
                    Point3 a = p1-p0;

                    for (size_t y=0; y <= numPixels; y++){
                        for (size_t x=0; x <= numPixels; x++){
                            size_t pos_x = xx*numPixels+x;
                            size_t pos_y = yy*numPixels+y;
                            //hits[pos_y][pos_x] = -1;

                            double v = (-b[1]*pos_x + b[0]*pos_y + b[1]*p0[0] -b[0]*p0[1]) / (b[0]*a[1]-b[1]*a[0]); // nach Formel berechnet, wenn 0<= v <= 1, dann okay
                            double k = (pos_x-p0[0] - v*a[0]) / (b[0]); // nach Formel berechnet, wenn 0<= k <= 1, dann okay
                            //v = trunc(v*1000) / 1000;
                            //k = trunc(k*1000) / 1000;
                            float delta = 0.02;
                            d2_out++;
                            if (v+delta > 0 and v-delta < 1){
                                if(k+delta > 0 and k-delta < 1){
                                    //hits[pos_y][pos_x] = -1;
                                    Point3 fcoord = coord_to_field(v, k, c, points);
                                    d2_in++;
                                    if (grid->contains(c, fcoord)){
                                        //hits[pos_y][pos_x] += 1;
                                        //values[pos_y][pos_x] += 1;
                                        d3_in++;
                                        std::vector<std::vector<Point3>> points_sl = compute_streamline(field, fcoord, c, grid);
                                        // TODO: in runge kutta abbrechen, wenn punkt nicht mehr in eigentlicher Zelle ist, sondern schon in einer anderen?
                                        compute_convolution(noise, points_sl, pos_x, pos_y, values, hits, index_cell, c, numPixels, points);

/*
                                        infoLog() << "fcoord:    " << fcoord << " x:" << pos_x << " y:" << pos_y << " v:" << v << " k:" << k <<  std::endl;

                                        Point3 pp0 = points[c.index(0)];
                                        Point3 pp1 = points[c.index(1)];
                                        Point3 pp2 = points[c.index(2)];
                                        Point3 aa = pp2-pp0;
                                        Point3 bb = pp1-pp0;

                                        double vv = (-bb[1]*fcoord[0] + bb[0]*fcoord[1] + bb[1]*pp0[0] -bb[0]*pp0[1]) / (bb[0]*aa[1]-bb[1]*aa[0]);
                                        double kk = (fcoord[0]-pp0[0] - vv*aa[0]) / (bb[0]);
                                        Point3 d3_d2 = coord_to_texture(vv, kk, index_cell, partition);
                                        infoLog() << "trans:     " << d3_d2 << " x:" << pos_x << " y:" << pos_y << " v:" << vv << " k:" << kk <<  std::endl;
*/
                                    } else{
                                        d3_out++;
                                    }
                                }
                            }

                        }
                    }
                }
            }

            infoLog() << "d2_in:" << d2_in << std::endl;
            infoLog() << "d2_out:" << d2_out-d2_in << std::endl;
            infoLog() << "d3_in:" << d3_in << std::endl;
            infoLog() << "d3_out:" << d3_out << std::endl;

            infoLog() << "Schritt 3: Berechne Farbwerte." << std::endl;
            for (size_t y = 0; y < anzPixelTexturey; y++) {
                for (size_t x = 0; x < anzPixelTexturex; x++) {
                    // Color = Summierter Farbwert durch Anzahl Treffer

                    //noise2D->set({0, 0, 1}, x, y);
                    float sum = values[y][x] / float(hits[y][x]);
                    lic->set({sum, sum, sum}, x, y);

                }
            }

            infoLog() << "Schritt 4: Erzeuge Texturkoordinaten." << std::endl;
            int start = 0;
            int end = 500;
            std::vector< Point3 > texCoords;
            std::vector< Point3 > vertices;
            for(int j = 0; j < numRows; j++ ){
                int cellRow = 0;

                for(int i = start; i < end ; i++){

                    Cell cell = grid->cell(i);

                    Point3 p0 = points[cell.index(0)];
                    Point3 p1 = points[cell.index(2)];
                    Point3 p2 = points[cell.index(1)];
/*
                    Point3 a = p1 - p0;
                    Point3 b = p2 - p0;
                    Point3 c = p2 - p1;

                    Point3 ls;
                    Point3 pl, pr, p;
                    if (norm(a) > norm(b) and norm(a) > norm(c)){
                        // a
                        ls = a;
                        pl = p0;
                        pr = p1;
                        p = p2;
                    } else if (norm(b) > norm(a) and norm(b) > norm(c)){
                        // b
                        ls = b;
                        pl = p0;
                        pr = p2;
                        p = p1;
                    } else{
                        // c
                        ls = c;
                        pl = p1;
                        pr = p2;
                        p = p0;
                    }

                    Point3 k = p-pl;
                    k = k / norm(k);
                    double kx = k[0] / 500;
                    kx = kx < 0 ? -kx : kx;
                    double ky = k[1] / numRows;
                    ky = ky < 0 ? -ky : ky;
*/
                    vertices.push_back( points[cell.index(0)] );
                    vertices.push_back( points[cell.index(2)] );
                    vertices.push_back( points[cell.index(1)] );

                    texCoords.push_back( Point3(cellRow/500.00 + 0.0002  ,j*partition,0) );//1
                    //infoLog() << texCoords.back() << std::endl;
                    texCoords.push_back( Point3((cellRow+1)/500.00 - 0.0002, j*partition,0) );//2
                    //infoLog() << texCoords.back() << std::endl;
                    texCoords.push_back( Point3((cellRow)/500.00 + 0.001, j*partition+partition - (0.134/numRows), 0) );//3
                    //infoLog() << texCoords.back() << std::endl;
                    cellRow += 1;
                }

                start = end;
                end = start + 500;

                if(end > numCells){
                    end = numCells;
                }
            }
/*
            std::vector<Point3> cube(4);
            cube[0] = Point3(0, 0, 5);                               cube[1] = Point3(anzPixelTexturex, 0, 5);
            cube[2] = Point3(anzPixelTexturex, anzPixelTexturey, 5); cube[3] = Point3(0, anzPixelTexturey, 5);
            std::vector<unsigned int> sides(4);
            sides[0] = 0; sides[1] = 1;
            sides[2] = 2; sides[3] = 3;
            std::vector<Point3> texCoords2D(4);
            texCoords2D[0] = Point3(0.0, 0.0, 1.0); texCoords2D[1] = Point3(1.0, 0.0, 1.0);
            texCoords2D[2] = Point3(1.0, 1.0, 1.0); texCoords2D[3] = Point3(0.0, 1.0, 1.0);
*/
            noiseTexture->add( Primitive::TRIANGLES ).setTexCoords( 0, texCoords ).setVertices(vertices);
            //noiseTexture2D->add(Primitive::QUADS).setTexCoords(0, texCoords2D).setVertices(cube, sides);;
            licTexture->add( Primitive::TRIANGLES ).setTexCoords( 0, texCoords ).setVertices(vertices);
        }
    };

    AlgorithmRegister<LIC> reg("VisPraktikum/LIC", "LIC mit Simplex Noise");
}
