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

        const size_t anzPixelTexture = 400;
        const size_t minHitsStreamLines = 1;

        const size_t L = 40; // wie viele Pixel maximal in Berechnung nutzen?
        const int B = 20; // wie viele Nachbarn zu {sx,sy} duerfen einen Wert bekommen?
        const size_t max_hits = 20;

    public:
        struct Options : public VisAlgorithm::Options {
            Options(fantom::Options::Control &control)
                    : VisAlgorithm::Options(control) {
                add<TensorFieldContinuous<3, Vector3> >("Field", "A 3D vector field");

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

        inline Point2 coord_to_field(double x, double y){
            size_t len = anzPixelTexture / 10;
            double fx = (x / len) * 2.0 - 10;
            double fy = (y / len) * 2.0 - 10;
            return {fx, fy};
        }

        inline Point2 coord_to_texture(Point3 coord){
            size_t len = anzPixelTexture / 10;
            float x = (coord[0] + 10) * len / 2.0;
            float y = (coord[1] + 10) * len / 2.0;
            return {x, y};
        }

        std::vector<std::vector<Point3>> compute_streamline(std::shared_ptr<const TensorFieldContinuous<3, Vector3> > &field, double x, double y) {
            std::vector<Point3> points = runge_kutta(field, 20, 1, {x, y, 0}, true); // foreward
            std::vector<Point3> points2 = runge_kutta(field, 20, 1, {x, y, 0}, false); // bachward
            return {points, points2};
        }

        inline std::vector<std::vector<size_t>> bresenham(size_t x0, size_t y0, size_t x1, size_t y1)
        {
            std::vector<std::vector<size_t>> points;

            int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
            int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
            int err = dx+dy, e2; /* error value e_xy */

            // push 1.

            while(true){
                points.push_back({x0,y0});
                if (x0==x1 && y0==y1) break;
                e2 = 2*err;
                if (e2 > dy) { err += dy; x0 += sx; } /* e_xy+e_x > 0 */
                if (e2 < dx) { err += dx; y0 += sy; } /* e_xy+e_y < 0 */
                //if (x0==x1 && y0==y1) break;
                // push naechsten, verhindert doppelte am Ende

            }
            return points;
        }

        void compute_convolution(std::shared_ptr<Texture> &noise, std::vector<std::vector<Point3>> &pointsStreamline, size_t sx, size_t sy, std::vector<std::vector<double>> &values, std::vector<std::vector<double>> &hits) {
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
                    Point2 p0 = coord_to_texture(pointsStreamline[fb][i]); // Punkt1 aus Stromlinie
                    Point2 p1 = coord_to_texture(pointsStreamline[fb][i+1]); // Punkt2 aus Stromlinie
                    // Coord sind in Texturcoord und erzeuge Liste mit getroffenen Pixeln der Geraden aus p1 und p2

                    //std::vector<std::vector<size_t>> pixel = bresenham(p1[0], p1[1], p2[0], p2[1]);
                    //Point3 col = bresenham(p1[0], p1[1], p2[0], p2[1]);

                    size_t x0 = p0[0]; size_t y0 = p0[1];
                    size_t x1 = p1[0]; size_t y1 = p1[1];

                    int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
                    int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
                    int err = dx+dy, e2; /* error value e_xy */

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
                size_t xpos = pixel_list[ne][0]>=anzPixelTexture ? anzPixelTexture-1 : pixel_list[ne][0];
                size_t ypos = pixel_list[ne][1]>=anzPixelTexture ? anzPixelTexture-1 : pixel_list[ne][1];

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
                    val += noise->get(pixel_list[pos + i][0], pixel_list[pos + i][1], 0).r();
                    anzahl++;
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

        std::vector<Point3>
        runge_kutta(std::shared_ptr<const TensorFieldContinuous<3, Vector3> > &field, size_t steps, double h,
                    Point3 start, bool vor) {
            auto evaluator = field->makeEvaluator();

            std::vector<Point3> points;
            //return points;
            Point3 point = {start[0], start[1], start[2]};
            points.push_back(point);

            for (size_t i = 0; i < steps; i++) {
                if (evaluator->reset(point)) {
                    // Berechne q1-q4
                    auto v = evaluator->value();
                    Point3 q1 = h * v;
                    if (vor == false){
                        q1 = -q1;
                    }
                    if (evaluator->reset(point + 0.5 * q1)) {
                        v = evaluator->value();
                        Point3 q2 = h * v;
                        if (vor ==false){
                            q2 = -q2;
                        }
                        if (evaluator->reset(point + 0.5 * q2)) {
                            v = evaluator->value();
                            Point3 q3 = h * v;
                            if (vor ==false){
                                q3 = -q3;
                            }
                            if (evaluator->reset(point + 0.5 * q3)) {
                                v = evaluator->value();
                                Point3 q4 = h * v;
                                if (vor ==false){
                                    q4 = -q4;
                                }
                                point = point + (1.0 / 6.0) * (q1 + 2 * q2 + 2 * q3 + q4);
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

        virtual void execute(const Algorithm::Options &options, const volatile bool &abortFlag) override {
            std::shared_ptr<const TensorFieldContinuous<3, Vector3> > field = options.get<TensorFieldContinuous<3, Vector3> >(
                    "Field");

            if (!field) {
                debugLog() << "Input Field not set!" << std::endl;
                return;
            }

            noiseTexture = getGraphics("noiseTexture").makePrimitive();
            licTexture = getGraphics("licTexture").makePrimitive();

            std::shared_ptr<Texture> noise = createTexture(anzPixelTexture, anzPixelTexture, true);
            std::shared_ptr<Texture> lic = createTexture(anzPixelTexture, anzPixelTexture, false);

            std::vector<std::vector<double>> values;
            std::vector<std::vector<double>> hits;

            // Init: summierte Farbwerte und ANzahl Treffer pro Pixel
            for (size_t y = 0; y < anzPixelTexture; y++) {
                std::vector<double> xval;
                std::vector<double> xhits;
                for (size_t x = 0; x < anzPixelTexture; x++) {
                    xval.push_back(0);
                    xhits.push_back(0);
                }
                values.push_back(xval);
                hits.push_back(xhits);
            }

            infoLog() << "Starte LIC" << std::endl;

            // LIC
            for (size_t y = 0; y < anzPixelTexture; ++y) {
                for (size_t x = 0; x < anzPixelTexture; ++x) {
                    Point2 fcoord = coord_to_field(x, y);
                    if (hits[y][x] >= minHitsStreamLines) {
                        continue;
                    }
                    std::vector<std::vector<Point3>> points = compute_streamline(field, fcoord[0], fcoord[1]);
                    compute_convolution(noise, points, x, y, values, hits);
                }
            }

            infoLog() << "Berechne Farbwerte" << std::endl;
            for (size_t y = 0; y < anzPixelTexture; y++) {
                for (size_t x = 0; x < anzPixelTexture; x++) {
                    // Color = Summierter Farbwert durch Anzahl Treffer
                    float sum = values[y][x] / float(hits[y][x]);
                    /*
                    if (hits[y][x] <= 1 and values[y][x] < 0.1){
                        sum = values[x][y] / double(hits[x][y]);
                    }
                    //double sum = values[y][x] / double(hits[y][x]);
                     */
                    lic->set({sum, sum, sum}, x, y);
                }
            }

            std::vector<Point3> cube(4);
            cube[0] = Point3(20, -10, 5);
            cube[1] = Point3(40, -10, 5);
            cube[2] = Point3(40, 10, 5);
            cube[3] = Point3(20, 10, 5);

            std::vector<Point3> cube2(4);
            cube2[0] = Point3(-10, -10, -0.5);
            cube2[1] = Point3(10, -10, -0.5);
            cube2[2] = Point3(10, 10, -0.5);
            cube2[3] = Point3(-10, 10, -0.5);

            std::vector<Point3> texCoords(4);
            texCoords[0] = Point3(0.0, 0.0, 1.0);
            texCoords[1] = Point3(1.0, 0.0, 1.0);
            texCoords[2] = Point3(1.0, 1.0, 1.0);
            texCoords[3] = Point3(0.0, 1.0, 1.0);

            std::vector<unsigned int> sides(4);
            sides[0] = 0;
            sides[1] = 1;
            sides[2] = 2;
            sides[3] = 3;

            noiseTexture->setTexture(0, *noise);
            noiseTexture->add(Primitive::QUADS).setTexCoords(0, texCoords).setVertices(cube, sides);

            licTexture->setTexture(0, *lic);
            licTexture->add(Primitive::QUADS).setTexCoords(0, texCoords).setVertices(cube2, sides);
        }
    };

    AlgorithmRegister<LIC> reg("VisPraktikum/LIC", "LIC mit Simplex Noise");
}
