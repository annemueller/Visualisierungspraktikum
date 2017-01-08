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

#define FASTFLOOR(x) ( ((x)>0) ? ((int)x) : (((int)x)-1) )

using namespace fantom;
namespace {
    class LIC : public VisAlgorithm {
        std::unique_ptr<Primitive> noiseTexture;
        std::unique_ptr<Primitive> licTexture;

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
            //std::shared_ptr< Texture > texture2D = makeTexture( resourcePath() + "/perlin.jpg" );
            for (size_t j = 0; j < texture->height(); ++j) {
                for (size_t i = 0; i < texture->width(); ++i) {
                    if (isNoise) {
                        float color = (mynoise(i, j) + 1.0) * 0.5; // siehe perlin.h/cpp
                        texture->set({color, color, color}, i, j);
                        //Color color = texture2D.get()->get(i, j);
                        //texture->set(color, i, j);
                    } else {
                        texture->set({1, 0, 0}, i, j);
                    }
                }
            }
            return texture;
        }


        std::vector<Point3> compute_streamline(std::shared_ptr<const TensorFieldContinuous<3, Vector3> > &field, double x, double y) {
            std::vector<Point3> points = runge_kutta(field, 50, 0.1, {x, y, 0}, true); // foreward
            std::vector<Point3> points2 = runge_kutta(field, 50, 0.1, {x, y, 0}, false); // bachward
            points.insert(points.end(), points2.begin(), points2.end());
            return points;
        }

        Point3 compute_convolution(std::shared_ptr<Texture> &noise, std::vector<Point3> points) {
            float sumr = 0;
            float sumg = 0;
            float sumb = 0;
            if (points.size() < 40){
                return {1,0,0};
            }
            int k = points.size() / 2.0;
            for (size_t i = 0; i < 40/*points.size()*/; i++) {
                Point3 p = points[k-20+i];
                sumr += noise->get(((p[0] + 10) * 40.) / 2., ((p[1] + 10) * 40.) / 2., 0).r();
                sumg += noise->get(((p[0] + 10) * 40.) / 2., ((p[1] + 10) * 40.) / 2., 0).g();
                sumb += noise->get(((p[0] + 10) * 40.) / 2., ((p[1] + 10) * 40.) / 2., 0).b();
            }
            sumr /= 40;//float(points.size());
            sumg /= 40;//float(points.size());
            sumb /= 40;//float(points.size());

            return {sumr, sumg, sumb};
        }

        std::vector<Point3>
        runge_kutta(std::shared_ptr<const TensorFieldContinuous<3, Vector3> > &field, size_t steps, double h,
                    Point3 start, bool vor) {
            // foreward
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

            std::shared_ptr<Texture> noise = createTexture(400, 400, true);
            std::shared_ptr<Texture> lic = createTexture(400, 400, false);

            // LIC
            for (size_t y = 0; y < lic->height(); ++y) {
                for (size_t x = 0; x < lic->width(); ++x) {
                    std::vector<Point3> points = compute_streamline(field, ((x / 40.) * 2.) - 10, ((y / 40.) * 2.) - 10);
                    Point3 sum = compute_convolution(noise, points);
                    lic->set({sum[0], sum[1], sum[2]}, x, y);
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
