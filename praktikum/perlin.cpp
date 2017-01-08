//
// Created by kathleen on 08.01.17.
//
#include "perlin.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <iostream>

#define FASTFLOOR(x) ( ((x)>0) ? ((int)x) : (((int)x)-1) )

std::vector<int> PermutationVektor() {
    std::vector<int> p;
    p.resize(256);
    std::iota(p.begin(), p.end(), 0);
    std::default_random_engine engine(200);
    std::shuffle(p.begin(), p.end(), engine);
    p.insert(p.end(), p.begin(), p.end());
    return p;
}

std::vector<int> p = PermutationVektor(); // global, einmal

float grad(int hash, float x, float y) {
    int h = hash & 7;      // Convert low 3 bits of hash code
    float u = h < 4 ? x : y;  // into 8 simple gradient directions,
    float v = h < 4 ? y : x;  // and compute the dot product with (x,y).
    return ((h & 1) ? -u : u) + ((h & 2) ? -2.0f * v : 2.0f * v);
}

float mynoise(float x, float y) {
    //std::vector<int> p = PermutationVektor();

#define F2 0.366025403f // F2 = 0.5*(sqrt(3.0)-1.0)
#define G2 0.211324865f // G2 = (3.0-Math.sqrt(3.0))/6.0

    float n0, n1, n2; //Beitrag der 3 Ecken des Simplex

    //Gitter verzehren
    float s = (x + y) * F2;
    float xs = x + s;
    float ys = y + s;
    int i = FASTFLOOR(xs);
    int j = FASTFLOOR(ys);


    float t = (float) (i + j) * G2;
    float X0 = i - t;
    float Y0 = j - t;
    float x0 = x - X0;
    float y0 = y - Y0;

    //Welcher Simplex?
    int i1, j1;
    if (x0 > y0)// lower triangle, XY order: (0,0)->(1,0)->(1,1)
    {
        i1 = 1;
        j1 = 0;
    } else// upper triangle, YX order: (0,0)->(0,1)->(1,1)
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
    float t0 = 0.5f - x0 * x0 - y0 * y0;
    if (t0 < 0.0f) {
        n0 = 0.0f;
    } else {
        t0 *= t0;
        n0 = t0 * t0 * grad(p[ii + p[jj]], x0, y0);
    }

    //Farbbeitrag von Ecke 1
    float t1 = 0.5f - x1 * x1 - y1 * y1;
    if (t1 < 0.0f) {
        n1 = 0.0f;
    } else {
        t1 *= t1;
        n1 = t1 * t1 * grad(p[ii + i1 + p[jj + j1]], x1, y1);
    }

    //Farbbeitrag von Ecke 2
    float t2 = 0.5f - x2 * x2 - y2 * y2;
    if (t2 < 0.0f) {
        n2 = 0.0f;
    } else {
        t2 *= t2;
        n2 = t2 * t2 * grad(p[ii + 1 + p[jj + 1]], x2, y2);
    }

    return 40.0f * (n0 + n1 + n2);
}

int main() {


    float t = mynoise(0.3, 0.23);

    std::cout << t << std::endl;

}
