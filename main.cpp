// SPDX-License-Identifier: MPL-2.0
/**
 * @file main.cpp
 * @brief edge drawing program
 * @author Adle Ben Salem
 *         Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2025-2026
 */

#include "edge_drawing.h"
#include "image.h"
#include "cmdLine.h"
#include "io_png.h"
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>

using uchar = unsigned char;
struct Color {
    uchar t[3];
    Color() { for(int i=0; i<3; i++) t[i]= uchar(rand()%256); }
};

/// Draw a disk centered at \a p in color image \a out.
/// Color and radius are kept constant inside the function.
void draw_disk(const Point& p, uchar* out, int w, int h) {
    const uchar c[3]={(uchar)128, (uchar)128, (uchar)128};
    const int RAD=4;
    for(int i=-RAD; i<=RAD; i++)
        for(int j=-RAD; j<=RAD; j++) {
            int x=p.x+j, y=p.y+i;
            if(0<=x && x<w && 0<=y && y<h && i*i+j*j<=RAD*RAD)
                std::copy_n(c, 3, out+3*(x+w*y));
        }
}

/// Return positive value different from \a iref at one 8-neighbor of \a p in
/// image \a trace if it exists, otherwise -1.
/// Used for T-junction detection: a pixel of \a trace contains the index of
/// the edge containg it (or -1 if none) and \a p is an endpoint of the edge
/// of index \a iref.
int check_tjunc(const Image<int>& trace, const Point& p, int iref) {
    for(int i=-1; i<=1; i++)
        for(int j=-1; j<=1; j++)
            if(trace.inside(p.x+j,p.y+i)) {
                int k = trace(p.x+j,p.y+i);
                if(k!=-1 && k!=iref)
                    return k;
            }
    return -1;
}

/// Compute gradient of image and store in polar. Gradient magnitude is
/// set to closest integer.
void grad(Image<float> I[3], size_t c, Image<int>& G, Image<float>& Theta) {
    G.reset(I[0].w, I[0].h); G.fill(0);
    Theta.reset(I[0].w, I[0].h); Theta.fill(0);
    if(c==1) // 2x2 scheme, like LSD
        for(int y=0; y+1<G.h; y++)
            for(int x=0; x+1<G.w; x++) {
                float c1 = I[0](x+1,y+1) - I[0](x,y);
                float c2 = I[0](x+1,y) - I[0](x,y+1);
                float gx = c1+c2, gy = c1-c2;
                G(x,y) = (int)std::round(0.5f*std::hypot(gx, gy));
                if(G(x,y)>0)
                    Theta(x,y) = std::atan2(gy, gx);
            }
    if(c==3) { // Di Zenzo color gradient
        const float norm=1/std::sqrt(3);
        for(int y=0; y+1<G.h; y++)
            for(int x=0; x+1<G.w; x++) {
                float dx[3], dy[3];
                for(size_t i=0; i<c; i++) {
                    float c1 = I[i](x+1,y+1) - I[i](x,y);
                    float c2 = I[i](x+1,y) - I[i](x,y+1);
                    dx[i] = c1+c2; dy[i] = c1-c2;
                }
                float gxx = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]; // u.u
                float gyy = dy[0]*dy[0] + dy[1]*dy[1] + dy[2]*dy[2]; // v.v
                float gxy = dx[0]*dy[0] + dx[1]*dy[1] + dx[2]*dy[2]; // u.v
                float theta = 0.5f * std::atan2(2*gxy, gxx-gyy);
                Theta(x,y) = theta;
                float v = std::sqrt((gxx+gyy +
                                     (gxx-gyy)*cos(2*theta) +
                                     2*gxy*sin(2*theta))*0.5f) * norm;
                G(x,y) = (int)std::round(v);
            }
    }
}

int main(int argc, char **argv) {
    int gradMin=6, anchorGap=2;
    int lengthMin=10;
    float sigma=1.0f;
    std::string tjunc;
    float lEpsNFA=0;
    int segLevel=0;

    CmdLine cmd;
    cmd.add( make_option('g', gradMin, "grad-min")
             .doc("Min gradient") );
    cmd.add( make_option('a', anchorGap, "anchor-gap")
             .doc("Min gap of gradient for anchor") );
    cmd.add( make_option('l', lengthMin, "length-min")
             .doc("Min length of edge segment") );
    cmd.add( make_option('s', sigma, "sigma")
             .doc("Sigma of Gaussian blur") );
    cmd.add( make_option('t', tjunc, "t-junctions")
             .doc("T-junctions output image") );
    cmd.add( make_option('e', lEpsNFA, "epsNFA")
             .doc("log10(NFA) for validation, normally 0 or negative") );
    cmd.add( make_option('S', segLevel, "segLevel")
             .doc("Segmentation level for valid edges") );
    try {
        cmd.process(argc, argv);
    } catch(const std::string& s) {
        std::cerr << "Error: " << s << std::endl;
        return 1;
    }
    if(argc!=3 && (argc!=2 || !cmd.used('t'))) {
        std::cerr << "Usage: "<<argv[0] << " [options] in.png [out.png]\n" <<cmd
                  << "out.png may be omitted if -t or --t-junctions was used\n"
                  << "NFA validation only with -e and/or -S\n"
                  << "segLevel:\n"
                  << " 0: validate only full edges\n"
                  << " 1: validate full edge if a segment is valid\n"
                  << " other: validate only segments" << std::endl;
        return 1;
    }
    if(gradMin < 1) {
        std::cerr << "The grad-min parameter must be 1 or more" <<std::endl;
        return 1;
    }
    if(sigma < 0) {
        std::cerr << "The sigma parameter must be 0 or more" <<std::endl;
        return 1;
    }

    float* im; size_t w, h, c;
    im = io_png_read_f32(argv[1], &w, &h, &c);
    if(! im) {
        std::cerr << "Error reading image file " << argv[1] << std::endl;
        return 1;
    }

    Image<float> channels[3];
    for(size_t i=0; i<c; i++) {
        channels[i].reset(w,h);
        channels[i].read(im+i, c);
        blurGaussian(channels[i], sigma);
    }
    free(im);

    Image<int> G; Image<float> Theta;
    grad(channels, c, G, Theta);

    ED ed(G, Theta, gradMin, anchorGap, lengthMin);
    if(cmd.used('e') || cmd.used('S'))
        ed.validateNFA(lEpsNFA, segLevel);

    uchar* out = 0;
    if(argc>2) {
        out = new uchar[3*w*h];
        std::fill_n(out, 3*w*h, 0);
    }
    const std::vector<std::vector<Point>>& seg = ed.edges;
    Color* mapCol=0; bool* keep=0; Image<int> trace; uchar* outT=0; // T-junc
    if(! tjunc.empty()) {
        mapCol = new Color[seg.size()];
        keep = new bool[seg.size()];
        std::fill_n(keep, seg.size(), false);
        trace.reset(w,h);
        trace.fill(-1);
        outT = new uchar[3*w*h];
        std::fill_n(outT, 3*w*h, 0);
    }
    for(size_t i=0; i<seg.size(); ++i) {
        const std::vector<Point>& segi = seg[i];
        Color c = (mapCol? mapCol[i]: Color());
        for(const Point& p : segi) {
            size_t idx = p.y*w + p.x;
            if(out)
                std::copy_n(c.t, 3, out+3*idx);
            if(outT)
                trace(p) = (int)i;
        }
    }

    if(outT) { // Handle T-junctions
        std::vector<Point> junctions;
        for(size_t i=0; i<seg.size(); ++i) { // Record T-junc and involved edges
            int j = check_tjunc(trace, seg[i].front(), i);
            if(j>=0) {
                junctions.push_back(seg[i].front());
                keep[i]=keep[j]=true;
            }
            j = check_tjunc(trace, seg[i].back(), i);
            if(j>=0) {
                junctions.push_back(seg[i].back());
                keep[i]=keep[j]=true;
            }
        }
        for(size_t i=0; i<seg.size(); ++i) // Draw involved edges
            if(keep[i]) {
                const std::vector<Point>& segi = seg[i];
                for(const Point& p : segi) {
                    size_t idx = p.y*w + p.x;
                    std::copy_n(mapCol[i].t, 3, outT+3*idx);
                }
            }
        delete [] keep;
        delete [] mapCol;
        for(size_t i=0; i<junctions.size(); i++) // Draw junctions
            draw_disk(junctions[i], outT, w, h);
    }

    if(out && io_png_write_u8(argv[2], out, w, h, 3)) {
        std::cerr << "Error writing file " << argv[2] << std::endl;
        return 1;
    }
    if(outT && io_png_write_u8(tjunc.c_str(), outT, w, h, 3)) {
        std::cerr << "Error writing file " << tjunc << std::endl;
        return 1;
    }
    delete [] out;
    delete [] outT;
    return 0;
}
