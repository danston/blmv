// Copyright Dmitry Anisimov danston@ymail.com (c) 2017.

// An example on how to compute blended coordinates.

// Local includes
#include "./coords/BlendedR2.hpp"

// Main function.
int main() {

    using namespace gbc;

    // Polygon.
    std::vector<VertexR2> poly(7);

    poly[0] = VertexR2(-0.542, -0.740);
    poly[1] = VertexR2(-0.066, -0.740);
    poly[2] = VertexR2(   0.0, -0.086);
    poly[3] = VertexR2( 0.066, -0.740);
    poly[4] = VertexR2( 0.542, -0.740);
    poly[5] = VertexR2( 0.406,  0.444);
    poly[6] = VertexR2(-0.406,  0.444);

    // Evaluation point.
    VertexR2 query(0.0, 0.2);

    // Storage for the computed blended coordinates.
    std::vector<double> b;

    // Compute blended coordinates.
    BlendedR2 blc(poly);

    blc.setContinuity(2);
    blc.compute(query, b);

    // Output the resulting coordinates.
    std::cout << "\nResult: ";
    for (size_t i = 0; i < b.size(); ++i) std::cout << b[i] << " ";
    std::cout << "\n\n";
}
