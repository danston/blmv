// Copyright Dmitry Anisimov danston@ymail.com (c) 2017.

// README:
/*

    Blended coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp
    5. Halfedge.hpp
    6. Face.hpp
    7. MeshR2.hpp
    8. TriangulatorR2.hpp
    9. TriangleCoordinatesR2.hpp

*/

#ifndef GBC_BLENDEDR2_HPP
#define GBC_BLENDEDR2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/Face.hpp"
#include "../extra/MeshR2.hpp"
#include "../extra/Halfedge.hpp"
#include "../extra/VertexR2.hpp"
#include "../extra/TriangulatorR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // Blended coordinates in R2.
    class BlendedR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        BlendedR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol), _continuity(1) {

            // Triangulate polygon.
            std::vector<VertexR2> tp;
            std::vector<Face> tf;

            TriangulatorR2 tri(_v);

            tri.setPlanarGraph(true);
            tri.allowSteinerPoints(false);

            tri.generate(tp, tf);

            assert(tp.size() != 0 && tf.size() != 0);

            _mesh.initialize(tp, tf);
            _mesh.createFaces();

            // Preinitialize some vectors.
            _ver.resize(6);
            _ind.resize(6);

            _s.resize(9);
            _r.resize(9);
            _A.resize(9);
            _D.resize(9);
            _t.resize(9);
            _w.resize(9);

            _w0.resize(4);
            _w1.resize(4);
            _w2.resize(4);
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "BlendedR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // D. Anisimov, D. Panozzo, and K. Hormann. Blended barycentric coordinates.
        // Computer Aided Geometric Design, to appear, 2017 (see Appendix A).
        void compute(const VertexR2 &p, std::vector<double> &b) {

            b.clear();
            b.resize(_v.size(), 0.0);

            // Boundary.
            // Comment this line, if you want the local behaviour of 
            // the blended coordinates for interior query points!
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            const int triInd = _mesh.findFace(p); // currently the slow algorithm is used
            assert(triInd != -1);

            const size_t tri = (size_t) triInd;

            std::vector<size_t> neighs;
            neighs.reserve(3);

            const size_t k = _mesh.getFaceNeighbours(tri, neighs); // these neighbours are obtained through halfedges

            switch (k) {
                case 1:
                    caseOneQuadrilateral(p, neighs, b);
                    break;
                case 2:
                    caseTwoQuadrilaterals(p, neighs, b);
                    break;
                case 3:
                    caseThreeQuadrilaterals(p, neighs, b);
                    break;
                default:
                    break;
            }
        }

        // Compute the coordinates at p using the internal storage from the VertexR2 class.
        inline void compute(VertexR2 &p) {
            compute(p, p.b());
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb) {

            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) compute(p[i], bb[i]);
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p) {

            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) compute(p[i], p[i].b());
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

        // Set continuity of the resulting blended coordinates.
        // C1 continuity = 1, C2 continuity = 2.
        // The default continuity is 1 set in the class constructor.
        inline void setContinuity(const size_t continuity) {

            assert(continuity == 1 || continuity == 2);
            _continuity = continuity;
        }

        // Get currently used continuity of the coordinates.
        inline size_t continuity() const {
            return _continuity;
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        // Triangle mesh.
        MeshR2 _mesh;

        // Continuity of the coordinates.
        size_t _continuity;

        // Fields.
        std::vector<size_t> _ind;
        std::vector<VertexR2> _ver;

        std::vector<VertexR2> _s;

        std::vector<double> _r;
        std::vector<double> _A;
        std::vector<double> _D;
        std::vector<double> _t;
        std::vector<double> _w;

        std::vector<double> _w0;
        std::vector<double> _w1;
        std::vector<double> _w2;

        // Compute blended coordinates for the case with one quadrilateral.
        void caseOneQuadrilateral(const VertexR2 &p,
                                  const std::vector<size_t> &nn,
                                  std::vector<double> &b) {

            // Get vertices of the quadrilateral.
            const std::vector<VertexR2> &tp = _mesh.vertices();
            const std::vector<Halfedge> &th = _mesh.halfedges();

            const size_t h0 = nn[0];

            // Build a quadrilateral.
            const Halfedge &he0 = th[h0];

            _ver[0] = tp[he0.dest];
            _ind[0] = he0.dest;

            const Halfedge &he1 = th[he0.next];

            _ver[1] = tp[he1.dest];
            _ind[1] = he1.dest;

            const Halfedge &he2 = th[he1.next];

            _ver[2] = tp[he2.dest];
            _ind[2] = he2.dest;
         
            const Halfedge &hea = th[he0.neigh];
            const Halfedge &he3 = th[hea.next];

            _ver[3] = tp[he3.dest];
            _ind[3] = he3.dest;

            // Compute blended coordinates with respect to the quadrilateral above.
            for (size_t i = 0; i < 4; ++i) {
                _s[i] = _ver[i] - p;
                _r[i] = _s[i].length();
            }

            for (size_t i = 0; i < 4; ++i) {
                const size_t ip = (i + 1) % 4;

                _A[i] = _s[i].crossProduct(_s[ip]);
                _D[i] = _s[i].scalarProduct(_s[ip]);

                assert(fabs(_r[i] * _r[ip] + _D[i]) > 0.0);
                _t[i] = _A[i] / (_r[i] * _r[ip] + _D[i]);
            }

            double W = 0.0;

            for (size_t i = 0; i < 4; ++i) {

                const size_t im = (i + 3) % 4;
                assert(fabs(_r[i]) > 0.0);

                _w[i] = (_t[im] + _t[i]) / _r[i];
                W += _w[i];
            }

            assert(fabs(W) > 0.0);
            const double invW = 1.0 / W;

            for (size_t i = 0; i < 4; ++i) b[_ind[i]] = _w[i] * invW;
        }

        // Compute blended coordinates for the case with two overlapping quadrilaterals.
        void caseTwoQuadrilaterals(const VertexR2 &p,
                                   const std::vector<size_t> &nn,
                                   std::vector<double> &b) {

            // Get vertices of two overlapping quadrilaterals that is
            // 5 vertices of the pentagon.
            const std::vector<VertexR2> &tp = _mesh.vertices();
            const std::vector<Halfedge> &th = _mesh.halfedges();

            size_t h0 = nn[0];
            size_t h1 = nn[1];

            if (th[th[h0].next].neigh != -1) {
                const size_t tmp = h1;
                h1 = h0;
                h0 = tmp;
            }

            // Build a pentagon.
            const Halfedge &he0 = th[h0];

            _ver[0] = tp[he0.dest];
            _ind[0] = he0.dest;

            const Halfedge &he1 = th[he0.next];

            _ver[1] = tp[he1.dest];
            _ind[1] = he1.dest;

            const Halfedge &he2 = th[he1.next];

            _ver[2] = tp[he2.dest];
            _ind[2] = he2.dest;

            const Halfedge &hea = th[he0.neigh];
            const Halfedge &he3 = th[hea.next];

            _ver[3] = tp[he3.dest];
            _ind[3] = he3.dest;

            const Halfedge &heb = th[th[h1].neigh];
            const Halfedge &he4 = th[heb.next];

            _ver[4] = tp[he4.dest];
            _ind[4] = he4.dest;

            // Compute blended coordinates with respect to the overlapping quadrilaterals above.
            for (size_t i = 0; i < 5; ++i) {
                _s[i] = _ver[i] - p;
                _r[i] = _s[i].length();
            }

            double Q = 0.0;

            for (size_t i = 0; i < 3; ++i) {
                const size_t ip = (i + 1) % 3;

                _A[i] = _s[i].crossProduct(_s[ip]);
                _D[i] = _s[i].scalarProduct(_s[ip]);

                Q += _A[i];

                assert(fabs(_r[i] * _r[ip] + _D[i]) > 0.0);
                _t[i] = _A[i] / (_r[i] * _r[ip] + _D[i]);
            }

            assert(fabs(Q) > 0.0);
            double invQ = 1.0 / Q;

            const double q_0 = _A[2] * invQ;
            const double q_1 = _A[1] * invQ;

            const double q0 = smoothStep(q_0);
            const double q1 = smoothStep(q_1);

            Q = q0 + q1;

            assert(fabs(Q) > 0.0);
            invQ = 1.0 / Q;

            const double w_0 = q1 * invQ;
            const double w_1 = q0 * invQ;

            _A[3] = _s[2].crossProduct(_s[3]);
            _A[4] = _s[3].crossProduct(_s[0]);
            _A[5] = _s[1].crossProduct(_s[4]);
            _A[6] = _s[4].crossProduct(_s[2]);

            _D[3] = _s[2].scalarProduct(_s[3]);
            _D[4] = _s[3].scalarProduct(_s[0]);
            _D[5] = _s[1].scalarProduct(_s[4]);
            _D[6] = _s[4].scalarProduct(_s[2]);

            assert(fabs(_r[2] * _r[3] + _D[3]) > 0.0);
            _t[3] = _A[3] / (_r[2] * _r[3] + _D[3]);
            assert(fabs(_r[3] * _r[0] + _D[4]) > 0.0);
            _t[4] = _A[4] / (_r[3] * _r[0] + _D[4]);
            assert(fabs(_r[1] * _r[4] + _D[5]) > 0.0);
            _t[5] = _A[5] / (_r[1] * _r[4] + _D[5]);
            assert(fabs(_r[4] * _r[2] + _D[6]) > 0.0);
            _t[6] = _A[6] / (_r[4] * _r[2] + _D[6]);

            assert(fabs(_r[0]) > 0.0);
            _w0[0] = (_t[4] + _t[0]) / _r[0];
            assert(fabs(_r[1]) > 0.0);
            _w0[1] = (_t[0] + _t[1]) / _r[1];
            assert(fabs(_r[2]) > 0.0);
            _w0[2] = (_t[1] + _t[3]) / _r[2];
            assert(fabs(_r[3]) > 0.0);
            _w0[3] = (_t[3] + _t[4]) / _r[3];

            double W0 = 0.0;
            W0 += _w0[0]; W0 += _w0[1]; W0 += _w0[2]; W0 += _w0[3];

            assert(fabs(W0) > 0.0);
            const double invW0 = 1.0 / W0;

            _w1[0] = (_t[2] + _t[0]) / _r[0];
            _w1[1] = (_t[0] + _t[5]) / _r[1];
            _w1[2] = (_t[6] + _t[2]) / _r[2];
            _w1[3] = (_t[5] + _t[6]) / _r[4];

            double W1 = 0.0;
            W1 += _w1[0]; W1 += _w1[1]; W1 += _w1[2]; W1 += _w1[3];

            assert(fabs(W1) > 0.0);
            const double invW1 = 1.0 / W1;

            for (size_t i = 0; i < 3; ++i)
                b[_ind[i]] = _w0[i] * invW0 * w_0 + _w1[i] * invW1 * w_1;

            b[_ind[3]] = _w0[3] * invW0 * w_0;
            b[_ind[4]] = _w1[3] * invW1 * w_1;
        }

        // Compute blended coordinates for the case with three overlapping quadrilaterals.
        void caseThreeQuadrilaterals(const VertexR2 &p,
                                     const std::vector<size_t> &nn,
                                     std::vector<double> &b) {

            // Get vertices of three overlapping quadrilaterals that is
            // 6 vertices of the hexagon.
            const std::vector<VertexR2> &tp = _mesh.vertices();
            const std::vector<Halfedge> &th = _mesh.halfedges();

            size_t h0 = nn[0];
            size_t h1 = nn[1];
            size_t h2 = nn[2];

            // Build a hexagon.
            const Halfedge &he0 = th[h0];

            _ver[0] = tp[he0.dest];
            _ind[0] = he0.dest;

            const Halfedge &he1 = th[he0.next];

            _ver[1] = tp[he1.dest];
            _ind[1] = he1.dest;

            const Halfedge &he2 = th[he1.next];

            _ver[2] = tp[he2.dest];
            _ind[2] = he2.dest;

            const Halfedge &hea = th[he0.neigh];
            const Halfedge &he3 = th[hea.next];

            _ver[3] = tp[he3.dest];
            _ind[3] = he3.dest;

            const Halfedge &heb = th[th[h1].neigh];
            const Halfedge &he4 = th[heb.next];

            _ver[4] = tp[he4.dest];
            _ind[4] = he4.dest;

            const Halfedge &hec = th[th[h2].neigh];
            const Halfedge &he5 = th[hec.next];

            _ver[5] = tp[he5.dest];
            _ind[5] = he5.dest;

            // Compute blended coordinates with respect to the overlapping quadrilaterals above.
            for (size_t i = 0; i < 6; ++i) {
                _s[i] = _ver[i] - p;
                _r[i] = _s[i].length();
            }

            double Q = 0.0;

            for (size_t i = 0; i < 3; ++i) {
                const size_t ip = (i + 1) % 3;

                _A[i] = _s[i].crossProduct(_s[ip]);
                _D[i] = _s[i].scalarProduct(_s[ip]);

                Q += _A[i];

                assert(fabs(_r[i] * _r[ip] + _D[i]) > 0.0);
                _t[i] = _A[i] / (_r[i] * _r[ip] + _D[i]);
            }

            assert(fabs(Q) > 0.0);
            double invQ = 1.0 / Q;

            const double q_0 = _A[2] * invQ;
            const double q_1 = _A[0] * invQ;
            const double q_2 = _A[1] * invQ;

            const double q0 = smoothStep(q_0);
            const double q1 = smoothStep(q_1);
            const double q2 = smoothStep(q_2);

            const double sigma_0 = q1 * q2;
            const double sigma_1 = q2 * q0;
            const double sigma_2 = q0 * q1;

            Q = sigma_0 + sigma_1 + sigma_2;

            assert(fabs(Q) > 0.0);
            invQ = 1.0 / Q;

            const double w_0 = sigma_0 * invQ;
            const double w_1 = sigma_1 * invQ;
            const double w_2 = sigma_2 * invQ;

            _A[3] = _s[2].crossProduct(_s[3]);
            _A[4] = _s[3].crossProduct(_s[0]);
            _A[5] = _s[0].crossProduct(_s[4]);
            _A[6] = _s[4].crossProduct(_s[1]);
            _A[7] = _s[1].crossProduct(_s[5]);
            _A[8] = _s[5].crossProduct(_s[2]);

            _D[3] = _s[2].scalarProduct(_s[3]);
            _D[4] = _s[3].scalarProduct(_s[0]);
            _D[5] = _s[0].scalarProduct(_s[4]);
            _D[6] = _s[4].scalarProduct(_s[1]);
            _D[7] = _s[1].scalarProduct(_s[5]);
            _D[8] = _s[5].scalarProduct(_s[2]);

            assert(fabs(_r[2] * _r[3] + _D[3]) > 0.0);
            _t[3] = _A[3] / (_r[2] * _r[3] + _D[3]);
            assert(fabs(_r[3] * _r[0] + _D[4]) > 0.0);
            _t[4] = _A[4] / (_r[3] * _r[0] + _D[4]);
            assert(fabs(_r[0] * _r[4] + _D[5]) > 0.0);
            _t[5] = _A[5] / (_r[0] * _r[4] + _D[5]);
            assert(fabs(_r[4] * _r[1] + _D[6]) > 0.0);
            _t[6] = _A[6] / (_r[4] * _r[1] + _D[6]);
            assert(fabs(_r[1] * _r[5] + _D[7]) > 0.0);
            _t[7] = _A[7] / (_r[1] * _r[5] + _D[7]);
            assert(fabs(_r[5] * _r[2] + _D[8]) > 0.0);
            _t[8] = _A[8] / (_r[5] * _r[2] + _D[8]);

            assert(fabs(_r[0]) > 0.0);
            _w0[0] = (_t[4] + _t[0]) / _r[0];
            assert(fabs(_r[1]) > 0.0);
            _w0[1] = (_t[0] + _t[1]) / _r[1];
            assert(fabs(_r[2]) > 0.0);
            _w0[2] = (_t[1] + _t[3]) / _r[2];
            assert(fabs(_r[3]) > 0.0);
            _w0[3] = (_t[3] + _t[4]) / _r[3];

            double W0 = 0.0;
            W0 += _w0[0]; W0 += _w0[1]; W0 += _w0[2]; W0 += _w0[3];

            assert(fabs(W0) > 0.0);
            const double invW0 = 1.0 / W0;

            assert(fabs(_r[0]) > 0.0);
            _w1[0] = (_t[2] + _t[5]) / _r[0];
            assert(fabs(_r[1]) > 0.0);
            _w1[1] = (_t[6] + _t[1]) / _r[1];
            assert(fabs(_r[2]) > 0.0);
            _w1[2] = (_t[1] + _t[2]) / _r[2];
            assert(fabs(_r[4]) > 0.0);
            _w1[3] = (_t[5] + _t[6]) / _r[4];

            double W1 = 0.0;
            W1 += _w1[0]; W1 += _w1[1]; W1 += _w1[2]; W1 += _w1[3];

            assert(fabs(W1) > 0.0);
            const double invW1 = 1.0 / W1;

            assert(fabs(_r[0]) > 0.0);
            _w2[0] = (_t[2] + _t[0]) / _r[0];
            assert(fabs(_r[1]) > 0.0);
            _w2[1] = (_t[0] + _t[7]) / _r[1];
            assert(fabs(_r[2]) > 0.0);
            _w2[2] = (_t[8] + _t[2]) / _r[2];
            assert(fabs(_r[5]) > 0.0);
            _w2[3] = (_t[7] + _t[8]) / _r[5];

            double W2 = 0.0;
            W2 += _w2[0]; W2 += _w2[1]; W2 += _w2[2]; W2 += _w2[3];

            assert(fabs(W2) > 0.0);
            const double invW2 = 1.0 / W2;

            for (size_t i = 0; i < 3; ++i)
                b[_ind[i]] = _w0[i] * invW0 * w_0 + _w1[i] * invW1 * w_1 + _w2[i] * invW2 * w_2;

            b[_ind[3]] = _w0[3] * invW0 * w_0;
            b[_ind[4]] = _w1[3] * invW1 * w_1;
            b[_ind[5]] = _w2[3] * invW2 * w_2;
        }

        // Compute C^1 or C^2 smooth step.
        double smoothStep(const double x) const {
            switch(_continuity) {
                case 1:
                    return (3.0 - 2.0 * x) * x * x;
                    break;
                case 2:
                    return ((6.0 * x - 15.0) * x + 10.0) * x * x * x;
                    break; 
                default:
                    return (3.0 - 2.0 * x) * x * x;
                    break;
            }
        }

        // Compute the signed area of a triangle.
        inline double triangleArea(const VertexR2 &v0, const VertexR2 &v1, const VertexR2 &v2) const {
            return 0.5 * (v1 - v0).crossProduct(v2 - v0);
        }
    };

} // namespace gbc

#endif // GBC_BLENDEDR2_HPP
