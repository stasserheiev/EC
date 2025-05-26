#include <chrono>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/integer/mod_inverse.hpp>
#include <iostream>
#include <windows.h>
#include <bcrypt.h>
using namespace std;
using namespace boost::multiprecision;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

#pragma comment(lib, "bcrypt.lib")

using bigint = cpp_int;

class Curve;
class Point;
class P_Curve;
class P_Point;

/*Affine coordinates*/
/*Elliptic curve y^2 = x^3 + ax + b (mod p)*/
class Curve {
public:    bigint a, b, p;

      Curve(const bigint& a_, const bigint& b_, const bigint& p_)
          : a(a_), b(b_), p(p_) {
      };

public: bool is_point(const bigint& x, const bigint& y) const {
    return (y * y) % p == (x * x * x + a * x + b) % p;
}
};

class Point {
private: bigint x, y;
public:

    const Curve& curve;

    Point(const bigint& x_, const bigint& y_, const Curve& curve_)
        : x(x_), y(y_), curve(curve_) {
        if (!curve.is_point(x, y)) {
            cout << "The point is not on the curve";
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& P)
    {
        os << "(" << std::hex << P.x << ", " << std::hex << P.y << ")";
        return os;
    }

    friend P_Point toProjective(Point P, bigint& z, P_Curve curve);
};


/*Projective coordinates Y^2*Z = X^3 + a*X*Z^2 + b*Z^3 */
class P_Curve {
public: bigint a, b, p;

      P_Curve(const bigint& a_, const bigint& b_, const bigint& p_)
          : a(a_), b(b_), p(p_) {
      };

      bool is_point(const bigint& x, const bigint& y, const bigint& z) const {
          return (y * y * z) % p == (x * x * x + a * x * z * z + b * z * z * z) % p;
      }
};



class P_Point {
protected:
    bigint x, y, z;
public:
    const P_Curve& curve;
    bool isinf;

public:
    P_Point Double();
    P_Point ScalarMultMontgomeri(bigint k);

    /*constructor: projective coordinates*/
    P_Point(const bigint& x_, const bigint& y_, const bigint& z_, const P_Curve& curve_)
        : x(x_), y(y_), z(z_), isinf(false), curve(curve_) {
        if (!curve.is_point(x, y, z)) {
            cout << "The point is not on the curve";
        }
        if (x == 0 and y == 1 and z == 0) {
            cout << "Point at Infinity";
            isinf = true;
        }
        else {
            isinf = false;
        }
    }

    P_Point(const P_Curve& curve_)
        : x(0), y(1), z(0), curve(curve_) {
    }

    P_Point(const bigint& x_, const bigint& y_, const P_Curve& curve_)
        : x(x_), y(y_), curve(curve_) {
    }

    bool isPointAtInfinity() const {
        return isinf == true || z == 0;
    }

    bool isPointOnCurve() const {
        if (curve.is_point(x, y, z)) {
            cout << "The Point is on the curve" << endl;
            return true;
        }
        else {
            cout << "The Point is not on the curve" << endl;
            return false;
        }
    }

    bigint getz() const {
        bigint z1;
        z1 = this->z;

        return z1;
    }

    P_Point& operator=(const P_Point& P) {
        const P_Curve& curve = P.curve;

        if (this != &P) {

            x = P.x;
            y = P.y;
            z = P.z;
        }
        return *this;
    }

    /*Add of Points*/
    friend    P_Point operator+(P_Point P, P_Point Q) {
        bigint u1, u2, v1, v2, u, v, w, a;
        bigint Sx, Sy, Sz;

        const P_Curve& curve = P.curve;

        if (P.isPointAtInfinity()) {
            return Q;
        }
        else if (Q.isPointAtInfinity()) {
            return P;
        }
        else {
            u1 = Q.y * P.z;
            u2 = P.y * Q.z;
            v1 = Q.x * P.z;
            v2 = P.x * Q.z;

            if (v1 == v2)
            {
                if (u1 != u2)
                {
                    return P_Point(curve);
                }

                if (u1 = u2)
                {
                    return P.Double();
                }
            }
            else
            {
                u = u1 - u2;
                v = v1 - v2;
                w = P.z * Q.z;
                a = u * u * w - v * v * v - 2 * v * v * v2;

                Sx = (v * a) % curve.p;
                Sy = (u * (v * v * v2 - a) - v * v * v * u2) % curve.p;
                Sz = (v * v * v * w) % curve.p;

                if (Sx < 0) Sx += curve.p;
                if (Sy < 0) Sy += curve.p;
                if (Sz < 0) Sz += curve.p;

                return P_Point(Sx, Sy, Sz, curve);
            }
        }

    };

    friend std::ostream& operator<<(std::ostream& os, const P_Point& P)
    {
        os << "(" << std::hex << P.x << ", " << std::hex << P.y << ", " << std::hex << P.z << ")";
        return os;
    }

    friend Point toAffine(P_Point P, Curve curve) {
        bigint x_aff;
        bigint y_aff;

        bigint z_inv = boost::integer::mod_inverse(P.z, curve.p);
        x_aff = (P.x * z_inv) % curve.p;
        y_aff = (P.y * z_inv) % curve.p;

        if (x_aff < 0) x_aff += curve.p;
        if (y_aff < 0) y_aff += curve.p;

        return Point(x_aff, y_aff, curve);
    }
};

P_Point toProjective(Point P, bigint& z, P_Curve curve) {
    bigint x_proj;
    bigint y_proj;

    x_proj = (P.x * z) % curve.p;
    y_proj = (P.y * z) % curve.p;

    if (x_proj < 0) x_proj += curve.p;
    if (y_proj < 0) y_proj += curve.p;

    return P_Point(x_proj, y_proj, z, curve);
}


/*Double Point*/
P_Point P_Point::Double() {
    bigint w, s, b1, h;
    bigint Dx, Dy, Dz;

    const P_Curve& curve = this->curve;

    if (isPointAtInfinity() || (y == 0)) {
        return P_Point(curve);
    }
    else {
        w = curve.a * z * z + 3 * x * x;
        s = y * z;
        b1 = x * y * s;
        h = w * w - 8 * b1;

        Dx = (2 * h * s) % curve.p;
        Dy = (w * (4 * b1 - h) - 8 * y * y * s * s) % curve.p;
        Dz = (8 * s * s * s) % curve.p;

        if (Dx < 0) Dx += curve.p;
        if (Dy < 0) Dy += curve.p;
        if (Dz < 0) Dz += curve.p;

        return P_Point(Dx, Dy, Dz, curve);
    }
}


P_Point P_Point::ScalarMultMontgomeri(bigint k) {
    P_Point R0(curve);
    P_Point R1 = *this;

    size_t bit_count = (k == 0) ? 1 : msb(k) + 1;

    for (size_t i = bit_count; i-- > 0;) {
        if (bit_test(k, i) == 0) {
            R1 = R0 + R1;
            R0 = R0.Double();
        }
        else {
            R0 = R0 + R1;
            R1 = R1.Double();
        }
    }
    return R0;
}

/*Generating of Random point*/
P_Point RandomPoint(P_Point& G, const P_Curve& curve, const bigint& n) {
    BCRYPT_ALG_HANDLE hAlg;
    if (BCryptOpenAlgorithmProvider(&hAlg, BCRYPT_RNG_ALGORITHM, nullptr, 0) != 0) {
        throw std::runtime_error("Cannot open RNG provider");
    }
    unsigned char buffer[32];
    if (BCryptGenRandom(hAlg, buffer, sizeof(buffer), 0) != 0) {
        BCryptCloseAlgorithmProvider(hAlg, 0);
        throw std::runtime_error("Cannot generate random bytes");
    }
    BCryptCloseAlgorithmProvider(hAlg, 0);

    bigint k = 0;
    for (size_t i = 0; i < sizeof(buffer); ++i) {
        k = (k << 8) | buffer[i];
    }

    k = k % (n - 1);
    k += 1;

    cout << "k=" << k << endl;

    return G.ScalarMultMontgomeri(k);
}

int main() {
 /*Curve P-384*/
    bigint a("0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc");
    bigint b("0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");
    bigint p("0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff");

    bigint n("0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551");

    /*G - generator*/
    bigint Gx("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296");
    bigint Gy("0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
    bigint Gz("0x1");

    Curve curve_aff(a, b, p);
    P_Curve curve(a, b, p);
    P_Point G(Gx, Gy, Gz, curve);
    P_Point O(curve);

    cout << "Point at Infinity O: " << O << endl;
    cout << "O + O = " << O + O << endl;
    cout << "Generator G: " << G << endl;

    auto t1 = high_resolution_clock::now();
    P_Point test = RandomPoint(G, curve, n);
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    cout << "Duration of RandomPoint() = " << ms_double.count() << endl;

    cout << "Random Point: " << test << endl;
    t1 = high_resolution_clock::now();
    P_Point test1 = test.ScalarMultMontgomeri(n);
    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    cout << "Duration of ScalarMultiplicationMontgomery = " << ms_double.count() << endl;

    cout << "Test (Should be Point at infinity): " << test1 << endl;

    
    P_Point test2 = test1 + test;
    
    cout << "Random Point: " << test2 << endl;

    t1 = high_resolution_clock::now();
    P_Point test3 = test + G;
    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    cout << "Duration of PointAdd = " << ms_double.count() << endl;


    cin.get();
    return 0;
}
