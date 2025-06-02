#include <chrono>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/integer/mod_inverse.hpp>
#include <iostream>
#include <tuple>
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

class Point;
class P_Curve;
class P_Point;

/*Affine coordinates*/
/*Elliptic curve y^2 = x^3 + ax + b (mod p)*/
/*Projective coordinates Y^2*Z = X^3 + a*X*Z^2 + b*Z^3 */

class P_Curve {
public: bigint a, b, p;

      P_Curve(const bigint& a_, const bigint& b_, const bigint& p_)
          : a(a_), b(b_), p(p_) {
      };

      bool is_point(const bigint& x, const bigint& y, const bigint& z) const {
          return (y * y * z) % p == ((x * x * x) % p + (a * x * z * z) % p + (b * z * z * z) % p) % p;
      }

      bool is_point(const bigint& x, const bigint& y) const {
          return (y * y) % p == ((x * x * x) % p + (a * x + b) % p) % p;
      }
};

class Point {
private: bigint x, y = 0;
public:

    const P_Curve& curve;

    Point(const bigint& x_, const bigint& y_, const P_Curve& curve_)
        : x(x_), y(y_), curve(curve_) {
        if (!curve.is_point(x, y)) {
            cout << "The point is not on the curve" << endl;
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& P)
    {
        os << "(" << std::hex << P.x << ", " << std::hex << P.y << ")";
        return os;
    }

    friend P_Point toProjective(Point P, bigint& z, P_Curve curve);

    bigint getx() const {
        bigint x1;
        x1 = this->x;

        return x1;
    }

    bool isPointOnCurveAff() const {
        if (curve.is_point(x, y)) {
            cout << "The Point is on the curve" << endl;
            return true;
        }
        else {
            cout << "The Point is not on the curve" << endl;
            return false;
        }
    }

    Point& operator=(const Point& P) {
        const P_Curve& curve = P.curve;

        if (this != &P) {

            x = P.x;
            y = P.y;
        }
        return *this;
    }

    //Point FindPoint(const bigint& a, const bigint& b, const bigint& p, const bigint& n);
};

class P_Point {
private:
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
            cout << "The point is not on the curve" << endl;
        }
        if (x == 0 and y == 1 and z == 0) {
            cout << "Point at Infinity" << endl;
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

    bool isPointOnCurvePr() const {
        if (curve.is_point(x, y, z)) {
            cout << "The Point is on the curve" << endl;
            return true;
        }
        else {
            cout << "The Point is not on the curve" << endl;
            return false;
        }
    }

    bigint getx() const {
        bigint x1;
        x1 = this->x;

        return x1;
    }

    bigint gety() const {
        bigint y1;
        y1 = this->y;

        return y1;
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

    friend Point toAffine(P_Point P, P_Curve curve);
};

Point toAffine(P_Point P, P_Curve curve) {
    bigint z_inv = boost::integer::mod_inverse(P.getz(), curve.p);
    bigint x, y;
    x = (P.getx() * z_inv) % curve.p;
    y = (P.gety() * z_inv) % curve.p;

    if (x < 0) x += curve.p;
    if (y < 0) y += curve.p;

    Point Paff(x, y, curve);
    return Paff;
}

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

    return G.ScalarMultMontgomeri(k);
}

/*lab 2 starts here!!!*/
bigint RandomInt(const bigint& n)
{
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
    bigint rand;

    for (size_t i = 0; i < sizeof(buffer); ++i) {
        rand = (rand << 8) | buffer[i];
    }

    rand = rand % n;
    rand += 1;
    return rand;
}

P_Point FindPoint(const P_Curve& curve, const bigint& n) {
    bigint x, y;
    bigint z("0x1");
    bigint y_sq;
    bigint test = 0;
    bigint t1 = 2;
    bigint t2 = 4;

    bigint pow1 = (curve.p - 1) * boost::integer::mod_inverse(t1, curve.p) % curve.p;
    bigint pow2 = (curve.p + 1) * boost::integer::mod_inverse(t2, curve.p) % curve.p;

    while (test != 1) {
        x = RandomInt(curve.p);
        y_sq = (x * x * x + curve.a * x + curve.b) % curve.p;
        test = powm(y_sq, pow1, curve.p) % curve.p;
        cout << "test=" << hex << test << endl;
    }

    y = powm(y_sq, pow2, curve.p);
    P_Point P(x, y, z, curve);
    return P;
}

vector<unsigned char> bigintToBytes(const bigint& num) {
    std::vector<unsigned char> bytes(32);
    bigint temp = num;
    for (int i = 31; i >= 0; --i) {
        bytes[i] = static_cast<unsigned char>(temp & 0xFF);
        temp >>= 8;
    }
    return bytes;
}

// SHA-256 using BCrypt
bigint computeSHA256(const vector<unsigned char>& input) {
    BCRYPT_ALG_HANDLE hAlg = nullptr;
    BCRYPT_HASH_HANDLE hHash = nullptr;
    NTSTATUS status;
    std::vector<unsigned char> hash(32);

    status = BCryptOpenAlgorithmProvider(&hAlg, BCRYPT_SHA256_ALGORITHM, nullptr, 0);
    if (!BCRYPT_SUCCESS(status)) {
        throw std::runtime_error("Cannot open SHA-256 provider");
    }

    status = BCryptCreateHash(hAlg, &hHash, nullptr, 0, nullptr, 0, 0);
    if (!BCRYPT_SUCCESS(status)) {
        BCryptCloseAlgorithmProvider(hAlg, 0);
        throw std::runtime_error("Cannot create hash object");
    }

    status = BCryptHashData(hHash, (PUCHAR)input.data(), (ULONG)input.size(), 0);
    if (!BCRYPT_SUCCESS(status)) {
        BCryptDestroyHash(hHash);
        BCryptCloseAlgorithmProvider(hAlg, 0);
        throw std::runtime_error("Cannot hash data");
    }

    status = BCryptFinishHash(hHash, hash.data(), (ULONG)hash.size(), 0);
    if (!BCRYPT_SUCCESS(status)) {
        BCryptDestroyHash(hHash);
        BCryptCloseAlgorithmProvider(hAlg, 0);
        throw std::runtime_error("Cannot finalize hash");
    }

    BCryptDestroyHash(hHash);
    BCryptCloseAlgorithmProvider(hAlg, 0);

    bigint result = 0;
    for (unsigned char b : hash) {
        result = (result << 8) | b;
    }

    return result;
}

class ECDH {
private:
    bigint dA;
    bigint dB;
public:
    tuple<P_Point, P_Point, P_Point> KeyExchange(P_Point& G, const P_Curve& curve, const bigint& n);
    tuple<P_Point, bigint, bigint> encryption(const bigint& M, P_Point P, P_Point QB, const P_Curve& curve, const bigint& n);
    bigint Decryption(P_Point QA, bigint Ck, bigint Cm, const P_Curve& curve, const bigint& n);
    tuple<bigint, bigint> Signature(const vector<unsigned char>& M, P_Point P, const P_Curve& curve, const bigint& n);
    void CheckSignature(const bigint& r, const bigint& s, P_Point QA, const vector<unsigned char>& M, P_Point P, const P_Curve& curve, const bigint& n);
};


tuple<P_Point, P_Point, P_Point> ECDH::KeyExchange(P_Point& G, const P_Curve& curve, const bigint& n) {
    P_Point P = FindPoint(curve, n);

    dB = RandomInt(n);
    P_Point QB = P.ScalarMultMontgomeri(dB);//QB - Bob's public key

    dA = RandomInt(n);
    P_Point QA = P.ScalarMultMontgomeri(dA);//QA - Alice's public key

    P_Point keyA = QB.ScalarMultMontgomeri(dA);//Alice creates key
    P_Point keyB = QA.ScalarMultMontgomeri(dB);//Bob creates key

    Point aff_keyA = toAffine(keyA, curve);
    Point aff_keyB = toAffine(keyB, curve);

    cout << "Alice's key = " << aff_keyA << endl;
    cout << "Bob's key = " << aff_keyB << endl;

    bigint keyA1 = aff_keyA.getx() % curve.p;
    bigint keyB1 = aff_keyB.getx() % curve.p;

    if (keyA1 == keyB1) {
        cout << "Keys are equal!" << endl;
    }
    else {
        cout << "Keys are NOT equal!" << endl;
    }

    return { P, QA, QB };
}

tuple<P_Point, bigint, bigint> ECDH::encryption(const bigint& M, P_Point P, P_Point QB, const P_Curve& curve, const bigint& n) {

    //Генерація k
    bigint k = RandomInt(n);

    //Зашифрування
    bigint Cm = M ^ k;//xor

    /*Ефемерна ключова пара Аліси*/
    bigint eA = RandomInt(n);
    P_Point QA = P.ScalarMultMontgomeri(eA);

    /*Секретний ключ*/
    P_Point S = QB.ScalarMultMontgomeri(eA);
    Point aff_S = toAffine(S, curve);
    bigint Sx = aff_S.getx();

    //Інкапсуляція ключа
    bigint Ck = k ^ Sx;

    return { QA, Ck, Cm };
}

bigint ECDH::Decryption(P_Point QA, bigint Ck, bigint Cm, const P_Curve& curve, const bigint& n) {
    /*Секретний ключ*/
    P_Point S = QA.ScalarMultMontgomeri(dB);
    Point aff_S = toAffine(S, curve);
    bigint Sx = aff_S.getx();

    //Декапсуляція ключа
    bigint k = (Ck ^ Sx) % curve.p;
    cout << "Key = " << k << endl;

    //Розшифровка
    bigint M = Cm ^ k;
    cout << "Deciphered M = " << M << endl;

    return M;
}

tuple<bigint, bigint> ECDH::Signature(const vector<unsigned char>& M, P_Point P, const P_Curve& curve, const bigint& n) {

    bigint h = computeSHA256(M);
    //Генерація k
    bigint r = 0;
    bigint s;

    while (r == 0) {
        bigint k = RandomInt(n);

        P_Point kP = P.ScalarMultMontgomeri(k);
        Point aff_kP = toAffine(kP, curve);
        bigint Sx = aff_kP.getx();

        r = Sx % n;

        s = boost::integer::mod_inverse(k, n) * (h + dA * r) % n;
    }

    return { r, s };
}

void ECDH::CheckSignature(const bigint& r, const bigint& s, P_Point QA, const vector<unsigned char>& M, P_Point P, const P_Curve& curve, const bigint& n) {
    bigint h = computeSHA256(M);

    bigint u1 = boost::integer::mod_inverse(s, n) * h % n;
    bigint u2 = boost::integer::mod_inverse(s, n) * r % n;

    P_Point R = P.ScalarMultMontgomeri(u1) + QA.ScalarMultMontgomeri(u2);
    Point aff_R = toAffine(R, curve);
    bigint Sx = aff_R.getx();

    bigint v = Sx % n;

    if (v == r) {
        cout << "Signature is correct!" << endl;
    }
    else {
        cout << "Signature is not correct!" << endl;
    }
}

int main() {
    /*Curve P-256*/
    bigint a("0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc");
    bigint b("0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");
    bigint p("0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff");

    bigint n("0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551");

    /*G - generator*/
    bigint Gx("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296");
    bigint Gy("0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
    bigint Gz("0x1");

    P_Curve curve(a, b, p);
    P_Point G(Gx, Gy, Gz, curve);
    P_Point O(curve);

    P_Point test_random = FindPoint(curve, n);
    cout << test_random << endl;
    test_random.isPointOnCurvePr();


    P_Point test_random2 = test_random.ScalarMultMontgomeri(n);
    cout << test_random2 << endl;

    bigint M("0x4fe34f9bf9e162bce37c056406376b315ececbb8ee7eb4a82e2fe1a737bf51f5");
    ECDH S;

    auto [P, QA, QB] = S.KeyExchange(G, curve, n);
    auto [QA1, Ck, Cm] = S.encryption(M, P, QB, curve, n);
    bigint M_decyphered = S.Decryption(QA1, Ck, Cm, curve, n);

    if (M == M_decyphered) {
        cout << "Message M and decyphered M are equal" << endl;
    }
    else {
        cout << "Messages are not equal" << endl;
    }

    vector<unsigned char> M1 = { 0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF };
    auto [r, s] = S.Signature(M1, P, curve, n);
    S.CheckSignature(r, s, QA, M1, P, curve, n);

    cin.get();
    return 0;
}
