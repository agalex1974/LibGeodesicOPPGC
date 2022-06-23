#pragma once

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
#ifndef _AFX
#include <windows.h>
#endif
#include <tchar.h>
#include <ios>
#include <math.h>
#include <stdint.h>
#include <memory>
#include <string>
#include <stdexcept>
#include <float.h>
#include "Arrays.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// General MACROS and TYPES //////////////////////////////////////////

#ifndef _AFX
#ifndef ASSERT
#define ASSERT(expr)
#endif
#ifndef VERIFY
#define VERIFY(expr) (expr)
#endif
#endif

#define DELETE_DATA(data) {if(data) delete data; data=NULL;}
#define COPY_STRING(dest, src) {dest=_T("");if(src!=_T("")&& !src.IsEmpty())dest=src;	}
#define COPY_STRING2(dest, src) {dest=_T("");if(src!=NULL && _tcslen(src)>0)dest=src;	}

// set maximum tolerance TOL
const double MACHINE_EPSILON = sqrt(DBL_EPSILON); // machine depedent
const double TOL = 1e-6;
#define MAX_INT16 65535
const double SQR_ROOT_OF_2 = sqrt(2.0);

typedef struct _VECTOR2D {
	double x, y;
}VECTOR2D;

// set pi and pi/180
const double M_PI = 3.1415926535897932384626433832795;
const double M_PI_180 = M_PI / 180.0;

#ifndef MAX_INT_16
#define MAX_INT_16 32765
#endif

#ifndef MIN_INT_16
#define MIN_INT_16 -32765
#endif

#ifndef SQR
#define SQR(x) ((x) * (x))
#endif

#ifndef CUBE
#define CUBE(x) ((x) * (x) * (x))
#endif

#ifndef MAX_INT_32
#define MAX_INT_32 INT_MAX
#endif
#ifndef MIN_INT_32
#define MIN_INT_32 INT_MIN
#endif

#define MAX_DOUBLE DBL_MAX;
#define MIN_DOUBLE DBL_MIN;
#define MAX_FLOAT FLT_MAX;
#define MIN_FLOAT FLT_MIN;

#define INCH2MM(inch) ( ( (inch) * 25.4) )

#define LOG_BASE_2(x) (log(x)/log(2.0))

#define LPCSTR LPCWSTR

/*inline int MsgBoxError(const std::string& str)
{
	return ::MessageBox(0, str.c_str(), ("Error", MB_OK | MB_ICONERROR);
}

inline int MsgBoxQuestion(const std::string& str)
{
	return ::MessageBox(0, str.c_str(), "Question", MB_OKCANCEL | MB_ICONQUESTION);
}

inline int MsgBoxInfo(const std::string& str)
{
	return ::MessageBox(0, str.c_str(), "Info", MB_OK | MB_ICONEXCLAMATION);
}*/

template<typename ... Args>
std::string string_format(const std::string& format, Args ... args)
{
	int size_s = snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
	if (size_s <= 0) { throw std::runtime_error("Error during formatting."); }
	auto size = static_cast<size_t>(size_s);
	auto buf = make_unique<char[]>(size);
	snprintf(buf.get(), size, format.c_str(), args ...);
	return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

// Closest Power of 2
inline long ClosestPower2(long x)
{
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	x++;
	return x;
}

#define __LOG2A(s) ((s &0xffffffff00000000) ? (32 +__LOG2B((s) >>32)): (__LOG2B(s)))
#define __LOG2B(s) ((s &0xffff0000)         ? (16 +__LOG2C((s) >>16)): (__LOG2C(s)))
#define __LOG2C(s) ((s &0xff00)             ? (8  +__LOG2D((s) >> 8)) : (__LOG2D(s)))
#define __LOG2D(s) ((s &0xf0)               ? (4  +__LOG2E((s) >> 4)) : (__LOG2E(s)))
#define __LOG2E(s) ((s &0xc)                ? (2  +__LOG2F((s) >> 2)) : (__LOG2F(s)))
#define __LOG2F(s) ((s &0x2)                ? (1)                  : (0))

#define LOG2_UINT64 __LOG2A
#define LOG2_UINT32 __LOG2B
#define LOG2_UINT16 __LOG2C
#define LOG2_UINT8  __LOG2D

// Next Power of 2
inline uint64_t next_power_of_2(uint64_t i)
{
#if defined(__GNUC__)
	return 1UL << (64 - __builtin_clzl(i - 1));
#else
	return 1ULL << (1i64 + LOG2_UINT64((long long)i - 1i64));
#endif
}

inline double GetRandomNum(double a, double b)
{
	return a + (b - a) * (rand() / double(RAND_MAX));
}

inline int GetRandomNum(int a, int b)
{
	return a + rand() % (b + 1 - a);
}

// We need to define when a scalar should be consider equal to
// zero, or equal to unit. This is achieved through the ZERO_INTERVAL definition.
inline double ZERO_INTERVAL(double x)
{
	if (fabs(x) <= TOL)
		return 0.0;
	else
		if (fabs(fabs(x) - 1) <= TOL)
			return 1.0;
		else
			return x;
}

inline void ZERO_INTERVAL2(double& x)
{
	if (fabs(x) <= TOL)
		x = 0.0;
	else
		if (fabs(fabs(x) - 1) <= TOL)
			x = 1.0;
}

#define ZERO_INTERVAL_POINTEX(p) \
{	if(fabs((p).x) <= TOL)\
		(p).x = 0.0;\
	else\
	if(fabs(fabs((p).x) - 1) <= TOL)\
		(p).x = 1.0;\
	if(fabs((p).y) <= TOL)\
		(p).y = 0.0;\
	else\
	if(fabs(fabs((p).y) - 1) <= TOL)\
		(p).y = 1.0;\
}

const double TOL_NUM = 1e-12;
#define ZERO_INTERVAL_POINTEX3D(p) \
{	if(fabs((p).x) <= TOL_NUM)\
		(p).x = 0.0;\
	else\
	if(fabs(fabs((p).x) - 1) <= TOL_NUM)\
		(p).x = 1.0;\
	if(fabs((p).y) <= TOL_NUM)\
		(p).y = 0.0;\
	else\
	if(fabs(fabs((p).y) - 1) <= TOL_NUM)\
		(p).y = 1.0;\
	if(fabs((p).z) <= TOL_NUM)\
		(p).z = 0.0;\
	else\
	if(fabs(fabs((p).z) - 1) <= TOL_NUM)\
		(p).z = 1.0;\
}

inline double ZERO_INTERVAL_0(double x)
{
	if (fabs(x) <= TOL)
		return 0.0;
	else
		return x;
}
inline double ZERO_INTERVAL_1(double x)
{
	if (fabs(fabs(x) - 1) <= TOL)
		return 1.0;
	else
		return x;
}

// check if p lies within [0,1]
#define IS_VALID_INTERVAL(p)((p).x>=0 && (p).x<=1 && (p).y>=0 && (p).y<=1)

inline double sind(double theta)
{
	return sin(theta * M_PI_180);
}

inline double cosd(double theta)
{
	return cos(theta * M_PI_180);
}

inline double InRadians(double degrees)
{
	return (degrees* M_PI_180);
}

inline double InDegrees(double radians)
{
	return (radians / M_PI_180);
}

// sinus and cosine of theta in degrees
#define SIND(X) sin((X) * M_PI_180)
#define COSD(X) cos((X) * M_PI_180)
#define RADIANS(X) ((X) * M_PI_180) // return the radians corresponding to X degrees
#define DEGREES(X) ((X) / M_PI_180) // return the degrees corresponding to X radians
#define IsDirSep(X) ((X) == _T('\\') || (X) == _T('/') )


// *****************************************************************
// ** CPointEx
// ** ~~~~~~~
// ** This is an extension of the MFC class CPoint.
// ** The x, y members are declared as double
// ** It provides i/o serialization
// *****************************************************************
//////////////////////////////////////////////////////////////////////
// MACROS of CPointEx ////////////////////////////////////////////////
#define NORM_EX(X) sqrt(SQR((X).x) + SQR((X).y))
#define SQR_NORM_EX(X) (SQR((X).x) + SQR((X).y)) // The square of norm A
#define NORMALIZE_EX(X) { double d = NORM_EX(X); ASSERT(d != 0);	\
	(X).x /= d; (X).y /= d; }
#define DISTANCE_EX(A,B) sqrt(SQR((A).x-(B).x) + SQR((A).y-(B).y)) // The distance between A, B
#define DISTANCEf(A,B) sqrt(SQR((float)(A).x-(float)(B).x) + SQR((float)(A).y-(float)(B).y)) // The distance between A, B
#define SQR_DISTANCE_EX(A,B) (SQR((A).x-(B).x) + SQR((A).y-(B).y)) // The square distance between A, B
#define DOT_EX(A,B) (((A).x)*((B).x) + ((A).y)*((B).y)) // Returns the dot product A * B

// The cross vector of two CPointEx structures is a vector (0,0,x)
// normal to the xy plane. So, we return only its third coordinate.
#define CROSS_EX(A,B) (((A).x)*((B).y) - ((A).y)*((B).x))

// Rotate a given point (p), throught an arbitrary axis (fixed)
// and by an angle theta.
#define ROTATE_EX(A,C,F) \
	{	\
		double c_theta = COSD(F);	double s_theta = SIND(F);	\
		double x = A.x;		double y = A.y;	\
		A.x = x * c_theta - y * s_theta - C.x * c_theta + C.y * s_theta + C.x;	\
		A.y = x * s_theta + y * c_theta - C.x * s_theta - C.y * c_theta + C.y;	\
	}
// the above process is equivalent to the following transformations.
// expanding these transformations we get the above relations.
//	A -= C;
//	A.Rotate(F);
//	A += C;

// Check if 3 planar points are colinear using maximum tolerance TOL
#define COLINEAR_EX(A,B,C) \
	(fabs(((A).x-(B).x)*((C).y-(B).y) - ((C).x-(B).x)*((A).y-(B).y)) <= TOL)

// Homogeneous Scaling of point p using:
// fixed = center of scaling
// s = scale factor along x,y-axis
#define SCALE_EX(C,A,S) \
	{	\
		A.x = ((A).x - (C).x) * (S) + (C).x; \
		A.y = ((A).y - (C).y) * (S) + (C).y; \
	}

// Differential Scaling of point p using:
// fixed = center of scaling
// s = scale factor along x-axis
// t = scale factor along y-axis
#define DIF_SCALE_EX(C, A, S, T)	\
	{	\
		A.x = ((A).x - (C).x) * (S) + (C).x; \
		A.y = ((A).y - (C).y) * (T) + (C).y; \
	}

class CPointEx
{
	// Constructors
public:
	CPointEx() // default constructor
	{
		x = 0;	y = 0;
	}

	CPointEx(double xx, double yy)
	{
		x = xx;	y = yy;
	}

	CPointEx(const CPointEx& point_ex)
	{
		x = point_ex.x;	y = point_ex.y;
	}

	// data members
public:
	double x, y;

	// implementation
	inline BOOL operator == (const CPointEx& point_ex) const
	{
		return (x == point_ex.x && y == point_ex.y);
	}

	inline BOOL operator != (const CPointEx& point_ex) const
	{
		return (x != point_ex.x || y != point_ex.y);
	}

	inline void operator += (const CPointEx& point_ex)
	{
		x += point_ex.x;
		y += point_ex.y;
	}

	inline void operator -= (const CPointEx& point_ex)
	{
		x -= point_ex.x;
		y -= point_ex.y;
	}

	inline void operator = (const CPointEx& point_ex)
	{
		x = point_ex.x;
		y = point_ex.y;
	}

	inline void operator *= (double m)
	{
		x *= m;
		y *= m;
	}

	inline void operator /= (double m)
	{
		ASSERT(m != 0);
		x /= m;
		y /= m;
	}

	inline double& operator [] (int i)
	{
		ASSERT(i >= 0 && i < 2);
		if (i == 0)
			return x;
		else
			return y;
	}

	inline const double& operator [] (int i) const
	{
		ASSERT(i >= 0 && i < 2);
		if (i == 0)
			return x;
		else
			return y;
	}

	inline BOOL IsZero() const
	{
		return (x == 0 && y == 0);
	}

	inline void Rotate(double phi)
	{
		double cos_phi = COSD(phi);
		double sin_phi = SIND(phi);

		double xx;
		double yy;
		xx = x * cos_phi - y * sin_phi;
		yy = x * sin_phi + y * cos_phi;

		x = xx;
		y = yy;
	}

	inline void Move(const CPointEx& d)
	{
		x += d.x; y += d.y;
	}

	inline void Rotate(const CPointEx& fixed, double theta)
	{
		double c_theta = COSD(theta);
		double s_theta = SIND(theta);
		double xx = x;
		double yy = y;

		x = xx * c_theta - yy * s_theta - fixed.x * c_theta + fixed.y * s_theta + fixed.x;
		y = xx * s_theta + yy * c_theta - fixed.x * s_theta - fixed.y * c_theta + fixed.y;
	}

	// use MACROS instead of these functions for better performance
	inline double norm() const
	{
		return (sqrt(x*x + y * y));
	}

	inline double sqr_norm() const
	{
		return (x*x + y * y);
	}

	inline void normalize()
	{
		double d = norm();
		ASSERT(d != 0);
		if (d == 0)
			return;
		x = x / d;
		y = y / d;
	}

	inline CPointEx GetMirror(const CPointEx& p1, const CPointEx& p2) const
	{
		ASSERT(p1 != p2);

		CPointEx n1 = (CPointEx)(p2 - p1);
		NORMALIZE_EX(n1);

		CPointEx mirrorPnt(*this);
		CPointEx p_p1 = (CPointEx)(mirrorPnt - p1);
		mirrorPnt += -2 * (p_p1 - (p_p1 * n1) * n1);
		return mirrorPnt;
	}

	// (x1,y1) and (x2,y2) are the endpoints of line
	inline double DistanceFromLine(const double x1, const double y1,
		const double x2, const double y2)
	{
		double L = sqrt(SQR(x2 - x1) + SQR(y2 - y1));
		if (L <= 1e-12)
			return -1; // invalid line segment
		else
			return (fabs((x - x1)*(y2 - y1) - (y - y1)*(x2 - x1)) / L);
	}

	// p1 and p2 are the endpoints of line
	inline double DistanceFromLine(const CPointEx& p1, const CPointEx& p2)
	{
		double L = DISTANCE_EX(p1, p2);
		if (L <= 1e-12)
			return -1; // invalid line segment
		else
			return (fabs((x - p1.x)*(p2.y - p1.y) - (y - p1.y)*(p2.x - p1.x)) / L);
	}

	// operators
	friend double operator * (const CPointEx& w, const CPointEx& v);
	friend CPointEx operator * (const CPointEx& p, double m);
	friend CPointEx operator / (const CPointEx& p, double m);
	friend CPointEx operator * (double m, const CPointEx& p);
	friend CPointEx operator + (const CPointEx& p1, const CPointEx& p2);
	friend CPointEx operator - (const CPointEx& p1, const CPointEx& p2);
	friend double distance(const CPointEx& p1, const CPointEx& p2);
	friend double sqr_distance(const CPointEx& p1, const CPointEx& p2);

	static void Rotate(CPointEx* p, int nSize, double phi)
	{
		double cos_phi = COSD(phi);
		double sin_phi = SIND(phi);

		for (int i = 0; i < nSize; i++)
		{
			double xx = p[i].x * cos_phi - p[i].y * sin_phi;
			double yy = p[i].x * sin_phi + p[i].y * cos_phi;

			p[i].x = xx;
			p[i].y = yy;
		}
	}

	static void Rotate(CPointEx* p, int nSize, double theta, const CPointEx& fixed)
	{
		double c_theta = COSD(theta);
		double s_theta = SIND(theta);
		for (int i = 0; i < nSize; i++)
		{
			double xx = p[i].x;
			double yy = p[i].y;

			p[i].x = xx * c_theta - yy * s_theta - fixed.x * c_theta + fixed.y * s_theta + fixed.x;
			p[i].y = xx * s_theta + yy * c_theta - fixed.x * s_theta - fixed.y * c_theta + fixed.y;
		}
	}

	static void Move(CPointEx* p, int nSize, const CPointEx& d)
	{
		for (int i = 0; i < nSize; i++)
			p[i] += d;
	}
};


typedef CPointEx POINTEX;
typedef POINTEX* LPPOINTEX;
typedef const POINTEX* LPCTPOINTEX;
#define no_point2D CPointEx(-999999.000091, -999999.000091)
#define NULL_POINTEX CPointEx(-1e9, -1e9)

double PolygonArea(const CPointEx* poly, int nVertices);

//  calculate the length of nPoints sequential CPointEx points.
inline double CalcLength(LPCTPOINTEX p, int nPoints)
{
	double d = 0;
	for (int i = 0; i < nPoints - 1; i++)
		d += distance(p[i], p[i + 1]);
	return d;
}

// *****************************************************************
// ** CPointEx3D
// ** ~~~~~~~
// ** This is an extension of the class CPointEx to support 3D points.
// **
// *****************************************************************
class CPointEx3D
{
	// data members
public:
	double x, y, z, w;

	// Constructors
public:
	CPointEx3D() // default constructor
	{
		x = 0;	y = 0; z = 0; w = 1;
	}

	CPointEx3D(double xx, double yy, double zz, double ww = 1.0)
	{
		x = xx;	y = yy; z = zz; w = ww;
	}

	CPointEx3D(const CPointEx3D& p)
	{
		x = p.x; y = p.y; z = p.z; w = p.w;
	}

	// implementation

	inline BOOL operator == (const CPointEx3D& p) const
	{
		return (x == p.x && y == p.y && z == p.z && w == p.w);
	}

	inline BOOL operator != (const CPointEx3D& p) const
	{
		return (x != p.x || y != p.y || z != p.z || w != p.w);
	}

	inline void operator += (const CPointEx3D& p)
	{
		x += p.x;
		y += p.y;
		z += p.z;
		w = p.w;
	}

	inline void operator -= (const CPointEx3D& p)
	{
		x -= p.x;
		y -= p.y;
		z -= p.z;
		w = p.w;
	}

	inline void operator = (const CPointEx3D& p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
		w = p.w;
	}

	inline void operator *= (double m)
	{
		x *= m;
		y *= m;
		z *= m;
	}

	inline void operator /= (double m)
	{
		ASSERT(m != 0);
		x /= m;
		y /= m;
		z /= m;
	}

	inline void SetPoint(double xx, double yy, double zz)
	{
		x = xx;		y = yy;		z = zz;
	}

	inline double& operator [] (int i)
	{
		ASSERT(i >= 0 && i < 3);
		if (i == 0)
			return x;
		else
			if (i == 1)
				return y;
			else
				return z;
	}

	inline const double& operator [] (int i) const
	{
		ASSERT(i >= 0 && i < 3);
		if (i == 0)
			return x;
		else
			if (i == 1)
				return y;
			else
				return z;
	}

	inline BOOL IsZero() const
	{
		return (x == 0 && y == 0 && z == 0);
	}

	inline double norm() const
	{
		return (sqrt(x*x + y * y + z * z));
	}

	inline double sqr_norm() const
	{
		return (x*x + y * y + z * z);
	}

	inline void normalize()
	{
		double d = norm();
#ifndef _DEBUG
		if (d < 1e-12)
			return;
#endif
		ASSERT(d != 0);
		x = x / d;
		y = y / d;
		z = z / d;
		w = 1;
	}

	// operators
	friend double dot(const CPointEx3D& w, const CPointEx3D& v);
	friend double operator * (const CPointEx3D& w, const CPointEx3D& v);
	friend CPointEx3D operator * (const CPointEx3D& p, double m);
	friend CPointEx3D operator * (double m, const CPointEx3D& p);
	friend	CPointEx3D operator / (const CPointEx3D& p, double m);
	friend CPointEx3D operator + (const CPointEx3D& p1, const CPointEx3D& p2);
	friend CPointEx3D operator - (const CPointEx3D& p1, const CPointEx3D& p2);
	friend double distance(const CPointEx3D& p1, const CPointEx3D& p2);
	friend double sqr_norm(const CPointEx3D& p);
	friend CPointEx3D cross(const CPointEx3D& p1, const CPointEx3D& p2);
	friend double gross(const CPointEx3D& x, const CPointEx3D& y, const CPointEx3D& z);
	friend CPointEx3D Isometric(const CPointEx3D& p, double theta, double phi);
	friend CPointEx3D Dimetric(const CPointEx3D& p);
};

typedef CVector<CPointEx3D> CPointEx3DVector;

inline double gross(const CPointEx3D& x,
	const CPointEx3D& y,
	const CPointEx3D& z)
{
	return (x.x*(y.y*z.z - y.z*z.y) - x.y*(y.x*z.z - y.z*z.x) + x.z*(y.x*z.y - z.x*y.y));
}


typedef CPointEx3D POINTEX3D;
typedef POINTEX3D* LPPOINTEX3D;
typedef const POINTEX3D* LPCTPOINTEX3D;

#define no_point CPointEx3D(-999999.000091, -999999.000091, -999999.000091)
#define no_point2D CPointEx(-999999.000091, -999999.000091)

//  calculate the length of nPoints sequential CPointEx3D points.
inline double CalcLength(LPCTPOINTEX3D p, int nPoints)
{
	double d = 0;
	for (int i = 0; i < nPoints - 2; i++)
		d += distance(p[i], p[i + 1]);
	return d;
}

// Returns the dot product w * v
inline double operator * (const CPointEx& w, const CPointEx& v)
{
	return (w.x * v.x + w.y * v.y);
}

// Returns  p(vector) * m(scalar)
inline CPointEx operator * (const CPointEx& p, double m)
{
	return CPointEx(p.x * m, p.y * m);
}

// Returns m(scalar) * p(vector)
inline CPointEx operator * (double m, const CPointEx& p)
{
	return CPointEx(p.x * m, p.y * m);
}

// Returns p(vector) / m(scalar)
inline CPointEx operator / (const CPointEx& p, double m)
{
	ASSERT(m != 0);
	return CPointEx(p.x / m, p.y / m);
}

// Returns p1 + p2
inline CPointEx operator + (const CPointEx& p1, const CPointEx& p2)
{
	return CPointEx(p1.x + p2.x, p1.y + p2.y);
}

// Returns p1 - p2
inline CPointEx operator - (const CPointEx& p1, const CPointEx& p2)
{
	return CPointEx(p1.x - p2.x, p1.y - p2.y);
}

// returns the distance between p1,p2
inline double distance(const CPointEx& p1, const CPointEx& p2)
{
	CPointEx p = p2 - p1;
	return p.norm();
}

// returns the square distance between p1,p2
inline double sqr_distance(const CPointEx& p1, const CPointEx& p2)
{
	CPointEx p = p2 - p1;
	return p.sqr_norm();
}

// returns the dot product w * v
inline double dot(const CPointEx3D& w, const CPointEx3D& v)
{
	return (w.x * v.x + w.y * v.y + w.z * v.z);
}

// returns the dot product w * v
inline double operator * (const CPointEx3D& w, const CPointEx3D& v)
{
	return (w.x * v.x + w.y * v.y + w.z * v.z);
}

// returns p * scalar
inline CPointEx3D operator * (const CPointEx3D& p, double m)
{
	return CPointEx3D(p.x * m, p.y * m, p.z * m);
}

// returns scalar * p (vector)
inline CPointEx3D operator * (double m, const CPointEx3D& p)
{
	return CPointEx3D(p.x * m, p.y * m, p.z * m);
}

// Returns p(vector) / m(scalar)
inline CPointEx3D operator / (const CPointEx3D& p, double m)
{
	ASSERT(m != 0);
	return CPointEx3D(p.x / m, p.y / m, p.z / m);
}

// returns p1+p2
inline CPointEx3D operator + (const CPointEx3D& p1, const CPointEx3D& p2)
{
	return CPointEx3D(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}

// returns p1-p2
inline CPointEx3D operator - (const CPointEx3D& p1, const CPointEx3D& p2)
{
	return CPointEx3D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

// returns the distance between p1,p2
inline double distance(const CPointEx3D& p1, const CPointEx3D& p2)
{
	CPointEx3D p = p2 - p1;
	return p.norm();
}

// returns the square distance between p1,p2
inline double sqr_distance(const CPointEx3D& p1, const CPointEx3D& p2)
{
	CPointEx3D p = p2 - p1;
	return p.sqr_norm();
}

// returns the square of the norm of p
inline double sqr_norm(const CPointEx3D& p)
{
	return (p.x * p.x + p.y * p.y + p.z * p.z);
}

// returns the cross product p1 x p2
inline CPointEx3D cross(const CPointEx3D& p1, const CPointEx3D& p2)
{
	return CPointEx3D(p1.y*p2.z - p2.y*p1.z, -(p1.x*p2.z - p2.x*p1.z), p1.x*p2.y - p2.x*p1.y);
}

inline CPointEx3D mult(const CPointEx3D& p, double** M)
{
	return CPointEx3D(
		p.x * M[0][0] + p.y * M[0][1] + p.z * M[0][2] + /*p.w **/ M[0][3],
		p.x * M[1][0] + p.y * M[1][1] + p.z * M[1][2] + /*p.w **/ M[1][3],
		p.x * M[2][0] + p.y * M[2][1] + p.z * M[2][2] + /*p.w **/ M[2][3]);
}

inline CPointEx3D mult2(const CPointEx3D& p, double** M)
{
	return CPointEx3D(
		p.x * M[0][0] + p.y * M[1][0] + p.z * M[2][0] + /*p.w **/ M[3][0],
		p.x * M[0][1] + p.y * M[1][1] + p.z * M[2][1] + /*p.w **/ M[3][1],
		p.x * M[0][2] + p.y * M[1][2] + p.z * M[2][2] + /*p.w **/ M[3][2]);
}

inline CPointEx3D rotateX(const CPointEx3D& p, double theta)
{
	if (theta == 0) return p;
	double R[4][4];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			R[i][j] = 0.0;
			R[i][i] = 1.0;
		}

	R[1][1] = cosd(theta);
	R[1][2] = sind(theta);
	R[2][1] = -R[1][2];
	R[2][2] = R[1][1];

	return CPointEx3D(
		p.x * R[0][0] + p.y * R[1][0] + p.z * R[2][0] + /*p.w **/ R[3][0],
		p.x * R[0][1] + p.y * R[1][1] + p.z * R[2][1] + /*p.w **/ R[3][1],
		p.x * R[0][2] + p.y * R[1][2] + p.z * R[2][2] + /*p.w **/ R[3][2]);

}

inline CPointEx3D rotateY(const CPointEx3D& p, double theta)
{
	if (theta == 0) return p;
	double R[4][4];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			R[i][j] = 0.0;
			R[i][i] = 1.0;
		}

	R[0][0] = cosd(theta);
	R[0][2] = -sind(theta);
	R[2][0] = -R[0][2];
	R[2][2] = R[0][0];

	return CPointEx3D(
		p.x * R[0][0] + p.y * R[1][0] + p.z * R[2][0] + /*p.w **/ R[3][0],
		p.x * R[0][1] + p.y * R[1][1] + p.z * R[2][1] + /*p.w **/ R[3][1],
		p.x * R[0][2] + p.y * R[1][2] + p.z * R[2][2] + /*p.w **/ R[3][2]);

}

inline CPointEx3D rotateZ(const CPointEx3D& p, double theta)
{
	if (theta == 0) return p;
	double R[4][4];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			R[i][j] = 0.0;
			R[i][i] = 1.0;
		}

	R[0][0] = cosd(theta);
	R[0][1] = sind(theta);
	R[1][0] = -R[0][1];
	R[1][1] = R[0][0];

	return CPointEx3D(
		p.x * R[0][0] + p.y * R[1][0] + p.z * R[2][0] + /*p.w **/ R[3][0],
		p.x * R[0][1] + p.y * R[1][1] + p.z * R[2][1] + /*p.w **/ R[3][1],
		p.x * R[0][2] + p.y * R[1][2] + p.z * R[2][2] + /*p.w **/ R[3][2]);

}

// This function computes the Metacenter of a triagle defined by
// the three points entered in the pointer 'points'
inline
CPointEx GetMetacenter(const CPointEx* points)
{
	ASSERT(points != NULL); // are you kidding?
	return (points[0] + points[1] + points[2]) / 3.0;
	/*	CPointEx N = 0.5 * (points[1] + points[2]);
		return (points[0] + (2./3.) * (N - points[0]));*/
}


/***************************************************************************************************/
inline BOOL IntersectsSegment(double x, double y,
	double x0, double y0,
	double x1, double y1,
	double h2)
{
	// This function determines if a line segment defined
	// by P0(x0,y0) and P1(x1,y1) is intersected by a square
	// with center x, y and having edge length h = sqrt(h2)

	double dx = x1 - x0;
	double dy = y1 - y0;
	double f2 = dx * dx + dy * dy;
	if (f2 <= TOL)
		return FALSE;
	double t = ((x - x0) * dx + (y - y0) * dy) / f2;

	if ((t >= 0.0) && (t <= 1.0))
	{
		double f1 = fabs(dx) + fabs(dy);
		double dbox2 = h2 * f1*f1 / f2;
		f1 = x0 - x + t * dx;
		f2 = y0 - y + t * dy;
		double d2 = f1 * f1 + f2 * f2;
		if (d2 < dbox2)
			return TRUE;
	}

	return FALSE;
}

inline BOOL Intersection(const CPointEx& P0, const CPointEx& P1,
	const CPointEx& Q0, const CPointEx& Q1,
	double& u, double& v)
{
	// This function determines if two line segments defined
	// by P(u) = P0 + u(P1-P0) and Q(v) = Q0 + v(Q1-Q0) are
	// intersecting. If they are parallel the function returns FALSE.
	// Else u and v define the point of intersection.
	// If 0<= u,v <=1 then the intersection point lies between
	// each line segment. Else, it lies at the extension of the lines.
	// The tolerance accuracy is determined by the constant TOL.

	double _x1s = P0.x, _y1s = P0.y;
	double _x1e = P1.x, _y1e = P1.y;
	double _x2s = Q0.x, _y2s = Q0.y;
	double _x2e = Q1.x, _y2e = Q1.y;
	double _D = (_x2e - _x2s) * (_y1s - _y1e) - (_y2e - _y2s) * (_x1s - _x1e);
	if (fabs(_D) >= TOL)
	{
		double Ds = (_x1s - _x2s) * (_y1s - _y1e) - (_y1s - _y2s) * (_x1s - _x1e);
		double Dt = (_x2e - _x2s) * (_y1s - _y2s) - (_y2e - _y2s) * (_x1s - _x2s);

		u = Dt / _D;
		if ((u < 0 - TOL) || (u > 1 + TOL))
			return FALSE;

		v = Ds / _D;
		if ((v < 0 - TOL) || (v > 1 + TOL))
			return FALSE;

		// Make sure the params are between [0,1]
		if (u < 0) u = 0;
		else
			if (u > 1) u = 1;
		if (v < 0) v = 0;
		else
			if (v > 1) v = 1;

		return TRUE;
	}
	else
	{
		// we have to check if the end of the first line 
		// is the beginning for the second and vice versa
		if (DISTANCE_EX(P0, Q0) <= TOL)
		{
			u = v = 0.0;
			return TRUE;
		}
		else
			if (DISTANCE_EX(P0, Q1) <= TOL)
			{
				u = 0.0;
				v = 1.0;
				return TRUE;
			}
			else
				if (DISTANCE_EX(P1, Q0) <= TOL)
				{
					u = 1.0;
					v = 0.0;
					return TRUE;
				}
				else
					if (DISTANCE_EX(P1, Q1) <= TOL)
					{
						u = 1.0;
						v = 1.0;
						return TRUE;
					}
	}

	return FALSE;
}

inline BOOL Intersection(const CPointEx& P0, const CPointEx& P1,
	const CPointEx& Q0, const CPointEx& Q1,
	CPointEx& P,
	double& u, double& v)
{
	// This function determines if two line segments defined
	// by P(u) = P0 + u(P1-P0) and Q(v) = Q0 + v(Q1-Q0) are
	// intersecting. If they are parallel the function returns FALSE.
	// Else u and v define the point of intersection.
	// If 0<= u,v <=1 then the intersection point lies between
	// each line segment. Else, it lies at the extension of the lines.
	// The tolerance accuracy is determined by the constant TOL.

	double _x1s = P0.x, _y1s = P0.y;
	double _x1e = P1.x, _y1e = P1.y;
	double _x2s = Q0.x, _y2s = Q0.y;
	double _x2e = Q1.x, _y2e = Q1.y;
	double _D = (_x2e - _x2s) * (_y1s - _y1e) - (_y2e - _y2s) * (_x1s - _x1e);
	if (fabs(_D) >= TOL)
	{
		double Ds = (_x1s - _x2s) * (_y1s - _y1e) - (_y1s - _y2s) * (_x1s - _x1e);
		double Dt = (_x2e - _x2s) * (_y1s - _y2s) - (_y2e - _y2s) * (_x1s - _x2s);

		u = Dt / _D;
		if ((u < 0 - TOL) || (u > 1 + TOL))
			return FALSE;

		v = Ds / _D;
		if ((v < 0 - TOL) || (v > 1 + TOL))
			return FALSE;

		// Make sure the params are between [0,1]
		if (u < 0) u = 0;
		else
			if (u > 1) u = 1;
		if (v < 0) v = 0;
		else
			if (v > 1) v = 1;

		P.x = P0.x + u * (P1.x - P0.x); // get the actual point of intersection
		P.y = P0.y + u * (P1.y - P0.y); // get the actual point of intersection
		return TRUE;
	}
	else
	{
		// we have to check if the end of the first line 
		// is the beginning for the second and vice versa
		if (DISTANCE_EX(P0, Q0) <= TOL)
		{
			P = P0;
			u = v = 0.0;
			return TRUE;
		}
		else
			if (DISTANCE_EX(P0, Q1) <= TOL)
			{
				P = P0;
				u = 0.0;
				v = 1.0;
				return TRUE;
			}
			else
				if (DISTANCE_EX(P1, Q0) <= TOL)
				{
					P = P1;
					u = 1.0;
					v = 0.0;
					return TRUE;
				}
				else
					if (DISTANCE_EX(P1, Q1) <= TOL)
					{
						P = P1;
						u = 1.0;
						v = 1.0;
						return TRUE;
					}
	}

	return FALSE;
}

enum EXPAND_TYPE {
	EXPAND_NONE = 0,
	EXPAND_FIRST = 1,
	EXPAND_SECOND = 2,
	EXPAND_BOTH = 4,
};

inline BOOL Intersection(const CPointEx& P0, const CPointEx& P1,
	const CPointEx& Q0, const CPointEx& Q1,
	CPointEx& P,
	double& u, double& v,
	EXPAND_TYPE ExpType)
{
	// This function determines if two line segment defined
	// by P(u) = P0 + u(P1-P0) and Q(v) = Q0 + v(Q1-Q0) are
	// intersecting. If they are parallel the function return FALSE.
	// Else u and v define the point of intersection.
	// If 0<= u,v <=1 then the intersection point lies between
	// each line segment. Else, it lies at the extension of the lines.
	// The tolerance accuracy is determined by the constant TOL.

	double x1s = P0.x, y1s = P0.y;
	double x1e = P1.x, y1e = P1.y;
	double x2s = Q0.x, y2s = Q0.y;
	double x2e = Q1.x, y2e = Q1.y;
	double D = (x2e - x2s) * (y1s - y1e) - (y2e - y2s) * (x1s - x1e);
	if (fabs(D) >= TOL)
	{
		double Ds = (x1s - x2s) * (y1s - y1e) - (y1s - y2s) * (x1s - x1e);
		double Dt = (x2e - x2s) * (y1s - y2s) - (y2e - y2s) * (x1s - x2s);

		u = Dt / D;

		BOOL bExp1 = (u < 0 || u > 1);
		BOOL bAllow = (ExpType & (EXPAND_FIRST | EXPAND_BOTH));
		if (bExp1 && !bAllow)
			return FALSE;

		v = Ds / D;
		BOOL bExp2 = (v < 0 || v > 1);
		bAllow = (ExpType & (EXPAND_SECOND | EXPAND_BOTH));
		if (bExp2 && !bAllow)
			return FALSE;

		P.x = P0.x + u * (P1.x - P0.x); // get the actual point of intersection
		P.y = P0.y + u * (P1.y - P0.y); // get the actual point of intersection

		return TRUE;
	}
	else
	{
		// we have to check if the end of the first line 
		// is the beginning for the second and vice versa
		if (DISTANCE_EX(P0, Q0) <= TOL)
		{
			P = P0;
			u = v = 0.0;
			return TRUE;
		}
		else
			if (DISTANCE_EX(P0, Q1) <= TOL)
			{
				P = P0;
				u = 0.0;
				v = 1.0;
				return TRUE;
			}
			else
				if (DISTANCE_EX(P1, Q0) <= TOL)
				{
					P = P1;
					u = 1.0;
					v = 0.0;
					return TRUE;
				}
				else
					if (DISTANCE_EX(P1, Q1) <= TOL)
					{
						P = P1;
						u = 1.0;
						v = 1.0;
						return TRUE;
					}
	}
	return FALSE;
}

template <class T>
inline void MCLS_SWAP(T& n1, T& n2)
{
	T temp = n1;
	n1 = n2;
	n2 = temp;
}

// get the projection of the vector (point-p0) to the vector (p1-p0)
// u = ||projection|| / ||p1-p0||
inline CPointEx GetProjection(const CPointEx& p0,
	const CPointEx& p1,
	const CPointEx& point,
	double & u)
{
	double L1 = DISTANCE_EX(p0, p1);
	ASSERT(L1 > 0);

	CPointEx n = p1 - p0;
	double norm = NORM_EX(n);
	ASSERT(norm != 0);

	n /= norm;

	double L2 = n * (point - p0);
	u = L2 / L1;

	return (p0 + L2 * n);
}

