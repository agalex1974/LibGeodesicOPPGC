#include "pch.h"
#include "MCLSEXST.H"

BOOL CRectEx::Intersects(const CPointEx& p1, const CPointEx& p2,
	int& nIntersections, CPointEx& pnt1, CPointEx& pnt2, BOOL bExpandLine) const
{
	double u, v;
	EXPAND_TYPE ExpType = bExpandLine ? EXPAND_SECOND : EXPAND_NONE;
	CPointEx q1, q2, q3, q4, P;

	q1 = TopLeft();
	q2 = TopRight();
	q3 = BottomRight();
	q4 = BottomLeft();

	nIntersections = 0;
	if (Intersection(q1, q2, p1, p2, P, u, v, ExpType))
	{
		pnt1 = P;
		nIntersections = 1;
	}

	if (Intersection(q2, q3, p1, p2, P, u, v, ExpType))
	{
		if (nIntersections == 1)
		{
			pnt2 = P;
			nIntersections = 2;
			return TRUE;
		}
		else
		{
			pnt1 = P;
			nIntersections = 1;
		}
	}

	if (Intersection(q3, q4, p1, p2, P, u, v, ExpType))
	{
		if (nIntersections == 1)
		{
			pnt2 = P;
			nIntersections = 2;
			return TRUE;
		}
		else
		{
			pnt1 = P;
			nIntersections = 1;
		}
	}

	if (Intersection(q4, q1, p1, p2, P, u, v, ExpType))
	{
		if (nIntersections == 1)
		{
			pnt2 = P;
			nIntersections = 2;
			return TRUE;
		}
		else
		{
			pnt1 = P;
			nIntersections = 1;
		}
	}

	return (nIntersections > 0);
}


void SetMatrValue(Matrix4x4* mat, int n, double val)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			(*mat)[i][j] = val;
}

void SetMatrDiagValue(Matrix4x4* mat, int n, double val)
{
	for (int i = 0; i < n; i++)
		(*mat)[i][i] = val;
}

void SetMatrValues(Matrix4x4* mat, int n, double diag, double el)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				(*mat)[i][j] = diag;
			else
				(*mat)[i][j] = el;
		}
}

//multiply M = M1 * M2 square matrices of dim n
void mult(Matrix4x4* C, Matrix4x4 A, Matrix4x4 B, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			(*C)[i][j] = 0.0;
			for (int k = 0; k < n; k++)
				(*C)[i][j] += A[i][k] * B[k][j];
		}
}

void SetMatrValue(double*** mat, int n, double val)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			(*mat)[i][j] = val;
}

void SetMatrDiagValue(double*** mat, int n, double val)
{
	for (int i = 0; i < n; i++)
		(*mat)[i][i] = val;
}

void SetMatrValues(double*** mat, int n, double diag, double el)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				(*mat)[i][j] = diag;
			else
				(*mat)[i][j] = el;
		}
}

double ** alloc_array(int m, int n)
{
	double** data;

	data = new double*[m];

	if (data == NULL)
	{
		::MessageBox(0, (LPCSTR)L"Unable to allocate memory. Abnormal program termination!",
			(LPCSTR)L"Memory Error!", MB_ICONERROR);
		exit(1);
	}

	for (int i = 0; i < m; ++i)
	{
		data[i] = new double[n];
		if (data[i] == NULL)
		{
			::MessageBox(0, (LPCSTR)L"Unable to allocate memory. Abnormal program termination!",
				(LPCSTR)L"Memory Error!", MB_ICONERROR);
			exit(1);
		}
	}

	return(data); //returns pointer to matrix
}

void disp_array(double **mem, int m)
{
	ASSERT(mem != NULL);
	ASSERT(*mem != NULL);
	int i;
	for (i = 0; i < m; i++)
		delete[] * (mem + i);

	delete[] mem;

	mem = NULL;
}

//multiply M = M1 * M2 square matrices of dim n
void mult(double*** C, double** A, double** B, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			(*C)[i][j] = 0.0;
			for (int k = 0; k < n; k++)
				(*C)[i][j] += A[i][k] * B[k][j];
		}
}

double** GetTransMatrix(double x0, double y0, double z0,
	double l, double m, double n,
	double dx, double dy, double dz,
	double theta)
{
	double** M = alloc_array(4, 4);

	// Compute Rotation-Translation matrix
	double** T = alloc_array(4, 4);
	SetMatrValues(&T, 4, 1.0, 0);

	double cos_theta = COSD(theta);
	double sin_theta = SIND(theta);

	// diagonal elements
	T[0][0] = l * l * (1 - cos_theta) + cos_theta;
	T[1][1] = m * m * (1 - cos_theta) + cos_theta;
	T[2][2] = n * n * (1 - cos_theta) + cos_theta;


	double a10_1 = m * l * (1 - cos_theta);
	double a10_2 = n * sin_theta;

	double a20_1 = n * l * (1 - cos_theta);
	double a20_2 = m * sin_theta;

	double a21_1 = m * n * (1 - cos_theta);
	double a21_2 = l * sin_theta;

	// lower triangle
	T[1][0] = a10_1 + a10_2;
	T[2][0] = a20_1 - a20_2;
	T[2][1] = a21_1 + a21_2;

	// upper triangle
	T[0][1] = a10_1 - a10_2;
	T[0][2] = a20_1 + a20_2;
	T[1][2] = a21_1 - a21_2;

	// last column
	T[0][3] = dx;
	T[1][3] = dy;
	T[2][3] = dz;


	// set translation
	double** T1 = alloc_array(4, 4);
	SetMatrValues(&T1, 4, 1.0, 0);
	T1[0][3] = -x0;
	T1[1][3] = -y0;
	T1[2][3] = -z0;

	// set translation
	double** T2 = alloc_array(4, 4);
	SetMatrValues(&T2, 4, 1.0, 0);
	T2[0][3] = x0;
	T2[1][3] = y0;
	T2[2][3] = z0;

	double** temp = alloc_array(4, 4); // temp array

	// construct overall trans. matrix M
	mult(&temp, T, T1, 4);
	mult(&M, T2, temp, 4);

	// free mem
	disp_array(temp, 4);
	disp_array(T, 4);
	disp_array(T1, 4);
	disp_array(T2, 4);

	return M;
}

void GetTransMatrix(double*** mat,
	double x0, double y0, double z0,
	double l, double m, double n,
	double dx, double dy, double dz,
	double theta)
{
	double cos_theta = COSD(theta);
	double sin_theta = SIND(theta);

	// diagonal elements
	(*mat)[0][0] = l * l * (1 - cos_theta) + cos_theta;
	(*mat)[1][1] = m * m * (1 - cos_theta) + cos_theta;
	(*mat)[2][2] = n * n * (1 - cos_theta) + cos_theta;

	double a10_1 = m * l * (1 - cos_theta);
	double a10_2 = n * sin_theta;

	double a20_1 = n * l * (1 - cos_theta);
	double a20_2 = m * sin_theta;

	double a21_1 = m * n * (1 - cos_theta);
	double a21_2 = l * sin_theta;

	// lower triangle
	(*mat)[1][0] = a10_1 + a10_2;
	(*mat)[2][0] = a20_1 - a20_2;
	(*mat)[2][1] = a21_1 + a21_2;

	// upper triangle
	(*mat)[0][1] = a10_1 - a10_2;
	(*mat)[0][2] = a20_1 + a20_2;
	(*mat)[1][2] = a21_1 - a21_2;

	// last column
	(*mat)[0][3] = -x0 * (*mat)[0][0] - y0 * (*mat)[0][1] - z0 * (*mat)[0][2] + dx + x0;
	(*mat)[1][3] = -x0 * (*mat)[1][0] - y0 * (*mat)[1][1] - z0 * (*mat)[1][2] + dy + y0;
	(*mat)[2][3] = -x0 * (*mat)[2][0] - y0 * (*mat)[2][1] - z0 * (*mat)[2][2] + dz + z0;

	// last row
	(*mat)[3][0] = 0.0;
	(*mat)[3][1] = 0.0;
	(*mat)[3][2] = 0.0;
	(*mat)[3][3] = 1.0;
}


void GetTransMatrix2(double*** mat,
	double x0, double y0, double z0,
	double l, double m, double n,
	double dx, double dy, double dz,
	double theta)
{
	// Compute Rotation-Translation matrix
	double** T = alloc_array(4, 4);
	SetMatrValues(&T, 4, 1.0, 0);

	double cos_theta = COSD(theta);
	double sin_theta = SIND(theta);

	// diagonal elements
	T[0][0] = l * l * (1 - cos_theta) + cos_theta;
	T[1][1] = m * m * (1 - cos_theta) + cos_theta;
	T[2][2] = n * n * (1 - cos_theta) + cos_theta;


	double a10_1 = m * l * (1 - cos_theta);
	double a10_2 = n * sin_theta;

	double a20_1 = n * l * (1 - cos_theta);
	double a20_2 = m * sin_theta;

	double a21_1 = m * n * (1 - cos_theta);
	double a21_2 = l * sin_theta;

	// lower triangle
	T[1][0] = a10_1 + a10_2;
	T[2][0] = a20_1 - a20_2;
	T[2][1] = a21_1 + a21_2;

	// upper triangle
	T[0][1] = a10_1 - a10_2;
	T[0][2] = a20_1 + a20_2;
	T[1][2] = a21_1 - a21_2;

	// last column
	T[0][3] = dx;
	T[1][3] = dy;
	T[2][3] = dz;


	// set translation
	double** T1 = alloc_array(4, 4);
	SetMatrValues(&T1, 4, 1.0, 0);
	T1[0][3] = -x0;
	T1[1][3] = -y0;
	T1[2][3] = -z0;

	// set translation
	double** T2 = alloc_array(4, 4);
	SetMatrValues(&T2, 4, 1.0, 0);
	T2[0][3] = x0;
	T2[1][3] = y0;
	T2[2][3] = z0;

	double** temp = alloc_array(4, 4); // temp array

	// construct overall trans. matrix M
	mult(&temp, T, T1, 4);
	mult(mat, T2, temp, 4);

	// free mem
	disp_array(temp, 4);
	disp_array(T, 4);
	disp_array(T1, 4);
	disp_array(T2, 4);
}

void GetTransMatrix(Matrix4x4* mat,
	double x0, double y0, double z0,
	double l, double m, double n,
	double dx, double dy, double dz,
	double theta)
{
	// Compute Rotation-Translation matrix
	Matrix4x4 T;
	SetMatrValues(&T, 4, 1.0, 0);

	double cos_theta = COSD(theta);
	double sin_theta = SIND(theta);

	double one_cos_theta = 1 - cos_theta;

	// diagonal elements
	T[0][0] = l * l * (one_cos_theta)+cos_theta;
	T[1][1] = m * m * (one_cos_theta)+cos_theta;
	T[2][2] = n * n * (one_cos_theta)+cos_theta;


	double a10_1 = m * l * (one_cos_theta);
	double a10_2 = n * sin_theta;

	double a20_1 = n * l * (one_cos_theta);
	double a20_2 = m * sin_theta;

	double a21_1 = m * n * (one_cos_theta);
	double a21_2 = l * sin_theta;

	// lower triangle
	T[1][0] = a10_1 + a10_2;
	T[2][0] = a20_1 - a20_2;
	T[2][1] = a21_1 + a21_2;

	// upper triangle
	T[0][1] = a10_1 - a10_2;
	T[0][2] = a20_1 + a20_2;
	T[1][2] = a21_1 - a21_2;

	// last column
	T[0][3] = dx;
	T[1][3] = dy;
	T[2][3] = dz;


	// set translation
	Matrix4x4 T1;
	SetMatrValues(&T1, 4, 1.0, 0);
	T1[0][3] = -x0;
	T1[1][3] = -y0;
	T1[2][3] = -z0;

	// set translation
	Matrix4x4 T2;
	SetMatrValues(&T2, 4, 1.0, 0);
	T2[0][3] = x0;
	T2[1][3] = y0;
	T2[2][3] = z0;

	Matrix4x4 temp; // temp array

	// construct overall trans. matrix M
	mult(&temp, T, T1, 4);
	mult(mat, T2, temp, 4);
}


void GetRotMatrix(Matrix4x4* mat,
	double x0, double y0, double z0,
	double l, double m, double n,
	double theta)
{
	// Compute Rotation-Translation matrix
	static Matrix4x4 T;
	SetMatrValues(&T, 4, 1.0, 0);

	double cos_theta = COSD(theta);
	double sin_theta = SIND(theta);

	double one_cos_theta = 1 - cos_theta;

	// diagonal elements
	T[0][0] = l * l * (one_cos_theta)+cos_theta;
	T[1][1] = m * m * (one_cos_theta)+cos_theta;
	T[2][2] = n * n * (one_cos_theta)+cos_theta;


	double a10_1 = m * l * (one_cos_theta);
	double a10_2 = n * sin_theta;

	double a20_1 = n * l * (one_cos_theta);
	double a20_2 = m * sin_theta;

	double a21_1 = m * n * (one_cos_theta);
	double a21_2 = l * sin_theta;

	// lower triangle
	T[1][0] = a10_1 + a10_2;
	T[2][0] = a20_1 - a20_2;
	T[2][1] = a21_1 + a21_2;

	// upper triangle
	T[0][1] = a10_1 - a10_2;
	T[0][2] = a20_1 + a20_2;
	T[1][2] = a21_1 - a21_2;

	// last column
	T[0][3] = 0;
	T[1][3] = 0;
	T[2][3] = 0;


	// set translation
	static Matrix4x4 T1;
	SetMatrValues(&T1, 4, 1.0, 0);
	T1[0][3] = -x0;
	T1[1][3] = -y0;
	T1[2][3] = -z0;

	// set translation
	static Matrix4x4 T2;
	SetMatrValues(&T2, 4, 1.0, 0);
	T2[0][3] = x0;
	T2[1][3] = y0;
	T2[2][3] = z0;

	static Matrix4x4 temp; // temp array

	// construct overall trans. matrix M
	mult(&temp, T, T1, 4);
	mult(mat, T2, temp, 4);
}

// set a comparison function for QuickSort; sort ascent
int SortAscentInt(const void *elem1, const void *elem2)
{
	int n1 = *((int*)elem1);
	int n2 = *((int*)elem2);

	if (n1 > n2)
		return 1;
	else
		if (n1 == n2)
			return 0;
		else
			return -1;
}

// set a comparison function for QuickSort; sort ascent
int SortAscentDoubleF(const void *elem1, const void *elem2)
{
	double n1 = *((double*)elem1);
	double n2 = *((double*)elem2);

	if (n1 > n2)
		return 1;
	else
		if (n1 == n2)
			return 0;
		else
			return -1;
}

void SortAscentDouble(double* data, int size)
{
	VERIFY(data != NULL);
	VERIFY(size > 0);
	/* Sort using Quicksort algorithm: */
	qsort((void *)data, (size_t)size, sizeof(double),
		(int(*)(const void*, const void*))SortAscentDoubleF);

}


// set a comparison function for QuickSort; sort descent
int SortDescentDoubleF(const void *elem1, const void *elem2)
{
	double n1 = *((double*)elem1);
	double n2 = *((double*)elem2);

	if (n1 < n2)
		return 1;
	else
		if (n1 == n2)
			return 0;
		else
			return -1;
}

void SortDescentDouble(double* data, int size)
{
	VERIFY(data != NULL);
	VERIFY(size > 0);
	/* Sort using Quicksort algorithm: */
	qsort((void *)data, (size_t)size, sizeof(double),
		(int(*)(const void*, const void*))SortDescentDoubleF);

}

