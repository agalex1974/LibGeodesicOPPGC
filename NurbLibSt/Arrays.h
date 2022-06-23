#ifndef __INTEL_COMPILER
#pragma once
#else
#pragma warning ( disable : 1885 )
#endif

/**********************************************************************************
 * Arrays library 
 * ~~~~~~~~~~~~~~~
 *
 * Purpose:
 *		To provide an efficient interface for communicating
		with the IMSL numerical library, which is written in FORTRAN.
		The problem occurs in dynamic allocation of NxM matrices
		since the interface functions of IMSL assume matrix data passed
		as vector arrays.
 *
 * Solution:
		A NxM matrix is stored as a vector of NxM elements.
 *
 * Drawback:
 *		A small extra computational burden is required
		in order to access an element of the array.
 *
 * Author:
		Philip N. Azariadis
 *
 * Date created: 6/8/98
 * Revision1: 25/5/99 (Phillip N. Azariadis)
 * Revision2: 7/7/99 - Added CVector::FlipDataOrder() (Phillip N. Azariadis)
 * Revision3: 4/8/99 - Added CVector::Add(T element) (Phillip N. Azariadis)
 * Revision4: 1/2/00 - Bug fixed at CMatrix::ReSize(..) (Phillip N. Azariadis)
 * Revision5: 10/3/00 - Bug fixed at CVector::operator = (..) (Phillip N. Azariadis)
 * Revision6: 15/6/01 - Added CVector::Add(CVector<>) (Phillip N. Azariadis)
 * Revision7: 20/6/01 - Added CVector::TrimLeft(int index) (Phillip N. Azariadis)
 * Revision7: 20/6/01 - Added CVector::TrimRight(int index) (Phillip N. Azariadis)
 * Revision7a: 12/12/01 - Fix CVector::TrimRight(int index) (Phillip N. Azariadis)
 * Revision7b: 3/05/02 - Added CMatrix::FreeData() (Should be called when SetData is used, and before calling again SetData. Phillip N. Azariadis)
 * Revision7b: 3/05/02 - Added CVector::FreeData() (Should be called when SetData is used, and before calling again SetData. Phillip N. Azariadis)
 * Revision8: 28/05/02 - Re-written CVector::ReSize(int NewSize, int GrowRate = 6) (More clear implementation. Phillip N. Azariadis)
 * Revision8: 28/05/02 - Bug fixed CVector::Add(CVector<>) (Phillip N. Azariadis)
 * Revision8: 28/05/02 - Bug fixed CVector::Add(T element) (Phillip N. Azariadis)
 * Revision9: 02/07/02 - Added CVector::T* GetDataAt(int index) (Phillip N. Azariadis)
 * Revision10: 26/08/02 - Added CVector::CVector(const CVector<T>& source) (Phillip N. Azariadis) (copy constructor)
 * Revision10: 26/08/02 - Added CMatrix::CMatrix(const CMatrix<T>& source) (Phillip N. Azariadis) (copy constructor)
 * Revision10: 26/08/02 - Re-written CMatrix::ReSize(int NewRows, int NewCols) (More clear and efficient implementation. Phillip N. Azariadis)
 * Revision11: 28/12/02 - Added CVector::IsEmpty() const  (Phillip N. Azariadis)
 * Revision11: 28/12/02 - Added CMatrix::IsEmpty() const  (Phillip N. Azariadis)
 * Revision11: 28/12/02 - Bug fixed CVector::Resize() (Phillip N. Azariadis)
 * Revision12: 16/06/03 - Added CVector::void operator += (T s) (Phillip N. Azariadis)
 * Revision12: 16/06/03 - Added CMatrix::void operator += (T s) (Phillip N. Azariadis)
 * Revision13: 28/08/03 - I replaced 'memcpy' with '= operator' for type safe memory copy - applications should define = operator properly (Phillip N. Azariadis)
 * Revision14a: 19/01/04 - Added CVector::operator T*()
 * Revision14b: 19/01/04 - Added CVector::operator const T*()
 * Revision14c: 19/01/04 - Added CMatrix::operator T*()
 * Revision14d: 19/01/04 - Added CMatrix::operator const T*()
 * Revision15: 15/04/04 -  void ReSize(int NewSize, int GrowRate = 6): When NewSize =0, the memory is not freed is max_size<1024 for speed gains
 * Revision16: 26/07/04 -  void CVector::Copy(const CVector<T>& source, int nStart, int nSize) // zero-based nStart; WARNING this is NOT TESTED!!!!
 * Revision17: 03/02/06 - Bug fixed CVector::TrimLeft(int index)
 * Revision18: 14/12/06 - Bug fixed CVector::IsEmpty() const 
 * Revision19: 01/12/09 - OMP support added with USE_OMP_AZA flag (still under investigation --> do not use for serious work)
 * Revision20: 02/05/13 - Added CVector::Sum()
 * Revision21a: 04/06/13 - Added CVector::PushBack()
 * Revision21b: 02/05/13 - Added CVector::PushFront
 * Revision22: 02/05/13 - OMP support is removed.
 * Revision23: 06/09/13 - Added CVector::DelAt(int n); warning this is slow and should not be used with large data; use CList with large data and insert/del operations
 * Revision24: 25/01/14 - Added CVector::Release() & CMatrix::Release()
 * Revision24: 28/09/14 - Added CVector::Contains(const T& x) const 
 *
 **********************************************************************************
 */


#ifndef _ARRAYS_HPP_
#define _ARRAYS_HPP_
#include <cmath>
#ifndef _INC_MEMORY
#include <memory.h>
#endif

#ifndef _AFX

#ifndef BOOL
typedef int                 BOOL;
//#define BOOL bool
#endif

#ifndef BYTE
#define BYTE short
#endif

#ifndef UINT
#define UINT unsigned long
#endif

#ifndef COLORREF
#define COLORREF unsigned long
#endif

#ifndef TRUE
#define TRUE true
#endif
#ifndef FALSE
#define FALSE false
#endif

#ifndef NULL
#define NULL 0
#endif

#ifndef _ASSERT
#define _ASSERT
#endif

#ifndef ASSERT
#define ASSERT _ASSERT
#endif

#endif // _AFX

const int USE_OMP_AZA_LMT = 0; // this is set only when you compile with OMP

template <class T>
class CVector
{
  public:
	CVector(int MaxSize = 0, int GrowRate = 6)				
	{ m_maxSize = MaxSize; m_maxAlloc = MaxSize; m_growRate = GrowRate; 
		m_pData = NULL; if( m_maxSize > 0) m_pData = new T[m_maxSize]; 
		if(m_growRate < 1) m_growRate = 1; m_bAutoDelete = TRUE;}

	CVector(const CVector<T>& source)
	{	m_maxSize = 0; m_maxAlloc = 0; m_growRate = 6; 
		m_pData = NULL; m_bAutoDelete = TRUE;
		if(source.m_pData == NULL || source.m_maxSize <=0){	return;}

		ReSize(source.m_maxSize, source.m_growRate); // need allocation
		
		_ASSERT(source.m_maxSize > 0);
		_ASSERT(m_pData != NULL); _ASSERT(source.m_pData != NULL);
		for(int i=0; i<m_maxSize; i++) m_pData[i] = source.m_pData[i];
	}

	virtual ~CVector()
	{ if(m_bAutoDelete && m_pData)  { delete [] m_pData;  m_pData = NULL; } }
	inline void Release() { ReSize(0); }

    inline void SetSize(int NewSize, int GrowRate = 6){ReSize(NewSize,GrowRate);} // for compatibility with MS "CArray"
    inline void ReSize(int NewSize, int GrowRate = 6)
    {
		m_growRate = GrowRate; if(m_growRate < 1) m_growRate = 1;
		
		if(m_pData != NULL) // we have already allocation
		{
			_ASSERT(m_pData);
			if(NewSize == 0) // destroy vector data
			{
				if(m_maxSize<=1024)
				{m_maxSize = 0; m_bAutoDelete = TRUE; } // do not destroy memory data but reset the array for speed gains
				else
				{
					delete [] m_pData;	m_pData = NULL;	m_bAutoDelete = TRUE;  
					m_maxSize = m_maxAlloc = 0; 
				}
			}
			else
			if(m_maxSize == NewSize ||
				NewSize <= m_maxAlloc) // already allocated
			{
				m_maxSize = NewSize; m_bAutoDelete = TRUE; 
			}
			else
			{
				// need reallocation
				delete [] m_pData; m_pData = NULL;
				m_maxAlloc = (m_maxSize = NewSize) + m_growRate;
				m_pData = new T[m_maxAlloc];
				_ASSERT(m_pData);
				m_bAutoDelete = TRUE;
			}
		}
		else // need reallocation
		{
			_ASSERT(m_pData == NULL);
			if(NewSize == 0)
				{ m_maxSize = m_maxAlloc = 0; m_bAutoDelete = TRUE;  return; }
			m_maxAlloc = (m_maxSize = NewSize) + m_growRate;
			m_pData = new T[m_maxAlloc];
			_ASSERT(m_pData);
			m_bAutoDelete = TRUE;
		}
    }
   
	inline T& operator [](int index)
      { _ASSERT(CheckIndex(index)); return m_pData[index]; }
	inline const T& operator [](int index) const
      { _ASSERT(CheckIndex(index)); return m_pData[index]; }

	inline T& operator [](unsigned index)
      { _ASSERT(CheckIndex(index)); return m_pData[index]; }
	inline const T& operator [](unsigned index) const
      { _ASSERT(CheckIndex(index)); return m_pData[index]; }

	inline T& operator [](long index)
      { _ASSERT(CheckIndex(index)); return m_pData[index]; }

	inline const T& operator [](long index) const
      { _ASSERT(CheckIndex(index)); return m_pData[index]; }
    
	inline operator T*() {return (T*)m_pData;}
	inline operator const T*() {return (T*) m_pData;}

	inline void Add(const T& element)
	{ 
		if(m_pData == NULL) // need allocation
		{
			m_maxAlloc = m_maxSize + m_growRate;
			m_pData = new T[m_maxAlloc];
			_ASSERT(m_pData);
			m_bAutoDelete = TRUE;
		}
		if(m_maxSize >= m_maxAlloc)
		{
			_ASSERT(m_pData);
			T* new_pData = new T[m_maxSize + m_growRate];
			_ASSERT(new_pData);
			for(int i=0; i<m_maxSize; i++) new_pData[i] = m_pData[i];
			delete[] m_pData; m_pData = new_pData; m_maxAlloc = m_maxSize + m_growRate; 
			m_maxSize++; m_bAutoDelete = TRUE;	
		}
		else
			m_maxSize++;
	  (*this)[m_maxSize-1] = element;
	}

	inline void PushBack(const T& element)
	{
		Add(element);
	}

	inline void PushFront(const T& element)
	{
		if(m_pData == NULL) // need allocation
		{
			m_maxAlloc = m_maxSize + m_growRate;
			m_pData = new T[m_maxAlloc];
			_ASSERT(m_pData);
			m_bAutoDelete = TRUE;
		}
		if(m_maxSize >= m_maxAlloc)
		{
			_ASSERT(m_pData);
			T* new_pData = new T[m_maxSize + m_growRate];
			_ASSERT(new_pData);

			for(int i=m_maxSize-1; i>=0; i--) new_pData[i+1] = m_pData[i];
			new_pData[0] = element;
			delete[] m_pData; m_pData = new_pData; m_maxAlloc = m_maxSize + m_growRate; 
			m_maxSize++; m_bAutoDelete = TRUE;	
		}
		else
		{
			for(int i=m_maxSize-1; i>=0; i--) m_pData[i+1] = m_pData[i];
			m_pData[0] = element;
			m_maxSize++;
		}
	}

	inline void Add(const CVector<T>& source)
	{ 
		if(m_pData == NULL) // need allocation
		{
			m_maxAlloc = m_maxSize + m_growRate;
			m_pData = new T[m_maxAlloc];
			_ASSERT(m_pData);
			m_bAutoDelete = TRUE;
		}
		if(source.m_pData == NULL){	ASSERT(FALSE);return;}
		int newSize = m_maxSize + source.m_maxSize;
		if(newSize >= m_maxAlloc){
		m_maxAlloc = newSize + m_growRate;
		T* new_pData = new T[m_maxAlloc];
		int i=0;
		for(i=0; i<m_maxSize; i++) new_pData[i] = m_pData[i]; // copy existing data
		for(i=0; i<source.m_maxSize; i++) new_pData[i+m_maxSize] = source.m_pData[i]; // copy source data

		delete[] m_pData; m_pData = new_pData;; 
		m_maxSize = newSize; 
		}
	  else
	  {
		for(int i=0; i<source.m_maxSize; i++) m_pData[i+m_maxSize] = source.m_pData[i]; // copy source data
		m_maxSize = newSize;
	  }
	  m_bAutoDelete = TRUE;
	}

	inline void DelAt(int n)
	{
		_ASSERT(CheckIndex(n));
		if(!(n>=0 && n<GetSize()) )
			return; // this will assert in debug

		T* new_pData = new T[m_maxAlloc];
		for(int i=0, i2=0; i<GetSize(); i++)
		{
			if(i==n)
				continue; // exclude n
			new_pData[i2++] = m_pData[i];
		}
		delete[] m_pData; m_pData = new_pData;; 
		m_maxSize--; 
	}



    inline int GetSize() const
      { return m_maxSize; }
    inline int CheckIndex(int index) const
      { return  (index >= 0 && index < m_maxSize) ? 1 : 0; }

	inline BOOL IsEmpty() const 
	{ 
#ifdef _DEBUG 
		if(m_maxAlloc == 0) _ASSERT(m_pData == NULL);
		if(m_maxSize > 0) _ASSERT(m_pData != NULL);
#endif
		return (m_maxSize==0);
	}

	inline T* GetData()
		{	return (T*) m_pData;	}

	inline const T* GetData() const 
		{	return (const T*) m_pData;	}

	inline T* GetDataAt(int index)
      { _ASSERT(CheckIndex(index)); return (T*) (m_pData+index); }

	inline void SetData(T* pData, int nSize, BOOL bAutoDelete = FALSE)
		{	_ASSERT(pData != NULL); if(m_pData != NULL) { delete [] m_pData; }
			m_pData = pData; m_maxSize = nSize;	m_bAutoDelete = bAutoDelete;
		}
	inline void FreeData() { m_pData = NULL; m_bAutoDelete = TRUE;m_maxSize = 0; }

	inline void SetValue(const T& val)
		{	_ASSERT(m_pData != NULL);
			for(int i = 0; i < m_maxSize; i++)m_pData[i] = val;
		}

	inline void SetVal(const T& val) { SetValue(val); }

	inline void Copy(const CVector<T>& source){*this = source;} // for compatibility with MS "CArray"

	inline void Copy(const CVector<T>& source, int nStart, int nSize) // zero-based nStart; WARNING this is NOT TESTED!!!!
	{
		if(nSize==0) return; // nothing to copy
		if(nStart<0) // copy entiry source
		{
			Copy(source);
			return;
		}

		if(nSize<0 || nStart+nSize>source.m_maxSize) nSize = source.m_maxSize-nStart; // copy all source starting from nStart

		ReSize(nSize);
		for(int i=0; i<nSize; i++)	m_pData[i] = source.m_pData[nStart+i]; // copy existing data

	} // for compatibility with MS "CArray"

	inline void operator = (const CVector<T>& source)
	{
		if(source.m_pData == NULL || source.m_maxSize ==0){	ReSize(0);return;}

		if(	m_pData == NULL || 
			m_maxSize != source.m_maxSize)
			ReSize(source.m_maxSize, source.m_growRate); // need allocation
		
		if(source.m_maxSize > 0)
		{
			for(int i=0; i<m_maxSize; i++) m_pData[i] = source.m_pData[i]; // copy existing data
		}
#ifdef _DEBUG
		else
			_ASSERT(m_pData == NULL);
#endif
	}

	inline void CopyData(const T* source, int size)
	{
		_ASSERT(source != NULL);
		_ASSERT(size > 0);
		if (m_pData != NULL)
		{
			delete[] m_pData;	m_pData = NULL;	
			m_maxSize = m_maxAlloc = 0; m_bAutoDelete = TRUE;  
		}
		ReSize(size);
		memcpy(m_pData, source, sizeof(T) * m_maxSize);
	}


	inline void operator *= (double s)
		{	_ASSERT(m_pData != NULL);
			for(int i = 0; i < m_maxSize; i++)	m_pData[i] *= s;
		}

	inline void operator += (const T& s)
		{	_ASSERT(m_pData != NULL);
			for(int i = 0; i < m_maxSize; i++)	m_pData[i] += (T)s;
		}

	inline void operator /= (double s)
		{	_ASSERT(m_pData != NULL);	_ASSERT(s != 0);
			for(int i = 0; i < m_maxSize; i++)	m_pData[i] /= s;
		}

	inline void FlipDataOrder()
		{	_ASSERT(m_pData != NULL);
			CVector<T> tmp(m_maxSize, m_growRate);
			for(int i = 0; i < m_maxSize; i++)tmp.m_pData[m_maxSize-1 - i] = m_pData[i];
			*this = tmp;
		}

	inline void TrimLeft(int index) // index is unit based; 1 <= index <= m_max_size
		{	_ASSERT(m_pData != NULL); _ASSERT(CheckIndex(index-1));
			_ASSERT(index != m_maxSize); // destroy data with a different way!
			if(index == 0) return;
			CVector<T> tmp(m_maxSize-index, m_growRate);
			for(int i = 0; i < m_maxSize-index; i++)tmp.m_pData[i] = m_pData[i+index]; // bug fix 3/02/2006
			*this = tmp;
		}

	inline void TrimRight(int index) // index is unit based; 1 <= index <= m_max_size
		{	_ASSERT(m_pData != NULL); _ASSERT(CheckIndex(index-1));
			_ASSERT(index != m_maxSize); // destroy data with a different way!
			_ASSERT(index != 0); // nothing to trim
			if(index == 0) return;
			m_maxSize -= index; // just reduce the size of the vector
			if(m_maxSize <= 0) ReSize(0,m_growRate); // delete array
		}

	inline T Sum() const {
			_ASSERT(m_pData != NULL);
			if(!m_pData ) return 0;
			T sum = m_pData[0];
			for(int i=1; i<m_maxSize; i++)
				sum += m_pData[i];
			return sum;
	}

	inline BOOL Contains(const T& x) const {
		if (!m_pData) return FALSE;
		for (int i = 0; i < m_maxSize; i++)
			if (m_pData[i] == x)
				return TRUE;
		return FALSE;
	}

  protected:
	BOOL m_bAutoDelete;
    int m_maxSize, m_maxAlloc, m_growRate;
    T* m_pData;
};

// define some types of vectors
typedef CVector<BYTE> CByteVector;
typedef CVector<BOOL> CBOOLVector;
#define CBoolVector CBOOLVector
typedef CVector<int> CIntVector; // –2,147,483,648 to 2,147,483,647
typedef CVector<unsigned int> CUIntVector; //0 to 4,294,967,295
typedef CVector<unsigned int> CUINTVector; //0 to 4,294,967,295;  for compatibility with old programms
typedef CVector<__int8> CInt8Vector; // -128 to 127
typedef CVector<unsigned __int8> CUInt8Vector; // 0 to 255
typedef CVector<unsigned __int16> CUInt16Vector; // 0 to 65535
typedef CVector<__int16> CInt16Vector; // –32,768 to 32,767
//#define CUINTVector CUIntVector
typedef CVector<long> CLongVector;
typedef CVector<float> CFloatVector;
typedef CVector<double> CDoubleVector;
typedef CVector<long double> CLDoubleVector;
typedef CVector<COLORREF> CColorVector;


// ------------------------------------------------------------------------------------------------
template <class T>
class CMatrix
{
  public:
    CMatrix()
     {	m_pData = NULL; m_bAutoDelete = TRUE; m_maxCols = m_maxRows = m_maxAlloc = 0; }

	CMatrix(const CMatrix<T>& source)
	{
		m_pData = NULL; m_bAutoDelete = TRUE; m_maxCols = m_maxRows = m_maxAlloc = 0;
		if(source.m_pData == NULL) return;

		_ASSERT(source.m_maxRows*source.m_maxCols>0);
		ReSize(source.m_maxRows, source.m_maxCols);
		_ASSERT(m_pData != NULL); _ASSERT(m_maxRows == source.m_maxRows);
		_ASSERT(m_maxCols == source.m_maxCols);
		for(int i=0; i<m_maxRows * m_maxCols; i++) m_pData[i] = source.m_pData[i]; // copy existing data
	}

    CMatrix(int MaxRows, int MaxCols)
     {
       m_maxCols = MaxCols;	m_maxRows = MaxRows;
	   m_maxAlloc = m_maxRows * m_maxCols;
       m_pData = new T[m_maxAlloc];
	   m_bAutoDelete = TRUE;
     }
    ~CMatrix()
      { if(m_bAutoDelete && m_pData) delete [] m_pData;  m_pData = NULL; }
	inline void Release() { ReSize(0,0); }
    inline void ReSize(int NewRows, int NewCols)
       {
		 int NewAlloc = NewRows * NewCols;
		 // if new size equals zero deallocate and return
		 if(NewAlloc <= 0)
		 {	 if(m_pData) delete [] m_pData;	m_pData = NULL;	m_bAutoDelete = TRUE;  
					m_maxCols = m_maxRows = m_maxAlloc = 0;return;
		 }

		 if(m_maxRows == NewRows && m_maxCols == NewCols)
			   return; // nothing to change

		 if(m_pData != NULL)
		 {
			 // check if we have already enough storage
			 if(NewAlloc<=m_maxAlloc)
			 {   m_maxCols = NewCols;	m_maxRows = NewRows;
				 m_bAutoDelete = TRUE;
				 return;
			 }
			 else // not sufficient; need reallocation
			 {	 delete [] m_pData;	m_pData = NULL;	m_bAutoDelete = TRUE; m_maxCols = m_maxRows = m_maxAlloc = 0;}
		 }

         _ASSERT(m_pData==NULL);

		 // we need to allocate new data
         m_maxCols = NewCols;
         m_maxRows = NewRows;
		 m_maxAlloc = NewAlloc;
		 _ASSERT(m_maxCols>0 && m_maxRows>0);
         m_pData = new T[m_maxAlloc]; _ASSERT(m_pData!=NULL);
		 m_bAutoDelete = TRUE;
       }

     inline T& operator ()(int row, int col)
      { _ASSERT(CheckRowCol(row,col)); return m_pData[row + m_maxRows * col]; }

     inline T operator () (int row, int col) const 
      { _ASSERT(CheckRowCol(row,col)); return m_pData[row + m_maxRows * col]; }

    inline int CheckRowCol(int row, int col) const
      { return  (row >= 0 && row < m_maxRows &&
                 col >= 0 && col < m_maxCols) ? 1 : 0; }
    inline int GetRows() const
      { return m_maxRows; }
    inline int GetCols() const
      { return m_maxCols; }

	inline T* GetData()
		{ return (T*) m_pData; }

	inline T* GetData() const
		{ return (T*) m_pData; }

	inline operator T*() {return (T*)m_pData;}
//	inline operator const T*() {return const (T*) m_pData;}

	inline T* GetData(int row, int col)
		{ _ASSERT(CheckRowCol(row,col)); return (T*) (m_pData+row + m_maxRows * col); }

	inline T* GetRow(int row)
		{ _ASSERT(row < m_maxRows); return (T*) (m_pData + row * m_maxCols); }

	inline T* GetRow(int row) const
		{ _ASSERT(row < m_maxRows); return (T*) (m_pData + row * m_maxCols); }

	inline void SetData(T* pData, int nRows, int nCols)
		{	_ASSERT(pData != NULL); if(m_pData != NULL) { delete [] m_pData; }
			m_pData = pData; m_maxRows = nRows; m_maxCols = nCols;m_maxAlloc = m_maxRows*m_maxCols;
			m_bAutoDelete = FALSE;
		}
	inline void FreeData() { m_pData = NULL; m_bAutoDelete = TRUE;m_maxCols = m_maxRows = m_maxAlloc = 0;}

	inline BOOL IsEmpty() const 
	{ 
#ifdef _DEBUG 
		if(m_maxRows*m_maxCols == 0) _ASSERT(m_pData == NULL);
		else
		if(m_maxRows*m_maxCols > 0) _ASSERT(m_pData != NULL);
		else
			_ASSERT(FALSE); // invalid m_maxRows*m_maxCols
#endif
		return (m_pData == NULL);
	}

	inline void operator = (const CMatrix<T>& source)
	{
		_ASSERT(source.m_pData != NULL);

		if(	m_pData == NULL || 
			m_maxRows != source.m_maxRows ||
			m_maxCols != source.m_maxCols )
				ReSize(source.m_maxRows, source.m_maxCols);
		_ASSERT(m_pData != NULL); _ASSERT(m_maxRows == source.m_maxRows);
		_ASSERT(m_maxCols == source.m_maxCols);

//		memcpy(m_pData, source.m_pData, sizeof(T) * m_maxRows * m_maxCols);
		for(int i=0; i<m_maxRows * m_maxCols; i++) m_pData[i] = source.m_pData[i]; // copy existing data
	}

	inline void CopyData(const T* source, int rows, int cols)
	{
		_ASSERT(source != NULL);
		_ASSERT(rows * cols > 0);
		if (m_pData != NULL)
		{
			delete[] m_pData;	m_pData = NULL;	m_bAutoDelete = TRUE; m_maxCols = m_maxRows = m_maxAlloc = 0;
		}
		ReSize(rows, cols);
		memcpy(m_pData, source, sizeof(T) * m_maxRows * m_maxCols);
	}


	inline void SetValue(const T& val)
		{	_ASSERT(m_pData != NULL);
			for(int i = 0; i < m_maxRows * m_maxCols; i++)	m_pData[i] = val;
		}

	inline void SetVal(const T& val) { SetValue(val); }

	inline void operator *= (double s)
		{	_ASSERT(m_pData != NULL);
			for(int i = 0; i < m_maxRows * m_maxCols; i++)	m_pData[i] *= s;
		}

	inline void operator /= (double s)
		{	_ASSERT(m_pData != NULL);	_ASSERT(s != 0);
			for(int i = 0; i < m_maxRows * m_maxCols; i++)m_pData[i] /= s;
		}

	inline void operator += (const T& s)
		{	_ASSERT(m_pData != NULL);
			for(int i = 0; i < m_maxRows * m_maxCols; i++)	m_pData[i] += (T)s;
		}

	inline void GetTranspose(CMatrix<T>& m) const
		{	_ASSERT(m_pData != NULL);
			m.ReSize(m_maxCols, m_maxRows);
			for(int i = 0; i < m_maxRows; i++)
				for(int j = 0; j < m_maxCols; j++)
					m(j,i) = (*this)(i,j);
		}

	inline CMatrix<T> GetTranspose() const
		{	_ASSERT(m_pData != NULL);
			CMatrix<T> m;
			m.ReSize(m_maxCols, m_maxRows);
			for(int i = 0; i < m_maxRows; i++)
				for(int j = 0; j < m_maxCols; j++)
					m(j,i) = (*this)(i,j);
			return m;
		}

	inline void Identity()
	{
		for(int i = 0; i < m_maxRows; i++)
		{
			for(int j = 0; j < m_maxCols; j++)
				(*this)(i,j) = 0;
			(*this)(i,i) = 1;
		}
	}

  protected:
	BOOL m_bAutoDelete;
    int m_maxRows;
    int m_maxCols;
	int m_maxAlloc;
    T* m_pData;
};


template<class T>  //multiply M1(mxn) * M2(n*l) -> res(mxl)
inline void Mult(const CMatrix<T>& matr1, const CMatrix<T>& matr2, CMatrix<T>& res)
{
	int m = matr1.GetRows();
	int n = matr1.GetCols();
	ASSERT(n == matr2.GetRows());
	int l = matr2.GetCols();
	ASSERT(res.GetRows()==m);
	ASSERT(res.GetCols()==l);
//	res.ReSize(m,l);
	for(int i = 0; i< m; i++)
		for(int j = 0; j < l; j++)
		{
			T data = 0;
			for(int k = 0; k < n; k++)
				data += (matr1(i,k) * matr2(k,j));
			res(i,j) = data; 
		}
}

template<class T>  //multiply M(nxm) * scalar a
inline void Mult(CMatrix<T>& matr, T a)
{
	int m = matr.GetRows();
	int n = matr.GetCols();
	if(m*n>USE_OMP_AZA_LMT)
	{
#ifdef USE_OMP_AZA
	#pragma omp parallel for
#endif
		for(int i = 0; i< m; i++)
			for(int j = 0; j < n; j++)
			{
				matr(i,j) *= a;
			}
	}
	else
	{
		for(int i = 0; i< m; i++)
			for(int j = 0; j < n; j++)
			{
				matr(i,j) *= a;
			}
	}
}

template<class T>  //add M1(nxn) + M2(n*n) -> res(nxn)
inline void Add(const CMatrix<T>& matr1, const CMatrix<T>& matr2, CMatrix<T>& res)
{
	int n = matr1.GetRows();
	ASSERT(n == matr2.GetRows());
	ASSERT(n == matr1.GetCols());
	ASSERT(n == matr2.GetCols());
	ASSERT(n == res.GetCols());
	ASSERT(n == res.GetRows());
	if(n>10*USE_OMP_AZA_LMT)
	{
#ifdef USE_OMP_AZA
	#pragma omp parallel for
#endif
		for(int i = 0; i< n; i++)
			for(int j = 0; j < n; j++)
			{
				res(i,j) = matr1(i,j)+matr2(i,j);
			}
	}
	else
	{
		for(int i = 0; i< n; i++)
			for(int j = 0; j < n; j++)
			{
				res(i,j) = matr1(i,j)+matr2(i,j);
			}
	}
}

template<class T>  //multiply M(mxn) * V(n*1) ->vec(mx1)
inline CVector<T> Mult(const CMatrix<T>& matr, const CVector<T>& vec)
{
	int m = matr.GetRows();
	int n = matr.GetCols();
	ASSERT(m > 0);
	CVector<T> res(m);
	for(int i = 0; i < m; i++)
    {
		T data = 0;
		for(int j = 0; j < n; j++)
			data += (matr(i,j) * vec[j]);
		res[i] = data;
    }

	return res;
}

template<class T>  //multiply M(mxn) * V(n*1) ->vec(mx1)
inline void Mult(const CMatrix<T>& matr, const CVector<T>& vec, CVector<T>& res)
{
	int m = matr.GetRows();
	int n = matr.GetCols();
	ASSERT(vec.GetSize() == n);
	ASSERT(res.GetSize() == m);
//	res.ReSize(m);
	for(int i = 0; i < m; i++)
    {
		T data = 0;
		for(int j = 0; j < n; j++)
			data += (matr(i,j) * vec[j]);
		res[i] = data;
    }
}

template<class T>  //multiply M(mxn) * N(n*l) -> Mat(mxl)
inline CMatrix<T> Mult(const CMatrix<T>& mat1, const CMatrix<T>& mat2)
{
	int m = mat1.GetRows();
	int n = mat1.GetCols();
	int l = mat2.GetCols();
	ASSERT(mat2.GetRows() == n);
	CMatrix<T> matr(m, l);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < l; j++)
		{
			T data = T(0);
			for (int k = 0; k < n; k++)
			{
				data += mat1(i, k) * mat2(k, j);
			}
			matr(i, j) = data;
		}
	}
	return matr;
}

// define some types of matrices
typedef CMatrix<BOOL> CBOOLMatrix; // 0 to 255
typedef CMatrix<__int8> CInt8Matrix; // -128 to 127
typedef CMatrix<unsigned __int8> CUInt8Matrix; // 0 to 255
typedef CMatrix<__int16> CInt16Matrix;// –32,768 to 32,767
typedef CMatrix<unsigned __int16> CUInt16Matrix; // 0 to 65535
typedef CMatrix<int> CIntMatrix;// –2,147,483,648 to 2,147,483,647
typedef CMatrix<unsigned int> CUIntMatrix;//0 to 4,294,967,295
typedef CMatrix<long> CLongMatrix;
typedef CMatrix<float> CFloatMatrix;
typedef CMatrix<double> CDoubleMatrix;
typedef CMatrix<long double> CLDoubleMatrix;
typedef CMatrix<COLORREF> CColorMatrix;
typedef CMatrix<BYTE> CByteMatrix;

#endif // _ARRAYS_HPP_


