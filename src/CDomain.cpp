#include "CDomain.hpp"

template <class T>
CDomain<T>::CDomain(int UID, CVector<3, int> size, CVector<3, int> origin_cell, CVector<3, T> length) :
		_UID(UID), _size(size), _origin_cell(origin_cell), _length(length)
{
}

template <class T>
CDomain<T>::CDomain(int UID, CVector<3, int> size) :
		_UID(UID), _size(size)
{
	_origin_cell = CVector<3, int>(0, 0, 0);
	_length = CVector<3, T>(0.05, 0.05, 0.05);
}

template <class T>
CDomain<T>::~CDomain()
{
}

template <class T>
CVector<3,int> CDomain<T>::getOrigin() const
{
	return _origin_cell;
}

template <class T>
CVector<3,int> CDomain<T>::getSize() const
{
	return _size;
}

template <class T>
int CDomain<T>::getUid() const
{
	return _UID;
}

template <class T>
CVector<3,T> CDomain<T>::getLength() const
{
	return _length;
}

template class CDomain<double>;
template class CDomain<float>;
