#pragma once
#include <assert.h>
// A container support random access, not thread safe,
// with default pre allocated memory, in stack
// "P" in CPArray means pre allocate
// designed to contain basic numerical types, like int, double, pointer
// not safe for types that require destructor

template<typename T, int preAllocateSize >
class CPArray {
public:
	CPArray() : 
		pMem(pPreAllocated)
	{};

	// deep copy constructor
	CPArray(const CPArray & arr) :
		pMem(pPreAllocated)
	{
		reserve(arr.capacity());
		mSize = 0;
		for (int i = 0; i < arr.size(); i++)
		{
			pMem[mSize] = arr[i];
			++mSize;;
		}
	};

	~CPArray() {
		if (pMem != pPreAllocated) {
			delete[] pMem;
		}
	}

	T& operator[](const int & i) {
		assert(i < mSize);
		return pMem[i];
	}

	const T& operator[](const int& i) const {
		assert(i < mSize);
		return pMem[i];
	}

	void push_back(const T & newMember) {
		if (mSize + 1 > mCapacity) {
			reserve(mCapacity + (int)(mCapacity / 2) + 3);
		}

		pMem[mSize] = newMember;
		++mSize;
	}

	void push_back(T && newMember) {
		if (mSize + 1 > mCapacity) {
			reserve(mCapacity + (int)(mCapacity / 2) + 3);
		}

		pMem[mSize] = std::move(newMember);
		++mSize;
	}

	// reserve memory, mention that CPArray cannot reserve memory with size less than preAllocateSize
	void reserve(const size_t & newCap) {
		if (mSize > newCap) {
			return;
		}
		if (newCap <= preAllocateSize) {
			if (pPreAllocated == pMem) {
				return;
			}
			else{
				std::move(pMem, pMem + mSize, pPreAllocated);
				delete[] pMem;
				pMem = pPreAllocated;
				mCapacity = preAllocateSize;
			}
		}
		else{
			T * newMem = new T[newCap];
			std::move(pMem, pMem + mSize, newMem);
			mCapacity = newCap;
			if (pMem != pPreAllocated) {
				delete[] pMem;
			}
			pMem = newMem;
		}

	}

	void clear() {
		mSize = 0;
	}

	void erase(size_t i) {
		for (size_t j = i; int(j) < int(mSize) - 1; ++j){
			pMem[j] = std::move(pMem[j+1]);
		}
		--mSize;
	}

	//erase form index i to i+n-1
	void eraseN(size_t i, size_t n) {
		for (size_t j = i; int(j) < int(mSize) - n; ++j) {
			pMem[j] = std::move(pMem[j + n]);
		}
		mSize -= n;
	}

	void insert(size_t i, const T & memberToInsert) {
		if (mSize + 1 > mCapacity) {
			reserve(mCapacity + (int)(mCapacity / 2) + 3);
		}
		for (size_t j = mSize; j > i + 1; --j) {
			pMem[j] = std::move(pMem[j - 1]);
		}
		pMem[i+1] = memberToInsert;
		++mSize;
	}
	void insert(size_t i, T && memberToInsert) {
		if (mSize + 1 > mCapacity) {
			reserve(mCapacity + (int)(mCapacity / 2) + 3);
		}
		for (size_t j = mSize; j > i + 1; --j) {
			pMem[j] = std::move(pMem[j-1]);
		}
		pMem[i+1] = std::move(memberToInsert);
		++mSize;
	}

	//insert n elemets after element with index i, only make rooms for those elements,
	//need to be initailized by you afterwards.
	void insertN(size_t i, int n) {	
		if (mSize + n > mCapacity) {
			reserve(mCapacity + (int)(mCapacity / 2) + n + 3);
		}
		for (size_t j = mSize + n - 1; j > n + i; --j) {
			pMem[j] = std::move(pMem[j - n]);
		}
		mSize += n;
	}

	T* begin() {
		return pMem;
	}
	T* end() {
		return pMem + mSize;
	}

	T& front() {
		return *pMem;
	}
	T& back() {
		return pMem[mSize - 1];
	}
	void pop_back() {
		--mSize;
	}

	bool has(const T & t) {
		for (size_t j = 0; j < mSize; ++j) {
			if (t == pMem[j]) {
				return true;
			}
		}
		return false;
	}

	size_t size() const { return mSize; }
	size_t capacity() const { return mCapacity; }
private:
	T pPreAllocated[preAllocateSize];
	T * pMem;
	size_t mSize = 0;
	size_t mCapacity = preAllocateSize;
};

// this Array does not allow dynamic memory allocation
template<typename T, int preAllocateSize >
class CPArrayStatic {
public:
	CPArrayStatic() 
	{};

	~CPArrayStatic() {
	}

	T& operator[](const int& i) {
		assert(i < mSize);
		return pPreAllocated[i];
	}

	bool push_back(const T& newMember) {
		if (mSize + 1 > mCapacity) {
			return false;
		}

		pPreAllocated[mSize] = newMember;
		++mSize;
		return true;
	}

	bool push_back(T&& newMember) {
		if (mSize + 1 > mCapacity) {
			false;
		}

		pPreAllocated[mSize] = std::move(newMember);
		++mSize;
		return true;
	}

	void clear() {
		mSize = 0;
	}

	void erase(size_t i) {
		for (size_t j = i; int(j) < int(mSize) - 1; ++j) {
			pPreAllocated[j] = std::move(pPreAllocated[j + 1]);
		}
		--mSize;
	}

	//erase form index i to i+n-1
	void eraseN(size_t i, size_t n) {
		for (size_t j = i; int(j) < int(mSize) - n; ++j) {
			pPreAllocated[j] = std::move(pPreAllocated[j + n]);
		}
		mSize -= n;
	}

	bool insert(size_t i, const T& memberToInsert) {
		if (mSize + 1 > mCapacity) {
			return false;
		}
		for (size_t j = mSize; j > i + 1; --j) {
			pPreAllocated[j] = std::move(pPreAllocated[j - 1]);
		}
		pPreAllocated[i + 1] = memberToInsert;
		++mSize;
		return true;
	}

	bool insert(size_t i, T&& memberToInsert) {
		if (mSize + 1 > mCapacity) {
			return false;
		}
		for (size_t j = mSize; j > i + 1; --j) {
			pPreAllocated[j] = std::move(pPreAllocated[j - 1]);
		}
		pPreAllocated[i + 1] = std::move(memberToInsert);
		++mSize;
		return true;
	}

	//insert n elemets after element with index i, only make rooms for those elements,
	//need to be initailized by you afterwards.
	bool insertN(size_t i, int n) {
		if (mSize + n > mCapacity) {
			return false;
		}
		for (size_t j = mSize + n - 1; j > n + i; --j) {
			pPreAllocated[j] = std::move(pPreAllocated[j - n]);
		}
		mSize += n;
		return true;
	}

	T* begin() {
		return pPreAllocated;
	}
	T* end() {
		return pPreAllocated + mSize;
	}

	T& front() {
		return *pPreAllocated;
	}
	T& back() {
		return pPreAllocated[mSize - 1];
	}
	T pop_back() {
		--mSize;
		return pPreAllocated[mSize];
	}

	bool empty() {
		return mSize == 0;
	}

	bool has(const T& t) {
		for (int j = 0; j < mSize; ++j) {
			if (t == pPreAllocated[j]) {
				return true;
			}
		}
		return false;
	}

	size_t size() { return mSize; }
	size_t capacity() { return mCapacity; }
private:
	alignas(16) T pPreAllocated[preAllocateSize];
	size_t mSize = 0;
	size_t mCapacity = preAllocateSize;
};

template <typename T, int maxSize>
class CircularArray
{
public:
	CircularArray() {};
	~CircularArray() {};

	int size() const { return mSize; }

	void clear(bool destruct = false) {
		while (mSize != 0)
		{
			popBack(destruct);
		}
	}

	bool has(const T & d) {
		for (int i = 0; i < size(); i++)
		{
			int actualIndex = i + mBeginIdx;
			if (actualIndex >= maxSize)
			{
				actualIndex -= maxSize;
			}
			if (d == mData[actualIndex])
			{
				return true;
			}
		}
		return false;
	}

	bool hasBackward(const T& d) {
		for (int i = size()-1; i >= 0 ; i--)
		{
			int actualIndex = i + mBeginIdx;
			if (actualIndex >= maxSize)
			{
				actualIndex -= maxSize;
			}
			if (d == mData[actualIndex])
			{
				return true;
			}
		}
		return false;
	}

	T& front() {
		if (!mSize)
		{
			assert(false && "Array is empty! ");
		}
		return mData[mBeginIdx];
	}

	T& back() {
		if (!mSize)
		{
			assert(false && "Array is empty! ");
		}
		return mData[mEndIdx];
	}

	void push_back(const T& newElement) {
		if (mSize == 0)
			// array is empty, in this case mBeginIdx == mEndIdx
		{
			mData[mEndIdx] = newElement;
			++mSize;
		}
		else
			// array is not empty
		{
			moveIdxForward(mEndIdx);
			if (mEndIdx != mBeginIdx)
				// end index hasn't caught up with begin index, increase the size
			{
				mData[mEndIdx] = newElement;
				++mSize;
			}
			else
				// end index has caught up with begin index, the array is full; overwrite the oldest element
				// the size won't increase in this case
			{
				moveIdxForward(mBeginIdx);
				mData[mEndIdx] = newElement;
			}
		}


	}

	// set destruct = true to explictly destruct the element
	void popFront(bool destruct = false) {

		if (mSize != 0)
		{
			--mSize;

			int mBeginIdxOld = mBeginIdx;
			// if its already empty, stay at where it is to make sure mEndIdx == mBeginIds
			if (mSize != 0) {
				moveIdxForward(mBeginIdx);
			}
			if (destruct)
			{
				mData[mBeginIdxOld].~T();
			}
			else {
			}
		}
		else
		{
			assert(false && "Trying to pop an empty array! ");
		}
	}

	// set destruct = true to explictly destruct the element
	void popBack(bool destruct = false) {

		if (mSize != 0)
		{
			--mSize;
			int mEndIdxOld = mEndIdx;
			// if its already empty, stay at where it is to make sure mEndIdx == mBeginIds
			if (mSize != 0) {
				moveIdxBackward(mEndIdx);
			}
			if (destruct)
			{
				mData[mEndIdxOld].~T();
			}
			else {
			}
		}
		else
		{
			assert(false && "Trying to pop an empty array! ");
		}
	}

	// count from the oldest element to newest element
	// 0th is the oldest, (size-1)th is the newest element
	// it won't detect the if it has exceeded the boundary
	// you will have to make sure idx < size()
	T& operator [] (int idx) {
		int actualIndex = idx + mBeginIdx;
		if (actualIndex >= maxSize)
		{
			actualIndex -= maxSize;
		}

		return mData[actualIndex];
	}

	const T& operator []  (int idx) const {
		int actualIndex = idx + mBeginIdx;
		if (actualIndex >= maxSize)
		{
			actualIndex -= maxSize;
		}

		return mData[actualIndex];
	}



protected:
	inline void moveIdxForward(int& idx) {
		++idx;
		if (idx >= maxSize) {
			idx = idx - maxSize;
		}
	}

	inline void moveIdxBackward(int& idx) {
		--idx;
		if (idx < 0) {
			idx = maxSize - 1;
		}
	}

	int mSize = 0;

	// if mBeginIdx == mEndIdx, the array can either be empty or has one element; in this case we need to rely on mSize to determine
	int mBeginIdx = 0; // the index of the latest element
	int mEndIdx = 0;   // the index of the newest element

	alignas(16) T mData[maxSize];
};

template <typename T, int maxSize>
inline std::ostream& operator<<(std::ostream& os, const CircularArray<T, maxSize>& arr)
{
	for (size_t i = 0; i < arr.size(); i++)
	{
		os << " " << arr[i];
	};
	return os;
}
