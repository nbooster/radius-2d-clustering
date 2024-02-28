#pragma once

#include <vector>
#include <iterator>
#include <cstddef>

/*
 * Each chunk now is 8 bytes long (from 24),
 * and the method can handle up to 2^32-2 ~ 4.3 billion input points
 * (if there is enough RAM of course).
 */

using Integer = std::uint_fast32_t;
constexpr auto maxValue = std::numeric_limits<Integer>::max();

class Chunk
{
	Integer startIndex_{ maxValue };
	Integer endIndex_{ maxValue };
	static inline std::vector<Integer>* indicesPtr_{ nullptr };

	struct Iterator
	{
		using iterator_category = std::forward_iterator_tag;
		using difference_type = Integer;
		using value_type = Integer;
		using pointer = value_type*;
		using reference = value_type&;

		explicit Iterator(const value_type& index) :
			index_{ index }
		{}

		const value_type& operator*() const noexcept
		{
			return index_;
		}

		Iterator& operator++() noexcept
		{
			index_ = (*indicesPtr_)[index_];
			return *this;
		}

		Iterator operator++(int) noexcept
		{
			Iterator tmp = *this;
			++(*this);
			return tmp;
		}

		friend bool operator== (const Iterator& a, const Iterator& b) = default;

	private:

		value_type index_;
	};

public:

	Chunk() = default;

	static void setIndices(std::vector<Integer>& indices) noexcept
	{
		indicesPtr_ = &indices;
	}

	bool isEmpty(void) const noexcept
	{
		return startIndex_ == maxValue;
	}

	void addIndex(const Integer index) noexcept
	{
		if (startIndex_ == maxValue)
		{
			startIndex_ = endIndex_ = index;
			return;
		}

		(*indicesPtr_)[endIndex_] = index;
		endIndex_ = index;
	}

	size_t size() const noexcept
	{
		size_t size = 0;

		for (auto currentIndex{ startIndex_ }; currentIndex != endIndex_; currentIndex = (*indicesPtr_)[currentIndex])
			++size;

		return size;
	}

	Iterator begin() const noexcept
	{
		return Iterator(startIndex_);
	}

	Iterator end() const noexcept
	{
		return Iterator(maxValue);
	}
};