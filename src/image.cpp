#include "image.hpp"
#include "common.hpp"

#include <algorithm>

template <typename T, std::size_t n> Histgram<T, n>::Histgram()
{
    std::fill(rs.begin(), rs.end(), 0);
    std::fill(gs.begin(), gs.end(), 0);
    std::fill(bs.begin(), bs.end(), 0);
}

template <typename T>
Image<T>::Image(std::size_t height, std::size_t width)
    : height_(height)
    , width_(width)
{
    data_.resize(height * width);
}

template <typename T> T* Image<T>::data() { return data_.data(); }

template <typename T> const T* Image<T>::data() const { return data_.data(); }

template <typename T> T& Image<T>::pixel(std::size_t y, std::size_t x) { return data_[width_ * y + x]; }

template <typename T> const T& Image<T>::pixel(std::size_t y, std::size_t x) const { return data_[width_ * y + x]; }

// explicit template

template class Histgram<value_t>;

template class Image<Color<value_t>>;
