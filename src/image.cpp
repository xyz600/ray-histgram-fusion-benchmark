#include "image.hpp"
#include "common.hpp"

#include <algorithm>
#include <memory>

template <typename T> void Color<T>::operator+=(const Color<T>& src)
{
    r += src.r;
    g += src.g;
    b += src.b;
}

template <typename T> void Color<T>::operator/=(const T divisor)
{
    r /= divisor;
    g /= divisor;
    b /= divisor;
}

template <typename T> Color<T>::Color() { r = g = b = 0; }

template <typename T, std::size_t n> Histgram<T, n>::Histgram() { std::fill(bins.begin(), bins.end(), 0); }

template <typename T>
Image<T>::Image(std::size_t height, std::size_t width)
    : height_(height)
    , width_(width)
    , data_(new T[height * width])
{
}

template <typename T> T* Image<T>::data() { return data_.get(); }

template <typename T> const T* Image<T>::data() const { return data_.get(); }

template <typename T> T& Image<T>::pixel(std::size_t y, std::size_t x) { return data_[width_ * y + x]; }

template <typename T> const T& Image<T>::pixel(std::size_t y, std::size_t x) const { return data_[width_ * y + x]; }

template <typename T> std::size_t Image<T>::height() const { return height_; }

template <typename T> std::size_t Image<T>::width() const { return width_; };

// explicit template

template class Color<value_t>;

template class Histgram<value_t>;

template class Image<Color<value_t>>;

template class Image<Histgram<value_t>>;