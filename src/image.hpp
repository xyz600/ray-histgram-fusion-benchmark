#pragma once

#include <array>
#include <memory>

template <typename T> struct Color
{
    using value_type = T;

    T r;
    T g;
    T b;

    Color();

    void operator+=(const Color& src);
    void operator/=(const T divisor);
};

template <typename T, std::size_t n = 20> struct Histgram
{
    using value_type = T;
    static constexpr std::size_t BinSizePerColor = n;
    static constexpr std::size_t BinSize = BinSizePerColor * 3;

    // r, g, b のヒストグラムが concat されている
    std::array<T, n * 3> bins;

    std::uint32_t nonzero_sample;

    Histgram();
};

template <typename T> class Image
{
public:
    Image(std::size_t height, std::size_t width);

    T* data();
    const T* data() const;

    T& pixel(std::size_t y, std::size_t x);
    const T& pixel(std::size_t y, std::size_t x) const;

    std::size_t height() const;

    std::size_t width() const;

private:
    std::unique_ptr<T[]> data_;

    std::size_t height_;

    std::size_t width_;
};