#include "common.hpp"
#include "image.hpp"

#include <chrono>
#include <cmath>
#include <iostream>

struct Config
{
    // 画像の高さ
    index_t height;

    // 画像の幅
    index_t width;

    // デノイズの判定最小単位(1辺 2 * patch_size + 1)
    index_t patch_size;

    // fusion する patch を探索する範囲(1辺 2 * search_range + 1)
    index_t search_range;

    // path tracing の pixel あたりサンプル数
    std::uint32_t number_of_sample;

    // patch 間の距離の threashold
    value_t kappa;

    // ヒストグラムの基準になる輝度。
    // ヒストグラムの各 bin のしきい値は以下の通り
    // [m / k, m * 2 / k, m * 3 / k, ..., m, s * m]
    value_t m;

    // ヒストグラムの最終 bin に含む係数
    value_t s;

    void setup()
    {
        std::cin >> height >> width >> number_of_sample;
        patch_size = 1;
        search_range = 5;
        kappa = 2.0;
    }
};

index_t calculate_histgram_index(value_t v, value_t m, value_t s, index_t bin_size)
{
    // threashold
    // [m * 1.0 / (bin_size - 1), m * 2.0 / (bin_size - 1), ..., m * (bin_size - 2) / (bin_size - 1), s * m]
    // x (bin_size - 1) -> [m, 2m, 3m, ..., (bin_size - 2)m, (bin_size - 1) * sm]
    // / m -> [1, 2, 3, ..., bin_size - 2, (bin_size - 1) * s]
    value_t normalized_v = std::pow(v, 1.0 / 2.2) * (bin_size - 1) / m;

    if (normalized_v < (bin_size - 2))
    {
        return static_cast<index_t>(normalized_v);
    }
    else if (normalized_v < (bin_size - 1) * s)
    {
        return bin_size - 2;
    }
    else
    {
        return bin_size - 1;
    }
}

void setup(Image<Color<value_t>>& input, Image<Histgram<value_t>>& color_histgram, const Config& config)
{
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);

    auto clamp = [](value_t lower, value_t v, value_t upper) { return std::min<>(std::max<>(lower, v), upper); };

    for (index_t y = 0; y < input.height(); y++)
    {
        for (index_t x = 0; x < input.width(); x++)
        {
            auto& hist = color_histgram.pixel(y, x);
            auto& pixel = input.pixel(y, x);

            for (std::uint32_t n = 0; n < config.number_of_sample; n++)
            {
                Color<value_t> c;
                std::cin >> c.r >> c.g >> c.b;

                pixel.r += c.r;
                pixel.g += c.g;
                pixel.b += c.b;

                c.r = clamp(0.0f, c.r, config.m * config.s);
                c.g = clamp(0.0f, c.g, config.m * config.s);
                c.b = clamp(0.0f, c.b, config.m * config.s);

                {
                    const auto index_r
                        = calculate_histgram_index(c.r, config.m, config.s, Histgram<value_t>::BinSizePerColor);
                    hist.bins[index_r]++;
                    if (index_r != 0)
                    {
                        hist.nonzero_sample++;
                    }
                }
                {
                    const auto index_g
                        = calculate_histgram_index(c.g, config.m, config.s, Histgram<value_t>::BinSizePerColor);
                    hist.bins[index_g + Histgram<value_t>::BinSizePerColor]++;
                    if (index_g != 0)
                    {
                        hist.nonzero_sample++;
                    }
                }
                {
                    const auto index_b
                        = calculate_histgram_index(c.b, config.m, config.s, Histgram<value_t>::BinSizePerColor);
                    hist.bins[index_b + 2 * Histgram<value_t>::BinSizePerColor]++;
                    if (index_b != 0)
                    {
                        hist.nonzero_sample++;
                    }
                }
            }
        }
    }
}

value_t xi_squared_distance(Histgram<value_t>& h1, Histgram<value_t>& h2)
{
    const auto kxy = h1.nonzero_sample + h2.nonzero_sample;
    value_t sum_all = 0;
    for (index_t i = 0; i < Histgram<value_t>::BinSize; i++)
    {
        const value_t sum = h1.bins[i] + h2.bins[i];
        if (sum > 0)
        {
            const value_t diff = h1.bins[i] - h2.bins[i];
            sum_all += (diff * diff) / sum;
        }
    }
    return sum_all / kxy;
}

void ray_histgram_fusion(Image<Color<value_t>>& input, Image<Histgram<value_t>>& color_histgram,
    Image<Color<value_t>>& output, Config& config)
{
    // initialize output buffer
    for (index_t y = 0; y < config.height; y++)
    {
        for (index_t x = 0; x < config.width; x++)
        {
            output.pixel(y, x).r = input.pixel(y, x).r;
            output.pixel(y, x).g = input.pixel(y, x).g;
            output.pixel(y, x).b = input.pixel(y, x).b;
        }
    }

    const auto max_range = std::max<>(config.patch_size, config.search_range);
    const auto max_range_all = 2 * max_range + 1;

    // pixel 間の xi-squared distance.
    // シンプルさと並列度を重視して、重複あり、自身との距離計算ありで実施
    std::unique_ptr<value_t[]> pixel_distance_table(
        new value_t[config.height * config.width * max_range_all * max_range_all]);

    // ピクセル間の xi-squared-distance を計算
    for (index_t y = 0; y < config.height; y++)
    {
        for (index_t x = 0; x < config.width; x++)
        {
            const auto index = y * config.width + x;

            for (index_t dy = -max_range; dy <= max_range; dy++)
            {
                for (index_t dx = -max_range; dx <= max_range; dx++)
                {
                    const auto d_index = (dy + max_range) * max_range_all + dx + max_range;
                    // pixel_distance_table[toIndex(y, x)][toIndex(dy, dx)] = distance(y, x, y + dy, x + dx);
                }
            }
        }
    }

    // filter_count[y][x] = 中心が (y ,x) の patchを重ねる回数
    std::unique_ptr<value_t[]> patch_count(new value_t[config.height * config.width]);

    const auto search_range_all = 2 * config.search_range + 1;
    const auto patch_range_all = 2 * config.patch_size + 1;

    // sub_buffer[(y, x, p_dy, p_dx)] = y, x の diff を計算
    // ピクセル間で独立に操作したいので、中心座標毎にバッファを別に持つ
    std::unique_ptr<Color<value_t>[]> sub_buffer(
        new Color<value_t>[config.height * config.width * patch_range_all * patch_range_all]);

    // patch 間距離を計算して、一定値以下なら output_buffer に足し込む
    for (index_t y = 0; y < config.height; y++)
    {
        for (index_t x = 0; x < config.width; x++)
        {
            const auto index = y * config.width + x;

            for (index_t s_dy = -config.search_range; s_dy <= config.search_range; s_dy++)
            {
                for (index_t s_dx = -config.search_range; s_dx <= config.search_range; s_dx++)
                {
                    // distance = distance(P(y, x), P(y + s_dy, x + s_dx))
                    // if (distance < kappa) {
                    // sub_buffer[toIndex(y, x)] += Patch(y + s_dy, x + s_dx)
                    // }

                    const auto c_y = y + s_dy;
                    const auto c_x = x + s_dx;
                    const auto c_index = y * config.width + x;

                    const auto s_index = (s_dy + config.search_range) * search_range_all + (s_dx + config.search_range);

                    value_t distance_sum = 0.0;
                    for (index_t p_dy = -config.patch_size; p_dy <= config.patch_size; p_dy++)
                    {
                        for (index_t p_dx = -config.patch_size; p_dx <= config.patch_size; p_dx++)
                        {
                            distance_sum += pixel_distance_table[c_index * max_range_all + s_index];
                        }
                    }
                    if (distance_sum <= config.kappa)
                    {
                        patch_count[index]++;
                        // Patch(y, x) += Patch(y + s_dy, x + s_dx);
                        for (index_t p_dy = -config.patch_size; p_dy <= config.patch_size; p_dy++)
                        {
                            for (index_t p_dx = -config.patch_size; p_dx <= config.patch_size; p_dx++)
                            {
                                const auto p_index
                                    = (p_dy + config.patch_size) * patch_range_all + (p_dx + config.patch_size);

                                const auto dst_index = (y + s_dy + p_dy) * config.width + (x + s_dx + p_dx);

                                auto& dst = sub_buffer[dst_index * (patch_range_all * patch_range_all) + p_index];
                                auto& src = input.pixel(y + s_dy + p_dy, x + s_dx + p_dx);
                                dst.r += src.r;
                                dst.g += src.g;
                                dst.b += src.b;
                            }
                        }
                    }
                }
            }
        }
    }

    // sub_buffer, patch_count から、最終的なバッファに足しこんでいく
    // patch 間距離を計算して、一定値以下なら output_buffer に足し込む
    for (index_t y = 0; y < config.height; y++)
    {
        for (index_t x = 0; x < config.width; x++)
        {
            for (index_t p_dy = -config.patch_size; p_dy <= config.patch_size; p_dy++)
            {
                for (index_t p_dx = -config.patch_size; p_dx <= config.patch_size; p_dx++)
                {
                    // output(y, x) += sub_buffer(y, x, p_dy, p_dx);
                    // count(y, x) += patch_count(y, x, p_dy, p_dx);
                }
            }
            // output(y, x) /= count(y, x)
        }
    }
}

int main()
{
    Config config;
    config.setup();

    Image<Color<value_t>> input(config.height, config.width);
    Image<Histgram<value_t>> color_histgram(config.height, config.width);
    Image<Color<value_t>> output(config.height, config.width);

    setup(input, color_histgram, config);

    const auto start = std::chrono::system_clock::now();

    ray_histgram_fusion(input, color_histgram, output, config);

    const auto end = std::chrono::system_clock::now();

    std::cout << "elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "[ms]" << std::endl;
}