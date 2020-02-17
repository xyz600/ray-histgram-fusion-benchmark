#include "common.hpp"
#include "image.hpp"

#include <cassert>
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
        kappa = 0.5;
        m = 7.5;
        s = 2;
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

    const auto half_padding = config.patch_size;

    for (index_t y = half_padding; y < config.height + half_padding; y++)
    {
        for (index_t x = half_padding; x < config.width + half_padding; x++)
        {
            auto& hist = color_histgram.pixel(y, x);
            auto& pixel = input.pixel(y, x);

            for (std::uint32_t n = 0; n < config.number_of_sample; n++)
            {
                Color<value_t> c;
                std::cin >> c.r >> c.g >> c.b;

                pixel += c;

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

    // padding 部分を埋める
    for (index_t y = 0; y < input.height(); y++)
    {
        for (index_t x = 0; x < input.width(); x++)
        {
            const auto cy = clamp(half_padding, y, config.height + half_padding);
            const auto cx = clamp(half_padding, x, config.width + half_padding);

            input.pixel(y, x) = input.pixel(cy, cx);
            color_histgram.pixel(y, x) = color_histgram.pixel(cy, cx);
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

void calculate_pixel_distance(
    value_t* pixel_distance_table, Image<Histgram<value_t>>& color_histgram, Config& config, const index_t max_range)
{
    const auto max_range_all = 2 * max_range + 1;
    const auto half_padding = config.patch_size;

    // ピクセル間の xi-squared-distance を計算
    for (index_t y = half_padding; y < config.height + half_padding; y++)
    {
        for (index_t x = half_padding; x < config.width + half_padding; x++)
        {
            const auto index = y * color_histgram.width() + x;

            for (index_t dy = 0; dy < max_range_all; dy++)
            {
                for (index_t dx = 0; dx < max_range_all; dx++)
                {
                    const auto d_index = dy * max_range_all + dx;

                    pixel_distance_table[index * max_range_all + d_index] = xi_squared_distance(
                        color_histgram.pixel(y, x), color_histgram.pixel(y + dy - max_range, x + dx - max_range));
                }
            }
        }
    }
}

value_t calculate_distance_patch(const value_t* pixel_distance_table, index_t y1, index_t x1, index_t y2, index_t x2,
    const index_t patch_size, const index_t max_range, const index_t stride, const index_t height)
{
    const auto max_range_all = 2 * max_range + 1;
    const auto dy = y2 - y1;
    const auto dx = x2 - x1;
    const auto d_index = (dy + max_range) * max_range_all + (dx + max_range);

    value_t distance_sum = 0.0;
    for (index_t p_dy = -patch_size; p_dy <= patch_size; p_dy++)
    {
        for (index_t p_dx = -patch_size; p_dx <= patch_size; p_dx++)
        {
            const auto p_index = (y1 + p_dy) * stride + x1 + p_dx;

            // TODO: d_index がループ内固定なので、gather すると持ってこれる
            distance_sum += pixel_distance_table[p_index * max_range_all * max_range_all + d_index];
        }
    }
    return distance_sum;
}

void add_patch(Color<value_t>* sub_buffer, index_t src_y, index_t src_x, index_t dst_y, index_t dst_x,
    Image<Color<value_t>>& input, const index_t patch_size)
{
    const index_t patch_range_all = patch_size * 2 + 1;

    for (index_t p_dy = -patch_size; p_dy <= patch_size; p_dy++)
    {
        for (index_t p_dx = -patch_size; p_dx <= patch_size; p_dx++)
        {
            const auto p_index = (p_dy + patch_size) * patch_range_all + (p_dx + patch_size);

            const auto dst_index = (dst_y + p_dy) * input.width() + (dst_x + p_dx);

            auto& dst = sub_buffer[dst_index * (patch_range_all * patch_range_all) + p_index];
            auto& src = input.pixel(src_y + p_dy, src_x + p_dx);
            dst += src;
        }
    }
}

void ray_histgram_fusion(Image<Color<value_t>>& input, Image<Histgram<value_t>>& color_histgram,
    Image<Color<value_t>>& output, Config& config)
{
    auto start1 = std::chrono::system_clock::now();

    // initialize output buffer
    for (index_t y = 0; y < output.height(); y++)
    {
        for (index_t x = 0; x < output.width(); x++)
        {
            auto& p = output.pixel(y, x);
            p.r = p.g = p.b = 0;
        }
    }

    const auto max_range = std::max<>(config.patch_size, config.search_range);
    const auto max_range_all = 2 * max_range + 1;

    const auto search_range_all = 2 * config.search_range + 1;
    const auto patch_range_all = 2 * config.patch_size + 1;
    const auto half_padding = config.patch_size;

    // pixel 間の xi-squared distance.
    // シンプルさと並列度を重視して、重複あり、自身との距離計算ありで実施
    std::unique_ptr<value_t[]> pixel_distance_table(
        new value_t[input.height() * input.width() * max_range_all * max_range_all]);

    // filter_count[y][x] = 中心が (y ,x) の patchを重ねる回数
    std::unique_ptr<value_t[]> patch_count(new value_t[input.height() * input.width()]);
    for (index_t i = 0; i < input.height() * input.width(); i++)
    {
        patch_count[i] = 0;
    }

    // sub_buffer[(y, x, p_dy, p_dx)] = y, x の diff を計算
    // ピクセル間で独立に操作したいので、中心座標毎にバッファを別に持つ
    std::unique_ptr<Color<value_t>[]> sub_buffer(
        new Color<value_t>[input.height() * input.width() * patch_range_all * patch_range_all]);

    auto end0 = std::chrono::system_clock::now();
    std::cerr << "elapsed time of initialization: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - start1).count() << "[ms]" << std::endl;

    calculate_pixel_distance(pixel_distance_table.get(), color_histgram, config, max_range);

    auto end1 = std::chrono::system_clock::now();
    std::cerr << "elapsed time of pixel distance: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - end0).count() << "[ms]" << std::endl;

    // patch 間距離を計算して、一定値以下なら output_buffer に足し込む
    for (index_t y = half_padding; y < config.height + half_padding; y++)
    {
        for (index_t x = half_padding; x < config.width + half_padding; x++)
        {
            const auto index = y * input.width() + x;

            const auto s_sy = std::max<>(half_padding, y - config.search_range);
            const auto s_ey = std::min<>(y + config.search_range, config.height - 1 + half_padding);
            const auto s_sx = std::max<>(half_padding, x - config.search_range);
            const auto s_ex = std::min<>(x + config.search_range, config.width - 1 + half_padding);

            // search range
            for (index_t s_y = s_sy; s_y <= s_ey; s_y++)
            {
                for (index_t s_x = s_sx; s_x <= s_ex; s_x++)
                {
                    const auto distance = calculate_distance_patch(pixel_distance_table.get(), y, x, s_y, s_x,
                        config.patch_size, max_range, color_histgram.width(), color_histgram.height());
                    if (distance < config.kappa)
                    {
                        patch_count[index]++;
                        add_patch(sub_buffer.get(), y, x, s_y, s_x, input, config.patch_size);
                    }
                }
            }
        }
    }

    auto end2 = std::chrono::system_clock::now();
    std::cerr << "elapsed time of patch distance: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end1).count() << "[ms]" << std::endl;

    // sub_buffer, patch_count から、最終的なバッファに足しこんでいく
    // patch 間距離を計算して、一定値以下なら output_buffer に足し込む
    for (index_t y = half_padding; y < config.height + half_padding; y++)
    {
        for (index_t x = half_padding; x < config.width + half_padding; x++)
        {
            const index_t index = y * input.width() + x;

            index_t count = 0;

            for (index_t p_dy = -config.patch_size; p_dy <= config.patch_size; p_dy++)
            {
                for (index_t p_dx = -config.patch_size; p_dx <= config.patch_size; p_dx++)
                {
                    const auto ny = y + p_dy;
                    const auto nx = x + p_dx;
                    const auto n_index = ny * input.width() + nx;
                    // 周辺 pixel の (y, x) 相当の場所を見て確認する
                    const auto p_index = (-p_dy + config.patch_size) * patch_range_all - p_dx + config.patch_size;

                    output.pixel(y, x) += sub_buffer[n_index * patch_range_all + p_index];
                    count += patch_count[n_index];
                }
            }
            output.pixel(y, x) /= count;
        }
    }

    auto end3 = std::chrono::system_clock::now();
    std::cerr << "elapsed time of accumulate buffer: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - end2).count() << "[ms]" << std::endl;
}

int main()
{
    Config config;
    config.setup();

    const auto padding = config.patch_size * 2;

    Image<Color<value_t>> input(config.height + padding, config.width + padding);
    Image<Histgram<value_t>> color_histgram(config.height + padding, config.width + padding);
    Image<Color<value_t>> output(config.height + padding, config.width + padding);

    setup(input, color_histgram, config);

    const auto start = std::chrono::system_clock::now();

    ray_histgram_fusion(input, color_histgram, output, config);

    const auto end = std::chrono::system_clock::now();

    std::cout << "elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "[ms]" << std::endl;

    return 0;
}