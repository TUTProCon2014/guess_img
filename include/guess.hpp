#pragma once

#include "../../inout/include/inout.hpp"
#include "../../utils/include/template.hpp"
#include "../../utils/include/types.hpp"
#include "../../utils/include/range.hpp"

#include <vector>
#include <set>
#include <array>
#include <deque>
#include <unordered_set>

namespace procon { namespace guess {

using namespace utils;

/**
Predicate fは、f(image1, image2, Direction::up) -> double を返す
doubleの`絶対値の値が小さい方`を優先します。

こんなふうに使う。g++4.9だとgeneric lambdaが使えるけど、g++4.8とかだと使えないしつらい。
Example:
------------
auto idxRC = guess(problem, [](Image const & p1,
                               Image const & p2,
                               Direction dir)
                               {
                                    // ... 比較関数
                                    // p1のdir方向にp2がどれだけ結合しやすいかをdoubleで返す。
                                    // 返す値は、別に[0, 1)じゃなくてもいい。
                                    // 負の数を返しても良いし、infinityを返しても良い
                               });
------------
*/
template <typename BinFunc>
std::vector<std::vector<ImageID>> guess(Problem const & problem, BinFunc f)
{
    auto remain = [&](){
        std::unordered_set<ImageID> dst;
        for(auto r : iota(problem.div_y()))
            for(auto c : iota(problem.div_x())){
                dst.insert(ImageID(r, c));
            }

        dst.erase(ImageID(0, 0));  // (0, 0)は原点として最初から使う
        return dst;
    }();


    // 画像(r, c) = originを中心にして、縦方向か横方向に画像を結合していく
    auto guess_oneline = [&](ImageID origin, bool isVerticalLine) -> std::vector<ImageID>
    {
        std::deque<ImageID> dst;
        dst.push_back(origin);      // 最初に原点がある
        for(auto t : iota(isVerticalLine ? problem.div_y()-1 : problem.div_x()-1)){
            Direction dir;
            ImageID mIdx;
            double min = std::numeric_limits<double>::infinity();

            // std::array<std::size_t, 2> ds = {0, dst.size()-1};    // {0 : 上(左), dst.size()-1 : 下(右)}
            for(auto dI : iota(2)){               // 結合した画像集合の上か下にくっつくはず
                Direction d = isVerticalLine
                    ? (dI == 0 ? Direction::up : Direction::down)
                    : (dI == 0 ? Direction::left : Direction::right);

                std::size_t tgtIdx = 0;
                if(dI == 1)
                    tgtIdx = dst.size() - 1;

                for(auto& idx : remain){    // 残っている画像の中から探す
                    const double v = std::abs(f(problem.get_element(dst[tgtIdx]),
                                                problem.get_element(idx),
                                                d));
                    if(min >= v){   // min == v == infのときは入れ替える
                        min = v;
                        dir = d;
                        mIdx = idx;
                    }
                }
            }
            if(dir == (isVerticalLine ? Direction::up : Direction::left))    // 先頭にくっつける
                dst.push_front(mIdx);
            else                        // 後ろにくっつける
                dst.push_back(mIdx);

            remain.erase(mIdx);         // 決定したので消す
        }

        std::vector<ImageID> vec(dst.begin(), dst.end());
        return vec;
    };


    // (0, 0)の画像に対して、まずは縦方向に結合し、
    // その後、横方向に結合していく
    // イメージ的には、
    // 1.      ↑                    (R1, C1)
    //      (0, 0)           =>     (0 , 0 )
    //         ↓                    (R2, C2)
    //
    // 2. この操作は 1. で求めた各要素に対して行う
    //    ← (R1, C1) →       => (R3, C3), (R1, C1), (R4, C4)
    //
    //
    std::vector<std::vector<ImageID>> dst;
    for(auto& o : guess_oneline(ImageID(0, 0), true))
        dst.push_back(guess_oneline(o, false));

    return dst;
}


/**
ある画像img1の方角directionに対して、画像img2がどの程度相関があるかを返します。
相関があるほど返す値は絶対値が小さくなります。
また、返す値は必ず正です。
*/
template <typename T, typename U
#ifdef SUPPORT_TEMPLATE_CONSTRAINTS
    , PROCON_TEMPLATE_CONSTRAINTS(is_image<T>() && is_image<U>())    // T, Uともに画像であるという制約
#endif
>
double diff_connection(T const & img1, U const & img2, Direction direction)
{
    double sum = 0;

    if(img1.height() != img2.height() || img1.width() != img2.width())
        return std::numeric_limits<double>::infinity();

    std::size_t r1 = 0, c1 = 0, r2 = 0, c2 = 0;
    switch(direction)
    {
      case Direction::right:
        c1 = img1.width() - 1;
        c2 = 0;
        break;

      case Direction::up:
        r1 = 0;
        r2 = img2.height() - 1;
        break;

      case Direction::left:
        c1 = 0;
        c2 =  img2.width() - 1;
        break;

      case Direction::down:
        r1 = img1.height() - 1;
        r2 = 0;
        break;
    }

    switch(direction)
    {
      case Direction::right:
      case Direction::left:
        for(std::size_t r = 0; r < img1.height(); ++r){
            auto p1 = img1.get_pixel(r, c1).vec();
            auto p2 = img2.get_pixel(r, c2).vec();
            for(std::size_t i = 0; i < 3; ++i)
                sum += std::abs(static_cast<float>(p1[i]) - static_cast<float>(p2[i]));
        }
        sum /= img1.height();
        break;

      case Direction::up:
      case Direction::down:
        for(std::size_t c = 0; c < img1.width(); ++c){
            auto p1 = img1.get_pixel(r1, c).vec();
            auto p2 = img2.get_pixel(r2, c).vec();
            for(std::size_t i = 0; i < 3; ++i)
                sum += std::abs(static_cast<float>(p1[i]) - static_cast<float>(p2[i]));
        }
        sum /= img1.width();
        break;
    }

    return sum;
}


// PROCON_DEF_STRUCT_FUNCTION(DiffConnection, diff_connection);

}} // namespace procon::guess