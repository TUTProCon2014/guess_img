#pragma once

#include "../../inout/include/inout.hpp"
#include "../../utils/include/template.hpp"
#include "../../utils/include/types.hpp"
#include "../../utils/include/range.hpp"

#include <vector>
#include <set>
#include <array>
#include <deque>

namespace procon { namespace guess {

using namespace utils;

/**
Predicate fは、f(image1, image2, Direction::up) -> double を返す
doubleの`絶対値の値が小さい方`を優先します。

こんなふうに使う。g++4.9だとgeneric lambdaが使えるけど、g++4.8とかだと使えないしつらい。
Example:
------------
auto idxRC = guess(problem, [](utils::ElementImage const & p1,
                               utils::ElementImage const & p2,
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
std::vector<std::vector<Index2D>> guess(utils::Problem const & problem, BinFunc f)
{
    auto remain = [&](){
        std::set<Index2D> dst;
        for(auto r : utils::iota(problem.div_y()))
            for(auto c : utils::iota(problem.div_x())){
                dst.insert(makeIndex2D(r, c));
            }

        Index2D idx; idx[0] = 0; idx[1] = 0;
        dst.erase(idx);  // (0, 0)は原点として最初から使う
        return dst;
    }();


    // 画像(r, c) = originを中心にして、縦方向か横方向に画像を結合していく
    auto guess_oneline = [&](Index2D origin, bool isVerticalLine){
        std::deque<Index2D> dst;
        dst.push_back(origin);      // 最初に原点がある
        for(auto t : utils::iota(isVerticalLine ? problem.div_y()-1 : problem.div_x()-1)){
            utils::Direction dir;
            Index2D mIdx;
            double min = std::numeric_limits<double>::infinity();

            // std::array<std::size_t, 2> ds = {0, dst.size()-1};    // {0 : 上(左), dst.size()-1 : 下(右)}
            for(auto dI : utils::iota(2)){               // 結合した画像集合の上か下にくっつくはず
                utils::Direction d = isVerticalLine
                    ? (dI == 0 ? utils::Direction::up : utils::Direction::down)
                    : (dI == 0 ? utils::Direction::left : utils::Direction::right);

                std::size_t tgtIdx = 0;
                if(dI == 1)
                    tgtIdx = dst.size() - 1;

                for(auto& idx : remain){    // 残っている画像の中から探す
                    const double v = std::abs(f(problem.get_element(dst[tgtIdx][0], dst[tgtIdx][1]),
                                                problem.get_element(idx[0], idx[1]),
                                                d));
                    if(min >= v){   // min == v == infのときは入れ替える
                        min = v;
                        dir = d;
                        mIdx = idx;
                    }
                }
            }
            if(dir == (isVerticalLine ? utils::Direction::up : utils::Direction::left))    // 先頭にくっつける
                dst.push_front(mIdx);
            else                        // 後ろにくっつける
                dst.push_back(mIdx);

            remain.erase(mIdx);         // 決定したので消す
        }

        std::vector<Index2D> vec(dst.begin(), dst.end());
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
    std::vector<std::vector<Index2D>> dst;
    for(auto& o : guess_oneline(makeIndex2D(0, 0), true))
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
    , PROCON_TEMPLATE_CONSTRAINTS(utils::is_image<T>() && utils::is_image<U>())    // T, Uともに画像であるという制約
#endif
>
double diff_connection(T const & img1, U const & img2, utils::Direction direction)
{
    double sum = 0;

    if(img1.height() != img2.height() || img1.width() != img2.width())
        return std::numeric_limits<double>::infinity();

    std::size_t r1 = 0, c1 = 0, r2 = 0, c2 = 0;
    switch(direction)
    {
      case utils::Direction::right:
        c1 = img1.width() - 1;
        c2 = 0;
        break;

      case utils::Direction::up:
        r1 = 0;
        r2 = img2.height() - 1;
        break;

      case utils::Direction::left:
        c1 = 0;
        c2 =  img2.width() - 1;
        break;

      case utils::Direction::down:
        r1 = img1.height() - 1;
        r2 = 0;
        break;
    }

    switch(direction)
    {
      case utils::Direction::right:
      case utils::Direction::left:
        for(std::size_t r = 0; r < img1.height(); ++r){
            auto p1 = img1.get_pixel(r, c1).vec();
            auto p2 = img2.get_pixel(r, c2).vec();
            for(std::size_t i = 0; i < 3; ++i)
                sum += std::abs(static_cast<float>(p1[i]) - static_cast<float>(p2[i]));
        }
        sum /= img1.height();
        break;

      case utils::Direction::up:
      case utils::Direction::down:
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