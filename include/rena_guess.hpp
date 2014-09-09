/*
 使い方
 1.test.cppにこのファイルをインクルードして、rena_guess::rena_guess()を呼び出す
 2.評価関数をrena_guess::diff_connection_rena()に変える

 レナお姉さん2に特化したプログラム
*/

#pragma once

#include "../../inout/include/inout.hpp"
#include "../../utils/include/template.hpp"
#include "../../utils/include/types.hpp"
#include "../../utils/include/range.hpp"
#include "../../utils/include/exception.hpp"

#include <vector>
#include <set>
#include <array>
#include <deque>

namespace procon { namespace rena_guess {

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
std::vector<std::vector<Index2D>> rena_guess(utils::Problem const & problem, BinFunc f)
{
    auto remain = [&](){
        std::set<Index2D> dst;
        for(auto r : utils::iota(problem.div_y()))
            for(auto c : utils::iota(problem.div_x())){
                dst.insert(makeIndex2D(r, c));
            }
        return dst;
    }();


    /// 画像originのdir方向に最適な画像を選び出す
    auto choose_best_one = [&](std::set<Index2D> const & remain, Index2D origin, utils::Direction dir, double *pPV)
    {
        Index2D mIdx;
        double min = std::numeric_limits<double>::infinity();

        for(auto const & idx: remain){
            const double v = std::abs(f(problem.get_element(origin[0], origin[1]),
                                        problem.get_element(idx[0], idx[1]),
                                        dir));

            if(min >= v){   // min == v == infのときは入れ替える
                min = v;
                mIdx = idx;
            }
        }

        if(pPV)
            *pPV = min;

        return mIdx;
    };


    // 画像リストdstを元にして、縦方向か横方向に画像を結合していき、
    // 最終的に、maxN個の画像のリストになるまで結合を進めます。
    // pTopN :  先頭に何個追加したかが格納される。
    // pIncPV : 評価関数の増加値が格納される
    auto guess_bidirectional =
    [&](std::set<Index2D> & remain, std::deque<Index2D> & dst, bool isVerticalLine,
        std::size_t maxN, std::size_t *pTopN, double *pIncPV)
    {
        PROCON_ENFORCE(dst.size() >= 1, "結合素材の画像が存在しません");

        double incPV = 0;
        std::size_t topN = 0;       // 先頭に何個追加したか

        while(!remain.empty() && dst.size() < maxN){
            auto dir = isVerticalLine ? utils::Direction::up : utils::Direction::left;
            double predValue;
            Index2D mIdx = choose_best_one(remain, dst[0], dir, &predValue);

            {
                const auto dir2 = isVerticalLine ? utils::Direction::down : utils::Direction::right;
                double predV2;
                Index2D mIdx2 = choose_best_one(remain, dst[dst.size()-1], dir2, &predV2);

                if(predV2 < predValue){
                    dir = dir2;
                    predValue = predV2;
                    mIdx = std::move(mIdx2);
                }
            }

            if(dir == (isVerticalLine ? utils::Direction::up : utils::Direction::left)){    // 先頭にくっつける
                dst.push_front(mIdx);
                ++topN;
            }else                        // 後ろにくっつける
                dst.push_back(mIdx);

            incPV += predValue;
            remain.erase(mIdx);         // 決定したので消す
        }

        if(pTopN)
            *pTopN = topN;
        if(pIncPV)
            *pIncPV = incPV;
    };


    // 画像リストdstを元にして、単方向に連結していき、最終的に、maxN個の画像のリストになるまで連結を進めます。
    auto guess_singlyLink =
    [&](std::set<Index2D> & remain, std::deque<Index2D> & dst, utils::Direction dir,
        std::size_t maxN, double *pIncPV)
    {
        PROCON_ENFORCE(dst.size() > 0, "結合素材の画像がありません");

        double incPV = 0;
        while(!remain.empty() && dst.size() < maxN){
            Index2D idx;

            if(dir == utils::Direction::up || dir == utils::Direction::left)
                idx = dst[0];
            else
                idx = dst[dst.size()-1];

            double predV = 0;
            auto bestIdx = choose_best_one(remain, idx, dir, &predV);
            remain.erase(bestIdx);

            if(idx == dst[0])
                dst.push_front(bestIdx);
            else
                dst.push_back(bestIdx);
            incPV += predV;
        }

        if(pIncPV)
            *pIncPV += incPV;
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
    std::size_t m_ind = 0;
    double min = std::numeric_limits<double>::infinity();

    //垂直の結合で一番評価関数の合計値が最小になるものを探す
    for(auto i: utils::iota(problem.div_y())){
        auto org = makeIndex2D(i, 0);
        remain.erase(org);
        std::deque<Index2D> list = { org };

        double v;
        guess_bidirectional(remain, list, true, problem.div_y(), nullptr, &v);

        if(v <= min){
            min = v;
            m_ind = i;
        }

        // remainに list に入っている分を戻す
        for(auto& e: list)
            remain.insert(std::move(e));
    }

    //最小の評価関数になる断片たちを垂直結合
    std::deque<Index2D> vert = { makeIndex2D(m_ind, 0) };
    guess_bidirectional(remain, vert, true, problem.div_y(), nullptr, nullptr);

    min = std::numeric_limits<double>::infinity();
    std::size_t limLN = 0;
    //水平の結合で一番評価関数の合計値が最小になるものの左連結数と右連結数を計算
    for(auto& o : vert){
        std::deque<Index2D> list = { o };
        double v = 0;
        std::size_t leftN = 0;
        guess_bidirectional(remain, list, false, problem.div_x(), &leftN, &v);
        if(v <= min){
            min = v;
            limLN = leftN;
        }

        for(auto& e: list)
            if(e != o)
                remain.insert(std::move(e));
    }

    std::vector<std::vector<Index2D>> dst;
    //求めた左連結数と右連結数で結合
    for(auto& o : vert){
        std::deque<Index2D> hlist = { o };
        guess_singlyLink(remain, hlist, utils::Direction::left, limLN+1, nullptr);
        guess_singlyLink(remain, hlist, utils::Direction::right, problem.div_x(), nullptr);

        dst.emplace_back(hlist.begin(), hlist.end());
    }
    
    return dst;
}


/**
 画素の差がかなり小さい画素の数を数え、その個数で
 diff_connection()の評価値を割ったものを評価値とする

 誤った画像の結合は、画素の差の分散が小さい誤った画像があるときに起こりやすいと考えた

 誤るときの画素の差の列の例
 正しい断片　   [ほぼ同じ][ほぼ同じ][ほぼ同じ][かなり異なる][かなり異なる]
 誤っている断片 [やや異なる][やや異なる][やや異なる][やや異なる][やや異なる]

 こういうときに画素の差を平均して正しい断片より誤っている断片が
 優れていることになってしまうことがある
*/

template <typename T, typename U
#ifdef SUPPORT_CONSTRAINTS
    , PROCON_TEMPLATE_CONSTRAINTS(utils::is_image<T>() && utils::is_image<U>())   // T, Uともに画像であるという制約
#endif
>
double diff_connection_rena(T const & img1, U const & img2, utils::Direction direction)
{
    double eval = 0;
    double num = 0;

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
            for(std::size_t i = 0; i < 3; ++i){
                auto v = std::abs(static_cast<float>(p1[i]) - static_cast<float>(p2[i]));
                eval += v;

                if(v < 7) ++num;
            }
        }
        eval /= img1.height();
        break;

      case utils::Direction::up:
      case utils::Direction::down:
        for(std::size_t c = 0; c < img1.width(); ++c){
            auto p1 = img1.get_pixel(r1, c).vec();
            auto p2 = img2.get_pixel(r2, c).vec();
            for(std::size_t i = 0; i < 3; ++i){
                auto v = std::abs(static_cast<float>(p1[i]) - static_cast<float>(p2[i]));
                eval += v;

                if(v < 7) ++num;
            }
        }
        eval /= img1.width();
        break;
    }

    return eval / (num*100 + 1);
}
}} 
