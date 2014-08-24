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

	//blim個後ろに断片を結合した後、flimだけ先頭に断片を結合
    auto guess_oneline_lim = [&](Index2D origin, bool isVerticalLine, int flim, int blim){
        std::deque<Index2D> dst;
        dst.push_back(origin);      // 最初に原点がある
        for(auto t : utils::iota(isVerticalLine ? problem.div_y()-1 : blim)){
            utils::Direction dir;
            Index2D mIdx;
            double min = std::numeric_limits<double>::infinity();

            // std::array<std::size_t, 2> ds = {0, dst.size()-1};    // {0 : 上(左), dst.size()-1 : 下(右)}
			utils::Direction d = utils::Direction::right;

			std::size_t tgtIdx = dst.size() - 1;

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
            if(dir == (isVerticalLine ? utils::Direction::up : utils::Direction::left))    // 先頭にくっつける
                dst.push_front(mIdx);
            else                        // 後ろにくっつける
                dst.push_back(mIdx);

            remain.erase(mIdx);         // 決定したので消す
        }


        for(auto t : utils::iota(isVerticalLine ? problem.div_y()-1 : flim)){
            utils::Direction dir;
            Index2D mIdx;
            double min = std::numeric_limits<double>::infinity();

            // std::array<std::size_t, 2> ds = {0, dst.size()-1};    // {0 : 上(左), dst.size()-1 : 下(右)}
			utils::Direction d = utils::Direction::left;

			std::size_t tgtIdx = 0;

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
            if(dir == (isVerticalLine ? utils::Direction::up : utils::Direction::left))    // 先頭にくっつける
                dst.push_front(mIdx);
            else                        // 後ろにくっつける
                dst.push_back(mIdx);

            remain.erase(mIdx);         // 決定したので消す
        }

        std::vector<Index2D> vec(dst.begin(), dst.end());
        return vec;
    };

	//remainの要素を消さずに評価関数の合計値を求める(最小の評価関数の合計値となる結合の前と後ろに結合した数がleftとrightに入る)
    auto guess_minline = [&](Index2D origin, bool isVerticalLine, int& left, int & right){
		auto rt = remain;
		double ret = 0; //直線に結合したときの、評価関数の合計値
		left = 0;
		right = 0;
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

                for(auto& idx : rt){    // 残っている画像の中から探す
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
			ret += min;

            if(dir == (isVerticalLine ? utils::Direction::up : utils::Direction::left)){    // 先頭にくっつける
                dst.push_front(mIdx);
				left++;
            }else{                        // 後ろにくっつける
                dst.push_back(mIdx);
				right++;
			}

            rt.erase(mIdx);         // 決定したので消す
        }

        std::vector<Index2D> vec(dst.begin(), dst.end());
        return ret;
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
	int m_ind = 0;
	int liml = 0;
	int limr = 0;
	int min = std::numeric_limits<int>::max();

	//垂直の結合で一番評価関数の合計値が最小になるものを探す
	for(int i=0; i<problem.div_y(); i++){
		int left, right;
		auto v = guess_minline(makeIndex2D(i, 0), true, left, right);
		if(v < min){
			min = v;
			m_ind = i;
		}
	}

	//最小の評価関数になる断片たちを垂直結合
	auto vert = guess_oneline(makeIndex2D(m_ind, 0), true);

	min = std::numeric_limits<int>::max();

	//水平の結合で一番評価関数の合計値が最小になるものの左連結数と右連結数を計算
    for(auto& o : vert){
		int left, right;
		auto v = guess_minline(o, false, left, right);
		if(v < min){
			min = v;
			liml = left;
			limr = right;
		}
	}

	//求めた左連結数と右連結数で結合
	for(auto& o : vert){
		dst.push_back(guess_oneline_lim(o, false, liml, limr));
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

#ifdef NOT_SUPPORT_CONSTEXPR
template <typename T, typename U>
#else
template <typename T, typename U,
    PROCON_TEMPLATE_CONSTRAINTS(utils::is_image<T>() && utils::is_image<U>())>   // T, Uともに画像であるという制約
#endif
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
            for(std::size_t i = 0; i < 3; ++i)
                eval += std::abs(static_cast<float>(p1[i]) - static_cast<float>(p2[i]));
        }
        eval /= img1.height();
        break;

      case utils::Direction::up:
      case utils::Direction::down:
        for(std::size_t c = 0; c < img1.width(); ++c){
            auto p1 = img1.get_pixel(r1, c).vec();
            auto p2 = img2.get_pixel(r2, c).vec();
            for(std::size_t i = 0; i < 3; ++i)
                eval += std::abs(static_cast<float>(p1[i]) - static_cast<float>(p2[i]));
        }
        eval /= img1.width();
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
				if(std::abs(static_cast<float>(p1[i]) - static_cast<float>(p2[i])) < 7){
					num++;
				}
			}
        }
        break;

      case utils::Direction::up:
      case utils::Direction::down:
        for(std::size_t c = 0; c < img1.width(); ++c){
            auto p1 = img1.get_pixel(r1, c).vec();
            auto p2 = img2.get_pixel(r2, c).vec();
            for(std::size_t i = 0; i < 3; ++i){
				if(std::abs(static_cast<float>(p1[i]) - static_cast<float>(p2[i])) < 7){
					num++;
				}
			}
        }
        break;
    }
	eval = eval / (num*100 + 1);

    return eval;
}
}} 
