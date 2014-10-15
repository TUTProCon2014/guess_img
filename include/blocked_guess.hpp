#pragma once

#include <deque>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <thread>
#include <future>

#include "../../utils/include/image.hpp"
#include "../../modify_guess_image/common.hpp"
#include "../../modify_guess_image/interactive_guess.hpp"
#include "../../utils/include/dwrite.hpp"


namespace procon { namespace blocked_guess {

using namespace utils;


int getLogExp2(int n)
{
    double x = n;
    x = std::frexp(x, &n);
    return n-1;
}


template <typename BinFunc>
modify::Group createGroup(Problem const & pb, BinFunc const & f, ImageID origin, unsigned int n,
                                        std::unordered_set<ImageID>& remain)
{
    if(n > pb.div_y())
        n = pb.div_y();

    std::deque<ImageID> list;
    list.emplace_back(origin);
    remain.erase(origin);
    while(list.size() < n){
        Direction dir = Direction::up;
        ImageID mIdx;
        double min = std::numeric_limits<double>::infinity();

        for(auto dI : iota(2)){
            Direction d = dI == 0 ? Direction::up : Direction::down;

            std::size_t tgtIdx = 0;
            if(dI == 1)
                tgtIdx = list.size() - 1;

            for(auto& id: remain){
                const double v = std::abs(f(list[tgtIdx],
                                            id,
                                            d));
                if(min >= v){
                    min = v;
                    dir = d;
                    mIdx = id;
                }
            }
        }

        if(dir == Direction::up)
            list.push_front(mIdx);
        else
            list.push_back(mIdx);

        remain.erase(mIdx);
    }

    writeln(list);

    modify::Group group;
    for(auto& e: list){
        group.emplace_back(e, std::array<std::ptrdiff_t, 2>(
                                {static_cast<std::ptrdiff_t>(group.size()), 0}));
    }


    std::array<std::ptrdiff_t, 2> fst = std::get<1>(group[0]);
    for(auto& e: group){
        auto& lst = std::get<1>(e);
        lst[0] -= fst[0];
        lst[1] -= fst[1];
    }

    return group;
}


template <typename BinFunc>
std::vector<std::vector<ImageID>> guess(Problem const & pb, BinFunc const & f)
{
    std::unordered_set<ImageID> remain;
    DividedImage::foreach(pb, [&](size_t i, size_t j){
        remain.emplace(i, j);
    });

    // all none
    modify::OptionalMap omp(pb.div_y(), std::vector<boost::optional<ImageID>>(pb.div_x()));


    std::vector<std::future<std::tuple<double, modify::ImgMap>>> ths;
    for(auto i: utils::iota(pb.div_x() / getLogExp2(pb.div_x()))){
        ths.emplace_back(std::async(
            std::launch::async,
            [&](size_t i){
                std::unordered_set<ImageID> rm = remain;
                auto gp = createGroup(pb, f, ImageID(0, i), getLogExp2(pb.div_x()) * 2, rm);
                return modify::position_bfs(&gp, &gp + 1, omp, rm, pb, f);
            },
            i));
    }


    std::tuple<double, modify::ImgMap> min;
    std::get<0>(min) = std::numeric_limits<double>::infinity();

    for(auto& e: ths){
        e.wait();
        auto res = e.get();

        if(std::get<0>(min) >= std::get<0>(res))
            min = res;
    }

    return std::get<1>(min);
}

}}