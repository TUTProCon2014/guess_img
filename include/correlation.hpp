#pragma once

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "../../utils/include/image.hpp"
#include "../../utils/include/types.hpp"
#include "../../utils/include/constants.hpp"


namespace procon { namespace guess {

struct Correlator
{
    Correlator(utils::Problem const & pb)
    {
        const size_t w = pb.width() / pb.div_x();
        const size_t h = pb.height() / pb.div_y();

        for(size_t i = 0; i < pb.div_y(); ++i)
            for(size_t j = 0; j < pb.div_x(); ++j){
                auto& img = pb.get_element(i, j);
                utils::Direction dirs[4] = {utils::Direction::right, utils::Direction::up, utils::Direction::left, utils::Direction::down};

                std::unordered_map<utils::Direction, std::vector<float>> pxsMap;

                for(auto dir: dirs){
                    if(dir == utils::Direction::up || dir == utils::Direction::down){
                        const auto l = dir == utils::Direction::up ? 0 : h-1;
                        std::vector<float> pxs; pxs.reserve(w);
                        for(size_t k = 0; k < w; ++k){
                            auto v = img.get_pixel(l, k).vec();

                            for(size_t pidx = 0; pidx < 3; ++pidx)
                                pxs.emplace_back(v[pidx]);
                        }
                        pxsMap.emplace(dir, std::move(pxs));
                    }
                    else if(dir == utils::Direction::right || dir == utils::Direction::left){
                        const auto l = dir == utils::Direction::left ? 0 : w-1;
                        std::vector<float> pxs; pxs.reserve(h);
                        for(size_t k = 0; k < h; ++k){
                            auto v = img.get_pixel(k, l).vec();

                            for(size_t pidx = 0; pidx < 3; ++pidx)
                                pxs.emplace_back(v[pidx]);
                        }
                        pxsMap.emplace(dir, std::move(pxs));
                    }
                }

                this->_memo.emplace(utils::ImageID(i, j), std::move(pxsMap));
            }
    }


    double operator()(utils::ImageID const & img1, utils::ImageID const & img2, utils::Direction dir) const
    {
        const auto dir2 = [&](){
            switch(dir){
                case utils::Direction::right: return utils::Direction::left;
                case utils::Direction::up:    return utils::Direction::down;
                case utils::Direction::left:  return utils::Direction::right;
                case utils::Direction::down:  return utils::Direction::up;
                default:
                    PROCON_ENFORCE(0, "Switch error");
                    return utils::Direction::right;
            }
        }();

        auto& px1 = _memo.at(img1).at(dir);
        auto p1 = &(px1[0]);
        auto p2 = &(_memo.at(img2).at(dir2)[0]);

        double sum = 0;
        const size_t n = px1.size();
        const auto e1 = p1 + n;
        while(p1 != e1){
            sum += std::abs(*p1 - *p2);
            ++p1; ++p2;
        }

        return sum / n;
    }


  private:
    std::unordered_map<utils::ImageID, std::unordered_map<utils::Direction, std::vector<float>>>
        _memo;
};

}}