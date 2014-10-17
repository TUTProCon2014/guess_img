#pragma once

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "../../utils/include/image.hpp"
#include "../../utils/include/types.hpp"
#include "../../utils/include/constants.hpp"


namespace procon { namespace guess {

class Ajust
{
public:
	static void ajust( std::vector<float>& pxs , int a )
	{
		std::vector<float> tpxs = pxs;
		int w = pxs.size() / a;

		//{‚PB|‚P ‚Ì‚¢‚¿‚ğ‘«‚µ‡‚í‚¹‚é
		for (size_t k = 0; k < w; ++k)
		{
			int km1 = k <= 0 ? k : k - 1;
			int kp1 = k >= w - 1 ? k : k + 1;

			for (int i = 0; i < a; i++)
			{
				pxs.at(k * 3 + i) += tpxs.at(km1 * 3 + i)*0.5 + tpxs.at(kp1 * 3 + i)*0.5;
			}
		}
	}

	static void ajust_s(std::vector<float>& pxs, int a)
	{
		std::vector<float> tpxs = pxs;
		int w = pxs.size() / a;

		// 1 , 0 , -1 ‚Ì‚¢‚¿‚ğ‘«‚µ‡‚í‚¹‚é
		for (size_t k = 0; k < w; ++k)
		{
			int km1 = k <= 0 ? k : k - 1;
			int kp1 = k >= w - 1 ? k : k + 1;

			for (int i = 0; i < a; i++)
			{
				pxs.at(k * 3 + i) = tpxs.at(km1 * 3 + i) - tpxs.at(kp1 * 3 + i);
			}
		}
	}
};

struct Correlator
{
	//‚·‚×‚Ä‚Ì‰æ‘œ•Ğ‚É‘Î‚µ‚ÄAã‰º¶‰E‚ÅŠÖ˜A•t‚¯‚Ä‚»‚Ì•ûŒü‚Pƒrƒbƒg‚Ì‰æ‘f”‚ğ“o˜^
    Correlator(utils::Problem const & pb)
    {
        const size_t w = pb.width() / pb.div_x();
        const size_t h = pb.height() / pb.div_y();

        for(size_t i = 0; i < pb.div_y(); ++i)
            for(size_t j = 0; j < pb.div_x(); ++j)
			{
				//‰æ‘œ‚Ìæ“¾
                auto& img = pb.get_element(i, j);

				//•ûŒü”z—ñ‚Ìì¬[‰EAãA¶A‰º]
                utils::Direction dirs[4] = {utils::Direction::right, utils::Direction::up, utils::Direction::left, utils::Direction::down};

				//•ûŒü‚Æ‰æ‘f”
                std::unordered_map<utils::Direction, std::vector<float>> pxsMap;
				std::unordered_map<utils::Direction, std::vector<float>> pxsMap_s;

				//•ûŒü‚Ì”‚¾‚¯ƒ‹[ƒv
                for(auto dir: dirs)
				{
					//ã‚©‰º‚È‚ç
                    if(dir == utils::Direction::up || dir == utils::Direction::down)
					{
                        const auto l = dir == utils::Direction::up ? 0 : h-1;
						std::vector<float> pxs, pxs_s; pxs_s.reserve(w); pxs.reserve(w);
						for(size_t k = 0; k < w; ++k)
						{
                            auto v = img.get_pixel(l, k).vec();
							for (size_t pidx = 0; pidx < 3; ++pidx)
							{
								pxs.emplace_back(v[pidx]);
								pxs_s.emplace_back(v[pidx]);
							}
						}
						Ajust::ajust_s(pxs_s, 3);
                        pxsMap.emplace(dir, std::move(pxs));
						pxsMap_s.emplace(dir, std::move(pxs_s));
					}
                    else if(dir == utils::Direction::right || dir == utils::Direction::left){
                        const auto l = dir == utils::Direction::left ? 0 : w-1;
						std::vector<float> pxs, pxs_s; pxs_s.reserve(h); pxs.reserve(h);
						for (size_t k = 0; k < h; ++k)
						{
							auto v = img.get_pixel(k, l).vec();
							for (size_t pidx = 0; pidx < 3; ++pidx)
							{
								pxs.emplace_back(v[pidx]);
								pxs_s.emplace_back(v[pidx]);
							}
						}
						Ajust::ajust_s(pxs_s, 3);
						pxsMap.emplace(dir, std::move(pxs));
						pxsMap_s.emplace(dir, std::move(pxs_s));
					}
                }

                this->_memo.emplace(utils::ImageID(i, j), std::move(pxsMap));
				this->_memo_s.emplace(utils::ImageID(i, j), std::move(pxsMap_s));
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


		double sum = 0;
		double var = 0;
		double ave = 0;

		{
			//•’Ê‚Ì
			auto& px1 = _memo.at(img1).at(dir);
			auto p1 = &(px1[0]);
			auto p2 = &(_memo.at(img2).at(dir2)[0]);

			const size_t n = px1.size();
			const auto e1 = p1 + n;

			//•½‹Ï
			while (p1 != e1){
				sum += std::abs(*p1 - *p2);
				++p1; ++p2;
			}

			p1 = &(px1[0]);
			p2 = &(_memo.at(img2).at(dir2)[0]);

			ave = sum / n;

			while (p1 != e1){
				auto x = std::abs(*p1 - *p2);

				var += std::abs((x - ave));

				++p1; ++p2;
			}

			var /= n;
		}

		double sums = 0;
		double vars = 0;
		double aves = 0;

		{
			//”÷•ª‚Ì
			auto& px1 = _memo_s.at(img1).at(dir);
			auto p1 = &(px1[0]);
			auto p2 = &(_memo_s.at(img2).at(dir2)[0]);

			const size_t n = px1.size();
			const auto e1 = p1 + n;

			//•½‹Ï
			while (p1 != e1){
				sums += std::abs(*p1 - *p2);
				++p1; ++p2;
			}

			p1 = &(px1[0]);
			p2 = &(_memo_s.at(img2).at(dir2)[0]);

			aves = sums / n;

			while (p1 != e1){
				auto x = std::abs(*p1 - *p2);
				vars += (x - aves)*(x - aves);
				++p1; ++p2;
			}

			vars /= n;
		}


        return ave * var * aves * vars;
    }


  private:
    std::unordered_map<utils::ImageID, std::unordered_map<utils::Direction, std::vector<float>>> _memo;
	std::unordered_map<utils::ImageID, std::unordered_map<utils::Direction, std::vector<float>>> _memo_s;
};

}}