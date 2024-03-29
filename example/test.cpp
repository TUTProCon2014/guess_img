#include "../include/guess.hpp"
#include "../include/pso_guess.hpp"
#include "../include/rena_guess.hpp"
#include "../include/bfs_guess.hpp"
#include "../include/blocked_guess.hpp"
#include "../include/correlation.hpp"
#include "../../inout/include/inout.hpp"
#include "../../utils/include/types.hpp"
#include "../../utils/include/dwrite.hpp"

#define GUESS_FUNC blocked_guess::guess

using namespace procon;

int main()
{
    auto p_opt = utils::Problem::get("img9.ppm");

    if(!p_opt)
        p_opt = inout::get_problem_from_test_server(9);

    if(p_opt){
        const utils::Problem& p = *p_opt;

        // 復元に使うための、2つの画像のくっつき度合いを返す関数
        auto pred = guess::Correlator(p);

        // 復元
        auto idxs = GUESS_FUNC(p, pred);

        // 復元できたインデックスを表示してみる
        for(auto& ee: idxs){
            for (auto& e : ee){
                auto i2 = e.get_index();
                std::cout << "(" << i2[0] << ", " <<i2[1] << "), ";
            }
            std::cout << std::endl;
        }

        auto dst = p.clone();

        // 画像を復元してみる
        {
            size_t r = 0;
            for(auto ee: idxs){
                size_t c = 0;
                for(auto e: ee){
                    p.get_element(e).cvMat().copyTo(dst.get_element(r, c).cvMat());
                    ++c;
                }
                ++r;
            }
        }

        cv::namedWindow("image1", cv::WINDOW_AUTOSIZE);
        
        // // ウィンドウ名でウィンドウを指定して，そこに画像を描画
        cv::imshow("image1", dst.cvMat());
        // // キー入力を（無限に）待つ
        cv::waitKey(0);
    }else
        std::cout << "死" << std::endl;

    return 0;
}
