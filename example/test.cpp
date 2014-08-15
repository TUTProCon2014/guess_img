#include "../include/guess.hpp"
#include "../../inout/include/inout.hpp"
#include "../../utils/include/types.hpp"

using namespace procon;

int main()
{
    auto p_opt =  inout::get_problem_from_test_server(8);
    if(p_opt){
        const utils::Problem& p = *p_opt;

        // 復元に使うための、2つの画像のくっつき度合いを返す関数
        guess::DiffConnection pred;

        // 復元
        auto idxs = guess::guess(p, pred);

        // 復元できたインデックスを表示してみる
        for(auto& ee: idxs){
            for(auto& e : ee)
                std::cout << "(" << e[0] << ", " << e[1] << "), ";
            std::cout << std::endl;
        }

        auto& src_img = p.cvMat();
        cv::Mat dst = src_img.clone();

        // 画像を復元してみる
        {
            size_t r = 0;
            for(auto ee: idxs){
                size_t c = 0;
                for(auto e: ee)
                {
                    auto elem = p.get_element(e[0], e[1]);
                    auto hh = p.height() / p.div_y();
                    auto ww = p.width() / p.div_x();

                    for(auto rr: utils::iota(hh))
                        for(auto cc: utils::iota(ww))
                            dst.at<cv::Vec3b>(r * hh + rr, c * ww + cc)
                                = elem.get_pixel(rr, cc).vec();
                    ++c;
                }
                ++r;
            }
        }

        cv::namedWindow("image1", cv::WINDOW_AUTOSIZE|cv::WINDOW_FREERATIO);
        
        // // ウィンドウ名でウィンドウを指定して，そこに画像を描画
        cv::imshow("image1", dst);
        // // キー入力を（無限に）待つ
        cv::waitKey(0);
    }else
        std::cout << "死" << std::endl;
}
