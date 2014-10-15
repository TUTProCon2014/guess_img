/*
  使い方
  test.cppにこのファイルをインクルードして、pso_guess::pso_guess()を呼び出す

  Particle Swarm Optimization(粒子群最適化)という
  メタヒューリスティクスにより画像復元を行う

  -PSOは粒子たちを探索空間内にばらまき、各粒子の位置は断片の
  ある配置を表現した数値を意味する
  -PSOははじめは探索空間内を全体的に探索し、最後の方は
  良い解がすでに見つかっている場所(良い解がありそうな場所)
  付近を重点探索する。
  -繰り返し回数、粒子数が多ければ多いほど性能は上がるが、
  速度が落ちる。(このプログラムじゃあんまり性能も上がらない...)
*/

#pragma once

#include <algorithm>
#include <vector>
#include <iostream>
#include <memory>
#include <cmath>
#include <limits>
#include <random>

namespace procon{ namespace pso_guess {

using namespace utils;

template <typename BinFunc>
class Particle{
    private:
        std::vector<double> _x;                         //位置ベクトル
        std::vector<size_t> _dis_x;                     //重複除去と離散化した位置ベクトル
        std::vector<double> _v;                         //速度ベクトル
        std::vector<double> _pbest;                     //personal best
        double _pvalue;                                 //個人の最高評価値
        int _dim;                                       //問題の次元
        BinFunc const *_f;                              //評価関数
        const Problem _problem;                         //問題情報
        std::mt19937 _rnd;                              //擬似乱数生成器
        std::uniform_real_distribution<double> _dist;   //一様分布生成器
        const double _c1 = 2.0;                         //移動係数(pbestへの近づきやすさ)
        const double _c2 = 2.0;                         //移動係数(gbestへの近づきやすさ)

    public:
        //粒子のランダム生成を行うコンストラクタ
        Particle(BinFunc f, Problem pro)
            : _f(f), _problem(pro), _rnd(std::random_device()()), _dist(0.0, 1.0)
        {
            _dim = _problem.div_x() * _problem.div_y(); //次元の計算

            for(int i=0; i < _dim; i++){
                //0から断片数までの範囲の値を生成
                double temp = _dist(_rnd) * _dim;
                _x.push_back(temp);
                _dis_x.push_back(0);
                //-断片数から断片数までの範囲の値を速度とする
                _v.push_back(_dist(_rnd) * 1.0 * _dim - _dim/2.0);
            }

            this->format(); //PSOの解を今回の問題の解に変換
            _pbest = _x;    //はじめは初期位置がこれまでの最善解となる
            _pvalue = std::numeric_limits<double>::max(); //pvalueを最大値で初期化
            this->calc_pvalue(); //pvalueの計算
        }

        //位置ベクトルを返す
        std::vector<double> x() const{
            return _x;
        }

        //速度ベクトルを返す
        std::vector<double> v() const{
            return _v;
        }

        //personal bestを返す
        std::vector<double> pbest() const{
            return _pbest;
        }

        //粒子個人の最高評価値を返す
        double pvalue() const{
            return _pvalue;
        }

        //PSOの解を今回の問題に変換(離散化+重複除去)
        //
        //PSOは連続探索空間を探索するアルゴリズム
        //今回以下のように断片の位置を番号付けする
        // | 1 2 3 |
        // | 4 5 6 |
        // | 7 8 9 |
        //PSOはそのまま実行すると粒子が同じ値を持つ可能性があるので
        //そのまま粒子の持つxを上記の番号とすることはできない
        //(重複したら同じ断片を２回使うことになる)
        //そこで0から次元までのベクトルのi番目を取り出すというときのiを粒子のもつxとする
        // [0 1 2 3 4 5 6 7 8]
        //取り出すごとにベクトルからその要素を消し去り、次の粒子のxは0から次元-1の値にする
        //こうすることで重複をなくせる
        void format(){
            for(int i=0; i < _dim; i++){
                double temp;
                temp = _x[i] * (_dim-i)*1.0/_dim; //重複除去
                if(temp >= _dim - i) temp = _dim - i - 1; //_x[i] = _dim - iのときのエラーを防ぐ
                _dis_x[i] = (size_t)floor(temp); //離散化
            }
        }

        //Indexの2次元配列に変換する
        std::vector<std::vector<ImageID>> make_indexv(){
            std::vector<std::vector<ImageID>> ret(_problem.div_y(), std::vector<ImageID>(_problem.div_x()));

            //0から次元までの数のベクトルを生成
            std::vector<int> remain;
            for(int i=0; i < _dim; i++){
                remain.push_back(i);
            }

            //index配列の生成
            for(int i=0; i < _dim; i++){
                int select = remain[_dis_x[i]]; //要素の選択

                remain.erase(std::remove(remain.begin(), remain.end(), remain[_dis_x[i]]), remain.end()); //選ばれた要素の削除

                Index2D idx;
                idx[1] = (size_t)select % _problem.div_x();     //x
                idx[0] = (size_t)(select / _problem.div_x());   //y
                size_t x = i % _problem.div_x();
                size_t y = (size_t)(i / _problem.div_x());
                ret[y][x] = ImageID(idx);
            }

            return ret;
        }

        //pvalueの計算とpbestの更新
        void calc_pvalue(){
            double val = 0;
            int dx[4] = {1, 0, -1, 0};
            int dy[4] = {0, -1, 0, 1};
            auto idxs = this->make_indexv(); //index2次元配列を受け取る

            //各断片について周りの断片との結合度を評価し、その合計を評価値とする
            for(int i=0; i < _problem.div_x(); i++){
                for(int j=0; j < _problem.div_y(); j++){
                    for(int k=0; k<4; k++){
                        size_t sx = i + dx[k];
                        size_t sy = j + dy[k];
                        if(sx >= 0 && sx < _problem.div_x() && sy >= 0 && sy < _problem.div_y()){
                            
                            double v = std::abs((*_f)(idxs[j][i],
                                                   idxs[sy][sx],
                                                  (utils::Direction)k));
                            val += v;
                        }
                    }
                }
            }

            //_pvalueの更新
            if(val < _pvalue){
                _pvalue = val;
                _pbest = _x;
            }
        }

        //Particleの移動
        void move(double w, std::vector<double> gbest){
            //乱数項の生成
            std::vector<double> escape; 
            for(int j=0; j<_dim; j++){
                escape.push_back(-_dim + _dist(_rnd) * 2.0 * _dim);
            }

            for(int i=0; i < _dim; i++){
                //粒子の速度の更新
                _v[i] = w * _v[i] + _c1 * _dist(_rnd) * (_pbest[i] - _x[i]) + _c2 * _dist(_rnd) * (gbest[i] - _x[i]) + 0.01 * escape[i];
                //粒子の速度が限界を超えるのを防ぐ
                if(_v[i] >= _dim/2.0) _v[i] = _dim/2.0;
                else if(_v[i] <= -_dim/2.0) _v[i] = -_dim/2.0;

                //粒子の位置の更新
                if(abs(_v[i]) > 10e-9){ //あまりに小さい数を足すと速度が遅くなることが以前あった
                    _x[i] = _x[i] + _v[i];

                    //粒子が探索領域外へ出ることを防ぐ(粒子は探索領域の壁で反射する)
                    if(_x[i] <= 0){
                        _x[i] = abs(_x[i]);
                        _v[i] *= -1;
                    }else if(_x[i] >= _dim){
                        _x[i] = _dim - abs(_x[i] - _dim);
                        _v[i] *= -1;
                    }
                }
            }

            this->format(); //解を問題に合わせる
            this->calc_pvalue(); //pbestの計算
        }
};

template <typename BinFunc>
std::vector<std::vector<ImageID>> pso_guess(utils::Problem const & problem, BinFunc const & f){
    const int p_num = 30;       //粒子数
    const int tmax = 300;       //イテレーション回数
    double w = 0.9;             //慣性項
    double gvalue;              //global best(最終的な解)
    std::vector<double> gbest;  //粒子みんなの最高評価値
    std::vector<std::vector<ImageID>> dst; //答えとなるインデックス2次元配列

    //粒子の生成
    std::vector<Particle<BinFunc>> p;
    for(int i=0; i < p_num; i++){
        Particle<BinFunc> t(f, problem);
        p.push_back(t);
    }

    //gbestの初期化
    gvalue = p[0].pvalue();
    gbest = p[0].pbest();
    dst = p[0].make_indexv();
    for(int i=1; i < p_num; i++){
        if(gvalue > p[i].pvalue()){
            gvalue = p[i].pvalue();
            gbest = p[i].pbest();
            dst = p[i].make_indexv();
        }
    }


    //粒子による探索
    for(int i=0; i < tmax; i++){
        //以下の1行で評価の遷移をみれる
        std::cout << i << " " << gvalue << std::endl;

        //本来のPSOでは以下の１行を入れるほうがいいはずだが、今回の探索空間ではないほうがよさそうかも
        w = 0.9 - i*1.0/tmax * 0.5;

        //粒子の移動 + pbestの更新
        for(int j=0; j < p_num; j++){
            p[j].move(w, gbest);
        }

        //gbestの更新
        for(int j=0; j < p_num; j++){
            if(gvalue > p[j].pvalue()){
                gvalue = p[j].pvalue();
                gbest = p[j].pbest();
                dst = p[j].make_indexv();
            }
        }
    }
    
    //最終評価
    std::cout << "final result : " << gvalue << std::endl;  

    return dst;
}

}}
