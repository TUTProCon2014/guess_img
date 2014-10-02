#pragma once

#include "../../inout/include/inout.hpp"
#include "../../utils/include/template.hpp"
#include "../../utils/include/types.hpp"
#include "../../utils/include/range.hpp"
#include "../../utils/include/dwrite.hpp"
#include "guess.hpp"

#include <vector>
#include <set>
#include <array>
#include <deque>
#include <deque>
#include <tuple>
#include <thread>

namespace procon { namespace bfs_guess {

using namespace utils;


inline ImageID convToImageID(std::size_t i, std::size_t w)
{
    return ImageID(i / w, i % w);
}


/**
T型のもつvalue()を評価し、その平均値と標準偏差を計算し、predの評価結果がtrueな要素をdstに入れます
その際、srcは破壊します。
*/
template <typename T, typename Pred>
void refine(std::deque<T>&& src, std::deque<T>& dst, Pred pred)
{
    double m = 0, v = 0;
    for(auto& e: src){
        const double vv = e.value();
        m += vv;
        v += vv * vv;
    }

    m /= src.size();
    v = v / src.size() - m * m;
    const double sd = std::sqrt(v);

    for(auto& e: src)
        if(pred(e.value(), m, sd))
            dst.push_back(std::move(e));
    src.clear();
}



template <typename State>
void bfs_guess_impl(std::deque<State>& state, std::size_t maxSize)
{
    std::size_t iterCount = 0;
    while(!state.empty() && !state[0].isEnd()){
        std::deque<State> dst1;
        for(auto& e: state){
            e.update(dst1);

            std::deque<State> dst2;
            refine(std::move(dst1), dst2, [](double v, double m, double sd){ return v <= m + sd; });
            dst1 = [](){ std::deque<State> ini; return ini; }();

            std::set<State> rbtree;
            for(auto& e: dst2)
                rbtree.insert(std::move(e));

            for(auto& e: rbtree){
                dst1.push_back(std::move(e));
                if(dst1.size() >= maxSize)
                    break;
            }
        }

        state = std::move(dst1);
    }
}



template <typename BinFunc>
struct State1st
{
    State1st(utils::Problem const *pb, BinFunc *func, std::deque<ImageID> const & idx, std::vector<bool> const & rem, double ev)
    : _pb(pb), _pred(func), _idx(idx), _remain(rem), _ev(ev) {}


    double value() const
    {
        return _ev;
    }


    std::deque<ImageID> const & index() const
    {
        return _idx;
    }


    bool operator==(State1st<BinFunc> const & rhs) const
    { return utils::equal(index(), rhs.index()); }
    bool operator!=(State1st<BinFunc> const & rhs) const
    { return !(*(this) == rhs); }
    bool operator<(State1st<BinFunc> const & rhs) const
    { return this->value() < rhs.value(); }


    bool isEnd() const { return _idx.size() == _pb->div_y(); }


    void update(std::deque<State1st> & q)
    {
        std::size_t i = 0;
        for(bool e: _remain){
            if(e){
                State1st<BinFunc> dupTop = *this;
                dupTop.insert(Direction::up, i);
                q.push_back(std::move(dupTop));

                State1st<BinFunc> dupBottom = *this;
                dupBottom.insert(Direction::down, i);
                q.push_back(std::move(dupBottom));
            }
            ++i;
        }
    }


    Problem const *_pb;
    BinFunc *_pred;
    std::deque<ImageID> _idx;
    std::vector<bool> _remain;
    double _ev;


    double pred_value(Direction dir, ImageID const & index)
    {
        ImageID tI;
        if(dir == Direction::up)
            tI = _idx[0];
        else
            tI = _idx[_idx.size()-1];

        return std::abs((*_pred)(_pb->get_element(tI),
                        _pb->get_element(index),
                        dir));
    }


    void insert(Direction dir, std::size_t i)
    {
        const ImageID index = convToImageID(i, _pb->div_x());

        _ev += pred_value(dir, index);

        if(dir == Direction::up)
            _idx.push_front(index);
        else
            _idx.push_back(index);

        _remain[i] = false;
    }
};


template <typename BinFunc>
struct State2nd
{
    explicit State2nd(State1st<BinFunc>&& state)
    : _1st(state), _idx({state.index()[0]}), _cntLN(0) {}


    double value() const { return _1st.value(); }


    bool operator==(State2nd<BinFunc> const & rhs) const
    {
        return _1st == rhs._1st && utils::equal(_idx, rhs._idx);
    }


    bool operator!=(State2nd<BinFunc> const & rhs) const
    {
        return !(*this == rhs);
    }

    bool operator<(State2nd<BinFunc> const & rhs) const
    { return this->value() < rhs.value(); }


    bool isEnd() const { return _idx.size() == _1st._pb->div_x(); }


    void update(std::deque<State2nd>& dst)
    {
        std::size_t i = 0;
        for(bool e: _1st._remain){
            if(e){
                auto dupLeft = *this;
                dupLeft.insert(Direction::left, i);
                dst.push_back(std::move(dupLeft));

                auto dupRight = *this;
                dupRight.insert(Direction::right, i);
                dst.push_back(std::move(dupRight));
            }
            ++i;
        }
    }


    std::deque<ImageID> const & index2nd()
    {
        return _idx;
    }


    std::deque<ImageID> const & index1st()
    {
        return _1st.index();
    }


    std::size_t insertLeftCount() const { return _cntLN; }


    State1st<BinFunc> _1st;
    std::deque<ImageID> _idx;
    std::size_t _cntLN;


    double pred_value(Direction dir, ImageID const & index)
    {
        auto tI = _idx[0];
        if(dir == Direction::right)
            tI = _idx[_idx.size() - 1];

        return std::abs((*_1st._pred)(_1st._pb->get_element(tI),
                                      _1st._pb->get_element(index),
                                      dir));
    }


    void insert(Direction dir, std::size_t i)
    {
        const ImageID index = convToImageID(i, _1st._pb->div_x());

        _1st._ev += pred_value(dir, index);

        if(dir == Direction::left){
            _idx.push_front(index);
            ++_cntLN;
        }
        else
            _idx.push_back(index);

        _1st._remain[i] = false;
    }
};


template <typename BinFunc>
struct State3rd
{
    explicit State3rd(State2nd<BinFunc>&& state)
    : _1st(state._1st), _idx({state._idx}), _cntLN(state._cntLN), _ctIdx(0) {
        for(auto i: utils::iota(static_cast<size_t>(1), _1st._pb->div_y()))
            _idx.push_back({_1st.index()[i]});
    }


    double value() const { return _1st.value(); }


    bool operator==(State3rd<BinFunc> const & rhs)
    {
        if(_cntLN != rhs._cntLN)
            return false;

        if(_idx.size() != rhs._idx.size())
            return false;

        for(auto i: utils::iota(_idx.size()))
            if(!utils::equal(_idx[i], rhs._idx[i]))
                return false;

        return true;
    }


    bool operator!=(State3rd<BinFunc> const & rhs)
    {
        return !(*this == rhs);
    }

    bool operator<(State3rd<BinFunc> const & rhs) const
    { return this->value() < rhs.value(); }


    std::vector<std::deque<ImageID>> const & index() const { return _idx; }


    void update(std::deque<State3rd<BinFunc>>& dst)
    {
        std::size_t i = 0;
        for(bool e: _1st._remain){
            if(e){
                State3rd<BinFunc> dup = *this;
                dup.insert(i);
                dst.push_back(std::move(dup));
            }
            ++i;
        }
    }


    bool isEnd() const {
        return _ctIdx == (_1st._pb->div_x() - 1) * (_1st._pb->div_y() - 1);
    }


    State1st<BinFunc> _1st;
    std::vector<std::deque<ImageID>> _idx;
    std::size_t _cntLN;
    std::size_t _ctIdx;


    void insert(std::size_t i){
        const ImageID index = convToImageID(i, _1st._pb->div_x());

        auto pos = nowPos();
        ImageID tIh, tIv;

        Direction dir = Direction::left;
        if(pos[1] > _cntLN){
            tIh = _idx[pos[0]].back();
            tIv = _idx[pos[0]-1][pos[1]];
            dir = Direction::right;
            _idx[pos[0]].push_back(index);
        }else{
            tIh = _idx[pos[0]].front();
            tIv = _idx[pos[0]-1][pos[1]];
            _idx[pos[0]].push_front(index);
        }

        _1st._remain[i] = false;

        _1st._ev += std::abs((*_1st._pred)(_1st._pb->get_element(tIh),
                                  _1st._pb->get_element(index),
                                  dir));

        _1st._ev += std::abs((*_1st._pred)(_1st._pb->get_element(tIv),
                                 _1st._pb->get_element(index),
                                 Direction::down));

        ++_ctIdx;
    }


    Index2D nowPos() const {
        const size_t r = 1 + _ctIdx / (_1st._pb->div_x() - 1);
        size_t c = _ctIdx % (_1st._pb->div_x() - 1);

        if(c >= _cntLN)
            ++c;
        else
            c = _cntLN - c - 1;

        return makeIndex2D(r, c);
    }
};


template <typename BinFunc>
std::vector<std::vector<ImageID>> bfs_guess(utils::Problem const & pb, BinFunc f)
{
    // stage1
    std::vector<bool> rem(pb.div_y() * pb.div_x(), true);

    std::deque<State1st<BinFunc>> state1;
    for(auto i: utils::iota(rem.size())){
        auto remdup = rem;
        remdup[i] = false;
        state1.emplace_back(&pb, &f, std::deque<ImageID>({ convToImageID(i, pb.div_x()) }), remdup, 0.0);
    }

    writeln("Stage1");
    bfs_guess_impl(state1, 128);

    // stage2
    std::deque<State2nd<BinFunc>> state2;
    for(State1st<BinFunc>& e: state1)
        state2.emplace_back(std::move(e));

    writeln("Stage2");
    bfs_guess_impl(state2, 128);

    // stage3
    std::deque<State3rd<BinFunc>> state3;
    for(State2nd<BinFunc>& e: state2)
        state3.emplace_back(std::move(e));

    writeln("Stage3");
    bfs_guess_impl(state3, 128);

    if(state3.empty())
        return guess::guess(pb, f);
    else{
        auto mat = state3[0].index();
        std::vector<std::vector<ImageID>> dst;
        for(auto& v: mat)
            dst.emplace_back(v.begin(), v.end());

        return dst;
    }
}


template <typename BinFunc>
std::vector<std::vector<ImageID>> bfs_guess_parallel(utils::Problem const & pb, BinFunc f)
{
    const std::size_t allTileN = pb.div_y() * pb.div_x();
    const std::size_t threadN = static_cast<std::size_t>(std::floor(allTileN / (allTileN >= 64 ? std::sqrt(pb.div_y()) : 1) / (allTileN >= 144 ? std::sqrt(pb.div_y()) : 1)));

    // マルチスレッドで`bfs_guess_impl`を呼ぶ
    auto parallel_guess_impl = [&pb](auto& states){
        std::vector<std::thread> ths;
        for(auto& e: states)
            ths.emplace_back([&](){ bfs_guess_impl(e, 4096 / pb.div_x() / pb.div_y()); });

        std::size_t cnt = 0;
        for (auto& e : ths){
            e.join();
            utils::writeln("end: ", ++cnt);
        }
    };


    std::vector<bool> rem(pb.div_x() * pb.div_y(), true);

    // stage1
    std::vector<std::deque<State1st<BinFunc>>> state1(threadN);
    for(auto i: utils::iota(state1.size())){
        auto remdup = rem;
        remdup[i] = false;
        state1[i].emplace_back(&pb, &f, std::deque<ImageID>({ convToImageID(i, pb.div_x()) }), remdup, 0.0);
    }


    utils::writeln("Stage1");
    parallel_guess_impl(state1);

    // stage2
    std::vector<std::deque<State2nd<BinFunc>>> state2; state2.reserve(state1.size());
    for(auto& q: state1){
        std::deque<State2nd<BinFunc>> qq;
        for(auto& e: q)
            qq.emplace_back(std::move(e));

        state2.emplace_back(std::move(qq));
    }

    writeln("Stage2");
    parallel_guess_impl(state2);

    // stage3
    std::vector<std::deque<State3rd<BinFunc>>> state3; state3.reserve(state2.size());
    for(auto& q: state2){
        std::deque<State3rd<BinFunc>> qq;
        for (auto& e : q){
            if (e.isEnd())
                qq.emplace_back(std::move(e));
        }

        state3.emplace_back(std::move(qq));
    }

    writeln("Stage3");
    parallel_guess_impl(state3);


    // もっとも良い結果の選択
    State3rd<BinFunc> *minState = nullptr;
    for(auto& q: state3)
        for(auto& e: q)
            if(minState == nullptr || e < *minState)
                minState = &e;

    if(minState == nullptr){
        std::cout << "empty !!!" << std::endl;
        return guess::guess(pb, f);
    }else{
        auto mat = minState->index();
        std::vector<std::vector<ImageID>> dst;
        for(auto& v: mat)
            dst.emplace_back(v.begin(), v.end());

        return dst;
    }
}


/**
ある画像img1に対して、方角directionに画像img2がどの程度相関があるかを返します。
相関があるほど返す値は絶対値が小さくなります。
また、返す値は必ず正です。
*/
template <typename T, typename U
#ifdef SUPPORT_TEMPLATE_CONSTRAINTS
    , PROCON_TEMPLATE_CONSTRAINTS(utils::is_image<T>() && utils::is_image<U>())    // T, Uともに画像であるという制約
#endif
>
double diff_connection(T const & img1, U const & img2, utils::Direction dir)
{
    const double value = guess::diff_connection(img1, img2, dir);

    if(dir == Direction::up || dir == Direction::down)
        return value * img1.width();
    else
        return value * img1.height();
}

}}