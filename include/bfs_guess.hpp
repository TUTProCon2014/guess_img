#pragma once

#include "../../inout/include/inout.hpp"
#include "../../utils/include/template.hpp"
#include "../../utils/include/types.hpp"
#include "../../utils/include/range.hpp"
#include "../../utils/include/dwrite.hpp"

#include <vector>
#include <set>
#include <array>
#include <deque>
#include <deque>
#include <tuple>

namespace procon { namespace bfs_guess {

using namespace utils;


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
void bfs_guess_impl(std::deque<State>& state)
{
    std::size_t iterCount = 0;
    while(!state.empty() && !state[0].isEnd()){
        writeln(std::cout, "iteration: ", iterCount++);

        std::deque<State> dst1;
        for(auto& e: state){
            e.update(dst1);

            std::deque<State> dst2;
            refine(std::move(dst1), dst2, [](double v, double m, double sd){ return v <= m + sd; });
            dst1 = [](){ std::deque<State> ini; return ini; }();

            std::sort(dst2.begin(), dst2.end(), [](State & a, State & b){ return a.value() < b.value(); });
            while(!dst2.empty() && dst1.size() < 128){
                if(dst1.empty() || dst1[dst1.size()-1] != dst2[0])
                    dst1.push_back(std::move(dst2[0]));
                dst2.pop_front();
            }
        }

        state = std::move(dst1);
    }
}



template <typename BinFunc>
struct State1st
{
    State1st(utils::Problem const *pb, BinFunc *func, std::deque<Index2D> const & idx, std::set<Index2D> const & rem, double ev)
    : _pb(pb), _pred(func), _idx(idx), _remain(rem), _ev(ev) {}


    double value() const
    {
        return _ev;
    }


    std::deque<Index2D> const & index() const
    {
        return _idx;
    }


    bool operator==(State1st<BinFunc> const & rhs) const
    {
        return utils::equal(index(), rhs.index());
    }


    bool operator!=(State1st<BinFunc> const & rhs) const
    {
        return !(*(this) == rhs);
    }


    bool isEnd() const { return _idx.size() == _pb->div_y(); }


    void update(std::deque<State1st> & q)
    {
        for(auto& e: _remain){
            State1st<BinFunc> dupTop = *this;
            dupTop.insert(Direction::up, e);
            q.push_back(std::move(dupTop));

            State1st<BinFunc> dupBottom = *this;
            dupBottom.insert(Direction::down, e);
            q.push_back(std::move(dupBottom));
        }
    }


    Problem const *_pb;
    BinFunc *_pred;
    std::deque<Index2D> _idx;
    std::set<Index2D> _remain;
    double _ev;


    double pred_value(Direction dir, Index2D const & index)
    {
        Index2D tI;
        if(dir == Direction::up)
            tI = _idx[0];
        else
            tI = _idx[_idx.size()-1];

        return std::abs((*_pred)(_pb->get_element(tI[0], tI[1]),
                        _pb->get_element(index[0], index[1]),
                        dir));
    }


    void insert(Direction dir, Index2D const & index)
    {
        _ev += pred_value(dir, index);

        if(dir == Direction::up)
            _idx.push_front(index);
        else
            _idx.push_back(index);
        _remain.erase(index);
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


    bool isEnd() const { return _idx.size() == _1st._pb->div_x(); }


    void update(std::deque<State2nd>& dst)
    {
        for(auto& e: _1st._remain){
            auto dupLeft = *this;
            dupLeft.insert(Direction::left, e);
            dst.push_back(std::move(dupLeft));

            auto dupRight = *this;
            dupRight.insert(Direction::right, e);
            dst.push_back(std::move(dupRight));
        }
    }


    std::deque<Index2D> const & index2nd()
    {
        return _idx;
    }


    std::deque<Index2D> const & index1st()
    {
        return _1st.index();
    }


    std::size_t insertLeftCount() const { return _cntLN; }


    State1st<BinFunc> _1st;
    std::deque<Index2D> _idx;
    std::size_t _cntLN;


    double pred_value(Direction dir, Index2D const & index)
    {
        auto tI = _idx[0];
        if(dir == Direction::right)
            tI = _idx[_idx.size() - 1];

        return std::abs((*_1st._pred)(_1st._pb->get_element(tI[0], tI[1]),
                                      _1st._pb->get_element(index[0], index[1]),
                                      dir));
    }


    void insert(Direction dir, Index2D const & index)
    {
        _1st._ev += pred_value(dir, index);

        if(dir == Direction::left){
            _idx.push_front(index);
            ++_cntLN;
        }
        else
            _idx.push_back(index);

        _1st._remain.erase(index);
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


    std::vector<std::deque<Index2D>> const & index() const { return _idx; }


    void update(std::deque<State3rd<BinFunc>>& dst)
    {
        for(auto& e: _1st._remain){
            State3rd<BinFunc> dup = *this;
            dup.insert(e);
            dst.push_back(std::move(dup));
        }
    }


    bool isEnd() const {
        return _ctIdx == (_1st._pb->div_x() - 1) * (_1st._pb->div_y() - 1);
    }


    State1st<BinFunc> _1st;
    std::vector<std::deque<Index2D>> _idx;
    std::size_t _cntLN;
    std::size_t _ctIdx;


    void insert(Index2D const & index){
        auto pos = nowPos();
        Index2D tIh, tIv;

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

        _1st._remain.erase(index);

        _1st._ev += std::abs((*_1st._pred)(_1st._pb->get_element(tIh[0], tIh[1]),
                                  _1st._pb->get_element(index[0], index[1]),
                                  dir));

        _1st._ev += std::abs((*_1st._pred)(_1st._pb->get_element(tIv[0], tIv[1]),
                                 _1st._pb->get_element(index[0], index[1]),
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
std::vector<std::vector<Index2D>> bfs_guess(utils::Problem const & pb, BinFunc f)
{
    // stage1
    std::set<Index2D> rem;
    for(auto i: utils::iota(pb.div_y()))
        for(auto j: utils::iota(pb.div_x()))
            rem.insert(makeIndex2D(i, j));

    std::deque<State1st<BinFunc>> state1;
    for(auto& e: rem){
        auto remdup = rem;
        remdup.erase(e);
        state1.emplace_back(&pb, &f, std::deque<Index2D>({ e }), remdup, 0.0);
    }

    writeln(std::cout, "Stage1");
    bfs_guess_impl(state1);

    // stage2
    std::deque<State2nd<BinFunc>> state2;
    for(State1st<BinFunc>& e: state1)
        state2.emplace_back(std::move(e));

    writeln(std::cout, "Stage2");
    bfs_guess_impl(state2);

    // stage3
    std::deque<State3rd<BinFunc>> state3;
    for(State2nd<BinFunc>& e: state2)
        state3.emplace_back(std::move(e));

    writeln(std::cout, "Stage3");
    bfs_guess_impl(state3);

    if(state3.empty())
        return guess::guess(pb, f);
    else{
        auto mat = state3[0].index();
        std::vector<std::vector<Index2D>> dst;
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