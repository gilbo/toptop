// +-------------------------------------------------------------------------
// | maybe.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2012
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the TopTop library.
// |
// |    TopTop is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    TopTop is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with TopTop.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

template<class T>
class Maybe {
    // please use these aliases instead of constructor functions
    template<class U>
    friend Maybe<U> nothing();
    template<class U>
    friend Maybe<U> just(const U &);
public: // this must be public so that there's a default value for a Maybe
    Maybe() : present(false) {}
private:
    Maybe(const T& just) : present(true), data(just) {}
    
public:
    Maybe(const Maybe<T> &cp) : present(cp.present), data(cp.data) {}
    //Maybe(Maybe<T>&& other) : present(other.present), data(other.data) {}
    ~Maybe() {} // data destructor called automatically
    
    const T& just() const { return data; }
          T& just()       { return data; }
    
    explicit operator bool() const {
        return present;
    }
private:
    bool present;
    T data;
};

template<class T>
Maybe<T> nothing() {
    return Maybe<T>();
}

template<class T>
Maybe<T> just(const T &data) {
    return Maybe<T>(data);
}

template<class Iterator, class Function>
inline
void for_maybe(
    Iterator begin,
    Iterator end,
    Function work
) {
    for( ; begin != end; ++begin) {
        if(*begin)
            work(begin->just());
    }
}

