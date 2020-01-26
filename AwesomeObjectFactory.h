//
// Created by s2 on 2020-01-25.
//

#ifndef TESTSET_AWESOMEOBJECTFACTORY_H
#define TESTSET_AWESOMEOBJECTFACTORY_H

#include <vector>
#include <map>
#include <stdexcept>
#include <iostream>

namespace icy {
    template<class T>
    class AwesomeObjectFactory;
}


template<class T>
class icy::AwesomeObjectFactory {

public:
    AwesomeObjectFactory(int initialSize);
    ~AwesomeObjectFactory();
    T* take();
    void release(T* obj);
    void releaseAll();
    void printout(); // for testing

private:
    std::vector<T*> available;           // any item on the stack is free to check out
    std::map<T*, bool> registry;       // keep track of all objects; true == 'taken'
};

template<class T>
icy::AwesomeObjectFactory<T>::AwesomeObjectFactory(int initialSize)
{
    for(int i=0;i<initialSize;i++) {
        T* obj = new T;
        registry[obj] = false;
        available.push_back(obj);
    }
}

template<class T>
icy::AwesomeObjectFactory<T>::~AwesomeObjectFactory()
{
    for(auto &x : registry)
        delete x.first;
}

template <class T>
T* icy::AwesomeObjectFactory<T>::take()
{
    T *obj;
    if(available.size() == 0) {
        obj = new T;
    }
    else {
        obj = available.back();
        available.pop_back();
    }

    registry[obj] = true;
    return obj;
}

template<class T>
void icy::AwesomeObjectFactory<T>::release(T* obj)
{
    auto it = registry.find(obj);
    if(it == registry.end())
        throw std::runtime_error("trying to release non-allocated object");
    else {
        registry[obj] = false;
        available.push_back(obj);
    }
}

template<class T>
void icy::AwesomeObjectFactory<T>::releaseAll()
{
    available.clear();
    for(auto &x : registry) {
        T* key = x.first;
        available.push_back(key);
        x.second = false;
    }
}

template<class T>
void icy::AwesomeObjectFactory<T>::printout()
{
    std::cout << "available.size " << available.size() << std::endl;
    std::cout << "registry: " << std::endl;
    for(auto const &x : registry) {
        T* key = x.first;
        bool val = x.second;
        std::cout << "key " << key << "; value " << val << std::endl;
    }
    std::cout << std::endl;
}

#endif //TESTSET_AWESOMEOBJECTFACTORY_H
