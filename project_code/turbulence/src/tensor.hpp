#include <iostream>
#include <vector>

#ifndef TENSOR
#define TENSOR

template <typename T>
class Tensor
{
    std::vector<T> components;

    public:
        Tensor();
        Tensor(const std::vector<T>& components);

        void push_component(const T& component);

        std::vector<T> get_components();
        T get_component(int index);

        // type specific
        Tensor<T> transpose();
        void add(Tensor<T>& tensor);
        void scale_by(double scaling_factor);
        double double_dot(Tensor<T>& tensor);
};



template <typename T>
Tensor<T>::Tensor()
{}

template <typename T>
Tensor<T>::Tensor(const std::vector<T>& components)
{
    this->components = components;
}



template <typename T>
void Tensor<T>::push_component(const T& comp)
{
    components.push_back(comp);
}



template <typename T>
std::vector<T> Tensor<T>::get_components()
{
    return components;
}

template <typename T>
T Tensor<T>::get_component(int index)
{
    return components[index];
}



// type specific implementations in cpp file



#endif