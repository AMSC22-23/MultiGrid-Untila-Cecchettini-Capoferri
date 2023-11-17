#ifndef CLASSES_H
#define CLASSES_H

template<typename T_m, typename T_v>
class Iteration{
    public:
        virtual Iteration(){

        }
    protected:
        T_m* m_A;
        T_v* m_b;
        size_t m_size;
};

#endif