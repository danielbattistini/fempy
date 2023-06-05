// the function using the function pointers:
#include "Riostream.h"

void somefunction(void (*fptr)(void*, int, int), void* context) {
    fptr(context, 17, 42);
}

void non_member(void*, int i0, int i1) {
    std::cout << "I don't need any context! i0=" << i0 << " i1=" << i1 << "\n";
}

struct foo {
    void member(int i0, int i1) {
        std::cout << "member function: this=" << this << " i0=" << i0 << " i1=" << i1 << "\n";
    }
};

void forwarder(void* context, int i0, int i1) {
    static_cast<foo*>(context)->member(i0, i1);
}

void pointerToMemberFunction() {
    somefunction(&non_member, nullptr);
    foo object;
    somefunction(&forwarder, &object);
}