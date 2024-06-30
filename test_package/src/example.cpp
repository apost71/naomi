#include "naomi.h"
#include <vector>
#include <string>

int main() {
    naomi();

    std::vector<std::string> vec;
    vec.push_back("test_package");

    naomi_print_vector(vec);
}
