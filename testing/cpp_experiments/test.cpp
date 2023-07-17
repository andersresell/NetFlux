#include <iostream>
#include <cassert>
#include <limits>
using namespace std;

#include <iostream>
#include <vector>
#include <algorithm>

struct MyType
{
    int value;
    // Add comparison operator for sorting
    bool operator<(const MyType &other) const
    {
        return value < other.value;
    }
};

int main()
{
    std::vector<MyType> guidingVector = {{4}, {2}, {1}, {3}};
    std::vector<int> otherVector1 = {40, 20, 10, 30};
    std::vector<char> otherVector2 = {'D', 'B', 'A', 'C'};

    // Create an index vector to store the original indices
    std::vector<size_t> indices(guidingVector.size());
    for (size_t i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    // Sort the indices based on the guiding vector
    std::sort(indices.begin(), indices.end(), [&guidingVector](size_t a, size_t b)
              { return guidingVector[a] < guidingVector[b]; });

    // Rearrange the other vectors using the sorted indices
    std::vector<int> sortedOtherVector1(guidingVector.size());
    std::vector<char> sortedOtherVector2(guidingVector.size());

    for (size_t i = 0; i < guidingVector.size(); ++i)
    {
        sortedOtherVector1[i] = otherVector1[indices[i]];
        sortedOtherVector2[i] = otherVector2[indices[i]];
    }

    // Output the sorted vectors
    for (const auto &element : guidingVector)
    {
        std::cout << element.value << " ";
    }
    std::cout << std::endl;

    for (const auto &element : sortedOtherVector1)
    {
        std::cout << element << " ";
    }
    std::cout << std::endl;

    for (const auto &element : sortedOtherVector2)
    {
        std::cout << element << " ";
    }
    std::cout << std::endl;

    return 0;
}
