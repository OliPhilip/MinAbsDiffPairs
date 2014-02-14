///
/// Finding all minimum absolute distance pairs in an integers array
///
/// LinkedIn discussion in Algorithms and Data Structures Development's group started by Ramkumar S.
///
/// All codes are done from scratch except Mersenne twister borrowed from
/// http://my.opera.com/metrallik/blog/2013/04/19/c-class-for-random-generation-with-mersenne-twister-method
///

#include "algo_datastruct.h"

#ifndef cf_chrono
    #define cf_chrono
    #include <chrono>
#endif
#ifndef cf_iostream
    #define cf_iostream
    #include <iostream>
#endif

using namespace std;
using namespace algo;
using namespace data_struct;
int main()
{
    clock_t cpu_time;
    string _na;
    uint32_t *ans, i, nb = 50000, maxloop = 8;

    for (i = 0; i < maxloop; ++i) {
        _na = "RandomArray ";
        cpu_time = clock();
        ans = RandomArray<uint32_t>(nb, 0xFFFFFFFF);
        cpu_time = clock() - cpu_time;
        cout << _na << cpu_time << endl;

        _na = "MinDiffPairs_radix ";
        cpu_time = clock();
        MinDiffPairs_radix(ans, nb);
        cpu_time = clock() - cpu_time;
        cout << _na << cpu_time << endl;

        _na = "MinDiffPairs_quick ";
        cpu_time = clock();
        MinDiffPairs_quick(ans, nb);
        cpu_time = clock() - cpu_time;
        cout << _na << cpu_time << endl;

        _na = "MinDiffPairs_hash ";
        cpu_time = clock();
        MinDiffPairs_hash(ans, nb);
        cpu_time = clock() - cpu_time;
        cout << _na << cpu_time << endl;

        free(ans);
    }

    return 0;
}
