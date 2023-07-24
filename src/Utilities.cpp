#include "../include/Utilities.hpp"

string Timer::get_elapsed_time() const
{
    using std::chrono::duration_cast;

    Time end_time = Clock::now();
    size_t total_miliseconds = duration_cast<Milliseconds>(end_time - start_time).count();
    size_t hours = total_miliseconds / (1000 * 60 * 60);
    size_t minutes = (total_miliseconds / (1000 * 60)) % 60;
    size_t seconds = (total_miliseconds / 1000) % 60;
    size_t milliseconds = total_miliseconds & 1000;

    using namespace std;
    constexpr size_t WIDTH = 4;
    std::stringstream ss;
    ss << "\n"
       << "    Hours:" << setw(WIDTH + 2) << hours << endl
       << "    Minutes:" << setw(WIDTH) << minutes << endl
       << "    Seconds:" << setw(WIDTH) << seconds << "." << milliseconds << endl;

    return ss.str();
};