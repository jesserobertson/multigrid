/*
    main.cpp (Multigrid)
    Jesse Robertson, 2010-05-20

    Driver routine for multigrid solver
*/

//#include <dispatch/dispatch.h> // GCD Queue library
#include <multigrid/multigrid.hpp>
#include "mosolov.hpp"
#include "mosolov_settings.hpp"

using namespace mgrid;
using namespace std;
using namespace blitz;

typedef TinyVector<double, 2>::T_vector ABTuple;

inline double critical_bingham(const double aspect) {
    return (2 + aspect - sqrt(4 + aspect*aspect + (2*pi - 4)*aspect))/(4-pi);
}

void calculate_flow(ABTuple aspectBinghamPair) {
    // Define settings
    MosolovSettings settings;
    settings.multigridSettings.aspectRatio = aspectBinghamPair(0);
    settings.binghamNumber = aspectBinghamPair(1);
    
    // Construct and solve problem
    std::auto_ptr<Mosolov> problem(new Mosolov(settings));
    problem->solve();
    problem->write(3); // Write out velocity, strain rate and residual
}

// void calculate_lists() {
//     // = Calculation settings =
//     // Specify aspect ratios to iterate over
//     static const double aspectRatios[] =
//         {2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10,
//          10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5,
//          18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25};
     
//     // Specify number of Bingham gradations per aspect ratio. The calculation
//     // will not include B=0 or B=B* since these are known analytically.
//     static const int nBingham = 32;
//     static const double binghamFrac[nBingham] =
//         {0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 
//          0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 
//          0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};
    
//     // = Parallel main loop =
//     // Get the default concurrent job queue for the application, and define
//     // a group of jobs which can be executed
//     dispatch_queue_t queue \
//         = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
//     dispatch_group_t group = dispatch_group_create();
    
//     // Generate a list of aspect-bingham pairs
//     vector<ABTuple> aspectBinghamPairs;
//     static const double binghamFracs[] = {0.8, 0.85, 0.9, 0.95};
//     foreach(double binghamFrac, binghamFracs) {
//         ABTuple abPair(18, binghamFrac*critical_bingham(18));
//         dispatch_group_async(group, queue, ^{calculate_flow(abPair);});
//     }
    
//     // Wait for threads to finish, then clean up
//     dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
//     dispatch_release(group);
// }

// void calculate_spec_pairs() {
//     vector<ABTuple> aspectBinghamPairs;
//     aspectBinghamPairs.push_back(ABTuple(5.9407, 0.2533));
//     aspectBinghamPairs.push_back(ABTuple(5.9407, 0.3149));
//     aspectBinghamPairs.push_back(ABTuple(4.5970, 0.0));
//     aspectBinghamPairs.push_back(ABTuple(8.1437, 0.2533));
//     aspectBinghamPairs.push_back(ABTuple(8.1437, 0.3149));
//     aspectBinghamPairs.push_back(ABTuple(6.7490, 0.0));
//     aspectBinghamPairs.push_back(ABTuple(7.0422, 0.2841));
//     aspectBinghamPairs.push_back(ABTuple(5.673, 0.0));
    
//     // Get the default concurrent job queue for the application, and define
//     // a group of jobs which can be executed
//     dispatch_queue_t queue \
//         = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
//     dispatch_group_t group = dispatch_group_create();
    
//     // Main loop
//     foreach(ABTuple abPair, aspectBinghamPairs)
//         dispatch_group_async(group, queue, ^{calculate_flow(abPair);});
//     dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
//     dispatch_release(group);
// }

int main() {
    // calculate_spec_pairs();
    calculate_flow(ABTuple(2, 0.24));
}