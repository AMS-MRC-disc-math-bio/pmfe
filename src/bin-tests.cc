// Copyright (c) 2015 Andrew Gainer-Dewar.

// This is a magic file to build the testing binary.
// All the work is done in catch.hpp.
// All test cases are defined in other files ("test-*.cc") to keep compilation zippy.
// There should never be any need to modify this file.

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#define BOOST_LOG_DYN_LINK 1 // Fix an issue with dynamic library loading
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

int main( int argc, char* const argv[] )
{
    // Setup
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::warning);

    // Run tests
    int result = Catch::Session().run( argc, argv );

    // Cleanup

    // Finish
    return result;
}
