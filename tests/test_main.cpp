//
// Created by huangcanjia on 25-11-13.
//

#include <gtest/gtest.h>
#include "test_environment_geogram.h"
#include "test_environment_spdlog.h"

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);

    testing::AddGlobalTestEnvironment(new SpdlogTestEnvironment);
    testing::AddGlobalTestEnvironment(new GeogramTestEnvironment);

    return RUN_ALL_TESTS();
}