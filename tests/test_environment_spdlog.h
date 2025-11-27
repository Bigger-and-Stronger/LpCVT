//
// Created by huangcanjia on 25-11-13.
//

#ifndef TEST_ENVIRONMENT_SPDLOG_H
#define TEST_ENVIRONMENT_SPDLOG_H

#include <gtest/gtest.h>
#include "utils/log.h"

class SpdlogTestEnvironment final : public testing::Environment {
public:
    void SetUp() override {
        spdlog::set_level(spdlog::level::trace);
    }
};

#endif //TEST_ENVIRONMENT_SPDLOG_H
