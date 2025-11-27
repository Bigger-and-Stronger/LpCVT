//
// Created by huangcanjia on 25-11-13.
//

#ifndef TEST_ENVIRONMENT_GEOGRAM_H
#define TEST_ENVIRONMENT_GEOGRAM_H

#include <geogram/basic/command_line_args.h>
#include <geogram/basic/common.h>
#include <geogram/basic/attributes.h>
#include <gtest/gtest.h>
#include "base/base.hpp"

class GeogramTestEnvironment final : public testing::Environment {
public:
    void SetUp() override {
        GEO::initialize(GEO::GEOGRAM_INSTALL_ALL);
        GEO::CmdLine::import_arg_group("standard");
        GEO::geo_register_attribute_type<hwpcb::euid>("euid");
    }
};

#endif //TEST_ENVIRONMENT_GEOGRAM_H
