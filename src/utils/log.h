//
// Created by huangcanjia on 25-11-27.
//

#ifndef LOG_H
#define LOG_H

#include <spdlog/spdlog.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/sinks/sink.h>

#ifdef ERROR
#undef ERROR
#endif

namespace LOG {
    namespace detail {
        class SourceLocation {
        public:
            constexpr SourceLocation(const char *fileName = __builtin_FILE(),
                                     const char *funcName = __builtin_FUNCTION(),
                                     std::uint32_t lineNum = __builtin_LINE()) noexcept
                    : _fileName(fileName), _funcName(funcName), _lineNum(lineNum) {
            }

            [[nodiscard]] constexpr const char *FileName() const noexcept {
                return _fileName;
            }

            [[nodiscard]] constexpr const char *FuncName() const noexcept {
                return _funcName;
            }

            [[nodiscard]] constexpr std::uint32_t LineNum() const noexcept {
                return _lineNum;
            }

        private:
            const char *_fileName;
            const char *_funcName;
            const std::uint32_t _lineNum;
        };

        constexpr auto GetLogSourceLocation(const SourceLocation &location) {
            return spdlog::source_loc{location.FileName(), static_cast<int>(location.LineNum()),
                                      location.FuncName()};
        }

    } // namespace detail

    // trace
    template<typename... Args>
    struct TRACE {
        explicit constexpr TRACE(fmt::format_string<Args...> fmt, Args &&...args, const detail::SourceLocation &location = {}) {
            spdlog::log(GetLogSourceLocation(location), spdlog::level::trace, fmt, std::forward<Args>(args)...);
        }
    };

    template<typename... Args>
    TRACE(fmt::format_string<Args...> fmt, Args &&...args) -> TRACE<Args...>;

    // debug
    template<typename... Args>
    struct DEBUG {
        explicit constexpr DEBUG(fmt::format_string<Args...> fmt, Args &&...args, const detail::SourceLocation &location = {}) {
            spdlog::log(GetLogSourceLocation(location), spdlog::level::debug, fmt, std::forward<Args>(args)...);
        }
    };

    template<typename... Args>
    DEBUG(fmt::format_string<Args...> fmt, Args &&...args) -> DEBUG<Args...>;

    // info
    template<typename... Args>
    struct INFO {
        explicit constexpr INFO(fmt::format_string<Args...> fmt, Args &&...args, const detail::SourceLocation &location = {}) {
            spdlog::log(GetLogSourceLocation(location), spdlog::level::info, fmt, std::forward<Args>(args)...);
        }
    };

    template<typename... Args>
    INFO(fmt::format_string<Args...> fmt, Args &&...args) -> INFO<Args...>;

    // warn
    template<typename... Args>
    struct WARN {
        explicit constexpr WARN(fmt::format_string<Args...> fmt, Args &&...args, const detail::SourceLocation &location = {}) {
            spdlog::log(GetLogSourceLocation(location), spdlog::level::warn, fmt, std::forward<Args>(args)...);
        }
    };

    template<typename... Args>
    WARN(fmt::format_string<Args...> fmt, Args &&...args) -> WARN<Args...>;

    // error
    template<typename... Args>
    struct ERROR {
        explicit constexpr ERROR(fmt::format_string<Args...> fmt, Args &&...args, const detail::SourceLocation &location = {}) {
            spdlog::log(GetLogSourceLocation(location), spdlog::level::err, fmt, std::forward<Args>(args)...);
        }
    };

    template<typename... Args>
    ERROR(fmt::format_string<Args...> fmt, Args &&...args) -> ERROR<Args...>;

    // CRITICAL
    template<typename... Args>
    struct CRITICAL {
        explicit constexpr CRITICAL(fmt::format_string<Args...> fmt, Args &&...args, const detail::SourceLocation &location = {}) {
            spdlog::log(GetLogSourceLocation(location), spdlog::level::critical, fmt, std::forward<Args>(args)...);
        }
    };

    template<typename... Args>
    CRITICAL(fmt::format_string<Args...> fmt, Args &&...args) -> CRITICAL<Args...>;
} // namespace LOG

#endif //LOG_H
