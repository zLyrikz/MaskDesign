include(FetchContent)
FetchContent_Declare(
    finite-diff
    GIT_REPOSITORY https://github.com/zfergus/finite-diff.git
    GIT_TAG v1.0.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(finite-diff)