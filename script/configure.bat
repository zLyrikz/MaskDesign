@echo off
setlocal

rem 设置目标目录
set "release_dir=..\bin\Release"
set "debug_dir=..\bin\Debug"

rem 创建目标目录（如果不存在）
if not exist "%release_dir%" mkdir "%release_dir%"
if not exist "%debug_dir%" mkdir "%debug_dir%"

rem 设置源文件夹路径
set "src_dir=..\build\_deps\libigl-src\external\mpfr\lib"
set "src_dirr=..\build\_deps\libigl-src\external\gmp\lib"

rem 移动 libgmp*.lib 文件到 Release 目录
copy "%src_dir%\libmpfr*.dll" "%release_dir%"
copy "%src_dirr%\libgmp*.dll" "%release_dir%"

rem 移动 libgmp*.lib 文件到 Debug 目录
copy "%src_dir%\libmpfr*.dll" "%debug_dir%"
copy "%src_dirr%\libgmp*.dll" "%debug_dir%"

echo Complete

endlocal


