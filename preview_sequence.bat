@echo off
setlocal EnableDelayedExpansion

set EXE=build\Release\AtmosphereBakePreview
set INDIR=.\out\am0_v2
set OUTDIR=.\out\preview_seq_v2

if not exist "%OUTDIR%" mkdir "%OUTDIR%"

for /L %%I in (0,1,190) do (
    set /A "ZI=%%I/2"
    set /A "ZF=(%%I %% 2)*5"

    set ZENITH=!ZI!.!ZF!

    set PAD=0000%%I
    set FRAME=!PAD:~-4!

    echo Rendering frame !FRAME!  sun-zenith=!ZENITH! deg

    "%EXE%" ^
      --in "%INDIR%" ^
      --out "%OUTDIR%\preview_!FRAME!.png" ^
      --width 1600 ^
      --height 800 ^
      --sun-zenith-deg !ZENITH! ^
      --sun-azimuth-deg 90 ^
      --exposure 0.2 ^
      --sky-mu 64 ^
      --sky-mu-s 64 ^
      --sky-nu 128

    if errorlevel 1 (
        echo Error: failed at frame !FRAME! ^(zenith=!ZENITH! deg^)
        exit /b 1
    )
)

echo Done.
endlocal