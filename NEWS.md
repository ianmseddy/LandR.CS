Known issues: https://github.com/ianmseddy/LandR.CS/issues

version 0.0.1
=============

## bug fixes

* `LandR.CS::own()` no longer uses `grep(..., sys.calls())` which takes many many minutes if `...` is a simList. Use `reproducible::.grepSysCalls` which is for this purpose.


