@echo off

::repeat runs
for /l %%r in (1,1,1) do (
::batches
for /l %%x in (1,1,1) do (
::size
for /l %%y in (1,1,3) do (
	::calculate a = 2^x
	set a=1
	for /l %%i in (1,1,%%x%%) do set /a a*=2
	::calculate b = 2^y
	set b=1
	for /l %%i in (1,1,%y%) do set /a b*=2

	::echo %a% batches and %b% size
	echo "x64/Release/LP.exe" benchmarks/%b% 
)
)
)
	::"x64/Release/LP.exe" benchmarks/%%y %x%
