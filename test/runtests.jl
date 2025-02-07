using Aqua
using TestItemRunner
using NetworkJumpProcesses

Aqua.test_all(NetworkJumpProcesses)
@run_package_tests verbose=true
