#!/bin/bash
#
# This script is intended to debug a segfault
# on Linux in CI environment where there is no
# way to run a program interactively ; it relies
# on GDB
#

echo "$(pwd)"

# check parameters
if [ -z "$1" ]
then
  echo "Usage: $0  <program command line>"
  exit 1
fi

# create a gdb command file to print backtrace
echo "backtrace
quit" > ./gdb_backtrace

# run faulty program
ulimit -c unlimited
"$@"
status="$?"

# run gdb to analyse core dump file
if [[ "$status" -eq 139 ]]
then
    gdb -q $1 core -x ./gdb_backtrace
fi

exit "$status"