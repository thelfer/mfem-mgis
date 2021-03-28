# /usr/bin/env bash

if command -v clang-format &> /dev/null
then
  clang-format -i $(find . -name "*xx")
fi
