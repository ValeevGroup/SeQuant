#!/usr/bin/env bash

cmake --build ../build --target srcc -- -j7 && ../build/srcc
