FROM debian:stable AS build

WORKDIR /src

COPY . .

RUN mkdir -p build \
    && ./ci/install_dependencies.sh

WORKDIR /src/build

RUN cmake -DENABLE_REGRESSION_TEST:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .. \
    && cmake --build . -- -j$(nproc)

# This allows for building a release without having to potentially re-run the tests.
FROM debian:stable-slim AS release

# Install runtime dependencies.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 \
    && apt-get clean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /src/build/src/CovidSim /usr/bin/.
COPY data /data/input

ENTRYPOINT ["/usr/bin/CovidSim"]

FROM build AS test

WORKDIR /src/build

RUN ctest --tests-regex regressiontest_UK_100th --verbose

# This ensures that the default behavior is to run the tests and then create a release
FROM release
