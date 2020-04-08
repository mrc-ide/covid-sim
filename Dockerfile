FROM debian:stable AS build

ARG COUNTRY=UK

WORKDIR /src

COPY . .

RUN mkdir -p build \
    && ./ci/install_dependencies.sh

WORKDIR /src/build

RUN cmake -DCOUNTRY=${COUNTRY} ../src \
    && make

# This allows for building a release without having to potentially re-run the tests.
FROM debian:stable-slim AS release

# Install runtime dependencies.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 \
    && apt-get clean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /src/build/SpatialSim /usr/bin/.
COPY data /data/input

ENTRYPOINT ["/usr/bin/SpatialSim"]

FROM build AS test

WORKDIR /src/tests

RUN pwd \
    && ls -l \
    && ./regressiontest_UK_100th.py

# This ensures that the default behavior is to run the tests and then create a release
FROM release
