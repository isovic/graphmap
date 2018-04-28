FROM ubuntu as builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/graphmap
COPY . .
RUN make && make install
ENTRYPOINT ["/usr/bin/graphmap"]

# Second stage -- smaller image (( 400MB before --> 117MB now ))
FROM ubuntu

RUN apt-get update && apt-get install -y --no-install-recommends \
        zlib1g \
        libgomp1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data
COPY --from=builder /usr/bin/graphmap /usr/bin/graphmap
ENTRYPOINT ["/usr/bin/graphmap"]

