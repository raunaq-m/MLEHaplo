FROM ubuntu:14.04

# install base packages
RUN apt-get update && apt-get install -y wget git g++ make zlib1g-dev bc

# Change workdir
WORKDIR /opt
RUN git clone https://github.com/raunaq-m/MLEHaplo.git

# Install multi-dsk
WORKDIR /opt/MLEHaplo/multi-dsk
RUN make k=64

# Install packages for perl scripts
RUN apt-get install -y curl
RUN curl -L https://cpanmin.us | perl - --sudo App::cpanminus
RUN cpanm Graph
RUN apt-get install -y libexpat1-dev
RUN cpanm XML::Parser
RUN cpanm Bio::SeqIO
# RUN cpanm --force Net::SSLeay
# RUN cpanm IO::Socket::SSL
# RUN cpanm LWP::Protocol::https
# RUN cpanm Bio::DB::GenBank
# RUN cpanm Bio::DB::GenPept
# RUN cpanm --force Bio::Perl


WORKDIR /opt/MLEHaplo

COPY Dockerfile /opt/
LABEL maintainer="li.weiling.1112@gmail.com"
LABEL author="raunaq.123@gmail.com"

