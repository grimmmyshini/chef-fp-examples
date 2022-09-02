FROM ubuntu:bionic

RUN apt -y update
RUN apt -y upgrade
RUN apt install -y build-essential
RUN apt install -y wget cmake gcc git python3 python3-pip curl

WORKDIR /code
RUN wget https://raw.githubusercontent.com/grimmmyshini/clad-fp-error-est-examples/main/setup.sh
RUN chmod +x setup.sh
RUN ./setup.sh -s
RUN ./setup.sh -r
RUN cat results.txt | curl -F 'sprunge=<-' http://sprunge.us
