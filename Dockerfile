FROM rootproject/root:6.32.04-ubuntu24.04


ENV GARFIELD_HOME=/opt/garfieldpp

WORKDIR /opt/garfield_build

RUN mkdir $GARFIELD_HOME && \
    git clone https://gitlab.cern.ch/garfield/garfieldpp.git && \
    cd garfieldpp && mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=${GARFIELD_HOME} .. && \
    make -j $(nproc) && make install && \
    cd ${GARFIELD_HOME} && rm -rf /opt/garfield_build

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${GARFIELD_HOME}/lib

CMD [ "bash" ]