FROM rootproject/root:6.32.04-ubuntu24.04


ENV ROOTSYS=/opt/root/
ENV GARFIELD_HOME=/opt/garfieldpp
ENV LD_LIBRARY_PATH=/lib:/lib64:/usr/lib:/opt/root/lib:$LD_LIBRARY_PATH

WORKDIR /opt/

RUN git clone https://gitlab.cern.ch/garfield/garfieldpp.git && \
    cd garfieldpp && mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=${GARFIELD_HOME}/install .. && \
    make -j $(nproc) && make install && \
    mkdir ${GARFIELD_HOME}/install/lib64/ && \
    ln -s ${GARFIELD_HOME}/install/lib/* ${GARFIELD_HOME}/install/lib64/ 

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${GARFIELD_HOME}/install/lib
ENV GARFIELD_IONDATA=${GARFIELD_HOME}/Data

CMD [ "bash" ]