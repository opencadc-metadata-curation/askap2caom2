FROM opencadc/astropy

RUN pip install astropy && pip install numpy && \
        pip install spherical-geometry

RUN pip install cadcdata && \
        pip install caom2 && \
        pip install caom2repo && \
        pip install caom2utils && \
        pip install PyYAML && pip install vos && \
        pip install caom2

RUN pip install pytest && pip install mock && pip install flake8 && \
        pip install funcsigs && pip install xml-compare && \
        pip install pytest-cov && pip install aenum && pip install future

WORKDIR /usr/src/app
RUN git clone https://github.com/opencadc-metadata-curation/caom2pipe.git && \
  pip install ./caom2pipe
  
RUN git clone https://github.com/opencadc-metadata-curation/askap2caom2.git && \
  pip install ./askap2caom2

COPY ./docker-entrypoint.sh ./

ENTRYPOINT ["./docker-entrypoint.sh"]

