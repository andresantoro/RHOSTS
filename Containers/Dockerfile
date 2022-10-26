# Dockerfile, Image, Container
FROM openjdk:slim
COPY --from=python:3.9 / /

ENV PATH "$PATH:/opt/jython/bin"


RUN wget https://repo1.maven.org/maven2/org/python/jython-installer/2.7.3/jython-installer-2.7.3.jar

RUN mkdir /opt/jython && mkdir /repo
RUN java -jar jython-installer-2.7.3.jar -s -d /opt/jython/


RUN pip install numpy>=1.16.0 scipy>=1.5.0 persim>=0.3.1 pybind11>=2.5.0 h5py>=2.10.0

RUN pip install --use-deprecated=legacy-resolver cechmate>=0.0.10

