FROM mdstudio/mdstudio_docker2:0.0.1

# Install requirements
RUN apt-get update -y && apt-get install swig libopenbabel-dev openbabel -y

COPY . /home/mdstudio/lie_structures

RUN chown mdstudio:mdstudio /home/mdstudio/lie_structures

WORKDIR /home/mdstudio/lie_structures

RUN pip install "https://github.com/cinfony/cinfony/tarball/master#egg=cinfony-1.2"

RUN pip install .

USER mdstudio

CMD ["bash", "entry_point_lie_structures.sh"]
