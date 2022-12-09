FROM dolfinx/lab

USER root
RUN pip install meshio h5py==3.6.0 pygmsh 

WORKDIR /home/app        

RUN git clone https://github.com/UCL/linearized-NSE-data-assimilation.git

WORKDIR /home/app/linearized-NSE-data-assimilation
