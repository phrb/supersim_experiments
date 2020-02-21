FROM debian
LABEL name="supersim"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential                                          \
        vim                                                      \
        htop                                                     \
        cmake                                                    \
        git                                                      \
        g++                                                      \
        python3                                                  \
        python3-pip                                              \
        python3-dev                                              \
        python3-tk                                               \
        libtclap-dev                                             \
        zlib1g-dev                                               \
        wget                                                     \
        python3-setuptools                                       \
        curl                                                     \
        ca-certificates                                          \
        r-base &&                                                \
        rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("rsm", "dplyr", "DiceKriging", "DiceDesign", "DiceOptim","randtoolbox", "future.apply"), repos="https://cran.rstudio.com")'

RUN pip3 install     \
        'numpy'      \
        'matplotlib' \
        'psutil'

RUN ln -s /usr/bin/python3 /usr/bin/python & \
    ln -s /usr/bin/pip3 /usr/bin/pip

ENV PATH /root/.local/bin\:$PATH

RUN git clone -b supersim_updated_install --single-branch https://github.com/phrb/supersim.git
RUN cd supersim/scripts && chmod +x installpy && ./installpy clone install

RUN git clone https://github.com/nicmcd/make-c-cpp ~/.makeccpp
RUN cd ~/.makeccpp && make

RUN cd supersim/scripts && chmod +x installcc && ./installcc clone build

# ENV http_proxy "http://web-proxy-pa.labs.hpecorp.net:8088/"
# ENV https_proxy "http://web-proxy-pa.labs.hpecorp.net:8088/"
