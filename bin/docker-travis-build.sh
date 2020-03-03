#!/bin/sh

# this script builds a 'Bionic' env docker image used by Travis-CI for ValeevGroup/SeQuant2 (aka SeQuant) project
#
# to run bash in the image: docker run -it sequant-travis-debug bash -l
# N.B. relevant locations:
#   - source dir: /home/travis/build/ValeevGroup/SeQuant (TRAVIS_BUILD_DIR env in Travis jobs)
#   - build dir: /home/travis/_build
#   - install dir: /home/travis/_install

# this is where in the container file system Travis-CI "starts"
export TRAVIS_BUILD_TOPDIR=/home/travis/build
export CURRENT_REV_SHA=`git rev-parse HEAD`

##############################################################
# make a script to download all prereqs and clone ValeevGroup/SeQuant2 repo
setup=setup.sh
cat > $setup << END
#!/bin/sh
curl -sSL "http://apt.llvm.org/llvm-snapshot.gpg.key" | apt-key add -
echo "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-9 main" | tee -a /etc/apt/sources.list > /dev/null
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB | sudo apt-key add -
echo "deb https://apt.repos.intel.com/tbb all main" | tee -a /etc/apt/sources.list.d/intel-tbb.list > /dev/null
apt-add-repository -y "ppa:ubuntu-toolchain-r/test"
apt-get -yq update >> ~/apt-get-update.log
apt-get -yq --no-install-suggests --no-install-recommends --force-yes install g++-8 g++-9 clang-8 clang-9 libblas-dev liblapack-dev liblapacke-dev gdb cmake cmake-data intel-tbb-2019.8-075
mkdir -p ${TRAVIS_BUILD_TOPDIR}
cd ${TRAVIS_BUILD_TOPDIR}
git clone git@github.com:ValeevGroup/SeQuant2 ${TRAVIS_BUILD_TOPDIR}/ValeevGroup/SeQuant
cd ${TRAVIS_BUILD_TOPDIR}/ValeevGroup/SeQuant
git checkout ${CURRENT_REV_SHA}
END
chmod +x $setup

##############################################################
# make a script to build all extra prereqs once in the container
build=build.sh
cat > $build << END
#!/bin/sh
cd /home/travis/_build
export BUILD_PREFIX=/home/travis/_build
export INSTALL_PREFIX=/home/travis/_install
export TRAVIS_BUILD_DIR=${TRAVIS_BUILD_TOPDIR}/ValeevGroup/SeQuant
\${TRAVIS_BUILD_DIR}/bin/build-linux.sh
END
chmod +x $build

##############################################################
# make Dockerfile
cat > Dockerfile << END
# Travis default 'Bionic' image
FROM travisci/ci-ubuntu-1804:packer-1577347966-74db3f91

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# create source, build, and install dirs
RUN mkdir -p /home/travis/_build
RUN mkdir -p /home/travis/_install

# Authorize SSH Host
RUN mkdir -p /root/.ssh && \
    chmod 0700 /root/.ssh && \
    ssh-keyscan github.com > /root/.ssh/known_hosts

# Add the keys and set permissions
ADD docker.id_rsa /root/.ssh/id_rsa
ADD docker.id_rsa.pub /root/.ssh/id_rsa.pub
RUN chmod 600 /root/.ssh/id_rsa && \
    chmod 600 /root/.ssh/id_rsa.pub

# install prereqs
ADD $setup /home/travis/_build/$setup
RUN /home/travis/_build/$setup

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# copy travis scripts + ssh key
ADD $build /home/travis/_build/$build
END

function clean_up {
  rm -f $setup $build Dockerfile
  exit
}

trap clean_up SIGHUP SIGINT SIGTERM

##############################################################
# build a dev image
docker build -t sequant-travis-debug .

##############################################################
# extra admin tasks, uncomment as needed

##############################################################
# done
clean_up
