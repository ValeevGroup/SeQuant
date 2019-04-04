#!/bin/sh

# this script builds a 'Trusty' env docker image used by Travis-CI for ValeevGroup/SeQuant2 (aka SeQuant) project
#
# to run bash in the image: docker run -it sequant-travis-debug bash -l
# N.B. relevant locations:
#   - source dir: /home/travis/build/ValeevGroup/SeQuant (TRAVIS_BUILD_DIR env in Travis jobs)
#   - build dir: /home/travis/_build
#   - install dir: /home/travis/_install

# this is where in the container file system Travis-CI "starts"
export TRAVIS_BUILD_TOPDIR=/home/travis/build

##############################################################
# make a script to download all prereqs and clone ValeevGroup/SeQuant2 repo
setup=setup.sh
cat > $setup << END
#!/bin/sh
curl -sSL "http://apt.llvm.org/llvm-snapshot.gpg.key" | apt-key add -
echo "deb http://apt.llvm.org/xenial/ llvm-toolchain-xenial-8 main" | tee -a /etc/apt/sources.list > /dev/null
apt-add-repository -y "ppa:ubuntu-toolchain-r/test"
apt-get -yq update >> ~/apt-get-update.log
apt-get -yq --no-install-suggests --no-install-recommends --force-yes install g++-8 clang-8 libc++-8-dev libc++abi-8-dev gdb cmake cmake-data
mkdir -p ${TRAVIS_BUILD_TOPDIR}
cd ${TRAVIS_BUILD_TOPDIR}
git clone git@github.com:ValeevGroup/SeQuant2 ${TRAVIS_BUILD_TOPDIR}/ValeevGroup/SeQuant
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
./build-linux.sh
END
chmod +x $build

##############################################################
# make Dockerfile
cat > Dockerfile << END
# Travis default 'Xenial' image
# for up-to-date info: https://docs.travis-ci.com/user/common-build-problems/#troubleshooting-locally-in-a-docker-image
FROM travisci/ci-sardonyx:packer-1541445940-e193d27

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
ADD build-boost-linux.sh /home/travis/_build/build-boost-linux.sh
ADD build-linux.sh /home/travis/_build/build-linux.sh
ADD build-rangev3-linux.sh /home/travis/_build/build-rangev3-linux.sh
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
