#!/bin/sh -e

# workaround for really bare-bones Archlinux containers:
if [ -x "$(command -v pacman)" ]; then
    pacman --noconfirm -Sy
    pacman --noconfirm -S grep gawk sed
fi

distro_id=$(grep '^ID=' /etc/os-release|awk -F = '{print $2}'|sed 's/\"//g')

case "$distro_id" in
    'fedora')
        dnf -y --refresh install cmake make g++ diffutils #wget git
        ;;

    'debian')
        apt-get update
        apt-get install -y cmake make g++ python3 #wget git
        ;;

    'arch')
        pacman --noconfirm -Syu
        pacman --noconfirm -S cmake make gcc diffutils python3 #wget git
        ;;

    'ubuntu')
        apt-get update
        apt-get install -y cmake make g++ python3 doxygen #wget
        ;;

    'centos'|'rhel')
        yum -y install epel-release
        yum -y install cmake3 make gcc-c++ python3 #wget git
        ln -s /bin/cmake3 /bin/cmake
        ;;

    'opensuse'|'opensuse-tumbleweed')
        zypper --non-interactive refresh
        zypper --non-interactive install cmake make gcc-c++ #wget tar git
        ;;

    'alpine')
        apk update
        apk add cmake make g++ bash python3 coreutils #git
        ;;
    *)
        echo "Unsupported distribution: $distro_id"
        exit 1
        ;;
esac

# Download CMake and add its location to the path
#wget https://cmake.org/files/v3.17/cmake-3.17.2-Linux-x86_64.tar.gz
#tar xf cmake-3.17.2-Linux-x86_64.tar.gz
#export PATH="`pwd`/cmake-3.17.2-Linux-x86_64/bin:$PATH"
