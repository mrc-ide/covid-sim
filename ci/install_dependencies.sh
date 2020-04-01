#!/bin/sh -e

# workaround for really bare-bones Archlinux containers:
if [ -x "$(command -v pacman)" ]; then
    pacman --noconfirm -Sy
    pacman --noconfirm -S grep gawk sed
fi

distro_id=$(grep '^ID=' /etc/os-release|awk -F = '{print $2}'|sed 's/\"//g')

case "$distro_id" in
    'fedora')
        dnf -y --refresh install cmake make g++ diffutils
        ;;

    'debian')
        apt-get update
        apt-get install -y cmake make g++ python3
        ;;

    'arch')
        pacman --noconfirm -Syu
        pacman --noconfirm -S cmake make gcc diffutils python3
        ;;

    'ubuntu')
        apt-get update
        apt-get install -y cmake make g++ python3
        ;;

    'centos'|'rhel')
        yum -y install epel-release
        yum -y install cmake3 make gcc-c++
        ln -s /bin/cmake3 /bin/cmake
        ;;

    'opensuse'|'opensuse-tumbleweed')
        zypper --non-interactive refresh
        zypper --non-interactive install cmake make gcc-c++
        ;;

    'alpine')
        apk update
        apk add cmake make g++ bash python3 coreutils
        ;;
    *)
        echo "Unsupported distribution: $distro_id"
        exit 1
        ;;
esac
