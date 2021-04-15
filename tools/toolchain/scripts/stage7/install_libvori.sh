#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

# shellcheck disable=SC1090

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libvori_ver="210412"
libvori_sha256="331886aea9d093d8c44b95a07fab13d47f101b1f94a0640d7d670eb722bf90ac"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libvori" ] && rm "${BUILDDIR}/setup_libvori"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_libvori:=__INSTALL__}" in
  __INSTALL__)
    echo "==================== Installing libvori ===================="
    pkg_install_dir="${INSTALLDIR}/libvori-${libvori_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libvori-${libvori_ver} is already installed, skipping it."
    else
      if [ -f libvori-${libvori_ver}.tar.gz ]; then
        echo "libvori-${libvori_ver}.tar.gz is found"
      else
        # shellcheck disable=SC2086
        download_pkg ${DOWNLOADER_FLAGS} ${libvori_sha256} \
          "https://www.cp2k.org/static/downloads/libvori-${libvori_ver}.tar.gz"
      fi

      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libvori-${libvori_ver} ] && rm -rf libvori-${libvori_ver}
      tar -xzf libvori-${libvori_ver}.tar.gz

      mkdir "libvori-${libvori_ver}/build"
      cd "libvori-${libvori_ver}/build"
      cmake \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE \
        .. > cmake.log 2>&1
      CMAKE_BUILD_PARALLEL_LEVEL="$(get_nprocs)" cmake --build . > build.log 2>&1
      CMAKE_BUILD_PARALLEL_LEVEL="$(get_nprocs)" cmake --build . --target test > test.log 2>&1
      CMAKE_BUILD_PARALLEL_LEVEL="$(get_nprocs)" cmake --build . --target install > install.log 2>&1

      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename "${SCRIPT_NAME}")"
    fi
    LIBVORI_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding libvori from system paths ===================="
    check_lib -lvori "libvori"
    add_lib_from_paths LIBVORI_LDFLAGS "libvori.*" "$LIB_PATHS"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking libvori to user paths ===================="
    pkg_install_dir="${with_libvori}"

    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    LIBVORI_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && LIBVORI_LIBDIR="${pkg_install_dir}/lib64"

    check_dir "${LIBVORI_LIBDIR}"
    LIBVORI_LDFLAGS="-L'${LIBVORI_LIBDIR}' -Wl,-rpath='${LIBVORI_LIBDIR}'"
    ;;
esac

if [ "$with_libvori" != "__DONTUSE__" ]; then
  LIBVORI_LIBS="-lvori -lstdc++"
  if [ "$with_libvori" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_libvori"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
export LIBVORI_LIBS="${LIBVORI_LIBS}"
export LIBVORI_ROOT="${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_libvori"
export LIBVORI_ROOT="${pkg_install_dir}"
export LIBVORI_VERSION=${libvori_ver}
export LIBVORI_LDFLAGS="${LIBVORI_LDFLAGS}"
export LIBVORI_LIBRARY="-lvori"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBVORI"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBVORI_LDFLAGS}"
export CP_LIBS="\${CP_LIBS} ${LIBVORI_LIBS}"
EOF
  cat "${BUILDDIR}/setup_libvori" >> "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libvori"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libvori"
