include (ExternalProject)

# Example call:
#
# ```cmake
# declare_datasource (
#   FILE pdb100d.ent.gz # build/data/pdb100d.ent.gz
#   URL ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/00/pdb100d.ent.gz # 16KiloByte
#   URL_HASH SHA256=c2b8f884568b07f58519966e256e2f3aa440508e8013bd10e0ee338e138e62a0)
# ```
#
# Options:
#
# declare_datasource (FILE <datasource name> URL <url1> [<url2>...] [URL_HASH <algo>=<hashValue>] [<option>...])
#
# FILE <datasource name> (required):
#    The name of the downloaded file (must be unique across all declared datasources).
#    This will put the file into the ``<build>/data/<datasource name>`` folder.
#
# URL <url1> [<url2>...] (required):
#    List of paths and/or URL(s) of the external datasource. When more than one URL is given, they are tried in
#    turn until one succeeds. A URL may be an ordinary path in the local file system (in which case it must be the only
#    URL provided) or any downloadable URL supported by the file(DOWNLOAD) command. A local filesystem path may refer to
#    either an existing directory or to an archive file, whereas a URL is expected to point to a file which can be
#    treated as an archive.
#
# URL_HASH <algo>=<hashValue> (optional):
#    Hash of the archive file to be downloaded. The argument should be of the form <algo>=<hashValue> where algo can be
#    any of the hashing algorithms supported by the file() command. Specifying this option is strongly recommended for
#    URL downloads, as it ensures the integrity of the downloaded content. It is also used as a check for a previously
#    downloaded file, allowing connection to the remote location to be avoided altogether if the local directory already
#    has a file from an earlier download that matches the specified hash.
#
# This uses under the hood ExternalProject's and you can pass any viable option of ExternalProject to this function and
# overwrite the default behaviour. See https://cmake.org/cmake/help/latest/module/ExternalProject.html for more
# information.
function (declare_datasource)
    set (options "")
    set (one_value_args FILE URL_HASH)
    set (multi_value_args URL)

    cmake_parse_arguments (ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    string (TOLOWER "datasource--${ARG_FILE}" datasource_name)

    # create data folder
    file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)

    ExternalProject_Add (
        "${datasource_name}"
        URL "${ARG_URL}"
        URL_HASH "${ARG_URL_HASH}"
        DOWNLOAD_NAME "${ARG_FILE}"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND
            ${CMAKE_COMMAND} -E create_symlink <DOWNLOADED_FILE> ${CMAKE_CURRENT_BINARY_DIR}/data/${ARG_FILE}
        TEST_COMMAND ""
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_datasources"
        DOWNLOAD_NO_EXTRACT TRUE # don't extract archive files like .tar.gz.
        ${ARG_UNPARSED_ARGUMENTS}
    )
endfunction ()

# Example call:
#
# ```cmake
# # my_app uses and needs the files RF00001.fa.gz, pdb100d.ent.gz, and GRCh38_latest_clinvar.vcf.gz.
# target_use_datasources (my_app FILES RF00001.fa.gz pdb100d.ent.gz GRCh38_latest_clinvar.vcf.gz)
# ```
#
# Options:
#
# target_use_datasources (<target> FILES <file1> [<file2>...])
#
# It declares that a <target> uses and depends on the following files.
#
# This also sets the build requirement that the files must be downloaded before the <target> will be build.
#
# The named <target> must have been created by a command such as add_executable() or add_library().
function (target_use_datasources target)
    set (options "")
    set (one_value_args)
    set (multi_value_args "FILES")

    cmake_parse_arguments (ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    foreach (filename ${ARG_FILES})
        string (TOLOWER "datasource--${filename}" datasource_name)
        add_dependencies ("${target}" "${datasource_name}")
    endforeach ()
endfunction ()
