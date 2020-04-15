add_subdirectory (async)
add_subdirectory (sync)

if (OPENSSL_FOUND)
    add_subdirectory (async-ssl)
    add_subdirectory (sync-ssl)
endif()
