if (OPENSSL_FOUND)
    add_executable (http-client-sync-ssl http_client_sync_ssl.cpp)
    target_link_libraries (http-client-sync-ssl  ${Boost_LIBRARIES} ${OPENSSL_LIBRARIES})
endif()
