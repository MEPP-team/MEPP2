if (OPENSSL_FOUND)
    add_executable (http-client-async-ssl http_client_async_ssl.cpp)
    target_link_libraries (http-client-async-ssl ${Boost_LIBRARIES} ${OPENSSL_LIBRARIES})
endif()
