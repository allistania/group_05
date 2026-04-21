cmake_minimum_required(VERSION 3.10)
project(Shem)

# Устанавливаем стандарт C++ (рекомендуется 17 из-за structured bindings)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# НАСТРОЙКА ВЫХОДНЫХ ДИРЕКТОРИЙ - чтобы все бинарные файлы были в одной папке
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Создаем папку для бинарных файлов
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Включаем все исходные файлы
add_executable(${PROJECT_NAME}
    main.cpp
    Feoktistov/FeokFunktion.cpp
    Medyakova/Output.cpp
    Batarin/parametrs.cpp
)

# Указываем пути для поиска заголовочных файлов
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)



