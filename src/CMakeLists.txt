add_executable(main main.f90)
target_link_libraries(main
    PRIVATE
        calculation
        integration
        parameters
        initialization
    )

add_fortran_library(calculation
    calculation.f90
    parameters.f90
    )

add_fortran_library(integration
    integration.f90
    calculation.f90
    parameters.f90
    )

add_fortran_library(parameters
    parameters.f90
    )

# add_fortran_library(transformation
#     transformation.f90
#     calculation.f90
#     parameters.f90
#     )

add_fortran_library(initialization
    initialization.f90
    calculation.f90
    parameters.f90
    )

# add_subdirectory(tests)
