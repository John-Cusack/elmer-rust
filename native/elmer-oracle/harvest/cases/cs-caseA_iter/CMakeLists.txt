INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(caseA_iter.sif caseA_iter.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO chan scalars.dat scalars.dat.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(cs-caseA_iter LABELS harvest)
