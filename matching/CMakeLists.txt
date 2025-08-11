set(ENGINE_FILE_LIST
        matchingcommand.cpp matchingcommand.h
        engine.cpp
        execution_tree.h
        execution_tree_node.h
        execution_tree_node.cpp
        execution_tree_generator.cpp
        execution_tree_generator.h
        execution_tree.cpp
        preprocessor.cpp
        preprocessor.h
        encoder.cpp
        encoder.h
        encoderSP.cpp
        encoderSP.h
        query_plan_generator.cpp
        query_plan_generator.h
        crqEnum.cpp
        crqEnum.h
        drqEnum.cpp
        drqEnum.h
        crqSampling.cpp
        crqSampling.h
        drqSampling.cpp
        drqSampling.h
)

add_executable(GPESCO.out ${ENGINE_FILE_LIST})
add_executable(Playground.out playground.cpp)
target_link_libraries(GPESCO.out graph utility nd)