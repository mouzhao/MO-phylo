set( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -msse3 -DRAXML_USE_LLALLOC -D_USE_PTHREADS -D_OPTIMIZED_FUNCTIONS -D__SIM_SSE3 -fno-builtin" )

add_executable( raxml_light axml.c  optimizeModel.c trash.c searchAlgo.c topologies.c treeIO.c models.c evaluatePartialGenericSpecial.c evaluateGenericSpecial.c newviewGenericSpecial.c makenewzGenericSpecial.c bipartitionList.c restartHashTable.c fastDNAparsimony.c randomTree.c lockless_allocator/ll_alloc.c mem_alloc.c recom.c)

target_link_libraries( raxml_light m pthread )