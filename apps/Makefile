# Enable this if the input graph has weights. Disable this if you want to use synthetically generated weights.
ifdef WEIGHTED
EDGEDATA = -DEDGEDATA
endif

# Enable this if number of vertices > 2^32 
ifdef LONGVERTEXCOUNT
INTV = -DLONG
endif



ifdef WORKFLAG
EDGEWORK = -DEDGEWORK
endif


ifdef IGNOREDEVIATION
IGNORE = -DIGNORE
endif

ifdef ONEAGGREGATIONDEVIATION
ONE = -DONE
endif

ifdef MULTIAGGREGATIONDEVIATION
MULTI = -DMULTI
endif

ifdef PRINTVALUEATTENTIME
TENTIME = -DTENTIME
endif

ifdef COUNTINGEDGE
COUNT = -DCOUNT
endif

ifdef DELTATIMECOMPUTE
DELTATIME = -DELTATIME
endif

ifdef TINYMODE
TINY = -DTINY
endif

ifdef SUBTIMEDEBUG
SUBTIME = -DSUBTIME
endif

INTE = -DEDGELONG

#compilers
PCC = g++
LDFLAGS = -L../lib/mimalloc/out/release -lmimalloc 
# PCFLAGS = -std=c++14 -g -fcilkplus -lcilkrts -O3 -DCILK $(INTV) $(INTE) $(EDGEDATA) $(EDGEWORK)
PCFLAGS = -std=c++14 -g -fcilkplus -lcilkrts -O3 -DCILK $(INTV) $(INTE) $(EDGEDATA) $(EDGEWORK) $(IGNORE) $(ONE) $(MULTI) $(TENTIME) $(COUNT) $(DELTATIME) $(TINY) $(SUBTIME)

# dependencies
COMMON=../core/common/binary_search.h ../core/common/bitsetscheduler.h ../core/common/blockRadixSort.h ../core/common/densebitset.h ../core/common/gettime.h ../core/common/index_map.h ../core/common/matrix.h ../core/common/maybe.h ../core/common/parallel.h ../core/common/parseCommandLine.h ../core/common/transpose.h ../core/common/quickSort.h ../core/common/rwlock.h ../core/common/sequence.h ../core/common/utils.h ../core/common/debugUtils.h

GRAPH= ../core/graph/graph.h ../core/graph/graphUtils.h ../core/graph/IO.h ../core/graph/vertex.h ../core/graph/vertexSubset.h

GRAPHBOLT=../core/graphBolt/AdaptiveExecutor.h ../core/graphBolt/GraphBoltEngine.h ../core/graphBolt/GraphBoltEngine_simple.h ../core/graphBolt/ingestor.h

OTHERS=../core/main.h

ALL=PageRank LabelPropagation COEM 

# make

all: $(ALL)

% : %.C $(COMMON) $(GRAPH) $(GRAPHBOLT) $(OTHERS)
	$(PCC) $(PCFLAGS) $< $(LDFLAGS) -o $@

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)

