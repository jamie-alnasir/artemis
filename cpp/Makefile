# Override CXX for native or cross-compilation; no CPU-specific flags by default.
CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic
LDFLAGS ?=
LDLIBS ?=
EXEEXT ?=
TARGET = artemis$(EXEEXT)
TEST_TARGET = test_artemis$(EXEEXT)

.PHONY: all test clean
all: $(TARGET)

$(TARGET): artemis.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) artemis.cpp $(LDFLAGS) $(LDLIBS) -o $@

$(TEST_TARGET): tests/test_artemis.cpp artemis.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests/test_artemis.cpp $(LDFLAGS) $(LDLIBS) -o $@

test: $(TEST_TARGET)
	./$(TEST_TARGET)

clean:
	$(RM) $(TARGET) $(TEST_TARGET)
