\
FC      = gfortran
FFLAGS  = -O2 -Wall -Wextra -fcheck=bounds -J build
LDLIBS  = -llapack -lblas

BUILD   = build
RESULTS = results

# Order matters: each module must be compiled after the modules it
# depends on ("use"s). fpm figures this out automatically; with a
# plain Makefile we spell it out explicitly.
SRC_ORDER = src/kinds_module.f90 \
            src/grid_module.f90 \
            src/io_module.f90 \
            src/poisson_module.f90 \
            src/transport_module.f90 \
            src/problem_definitions_module.f90

OBJS = $(patsubst src/%.f90,$(BUILD)/%.o,$(SRC_ORDER))

.PHONY: all clean test run-poisson run-transport dirs

all: dirs poisson transport

dirs:
	mkdir -p $(BUILD) $(RESULTS)

$(BUILD)/%.o: src/%.f90 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

poisson: dirs $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) app/main_poisson.f90 -o $(BUILD)/poisson $(LDLIBS)

transport: dirs $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) app/main_transport.f90 -o $(BUILD)/transport $(LDLIBS)

test: dirs $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) test/check_poisson.f90 -o $(BUILD)/check_poisson $(LDLIBS)
	./$(BUILD)/check_poisson

run-poisson: poisson
	./$(BUILD)/poisson

run-transport: transport
	./$(BUILD)/transport

clean:
	rm -rf $(BUILD) $(RESULTS)/*.dat
