.PHONY: all scenarios ICU_durations

all: scenarios ICU_durations

scenarios:
	cd scenarios && cargo build --release

ICU_durations:
	cd ICU_durations && cargo build --release
