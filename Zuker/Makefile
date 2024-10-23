NAME := zuker

SRCDIR = src
OBJDIR = obj

SRC := main.cpp \
		zuker.cpp
OBJ := $(SRC:.cpp=.o)

SRC := $(addprefix $(SRCDIR)/, $(SRC))
OBJ := $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

CXX := c++
CXXFLAGS := -Wall -Wextra -Werror -std=c++20 -O2 -Iinclude -MMD -MP
CXXFLAGS += -fsanitize=address,undefined,bounds

all: $(NAME)
.PHONY: all

debug: CXXFLAGS += -DDEBUG
debug: re
.PHONY: debug

$(NAME): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(NAME)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(OBJDIR) $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) -r $(OBJDIR)
.PHONY: clean

distclean: clean
	$(RM) $(NAME)
.PHONY: distclean

re: distclean all
.PHONY: re
