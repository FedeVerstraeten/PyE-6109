TARGET = TP.pdf
ENGINE = pdflatex
OUT_FORMAT = latex
IN_FORMAT = markdown
HEADER = header.tex
TEMPLATE = fiuba-template.tex

.PHONY: all clean

all: $(TARGET)

$(TARGET): *.yaml *.md
	pandoc --latex-engine=$(ENGINE) -H $(HEADER) -f $(IN_FORMAT) -t $(OUT_FORMAT) --template $(TEMPLATE) -s $^ -o $(TARGET)

clean:
	-@rm -f $(TARGET)
