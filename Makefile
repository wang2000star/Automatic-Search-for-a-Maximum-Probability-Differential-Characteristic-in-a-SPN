TAR=main.exe
.PHONY: all clean

all:$(TAR)

$(TAR):  fenzhidingjie.o main.o
	gcc -O2 -o $@ $^

%.o : %.c
	gcc -O2 -o $@ -c $<

clean:
	del *.o
	del $(TAR)