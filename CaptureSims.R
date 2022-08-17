rm(list = ls())
library(data.table)
library(ggplot2)
library(vecsets)
library(plyr)
library(ggpubr)


FRAGS=10000 # Number of unique fragments in library
End=0.1 # % endogenous
TARGET=round(FRAGS*End, 0)
Lpcr=1  # PCR Rounds after Library prep before Capture
CAPpcr=2 # PCR Rounds after 1st round of Capture, before 2nd round
STEP=100 # Step size
MAX=10000 # Max seq effort
SPEC1=0.1 # Capture specificty first round (prob of capturing off-target read)
SPEC2=0.1 # Capture specificty second round (prob of capturing off-target read)

EFF1=1 # Capture efficiency first round (prob of capturing on-target read)
EFF2=1 # Capture efficiency second round (prob of capturing on-target read)
LengthOut=20000
res=data.table(n=numeric(LengthOut), unique=numeric(LengthOut), uniqueOnTarget=numeric(LengthOut), Type=factor(LengthOut), Run=numeric(LengthOut))
Y=1

ldply(1:50, function(x) {
  print(x)
  MAX <<- MAX
  STEP <<- STEP
  CAPpcr <<- CAPpcr
  Lpcr <<- Lpcr
  EFF1 <<-EFF1
  EFF2 <<- EFF2
  SPEC1 <<- SPEC1
  SPEC2 <<- SPEC2
  TARGET <<- TARGET
  FRAGS <<- FRAGS
  
  #Cycles=MAX/FRAGS
  Library=rep(1:FRAGS, Lpcr) 
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Library)>STEP) {
    N=N+STEP
    Sub=sample(Library, STEP)
    Library=vsetdiff(Library, Sub)
    Seq=c(Seq, Sub)
    Uni=length(unique(Seq))
    TargetUni=length(unique(Seq[Seq<=TARGET]))
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Shot", Run=x)
    Nadd=nrow(New)
    res[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }

  Library=rep(1:FRAGS, 2^Lpcr) #minimal PCRs before going on to capture
  Hyb <- rep(NA_real_, length(Library))
  capON <- sample(c(T,F), length(Library),replace=T,prob=c(EFF1,1-EFF1))
  capOFF <- sample(c(T,F), length(Library),replace=T,prob=c(SPEC1,1-SPEC1))
  Hyb=na.exclude(ifelse( (Library<TARGET & capON==T) | (Library>=TARGET & capOFF==T) , Library , Hyb))
  Cycles=MAX/length(Hyb) 
  Capture=rep(Hyb, Cycles)
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Capture)>STEP) {
    N=N+STEP
    Sub=sample(Capture, STEP)
    Capture=vsetdiff(Capture, Sub)
    Seq=c(Seq, Sub)
    Uni=length(unique(Seq))
    TargetUni=length(unique(Seq[Seq<=TARGET]))
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Capture1", Run=x)
    Nadd=nrow(New)
    res[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }

  Library=rep(1:FRAGS, 2^Lpcr)
  Hyb <- rep(NA_real_, length(Library))
  capON <- sample(c(T,F), length(Library),replace=T,prob=c(EFF1,1-EFF1))
  capOFF <- sample(c(T,F), length(Library),replace=T,prob=c(SPEC1,1-SPEC1))
  Hyb=na.exclude(ifelse( (Library<TARGET & capON==T) | (Library>=TARGET & capOFF==T) , Library , Hyb))
  Capture=rep(Hyb, 2^CAPpcr)
  Hyb2 <- rep(NA_real_, length(Capture))
  capON <- sample(c(T,F), length(Capture),replace=T,prob=c(EFF2,1-EFF2))
  capOFF <- sample(c(T,F), length(Capture),replace=T,prob=c(SPEC2,1-SPEC2))
  Hyb2=na.exclude(ifelse( (Capture<TARGET & capON==T) | (Capture>=TARGET & capOFF==T) , Capture , Hyb2))
  Cycles=MAX/length(Hyb2) 
  Capture2=rep(Hyb2, Cycles)
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Capture2)>STEP) {
    N=N+STEP
    Sub=sample(Capture2, STEP)
    Capture2=vsetdiff(Capture2, Sub)
    Seq=c(Seq, Sub)
    Uni=length(unique(Seq))
    TargetUni=length(unique(Seq[Seq<=TARGET]))
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Capture2", Run=x)
    Nadd=nrow(New)
    res[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }
})
res=res[Type!="20000"]

EFF1=0.5 # Capture efficiency first round (prob of capturing on-target read)
EFF2=0.5 # Capture efficiency second round (prob of capturing on-target read)
res2=data.table(n=numeric(LengthOut), unique=numeric(LengthOut), uniqueOnTarget=numeric(LengthOut), Type=factor(LengthOut), Run=numeric(LengthOut))
Y=1
ldply(1:50, function(x) {
  print(x)
  MAX <<- MAX
  STEP <<- STEP
  CAPpcr <<- CAPpcr
  Lpcr <<- Lpcr
  EFF1 <<-EFF1
  EFF2 <<- EFF2
  SPEC1 <<- SPEC1
  SPEC2 <<- SPEC2
  TARGET <<- TARGET
  FRAGS <<- FRAGS
  
  #Cycles=MAX/FRAGS
  Library=rep(1:FRAGS, Lpcr) 
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Library)>STEP) {
    N=N+STEP
    Sub=sample(Library, STEP)
    Library=vsetdiff(Library, Sub)
    Seq=c(Seq, Sub)
    Uni=length(unique(Seq))
    TargetUni=length(unique(Seq[Seq<=TARGET]))
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Shot", Run=x)
    Nadd=nrow(New)
    res2[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }
  
  Library=rep(1:FRAGS, 2^Lpcr) #minimal PCRs before going on to capture
  Hyb <- rep(NA_real_, length(Library))
  capON <- sample(c(T,F), length(Library),replace=T,prob=c(EFF1,1-EFF1))
  capOFF <- sample(c(T,F), length(Library),replace=T,prob=c(SPEC1,1-SPEC1))
  Hyb=na.exclude(ifelse( (Library<TARGET & capON==T) | (Library>=TARGET & capOFF==T) , Library , Hyb))
  Cycles=MAX/length(Hyb) 
  Capture=rep(Hyb, Cycles)
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Capture)>STEP) {
    N=N+STEP
    Sub=sample(Capture, STEP)
    Capture=vsetdiff(Capture, Sub)
    Seq=c(Seq, Sub)
    Uni=length(unique(Seq))
    TargetUni=length(unique(Seq[Seq<=TARGET]))
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Capture1", Run=x)
    Nadd=nrow(New)
    res2[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }
  
  Library=rep(1:FRAGS, 2^Lpcr)
  Hyb <- rep(NA_real_, length(Library))
  capON <- sample(c(T,F), length(Library),replace=T,prob=c(EFF1,1-EFF1))
  capOFF <- sample(c(T,F), length(Library),replace=T,prob=c(SPEC1,1-SPEC1))
  Hyb=na.exclude(ifelse( (Library<TARGET & capON==T) | (Library>=TARGET & capOFF==T) , Library , Hyb))
  Capture=rep(Hyb, 2^CAPpcr)
  Hyb2 <- rep(NA_real_, length(Capture))
  capON <- sample(c(T,F), length(Capture),replace=T,prob=c(EFF2,1-EFF2))
  capOFF <- sample(c(T,F), length(Capture),replace=T,prob=c(SPEC2,1-SPEC2))
  Hyb2=na.exclude(ifelse( (Capture<TARGET & capON==T) | (Capture>=TARGET & capOFF==T) , Capture , Hyb2))
  Cycles=MAX/length(Hyb2) 
  Capture2=rep(Hyb2, Cycles)
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Capture2)>STEP) {
    N=N+STEP
    Sub=sample(Capture2, STEP)
    Capture2=vsetdiff(Capture2, Sub)
    Seq=c(Seq, Sub)
    Uni=length(unique(Seq))
    TargetUni=length(unique(Seq[Seq<=TARGET]))
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Capture2", Run=x)
    Nadd=nrow(New)
    res2[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }
})
res2=res2[Type!="20000"]

EFF1=0.25 # Capture efficiency first round (prob of capturing on-target read)
EFF2=0.25 # Capture efficiency second round (prob of capturing on-target read)
res3=data.table(n=numeric(LengthOut), unique=numeric(LengthOut), uniqueOnTarget=numeric(LengthOut), Type=factor(LengthOut), Run=numeric(LengthOut))
Y=1
ldply(1:50, function(x) {
  print(x)
  MAX <<- MAX
  STEP <<- STEP
  CAPpcr <<- CAPpcr
  Lpcr <<- Lpcr
  EFF1 <<-EFF1
  EFF2 <<- EFF2
  SPEC1 <<- SPEC1
  SPEC2 <<- SPEC2
  TARGET <<- TARGET
  FRAGS <<- FRAGS
  
  #Cycles=MAX/FRAGS
  Library=rep(1:FRAGS, Lpcr) 
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Library)>STEP) {
    N=N+STEP
    Sub=sample(Library, STEP)
    Library=vsetdiff(Library, Sub)
    Seq=c(Seq, Sub)
    Uni=length(unique(Seq))
    TargetUni=length(unique(Seq[Seq<=TARGET]))
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Shot", Run=x)
    Nadd=nrow(New)
    res3[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }
  
  Library=rep(1:FRAGS, 2^Lpcr) #minimal PCRs before going on to capture
  Hyb <- rep(NA_real_, length(Library))
  capON <- sample(c(T,F), length(Library),replace=T,prob=c(EFF1,1-EFF1))
  capOFF <- sample(c(T,F), length(Library),replace=T,prob=c(SPEC1,1-SPEC1))
  Hyb=na.exclude(ifelse( (Library<TARGET & capON==T) | (Library>=TARGET & capOFF==T) , Library , Hyb))
  Cycles=MAX/length(Hyb) 
  Capture=rep(Hyb, Cycles)
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Capture)>STEP) {
    N=N+STEP
    Sub=sample(Capture, STEP)
    Capture=vsetdiff(Capture, Sub)
    Seq=c(Seq, Sub)
    Uni=length(unique(Seq))
    TargetUni=length(unique(Seq[Seq<=TARGET]))
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Capture1", Run=x)
    Nadd=nrow(New)
    res3[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }
  
  Library=rep(1:FRAGS, 2^Lpcr)
  Hyb <- rep(NA_real_, length(Library))
  capON <- sample(c(T,F), length(Library),replace=T,prob=c(EFF1,1-EFF1))
  capOFF <- sample(c(T,F), length(Library),replace=T,prob=c(SPEC1,1-SPEC1))
  Hyb=na.exclude(ifelse( (Library<TARGET & capON==T) | (Library>=TARGET & capOFF==T) , Library , Hyb))
  Capture=rep(Hyb, 2^CAPpcr)
  Hyb2 <- rep(NA_real_, length(Capture))
  capON <- sample(c(T,F), length(Capture),replace=T,prob=c(EFF2,1-EFF2))
  capOFF <- sample(c(T,F), length(Capture),replace=T,prob=c(SPEC2,1-SPEC2))
  Hyb2=na.exclude(ifelse( (Capture<TARGET & capON==T) | (Capture>=TARGET & capOFF==T) , Capture , Hyb2))
  Cycles=MAX/length(Hyb2) 
  Capture2=rep(Hyb2, Cycles)
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Capture2)>STEP) {
    N=N+STEP
    Sub=sample(Capture2, STEP)
    Capture2=vsetdiff(Capture2, Sub)
    Seq=c(Seq, Sub)
    Uni=length(unique(Seq))
    TargetUni=length(unique(Seq[Seq<=TARGET]))
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Capture2", Run=x)
    Nadd=nrow(New)
    res3[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }
})
res3=res3[Type!="20000"]



p1=ggplot(res, aes(x=n, y=uniqueOnTarget, col=as.factor(Type), group=interaction(Run, Type))) +
  xlab("Number of Reads Sequenced") + ylab("Unique Reads on Target") +
  annotate("rect", xmin = 100, xmax = 5000, ymin = 0, ymax = 1000, fill = "yellow", alpha = 0.2) +
  #annotate("rect", xmin = 1600, xmax = 6500, ymin = 0, ymax = 1000, fill = "green", alpha = 0.2) +
  #annotate("rect", xmin = 6600, xmax = 10000, ymin = 0, ymax = 1000, fill = "red", alpha = 0.2) +
  annotate("label", label="Capture Efficiency = 100%", x=7500, y=50, label.size=NA, fill="white", size=5) +
  geom_line(alpha=0.5)  +
  guides(color = guide_legend(override.aes = list(lwd = 3, alpha=1) ) ) +
  scale_color_viridis_d(labels=c("Shotgun","Capture\n1 round","Capture\n2 rounds")) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=15),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_blank())
p2=ggplot(res2, aes(x=n, y=uniqueOnTarget, col=as.factor(Type), group=interaction(Run, Type))) +
  xlab("Number of Reads Sequenced") + ylab("Unique Reads on Target") +
  annotate("rect", xmin = 100, xmax = 4500, ymin = 0, ymax = 1000, fill = "yellow", alpha = 0.2) +
  annotate("rect", xmin = 4600, xmax = 7400, ymin = 0, ymax = 1000, fill = "mediumturquoise", alpha = 0.2) +
  annotate("rect", xmin = 7500, xmax = 10000, ymin = 0, ymax = 1000, fill = "purple", alpha = 0.2) +
  annotate("label", label="Capture Efficiency = 50%", x=7500, y=50, label.size=NA, fill="white", size=5) +
  geom_line(alpha=0.5)  +
  guides(color = guide_legend(override.aes = list(lwd = 3, alpha=1) ) ) +
  scale_color_viridis_d(labels=c("Shotgun","Capture\n1 round","Capture\n2 rounds")) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=15),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_blank())
p3=ggplot(res3, aes(x=n, y=uniqueOnTarget, col=as.factor(Type), group=interaction(Run, Type))) +
  xlab("Number of Reads Sequenced") + ylab("Unique Reads on Target") +
  annotate("rect", xmin = 100, xmax = 1500, ymin = 0, ymax = 1000, fill = "yellow", alpha = 0.2) +
  annotate("rect", xmin = 1600, xmax = 4000, ymin = 0, ymax = 1000, fill = "mediumturquoise", alpha = 0.2) +
  annotate("rect", xmin = 4100, xmax = 10000, ymin = 0, ymax = 1000, fill = "purple", alpha = 0.2) +
  annotate("label", label="Capture Efficiency = 25%", x=7500, y=50, label.size=NA, fill="white", size=5) +
  geom_line(alpha=0.5) +
  guides(color = guide_legend(override.aes = list(lwd = 3, alpha=1) ) ) +
  scale_color_viridis_d(labels=c("Shotgun","Capture\n1 round","Capture\n2 rounds")) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=15),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_blank())

fig=ggarrange(p3 + rremove("ylab") + rremove("xlab"), p2 + rremove("ylab") + rremove("xlab"), p1+ rremove("ylab") + rremove("xlab"),
          nrow=3, common.legend = TRUE, legend="top")

annotate_figure(fig, left = text_grob("Unique reads on-target", rot=90),
                bottom = text_grob("Number of reads sequenced"))
