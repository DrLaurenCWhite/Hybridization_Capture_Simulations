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
CAPpcr=1 # PCR Rounds after 1st round of Capture, before 2nd round
STEP=100 # Step size
SPEC1=0.1 # Capture specificty first round (prob of capturing off-target read)
SPEC2=0.1 # Capture specificty second round (prob of capturing off-target read)
LengthOut=20000
MAX=LengthOut #Max sequencing effort

#100% Capture efficiency
EFF1=1 # Capture efficiency first round (prob of capturing on-target read)
EFF2=1 # Capture efficiency second round (prob of capturing on-target read)
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
  
  #Shotgun
  Library=rep(1:FRAGS, 2^Lpcr) 
  Seq=c(NULL)
  N=0
  while (N<MAX & length(Library)>STEP) {
    N=N+STEP
    Sub=sample(Library, STEP) #Sequence the library by sampling fragments
    Library=vsetdiff(Library, Sub) #remove the fragments already sequenced
    Seq=c(Seq, Sub) #Add sequenced fragments to the sequenced list
    Uni=length(unique(Seq)) #Count how many sequenced fragments are unique
    TargetUni=length(unique(Seq[Seq<=TARGET])) #Count how many unqiue reads match the target
    New=data.table(n=N, unique=Uni, uniqueOnTarget=TargetUni, Type="Shot", Run=x) #Tabulate results
    Nadd=nrow(New)
    res[Y:(Y+Nadd-1)] <<- New
    Y <<- Y+Nadd
  }

  #One round of capture
  Library=rep(1:FRAGS, 2^Lpcr) 
  Hyb <- rep(NA_real_, length(Library))
  capON <- sample(c(T,F), length(Library),replace=T,prob=c(EFF1,1-EFF1))
  capOFF <- sample(c(T,F), length(Library),replace=T,prob=c(SPEC1,1-SPEC1))
  Hyb=na.exclude(ifelse( (Library<TARGET & capON==T) | (Library>=TARGET & capOFF==T) , Library , Hyb))
  Cycles=ifelse(length(Hyb)<=MAX, ceiling(MAX/length(Hyb)), 1) 
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

  #Two rounds of capture
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
  Cycles=ifelse(length(Hyb2)<=MAX, ceiling(MAX/length(Hyb2)), 1)
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

#50% Capture efficiency
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
  
  #Shotgun
  Library=rep(1:FRAGS, 2^Lpcr) 
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
  
  #One round of capture
  Library=rep(1:FRAGS, 2^Lpcr) #minimal PCRs before going on to capture
  Hyb <- rep(NA_real_, length(Library))
  capON <- sample(c(T,F), length(Library),replace=T,prob=c(EFF1,1-EFF1))
  capOFF <- sample(c(T,F), length(Library),replace=T,prob=c(SPEC1,1-SPEC1))
  Hyb=na.exclude(ifelse( (Library<TARGET & capON==T) | (Library>=TARGET & capOFF==T) , Library , Hyb))
  Cycles=ifelse(length(Hyb)<=MAX, ceiling(MAX/length(Hyb)), 1) 
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
  
  #Two rounds of capture
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
  Cycles=ifelse(length(Hyb2)<=MAX, ceiling(MAX/length(Hyb2)), 1)
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


#25% Capture efficiency
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
  
  #Shotgun
  Library=rep(1:FRAGS, 2^Lpcr) 
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
  
  #One round of capture
  Library=rep(1:FRAGS, 2^Lpcr) #minimal PCRs before going on to capture
  Hyb <- rep(NA_real_, length(Library))
  capON <- sample(c(T,F), length(Library),replace=T,prob=c(EFF1,1-EFF1))
  capOFF <- sample(c(T,F), length(Library),replace=T,prob=c(SPEC1,1-SPEC1))
  Hyb=na.exclude(ifelse( (Library<TARGET & capON==T) | (Library>=TARGET & capOFF==T) , Library , Hyb))
  Cycles=ifelse(length(Hyb)<=MAX, ceiling(MAX/length(Hyb)), 1)
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
  
  #Two rounds of capture
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
  Cycles=ifelse(length(Hyb2)<=MAX, ceiling(MAX/length(Hyb2)), 1)
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
  annotate("rect", xmin = 100, xmax = 10000, ymin = 0, ymax = TARGET, fill = "yellow", alpha = 0.2) +
#  annotate("rect", xmin = 9100, xmax = 10000, ymin = 0, ymax = TARGET, fill = "mediumturquoise", alpha = 0.2) +
  #annotate("rect", xmin = 6600, xmax = 10000, ymin = 0, ymax = TARGET, fill = "red", alpha = 0.2) +
  annotate("label", label="Capture Efficiency = 100%", x=17000, y=50, label.size=NA, fill="white", size=5) +
  geom_line(alpha=0.5)  + ylim(c(0, TARGET)) +
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
  annotate("rect", xmin = 100, xmax = 3000, ymin = 0, ymax = TARGET, fill = "yellow", alpha = 0.2) +
  annotate("rect", xmin = 3100, xmax = 10000, ymin = 0, ymax = TARGET, fill = "mediumturquoise", alpha = 0.2) +
  annotate("rect", xmin = 10100, xmax = 20000, ymin = 0, ymax = TARGET, fill = "purple", alpha = 0.2) +
  annotate("label", label="Capture Efficiency = 50%", x=17000, y=50, label.size=NA, fill="white", size=5) +
  geom_line(alpha=0.5)  + ylim(c(0, TARGET)) +
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
  annotate("rect", xmin = 100, xmax =900, ymin = 0, ymax = TARGET, fill = "yellow", alpha = 0.2) +
  annotate("rect", xmin = 1000, xmax = 4000, ymin = 0, ymax = TARGET, fill = "mediumturquoise", alpha = 0.2) +
  annotate("rect", xmin = 4100, xmax = 20000, ymin = 0, ymax = TARGET, fill = "purple", alpha = 0.2) +
  annotate("label", label="Capture Efficiency = 25%", x=17000, y=50, label.size=NA, fill="white", size=5) +
  geom_line(alpha=0.5) + ylim(c(0, TARGET)) +
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

