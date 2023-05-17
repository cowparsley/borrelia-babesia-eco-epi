library(RMark)

# rm(list = ls())

# read in all the data, BI and CT combined
setwd("C:/Users/Ben Adams/OneDrive - University of Bath/borrelia communities/multistrain/babesiaBorrelia/semidiscrete/data/mouseInfectionsMSM")
data = read.csv("./mice_data_2014-2016.csv")

# extract the data for sessions 1 - 7 only, omit session 8 as not consistently sampled
data.CT = data[data$State == "CT" & data$Session %in% 1:7, ]  # data for CT region
data.BI = data[data$State == "RI" & data$Session %in% 1:7, ]  # data for BI region

data = data.CT  # use the CT data
# data = data.BI  # use the BI data
# add a column to the dataframe summarising the state as a single integer
# 1: no infection, 2: Bb only, 3: Bab only, 4: coinfection
data$infState = (data$Bab_micro == 'Y')*2 + (data$Bor_burg == 'Y') + 1

# add a column that combines year and tag number so that all mice are uniquely tagged
# some tags appear in more than one year, but no mice were sampled across years
data$year_tag_num = paste0(data$Year,'.',data$tag_num)

# get the unique tags of each mouse
unique.tag.num = unique(data$year_tag_num)

# sort the data frame so that unique mice are grouped and in session order
ordered.data <- data[order(data$year_tag_num, data$Session),]

# get the number of sessions i.e. the unique number of sessions recorded in the data (could be less than 7 for some subsets)
num.time.steps = length(unique(data$Session))
# get the number of unique mice (some were sampled in multiple sessions)
num.mice = length(unique(data$tag_num))
# set up a matrix of observations times - each row is unique mouse, column is session, initially populate with 0s
obs.time = matrix(0, nrow=num.mice, ncol=num.time.steps)
# set up a matrix of infection states - each row is unique mouse, column is session, initially populate with 0s (state not observed)
Y = matrix(0, nrow=num.mice, ncol=num.time.steps)
# get the unique tags of each mouse
unique.tag.num = unique(data$tag_num)

# each row of the state matrix corresponds to a mouse, each column is its state at the 7 observation times
for(m in 1:num.mice) {
  # for each mouse get all the observations and remove any duplicates
  data.m = ordered.data[ordered.data$tag_num == unique.tag.num[m], ]
  data.m = data.m[!duplicated(data.m$Session),]
  # get the sessions in which it was observed, and the observed states
  sessions = data.m$Session
  states = data.m$infState
  # record the state when observed - default when not observed is 0, uninfected = 1, Bb only = 2, Bab only = 3, Both = 4
  # each row of Y is now the observed state history for an individual mouse
  Y[m,sessions] = states
}

# transform the state matrix into a list of strings that can be imported into RMark
# here apply transforms the elements in each row of Y to a single string
# we store the result in column ch of a dataframe with additional columns frequency (1 for eaverything)
mice.df=as.data.frame( apply(Y, 1, paste, collapse = "") ); names(mice.df) = "ch"; mice.df$freq = 1; 
head(mice.df)

# preliminary analysis: show possible transitions 
find.possible.transitions(mice.df$ch)


# show transition pairs for same data
tp <- transition.pairs(mice.df$ch)
# display, and get the total number of transitions
tp
sum(tp, na.rm = TRUE)

# process the data to get it into the right format for mark
mice.processed = process.data(mice.df, model = "Multistrata") 

# construct the design data
# we specify that transitions that are not observed/rarely observed, as shown in tp above, should be calculated by subtraction 
#mice.ddl = make.design.data(mice.processed, parameters= list(Psi=list(subtract.stratum=c("2","4","2", "1"))))  # use this for BI
mice.ddl = make.design.data(mice.processed, parameters= list(Psi=list(subtract.stratum=c("2","2","2", "2"))))  # use this for CT


# Define models for S (survival probability), formula refers to GLM for this parameter set
S.dot = list(formula=~1)        # survival same for all states, independent of time
# S.stratum=list(formula=~stratum)  # survival different for each state (stratum), independent of time

# Define models for p (observation probability)
p.dot = list(formula=~1)       # observation probability same for all states, independent of time
# p.stratum=list(formula=~stratum)  # observation probability depends on state, independent of time

# Define model for Psi (state transition probabilities)
Psi.dot = list(formula=~-1+stratum:tostratum) # independent of time, using ~ stratum gives same result for Psi, but different betas, see make.mark.model

# Run model with S, p independent of stratum and time, Psi independent of time
model.dot = mark(data=mice.processed, ddl=mice.ddl, model.parameters = list(S = S.dot, p = p.dot, Psi = Psi.dot), output = FALSE)

model.dot.real = model.dot$results$real   # dataframe of the simplified real parameters

summary(model.dot)
model.dot.real


# get the transition probabilities between two observations
Psilist=get.real(model.dot,"Psi",vcv=TRUE)
Psivalues=Psilist$estimates

# get the estimated at a particular time point - here any time is OK as transition rates independent of time
pmat.mark = TransitionMatrix(Psivalues[Psivalues$time==2 & Psivalues$age == 1,],vcv.real=Psilist$vcv.real)


