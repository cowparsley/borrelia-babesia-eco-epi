library(msm)

# rm(list = ls())

# read in all the data, BI and CT combined
setwd("C:/Users/Ben Adams/OneDrive - University of Bath/borrelia communities/multistrain/babesiaBorrelia/semidiscrete/data/mouseInfectionsMSM")
data = read.csv("./mice_data_2014-2016.csv")
#year = strftime(as.Date(data$date, "%m/%d/%Y"), format = "%Y")
#data = data[year == "2015" & data$Session %in% 1:7, ]

# extract the CT data for sessions 1 - 7 only, omit session 8 as not consistently sampled
# data = data[data$State == "CT" & data$Session %in% 1:7, ]
# extract the BI data for sessions 1 - 7 only, omit session 8 as not consistently sampled
# data = data[data$State == "RI" & data$Session %in% 1:7, ]
# extract the combined CT and Bi data
data = data[data$State %in% c("CT","RI") & data$Session %in% 1:7, ]

# add a column to the dataframe summarising the state as a single integer
data$infState = (data$Bab_micro == 'Y')*2 + (data$Bor_burg == 'Y') + 1
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
# set up a matrix of infection states - each row is unique mouse, column is session, initially populate with -1s (state not observed)
Y = matrix(-1, nrow=num.mice, ncol=num.time.steps)
# get the unique tags of each mouse
unique.tag.num = unique(data$tag_num)

for(m in 1:num.mice) {
    # for each mouse get all the observations and remove any duplicates
    data.m = ordered.data[ordered.data$tag_num == unique.tag.num[m], ]
    data.m = data.m[!duplicated(data.m$Session),]
	  # get the sessions in which it was observed, and the observed states
    sessions = data.m$Session
    states = data.m$infState
	  # record the observation times 
    obs.time[m, 1:length(sessions)] = sessions
    # record the state when observed - default when not observed is -1, uninfected = 1, Bb only = 2, Bab only = 3, Both = 4
    # each row of Y is now the observed state history for an individual mouse
    Y[m,sessions] = states
}

# construct frequency table of consecutive states
statetable.msm(infState, year_tag_num, data = ordered.data)

# set up state transition matrix, including a fifth state: Dead
# diagonal entries are dummies as calculated internally, other 0 entries are fixed at 0 throughout
# other values are initial estimates
Q <- rbind( c(0, 1, 1, 1),
            c(1, 0, 1, 1),
            c(1, 1, 0, 1),
            c(1, 1, 1, 0))
  
# get a rough estimate of the state transition rates  
Q.crude <- crudeinits.msm(infState ~ Session, year_tag_num, data=ordered.data, qmatrix=Q)

# compute max likelihood estimates of transition intensities
infState.msm <- msm( infState ~ Session, subject = year_tag_num, data = ordered.data, qmatrix = Q.crude)

# extract the transition probabilities in the time interval between two sessions
pmat.msm <- pmatrix.msm(infState.msm, t=1, ci = "norm")

