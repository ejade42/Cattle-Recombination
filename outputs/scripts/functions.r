## Contains useful functions I use throughout the rest of the scripts, often just for diagnostic printing purposes.


## Python-like f strings
library(glue)   ## v1.6.2 
fstring <- function(string) {as.character(glue(string))}
fprint  <- function(string) {print(as.character(glue(string)), quote=F)}


## conversion between vector and string (again, why can't R simply be like Python)
string_to_vector  <- function(string) {as.numeric(unlist(strsplit(string, split = " ")))}
vector_to_string  <- function(vector) {paste(vector, collapse = " ")}


## print a blank line (for some bizarre reason R doesn't do this automatically when using print() with no arguments)
print_blank = function() {print("", quote = F)}


## time processing functions
get_time_diff_min <- function(time1, time2) {
    round(as.numeric(difftime(time2, time1, units = "min")), 2)
}

mins_to_hours <- function(time_mins) {
    num_hours <- floor(time_mins / 60)
    num_mins  <- time_mins %% 60
    if (num_hours == 1) {hour_message <- "hour"}   else {hour_message <- "hours"}
    if (num_mins  == 1) {mins_message <- "minute"} else {mins_message <- "minutes"}
    message <- paste(num_hours, hour_message, num_mins, mins_message)
    message
}

print_end_message <- function(start_time, done = T) {
    if (done == T) {fprint("Done.")}
    end_time <- Sys.time()
    fprint(paste("End time:", round(end_time)))
    fprint(paste("Total time elapsed:", get_time_diff_min(start_time, end_time), "minutes"))
}

