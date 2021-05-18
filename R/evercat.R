#General shared R utilities

evercat<-setRefClass("evercat",
                   fields = list(
                     previous.message.length = "numeric"
                   ),
                   methods = list
                   (
                     initialize=function()
                     {
                       previous.message.length<<-0
                     },
                     
                     # Persistent program print message, evercat
                     cat=function(message,new=FALSE,end=FALSE)
                     {
                       if(new | previous.message.length<1){
                        base::cat(message)
                       }
                       else {
                         base::cat(paste0(rep("\b", previous.message.length),collapse=""),message)
                       }
                       
                       if(end){
                         previous.message.length<<-0
                       } else {
                         previous.message.length<<-nchar(message)
                       }
                     }
                   )
)

# shruC<-setRefClass("shru",
#                    fields = list(
#                      evercat.previous.message.length = "numeric"
#                    ),
#                    methods = list
#                    (
#                      initialize=function()
#                      {
#                        evercat.previous.message.length<<-0
#                      },
#                      
#                      # Persistent program print message, evercat
#                      evercat=function(message,new=FALSE,end=FALSE)
#                      {
#                        if(new){
#                          evercat.previous.message.length<<-0
#                          cat("\n")
#                          cat(message)
#                        }
#                        else {
#                          cat(
#                            paste(rep("\b", evercat.previous.message.length+1),collapse=""),
#                            message)
#                        }
#                        
#                        if(end){
#                          evercat.previous.message.length<<-0
#                          cat("\n")
#                        } else {
#                          evercat.previous.message.length<<-nchar(message)
#                        }
#                        
#                      }
#                    )
# )
# 
# shru<-shruC()

# shru$evercat.previous.message.length=0
# shru$evercat<-function(message,new=FALSE,end=FALSE){
#   
#   if(end){
#     shru$evercat.previous.message.length<-0
#     cat("\n")
#   } else {
#     
#     if(new){
#       shru$evercat.previous.message.lengt<-0
#       cat("\n")
#       cat(message)
#     }
#     else {
#       cat(
#         paste(rep("\b", shru$evercat.previous.message.length), collapse = ""),
#         message
#       )
#     }
#    
#     shru$evercat.previous.message.length<-nchar(message)
#   }
#   
# }