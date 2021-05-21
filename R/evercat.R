#General shared R utilities

evercatClass<-setRefClass("evercat",
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
                     cat=function(message,append=F,new=F,end=F)
                     {
                       if(new | append | previous.message.length<1){
                          base::cat(message)
                       } else {
                          base::cat(rep("\b", previous.message.length),collapse="")
                          base::cat(message)
                       }
                       
                       if(end){
                         previous.message.length<<-0
                       } else if(append) {
                         previous.message.length<<-previous.message.length+nchar(message, type ="chars")
                       } else {
                         previous.message.length<<-nchar(message, type ="chars")
                       }
                     }
                   )
)

ever<-evercatClass()
