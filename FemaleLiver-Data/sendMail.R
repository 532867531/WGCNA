goMail=function(thebody){
    # MyEmail.R
    library(mailR)
    sender <- "xtlyszj@yeah.net"
    recipients <- c("532867531@qq.com")
    aletter=send.mail(from = sender,
                      to = recipients,
                      subject = paste("Program Done.",Sys.timezone(),Sys.time(),thebody,sep = "_"),
                      body = "My program is finished.",
                      smtp = list(host.name = "smtp.yeah.net", port = 25,
                                  user.name = "xtlyszj@yeah.net",
                                  passwd = "532867531asdfg", ssl = FALSE),
                      authenticate = TRUE,
                      send = FALSE)
    r=aletter$send()
}