options(stringsAsFactors=FALSE)

pixel2um  <- 0.293

p2toum    <- pixel2um^2

p2tomm    <- pixel2um^2/(1000^2)

getConstantBoundingBox <- function(cdid, borderPad = 20){
    if(cdid == "mel_19"){
        bbFOV <<- list(X0=1,Y0=-2298,X1=3640,Y1=1) 
    } else {
        bbFOV <<- list(X0=1,Y0=-3375,X1=5363,Y1=1)
    }
    #bb        <<- list(X0 = bbFOV$X0 + borderPad/pixel2um,
    #                   X1 = bbFOV$X1 - borderPad/pixel2um,
    #                   Y0 = bbFOV$Y0 + borderPad/pixel2um,
    #                   Y1 = bbFOV$Y1 - borderPad/pixel2um)
    #bb
}
