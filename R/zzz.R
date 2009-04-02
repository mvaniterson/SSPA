
.onAttach <- function(lib, pkg) {

  if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI == "Rgui"){
  	addVigs2WinMenu("SSPA")
  }

}

