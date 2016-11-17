# This ensures that a certain set of CMake features is available to
# your build.

cmake_minimum_required(VERSION 3.4.1)

# set the native library language
enable_language(C)


# Specifies a library name  specifies whether the library is STATIC or
# SHARED  and provides relative paths to the source code. You can
# define multiple libraries by adding multiple add.library() commands
# and CMake builds them for you. When you build your app  Gradle
# automatically packages shared libraries with your APK.



add_library(native-schur-lib SHARED FillCommandTab.c  FillCommandTab.h FrontPage.c  FrontPage.h  Getl.c  Getl.h  Makefile.am  Makefile.in  ReadWrite.c  ReadWrite.h Scanck.c  Scanck.h  ansiCodes.txt  apropos.c  apropos.h  bignums.c  bignums.h  branch.c  branch.h define.h  dim.h  dimensions.c  dimensions.h  g.c  g.h  gr.c  gr.h  hivesLRcoef.c  hivesLRcoef.h init.c  init.h  label.c  label.h  loadcomList.c  loadcomList.h  m.c  m.h  mymalloc.c  mymalloc.h outerskew.c  outerskew.h  outerskew_inc.c  q1.c  q1.h  q2.c  q2.h  r.c  r.h  rib_to_s.c rib_to_s.h  s.c  s.h  s0.c  s0.h  s1.c  s1.h  s2.c  s2.h  s3.c  s3.h  s4.c  s4.h  s5.c  s5.h  s6.c s6.h  s7.c  s7.h  s8.c  s8.h  savecomList.c  savecomList.h  schur.c  sets_mgmt.c  sets_mgmt.h skew.c  skew.h  standard.h  tableaux.c  tableaux.h  type.h  utils.c  utils.h  var.h  write.c write.h) 
# Specifies a path to native header files.
include_directories(src/main/cpp/schurenativelib/)

