cmake_minimum_required(VERSION 2.8)

project(ROBUSTBEM)

set(DIRECTORIES src examples)
foreach(directory ${DIRECTORIES})
	add_subdirectory(${directory})
endforeach(directory)




