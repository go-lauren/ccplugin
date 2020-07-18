
target_sources( ${PROJECT_NAME}
	PRIVATE
		${CMAKE_CURRENT_LIST_DIR}/ActionA.h
		${CMAKE_CURRENT_LIST_DIR}/qPanels.h
		${CMAKE_CURRENT_LIST_DIR}/fParamsDlg.h
		${CMAKE_CURRENT_LIST_DIR}/facade.h
		${CMAKE_CURRENT_LIST_DIR}/tools.h
)

target_include_directories( ${PROJECT_NAME}
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}
)
