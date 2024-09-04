# ##############################################################################
# CMakeLists.txt
#
# Copyright (C) 2024 Florian Kurpicz <florian@kurpicz.org>
#
# pasta::bit_vector is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# pasta::bit_vector is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# pasta::bit_vector.  If not, see <http://www.gnu.org/licenses/>.
#
# ##############################################################################


# find Doxygen
find_package(Doxygen)

if (DOXYGEN_FOUND)
  if (DOXYGEN_VERSION VERSION_LESS "1.10.0")
    message(
      WARNING
        "Doxygen must be version 1.10.0 or newer. Documentation may not be displayed correctly and CI may "
        "not pass even if checks pass locally."
    )
  endif ()
  FetchContent_Declare(
    doxygen-awesome-css
    GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css.git
    GIT_TAG v2.3.3
  )
  FetchContent_MakeAvailable(doxygen-awesome-css)

  # required variables for setting up doxygen-awesome-css
  set(DOXYGEN_GENERATE_TREEVIEW YES)
  set(DOXYGEN_DISABLE_INDEX NO)
  set(DOXYGEN_FULL_SIDEBAR NO)
  set(DOXYGEN_HTML_EXTRA_STYLESHEET
      ${doxygen-awesome-css_SOURCE_DIR}/doxygen-awesome.css
      ${doxygen-awesome-css_SOURCE_DIR}/doxygen-awesome-sidebar-only.css
      ${doxygen-awesome-css_SOURCE_DIR}/doxygen-awesome-sidebar-only-darkmode-toggle.css
      ${CMAKE_CURRENT_LIST_DIR}/custom_doxygen/custom.css
  )
  set(DOXYGEN_HTML_EXTRA_FILES
      ${doxygen-awesome-css_SOURCE_DIR}/doxygen-awesome-darkmode-toggle.js
      ${doxygen-awesome-css_SOURCE_DIR}/doxygen-awesome-fragment-copy-button.js
      ${doxygen-awesome-css_SOURCE_DIR}/doxygen-awesome-paragraph-link.js
  )
  set(DOXYGEN_HTML_HEADER ${CMAKE_CURRENT_LIST_DIR}/custom_doxygen/header.html)
  set(DOXYGEN_HTML_COLORSTYLE LIGHT)

  set(DOXYGEN_PROJECT_LOGO "${CMAKE_CURRENT_LIST_DIR}/images/logo_pasta_bit_vector.svg")
  set(DOXYGEN_IMAGE_PATH "${CMAKE_CURRENT_LIST_DIR}/images")
  set(DOXYGEN_GENERATE_TODOLIST NO)
  set(DOXYGEN_EXTRACT_STATIC YES)
  set(DOXYGEN_CASE_SENSE_NAMES YES)
  set(DOXYGEN_BUILTIN_STL_SUPPORT YES)
  set(DOXYGEN_RECURSIVE YES)
  set(DOXYGEN_STRIP_CODE_COMMENTS NO)
  set(DOXYGEN_GENERATE_LATEX YES)
  set(DOXYGEN_TAB_SIZE 2)
  set(DOXYGEN_FILE_PATTERNS
      *.cpp
      *.hpp
      *.dox
      *.md
  )
  set(DOXYGEN_WARN_AS_ERROR FAIL_ON_WARNINGS)
  set(DOXYGEN_CITE_BIB_FILES ${CMAKE_CURRENT_LIST_DIR}/literature.bib)

  doxygen_add_docs(
    docs ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/../include/
    COMMENT "Generating Documentation"
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
  )
else ()
  message(STATUS "Doxygen not found, not building docs")
endif ()
