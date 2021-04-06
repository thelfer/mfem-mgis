if(MFEMMGIS_HAVE_PANDOC)
  function(mfemmgis_pandoc_generate_html_page file)
    ADD_CUSTOM_COMMAND(
      OUTPUT    ${CMAKE_CURRENT_BINARY_DIR}/${file}.html
      DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/${file}.md
      DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/template.html
      DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/css/main.css
      COMMAND   ${PANDOC}
      ARGS      -f markdown-markdown_in_html_blocks+tex_math_single_backslash
      ARGS      --filter pandoc-crossref
      ARGS      --bibliography=${CMAKE_CURRENT_SOURCE_DIR}/bibliography.bib --citeproc
      ARGS      --template=${CMAKE_CURRENT_SOURCE_DIR}/template.html --css=${CMAKE_CURRENT_SOURCE_DIR}/css/main.css
      ARGS      --mathjax
      ARGS      --highlight-style=tango
      ARGS      --email-obfuscation=javascript
      ARGS       --default-image-extension=svg
      ARGS      ${ARGN}
      ARGS      ${CMAKE_CURRENT_SOURCE_DIR}/${file}.md -o ${file}.html)
    add_custom_target(${file}-html ALL DEPENDS ${file}.html)
    add_dependencies(website ${file}-html)
    if(MFEMMGIS_APPEND_SUFFIX)
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${file}.html
        DESTINATION share/doc/mfem-mgis-${MFEMMGIS_SUFFIX}/web
        COMPONENT website)
    else(MFEMMGIS_APPEND_SUFFIX)
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${file}.html
        DESTINATION share/doc/mfem-mgis/web
        COMPONENT website)
    endif(MFEMMGIS_APPEND_SUFFIX)
  endfunction(mfemmgis_pandoc_generate_html_page)
endif(MFEMMGIS_HAVE_PANDOC)

add_custom_target(website)
add_dependencies(doc website)

if(MFEMMGIS_HAVE_PANDOC)
  mfemmgis_pandoc_generate_html_page(index)
else()
endif(MFEMMGIS_HAVE_PANDOC)

if(MFEMMGIS_APPEND_SUFFIX)
  install(DIRECTORY css img js
    DESTINATION share/doc/mfem-mgis-${MFEMMGIS_SUFFIX}/web
    COMPONENT website
    PATTERN ".svn" EXCLUDE)
else(MFEMMGIS_APPEND_SUFFIX)
  install(DIRECTORY css img js
    DESTINATION share/doc/mfem-mgis/web
    COMPONENT website
    PATTERN ".svn" EXCLUDE)
endif(MFEMMGIS_APPEND_SUFFIX)
