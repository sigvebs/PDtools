#eval(headers.path    = $$OUT_PWD)
#eval(headers.files   += $$HEADERS)
#eval(INSTALLS       += headers)

INSTALL_PREFIX = $$OUT_PWD

for(header, HEADERS) {
    path = $${INSTALL_PREFIX}/$${dirname(header)}
    eval(headers_$${path}.files += $$header)
    eval(headers_$${path}.path = $$path)
    eval(INSTALLS *= headers_$${path})
}
