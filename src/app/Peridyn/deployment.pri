unix {
    isEmpty(DESTDIR) {
        target.path = /usr/local/bin
        CONFIG(sbsPGP) {
            target.path = /home/sigve/usr/bin/
        }
        export(target.path)
    } else {
        target.path = $$DESTDIR
        export(target.path)
    }
    INSTALLS += target
}

export(INSTALLS)
