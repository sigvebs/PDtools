unix {
#    isEmpty(target.path) {
##        target.path = /usr/local/bin
#        target.path = /home/sigve/Documents/SLETT/bin
#        export(target.path)
#    }

    isEmpty(DESTDIR) {
#        target.path = /usr/local/bin
        target.path = /home/sigve/Documents/SLETT/bin
        export(target.path)
    } else {
        target.path = $$DESTDIR
        export(target.path)
    }
    INSTALLS += target
}

message(DESTDIR: $$DESTDIR)
export(INSTALLS)
