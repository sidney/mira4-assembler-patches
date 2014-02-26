

#ifdef PUBLICQUIET
#define MLV_QUIET
#else
#define MLV_QUIET "(debugmsg) "
#endif

#ifdef TIMERESTRICTED
#define MLV_TIME "(time restricted) "
#else
#define MLV_TIME
#endif


#define GITBUILD "4.0_0+gea25053"
#define MIRALIBVERSION GITBUILD
#define MIRAVERSION MIRALIBVERSION " " MLV_QUIET MLV_TIME""
