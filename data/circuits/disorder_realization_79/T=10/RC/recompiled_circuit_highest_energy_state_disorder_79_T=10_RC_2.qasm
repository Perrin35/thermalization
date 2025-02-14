OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.41843721) q[0];
sx q[0];
rz(-0.96324459) q[0];
sx q[0];
rz(0.20382398) q[0];
rz(-2.0889497) q[1];
sx q[1];
rz(-0.79336762) q[1];
sx q[1];
rz(-1.5348943) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8730072) q[0];
sx q[0];
rz(-2.632336) q[0];
sx q[0];
rz(-1.0727706) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7947172) q[2];
sx q[2];
rz(-2.2636608) q[2];
sx q[2];
rz(-1.3077867) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4272449) q[1];
sx q[1];
rz(-2.3934919) q[1];
sx q[1];
rz(-1.0393049) q[1];
x q[2];
rz(-0.20333692) q[3];
sx q[3];
rz(-1.4383738) q[3];
sx q[3];
rz(1.9258969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7242929) q[2];
sx q[2];
rz(-0.55880004) q[2];
sx q[2];
rz(-2.144045) q[2];
rz(2.3857462) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(1.0409748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7928829) q[0];
sx q[0];
rz(-0.097276874) q[0];
sx q[0];
rz(-1.2874228) q[0];
rz(-0.48311326) q[1];
sx q[1];
rz(-1.0419507) q[1];
sx q[1];
rz(1.2189254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.481015) q[0];
sx q[0];
rz(-1.3754551) q[0];
sx q[0];
rz(-0.16150772) q[0];
rz(2.8701209) q[2];
sx q[2];
rz(-2.4515332) q[2];
sx q[2];
rz(0.59218237) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.87172958) q[1];
sx q[1];
rz(-0.31530373) q[1];
sx q[1];
rz(2.053715) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2628701) q[3];
sx q[3];
rz(-1.2936397) q[3];
sx q[3];
rz(-0.73776484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37811849) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(-2.5197869) q[2];
rz(-2.5032737) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(-0.84426713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3062375) q[0];
sx q[0];
rz(-2.4995646) q[0];
sx q[0];
rz(-0.21632347) q[0];
rz(0.59858876) q[1];
sx q[1];
rz(-1.8363771) q[1];
sx q[1];
rz(0.52282202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3591839) q[0];
sx q[0];
rz(-2.6893927) q[0];
sx q[0];
rz(2.3308999) q[0];
rz(-0.39117809) q[2];
sx q[2];
rz(-1.4257981) q[2];
sx q[2];
rz(0.36996335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2017858) q[1];
sx q[1];
rz(-2.5720189) q[1];
sx q[1];
rz(-0.42167432) q[1];
rz(-pi) q[2];
x q[2];
rz(2.834972) q[3];
sx q[3];
rz(-1.2546179) q[3];
sx q[3];
rz(-3.0512323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.181695) q[2];
sx q[2];
rz(-1.2744224) q[2];
sx q[2];
rz(1.6240906) q[2];
rz(-0.46487871) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(1.5215065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9972123) q[0];
sx q[0];
rz(-0.388044) q[0];
sx q[0];
rz(1.602518) q[0];
rz(-2.1082361) q[1];
sx q[1];
rz(-0.76834232) q[1];
sx q[1];
rz(-1.8446406) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015681277) q[0];
sx q[0];
rz(-2.3871579) q[0];
sx q[0];
rz(1.2638646) q[0];
x q[1];
rz(0.87501915) q[2];
sx q[2];
rz(-1.6187374) q[2];
sx q[2];
rz(-1.9557916) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9819543) q[1];
sx q[1];
rz(-0.57483638) q[1];
sx q[1];
rz(-0.82623085) q[1];
x q[2];
rz(-0.4850895) q[3];
sx q[3];
rz(-0.61044932) q[3];
sx q[3];
rz(-1.7318673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3492655) q[2];
sx q[2];
rz(-0.63750625) q[2];
sx q[2];
rz(1.0456592) q[2];
rz(2.9271434) q[3];
sx q[3];
rz(-0.72434536) q[3];
sx q[3];
rz(1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8056718) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(0.51934284) q[0];
rz(0.94201159) q[1];
sx q[1];
rz(-2.8603034) q[1];
sx q[1];
rz(-1.6090144) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.596622) q[0];
sx q[0];
rz(-0.76966296) q[0];
sx q[0];
rz(-0.53588698) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66560676) q[2];
sx q[2];
rz(-0.22826787) q[2];
sx q[2];
rz(0.98787243) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11621257) q[1];
sx q[1];
rz(-1.0702225) q[1];
sx q[1];
rz(0.0059555014) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4761837) q[3];
sx q[3];
rz(-1.4615371) q[3];
sx q[3];
rz(-3.008568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.64451009) q[2];
sx q[2];
rz(-0.77748674) q[2];
sx q[2];
rz(0.28437781) q[2];
rz(-0.55822462) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.072642) q[0];
sx q[0];
rz(-0.41256368) q[0];
sx q[0];
rz(0.64055842) q[0];
rz(1.5156281) q[1];
sx q[1];
rz(-1.1311572) q[1];
sx q[1];
rz(-2.9277149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17992526) q[0];
sx q[0];
rz(-2.7583987) q[0];
sx q[0];
rz(-2.2201594) q[0];
x q[1];
rz(2.7340545) q[2];
sx q[2];
rz(-0.78754163) q[2];
sx q[2];
rz(-2.119273) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37376172) q[1];
sx q[1];
rz(-1.5015748) q[1];
sx q[1];
rz(2.2009785) q[1];
rz(1.6951896) q[3];
sx q[3];
rz(-2.6181707) q[3];
sx q[3];
rz(-2.058213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3512257) q[2];
sx q[2];
rz(-2.6595317) q[2];
sx q[2];
rz(-0.17769979) q[2];
rz(1.7443582) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(0.41745225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0419256) q[0];
sx q[0];
rz(-2.5302027) q[0];
sx q[0];
rz(2.0360816) q[0];
rz(0.28327495) q[1];
sx q[1];
rz(-0.53422821) q[1];
sx q[1];
rz(1.6800605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38835907) q[0];
sx q[0];
rz(-0.93030158) q[0];
sx q[0];
rz(-2.4413013) q[0];
rz(-0.5080451) q[2];
sx q[2];
rz(-2.1461502) q[2];
sx q[2];
rz(-2.3066556) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0431598) q[1];
sx q[1];
rz(-2.5665356) q[1];
sx q[1];
rz(2.9270323) q[1];
x q[2];
rz(0.47847943) q[3];
sx q[3];
rz(-0.63204403) q[3];
sx q[3];
rz(1.6242336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0854411) q[2];
sx q[2];
rz(-1.5957007) q[2];
sx q[2];
rz(-2.936787) q[2];
rz(1.0497931) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(-0.36627305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3664704) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(3.108016) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-2.3971403) q[1];
sx q[1];
rz(-1.9974476) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0856253) q[0];
sx q[0];
rz(-2.6542943) q[0];
sx q[0];
rz(3.0874599) q[0];
rz(-pi) q[1];
rz(-3.0424007) q[2];
sx q[2];
rz(-1.8934403) q[2];
sx q[2];
rz(-1.1485554) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7209789) q[1];
sx q[1];
rz(-2.3434964) q[1];
sx q[1];
rz(-2.6047203) q[1];
rz(-0.28067067) q[3];
sx q[3];
rz(-1.9037399) q[3];
sx q[3];
rz(1.6368293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(-1.2516652) q[2];
rz(1.7898611) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(2.1613817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1397454) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(1.2563323) q[0];
rz(0.095257692) q[1];
sx q[1];
rz(-0.88905159) q[1];
sx q[1];
rz(-0.5307861) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7266709) q[0];
sx q[0];
rz(-1.4106531) q[0];
sx q[0];
rz(1.3199575) q[0];
rz(1.6094535) q[2];
sx q[2];
rz(-1.2529904) q[2];
sx q[2];
rz(3.0927049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7698313) q[1];
sx q[1];
rz(-2.8690845) q[1];
sx q[1];
rz(2.647577) q[1];
rz(-2.9050499) q[3];
sx q[3];
rz(-2.6841087) q[3];
sx q[3];
rz(1.903423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8119767) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(-0.54083332) q[2];
rz(2.3234308) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(-0.63280672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.624991) q[0];
sx q[0];
rz(-1.3636959) q[0];
sx q[0];
rz(-2.7428108) q[0];
rz(-1.7575691) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(-0.049364518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.787034) q[0];
sx q[0];
rz(-1.6794717) q[0];
sx q[0];
rz(3.0804481) q[0];
rz(-pi) q[1];
rz(-1.9513543) q[2];
sx q[2];
rz(-2.8597288) q[2];
sx q[2];
rz(-2.6388002) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80222622) q[1];
sx q[1];
rz(-2.5605533) q[1];
sx q[1];
rz(1.4895579) q[1];
rz(-pi) q[2];
rz(-2.8667198) q[3];
sx q[3];
rz(-2.4913553) q[3];
sx q[3];
rz(-2.6848328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.12386879) q[2];
sx q[2];
rz(-1.3088635) q[2];
sx q[2];
rz(-1.5105985) q[2];
rz(2.2137568) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(-1.235289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8265726) q[0];
sx q[0];
rz(-0.44354225) q[0];
sx q[0];
rz(-1.0916239) q[0];
rz(2.0768968) q[1];
sx q[1];
rz(-0.5263435) q[1];
sx q[1];
rz(1.975504) q[1];
rz(-0.76416107) q[2];
sx q[2];
rz(-2.4780826) q[2];
sx q[2];
rz(1.6586951) q[2];
rz(0.31872411) q[3];
sx q[3];
rz(-1.0464077) q[3];
sx q[3];
rz(2.9628021) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
