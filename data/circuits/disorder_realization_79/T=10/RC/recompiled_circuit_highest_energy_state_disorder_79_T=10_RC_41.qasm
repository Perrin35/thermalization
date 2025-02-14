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
rz(-2.7231554) q[0];
sx q[0];
rz(-2.1783481) q[0];
sx q[0];
rz(-0.20382398) q[0];
rz(-2.0889497) q[1];
sx q[1];
rz(-0.79336762) q[1];
sx q[1];
rz(-1.5348943) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2685854) q[0];
sx q[0];
rz(-2.632336) q[0];
sx q[0];
rz(1.0727706) q[0];
rz(-1.1821177) q[2];
sx q[2];
rz(-0.76180327) q[2];
sx q[2];
rz(2.3488059) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2633367) q[1];
sx q[1];
rz(-1.2188101) q[1];
sx q[1];
rz(-0.89603591) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4356218) q[3];
sx q[3];
rz(-1.7723284) q[3];
sx q[3];
rz(0.38231787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4172998) q[2];
sx q[2];
rz(-0.55880004) q[2];
sx q[2];
rz(0.99754769) q[2];
rz(2.3857462) q[3];
sx q[3];
rz(-1.3563124) q[3];
sx q[3];
rz(-1.0409748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34870979) q[0];
sx q[0];
rz(-0.097276874) q[0];
sx q[0];
rz(-1.8541699) q[0];
rz(2.6584794) q[1];
sx q[1];
rz(-1.0419507) q[1];
sx q[1];
rz(1.2189254) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.481015) q[0];
sx q[0];
rz(-1.7661375) q[0];
sx q[0];
rz(2.9800849) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67181113) q[2];
sx q[2];
rz(-1.3992568) q[2];
sx q[2];
rz(-1.951527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.375735) q[1];
sx q[1];
rz(-1.2925783) q[1];
sx q[1];
rz(-0.15032676) q[1];
x q[2];
rz(-0.87872259) q[3];
sx q[3];
rz(-1.847953) q[3];
sx q[3];
rz(-2.4038278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7634742) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(-2.5197869) q[2];
rz(-0.63831896) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(-2.2973255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3062375) q[0];
sx q[0];
rz(-2.4995646) q[0];
sx q[0];
rz(2.9252692) q[0];
rz(-0.59858876) q[1];
sx q[1];
rz(-1.8363771) q[1];
sx q[1];
rz(-0.52282202) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.110958) q[0];
sx q[0];
rz(-1.2485663) q[0];
sx q[0];
rz(2.8186174) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39117809) q[2];
sx q[2];
rz(-1.4257981) q[2];
sx q[2];
rz(-2.7716293) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93980689) q[1];
sx q[1];
rz(-2.5720189) q[1];
sx q[1];
rz(-0.42167432) q[1];
x q[2];
rz(1.9013635) q[3];
sx q[3];
rz(-1.8617612) q[3];
sx q[3];
rz(1.3823079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9598976) q[2];
sx q[2];
rz(-1.8671702) q[2];
sx q[2];
rz(-1.5175021) q[2];
rz(0.46487871) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(-1.5215065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14438039) q[0];
sx q[0];
rz(-2.7535487) q[0];
sx q[0];
rz(1.602518) q[0];
rz(1.0333565) q[1];
sx q[1];
rz(-0.76834232) q[1];
sx q[1];
rz(1.296952) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015681277) q[0];
sx q[0];
rz(-0.75443479) q[0];
sx q[0];
rz(1.8777281) q[0];
rz(-pi) q[1];
rz(-0.06242604) q[2];
sx q[2];
rz(-2.2656144) q[2];
sx q[2];
rz(0.42499396) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0721283) q[1];
sx q[1];
rz(-1.9481244) q[1];
sx q[1];
rz(-2.0153785) q[1];
x q[2];
rz(-0.55418535) q[3];
sx q[3];
rz(-1.8413787) q[3];
sx q[3];
rz(2.8949646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79232717) q[2];
sx q[2];
rz(-0.63750625) q[2];
sx q[2];
rz(1.0456592) q[2];
rz(0.21444923) q[3];
sx q[3];
rz(-0.72434536) q[3];
sx q[3];
rz(-1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8056718) q[0];
sx q[0];
rz(-1.8593973) q[0];
sx q[0];
rz(2.6222498) q[0];
rz(-0.94201159) q[1];
sx q[1];
rz(-2.8603034) q[1];
sx q[1];
rz(1.6090144) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62270861) q[0];
sx q[0];
rz(-1.9340589) q[0];
sx q[0];
rz(0.69464442) q[0];
rz(-0.18073323) q[2];
sx q[2];
rz(-1.4305947) q[2];
sx q[2];
rz(-3.0716346) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.457442) q[1];
sx q[1];
rz(-1.5760211) q[1];
sx q[1];
rz(-1.0702151) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10974613) q[3];
sx q[3];
rz(-1.4767495) q[3];
sx q[3];
rz(-1.7141683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4970826) q[2];
sx q[2];
rz(-2.3641059) q[2];
sx q[2];
rz(0.28437781) q[2];
rz(0.55822462) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(-0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.072642) q[0];
sx q[0];
rz(-2.729029) q[0];
sx q[0];
rz(2.5010342) q[0];
rz(-1.5156281) q[1];
sx q[1];
rz(-1.1311572) q[1];
sx q[1];
rz(-0.21387771) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77737264) q[0];
sx q[0];
rz(-1.7988482) q[0];
sx q[0];
rz(-1.8814726) q[0];
rz(0.40753813) q[2];
sx q[2];
rz(-0.78754163) q[2];
sx q[2];
rz(2.119273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2474413) q[1];
sx q[1];
rz(-0.94235984) q[1];
sx q[1];
rz(-3.0559866) q[1];
x q[2];
rz(1.446403) q[3];
sx q[3];
rz(-2.6181707) q[3];
sx q[3];
rz(-1.0833797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(2.9638929) q[2];
rz(1.3972345) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0419256) q[0];
sx q[0];
rz(-0.61138994) q[0];
sx q[0];
rz(1.1055111) q[0];
rz(-0.28327495) q[1];
sx q[1];
rz(-2.6073644) q[1];
sx q[1];
rz(-1.4615321) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7992226) q[0];
sx q[0];
rz(-0.91081753) q[0];
sx q[0];
rz(-0.8578542) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2093532) q[2];
sx q[2];
rz(-1.1503714) q[2];
sx q[2];
rz(-1.0300385) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9859559) q[1];
sx q[1];
rz(-2.1310622) q[1];
sx q[1];
rz(1.4336647) q[1];
x q[2];
rz(-1.8959778) q[3];
sx q[3];
rz(-1.0187314) q[3];
sx q[3];
rz(-2.0887041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0854411) q[2];
sx q[2];
rz(-1.5957007) q[2];
sx q[2];
rz(0.20480569) q[2];
rz(-2.0917995) q[3];
sx q[3];
rz(-0.27725163) q[3];
sx q[3];
rz(-2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751223) q[0];
sx q[0];
rz(-0.46638745) q[0];
sx q[0];
rz(-0.033576641) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-2.3971403) q[1];
sx q[1];
rz(1.144145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0559674) q[0];
sx q[0];
rz(-0.48729839) q[0];
sx q[0];
rz(-3.0874599) q[0];
x q[1];
rz(-0.099191908) q[2];
sx q[2];
rz(-1.2481523) q[2];
sx q[2];
rz(1.9930372) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.28539666) q[1];
sx q[1];
rz(-0.90803972) q[1];
sx q[1];
rz(-2.0539356) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2461478) q[3];
sx q[3];
rz(-0.43206462) q[3];
sx q[3];
rz(-0.91401446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3860151) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(-1.2516652) q[2];
rz(1.3517316) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(0.98021093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1397454) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(1.2563323) q[0];
rz(-3.046335) q[1];
sx q[1];
rz(-2.2525411) q[1];
sx q[1];
rz(-2.6108066) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40079257) q[0];
sx q[0];
rz(-2.8449028) q[0];
sx q[0];
rz(0.99389561) q[0];
rz(2.8235648) q[2];
sx q[2];
rz(-1.5340759) q[2];
sx q[2];
rz(-1.6317692) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4640749) q[1];
sx q[1];
rz(-1.6987659) q[1];
sx q[1];
rz(0.24125464) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9050499) q[3];
sx q[3];
rz(-0.45748392) q[3];
sx q[3];
rz(1.2381697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.329616) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(2.6007593) q[2];
rz(-2.3234308) q[3];
sx q[3];
rz(-1.4427789) q[3];
sx q[3];
rz(-0.63280672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5166017) q[0];
sx q[0];
rz(-1.3636959) q[0];
sx q[0];
rz(0.3987819) q[0];
rz(1.3840236) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(-0.049364518) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.931995) q[0];
sx q[0];
rz(-1.6315797) q[0];
sx q[0];
rz(1.6796735) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8334421) q[2];
sx q[2];
rz(-1.6742953) q[2];
sx q[2];
rz(-1.4348794) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4409742) q[1];
sx q[1];
rz(-1.6153533) q[1];
sx q[1];
rz(-0.99127165) q[1];
rz(-pi) q[2];
rz(0.63189854) q[3];
sx q[3];
rz(-1.4057341) q[3];
sx q[3];
rz(0.89323211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12386879) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(1.6309942) q[2];
rz(-2.2137568) q[3];
sx q[3];
rz(-1.8001013) q[3];
sx q[3];
rz(-1.235289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31502003) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(-2.0768968) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(-0.51382463) q[2];
sx q[2];
rz(-2.0110301) q[2];
sx q[2];
rz(-0.55883519) q[2];
rz(1.0743027) q[3];
sx q[3];
rz(-0.60582325) q[3];
sx q[3];
rz(-0.7614991) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
