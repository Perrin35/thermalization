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
rz(1.0526429) q[1];
sx q[1];
rz(3.9349603) q[1];
sx q[1];
rz(10.959672) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14102916) q[0];
sx q[0];
rz(-1.805843) q[0];
sx q[0];
rz(-1.1147333) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1821177) q[2];
sx q[2];
rz(-2.3797894) q[2];
sx q[2];
rz(-2.3488059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1032627) q[1];
sx q[1];
rz(-2.1974149) q[1];
sx q[1];
rz(-0.4396529) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9382557) q[3];
sx q[3];
rz(-1.4383738) q[3];
sx q[3];
rz(-1.2156957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4172998) q[2];
sx q[2];
rz(-2.5827926) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7928829) q[0];
sx q[0];
rz(-3.0443158) q[0];
sx q[0];
rz(-1.2874228) q[0];
rz(-0.48311326) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(-1.2189254) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12139509) q[0];
sx q[0];
rz(-1.4123865) q[0];
sx q[0];
rz(1.3729457) q[0];
rz(-pi) q[1];
rz(-2.4697815) q[2];
sx q[2];
rz(-1.3992568) q[2];
sx q[2];
rz(1.951527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2698631) q[1];
sx q[1];
rz(-0.31530373) q[1];
sx q[1];
rz(1.0878776) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.990149) q[3];
sx q[3];
rz(-2.4046728) q[3];
sx q[3];
rz(1.1518971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37811849) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(-2.5197869) q[2];
rz(-2.5032737) q[3];
sx q[3];
rz(-0.74898762) q[3];
sx q[3];
rz(-2.2973255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8353552) q[0];
sx q[0];
rz(-0.64202809) q[0];
sx q[0];
rz(-0.21632347) q[0];
rz(2.5430039) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(-2.6187706) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.110958) q[0];
sx q[0];
rz(-1.8930264) q[0];
sx q[0];
rz(-0.3229753) q[0];
rz(0.39117809) q[2];
sx q[2];
rz(-1.7157946) q[2];
sx q[2];
rz(0.36996335) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4292517) q[1];
sx q[1];
rz(-1.05637) q[1];
sx q[1];
rz(-1.8271258) q[1];
x q[2];
rz(0.30662068) q[3];
sx q[3];
rz(-1.2546179) q[3];
sx q[3];
rz(-0.090360377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9598976) q[2];
sx q[2];
rz(-1.8671702) q[2];
sx q[2];
rz(-1.5175021) q[2];
rz(-2.6767139) q[3];
sx q[3];
rz(-2.2113694) q[3];
sx q[3];
rz(-1.6200861) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14438039) q[0];
sx q[0];
rz(-0.388044) q[0];
sx q[0];
rz(-1.602518) q[0];
rz(2.1082361) q[1];
sx q[1];
rz(-0.76834232) q[1];
sx q[1];
rz(-1.296952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.715623) q[0];
sx q[0];
rz(-0.85943009) q[0];
sx q[0];
rz(-0.27669669) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4960852) q[2];
sx q[2];
rz(-0.69715188) q[2];
sx q[2];
rz(2.8139204) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.32758157) q[1];
sx q[1];
rz(-1.9821229) q[1];
sx q[1];
rz(-0.41366215) q[1];
x q[2];
rz(0.4850895) q[3];
sx q[3];
rz(-2.5311433) q[3];
sx q[3];
rz(-1.7318673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79232717) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(2.0959334) q[2];
rz(2.9271434) q[3];
sx q[3];
rz(-2.4172473) q[3];
sx q[3];
rz(-1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8056718) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(0.51934284) q[0];
rz(-0.94201159) q[1];
sx q[1];
rz(-0.28128925) q[1];
sx q[1];
rz(-1.6090144) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.518884) q[0];
sx q[0];
rz(-1.9340589) q[0];
sx q[0];
rz(0.69464442) q[0];
rz(-pi) q[1];
rz(-0.66560676) q[2];
sx q[2];
rz(-2.9133248) q[2];
sx q[2];
rz(-0.98787243) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11621257) q[1];
sx q[1];
rz(-2.0713701) q[1];
sx q[1];
rz(3.1356372) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4970826) q[2];
sx q[2];
rz(-2.3641059) q[2];
sx q[2];
rz(0.28437781) q[2];
rz(-0.55822462) q[3];
sx q[3];
rz(-1.3642718) q[3];
sx q[3];
rz(-0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068950653) q[0];
sx q[0];
rz(-2.729029) q[0];
sx q[0];
rz(2.5010342) q[0];
rz(1.6259646) q[1];
sx q[1];
rz(-2.0104355) q[1];
sx q[1];
rz(-2.9277149) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2757077) q[0];
sx q[0];
rz(-1.8731706) q[0];
sx q[0];
rz(2.9024966) q[0];
x q[1];
rz(-0.40753813) q[2];
sx q[2];
rz(-2.354051) q[2];
sx q[2];
rz(2.119273) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7678309) q[1];
sx q[1];
rz(-1.5015748) q[1];
sx q[1];
rz(0.94061416) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0701105) q[3];
sx q[3];
rz(-2.0897647) q[3];
sx q[3];
rz(-1.9148358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(-2.9638929) q[2];
rz(-1.3972345) q[3];
sx q[3];
rz(-1.1763562) q[3];
sx q[3];
rz(-0.41745225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099667065) q[0];
sx q[0];
rz(-0.61138994) q[0];
sx q[0];
rz(2.0360816) q[0];
rz(-0.28327495) q[1];
sx q[1];
rz(-2.6073644) q[1];
sx q[1];
rz(-1.4615321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7992226) q[0];
sx q[0];
rz(-2.2307751) q[0];
sx q[0];
rz(-0.8578542) q[0];
rz(2.2143557) q[2];
sx q[2];
rz(-0.74802784) q[2];
sx q[2];
rz(-0.03803703) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9859559) q[1];
sx q[1];
rz(-1.0105304) q[1];
sx q[1];
rz(1.7079279) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5652065) q[3];
sx q[3];
rz(-1.8462876) q[3];
sx q[3];
rz(2.798693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0561515) q[2];
sx q[2];
rz(-1.5957007) q[2];
sx q[2];
rz(-0.20480569) q[2];
rz(-2.0917995) q[3];
sx q[3];
rz(-0.27725163) q[3];
sx q[3];
rz(-2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3664704) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(0.033576641) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-0.74445236) q[1];
sx q[1];
rz(-1.144145) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9947203) q[0];
sx q[0];
rz(-2.0573186) q[0];
sx q[0];
rz(1.5994607) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0424007) q[2];
sx q[2];
rz(-1.2481523) q[2];
sx q[2];
rz(1.9930372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.42061372) q[1];
sx q[1];
rz(-2.3434964) q[1];
sx q[1];
rz(0.53687232) q[1];
rz(-pi) q[2];
rz(1.9162634) q[3];
sx q[3];
rz(-1.8356712) q[3];
sx q[3];
rz(3.1136852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(1.8899274) q[2];
rz(-1.3517316) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(-0.98021093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0018472483) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(-1.8852604) q[0];
rz(0.095257692) q[1];
sx q[1];
rz(-0.88905159) q[1];
sx q[1];
rz(2.6108066) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4149218) q[0];
sx q[0];
rz(-1.4106531) q[0];
sx q[0];
rz(1.3199575) q[0];
rz(2.8235648) q[2];
sx q[2];
rz(-1.5340759) q[2];
sx q[2];
rz(-1.6317692) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2797045) q[1];
sx q[1];
rz(-1.3315531) q[1];
sx q[1];
rz(1.4390535) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6952031) q[3];
sx q[3];
rz(-1.4671031) q[3];
sx q[3];
rz(0.54561347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.329616) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(-2.6007593) q[2];
rz(2.3234308) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5166017) q[0];
sx q[0];
rz(-1.7778968) q[0];
sx q[0];
rz(2.7428108) q[0];
rz(1.7575691) q[1];
sx q[1];
rz(-2.3868491) q[1];
sx q[1];
rz(-0.049364518) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20959768) q[0];
sx q[0];
rz(-1.5100129) q[0];
sx q[0];
rz(1.6796735) q[0];
x q[1];
rz(-0.10714679) q[2];
sx q[2];
rz(-1.8320036) q[2];
sx q[2];
rz(-0.10814737) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80222622) q[1];
sx q[1];
rz(-2.5605533) q[1];
sx q[1];
rz(1.4895579) q[1];
x q[2];
rz(-1.3672164) q[3];
sx q[3];
rz(-2.1927811) q[3];
sx q[3];
rz(-0.79727117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0177239) q[2];
sx q[2];
rz(-1.3088635) q[2];
sx q[2];
rz(1.5105985) q[2];
rz(0.92783582) q[3];
sx q[3];
rz(-1.8001013) q[3];
sx q[3];
rz(-1.235289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(1.0646959) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(2.627768) q[2];
sx q[2];
rz(-2.0110301) q[2];
sx q[2];
rz(-0.55883519) q[2];
rz(-1.0743027) q[3];
sx q[3];
rz(-2.5357694) q[3];
sx q[3];
rz(2.3800935) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
