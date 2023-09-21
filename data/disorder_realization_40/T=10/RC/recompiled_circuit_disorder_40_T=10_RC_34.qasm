OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5326795) q[0];
sx q[0];
rz(-2.764954) q[0];
sx q[0];
rz(-0.11178804) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(4.7987727) q[1];
sx q[1];
rz(6.12943) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.789334) q[0];
sx q[0];
rz(-0.52868045) q[0];
sx q[0];
rz(-2.5369011) q[0];
rz(-pi) q[1];
rz(-1.2674238) q[2];
sx q[2];
rz(-0.25250013) q[2];
sx q[2];
rz(-1.4361824) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.598688) q[1];
sx q[1];
rz(-2.821327) q[1];
sx q[1];
rz(-1.6673052) q[1];
rz(0.33711707) q[3];
sx q[3];
rz(-1.9883336) q[3];
sx q[3];
rz(-1.9330213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.443632) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(1.4367746) q[2];
rz(-0.73389655) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(-2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2519418) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(2.2170128) q[0];
rz(-0.997116) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(-1.3234214) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2307777) q[0];
sx q[0];
rz(-1.1417023) q[0];
sx q[0];
rz(0.50165117) q[0];
x q[1];
rz(-1.6152482) q[2];
sx q[2];
rz(-1.9200846) q[2];
sx q[2];
rz(2.8291707) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9604608) q[1];
sx q[1];
rz(-0.84003969) q[1];
sx q[1];
rz(-1.9366656) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1341392) q[3];
sx q[3];
rz(-0.65348071) q[3];
sx q[3];
rz(-0.37936488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.20415846) q[2];
sx q[2];
rz(-1.5519451) q[2];
sx q[2];
rz(-0.75817529) q[2];
rz(-0.6289064) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(1.988407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4085098) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(-0.68840233) q[0];
rz(0.06772659) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(0.53007954) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90572689) q[0];
sx q[0];
rz(-1.6775963) q[0];
sx q[0];
rz(1.8624767) q[0];
rz(-pi) q[1];
rz(2.7128503) q[2];
sx q[2];
rz(-2.0566166) q[2];
sx q[2];
rz(0.57810099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59280076) q[1];
sx q[1];
rz(-0.93344102) q[1];
sx q[1];
rz(0.42216502) q[1];
rz(-1.1881371) q[3];
sx q[3];
rz(-0.2573765) q[3];
sx q[3];
rz(2.5352258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80785859) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(-0.90144908) q[2];
rz(0.83550134) q[3];
sx q[3];
rz(-1.6136026) q[3];
sx q[3];
rz(1.3114595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9505342) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(2.3572671) q[0];
rz(3.0803608) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(-3.004946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89350677) q[0];
sx q[0];
rz(-0.49960217) q[0];
sx q[0];
rz(-1.8462371) q[0];
rz(-pi) q[1];
x q[1];
rz(1.253445) q[2];
sx q[2];
rz(-2.0145406) q[2];
sx q[2];
rz(1.2197942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.54388753) q[1];
sx q[1];
rz(-1.6391616) q[1];
sx q[1];
rz(-0.023936546) q[1];
rz(-3.0487719) q[3];
sx q[3];
rz(-0.99675677) q[3];
sx q[3];
rz(1.012158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1293929) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(0.56048918) q[2];
rz(0.012332049) q[3];
sx q[3];
rz(-2.2380232) q[3];
sx q[3];
rz(1.0906609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085768) q[0];
sx q[0];
rz(-0.60537678) q[0];
sx q[0];
rz(0.82114712) q[0];
rz(-0.87617809) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(-1.7339773) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53769892) q[0];
sx q[0];
rz(-1.961381) q[0];
sx q[0];
rz(-1.1917398) q[0];
rz(-2.1528835) q[2];
sx q[2];
rz(-1.8277797) q[2];
sx q[2];
rz(-0.35494057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.69406063) q[1];
sx q[1];
rz(-0.58528712) q[1];
sx q[1];
rz(-1.8156169) q[1];
rz(-1.0293343) q[3];
sx q[3];
rz(-0.41356219) q[3];
sx q[3];
rz(-3.001861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6901107) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(-0.042479854) q[2];
rz(-0.63043198) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(-2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(0.4831627) q[0];
rz(-2.0893611) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(-0.56484708) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0676346) q[0];
sx q[0];
rz(-1.7059776) q[0];
sx q[0];
rz(-1.0413175) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9916612) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(1.3644621) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.40183345) q[1];
sx q[1];
rz(-1.9976915) q[1];
sx q[1];
rz(0.065211936) q[1];
rz(-2.5552093) q[3];
sx q[3];
rz(-2.8660503) q[3];
sx q[3];
rz(-1.1875718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.57006449) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(1.8035536) q[2];
rz(1.3048874) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(-2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17334443) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(-0.54779732) q[0];
rz(-2.3563747) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(-2.8731667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7086605) q[0];
sx q[0];
rz(-1.6527358) q[0];
sx q[0];
rz(-2.8952778) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3586876) q[2];
sx q[2];
rz(-2.1856538) q[2];
sx q[2];
rz(2.8564786) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0722326) q[1];
sx q[1];
rz(-1.6688804) q[1];
sx q[1];
rz(-0.98274883) q[1];
x q[2];
rz(-0.3018474) q[3];
sx q[3];
rz(-1.0511304) q[3];
sx q[3];
rz(0.93320751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.61775529) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(-2.7653149) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58586621) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(-2.1222173) q[0];
rz(-2.2881919) q[1];
sx q[1];
rz(-1.9995721) q[1];
sx q[1];
rz(0.44874915) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73483) q[0];
sx q[0];
rz(-1.3376298) q[0];
sx q[0];
rz(1.083311) q[0];
rz(-pi) q[1];
rz(-0.86782311) q[2];
sx q[2];
rz(-2.7542369) q[2];
sx q[2];
rz(-0.58434904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1467421) q[1];
sx q[1];
rz(-2.5454306) q[1];
sx q[1];
rz(-0.46390987) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87047808) q[3];
sx q[3];
rz(-2.212489) q[3];
sx q[3];
rz(0.93922797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2723508) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(0.31420079) q[2];
rz(-0.82434404) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(0.79469386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0712414) q[0];
sx q[0];
rz(-3.0817139) q[0];
sx q[0];
rz(1.2605793) q[0];
rz(0.64385995) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(3.1226645) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8185794) q[0];
sx q[0];
rz(-1.5529263) q[0];
sx q[0];
rz(0.26364003) q[0];
x q[1];
rz(2.9855965) q[2];
sx q[2];
rz(-1.4737355) q[2];
sx q[2];
rz(2.7023466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70871204) q[1];
sx q[1];
rz(-1.6212665) q[1];
sx q[1];
rz(-2.9570079) q[1];
rz(-pi) q[2];
rz(-2.4631259) q[3];
sx q[3];
rz(-0.36558357) q[3];
sx q[3];
rz(1.4640704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5921322) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(2.4712759) q[2];
rz(0.51236764) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(1.4204773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55384127) q[0];
sx q[0];
rz(-1.2441664) q[0];
sx q[0];
rz(2.642139) q[0];
rz(-1.5746501) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(-2.0589028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2912746) q[0];
sx q[0];
rz(-0.58915888) q[0];
sx q[0];
rz(-3.0535112) q[0];
x q[1];
rz(0.81929368) q[2];
sx q[2];
rz(-1.594992) q[2];
sx q[2];
rz(1.0011315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6607099) q[1];
sx q[1];
rz(-2.1143267) q[1];
sx q[1];
rz(-1.5927614) q[1];
x q[2];
rz(1.8977676) q[3];
sx q[3];
rz(-1.9332814) q[3];
sx q[3];
rz(-2.8231951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3988951) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(2.6297074) q[2];
rz(2.7486457) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(-2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95505161) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(2.4004249) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(-1.0319866) q[2];
sx q[2];
rz(-2.451755) q[2];
sx q[2];
rz(-2.3103726) q[2];
rz(-0.74647222) q[3];
sx q[3];
rz(-1.6975879) q[3];
sx q[3];
rz(-3.0293037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
