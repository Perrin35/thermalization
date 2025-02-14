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
rz(-1.7186681) q[0];
sx q[0];
rz(-1.0942425) q[0];
sx q[0];
rz(-2.8835468) q[0];
rz(-0.99335042) q[1];
sx q[1];
rz(-1.8408096) q[1];
sx q[1];
rz(2.8859477) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6171744) q[0];
sx q[0];
rz(-0.99992311) q[0];
sx q[0];
rz(0.85321315) q[0];
rz(-pi) q[1];
rz(2.9886888) q[2];
sx q[2];
rz(-0.39275482) q[2];
sx q[2];
rz(-3.1043579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8430184) q[1];
sx q[1];
rz(-2.9461423) q[1];
sx q[1];
rz(-1.2091544) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7624038) q[3];
sx q[3];
rz(-2.6122836) q[3];
sx q[3];
rz(2.3922753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.51097441) q[2];
sx q[2];
rz(-1.3773842) q[2];
sx q[2];
rz(-2.0261436) q[2];
rz(-0.014178064) q[3];
sx q[3];
rz(-1.7970128) q[3];
sx q[3];
rz(-2.7052963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9344591) q[0];
sx q[0];
rz(-0.62196982) q[0];
sx q[0];
rz(-1.2472664) q[0];
rz(2.6994052) q[1];
sx q[1];
rz(-1.7555883) q[1];
sx q[1];
rz(-0.8173379) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2674448) q[0];
sx q[0];
rz(-2.5782881) q[0];
sx q[0];
rz(-1.9356217) q[0];
rz(2.0109315) q[2];
sx q[2];
rz(-1.3190184) q[2];
sx q[2];
rz(1.8735069) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8838447) q[1];
sx q[1];
rz(-0.93755975) q[1];
sx q[1];
rz(1.7764938) q[1];
rz(-2.2226187) q[3];
sx q[3];
rz(-2.2187244) q[3];
sx q[3];
rz(1.0046665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77357972) q[2];
sx q[2];
rz(-2.7648338) q[2];
sx q[2];
rz(2.9212941) q[2];
rz(1.1822654) q[3];
sx q[3];
rz(-1.1743098) q[3];
sx q[3];
rz(-0.063145414) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0910864) q[0];
sx q[0];
rz(-1.3171221) q[0];
sx q[0];
rz(-2.0029946) q[0];
rz(2.7298722) q[1];
sx q[1];
rz(-2.3717334) q[1];
sx q[1];
rz(-1.0134816) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0419755) q[0];
sx q[0];
rz(-1.8238153) q[0];
sx q[0];
rz(1.5554886) q[0];
x q[1];
rz(-1.1297497) q[2];
sx q[2];
rz(-1.6905606) q[2];
sx q[2];
rz(-0.83415595) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0358082) q[1];
sx q[1];
rz(-0.86584751) q[1];
sx q[1];
rz(0.99143274) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7168379) q[3];
sx q[3];
rz(-0.73196326) q[3];
sx q[3];
rz(0.41106658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1608405) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(1.8864924) q[2];
rz(-0.27215019) q[3];
sx q[3];
rz(-2.3641219) q[3];
sx q[3];
rz(0.14061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18144064) q[0];
sx q[0];
rz(-2.20521) q[0];
sx q[0];
rz(2.5312359) q[0];
rz(-2.4329674) q[1];
sx q[1];
rz(-1.0162063) q[1];
sx q[1];
rz(1.5707387) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40442586) q[0];
sx q[0];
rz(-0.20658399) q[0];
sx q[0];
rz(2.4086359) q[0];
x q[1];
rz(1.4518634) q[2];
sx q[2];
rz(-1.2858675) q[2];
sx q[2];
rz(-0.97078427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8766986) q[1];
sx q[1];
rz(-1.8421409) q[1];
sx q[1];
rz(1.62074) q[1];
x q[2];
rz(2.3409178) q[3];
sx q[3];
rz(-1.8599556) q[3];
sx q[3];
rz(-1.9137933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9307956) q[2];
sx q[2];
rz(-1.0127298) q[2];
sx q[2];
rz(0.55541682) q[2];
rz(-0.72426116) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(0.65653062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.637218) q[0];
sx q[0];
rz(-1.9509622) q[0];
sx q[0];
rz(2.7727238) q[0];
rz(0.24636191) q[1];
sx q[1];
rz(-1.8135704) q[1];
sx q[1];
rz(1.7074283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1911538) q[0];
sx q[0];
rz(-2.3847347) q[0];
sx q[0];
rz(0.0075154742) q[0];
x q[1];
rz(-2.6961961) q[2];
sx q[2];
rz(-1.7295618) q[2];
sx q[2];
rz(-2.7743055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7692657) q[1];
sx q[1];
rz(-2.7397836) q[1];
sx q[1];
rz(-3.1308082) q[1];
x q[2];
rz(3.0415972) q[3];
sx q[3];
rz(-1.7971562) q[3];
sx q[3];
rz(-0.55734998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4904334) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(-1.0820214) q[2];
rz(-2.3467482) q[3];
sx q[3];
rz(-1.669603) q[3];
sx q[3];
rz(1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27570462) q[0];
sx q[0];
rz(-0.10144932) q[0];
sx q[0];
rz(-2.8979982) q[0];
rz(-0.98681915) q[1];
sx q[1];
rz(-2.3934264) q[1];
sx q[1];
rz(-2.2056244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3098738) q[0];
sx q[0];
rz(-1.1748992) q[0];
sx q[0];
rz(2.7926366) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7924839) q[2];
sx q[2];
rz(-2.287068) q[2];
sx q[2];
rz(-2.8385065) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82560357) q[1];
sx q[1];
rz(-0.88615075) q[1];
sx q[1];
rz(2.0988093) q[1];
x q[2];
rz(0.021401568) q[3];
sx q[3];
rz(-2.170522) q[3];
sx q[3];
rz(-0.76364005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4794856) q[2];
sx q[2];
rz(-1.138849) q[2];
sx q[2];
rz(2.7916419) q[2];
rz(-0.55772603) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(0.85132712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6996985) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(-2.8398474) q[0];
rz(2.8217577) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(-1.1357657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0253191) q[0];
sx q[0];
rz(-0.92705446) q[0];
sx q[0];
rz(-1.6380861) q[0];
rz(2.202353) q[2];
sx q[2];
rz(-1.9983091) q[2];
sx q[2];
rz(-1.3434501) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3150683) q[1];
sx q[1];
rz(-0.60501912) q[1];
sx q[1];
rz(1.4304763) q[1];
rz(-pi) q[2];
rz(0.11140996) q[3];
sx q[3];
rz(-0.79395959) q[3];
sx q[3];
rz(0.68796989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7319506) q[2];
sx q[2];
rz(-0.66764098) q[2];
sx q[2];
rz(0.14275924) q[2];
rz(-1.3879294) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(2.4587542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9614354) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(-0.55602443) q[0];
rz(-2.9122638) q[1];
sx q[1];
rz(-1.8792968) q[1];
sx q[1];
rz(-1.3831327) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97998226) q[0];
sx q[0];
rz(-1.605827) q[0];
sx q[0];
rz(-0.73765124) q[0];
rz(3.0223614) q[2];
sx q[2];
rz(-2.6802353) q[2];
sx q[2];
rz(-2.2329494) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.759093) q[1];
sx q[1];
rz(-1.9803932) q[1];
sx q[1];
rz(0.31633693) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69139685) q[3];
sx q[3];
rz(-0.56400877) q[3];
sx q[3];
rz(2.9897652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6250299) q[2];
sx q[2];
rz(-0.27250686) q[2];
sx q[2];
rz(1.2410835) q[2];
rz(-3.0789913) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(-2.3144498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1271707) q[0];
sx q[0];
rz(-1.7272471) q[0];
sx q[0];
rz(2.993809) q[0];
rz(-2.6328909) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(-2.7630189) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050178278) q[0];
sx q[0];
rz(-1.1568406) q[0];
sx q[0];
rz(-1.3148091) q[0];
x q[1];
rz(2.1630387) q[2];
sx q[2];
rz(-1.8752974) q[2];
sx q[2];
rz(-1.026012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7853773) q[1];
sx q[1];
rz(-1.2357724) q[1];
sx q[1];
rz(0.32932333) q[1];
rz(-pi) q[2];
rz(-1.6949953) q[3];
sx q[3];
rz(-0.81362766) q[3];
sx q[3];
rz(-2.1424993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1754237) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(-1.2790722) q[2];
rz(-1.7133948) q[3];
sx q[3];
rz(-1.0019852) q[3];
sx q[3];
rz(-0.14555791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7487504) q[0];
sx q[0];
rz(-2.5635283) q[0];
sx q[0];
rz(0.064099126) q[0];
rz(2.6129258) q[1];
sx q[1];
rz(-1.8730947) q[1];
sx q[1];
rz(-2.166523) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.753016) q[0];
sx q[0];
rz(-2.574769) q[0];
sx q[0];
rz(0.50811572) q[0];
rz(-pi) q[1];
rz(-2.1553173) q[2];
sx q[2];
rz(-1.0800252) q[2];
sx q[2];
rz(1.402439) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31723695) q[1];
sx q[1];
rz(-1.1350766) q[1];
sx q[1];
rz(1.812326) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3417792) q[3];
sx q[3];
rz(-1.3327362) q[3];
sx q[3];
rz(0.89512596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9958682) q[2];
sx q[2];
rz(-2.4805562) q[2];
sx q[2];
rz(-0.60689849) q[2];
rz(0.25035826) q[3];
sx q[3];
rz(-1.4162049) q[3];
sx q[3];
rz(-2.6355766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.09457) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(0.49867123) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(-0.86317369) q[2];
sx q[2];
rz(-2.045608) q[2];
sx q[2];
rz(1.506293) q[2];
rz(0.12228431) q[3];
sx q[3];
rz(-1.3663843) q[3];
sx q[3];
rz(-2.1257675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
