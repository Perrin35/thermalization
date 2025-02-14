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
rz(1.2394387) q[0];
sx q[0];
rz(-1.8128938) q[0];
sx q[0];
rz(-0.21188307) q[0];
rz(1.7243241) q[1];
sx q[1];
rz(2.6091726) q[1];
sx q[1];
rz(6.6611023) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1013452) q[0];
sx q[0];
rz(-1.57124) q[0];
sx q[0];
rz(1.5612649) q[0];
x q[1];
rz(-1.2075519) q[2];
sx q[2];
rz(-0.9245199) q[2];
sx q[2];
rz(0.56528795) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18641414) q[1];
sx q[1];
rz(-0.70611533) q[1];
sx q[1];
rz(2.7104055) q[1];
x q[2];
rz(2.3842281) q[3];
sx q[3];
rz(-0.93776449) q[3];
sx q[3];
rz(-2.4773134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2892896) q[2];
sx q[2];
rz(-1.0785582) q[2];
sx q[2];
rz(3.0618073) q[2];
rz(0.96528178) q[3];
sx q[3];
rz(-1.691317) q[3];
sx q[3];
rz(-2.4108346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47736436) q[0];
sx q[0];
rz(-2.1381162) q[0];
sx q[0];
rz(-0.088951237) q[0];
rz(-1.7695919) q[1];
sx q[1];
rz(-1.7141432) q[1];
sx q[1];
rz(-2.037183) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7331284) q[0];
sx q[0];
rz(-2.1845572) q[0];
sx q[0];
rz(-1.2888786) q[0];
rz(2.6501765) q[2];
sx q[2];
rz(-0.84553908) q[2];
sx q[2];
rz(0.27057901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1762878) q[1];
sx q[1];
rz(-1.4215648) q[1];
sx q[1];
rz(1.1231827) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4711963) q[3];
sx q[3];
rz(-1.0336813) q[3];
sx q[3];
rz(-0.84505862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.264512) q[2];
sx q[2];
rz(-1.1254213) q[2];
sx q[2];
rz(1.6748927) q[2];
rz(-3.1125715) q[3];
sx q[3];
rz(-2.1252188) q[3];
sx q[3];
rz(0.69028729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90917176) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(-0.27798852) q[0];
rz(1.7104644) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(2.7412282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5059197) q[0];
sx q[0];
rz(-0.95209661) q[0];
sx q[0];
rz(2.1257504) q[0];
rz(-pi) q[1];
rz(1.8258926) q[2];
sx q[2];
rz(-0.9005024) q[2];
sx q[2];
rz(2.3885661) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4440205) q[1];
sx q[1];
rz(-2.6409147) q[1];
sx q[1];
rz(-0.33659192) q[1];
rz(-2.4147846) q[3];
sx q[3];
rz(-0.44379674) q[3];
sx q[3];
rz(-0.82647317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4572738) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(-0.45480967) q[2];
rz(1.8999892) q[3];
sx q[3];
rz(-0.95507115) q[3];
sx q[3];
rz(2.4212867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860745) q[0];
sx q[0];
rz(-1.8121413) q[0];
sx q[0];
rz(3.1410826) q[0];
rz(-0.60091248) q[1];
sx q[1];
rz(-0.83952236) q[1];
sx q[1];
rz(0.14437637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91286425) q[0];
sx q[0];
rz(-1.0236386) q[0];
sx q[0];
rz(-2.8990101) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52835502) q[2];
sx q[2];
rz(-0.29137416) q[2];
sx q[2];
rz(1.0519415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6295602) q[1];
sx q[1];
rz(-1.343439) q[1];
sx q[1];
rz(-2.4870706) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17474971) q[3];
sx q[3];
rz(-2.5876382) q[3];
sx q[3];
rz(0.21077158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1986177) q[2];
sx q[2];
rz(-1.696442) q[2];
sx q[2];
rz(-0.96251881) q[2];
rz(1.5396384) q[3];
sx q[3];
rz(-1.4055777) q[3];
sx q[3];
rz(-0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9050423) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(1.9566253) q[0];
rz(0.22625893) q[1];
sx q[1];
rz(-2.2626651) q[1];
sx q[1];
rz(2.102898) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0112891) q[0];
sx q[0];
rz(-1.652392) q[0];
sx q[0];
rz(2.5528583) q[0];
rz(0.73815411) q[2];
sx q[2];
rz(-1.0760436) q[2];
sx q[2];
rz(0.02502266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86335582) q[1];
sx q[1];
rz(-2.0788631) q[1];
sx q[1];
rz(-2.3418576) q[1];
rz(-pi) q[2];
rz(2.0108972) q[3];
sx q[3];
rz(-1.6301042) q[3];
sx q[3];
rz(-1.9594994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.431939) q[2];
sx q[2];
rz(-0.69810549) q[2];
sx q[2];
rz(-2.8254438) q[2];
rz(1.6759253) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50452152) q[0];
sx q[0];
rz(-0.9032473) q[0];
sx q[0];
rz(-1.1035408) q[0];
rz(2.6612813) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(-0.59741098) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99973122) q[0];
sx q[0];
rz(-1.7787361) q[0];
sx q[0];
rz(-0.19241649) q[0];
rz(-0.2603724) q[2];
sx q[2];
rz(-2.3790759) q[2];
sx q[2];
rz(-2.0338361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5525517) q[1];
sx q[1];
rz(-1.8218166) q[1];
sx q[1];
rz(1.7415206) q[1];
x q[2];
rz(-1.1998981) q[3];
sx q[3];
rz(-1.4850332) q[3];
sx q[3];
rz(-0.24468064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9395113) q[2];
sx q[2];
rz(-1.6263522) q[2];
sx q[2];
rz(-1.0763947) q[2];
rz(1.1427897) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(2.6514261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48653212) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(-3.1031188) q[0];
rz(-3.0768652) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(-0.23385349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9486987) q[0];
sx q[0];
rz(-1.3854376) q[0];
sx q[0];
rz(-0.22484397) q[0];
rz(-pi) q[1];
rz(1.2042768) q[2];
sx q[2];
rz(-1.8631753) q[2];
sx q[2];
rz(2.2077843) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95215248) q[1];
sx q[1];
rz(-0.92467151) q[1];
sx q[1];
rz(-1.6829856) q[1];
x q[2];
rz(2.7247808) q[3];
sx q[3];
rz(-2.9997065) q[3];
sx q[3];
rz(-0.2896504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48876277) q[2];
sx q[2];
rz(-1.5831999) q[2];
sx q[2];
rz(2.5433507) q[2];
rz(-0.11387842) q[3];
sx q[3];
rz(-1.7555534) q[3];
sx q[3];
rz(2.2940476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4050196) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(-2.8952428) q[0];
rz(-1.2975289) q[1];
sx q[1];
rz(-2.0228701) q[1];
sx q[1];
rz(-0.49682239) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637359) q[0];
sx q[0];
rz(-1.1564768) q[0];
sx q[0];
rz(0.68092771) q[0];
x q[1];
rz(0.17077568) q[2];
sx q[2];
rz(-0.71238067) q[2];
sx q[2];
rz(-2.8757489) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5766746) q[1];
sx q[1];
rz(-1.0905301) q[1];
sx q[1];
rz(2.8757921) q[1];
x q[2];
rz(-1.854135) q[3];
sx q[3];
rz(-2.6322798) q[3];
sx q[3];
rz(0.71374245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7184427) q[2];
sx q[2];
rz(-0.53360525) q[2];
sx q[2];
rz(-1.144484) q[2];
rz(-1.1540958) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(2.9437039) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941187) q[0];
sx q[0];
rz(-2.9070774) q[0];
sx q[0];
rz(-0.06614729) q[0];
rz(1.9937531) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(-2.537421) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2691374) q[0];
sx q[0];
rz(-2.7949484) q[0];
sx q[0];
rz(2.4983062) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2653052) q[2];
sx q[2];
rz(-2.5992706) q[2];
sx q[2];
rz(0.10077439) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0688404) q[1];
sx q[1];
rz(-1.3786331) q[1];
sx q[1];
rz(0.56984624) q[1];
rz(-2.3202032) q[3];
sx q[3];
rz(-0.8474955) q[3];
sx q[3];
rz(2.5625474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8784647) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(-0.43295941) q[2];
rz(1.912502) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(1.8142726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94806725) q[0];
sx q[0];
rz(-2.0663517) q[0];
sx q[0];
rz(-1.2731592) q[0];
rz(-0.46514568) q[1];
sx q[1];
rz(-1.7756614) q[1];
sx q[1];
rz(0.21496162) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0996159) q[0];
sx q[0];
rz(-0.93659217) q[0];
sx q[0];
rz(1.2553822) q[0];
rz(-pi) q[1];
rz(-2.4768171) q[2];
sx q[2];
rz(-2.4996346) q[2];
sx q[2];
rz(2.4017815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7219639) q[1];
sx q[1];
rz(-0.89744324) q[1];
sx q[1];
rz(-0.90760214) q[1];
rz(-pi) q[2];
rz(-0.39542012) q[3];
sx q[3];
rz(-2.5428465) q[3];
sx q[3];
rz(2.0778164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8421858) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(2.1827533) q[2];
rz(-1.1622608) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(-1.4452665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5398298) q[0];
sx q[0];
rz(-1.5400664) q[0];
sx q[0];
rz(-1.6590317) q[0];
rz(-0.70855793) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(0.49025771) q[2];
sx q[2];
rz(-2.5210862) q[2];
sx q[2];
rz(1.6747337) q[2];
rz(1.3832573) q[3];
sx q[3];
rz(-0.85761999) q[3];
sx q[3];
rz(1.7354538) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
