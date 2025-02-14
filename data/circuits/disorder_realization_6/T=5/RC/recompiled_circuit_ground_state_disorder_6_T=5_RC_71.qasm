OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6239983) q[0];
sx q[0];
rz(2.6156293) q[0];
sx q[0];
rz(9.2064657) q[0];
rz(1.4766308) q[1];
sx q[1];
rz(-2.6638439) q[1];
sx q[1];
rz(-0.55396095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4628578) q[0];
sx q[0];
rz(-2.4241872) q[0];
sx q[0];
rz(-0.82490246) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5289107) q[2];
sx q[2];
rz(-1.6220835) q[2];
sx q[2];
rz(-1.9330658) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.40383717) q[1];
sx q[1];
rz(-0.74450508) q[1];
sx q[1];
rz(1.2584524) q[1];
rz(-pi) q[2];
rz(0.64726909) q[3];
sx q[3];
rz(-2.4993901) q[3];
sx q[3];
rz(2.6489779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7654984) q[2];
sx q[2];
rz(-0.29355294) q[2];
sx q[2];
rz(-1.8739088) q[2];
rz(-0.82967657) q[3];
sx q[3];
rz(-1.5209578) q[3];
sx q[3];
rz(-0.42384306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053452881) q[0];
sx q[0];
rz(-1.8064878) q[0];
sx q[0];
rz(1.394519) q[0];
rz(2.246619) q[1];
sx q[1];
rz(-1.0286237) q[1];
sx q[1];
rz(-1.5590394) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5916067) q[0];
sx q[0];
rz(-2.215452) q[0];
sx q[0];
rz(-2.2695973) q[0];
rz(-pi) q[1];
rz(-0.87662794) q[2];
sx q[2];
rz(-2.6149984) q[2];
sx q[2];
rz(1.9283629) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6239418) q[1];
sx q[1];
rz(-1.2576767) q[1];
sx q[1];
rz(-2.4409358) q[1];
x q[2];
rz(-1.1432173) q[3];
sx q[3];
rz(-1.2914011) q[3];
sx q[3];
rz(-2.2277149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6643657) q[2];
sx q[2];
rz(-1.5218647) q[2];
sx q[2];
rz(-0.030755432) q[2];
rz(2.6886046) q[3];
sx q[3];
rz(-2.9121297) q[3];
sx q[3];
rz(-0.17791137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7396616) q[0];
sx q[0];
rz(-3.0406096) q[0];
sx q[0];
rz(2.3133551) q[0];
rz(3.0939057) q[1];
sx q[1];
rz(-2.2779155) q[1];
sx q[1];
rz(-1.2275068) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9102893) q[0];
sx q[0];
rz(-2.4132015) q[0];
sx q[0];
rz(-0.69407082) q[0];
x q[1];
rz(-1.0989855) q[2];
sx q[2];
rz(-2.6919439) q[2];
sx q[2];
rz(-1.4020021) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3674254) q[1];
sx q[1];
rz(-1.3102505) q[1];
sx q[1];
rz(-2.0266286) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0892598) q[3];
sx q[3];
rz(-2.298758) q[3];
sx q[3];
rz(-0.46848255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0487655) q[2];
sx q[2];
rz(-0.91492492) q[2];
sx q[2];
rz(-1.1009334) q[2];
rz(-0.01384211) q[3];
sx q[3];
rz(-1.3661386) q[3];
sx q[3];
rz(2.3069416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5562627) q[0];
sx q[0];
rz(-1.4034554) q[0];
sx q[0];
rz(-0.334326) q[0];
rz(0.73257929) q[1];
sx q[1];
rz(-0.94894797) q[1];
sx q[1];
rz(-0.11071959) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3601013) q[0];
sx q[0];
rz(-1.5859787) q[0];
sx q[0];
rz(0.76155565) q[0];
rz(1.6640856) q[2];
sx q[2];
rz(-1.7356725) q[2];
sx q[2];
rz(0.071135757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9690937) q[1];
sx q[1];
rz(-2.2374472) q[1];
sx q[1];
rz(0.73760017) q[1];
rz(-pi) q[2];
rz(-0.7752876) q[3];
sx q[3];
rz(-1.9248157) q[3];
sx q[3];
rz(0.86590761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.683814) q[2];
sx q[2];
rz(-0.28202287) q[2];
sx q[2];
rz(1.7337743) q[2];
rz(-1.1553361) q[3];
sx q[3];
rz(-1.9852394) q[3];
sx q[3];
rz(2.9467764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4487576) q[0];
sx q[0];
rz(-2.223707) q[0];
sx q[0];
rz(0.99739972) q[0];
rz(-1.6150486) q[1];
sx q[1];
rz(-0.63957447) q[1];
sx q[1];
rz(-1.4195199) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2364295) q[0];
sx q[0];
rz(-1.7337203) q[0];
sx q[0];
rz(1.7927367) q[0];
x q[1];
rz(-2.3676374) q[2];
sx q[2];
rz(-1.6822691) q[2];
sx q[2];
rz(-0.74711266) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7035326) q[1];
sx q[1];
rz(-2.4171481) q[1];
sx q[1];
rz(2.0634335) q[1];
x q[2];
rz(-0.12357281) q[3];
sx q[3];
rz(-1.6032748) q[3];
sx q[3];
rz(-0.44636727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.851696) q[2];
sx q[2];
rz(-0.09446129) q[2];
sx q[2];
rz(0.32290253) q[2];
rz(2.0276535) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(2.7285301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36393976) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(-2.3349578) q[0];
rz(0.58397645) q[1];
sx q[1];
rz(-1.1176502) q[1];
sx q[1];
rz(-1.9516099) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7936913) q[0];
sx q[0];
rz(-2.7712203) q[0];
sx q[0];
rz(1.9166975) q[0];
rz(-1.3596542) q[2];
sx q[2];
rz(-1.0462772) q[2];
sx q[2];
rz(0.013583029) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.62184238) q[1];
sx q[1];
rz(-1.6222266) q[1];
sx q[1];
rz(1.7396881) q[1];
x q[2];
rz(-0.1889008) q[3];
sx q[3];
rz(-2.402225) q[3];
sx q[3];
rz(0.96773558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7334062) q[2];
sx q[2];
rz(-1.5032282) q[2];
sx q[2];
rz(0.1850941) q[2];
rz(1.5271651) q[3];
sx q[3];
rz(-1.7388758) q[3];
sx q[3];
rz(2.4534498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1899034) q[0];
sx q[0];
rz(-2.1898495) q[0];
sx q[0];
rz(-0.21251799) q[0];
rz(-2.3566133) q[1];
sx q[1];
rz(-1.4947944) q[1];
sx q[1];
rz(2.3627538) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7424282) q[0];
sx q[0];
rz(-1.0013097) q[0];
sx q[0];
rz(-0.67230255) q[0];
rz(-pi) q[1];
rz(-1.1981702) q[2];
sx q[2];
rz(-1.0154795) q[2];
sx q[2];
rz(-0.057387847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1800715) q[1];
sx q[1];
rz(-2.8043282) q[1];
sx q[1];
rz(2.2757761) q[1];
x q[2];
rz(2.2810197) q[3];
sx q[3];
rz(-1.9509776) q[3];
sx q[3];
rz(3.1050888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6276041) q[2];
sx q[2];
rz(-0.49391654) q[2];
sx q[2];
rz(-2.7331875) q[2];
rz(2.4397395) q[3];
sx q[3];
rz(-2.2097094) q[3];
sx q[3];
rz(-1.2592038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34847611) q[0];
sx q[0];
rz(-1.9261253) q[0];
sx q[0];
rz(0.066019639) q[0];
rz(-1.5090212) q[1];
sx q[1];
rz(-1.0450109) q[1];
sx q[1];
rz(-0.95796934) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1692266) q[0];
sx q[0];
rz(-1.1145076) q[0];
sx q[0];
rz(1.8877939) q[0];
x q[1];
rz(-1.0220549) q[2];
sx q[2];
rz(-0.83344747) q[2];
sx q[2];
rz(1.1481783) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0762679) q[1];
sx q[1];
rz(-2.3889414) q[1];
sx q[1];
rz(0.46588834) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1470966) q[3];
sx q[3];
rz(-0.18261431) q[3];
sx q[3];
rz(1.2507358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3079188) q[2];
sx q[2];
rz(-1.6180399) q[2];
sx q[2];
rz(-1.4109122) q[2];
rz(-2.446567) q[3];
sx q[3];
rz(-1.4796939) q[3];
sx q[3];
rz(-0.01785774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0787635) q[0];
sx q[0];
rz(-1.6595027) q[0];
sx q[0];
rz(-0.98989809) q[0];
rz(-0.46317378) q[1];
sx q[1];
rz(-1.6894692) q[1];
sx q[1];
rz(-1.90082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2808229) q[0];
sx q[0];
rz(-1.6284571) q[0];
sx q[0];
rz(-1.2755544) q[0];
x q[1];
rz(-1.1282721) q[2];
sx q[2];
rz(-1.754369) q[2];
sx q[2];
rz(1.5300446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5338143) q[1];
sx q[1];
rz(-1.3097714) q[1];
sx q[1];
rz(-2.8487474) q[1];
rz(-2.1677446) q[3];
sx q[3];
rz(-2.429109) q[3];
sx q[3];
rz(-1.9402011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8076294) q[2];
sx q[2];
rz(-1.7619851) q[2];
sx q[2];
rz(-0.71869746) q[2];
rz(-1.1527609) q[3];
sx q[3];
rz(-1.4216239) q[3];
sx q[3];
rz(1.0144455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5929247) q[0];
sx q[0];
rz(-0.34807006) q[0];
sx q[0];
rz(-2.3396709) q[0];
rz(-2.0536664) q[1];
sx q[1];
rz(-1.3366924) q[1];
sx q[1];
rz(2.4235639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.021026) q[0];
sx q[0];
rz(-0.88827288) q[0];
sx q[0];
rz(-1.1521856) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0775349) q[2];
sx q[2];
rz(-2.0414957) q[2];
sx q[2];
rz(0.59322178) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5637569) q[1];
sx q[1];
rz(-0.33591649) q[1];
sx q[1];
rz(1.4116686) q[1];
rz(-pi) q[2];
rz(-2.2949785) q[3];
sx q[3];
rz(-0.9808971) q[3];
sx q[3];
rz(-2.5835832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3512909) q[2];
sx q[2];
rz(-0.21778926) q[2];
sx q[2];
rz(-2.6289319) q[2];
rz(-0.74448186) q[3];
sx q[3];
rz(-1.0267886) q[3];
sx q[3];
rz(0.7640394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9164593) q[0];
sx q[0];
rz(-1.8712578) q[0];
sx q[0];
rz(1.901392) q[0];
rz(2.4254639) q[1];
sx q[1];
rz(-0.54812535) q[1];
sx q[1];
rz(-2.5352238) q[1];
rz(1.7960346) q[2];
sx q[2];
rz(-1.4627395) q[2];
sx q[2];
rz(1.1278056) q[2];
rz(-1.8567139) q[3];
sx q[3];
rz(-0.34964041) q[3];
sx q[3];
rz(1.6007363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
