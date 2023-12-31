OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(5.3263721) q[0];
sx q[0];
rz(10.858067) q[0];
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(1.1448316) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24554907) q[0];
sx q[0];
rz(-2.890812) q[0];
sx q[0];
rz(1.301469) q[0];
x q[1];
rz(-1.990591) q[2];
sx q[2];
rz(-1.8325873) q[2];
sx q[2];
rz(-0.37444886) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5073587) q[1];
sx q[1];
rz(-0.45295742) q[1];
sx q[1];
rz(1.5276315) q[1];
rz(-pi) q[2];
rz(0.54361312) q[3];
sx q[3];
rz(-1.8279148) q[3];
sx q[3];
rz(-1.6888113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(-1.9072745) q[2];
rz(2.0862789) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-0.82988513) q[0];
rz(2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-2.2944962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0711489) q[0];
sx q[0];
rz(-2.0111472) q[0];
sx q[0];
rz(-1.1308934) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29909889) q[2];
sx q[2];
rz(-2.3397589) q[2];
sx q[2];
rz(2.9849844) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9287195) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(-2.5961155) q[1];
x q[2];
rz(-2.794572) q[3];
sx q[3];
rz(-1.1565398) q[3];
sx q[3];
rz(-1.8300717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0236686) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(-2.4222597) q[2];
rz(1.8524648) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2704724) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(0.57139325) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(2.6142696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29887154) q[0];
sx q[0];
rz(-1.5776331) q[0];
sx q[0];
rz(-2.8323035) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6730404) q[2];
sx q[2];
rz(-1.3125784) q[2];
sx q[2];
rz(0.25707993) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5273683) q[1];
sx q[1];
rz(-0.80263153) q[1];
sx q[1];
rz(-2.7781092) q[1];
rz(-pi) q[2];
rz(1.2445883) q[3];
sx q[3];
rz(-0.86498125) q[3];
sx q[3];
rz(-2.0446442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8335235) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(2.4148338) q[2];
rz(2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(-1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(1.0035275) q[0];
rz(0.040680496) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(-0.8262659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59506455) q[0];
sx q[0];
rz(-1.0142769) q[0];
sx q[0];
rz(-0.80842774) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.063886558) q[2];
sx q[2];
rz(-2.4404844) q[2];
sx q[2];
rz(-2.137616) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1288209) q[1];
sx q[1];
rz(-1.9073309) q[1];
sx q[1];
rz(-3.0613042) q[1];
x q[2];
rz(-2.5097333) q[3];
sx q[3];
rz(-1.8903036) q[3];
sx q[3];
rz(-2.7416122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(-2.2272002) q[2];
rz(3.0363723) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(-0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2545664) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(2.0078833) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(2.5240135) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3591946) q[0];
sx q[0];
rz(-0.054464666) q[0];
sx q[0];
rz(-0.99476238) q[0];
rz(-2.0134301) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(-2.9885459) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4434112) q[1];
sx q[1];
rz(-0.65084208) q[1];
sx q[1];
rz(-1.0990259) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0693552) q[3];
sx q[3];
rz(-2.9655955) q[3];
sx q[3];
rz(-2.3776059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79409838) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(2.9079672) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96930209) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(-2.2604997) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(2.9615013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0463379) q[0];
sx q[0];
rz(-1.006608) q[0];
sx q[0];
rz(2.6020223) q[0];
rz(-pi) q[1];
rz(2.1744556) q[2];
sx q[2];
rz(-1.7269616) q[2];
sx q[2];
rz(-2.9784163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5345726) q[1];
sx q[1];
rz(-1.8587451) q[1];
sx q[1];
rz(1.401591) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0635707) q[3];
sx q[3];
rz(-1.587095) q[3];
sx q[3];
rz(-1.8340045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2580516) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(1.9656666) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(-1.9765967) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(-0.040239008) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9366074) q[0];
sx q[0];
rz(-1.6523223) q[0];
sx q[0];
rz(-1.8100912) q[0];
rz(1.6778498) q[2];
sx q[2];
rz(-2.7333626) q[2];
sx q[2];
rz(1.5232616) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8602627) q[1];
sx q[1];
rz(-1.4308235) q[1];
sx q[1];
rz(0.94957385) q[1];
rz(-pi) q[2];
rz(0.76544806) q[3];
sx q[3];
rz(-2.7141889) q[3];
sx q[3];
rz(0.039507341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0030901) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(0.17383943) q[2];
rz(1.7447757) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(-0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.7370976) q[0];
rz(1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.6361902) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36639402) q[0];
sx q[0];
rz(-1.1530071) q[0];
sx q[0];
rz(1.226107) q[0];
x q[1];
rz(-2.3178188) q[2];
sx q[2];
rz(-2.0588377) q[2];
sx q[2];
rz(1.4866231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1464403) q[1];
sx q[1];
rz(-2.4447828) q[1];
sx q[1];
rz(-0.050683024) q[1];
x q[2];
rz(2.1721341) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(-0.45575842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3796842) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(-0.15052477) q[2];
rz(-1.6020417) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-2.495893) q[0];
rz(-0.49939108) q[1];
sx q[1];
rz(-1.7130518) q[1];
sx q[1];
rz(1.2649149) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7495959) q[0];
sx q[0];
rz(-0.8526593) q[0];
sx q[0];
rz(2.0977661) q[0];
rz(0.55391295) q[2];
sx q[2];
rz(-1.1914807) q[2];
sx q[2];
rz(-1.6162789) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7449194) q[1];
sx q[1];
rz(-1.5349689) q[1];
sx q[1];
rz(-0.14573914) q[1];
rz(-pi) q[2];
rz(-1.3690788) q[3];
sx q[3];
rz(-1.1551876) q[3];
sx q[3];
rz(-2.6430074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.634793) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(2.6169422) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6181347) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(-0.24542228) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(-1.6773178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9058162) q[0];
sx q[0];
rz(-2.896744) q[0];
sx q[0];
rz(0.4723627) q[0];
rz(0.98430888) q[2];
sx q[2];
rz(-1.7585635) q[2];
sx q[2];
rz(0.078660065) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.163584) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(2.7917557) q[1];
rz(-pi) q[2];
rz(-1.1226095) q[3];
sx q[3];
rz(-2.9467694) q[3];
sx q[3];
rz(2.4217055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-2.396092) q[2];
rz(1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-1.8650919) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(-1.8854234) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(-0.55143572) q[2];
sx q[2];
rz(-0.79179344) q[2];
sx q[2];
rz(1.0168016) q[2];
rz(-1.9477378) q[3];
sx q[3];
rz(-2.1115163) q[3];
sx q[3];
rz(1.9797309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
