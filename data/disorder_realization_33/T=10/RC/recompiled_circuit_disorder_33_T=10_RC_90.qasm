OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(6.0072748) q[0];
sx q[0];
rz(10.732565) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6558134) q[0];
sx q[0];
rz(-1.2278779) q[0];
sx q[0];
rz(-1.2113843) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70648944) q[2];
sx q[2];
rz(-0.90105614) q[2];
sx q[2];
rz(2.0073839) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5174487) q[1];
sx q[1];
rz(-1.2341208) q[1];
sx q[1];
rz(-1.3144073) q[1];
x q[2];
rz(-2.9458463) q[3];
sx q[3];
rz(-1.9276852) q[3];
sx q[3];
rz(1.9407879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2661665) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(2.0092633) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(2.1291389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9448626) q[0];
sx q[0];
rz(-2.9319627) q[0];
sx q[0];
rz(-2.9557513) q[0];
rz(-2.5813685) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(-0.21683189) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2730392) q[0];
sx q[0];
rz(-0.93398636) q[0];
sx q[0];
rz(-2.7815232) q[0];
rz(-pi) q[1];
rz(-1.0160604) q[2];
sx q[2];
rz(-1.1604571) q[2];
sx q[2];
rz(-1.0085269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20570457) q[1];
sx q[1];
rz(-1.3591896) q[1];
sx q[1];
rz(2.2536224) q[1];
rz(-pi) q[2];
rz(2.5012245) q[3];
sx q[3];
rz(-0.98207563) q[3];
sx q[3];
rz(2.0224188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8314787) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(-1.8537834) q[2];
rz(-0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(2.8365703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-2.537354) q[0];
rz(-1.8151981) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-0.93260971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9661449) q[0];
sx q[0];
rz(-1.6245337) q[0];
sx q[0];
rz(-1.2530112) q[0];
rz(-pi) q[1];
rz(1.7980174) q[2];
sx q[2];
rz(-1.7482687) q[2];
sx q[2];
rz(2.9858659) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0886791) q[1];
sx q[1];
rz(-2.1267849) q[1];
sx q[1];
rz(-1.710379) q[1];
rz(-pi) q[2];
rz(-0.44585769) q[3];
sx q[3];
rz(-2.0816396) q[3];
sx q[3];
rz(-1.0872935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9937667) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(-2.0489342) q[2];
rz(-2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7820691) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(2.6413667) q[0];
rz(0.80530986) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(1.4979699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26101199) q[0];
sx q[0];
rz(-1.2769165) q[0];
sx q[0];
rz(-2.0902993) q[0];
x q[1];
rz(1.3802714) q[2];
sx q[2];
rz(-1.5152144) q[2];
sx q[2];
rz(-1.7855102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8343463) q[1];
sx q[1];
rz(-1.3386968) q[1];
sx q[1];
rz(-2.8744065) q[1];
rz(-0.46521503) q[3];
sx q[3];
rz(-1.1606996) q[3];
sx q[3];
rz(1.0614392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(0.70181075) q[2];
rz(-2.3102405) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9005301) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(0.81533122) q[0];
rz(1.6197846) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(-1.048208) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5532903) q[0];
sx q[0];
rz(-0.36670812) q[0];
sx q[0];
rz(-0.81272965) q[0];
rz(-pi) q[1];
rz(0.65152119) q[2];
sx q[2];
rz(-1.5805575) q[2];
sx q[2];
rz(-2.5196645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2506927) q[1];
sx q[1];
rz(-0.39784583) q[1];
sx q[1];
rz(0.65244168) q[1];
x q[2];
rz(2.1225554) q[3];
sx q[3];
rz(-0.69283797) q[3];
sx q[3];
rz(0.90941959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(2.0416416) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0681756) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(0.90240479) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(-0.12983233) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31045612) q[0];
sx q[0];
rz(-1.2015011) q[0];
sx q[0];
rz(-2.5184758) q[0];
rz(-pi) q[1];
rz(2.1075222) q[2];
sx q[2];
rz(-0.64646361) q[2];
sx q[2];
rz(-1.549364) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1631158) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(-2.5549868) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3744266) q[3];
sx q[3];
rz(-1.0305627) q[3];
sx q[3];
rz(1.2353209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3123902) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(-0.20425805) q[2];
rz(1.2060818) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234574) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(1.4468505) q[0];
rz(1.2591259) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(0.68626219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26353729) q[0];
sx q[0];
rz(-1.1362846) q[0];
sx q[0];
rz(-2.6146019) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5574066) q[2];
sx q[2];
rz(-0.91911941) q[2];
sx q[2];
rz(-2.3551031) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4916617) q[1];
sx q[1];
rz(-1.4964536) q[1];
sx q[1];
rz(1.3851623) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2495743) q[3];
sx q[3];
rz(-1.9739082) q[3];
sx q[3];
rz(1.5954799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4454322) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(2.5799675) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6034265) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(1.6954533) q[0];
rz(2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(1.6400281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7555996) q[0];
sx q[0];
rz(-2.0634683) q[0];
sx q[0];
rz(-0.022105769) q[0];
x q[1];
rz(-0.026560606) q[2];
sx q[2];
rz(-1.5090669) q[2];
sx q[2];
rz(-0.82690566) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6540263) q[1];
sx q[1];
rz(-1.8974202) q[1];
sx q[1];
rz(1.6649654) q[1];
rz(-pi) q[2];
rz(2.5128965) q[3];
sx q[3];
rz(-2.0572212) q[3];
sx q[3];
rz(-1.4592255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35187307) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(-1.8224576) q[2];
rz(1.2119279) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(-2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(0.38326344) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4708913) q[0];
sx q[0];
rz(-2.0235217) q[0];
sx q[0];
rz(0.41505138) q[0];
rz(-pi) q[1];
rz(-0.36231626) q[2];
sx q[2];
rz(-10*pi/13) q[2];
sx q[2];
rz(2.97646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5896776) q[1];
sx q[1];
rz(-2.3304686) q[1];
sx q[1];
rz(-0.07304904) q[1];
rz(1.5185235) q[3];
sx q[3];
rz(-1.6802603) q[3];
sx q[3];
rz(2.0930406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(-1.8593672) q[2];
rz(-1.6451689) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(-2.0857247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4984109) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(-0.19432755) q[0];
rz(2.1037897) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(2.1077572) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72458306) q[0];
sx q[0];
rz(-1.4405182) q[0];
sx q[0];
rz(0.91086046) q[0];
rz(-0.86976544) q[2];
sx q[2];
rz(-1.1681721) q[2];
sx q[2];
rz(2.082777) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49627134) q[1];
sx q[1];
rz(-1.6722582) q[1];
sx q[1];
rz(-2.1437777) q[1];
x q[2];
rz(-0.25063534) q[3];
sx q[3];
rz(-2.4126629) q[3];
sx q[3];
rz(-0.38754101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0620492) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(-0.27030269) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939659) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.7059965) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(-1.5244665) q[2];
sx q[2];
rz(-2.5382858) q[2];
sx q[2];
rz(-0.45679191) q[2];
rz(1.2407606) q[3];
sx q[3];
rz(-1.6374554) q[3];
sx q[3];
rz(-1.0367254) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
