OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.090341181) q[0];
sx q[0];
rz(-1.0426961) q[0];
sx q[0];
rz(-2.8629177) q[0];
rz(-1.9947808) q[1];
sx q[1];
rz(-2.6169701) q[1];
sx q[1];
rz(1.2710748) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81084736) q[0];
sx q[0];
rz(-0.26118229) q[0];
sx q[0];
rz(-3.0708341) q[0];
rz(2.4248883) q[2];
sx q[2];
rz(-0.19858805) q[2];
sx q[2];
rz(-2.8633397) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6587352) q[1];
sx q[1];
rz(-2.4244011) q[1];
sx q[1];
rz(-1.6916005) q[1];
x q[2];
rz(-0.99156143) q[3];
sx q[3];
rz(-2.09822) q[3];
sx q[3];
rz(0.19471951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.110454) q[2];
sx q[2];
rz(-1.217507) q[2];
sx q[2];
rz(-2.8793867) q[2];
rz(-2.0860705) q[3];
sx q[3];
rz(-1.2473829) q[3];
sx q[3];
rz(3.054255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-1.5259041) q[0];
sx q[0];
rz(-0.50351244) q[0];
sx q[0];
rz(0.18158922) q[0];
rz(-1.4008105) q[1];
sx q[1];
rz(-1.9488275) q[1];
sx q[1];
rz(-0.84091944) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5530295) q[0];
sx q[0];
rz(-0.32051495) q[0];
sx q[0];
rz(-3.0856087) q[0];
x q[1];
rz(2.1757751) q[2];
sx q[2];
rz(-2.4060095) q[2];
sx q[2];
rz(-0.090026131) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5121756) q[1];
sx q[1];
rz(-1.4025153) q[1];
sx q[1];
rz(0.03668935) q[1];
rz(-pi) q[2];
rz(0.36786973) q[3];
sx q[3];
rz(-2.1506393) q[3];
sx q[3];
rz(0.90880064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.759364) q[2];
sx q[2];
rz(-0.21491924) q[2];
sx q[2];
rz(-2.0767029) q[2];
rz(3.0801638) q[3];
sx q[3];
rz(-1.8062785) q[3];
sx q[3];
rz(-1.5016865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8432124) q[0];
sx q[0];
rz(-1.9623663) q[0];
sx q[0];
rz(-2.9425353) q[0];
rz(2.3152323) q[1];
sx q[1];
rz(-2.2233456) q[1];
sx q[1];
rz(-2.3645649) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13764405) q[0];
sx q[0];
rz(-0.17034082) q[0];
sx q[0];
rz(1.9714703) q[0];
x q[1];
rz(-0.42962809) q[2];
sx q[2];
rz(-1.0681392) q[2];
sx q[2];
rz(-0.70870295) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2232845) q[1];
sx q[1];
rz(-1.3121843) q[1];
sx q[1];
rz(-1.0041987) q[1];
rz(2.8641508) q[3];
sx q[3];
rz(-2.146477) q[3];
sx q[3];
rz(-2.4572008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0394773) q[2];
sx q[2];
rz(-1.1481043) q[2];
sx q[2];
rz(-0.51304212) q[2];
rz(-0.78220621) q[3];
sx q[3];
rz(-1.2057883) q[3];
sx q[3];
rz(-1.7647083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9272598) q[0];
sx q[0];
rz(-1.5434649) q[0];
sx q[0];
rz(1.6511035) q[0];
rz(0.53228846) q[1];
sx q[1];
rz(-0.63096255) q[1];
sx q[1];
rz(1.5375563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69025264) q[0];
sx q[0];
rz(-1.4079388) q[0];
sx q[0];
rz(-1.4630408) q[0];
rz(0.093393383) q[2];
sx q[2];
rz(-1.0265145) q[2];
sx q[2];
rz(1.5285847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.93354952) q[1];
sx q[1];
rz(-1.4999823) q[1];
sx q[1];
rz(-2.2068702) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7521637) q[3];
sx q[3];
rz(-2.5797322) q[3];
sx q[3];
rz(-1.495468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6060467) q[2];
sx q[2];
rz(-2.594785) q[2];
sx q[2];
rz(3.1287076) q[2];
rz(1.9384725) q[3];
sx q[3];
rz(-1.4667526) q[3];
sx q[3];
rz(-2.119868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1383698) q[0];
sx q[0];
rz(-0.62507451) q[0];
sx q[0];
rz(1.7417997) q[0];
rz(-2.7895582) q[1];
sx q[1];
rz(-2.0875918) q[1];
sx q[1];
rz(-2.3037516) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1803446) q[0];
sx q[0];
rz(-0.77416468) q[0];
sx q[0];
rz(-1.3526239) q[0];
rz(-pi) q[1];
rz(0.61428563) q[2];
sx q[2];
rz(-1.0453684) q[2];
sx q[2];
rz(-3.0802259) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.067474631) q[1];
sx q[1];
rz(-0.58093853) q[1];
sx q[1];
rz(0.23260637) q[1];
x q[2];
rz(1.4170333) q[3];
sx q[3];
rz(-1.7332675) q[3];
sx q[3];
rz(-2.3547309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1167404) q[2];
sx q[2];
rz(-0.96472538) q[2];
sx q[2];
rz(-1.1836729) q[2];
rz(-2.8708894) q[3];
sx q[3];
rz(-2.201122) q[3];
sx q[3];
rz(-1.6089599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1149513) q[0];
sx q[0];
rz(-2.5178435) q[0];
sx q[0];
rz(0.57914105) q[0];
rz(-1.2222414) q[1];
sx q[1];
rz(-1.8341583) q[1];
sx q[1];
rz(-1.1096257) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0298808) q[0];
sx q[0];
rz(-1.7222026) q[0];
sx q[0];
rz(-2.9747308) q[0];
rz(-pi) q[1];
rz(-1.6272306) q[2];
sx q[2];
rz(-1.2662958) q[2];
sx q[2];
rz(-2.3431108) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5612003) q[1];
sx q[1];
rz(-1.322974) q[1];
sx q[1];
rz(0.3419673) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6251276) q[3];
sx q[3];
rz(-2.0544009) q[3];
sx q[3];
rz(2.2276239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9455202) q[2];
sx q[2];
rz(-1.2976982) q[2];
sx q[2];
rz(0.35745364) q[2];
rz(-1.9516021) q[3];
sx q[3];
rz(-1.9414732) q[3];
sx q[3];
rz(1.3397217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5424407) q[0];
sx q[0];
rz(-1.0170794) q[0];
sx q[0];
rz(2.3495112) q[0];
rz(0.7041086) q[1];
sx q[1];
rz(-2.6857565) q[1];
sx q[1];
rz(-1.5584996) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7575551) q[0];
sx q[0];
rz(-0.4641986) q[0];
sx q[0];
rz(-0.0065825855) q[0];
rz(-pi) q[1];
rz(2.3249945) q[2];
sx q[2];
rz(-2.3629521) q[2];
sx q[2];
rz(-1.9684079) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.577282) q[1];
sx q[1];
rz(-2.7140712) q[1];
sx q[1];
rz(3.0583818) q[1];
x q[2];
rz(0.24105234) q[3];
sx q[3];
rz(-0.76605443) q[3];
sx q[3];
rz(1.1055921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11129561) q[2];
sx q[2];
rz(-0.78812495) q[2];
sx q[2];
rz(-3.1374068) q[2];
rz(-2.8748416) q[3];
sx q[3];
rz(-1.5240074) q[3];
sx q[3];
rz(-0.72223103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58077145) q[0];
sx q[0];
rz(-2.4120791) q[0];
sx q[0];
rz(-0.15604493) q[0];
rz(1.5566298) q[1];
sx q[1];
rz(-2.2792351) q[1];
sx q[1];
rz(2.5468266) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0069591) q[0];
sx q[0];
rz(-1.3158321) q[0];
sx q[0];
rz(-2.6322271) q[0];
x q[1];
rz(-0.34551106) q[2];
sx q[2];
rz(-2.4191769) q[2];
sx q[2];
rz(-1.9986716) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3799499) q[1];
sx q[1];
rz(-2.159286) q[1];
sx q[1];
rz(2.2698927) q[1];
rz(-0.90601633) q[3];
sx q[3];
rz(-0.95905868) q[3];
sx q[3];
rz(-2.8869528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5742089) q[2];
sx q[2];
rz(-1.5625861) q[2];
sx q[2];
rz(-2.1136005) q[2];
rz(1.7993641) q[3];
sx q[3];
rz(-0.65516156) q[3];
sx q[3];
rz(-1.8067693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26355711) q[0];
sx q[0];
rz(-1.6850543) q[0];
sx q[0];
rz(-0.79826075) q[0];
rz(-0.80052605) q[1];
sx q[1];
rz(-2.6526178) q[1];
sx q[1];
rz(2.5953603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9437708) q[0];
sx q[0];
rz(-1.4620302) q[0];
sx q[0];
rz(2.1740211) q[0];
x q[1];
rz(2.5419919) q[2];
sx q[2];
rz(-2.4614459) q[2];
sx q[2];
rz(1.6161275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7475976) q[1];
sx q[1];
rz(-2.0556652) q[1];
sx q[1];
rz(-2.5970244) q[1];
x q[2];
rz(0.30719884) q[3];
sx q[3];
rz(-2.9044378) q[3];
sx q[3];
rz(-3.0320771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3937248) q[2];
sx q[2];
rz(-1.0900898) q[2];
sx q[2];
rz(3.0785576) q[2];
rz(-0.6978327) q[3];
sx q[3];
rz(-0.2581667) q[3];
sx q[3];
rz(0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49552712) q[0];
sx q[0];
rz(-0.87764469) q[0];
sx q[0];
rz(2.7255507) q[0];
rz(-1.7795732) q[1];
sx q[1];
rz(-1.1241309) q[1];
sx q[1];
rz(1.2927885) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36609367) q[0];
sx q[0];
rz(-1.3028569) q[0];
sx q[0];
rz(1.3331205) q[0];
x q[1];
rz(-2.6854158) q[2];
sx q[2];
rz(-1.8229228) q[2];
sx q[2];
rz(1.1751564) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.67329183) q[1];
sx q[1];
rz(-2.3929962) q[1];
sx q[1];
rz(0.8573726) q[1];
rz(-pi) q[2];
rz(-1.9927916) q[3];
sx q[3];
rz(-1.3316324) q[3];
sx q[3];
rz(-1.07475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4585939) q[2];
sx q[2];
rz(-1.3364044) q[2];
sx q[2];
rz(2.2340753) q[2];
rz(1.8806184) q[3];
sx q[3];
rz(-0.57456273) q[3];
sx q[3];
rz(-1.6921836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4954729) q[0];
sx q[0];
rz(-1.6874122) q[0];
sx q[0];
rz(0.73943403) q[0];
rz(2.6932035) q[1];
sx q[1];
rz(-1.3403475) q[1];
sx q[1];
rz(1.1833804) q[1];
rz(0.35198718) q[2];
sx q[2];
rz(-1.2766196) q[2];
sx q[2];
rz(2.5310493) q[2];
rz(-1.7200593) q[3];
sx q[3];
rz(-2.4076318) q[3];
sx q[3];
rz(-2.9416495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
