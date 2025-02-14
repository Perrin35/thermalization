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
rz(0.49405721) q[0];
sx q[0];
rz(2.9394975) q[0];
sx q[0];
rz(9.1261368) q[0];
rz(-0.59386295) q[1];
sx q[1];
rz(-1.2358103) q[1];
sx q[1];
rz(1.2716582) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8182871) q[0];
sx q[0];
rz(-3.0491017) q[0];
sx q[0];
rz(0.166786) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3321882) q[2];
sx q[2];
rz(-2.9706741) q[2];
sx q[2];
rz(-1.6910291) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.16954452) q[1];
sx q[1];
rz(-1.3005592) q[1];
sx q[1];
rz(1.0098021) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9416973) q[3];
sx q[3];
rz(-1.8270396) q[3];
sx q[3];
rz(-0.56396644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1539845) q[2];
sx q[2];
rz(-2.5488904) q[2];
sx q[2];
rz(-2.985037) q[2];
rz(2.9708059) q[3];
sx q[3];
rz(-2.4510866) q[3];
sx q[3];
rz(-1.7295711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.60468948) q[0];
sx q[0];
rz(-1.5770183) q[0];
sx q[0];
rz(-0.94737303) q[0];
rz(-1.0700215) q[1];
sx q[1];
rz(-1.0173631) q[1];
sx q[1];
rz(-2.1224461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6603834) q[0];
sx q[0];
rz(-1.0374746) q[0];
sx q[0];
rz(0.18444277) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48430125) q[2];
sx q[2];
rz(-1.2107841) q[2];
sx q[2];
rz(-1.1126435) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6883057) q[1];
sx q[1];
rz(-0.46436342) q[1];
sx q[1];
rz(2.7771441) q[1];
x q[2];
rz(-0.74599539) q[3];
sx q[3];
rz(-0.87019414) q[3];
sx q[3];
rz(-1.7563411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77000874) q[2];
sx q[2];
rz(-1.2086478) q[2];
sx q[2];
rz(-0.51119512) q[2];
rz(-2.8080071) q[3];
sx q[3];
rz(-2.249692) q[3];
sx q[3];
rz(0.79264486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12673512) q[0];
sx q[0];
rz(-2.6124586) q[0];
sx q[0];
rz(-1.0134617) q[0];
rz(-0.86499372) q[1];
sx q[1];
rz(-1.8987055) q[1];
sx q[1];
rz(-1.6046883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44492267) q[0];
sx q[0];
rz(-1.4620263) q[0];
sx q[0];
rz(-0.12420267) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3853677) q[2];
sx q[2];
rz(-0.87487223) q[2];
sx q[2];
rz(2.8967901) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.930248) q[1];
sx q[1];
rz(-2.4635163) q[1];
sx q[1];
rz(-0.67386595) q[1];
rz(0.38648326) q[3];
sx q[3];
rz(-2.1162652) q[3];
sx q[3];
rz(0.89279934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.65261877) q[2];
sx q[2];
rz(-1.3451312) q[2];
sx q[2];
rz(2.6441003) q[2];
rz(0.8688212) q[3];
sx q[3];
rz(-2.7610064) q[3];
sx q[3];
rz(0.07180056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82197613) q[0];
sx q[0];
rz(-0.0030586547) q[0];
sx q[0];
rz(-0.36913607) q[0];
rz(-2.4311851) q[1];
sx q[1];
rz(-2.0901168) q[1];
sx q[1];
rz(-1.2999473) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1064736) q[0];
sx q[0];
rz(-1.4958315) q[0];
sx q[0];
rz(-2.3142561) q[0];
x q[1];
rz(-1.4163912) q[2];
sx q[2];
rz(-2.0090299) q[2];
sx q[2];
rz(-1.1296425) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.877546) q[1];
sx q[1];
rz(-2.6461671) q[1];
sx q[1];
rz(-2.3645414) q[1];
rz(1.348576) q[3];
sx q[3];
rz(-0.45688964) q[3];
sx q[3];
rz(0.27049388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9909782) q[2];
sx q[2];
rz(-1.6852448) q[2];
sx q[2];
rz(2.4197104) q[2];
rz(2.7041096) q[3];
sx q[3];
rz(-0.36967725) q[3];
sx q[3];
rz(0.30445254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7337912) q[0];
sx q[0];
rz(-0.79923874) q[0];
sx q[0];
rz(-2.9275628) q[0];
rz(0.80263823) q[1];
sx q[1];
rz(-1.1624348) q[1];
sx q[1];
rz(-2.7208557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9393552) q[0];
sx q[0];
rz(-1.4927683) q[0];
sx q[0];
rz(1.3952888) q[0];
rz(1.8123295) q[2];
sx q[2];
rz(-1.0189172) q[2];
sx q[2];
rz(-0.59029366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.72374684) q[1];
sx q[1];
rz(-1.6366948) q[1];
sx q[1];
rz(1.4562277) q[1];
rz(-pi) q[2];
rz(-1.7047799) q[3];
sx q[3];
rz(-2.4089536) q[3];
sx q[3];
rz(2.028038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3558827) q[2];
sx q[2];
rz(-0.37675253) q[2];
sx q[2];
rz(-2.9316087) q[2];
rz(1.0979794) q[3];
sx q[3];
rz(-0.72303253) q[3];
sx q[3];
rz(-2.9155904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2881154) q[0];
sx q[0];
rz(-2.2324039) q[0];
sx q[0];
rz(-0.01532456) q[0];
rz(2.3769936) q[1];
sx q[1];
rz(-1.2654283) q[1];
sx q[1];
rz(-2.9008124) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.044996) q[0];
sx q[0];
rz(-1.2701663) q[0];
sx q[0];
rz(-2.0040254) q[0];
rz(-pi) q[1];
rz(0.44220512) q[2];
sx q[2];
rz(-0.54536146) q[2];
sx q[2];
rz(2.9186625) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7635212) q[1];
sx q[1];
rz(-1.3978036) q[1];
sx q[1];
rz(0.066336169) q[1];
rz(2.7731259) q[3];
sx q[3];
rz(-1.0771015) q[3];
sx q[3];
rz(-1.5696978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2601872) q[2];
sx q[2];
rz(-2.0892961) q[2];
sx q[2];
rz(-0.50797272) q[2];
rz(2.9925665) q[3];
sx q[3];
rz(-2.7601354) q[3];
sx q[3];
rz(1.6833359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43450612) q[0];
sx q[0];
rz(-0.91489804) q[0];
sx q[0];
rz(-1.1900505) q[0];
rz(3.1022364) q[1];
sx q[1];
rz(-2.1085565) q[1];
sx q[1];
rz(2.8796223) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16548115) q[0];
sx q[0];
rz(-1.9046749) q[0];
sx q[0];
rz(1.5206485) q[0];
x q[1];
rz(0.33032051) q[2];
sx q[2];
rz(-1.9322504) q[2];
sx q[2];
rz(0.20955958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90967623) q[1];
sx q[1];
rz(-1.2753829) q[1];
sx q[1];
rz(2.7353103) q[1];
x q[2];
rz(1.3218325) q[3];
sx q[3];
rz(-2.1002703) q[3];
sx q[3];
rz(-1.7789237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1956341) q[2];
sx q[2];
rz(-0.28886637) q[2];
sx q[2];
rz(1.757901) q[2];
rz(1.1727928) q[3];
sx q[3];
rz(-1.5044418) q[3];
sx q[3];
rz(2.8891532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4381572) q[0];
sx q[0];
rz(-2.1717838) q[0];
sx q[0];
rz(0.54617149) q[0];
rz(-0.81900412) q[1];
sx q[1];
rz(-1.9592229) q[1];
sx q[1];
rz(-0.46772734) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0696164) q[0];
sx q[0];
rz(-1.8307951) q[0];
sx q[0];
rz(-3.0896356) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7502656) q[2];
sx q[2];
rz(-1.9576009) q[2];
sx q[2];
rz(-2.4387173) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7744171) q[1];
sx q[1];
rz(-2.0655334) q[1];
sx q[1];
rz(-2.6082188) q[1];
rz(-pi) q[2];
rz(-2.9727861) q[3];
sx q[3];
rz(-1.3309817) q[3];
sx q[3];
rz(0.60926357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8972682) q[2];
sx q[2];
rz(-1.6017598) q[2];
sx q[2];
rz(2.3648105) q[2];
rz(1.6952093) q[3];
sx q[3];
rz(-1.8518238) q[3];
sx q[3];
rz(-2.9608534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.986213) q[0];
sx q[0];
rz(-1.7387094) q[0];
sx q[0];
rz(-0.078393161) q[0];
rz(-0.29131237) q[1];
sx q[1];
rz(-2.4261609) q[1];
sx q[1];
rz(-1.327347) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3122897) q[0];
sx q[0];
rz(-0.96594496) q[0];
sx q[0];
rz(0.74440794) q[0];
rz(0.38703309) q[2];
sx q[2];
rz(-3.0245028) q[2];
sx q[2];
rz(2.4486604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9207245) q[1];
sx q[1];
rz(-0.70573843) q[1];
sx q[1];
rz(-2.7138813) q[1];
rz(2.6895608) q[3];
sx q[3];
rz(-2.8252183) q[3];
sx q[3];
rz(-1.3133501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.766481) q[2];
sx q[2];
rz(-1.8152619) q[2];
sx q[2];
rz(0.43409902) q[2];
rz(0.33657524) q[3];
sx q[3];
rz(-1.7124636) q[3];
sx q[3];
rz(2.7460639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3340988) q[0];
sx q[0];
rz(-2.8039126) q[0];
sx q[0];
rz(-0.39921528) q[0];
rz(0.68398625) q[1];
sx q[1];
rz(-1.8268879) q[1];
sx q[1];
rz(-1.7302053) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1841662) q[0];
sx q[0];
rz(-1.7475351) q[0];
sx q[0];
rz(0.33542963) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67839854) q[2];
sx q[2];
rz(-1.7400555) q[2];
sx q[2];
rz(2.9321743) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.70727247) q[1];
sx q[1];
rz(-2.1967683) q[1];
sx q[1];
rz(-2.0213338) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3606922) q[3];
sx q[3];
rz(-2.0640496) q[3];
sx q[3];
rz(0.79003143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6902249) q[2];
sx q[2];
rz(-1.6121612) q[2];
sx q[2];
rz(2.8239047) q[2];
rz(-2.132527) q[3];
sx q[3];
rz(-1.7694446) q[3];
sx q[3];
rz(-2.7528929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4036564) q[0];
sx q[0];
rz(-1.6725578) q[0];
sx q[0];
rz(-2.9148711) q[0];
rz(-2.8059088) q[1];
sx q[1];
rz(-0.93135584) q[1];
sx q[1];
rz(2.9109536) q[1];
rz(2.4776501) q[2];
sx q[2];
rz(-1.5879122) q[2];
sx q[2];
rz(1.5664311) q[2];
rz(-0.54343358) q[3];
sx q[3];
rz(-2.4590878) q[3];
sx q[3];
rz(1.2868846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
