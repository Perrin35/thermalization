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
rz(3.1004768) q[0];
sx q[0];
rz(-1.2007204) q[0];
sx q[0];
rz(-1.6706985) q[0];
rz(2.7566551) q[1];
sx q[1];
rz(-2.5632783) q[1];
sx q[1];
rz(-2.2338423) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0931514) q[0];
sx q[0];
rz(-1.0142097) q[0];
sx q[0];
rz(1.0445892) q[0];
rz(-pi) q[1];
rz(-0.29405354) q[2];
sx q[2];
rz(-1.283125) q[2];
sx q[2];
rz(1.2534382) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.585184) q[1];
sx q[1];
rz(-2.6165462) q[1];
sx q[1];
rz(0.37918143) q[1];
rz(-pi) q[2];
rz(-1.5159881) q[3];
sx q[3];
rz(-1.6920838) q[3];
sx q[3];
rz(0.19586043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.82751983) q[2];
sx q[2];
rz(-1.6562409) q[2];
sx q[2];
rz(-0.32076389) q[2];
rz(3.1406506) q[3];
sx q[3];
rz(-2.2195246) q[3];
sx q[3];
rz(0.48085406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8568521) q[0];
sx q[0];
rz(-2.4212403) q[0];
sx q[0];
rz(2.5208933) q[0];
rz(1.2868098) q[1];
sx q[1];
rz(-1.4871253) q[1];
sx q[1];
rz(0.99498814) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5704203) q[0];
sx q[0];
rz(-1.3745527) q[0];
sx q[0];
rz(-1.061383) q[0];
rz(-2.0109977) q[2];
sx q[2];
rz(-1.2792247) q[2];
sx q[2];
rz(-2.7879813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8825922) q[1];
sx q[1];
rz(-2.0152806) q[1];
sx q[1];
rz(0.20436546) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0424588) q[3];
sx q[3];
rz(-0.83513672) q[3];
sx q[3];
rz(2.5723815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93028012) q[2];
sx q[2];
rz(-1.899753) q[2];
sx q[2];
rz(2.9105183) q[2];
rz(1.5288345) q[3];
sx q[3];
rz(-1.684618) q[3];
sx q[3];
rz(-0.55135623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424161) q[0];
sx q[0];
rz(-1.9743974) q[0];
sx q[0];
rz(-2.7447847) q[0];
rz(-1.9467758) q[1];
sx q[1];
rz(-1.8794941) q[1];
sx q[1];
rz(1.6669115) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68434381) q[0];
sx q[0];
rz(-1.5822344) q[0];
sx q[0];
rz(1.4829163) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54946396) q[2];
sx q[2];
rz(-1.1389175) q[2];
sx q[2];
rz(-2.9058676) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.296862) q[1];
sx q[1];
rz(-2.3760188) q[1];
sx q[1];
rz(2.0716448) q[1];
x q[2];
rz(-1.9649404) q[3];
sx q[3];
rz(-1.5817465) q[3];
sx q[3];
rz(-1.2436858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2954703) q[2];
sx q[2];
rz(-0.84540558) q[2];
sx q[2];
rz(-1.6579312) q[2];
rz(0.96019238) q[3];
sx q[3];
rz(-0.63513297) q[3];
sx q[3];
rz(-0.6944164) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33261499) q[0];
sx q[0];
rz(-1.3354477) q[0];
sx q[0];
rz(2.4131925) q[0];
rz(-1.2737466) q[1];
sx q[1];
rz(-1.1796917) q[1];
sx q[1];
rz(-1.5625928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.584986) q[0];
sx q[0];
rz(-1.5674563) q[0];
sx q[0];
rz(-1.5482315) q[0];
x q[1];
rz(3.098084) q[2];
sx q[2];
rz(-1.0232506) q[2];
sx q[2];
rz(-2.7491153) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9409318) q[1];
sx q[1];
rz(-0.25487772) q[1];
sx q[1];
rz(-0.55995877) q[1];
rz(-1.9779491) q[3];
sx q[3];
rz(-2.4302357) q[3];
sx q[3];
rz(0.16669455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.031521156) q[2];
sx q[2];
rz(-1.4036125) q[2];
sx q[2];
rz(2.9384379) q[2];
rz(-1.3567694) q[3];
sx q[3];
rz(-1.2099096) q[3];
sx q[3];
rz(-1.8206966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88306952) q[0];
sx q[0];
rz(-1.8061545) q[0];
sx q[0];
rz(1.6254599) q[0];
rz(2.7698703) q[1];
sx q[1];
rz(-2.5860791) q[1];
sx q[1];
rz(2.1518167) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55732841) q[0];
sx q[0];
rz(-0.66270486) q[0];
sx q[0];
rz(1.2602379) q[0];
rz(2.3690574) q[2];
sx q[2];
rz(-1.4896817) q[2];
sx q[2];
rz(0.81743956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.02640662) q[1];
sx q[1];
rz(-2.7640928) q[1];
sx q[1];
rz(-1.0708234) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6521381) q[3];
sx q[3];
rz(-1.8453997) q[3];
sx q[3];
rz(-2.5413562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61701361) q[2];
sx q[2];
rz(-2.9144574) q[2];
sx q[2];
rz(3.0830834) q[2];
rz(1.2837563) q[3];
sx q[3];
rz(-0.82710281) q[3];
sx q[3];
rz(1.9730998) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7002895) q[0];
sx q[0];
rz(-1.0463725) q[0];
sx q[0];
rz(0.48536479) q[0];
rz(-0.46733388) q[1];
sx q[1];
rz(-2.554437) q[1];
sx q[1];
rz(2.4553518) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29083911) q[0];
sx q[0];
rz(-2.3510482) q[0];
sx q[0];
rz(-2.683862) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21659084) q[2];
sx q[2];
rz(-2.1859043) q[2];
sx q[2];
rz(-1.5768346) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8778692) q[1];
sx q[1];
rz(-2.4179966) q[1];
sx q[1];
rz(-0.76189009) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83250203) q[3];
sx q[3];
rz(-2.0568536) q[3];
sx q[3];
rz(1.6532236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.301959) q[2];
sx q[2];
rz(-2.0531824) q[2];
sx q[2];
rz(2.9097617) q[2];
rz(-2.2202282) q[3];
sx q[3];
rz(-1.6145584) q[3];
sx q[3];
rz(-0.025402633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(2.5068518) q[0];
sx q[0];
rz(-2.1379819) q[0];
sx q[0];
rz(0.64742175) q[0];
rz(-2.1000775) q[1];
sx q[1];
rz(-1.8687277) q[1];
sx q[1];
rz(1.1127068) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1739832) q[0];
sx q[0];
rz(-0.69794535) q[0];
sx q[0];
rz(2.3036454) q[0];
rz(1.7137063) q[2];
sx q[2];
rz(-1.630703) q[2];
sx q[2];
rz(-0.026345677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0354516) q[1];
sx q[1];
rz(-1.7454495) q[1];
sx q[1];
rz(-2.5779252) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1313354) q[3];
sx q[3];
rz(-1.105765) q[3];
sx q[3];
rz(0.51562515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.57704321) q[2];
sx q[2];
rz(-2.3304522) q[2];
sx q[2];
rz(0.88456336) q[2];
rz(2.1029643) q[3];
sx q[3];
rz(-2.0197155) q[3];
sx q[3];
rz(-1.60892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3603947) q[0];
sx q[0];
rz(-0.21112694) q[0];
sx q[0];
rz(-1.1398844) q[0];
rz(-2.3431011) q[1];
sx q[1];
rz(-1.7059098) q[1];
sx q[1];
rz(3.0671157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89635623) q[0];
sx q[0];
rz(-0.6707317) q[0];
sx q[0];
rz(1.9696495) q[0];
rz(2.3788667) q[2];
sx q[2];
rz(-1.8755091) q[2];
sx q[2];
rz(0.21108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.18596126) q[1];
sx q[1];
rz(-2.4096386) q[1];
sx q[1];
rz(0.47853985) q[1];
rz(-1.2182323) q[3];
sx q[3];
rz(-1.2855913) q[3];
sx q[3];
rz(-0.5087626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8188339) q[2];
sx q[2];
rz(-2.5311311) q[2];
sx q[2];
rz(-0.16767137) q[2];
rz(0.75336114) q[3];
sx q[3];
rz(-1.2246917) q[3];
sx q[3];
rz(1.0838449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3566256) q[0];
sx q[0];
rz(-1.4545119) q[0];
sx q[0];
rz(2.0941358) q[0];
rz(0.52182237) q[1];
sx q[1];
rz(-1.1818886) q[1];
sx q[1];
rz(-0.80400115) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2521236) q[0];
sx q[0];
rz(-1.7458545) q[0];
sx q[0];
rz(1.996422) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1267407) q[2];
sx q[2];
rz(-1.7202411) q[2];
sx q[2];
rz(-0.99260703) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40506002) q[1];
sx q[1];
rz(-1.1195464) q[1];
sx q[1];
rz(1.8644144) q[1];
x q[2];
rz(1.9639765) q[3];
sx q[3];
rz(-0.90709527) q[3];
sx q[3];
rz(-2.2507325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1890586) q[2];
sx q[2];
rz(-1.4905832) q[2];
sx q[2];
rz(-0.35183364) q[2];
rz(-0.072988836) q[3];
sx q[3];
rz(-1.9609541) q[3];
sx q[3];
rz(0.73968691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9214582) q[0];
sx q[0];
rz(-0.84311068) q[0];
sx q[0];
rz(-1.5699009) q[0];
rz(-0.72228557) q[1];
sx q[1];
rz(-1.1415569) q[1];
sx q[1];
rz(-1.7438181) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1565052) q[0];
sx q[0];
rz(-2.6596375) q[0];
sx q[0];
rz(-2.5768649) q[0];
rz(-pi) q[1];
rz(0.70571122) q[2];
sx q[2];
rz(-1.3967665) q[2];
sx q[2];
rz(1.1104012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.074894431) q[1];
sx q[1];
rz(-0.5579307) q[1];
sx q[1];
rz(2.9900223) q[1];
rz(1.8595376) q[3];
sx q[3];
rz(-2.7061314) q[3];
sx q[3];
rz(-0.2138775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.34652823) q[2];
sx q[2];
rz(-1.3115735) q[2];
sx q[2];
rz(-1.2947003) q[2];
rz(-1.0558225) q[3];
sx q[3];
rz(-2.081223) q[3];
sx q[3];
rz(1.0052217) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1718564) q[0];
sx q[0];
rz(-0.91005253) q[0];
sx q[0];
rz(-1.4685312) q[0];
rz(2.3079474) q[1];
sx q[1];
rz(-1.2794762) q[1];
sx q[1];
rz(1.5247482) q[1];
rz(1.1396617) q[2];
sx q[2];
rz(-0.79513245) q[2];
sx q[2];
rz(0.68383696) q[2];
rz(2.3232719) q[3];
sx q[3];
rz(-2.5723895) q[3];
sx q[3];
rz(1.7583917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
