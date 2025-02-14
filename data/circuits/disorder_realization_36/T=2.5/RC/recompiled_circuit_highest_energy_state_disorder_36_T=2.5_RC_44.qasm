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
rz(-0.041115887) q[0];
sx q[0];
rz(-1.9408722) q[0];
sx q[0];
rz(1.6706985) q[0];
rz(-0.38493758) q[1];
sx q[1];
rz(-0.5783143) q[1];
sx q[1];
rz(-0.90775031) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9616763) q[0];
sx q[0];
rz(-1.130234) q[0];
sx q[0];
rz(-2.5178686) q[0];
rz(-pi) q[1];
rz(-1.8706246) q[2];
sx q[2];
rz(-1.289164) q[2];
sx q[2];
rz(-2.9099438) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1536394) q[1];
sx q[1];
rz(-2.0551632) q[1];
sx q[1];
rz(-1.7820249) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5159881) q[3];
sx q[3];
rz(-1.4495088) q[3];
sx q[3];
rz(-2.9457322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3140728) q[2];
sx q[2];
rz(-1.4853518) q[2];
sx q[2];
rz(-2.8208288) q[2];
rz(-0.00094207923) q[3];
sx q[3];
rz(-0.92206803) q[3];
sx q[3];
rz(-0.48085406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8568521) q[0];
sx q[0];
rz(-0.72035235) q[0];
sx q[0];
rz(-2.5208933) q[0];
rz(1.8547828) q[1];
sx q[1];
rz(-1.6544673) q[1];
sx q[1];
rz(0.99498814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8061229) q[0];
sx q[0];
rz(-2.5988082) q[0];
sx q[0];
rz(-1.1837028) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1844145) q[2];
sx q[2];
rz(-0.52268302) q[2];
sx q[2];
rz(1.3764639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7409119) q[1];
sx q[1];
rz(-1.3865292) q[1];
sx q[1];
rz(1.1180941) q[1];
rz(-1.0424588) q[3];
sx q[3];
rz(-2.3064559) q[3];
sx q[3];
rz(2.5723815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93028012) q[2];
sx q[2];
rz(-1.899753) q[2];
sx q[2];
rz(-2.9105183) q[2];
rz(-1.6127582) q[3];
sx q[3];
rz(-1.4569747) q[3];
sx q[3];
rz(-2.5902364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.424161) q[0];
sx q[0];
rz(-1.9743974) q[0];
sx q[0];
rz(-0.396808) q[0];
rz(-1.1948168) q[1];
sx q[1];
rz(-1.8794941) q[1];
sx q[1];
rz(-1.6669115) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4572488) q[0];
sx q[0];
rz(-1.5593582) q[0];
sx q[0];
rz(-1.4829163) q[0];
rz(-0.72309317) q[2];
sx q[2];
rz(-2.4567219) q[2];
sx q[2];
rz(-1.2069281) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.296862) q[1];
sx q[1];
rz(-0.76557388) q[1];
sx q[1];
rz(2.0716448) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.01185939) q[3];
sx q[3];
rz(-1.9649155) q[3];
sx q[3];
rz(-0.32255641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.84612238) q[2];
sx q[2];
rz(-0.84540558) q[2];
sx q[2];
rz(1.4836614) q[2];
rz(2.1814003) q[3];
sx q[3];
rz(-2.5064597) q[3];
sx q[3];
rz(-0.6944164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8089777) q[0];
sx q[0];
rz(-1.806145) q[0];
sx q[0];
rz(-2.4131925) q[0];
rz(1.2737466) q[1];
sx q[1];
rz(-1.1796917) q[1];
sx q[1];
rz(1.5625928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1274783) q[0];
sx q[0];
rz(-1.593361) q[0];
sx q[0];
rz(-3.1382518) q[0];
rz(0.043508675) q[2];
sx q[2];
rz(-1.0232506) q[2];
sx q[2];
rz(-0.39247733) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3660903) q[1];
sx q[1];
rz(-1.7860768) q[1];
sx q[1];
rz(-1.4332814) q[1];
x q[2];
rz(-0.32890851) q[3];
sx q[3];
rz(-2.2136627) q[3];
sx q[3];
rz(-2.7907284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1100715) q[2];
sx q[2];
rz(-1.7379802) q[2];
sx q[2];
rz(2.9384379) q[2];
rz(1.3567694) q[3];
sx q[3];
rz(-1.9316831) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88306952) q[0];
sx q[0];
rz(-1.3354381) q[0];
sx q[0];
rz(-1.6254599) q[0];
rz(2.7698703) q[1];
sx q[1];
rz(-2.5860791) q[1];
sx q[1];
rz(2.1518167) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1976372) q[0];
sx q[0];
rz(-2.1966874) q[0];
sx q[0];
rz(0.23412378) q[0];
x q[1];
rz(1.4577652) q[2];
sx q[2];
rz(-2.340132) q[2];
sx q[2];
rz(-0.67455268) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0142611) q[1];
sx q[1];
rz(-1.3931573) q[1];
sx q[1];
rz(1.9056715) q[1];
rz(0.53984657) q[3];
sx q[3];
rz(-2.5858736) q[3];
sx q[3];
rz(-2.641895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.524579) q[2];
sx q[2];
rz(-2.9144574) q[2];
sx q[2];
rz(3.0830834) q[2];
rz(1.8578364) q[3];
sx q[3];
rz(-2.3144898) q[3];
sx q[3];
rz(-1.1684928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7002895) q[0];
sx q[0];
rz(-2.0952201) q[0];
sx q[0];
rz(-0.48536479) q[0];
rz(2.6742588) q[1];
sx q[1];
rz(-0.58715564) q[1];
sx q[1];
rz(-2.4553518) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29083911) q[0];
sx q[0];
rz(-2.3510482) q[0];
sx q[0];
rz(-2.683862) q[0];
rz(-pi) q[1];
rz(2.9250018) q[2];
sx q[2];
rz(-0.95568839) q[2];
sx q[2];
rz(1.5647581) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.972922) q[1];
sx q[1];
rz(-1.0712364) q[1];
sx q[1];
rz(-2.1184177) q[1];
rz(-pi) q[2];
rz(-0.62028168) q[3];
sx q[3];
rz(-0.93343319) q[3];
sx q[3];
rz(-0.31951527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83963362) q[2];
sx q[2];
rz(-2.0531824) q[2];
sx q[2];
rz(-2.9097617) q[2];
rz(-2.2202282) q[3];
sx q[3];
rz(-1.5270343) q[3];
sx q[3];
rz(-3.11619) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5068518) q[0];
sx q[0];
rz(-2.1379819) q[0];
sx q[0];
rz(0.64742175) q[0];
rz(-1.0415152) q[1];
sx q[1];
rz(-1.8687277) q[1];
sx q[1];
rz(-1.1127068) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96760945) q[0];
sx q[0];
rz(-0.69794535) q[0];
sx q[0];
rz(0.83794727) q[0];
rz(-1.4278863) q[2];
sx q[2];
rz(-1.5108897) q[2];
sx q[2];
rz(0.026345677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0354516) q[1];
sx q[1];
rz(-1.3961432) q[1];
sx q[1];
rz(-0.56366745) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0358488) q[3];
sx q[3];
rz(-1.5616284) q[3];
sx q[3];
rz(-2.0818215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5645494) q[2];
sx q[2];
rz(-2.3304522) q[2];
sx q[2];
rz(0.88456336) q[2];
rz(1.0386284) q[3];
sx q[3];
rz(-1.1218772) q[3];
sx q[3];
rz(-1.60892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7811979) q[0];
sx q[0];
rz(-0.21112694) q[0];
sx q[0];
rz(-2.0017083) q[0];
rz(-0.79849157) q[1];
sx q[1];
rz(-1.4356828) q[1];
sx q[1];
rz(-0.07447695) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1482753) q[0];
sx q[0];
rz(-1.3269985) q[0];
sx q[0];
rz(-2.2021342) q[0];
rz(0.76272599) q[2];
sx q[2];
rz(-1.8755091) q[2];
sx q[2];
rz(2.9305127) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18596126) q[1];
sx q[1];
rz(-0.73195405) q[1];
sx q[1];
rz(0.47853985) q[1];
x q[2];
rz(2.8387855) q[3];
sx q[3];
rz(-1.2330556) q[3];
sx q[3];
rz(-2.1827121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3227587) q[2];
sx q[2];
rz(-0.61046159) q[2];
sx q[2];
rz(-2.9739213) q[2];
rz(-2.3882315) q[3];
sx q[3];
rz(-1.916901) q[3];
sx q[3];
rz(2.0577478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-1.9597041) q[1];
sx q[1];
rz(-2.3375915) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3814731) q[0];
sx q[0];
rz(-1.1520885) q[0];
sx q[0];
rz(-0.19180723) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1754403) q[2];
sx q[2];
rz(-1.0217624) q[2];
sx q[2];
rz(2.6556478) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2968353) q[1];
sx q[1];
rz(-1.8342819) q[1];
sx q[1];
rz(-2.6729463) q[1];
rz(-2.439043) q[3];
sx q[3];
rz(-1.2642198) q[3];
sx q[3];
rz(0.93010139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1890586) q[2];
sx q[2];
rz(-1.4905832) q[2];
sx q[2];
rz(0.35183364) q[2];
rz(-3.0686038) q[3];
sx q[3];
rz(-1.1806386) q[3];
sx q[3];
rz(-2.4019057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9214582) q[0];
sx q[0];
rz(-2.298482) q[0];
sx q[0];
rz(-1.5716918) q[0];
rz(2.4193071) q[1];
sx q[1];
rz(-2.0000358) q[1];
sx q[1];
rz(1.7438181) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0972526) q[0];
sx q[0];
rz(-1.821479) q[0];
sx q[0];
rz(0.41608019) q[0];
rz(-2.8768853) q[2];
sx q[2];
rz(-2.4183344) q[2];
sx q[2];
rz(0.25991752) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8885615) q[1];
sx q[1];
rz(-1.0200046) q[1];
sx q[1];
rz(-1.4768449) q[1];
rz(-pi) q[2];
rz(3.0098823) q[3];
sx q[3];
rz(-1.1545106) q[3];
sx q[3];
rz(-0.10271969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.34652823) q[2];
sx q[2];
rz(-1.3115735) q[2];
sx q[2];
rz(1.8468924) q[2];
rz(-2.0857701) q[3];
sx q[3];
rz(-1.0603696) q[3];
sx q[3];
rz(-2.136371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(1.1718564) q[0];
sx q[0];
rz(-0.91005253) q[0];
sx q[0];
rz(-1.4685312) q[0];
rz(-2.3079474) q[1];
sx q[1];
rz(-1.8621164) q[1];
sx q[1];
rz(-1.6168445) q[1];
rz(-2.7387754) q[2];
sx q[2];
rz(-0.86502148) q[2];
sx q[2];
rz(-1.8765052) q[2];
rz(0.81832073) q[3];
sx q[3];
rz(-0.56920316) q[3];
sx q[3];
rz(-1.3832009) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
