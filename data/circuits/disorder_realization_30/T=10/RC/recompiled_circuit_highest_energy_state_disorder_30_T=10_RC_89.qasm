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
rz(-2.2089145) q[0];
sx q[0];
rz(2.7490766) q[0];
sx q[0];
rz(8.3277185) q[0];
rz(-1.8707844) q[1];
sx q[1];
rz(-1.9440117) q[1];
sx q[1];
rz(-1.4374011) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28406806) q[0];
sx q[0];
rz(-1.9211313) q[0];
sx q[0];
rz(0.68035462) q[0];
rz(1.497606) q[2];
sx q[2];
rz(-1.8022259) q[2];
sx q[2];
rz(2.481593) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.7858251) q[1];
sx q[1];
rz(-2.4406129) q[1];
sx q[1];
rz(2.088083) q[1];
rz(0.64793628) q[3];
sx q[3];
rz(-1.5625232) q[3];
sx q[3];
rz(3.0870761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2854332) q[2];
sx q[2];
rz(-0.23468748) q[2];
sx q[2];
rz(-2.6402546) q[2];
rz(-1.3242138) q[3];
sx q[3];
rz(-2.5960077) q[3];
sx q[3];
rz(-1.6747624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.685574) q[0];
sx q[0];
rz(-2.7360003) q[0];
sx q[0];
rz(1.4612041) q[0];
rz(0.17924084) q[1];
sx q[1];
rz(-2.3724809) q[1];
sx q[1];
rz(2.5040212) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44762627) q[0];
sx q[0];
rz(-1.3424831) q[0];
sx q[0];
rz(-2.8986321) q[0];
rz(-pi) q[1];
rz(0.045863002) q[2];
sx q[2];
rz(-1.8872607) q[2];
sx q[2];
rz(-0.72039159) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.04794807) q[1];
sx q[1];
rz(-2.1971697) q[1];
sx q[1];
rz(2.6996524) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6735503) q[3];
sx q[3];
rz(-0.7193102) q[3];
sx q[3];
rz(-1.5507789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3011938) q[2];
sx q[2];
rz(-1.0777148) q[2];
sx q[2];
rz(-0.044142874) q[2];
rz(0.61331493) q[3];
sx q[3];
rz(-1.3215348) q[3];
sx q[3];
rz(0.80673748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.089207) q[0];
sx q[0];
rz(-3.0969924) q[0];
sx q[0];
rz(1.3432304) q[0];
rz(-1.4187468) q[1];
sx q[1];
rz(-0.90444618) q[1];
sx q[1];
rz(0.80352965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1508458) q[0];
sx q[0];
rz(-0.68239337) q[0];
sx q[0];
rz(0.093390246) q[0];
x q[1];
rz(-2.5577689) q[2];
sx q[2];
rz(-1.1808504) q[2];
sx q[2];
rz(-2.4175298) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5403998) q[1];
sx q[1];
rz(-1.6777039) q[1];
sx q[1];
rz(-1.8501297) q[1];
x q[2];
rz(-0.50602956) q[3];
sx q[3];
rz(-1.2057736) q[3];
sx q[3];
rz(-0.50382352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.558305) q[2];
sx q[2];
rz(-0.95197695) q[2];
sx q[2];
rz(-1.0179016) q[2];
rz(-1.2244276) q[3];
sx q[3];
rz(-1.8095576) q[3];
sx q[3];
rz(-1.1289736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1430436) q[0];
sx q[0];
rz(-0.7989378) q[0];
sx q[0];
rz(0.96027389) q[0];
rz(0.1943365) q[1];
sx q[1];
rz(-1.7002218) q[1];
sx q[1];
rz(-3.1365373) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1086297) q[0];
sx q[0];
rz(-0.89187183) q[0];
sx q[0];
rz(-2.4358948) q[0];
rz(-pi) q[1];
rz(0.20273613) q[2];
sx q[2];
rz(-2.954287) q[2];
sx q[2];
rz(-0.024307131) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21987629) q[1];
sx q[1];
rz(-1.3520693) q[1];
sx q[1];
rz(-1.4552067) q[1];
rz(-pi) q[2];
rz(-1.051733) q[3];
sx q[3];
rz(-1.4421652) q[3];
sx q[3];
rz(2.7277456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.888835) q[2];
sx q[2];
rz(-1.5450666) q[2];
sx q[2];
rz(-1.4463536) q[2];
rz(1.3085922) q[3];
sx q[3];
rz(-1.7596217) q[3];
sx q[3];
rz(0.60477177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5204891) q[0];
sx q[0];
rz(-2.4172754) q[0];
sx q[0];
rz(-2.5526168) q[0];
rz(3.1336054) q[1];
sx q[1];
rz(-1.1050478) q[1];
sx q[1];
rz(-2.0569107) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046291489) q[0];
sx q[0];
rz(-1.3372375) q[0];
sx q[0];
rz(0.26136847) q[0];
x q[1];
rz(0.27954526) q[2];
sx q[2];
rz(-1.5790006) q[2];
sx q[2];
rz(-1.4943124) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.95653043) q[1];
sx q[1];
rz(-1.9635728) q[1];
sx q[1];
rz(-2.3529733) q[1];
x q[2];
rz(2.233423) q[3];
sx q[3];
rz(-1.623933) q[3];
sx q[3];
rz(-0.65180627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40953723) q[2];
sx q[2];
rz(-1.0169225) q[2];
sx q[2];
rz(1.766905) q[2];
rz(0.095666766) q[3];
sx q[3];
rz(-1.5547662) q[3];
sx q[3];
rz(0.63053757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0365527) q[0];
sx q[0];
rz(-2.6281272) q[0];
sx q[0];
rz(-1.3404982) q[0];
rz(-2.4758677) q[1];
sx q[1];
rz(-1.7981073) q[1];
sx q[1];
rz(-0.7241157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.501717) q[0];
sx q[0];
rz(-2.6611009) q[0];
sx q[0];
rz(2.2697422) q[0];
rz(-pi) q[1];
rz(-0.55990368) q[2];
sx q[2];
rz(-2.0717868) q[2];
sx q[2];
rz(-0.5463813) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16038475) q[1];
sx q[1];
rz(-2.0295459) q[1];
sx q[1];
rz(-1.4190516) q[1];
rz(3.136803) q[3];
sx q[3];
rz(-1.6671204) q[3];
sx q[3];
rz(-0.19118689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93806997) q[2];
sx q[2];
rz(-1.4504434) q[2];
sx q[2];
rz(-3.0435437) q[2];
rz(2.8010662) q[3];
sx q[3];
rz(-2.2398658) q[3];
sx q[3];
rz(0.19671973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.606474) q[0];
sx q[0];
rz(-2.4381194) q[0];
sx q[0];
rz(0.055543609) q[0];
rz(-0.37560383) q[1];
sx q[1];
rz(-2.3665078) q[1];
sx q[1];
rz(1.6966049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8880896) q[0];
sx q[0];
rz(-1.4652592) q[0];
sx q[0];
rz(2.3679033) q[0];
rz(-pi) q[1];
rz(-1.8246006) q[2];
sx q[2];
rz(-1.6217167) q[2];
sx q[2];
rz(-1.6737991) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0241131) q[1];
sx q[1];
rz(-1.584519) q[1];
sx q[1];
rz(-0.3776457) q[1];
rz(-1.1678004) q[3];
sx q[3];
rz(-1.1366362) q[3];
sx q[3];
rz(-0.94769615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.043341788) q[2];
sx q[2];
rz(-1.8378374) q[2];
sx q[2];
rz(-1.8275758) q[2];
rz(3.1116389) q[3];
sx q[3];
rz(-1.6311389) q[3];
sx q[3];
rz(0.31908527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7444262) q[0];
sx q[0];
rz(-2.501896) q[0];
sx q[0];
rz(-1.2514914) q[0];
rz(0.80998069) q[1];
sx q[1];
rz(-2.4046343) q[1];
sx q[1];
rz(1.4553778) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3446218) q[0];
sx q[0];
rz(-0.83528548) q[0];
sx q[0];
rz(1.263144) q[0];
rz(-pi) q[1];
rz(-0.52233861) q[2];
sx q[2];
rz(-1.6185158) q[2];
sx q[2];
rz(-3.0850956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67564671) q[1];
sx q[1];
rz(-0.94871229) q[1];
sx q[1];
rz(-0.81333604) q[1];
rz(2.3076644) q[3];
sx q[3];
rz(-2.0788611) q[3];
sx q[3];
rz(-2.4968387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84888419) q[2];
sx q[2];
rz(-1.3564738) q[2];
sx q[2];
rz(-1.9902309) q[2];
rz(3.1274318) q[3];
sx q[3];
rz(-1.3580946) q[3];
sx q[3];
rz(-1.2667228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8779811) q[0];
sx q[0];
rz(-2.4284555) q[0];
sx q[0];
rz(-0.94733316) q[0];
rz(-2.5454648) q[1];
sx q[1];
rz(-1.4214186) q[1];
sx q[1];
rz(1.5790342) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57347711) q[0];
sx q[0];
rz(-1.7422424) q[0];
sx q[0];
rz(0.72400064) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7315699) q[2];
sx q[2];
rz(-2.0410182) q[2];
sx q[2];
rz(1.1755113) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4217675) q[1];
sx q[1];
rz(-0.9523069) q[1];
sx q[1];
rz(-2.5773125) q[1];
rz(-2.4189255) q[3];
sx q[3];
rz(-2.8407466) q[3];
sx q[3];
rz(2.2900555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24937135) q[2];
sx q[2];
rz(-2.5669079) q[2];
sx q[2];
rz(2.6577244) q[2];
rz(-2.5441235) q[3];
sx q[3];
rz(-1.4474844) q[3];
sx q[3];
rz(-2.5634403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4018965) q[0];
sx q[0];
rz(-1.551349) q[0];
sx q[0];
rz(-0.36648146) q[0];
rz(-1.81555) q[1];
sx q[1];
rz(-2.1791024) q[1];
sx q[1];
rz(-2.0749626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90049998) q[0];
sx q[0];
rz(-1.0783195) q[0];
sx q[0];
rz(0.89333138) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5353955) q[2];
sx q[2];
rz(-1.5328888) q[2];
sx q[2];
rz(2.4647922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6910256) q[1];
sx q[1];
rz(-2.1524309) q[1];
sx q[1];
rz(-0.29009351) q[1];
rz(-pi) q[2];
x q[2];
rz(2.048038) q[3];
sx q[3];
rz(-1.541409) q[3];
sx q[3];
rz(-3.0813956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8804973) q[2];
sx q[2];
rz(-0.77793056) q[2];
sx q[2];
rz(2.0757389) q[2];
rz(-2.8737658) q[3];
sx q[3];
rz(-1.6949751) q[3];
sx q[3];
rz(2.6563787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5601226) q[0];
sx q[0];
rz(-1.0030092) q[0];
sx q[0];
rz(0.21448294) q[0];
rz(-0.32196925) q[1];
sx q[1];
rz(-1.8245158) q[1];
sx q[1];
rz(1.239924) q[1];
rz(2.8602045) q[2];
sx q[2];
rz(-1.5783969) q[2];
sx q[2];
rz(1.6498298) q[2];
rz(-0.81372502) q[3];
sx q[3];
rz(-2.57434) q[3];
sx q[3];
rz(2.5634585) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
