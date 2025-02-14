OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71894574) q[0];
sx q[0];
rz(-2.5526241) q[0];
sx q[0];
rz(3.0229521) q[0];
rz(2.150382) q[1];
sx q[1];
rz(-1.041643) q[1];
sx q[1];
rz(0.51377327) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27072752) q[0];
sx q[0];
rz(-2.0692192) q[0];
sx q[0];
rz(0.99683572) q[0];
rz(1.6330326) q[2];
sx q[2];
rz(-2.0449491) q[2];
sx q[2];
rz(-1.4197503) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7440255) q[1];
sx q[1];
rz(-1.7850707) q[1];
sx q[1];
rz(-2.812723) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.051187201) q[3];
sx q[3];
rz(-0.59268206) q[3];
sx q[3];
rz(0.41646233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4673246) q[2];
sx q[2];
rz(-1.6215723) q[2];
sx q[2];
rz(-0.33217126) q[2];
rz(-0.25003555) q[3];
sx q[3];
rz(-1.9146405) q[3];
sx q[3];
rz(-1.8488319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2501204) q[0];
sx q[0];
rz(-1.8880867) q[0];
sx q[0];
rz(-0.44723311) q[0];
rz(-2.1121292) q[1];
sx q[1];
rz(-2.5599458) q[1];
sx q[1];
rz(3.0404125) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8120904) q[0];
sx q[0];
rz(-2.0116099) q[0];
sx q[0];
rz(-2.0157218) q[0];
rz(0.20002983) q[2];
sx q[2];
rz(-2.1406271) q[2];
sx q[2];
rz(2.9404158) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.14462073) q[1];
sx q[1];
rz(-2.640814) q[1];
sx q[1];
rz(2.5624496) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93191913) q[3];
sx q[3];
rz(-1.4804236) q[3];
sx q[3];
rz(-2.3109009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0826565) q[2];
sx q[2];
rz(-2.3857748) q[2];
sx q[2];
rz(-1.6015046) q[2];
rz(-0.61795175) q[3];
sx q[3];
rz(-0.58303419) q[3];
sx q[3];
rz(-1.9790953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59442941) q[0];
sx q[0];
rz(-1.0048486) q[0];
sx q[0];
rz(-0.51602236) q[0];
rz(1.0524606) q[1];
sx q[1];
rz(-2.2155589) q[1];
sx q[1];
rz(-1.1722391) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7965423) q[0];
sx q[0];
rz(-1.5597412) q[0];
sx q[0];
rz(0.42406908) q[0];
rz(-1.6555696) q[2];
sx q[2];
rz(-2.4293278) q[2];
sx q[2];
rz(0.9916457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7683923) q[1];
sx q[1];
rz(-2.6127763) q[1];
sx q[1];
rz(1.9463198) q[1];
rz(-pi) q[2];
rz(-1.5851791) q[3];
sx q[3];
rz(-0.77061117) q[3];
sx q[3];
rz(1.6028459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27511328) q[2];
sx q[2];
rz(-1.8786414) q[2];
sx q[2];
rz(-2.4647554) q[2];
rz(0.46946851) q[3];
sx q[3];
rz(-2.2468086) q[3];
sx q[3];
rz(1.2137871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4645828) q[0];
sx q[0];
rz(-1.0053758) q[0];
sx q[0];
rz(1.9418035) q[0];
rz(-0.59263539) q[1];
sx q[1];
rz(-1.0721782) q[1];
sx q[1];
rz(-1.6823654) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71818411) q[0];
sx q[0];
rz(-0.99331021) q[0];
sx q[0];
rz(-1.8052438) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8267385) q[2];
sx q[2];
rz(-2.3530966) q[2];
sx q[2];
rz(2.7389604) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45294083) q[1];
sx q[1];
rz(-2.4600852) q[1];
sx q[1];
rz(1.9863434) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6871485) q[3];
sx q[3];
rz(-2.4049149) q[3];
sx q[3];
rz(-0.018679623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.49372855) q[2];
sx q[2];
rz(-2.04546) q[2];
sx q[2];
rz(2.4596821) q[2];
rz(2.72825) q[3];
sx q[3];
rz(-1.5324493) q[3];
sx q[3];
rz(-1.9857064) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2979564) q[0];
sx q[0];
rz(-1.6652668) q[0];
sx q[0];
rz(0.27808878) q[0];
rz(0.9043215) q[1];
sx q[1];
rz(-1.4537289) q[1];
sx q[1];
rz(-0.98731891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073719115) q[0];
sx q[0];
rz(-0.6530531) q[0];
sx q[0];
rz(1.1730173) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9402307) q[2];
sx q[2];
rz(-0.1397396) q[2];
sx q[2];
rz(-2.1167378) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.029475538) q[1];
sx q[1];
rz(-1.6320845) q[1];
sx q[1];
rz(-1.1809096) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.868949) q[3];
sx q[3];
rz(-1.050625) q[3];
sx q[3];
rz(-0.70462728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.014293369) q[2];
sx q[2];
rz(-0.96724808) q[2];
sx q[2];
rz(-2.3822752) q[2];
rz(1.4826108) q[3];
sx q[3];
rz(-1.9448152) q[3];
sx q[3];
rz(-2.7827941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229014) q[0];
sx q[0];
rz(-2.3221115) q[0];
sx q[0];
rz(2.4500093) q[0];
rz(2.2870731) q[1];
sx q[1];
rz(-2.4311192) q[1];
sx q[1];
rz(2.7202594) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79165073) q[0];
sx q[0];
rz(-1.6393108) q[0];
sx q[0];
rz(-2.0035901) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3106903) q[2];
sx q[2];
rz(-0.38543265) q[2];
sx q[2];
rz(1.3192434) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9142575) q[1];
sx q[1];
rz(-2.7707639) q[1];
sx q[1];
rz(2.9079958) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60507994) q[3];
sx q[3];
rz(-1.7972094) q[3];
sx q[3];
rz(0.070158557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2104346) q[2];
sx q[2];
rz(-1.6633818) q[2];
sx q[2];
rz(0.20273905) q[2];
rz(1.5431131) q[3];
sx q[3];
rz(-2.8212382) q[3];
sx q[3];
rz(-1.4690442) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28432524) q[0];
sx q[0];
rz(-1.7131282) q[0];
sx q[0];
rz(0.39886928) q[0];
rz(-1.9786037) q[1];
sx q[1];
rz(-0.82610026) q[1];
sx q[1];
rz(2.861048) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4568854) q[0];
sx q[0];
rz(-2.5734757) q[0];
sx q[0];
rz(2.3756308) q[0];
rz(-pi) q[1];
rz(-2.7107096) q[2];
sx q[2];
rz(-2.2019349) q[2];
sx q[2];
rz(0.43090303) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.043634) q[1];
sx q[1];
rz(-0.53501883) q[1];
sx q[1];
rz(-2.19997) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3451398) q[3];
sx q[3];
rz(-1.1304677) q[3];
sx q[3];
rz(-1.8673107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3205388) q[2];
sx q[2];
rz(-2.1619449) q[2];
sx q[2];
rz(1.6542356) q[2];
rz(-0.94333831) q[3];
sx q[3];
rz(-2.2127559) q[3];
sx q[3];
rz(1.1613891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5281552) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(-0.50233895) q[0];
rz(0.54652864) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(2.6796403) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0620856) q[0];
sx q[0];
rz(-1.2794198) q[0];
sx q[0];
rz(1.0568922) q[0];
x q[1];
rz(-0.047872825) q[2];
sx q[2];
rz(-0.92000735) q[2];
sx q[2];
rz(1.6597468) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5344055) q[1];
sx q[1];
rz(-2.7173373) q[1];
sx q[1];
rz(1.3035745) q[1];
x q[2];
rz(0.021691247) q[3];
sx q[3];
rz(-2.6110161) q[3];
sx q[3];
rz(-2.2013045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97303566) q[2];
sx q[2];
rz(-0.23739561) q[2];
sx q[2];
rz(-2.9711704) q[2];
rz(2.6863344) q[3];
sx q[3];
rz(-1.5452496) q[3];
sx q[3];
rz(-0.71793238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1667204) q[0];
sx q[0];
rz(-1.6396739) q[0];
sx q[0];
rz(0.18044743) q[0];
rz(-2.6249053) q[1];
sx q[1];
rz(-1.5645212) q[1];
sx q[1];
rz(2.7076941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2045) q[0];
sx q[0];
rz(-0.70815101) q[0];
sx q[0];
rz(2.2592179) q[0];
rz(-pi) q[1];
rz(-2.216748) q[2];
sx q[2];
rz(-0.5415104) q[2];
sx q[2];
rz(-0.72586593) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89131195) q[1];
sx q[1];
rz(-1.6906594) q[1];
sx q[1];
rz(-2.2054172) q[1];
x q[2];
rz(-2.2614165) q[3];
sx q[3];
rz(-1.7978906) q[3];
sx q[3];
rz(0.86798641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5608998) q[2];
sx q[2];
rz(-2.7036724) q[2];
sx q[2];
rz(-1.8221347) q[2];
rz(-0.27160078) q[3];
sx q[3];
rz(-1.8235794) q[3];
sx q[3];
rz(-2.0297089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.44529799) q[0];
sx q[0];
rz(-0.91067186) q[0];
sx q[0];
rz(-1.4060422) q[0];
rz(-1.7156853) q[1];
sx q[1];
rz(-2.0366663) q[1];
sx q[1];
rz(-1.8941194) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1734764) q[0];
sx q[0];
rz(-2.0531468) q[0];
sx q[0];
rz(2.8456057) q[0];
x q[1];
rz(1.1389473) q[2];
sx q[2];
rz(-3.0122535) q[2];
sx q[2];
rz(2.6292588) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5312219) q[1];
sx q[1];
rz(-0.34451133) q[1];
sx q[1];
rz(-0.12126874) q[1];
rz(-2.5188451) q[3];
sx q[3];
rz(-1.0120262) q[3];
sx q[3];
rz(-1.2947242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4282816) q[2];
sx q[2];
rz(-2.9094628) q[2];
sx q[2];
rz(0.10078079) q[2];
rz(3.1320324) q[3];
sx q[3];
rz(-1.5848426) q[3];
sx q[3];
rz(-1.5338219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1495001) q[0];
sx q[0];
rz(-1.7763573) q[0];
sx q[0];
rz(-1.1471164) q[0];
rz(2.5485582) q[1];
sx q[1];
rz(-1.4321764) q[1];
sx q[1];
rz(2.166116) q[1];
rz(1.4040074) q[2];
sx q[2];
rz(-2.7381512) q[2];
sx q[2];
rz(2.9742407) q[2];
rz(1.5884052) q[3];
sx q[3];
rz(-1.3497769) q[3];
sx q[3];
rz(-2.6365437) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
