OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.934259) q[0];
sx q[0];
rz(-0.59036314) q[0];
sx q[0];
rz(-2.7705749) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(1.376027) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8126412) q[0];
sx q[0];
rz(-1.1321804) q[0];
sx q[0];
rz(-2.2493275) q[0];
x q[1];
rz(0.71557578) q[2];
sx q[2];
rz(-1.9543497) q[2];
sx q[2];
rz(0.54563145) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7264759) q[1];
sx q[1];
rz(-1.8980025) q[1];
sx q[1];
rz(0.16023689) q[1];
rz(-pi) q[2];
x q[2];
rz(1.514643) q[3];
sx q[3];
rz(-0.75222844) q[3];
sx q[3];
rz(2.5889531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(2.1526745) q[2];
rz(2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(0.83077103) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(-1.9616615) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(-2.4172799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82942078) q[0];
sx q[0];
rz(-2.2342626) q[0];
sx q[0];
rz(2.7396766) q[0];
rz(-pi) q[1];
rz(-2.9421147) q[2];
sx q[2];
rz(-1.4885159) q[2];
sx q[2];
rz(0.85180887) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1867379) q[1];
sx q[1];
rz(-1.2192982) q[1];
sx q[1];
rz(2.9224706) q[1];
x q[2];
rz(1.3734666) q[3];
sx q[3];
rz(-1.8036246) q[3];
sx q[3];
rz(0.80430921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(-0.36188564) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(-0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996465) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(1.746159) q[0];
rz(-0.46229258) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(-1.9225072) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65788668) q[0];
sx q[0];
rz(-1.833796) q[0];
sx q[0];
rz(-1.9348683) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8457882) q[2];
sx q[2];
rz(-0.89106262) q[2];
sx q[2];
rz(2.5909397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.29875007) q[1];
sx q[1];
rz(-1.0423653) q[1];
sx q[1];
rz(0.28038402) q[1];
rz(-pi) q[2];
rz(-3.1261256) q[3];
sx q[3];
rz(-1.8387715) q[3];
sx q[3];
rz(-2.1249352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(-2.5727663) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(0.7154243) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(0.82675654) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8985734) q[0];
sx q[0];
rz(-1.9802226) q[0];
sx q[0];
rz(-2.5893674) q[0];
x q[1];
rz(2.3739359) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(-1.8061639) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3863694) q[1];
sx q[1];
rz(-1.9519023) q[1];
sx q[1];
rz(-1.8783046) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2590253) q[3];
sx q[3];
rz(-2.2707553) q[3];
sx q[3];
rz(1.242897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3146661) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(1.1258874) q[0];
rz(2.2391438) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(0.28516969) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171705) q[0];
sx q[0];
rz(-1.5283913) q[0];
sx q[0];
rz(1.0958584) q[0];
x q[1];
rz(-0.4816149) q[2];
sx q[2];
rz(-2.4184605) q[2];
sx q[2];
rz(2.8358104) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3030745) q[1];
sx q[1];
rz(-2.8285366) q[1];
sx q[1];
rz(1.5721553) q[1];
rz(-pi) q[2];
rz(2.9167446) q[3];
sx q[3];
rz(-0.41089155) q[3];
sx q[3];
rz(-1.9083244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29331648) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(-0.53945333) q[2];
rz(-0.30682492) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(-0.69333386) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(-1.7745811) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.052664) q[0];
sx q[0];
rz(-3.0658709) q[0];
sx q[0];
rz(1.1195539) q[0];
rz(-2.1742646) q[2];
sx q[2];
rz(-2.2852995) q[2];
sx q[2];
rz(-2.8514903) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1937618) q[1];
sx q[1];
rz(-1.281922) q[1];
sx q[1];
rz(-3.0572901) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4671765) q[3];
sx q[3];
rz(-1.9220256) q[3];
sx q[3];
rz(-2.2675089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-2.7381251) q[2];
rz(0.48163313) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(0.51923716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.24213174) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(0.8738628) q[0];
rz(0.44772398) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(-1.1706932) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3160352) q[0];
sx q[0];
rz(-3.1177525) q[0];
sx q[0];
rz(2.8799879) q[0];
rz(-0.41202338) q[2];
sx q[2];
rz(-3.014866) q[2];
sx q[2];
rz(0.7574946) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16329855) q[1];
sx q[1];
rz(-0.98493176) q[1];
sx q[1];
rz(2.3818124) q[1];
x q[2];
rz(-2.6074334) q[3];
sx q[3];
rz(-1.8338406) q[3];
sx q[3];
rz(-1.1155038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(-0.61075413) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(-0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24191813) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(2.8444667) q[0];
rz(1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6883811) q[0];
sx q[0];
rz(-1.6717981) q[0];
sx q[0];
rz(1.7779011) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6693194) q[2];
sx q[2];
rz(-1.4204645) q[2];
sx q[2];
rz(2.3538102) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4375953) q[1];
sx q[1];
rz(-2.5181209) q[1];
sx q[1];
rz(0.61203476) q[1];
rz(2.5313247) q[3];
sx q[3];
rz(-2.2329997) q[3];
sx q[3];
rz(0.83412795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(0.45483744) q[2];
rz(2.440195) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(-1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69944537) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(3.033175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.425697) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(2.2349368) q[0];
rz(-1.3181395) q[2];
sx q[2];
rz(-2.9571819) q[2];
sx q[2];
rz(-0.84469634) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.11549599) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(0.91067578) q[1];
x q[2];
rz(-2.9547144) q[3];
sx q[3];
rz(-0.87516057) q[3];
sx q[3];
rz(-0.38284341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0124399) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-0.33995315) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(2.4627731) q[0];
rz(-0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(0.055158786) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.508651) q[0];
sx q[0];
rz(-1.635449) q[0];
sx q[0];
rz(0.37737329) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5818985) q[2];
sx q[2];
rz(-1.392138) q[2];
sx q[2];
rz(-0.44911227) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0810869) q[1];
sx q[1];
rz(-0.85716893) q[1];
sx q[1];
rz(-0.75018261) q[1];
rz(-2.0268029) q[3];
sx q[3];
rz(-1.6037233) q[3];
sx q[3];
rz(-1.959521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1577592) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(3.0040719) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(-0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(2.9329119) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(-2.8293777) q[2];
sx q[2];
rz(-1.4785462) q[2];
sx q[2];
rz(-1.2350456) q[2];
rz(1.5260074) q[3];
sx q[3];
rz(-1.3309892) q[3];
sx q[3];
rz(0.25467024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
