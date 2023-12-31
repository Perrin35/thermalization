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
rz(0.37101775) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(2.5420904) q[1];
sx q[1];
rz(11.190344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24297548) q[0];
sx q[0];
rz(-0.78865047) q[0];
sx q[0];
rz(2.2126161) q[0];
rz(2.5901428) q[2];
sx q[2];
rz(-0.79556888) q[2];
sx q[2];
rz(-1.431682) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4151167) q[1];
sx q[1];
rz(-1.2435902) q[1];
sx q[1];
rz(-0.16023689) q[1];
rz(2.3222378) q[3];
sx q[3];
rz(-1.5324394) q[3];
sx q[3];
rz(2.0824144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20724021) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(-1.9616615) q[0];
rz(0.99769366) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(0.72431272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99740072) q[0];
sx q[0];
rz(-1.2574982) q[0];
sx q[0];
rz(-2.2749167) q[0];
rz(-pi) q[1];
rz(2.9421147) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(-2.2897838) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69246768) q[1];
sx q[1];
rz(-1.3652703) q[1];
sx q[1];
rz(1.9301901) q[1];
rz(-pi) q[2];
rz(2.904326) q[3];
sx q[3];
rz(-1.7627343) q[3];
sx q[3];
rz(-2.3290079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90536845) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-2.779707) q[2];
rz(-3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(-0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996465) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(1.746159) q[0];
rz(0.46229258) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(1.9225072) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65788668) q[0];
sx q[0];
rz(-1.833796) q[0];
sx q[0];
rz(-1.2067243) q[0];
rz(-pi) q[1];
rz(-2.2723324) q[2];
sx q[2];
rz(-1.3420891) q[2];
sx q[2];
rz(1.9321835) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0137274) q[1];
sx q[1];
rz(-1.8121108) q[1];
sx q[1];
rz(-2.1167386) q[1];
rz(-pi) q[2];
rz(1.3027906) q[3];
sx q[3];
rz(-1.5558814) q[3];
sx q[3];
rz(-2.5915495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7200155) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(-0.56882632) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(-0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.22359426) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(-2.4261684) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(0.82675654) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2430192) q[0];
sx q[0];
rz(-1.1613701) q[0];
sx q[0];
rz(0.55222521) q[0];
x q[1];
rz(0.13917285) q[2];
sx q[2];
rz(-2.3690802) q[2];
sx q[2];
rz(-3.006209) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0485059) q[1];
sx q[1];
rz(-2.6566681) q[1];
sx q[1];
rz(-0.64694689) q[1];
x q[2];
rz(-2.4171962) q[3];
sx q[3];
rz(-1.3339692) q[3];
sx q[3];
rz(0.5326007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(-2.2391438) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(2.856423) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641159) q[0];
sx q[0];
rz(-2.6649094) q[0];
sx q[0];
rz(1.6633196) q[0];
x q[1];
rz(-2.6599777) q[2];
sx q[2];
rz(-2.4184605) q[2];
sx q[2];
rz(-2.8358104) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4080216) q[1];
sx q[1];
rz(-1.5703778) q[1];
sx q[1];
rz(1.2577406) q[1];
rz(0.22484803) q[3];
sx q[3];
rz(-0.41089155) q[3];
sx q[3];
rz(1.9083244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8482762) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(-0.53945333) q[2];
rz(0.30682492) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25813112) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(2.9845797) q[0];
rz(-0.69333386) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(-1.3670115) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54129823) q[0];
sx q[0];
rz(-1.6389264) q[0];
sx q[0];
rz(-3.1085204) q[0];
x q[1];
rz(0.579367) q[2];
sx q[2];
rz(-0.89951347) q[2];
sx q[2];
rz(0.52057779) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4945592) q[1];
sx q[1];
rz(-1.4899947) q[1];
sx q[1];
rz(-1.2809491) q[1];
rz(-pi) q[2];
rz(0.27512392) q[3];
sx q[3];
rz(-0.36558662) q[3];
sx q[3];
rz(2.5610353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(2.6599595) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994609) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(2.2677299) q[0];
rz(-2.6938687) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(1.9708995) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3160352) q[0];
sx q[0];
rz(-0.023840126) q[0];
sx q[0];
rz(2.8799879) q[0];
x q[1];
rz(0.41202338) q[2];
sx q[2];
rz(-3.014866) q[2];
sx q[2];
rz(2.3840981) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2611744) q[1];
sx q[1];
rz(-2.2194127) q[1];
sx q[1];
rz(-2.3748114) q[1];
x q[2];
rz(0.53415926) q[3];
sx q[3];
rz(-1.307752) q[3];
sx q[3];
rz(-2.0260889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0447023) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(-2.6664873) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8996745) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(0.29712594) q[0];
rz(-1.3946474) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(-2.4954605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4532115) q[0];
sx q[0];
rz(-1.6717981) q[0];
sx q[0];
rz(-1.7779011) q[0];
rz(-pi) q[1];
rz(-1.4023151) q[2];
sx q[2];
rz(-1.104276) q[2];
sx q[2];
rz(-0.70665765) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4375953) q[1];
sx q[1];
rz(-0.62347177) q[1];
sx q[1];
rz(-0.61203476) q[1];
rz(-pi) q[2];
rz(-2.2046702) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(-3.1261409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0769161) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(-2.440195) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(-1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69944537) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(0.79750693) q[0];
rz(2.6240255) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(0.10841766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71589564) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(2.2349368) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3181395) q[2];
sx q[2];
rz(-2.9571819) q[2];
sx q[2];
rz(0.84469634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0260967) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(-0.91067578) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2750862) q[3];
sx q[3];
rz(-1.7139072) q[3];
sx q[3];
rz(-2.0742311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(0.33995315) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(0.72559124) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323935) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-3.0864339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.178135) q[0];
sx q[0];
rz(-1.9473416) q[0];
sx q[0];
rz(-1.5012653) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5818985) q[2];
sx q[2];
rz(-1.392138) q[2];
sx q[2];
rz(2.6924804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0835417) q[1];
sx q[1];
rz(-2.1122879) q[1];
sx q[1];
rz(0.70152775) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.03667128) q[3];
sx q[3];
rz(-1.1150556) q[3];
sx q[3];
rz(2.7367221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(3.0040719) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4162083) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(-0.2086808) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(2.8490534) q[2];
sx q[2];
rz(-2.8164622) q[2];
sx q[2];
rz(-2.5278317) q[2];
rz(-0.24003868) q[3];
sx q[3];
rz(-1.6143027) q[3];
sx q[3];
rz(-1.3054813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
