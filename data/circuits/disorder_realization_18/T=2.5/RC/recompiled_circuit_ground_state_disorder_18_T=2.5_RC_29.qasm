OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.48489269) q[0];
sx q[0];
rz(-3.0527053) q[0];
sx q[0];
rz(-2.3645997) q[0];
rz(-2.4465893) q[1];
sx q[1];
rz(-0.13091317) q[1];
sx q[1];
rz(-0.34832365) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.191984) q[0];
sx q[0];
rz(-0.79553793) q[0];
sx q[0];
rz(2.863671) q[0];
rz(2.8856002) q[2];
sx q[2];
rz(-0.29689327) q[2];
sx q[2];
rz(0.96742899) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.45280734) q[1];
sx q[1];
rz(-1.0579351) q[1];
sx q[1];
rz(2.5393583) q[1];
rz(-1.6197085) q[3];
sx q[3];
rz(-1.9659316) q[3];
sx q[3];
rz(1.5103769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4942887) q[2];
sx q[2];
rz(-2.5755197) q[2];
sx q[2];
rz(1.0757793) q[2];
rz(-3.0268269) q[3];
sx q[3];
rz(-1.4817295) q[3];
sx q[3];
rz(-2.4858294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39605108) q[0];
sx q[0];
rz(-0.93463722) q[0];
sx q[0];
rz(2.4387687) q[0];
rz(-1.5952236) q[1];
sx q[1];
rz(-0.74235761) q[1];
sx q[1];
rz(2.5667618) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0718577) q[0];
sx q[0];
rz(-1.1848579) q[0];
sx q[0];
rz(-1.2161847) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61274075) q[2];
sx q[2];
rz(-1.1270971) q[2];
sx q[2];
rz(2.8982364) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0965754) q[1];
sx q[1];
rz(-2.5019766) q[1];
sx q[1];
rz(-2.9208697) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5089056) q[3];
sx q[3];
rz(-1.1345769) q[3];
sx q[3];
rz(-1.3496646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90291643) q[2];
sx q[2];
rz(-1.6409589) q[2];
sx q[2];
rz(0.73671663) q[2];
rz(0.77218562) q[3];
sx q[3];
rz(-0.18852791) q[3];
sx q[3];
rz(-2.7948715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18948983) q[0];
sx q[0];
rz(-2.0661856) q[0];
sx q[0];
rz(-1.9814251) q[0];
rz(0.53933764) q[1];
sx q[1];
rz(-2.5181006) q[1];
sx q[1];
rz(-0.070405237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5146014) q[0];
sx q[0];
rz(-0.87341269) q[0];
sx q[0];
rz(2.227289) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16164367) q[2];
sx q[2];
rz(-1.3934027) q[2];
sx q[2];
rz(-0.72189513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9765832) q[1];
sx q[1];
rz(-1.1900468) q[1];
sx q[1];
rz(1.5006719) q[1];
x q[2];
rz(0.65053026) q[3];
sx q[3];
rz(-1.1766947) q[3];
sx q[3];
rz(-0.89571834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77014852) q[2];
sx q[2];
rz(-1.7152717) q[2];
sx q[2];
rz(-2.9798129) q[2];
rz(2.3700304) q[3];
sx q[3];
rz(-3.0050889) q[3];
sx q[3];
rz(-1.0426883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55018392) q[0];
sx q[0];
rz(-0.78224459) q[0];
sx q[0];
rz(-0.25006205) q[0];
rz(-2.1024044) q[1];
sx q[1];
rz(-0.74165529) q[1];
sx q[1];
rz(2.340462) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.259819) q[0];
sx q[0];
rz(-1.7823813) q[0];
sx q[0];
rz(1.5552646) q[0];
rz(0.74015871) q[2];
sx q[2];
rz(-2.7974786) q[2];
sx q[2];
rz(2.0933852) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0274356) q[1];
sx q[1];
rz(-2.0832169) q[1];
sx q[1];
rz(1.6504662) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2053554) q[3];
sx q[3];
rz(-1.1574928) q[3];
sx q[3];
rz(2.6843021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4074126) q[2];
sx q[2];
rz(-1.4449747) q[2];
sx q[2];
rz(0.59047353) q[2];
rz(0.34481314) q[3];
sx q[3];
rz(-0.36546388) q[3];
sx q[3];
rz(-1.6706246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.19584) q[0];
sx q[0];
rz(-1.4386289) q[0];
sx q[0];
rz(1.6666743) q[0];
rz(-1.6189812) q[1];
sx q[1];
rz(-1.5432065) q[1];
sx q[1];
rz(0.408907) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98981114) q[0];
sx q[0];
rz(-1.5792002) q[0];
sx q[0];
rz(0.046812917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7402503) q[2];
sx q[2];
rz(-2.231488) q[2];
sx q[2];
rz(1.5398538) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6521367) q[1];
sx q[1];
rz(-2.5083275) q[1];
sx q[1];
rz(0.79852028) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4765396) q[3];
sx q[3];
rz(-0.89897147) q[3];
sx q[3];
rz(-0.49956027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54395479) q[2];
sx q[2];
rz(-0.13263098) q[2];
sx q[2];
rz(2.4908716) q[2];
rz(1.6808274) q[3];
sx q[3];
rz(-1.4131578) q[3];
sx q[3];
rz(0.4167324) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2862709) q[0];
sx q[0];
rz(-2.0725508) q[0];
sx q[0];
rz(-1.8977813) q[0];
rz(-2.6142201) q[1];
sx q[1];
rz(-1.4445883) q[1];
sx q[1];
rz(1.7605555) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52079372) q[0];
sx q[0];
rz(-0.24220151) q[0];
sx q[0];
rz(3.0907756) q[0];
x q[1];
rz(2.5041298) q[2];
sx q[2];
rz(-0.12345498) q[2];
sx q[2];
rz(-3.086498) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3917249) q[1];
sx q[1];
rz(-1.0942182) q[1];
sx q[1];
rz(-2.1471572) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45257537) q[3];
sx q[3];
rz(-2.0618478) q[3];
sx q[3];
rz(2.0752751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1018299) q[2];
sx q[2];
rz(-1.7888174) q[2];
sx q[2];
rz(1.9977894) q[2];
rz(0.77491289) q[3];
sx q[3];
rz(-1.9655922) q[3];
sx q[3];
rz(-2.9162143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5483122) q[0];
sx q[0];
rz(-0.53638023) q[0];
sx q[0];
rz(2.8450052) q[0];
rz(-0.93447724) q[1];
sx q[1];
rz(-2.2521033) q[1];
sx q[1];
rz(0.26781043) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6952301) q[0];
sx q[0];
rz(-0.3033378) q[0];
sx q[0];
rz(-2.3609451) q[0];
rz(-0.26508649) q[2];
sx q[2];
rz(-0.62787442) q[2];
sx q[2];
rz(1.6123079) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2909524) q[1];
sx q[1];
rz(-1.87504) q[1];
sx q[1];
rz(0.67493446) q[1];
x q[2];
rz(-1.7095836) q[3];
sx q[3];
rz(-1.7273081) q[3];
sx q[3];
rz(-1.9902347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3725793) q[2];
sx q[2];
rz(-2.09435) q[2];
sx q[2];
rz(0.11865842) q[2];
rz(2.3246121) q[3];
sx q[3];
rz(-1.1840362) q[3];
sx q[3];
rz(0.55678862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939222) q[0];
sx q[0];
rz(-0.477328) q[0];
sx q[0];
rz(-1.073786) q[0];
rz(-2.1444164) q[1];
sx q[1];
rz(-1.8938096) q[1];
sx q[1];
rz(1.006385) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2074315) q[0];
sx q[0];
rz(-1.6994358) q[0];
sx q[0];
rz(-3.1324803) q[0];
x q[1];
rz(0.7147737) q[2];
sx q[2];
rz(-0.48764418) q[2];
sx q[2];
rz(2.2478916) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80015228) q[1];
sx q[1];
rz(-1.7625215) q[1];
sx q[1];
rz(0.33445332) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48000042) q[3];
sx q[3];
rz(-2.1566803) q[3];
sx q[3];
rz(-0.65174499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61768475) q[2];
sx q[2];
rz(-1.074581) q[2];
sx q[2];
rz(-2.1802444) q[2];
rz(0.62054595) q[3];
sx q[3];
rz(-2.4397662) q[3];
sx q[3];
rz(-2.3294241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-2.0945702) q[0];
sx q[0];
rz(-0.71806878) q[0];
sx q[0];
rz(-2.4156003) q[0];
rz(-2.831931) q[1];
sx q[1];
rz(-2.1084712) q[1];
sx q[1];
rz(-0.18212254) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.108902) q[0];
sx q[0];
rz(-1.6918403) q[0];
sx q[0];
rz(2.3897479) q[0];
rz(2.0810037) q[2];
sx q[2];
rz(-2.0296718) q[2];
sx q[2];
rz(1.6422249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2606874) q[1];
sx q[1];
rz(-1.9753078) q[1];
sx q[1];
rz(-1.4244367) q[1];
rz(-pi) q[2];
rz(2.612669) q[3];
sx q[3];
rz(-2.4327603) q[3];
sx q[3];
rz(-0.20656997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2187659) q[2];
sx q[2];
rz(-1.7106235) q[2];
sx q[2];
rz(1.9975086) q[2];
rz(2.4036582) q[3];
sx q[3];
rz(-2.4149371) q[3];
sx q[3];
rz(1.1206333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69796193) q[0];
sx q[0];
rz(-1.8101036) q[0];
sx q[0];
rz(-2.4586082) q[0];
rz(-1.5890652) q[1];
sx q[1];
rz(-0.88337675) q[1];
sx q[1];
rz(2.24486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98600799) q[0];
sx q[0];
rz(-1.141662) q[0];
sx q[0];
rz(-1.3567423) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3405587) q[2];
sx q[2];
rz(-0.55998324) q[2];
sx q[2];
rz(-2.3513133) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7423279) q[1];
sx q[1];
rz(-0.15769449) q[1];
sx q[1];
rz(-2.5460973) q[1];
rz(-2.6399778) q[3];
sx q[3];
rz(-1.1715495) q[3];
sx q[3];
rz(-0.58540067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5985742) q[2];
sx q[2];
rz(-1.908058) q[2];
sx q[2];
rz(-0.42178556) q[2];
rz(2.5770523) q[3];
sx q[3];
rz(-0.83344236) q[3];
sx q[3];
rz(-0.88258266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17557872) q[0];
sx q[0];
rz(-1.6015552) q[0];
sx q[0];
rz(1.8764499) q[0];
rz(1.235678) q[1];
sx q[1];
rz(-2.0912248) q[1];
sx q[1];
rz(0.32774027) q[1];
rz(-2.7719978) q[2];
sx q[2];
rz(-1.7679749) q[2];
sx q[2];
rz(1.5883636) q[2];
rz(-2.0717703) q[3];
sx q[3];
rz(-1.9961431) q[3];
sx q[3];
rz(2.3186984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
