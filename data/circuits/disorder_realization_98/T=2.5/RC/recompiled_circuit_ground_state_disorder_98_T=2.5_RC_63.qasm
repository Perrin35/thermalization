OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.15220517) q[0];
sx q[0];
rz(2.9289991) q[0];
sx q[0];
rz(8.6906035) q[0];
rz(1.2165767) q[1];
sx q[1];
rz(-2.7956378) q[1];
sx q[1];
rz(0.024756519) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5261794) q[0];
sx q[0];
rz(-0.075881474) q[0];
sx q[0];
rz(-0.60657255) q[0];
rz(1.124618) q[2];
sx q[2];
rz(-2.8834448) q[2];
sx q[2];
rz(0.1218957) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5879257) q[1];
sx q[1];
rz(-0.30372639) q[1];
sx q[1];
rz(-0.17596735) q[1];
x q[2];
rz(1.3583899) q[3];
sx q[3];
rz(-1.5161361) q[3];
sx q[3];
rz(-2.2509991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88432246) q[2];
sx q[2];
rz(-1.4929644) q[2];
sx q[2];
rz(0.30882588) q[2];
rz(0.8257927) q[3];
sx q[3];
rz(-2.1335996) q[3];
sx q[3];
rz(2.3776313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17077133) q[0];
sx q[0];
rz(-1.9124557) q[0];
sx q[0];
rz(-2.1139076) q[0];
rz(0.81614256) q[1];
sx q[1];
rz(-2.0232537) q[1];
sx q[1];
rz(0.81248409) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4554169) q[0];
sx q[0];
rz(-0.38639613) q[0];
sx q[0];
rz(-1.1375582) q[0];
rz(-pi) q[1];
rz(1.5565926) q[2];
sx q[2];
rz(-1.9441295) q[2];
sx q[2];
rz(-2.4844784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5066487) q[1];
sx q[1];
rz(-2.6104984) q[1];
sx q[1];
rz(-1.8455161) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3489209) q[3];
sx q[3];
rz(-1.9066208) q[3];
sx q[3];
rz(-2.9909657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5673148) q[2];
sx q[2];
rz(-1.4894166) q[2];
sx q[2];
rz(2.2226649) q[2];
rz(-0.53141665) q[3];
sx q[3];
rz(-1.0454949) q[3];
sx q[3];
rz(-0.04118583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(-2.4384005) q[0];
sx q[0];
rz(-1.0696609) q[0];
sx q[0];
rz(2.0042787) q[0];
rz(-1.7510341) q[1];
sx q[1];
rz(-0.5568234) q[1];
sx q[1];
rz(2.4086319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.304564) q[0];
sx q[0];
rz(-0.51690716) q[0];
sx q[0];
rz(0.62249534) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3172313) q[2];
sx q[2];
rz(-2.3869053) q[2];
sx q[2];
rz(-3.1122308) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3324686) q[1];
sx q[1];
rz(-0.66624236) q[1];
sx q[1];
rz(0.82102832) q[1];
x q[2];
rz(-3.1350144) q[3];
sx q[3];
rz(-2.0195144) q[3];
sx q[3];
rz(-2.6980163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1375894) q[2];
sx q[2];
rz(-1.7294451) q[2];
sx q[2];
rz(-2.1027193) q[2];
rz(2.1393356) q[3];
sx q[3];
rz(-1.8398617) q[3];
sx q[3];
rz(2.8514298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.20659474) q[0];
sx q[0];
rz(-2.0552141) q[0];
sx q[0];
rz(2.2344053) q[0];
rz(-2.9838003) q[1];
sx q[1];
rz(-1.0148427) q[1];
sx q[1];
rz(-1.2812322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84128252) q[0];
sx q[0];
rz(-0.96409809) q[0];
sx q[0];
rz(-0.9228031) q[0];
rz(0.8550941) q[2];
sx q[2];
rz(-1.2629328) q[2];
sx q[2];
rz(1.5941053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43739) q[1];
sx q[1];
rz(-1.1466195) q[1];
sx q[1];
rz(0.68628879) q[1];
rz(-pi) q[2];
rz(0.82366039) q[3];
sx q[3];
rz(-2.3170217) q[3];
sx q[3];
rz(-2.0093105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8944051) q[2];
sx q[2];
rz(-2.2053568) q[2];
sx q[2];
rz(1.8518764) q[2];
rz(-1.7973409) q[3];
sx q[3];
rz(-1.684609) q[3];
sx q[3];
rz(-2.4414506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46071389) q[0];
sx q[0];
rz(-2.6432156) q[0];
sx q[0];
rz(-2.1233249) q[0];
rz(-2.0385888) q[1];
sx q[1];
rz(-0.7119199) q[1];
sx q[1];
rz(-1.8011372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.074891) q[0];
sx q[0];
rz(-1.929787) q[0];
sx q[0];
rz(-0.51407459) q[0];
rz(-pi) q[1];
rz(0.37968882) q[2];
sx q[2];
rz(-1.4799397) q[2];
sx q[2];
rz(1.0041217) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83507043) q[1];
sx q[1];
rz(-1.7002652) q[1];
sx q[1];
rz(-3.0135327) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0839447) q[3];
sx q[3];
rz(-2.1103922) q[3];
sx q[3];
rz(-3.0209013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5420142) q[2];
sx q[2];
rz(-0.88853637) q[2];
sx q[2];
rz(-1.0020024) q[2];
rz(2.4501948) q[3];
sx q[3];
rz(-1.9611497) q[3];
sx q[3];
rz(1.6208167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6028676) q[0];
sx q[0];
rz(-1.069101) q[0];
sx q[0];
rz(-2.5163203) q[0];
rz(-1.9484776) q[1];
sx q[1];
rz(-2.4517877) q[1];
sx q[1];
rz(0.56732059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0085908169) q[0];
sx q[0];
rz(-0.80256217) q[0];
sx q[0];
rz(0.24548291) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3151933) q[2];
sx q[2];
rz(-1.5763088) q[2];
sx q[2];
rz(-0.369584) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3750227) q[1];
sx q[1];
rz(-0.61900292) q[1];
sx q[1];
rz(-1.5260215) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17622275) q[3];
sx q[3];
rz(-1.7473975) q[3];
sx q[3];
rz(0.11208243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.304473) q[2];
sx q[2];
rz(-1.8581055) q[2];
sx q[2];
rz(-1.5597957) q[2];
rz(-0.61257735) q[3];
sx q[3];
rz(-1.9809096) q[3];
sx q[3];
rz(-0.15577236) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0953858) q[0];
sx q[0];
rz(-1.950773) q[0];
sx q[0];
rz(2.0147391) q[0];
rz(1.2696179) q[1];
sx q[1];
rz(-1.1879299) q[1];
sx q[1];
rz(-0.52106214) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4241735) q[0];
sx q[0];
rz(-1.9195286) q[0];
sx q[0];
rz(2.2762389) q[0];
rz(-pi) q[1];
rz(1.4748739) q[2];
sx q[2];
rz(-2.811199) q[2];
sx q[2];
rz(-1.0164193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.46494486) q[1];
sx q[1];
rz(-1.0527234) q[1];
sx q[1];
rz(2.2482616) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6005706) q[3];
sx q[3];
rz(-1.3329901) q[3];
sx q[3];
rz(-1.0353119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0236987) q[2];
sx q[2];
rz(-1.3641027) q[2];
sx q[2];
rz(1.2269616) q[2];
rz(-1.4620694) q[3];
sx q[3];
rz(-1.5173802) q[3];
sx q[3];
rz(0.44155651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11314497) q[0];
sx q[0];
rz(-1.0297091) q[0];
sx q[0];
rz(0.5471158) q[0];
rz(-0.088134915) q[1];
sx q[1];
rz(-1.3476177) q[1];
sx q[1];
rz(2.7190582) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9402661) q[0];
sx q[0];
rz(-0.28183386) q[0];
sx q[0];
rz(-1.4920217) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52253047) q[2];
sx q[2];
rz(-1.8944507) q[2];
sx q[2];
rz(0.98060545) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9692058) q[1];
sx q[1];
rz(-1.5017461) q[1];
sx q[1];
rz(0.056180908) q[1];
x q[2];
rz(3.001207) q[3];
sx q[3];
rz(-1.3891798) q[3];
sx q[3];
rz(-0.80206213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9162468) q[2];
sx q[2];
rz(-1.5044745) q[2];
sx q[2];
rz(2.8543191) q[2];
rz(-0.54287994) q[3];
sx q[3];
rz(-2.3714239) q[3];
sx q[3];
rz(0.058102593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7420237) q[0];
sx q[0];
rz(-0.68411198) q[0];
sx q[0];
rz(0.77731079) q[0];
rz(-0.9264535) q[1];
sx q[1];
rz(-1.9994206) q[1];
sx q[1];
rz(-1.3949589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1374717) q[0];
sx q[0];
rz(-0.51252675) q[0];
sx q[0];
rz(-0.88051535) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59916536) q[2];
sx q[2];
rz(-1.3627421) q[2];
sx q[2];
rz(0.64917246) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0277176) q[1];
sx q[1];
rz(-1.1166971) q[1];
sx q[1];
rz(0.36832955) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93392196) q[3];
sx q[3];
rz(-1.0524629) q[3];
sx q[3];
rz(2.3766488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42561439) q[2];
sx q[2];
rz(-1.7691111) q[2];
sx q[2];
rz(-1.0905637) q[2];
rz(2.1770832) q[3];
sx q[3];
rz(-2.1301853) q[3];
sx q[3];
rz(1.2151659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5665117) q[0];
sx q[0];
rz(-0.33841857) q[0];
sx q[0];
rz(1.5100719) q[0];
rz(-3.1072726) q[1];
sx q[1];
rz(-1.8020554) q[1];
sx q[1];
rz(0.43201772) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8513819) q[0];
sx q[0];
rz(-2.599596) q[0];
sx q[0];
rz(-1.4567514) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.081771227) q[2];
sx q[2];
rz(-1.2207979) q[2];
sx q[2];
rz(3.0584665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0610785) q[1];
sx q[1];
rz(-0.17560683) q[1];
sx q[1];
rz(-2.5422515) q[1];
rz(-pi) q[2];
rz(-1.179078) q[3];
sx q[3];
rz(-2.5234902) q[3];
sx q[3];
rz(-0.51392344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4711275) q[2];
sx q[2];
rz(-1.9894783) q[2];
sx q[2];
rz(2.067789) q[2];
rz(-0.18887575) q[3];
sx q[3];
rz(-2.5329068) q[3];
sx q[3];
rz(-0.3824189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043924532) q[0];
sx q[0];
rz(-2.8207939) q[0];
sx q[0];
rz(2.4378142) q[0];
rz(-1.7882998) q[1];
sx q[1];
rz(-0.67519338) q[1];
sx q[1];
rz(2.304362) q[1];
rz(0.64622579) q[2];
sx q[2];
rz(-2.2868509) q[2];
sx q[2];
rz(1.7522191) q[2];
rz(1.5021642) q[3];
sx q[3];
rz(-1.1123688) q[3];
sx q[3];
rz(0.49751626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
