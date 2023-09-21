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
rz(2.7603005) q[1];
sx q[1];
rz(-2.5420904) q[1];
sx q[1];
rz(-1.376027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24297548) q[0];
sx q[0];
rz(-2.3529422) q[0];
sx q[0];
rz(-2.2126161) q[0];
x q[1];
rz(0.55144989) q[2];
sx q[2];
rz(-0.79556888) q[2];
sx q[2];
rz(-1.7099107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.260533) q[1];
sx q[1];
rz(-0.36306371) q[1];
sx q[1];
rz(2.0102324) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.052470603) q[3];
sx q[3];
rz(-0.82004181) q[3];
sx q[3];
rz(-2.6657871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0573037) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(0.75254285) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(-2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.1799312) q[0];
rz(-2.143899) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(-0.72431272) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1441919) q[0];
sx q[0];
rz(-1.8840944) q[0];
sx q[0];
rz(0.86667592) q[0];
x q[1];
rz(-0.39436491) q[2];
sx q[2];
rz(-2.9260203) q[2];
sx q[2];
rz(-0.33286103) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.38072941) q[1];
sx q[1];
rz(-2.7298096) q[1];
sx q[1];
rz(-2.1058583) q[1];
rz(-pi) q[2];
x q[2];
rz(1.768126) q[3];
sx q[3];
rz(-1.8036246) q[3];
sx q[3];
rz(-0.80430921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-0.36188564) q[2];
rz(-3.0055255) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996465) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.3954337) q[0];
rz(0.46229258) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(1.2190855) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1299382) q[0];
sx q[0];
rz(-1.219795) q[0];
sx q[0];
rz(-0.28052335) q[0];
rz(-pi) q[1];
rz(-1.9169541) q[2];
sx q[2];
rz(-2.4097754) q[2];
sx q[2];
rz(3.0423622) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1278652) q[1];
sx q[1];
rz(-1.3294819) q[1];
sx q[1];
rz(2.1167386) q[1];
rz(3.1261256) q[3];
sx q[3];
rz(-1.3028212) q[3];
sx q[3];
rz(1.0166575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42157713) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(-1.6960309) q[2];
rz(-2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(3.0310757) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22359426) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(2.6233327) q[0];
rz(2.4261684) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(-0.82675654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2430192) q[0];
sx q[0];
rz(-1.9802226) q[0];
sx q[0];
rz(-0.55222521) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7051758) q[2];
sx q[2];
rz(-2.3339286) q[2];
sx q[2];
rz(0.32854167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0485059) q[1];
sx q[1];
rz(-2.6566681) q[1];
sx q[1];
rz(0.64694689) q[1];
rz(-1.8825674) q[3];
sx q[3];
rz(-2.2707553) q[3];
sx q[3];
rz(1.242897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15239079) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(3.0692696) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(-1.4962083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82692659) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(2.2391438) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(-2.856423) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0774768) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(1.6633196) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4816149) q[2];
sx q[2];
rz(-0.72313213) q[2];
sx q[2];
rz(-0.30578223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83708977) q[1];
sx q[1];
rz(-1.883852) q[1];
sx q[1];
rz(3.1411527) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6676335) q[3];
sx q[3];
rz(-1.1708461) q[3];
sx q[3];
rz(-1.6638343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(0.53945333) q[2];
rz(-0.30682492) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(-0.45421281) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(0.15701292) q[0];
rz(0.69333386) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(-1.7745811) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54129823) q[0];
sx q[0];
rz(-1.5026662) q[0];
sx q[0];
rz(3.1085204) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81142601) q[2];
sx q[2];
rz(-2.0138513) q[2];
sx q[2];
rz(-1.4366988) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4945592) q[1];
sx q[1];
rz(-1.651598) q[1];
sx q[1];
rz(-1.2809491) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4671765) q[3];
sx q[3];
rz(-1.9220256) q[3];
sx q[3];
rz(2.2675089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-2.7381251) q[2];
rz(2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(0.51923716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.8994609) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(-0.44772398) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(1.1706932) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3160352) q[0];
sx q[0];
rz(-3.1177525) q[0];
sx q[0];
rz(0.26160474) q[0];
rz(-0.11622073) q[2];
sx q[2];
rz(-1.5201609) q[2];
sx q[2];
rz(1.2223787) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2504187) q[1];
sx q[1];
rz(-2.182057) q[1];
sx q[1];
rz(-2.3120018) q[1];
rz(1.2675769) q[3];
sx q[3];
rz(-2.0847287) q[3];
sx q[3];
rz(0.60790387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(2.5308385) q[2];
rz(-0.47510535) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(-2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24191813) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(-2.8444667) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883811) q[0];
sx q[0];
rz(-1.4697945) q[0];
sx q[0];
rz(-1.3636916) q[0];
rz(-2.8201639) q[2];
sx q[2];
rz(-2.6476963) q[2];
sx q[2];
rz(-1.0682046) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1505193) q[1];
sx q[1];
rz(-1.0725613) q[1];
sx q[1];
rz(-1.1789765) q[1];
rz(0.93692245) q[3];
sx q[3];
rz(-0.86808944) q[3];
sx q[3];
rz(3.1261409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-0.45483744) q[2];
rz(-0.70139766) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(-2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69944537) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(0.79750693) q[0];
rz(-0.51756716) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(-0.10841766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90005504) q[0];
sx q[0];
rz(-0.90796048) q[0];
sx q[0];
rz(3.0682949) q[0];
rz(-pi) q[1];
rz(3.0949981) q[2];
sx q[2];
rz(-1.7492883) q[2];
sx q[2];
rz(2.5537234) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.11549599) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(-2.2309169) q[1];
x q[2];
rz(0.86650642) q[3];
sx q[3];
rz(-1.7139072) q[3];
sx q[3];
rz(2.0742311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-2.8016395) q[2];
rz(-2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(-0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6091992) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(-0.055158786) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96345761) q[0];
sx q[0];
rz(-1.9473416) q[0];
sx q[0];
rz(1.5012653) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8137384) q[2];
sx q[2];
rz(-0.58460669) q[2];
sx q[2];
rz(-2.2962928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0810869) q[1];
sx q[1];
rz(-2.2844237) q[1];
sx q[1];
rz(2.39141) q[1];
x q[2];
rz(3.1049214) q[3];
sx q[3];
rz(-1.1150556) q[3];
sx q[3];
rz(-0.40487056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(2.6514163) q[2];
rz(3.0040719) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(1.6677042) q[2];
sx q[2];
rz(-1.8816392) q[2];
sx q[2];
rz(0.30602602) q[2];
rz(-2.901554) q[3];
sx q[3];
rz(-1.5272899) q[3];
sx q[3];
rz(1.8361113) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];