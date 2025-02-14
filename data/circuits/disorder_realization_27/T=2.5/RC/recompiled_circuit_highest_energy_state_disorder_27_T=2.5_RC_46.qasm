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
rz(-1.6838411) q[0];
sx q[0];
rz(-0.10310752) q[0];
sx q[0];
rz(-2.4007894) q[0];
rz(-0.15751547) q[1];
sx q[1];
rz(-2.4060251) q[1];
sx q[1];
rz(-0.29677376) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4045217) q[0];
sx q[0];
rz(-1.1566136) q[0];
sx q[0];
rz(-0.41612723) q[0];
rz(-pi) q[1];
rz(2.4236083) q[2];
sx q[2];
rz(-1.3870001) q[2];
sx q[2];
rz(1.913547) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41488838) q[1];
sx q[1];
rz(-1.0203531) q[1];
sx q[1];
rz(-1.5797944) q[1];
x q[2];
rz(-1.2052676) q[3];
sx q[3];
rz(-2.6726279) q[3];
sx q[3];
rz(-0.381857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2510117) q[2];
sx q[2];
rz(-3.0457532) q[2];
sx q[2];
rz(-1.4061692) q[2];
rz(-2.3928394) q[3];
sx q[3];
rz(-0.65776062) q[3];
sx q[3];
rz(2.9296618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1521848) q[0];
sx q[0];
rz(-2.2693372) q[0];
sx q[0];
rz(2.8023791) q[0];
rz(2.5665414) q[1];
sx q[1];
rz(-1.1851858) q[1];
sx q[1];
rz(-0.38160479) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4313916) q[0];
sx q[0];
rz(-1.9862735) q[0];
sx q[0];
rz(-2.8870967) q[0];
x q[1];
rz(-1.9100045) q[2];
sx q[2];
rz(-1.8198065) q[2];
sx q[2];
rz(-0.23819645) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72605342) q[1];
sx q[1];
rz(-1.171863) q[1];
sx q[1];
rz(3.1311324) q[1];
x q[2];
rz(-0.53839243) q[3];
sx q[3];
rz(-0.88148553) q[3];
sx q[3];
rz(1.6834843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7316458) q[2];
sx q[2];
rz(-2.4531328) q[2];
sx q[2];
rz(2.6551969) q[2];
rz(-0.54388034) q[3];
sx q[3];
rz(-2.0982274) q[3];
sx q[3];
rz(2.9774408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.630702) q[0];
sx q[0];
rz(-0.4158026) q[0];
sx q[0];
rz(-2.1422332) q[0];
rz(2.4103145) q[1];
sx q[1];
rz(-0.46842289) q[1];
sx q[1];
rz(-2.9670002) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098487735) q[0];
sx q[0];
rz(-0.44666651) q[0];
sx q[0];
rz(1.5652324) q[0];
rz(-pi) q[1];
rz(0.83589696) q[2];
sx q[2];
rz(-0.23698254) q[2];
sx q[2];
rz(-1.8191172) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5970351) q[1];
sx q[1];
rz(-1.4824776) q[1];
sx q[1];
rz(3.0889325) q[1];
rz(1.1483795) q[3];
sx q[3];
rz(-2.6452521) q[3];
sx q[3];
rz(-1.5143192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.752855) q[2];
sx q[2];
rz(-2.2896705) q[2];
sx q[2];
rz(1.8641776) q[2];
rz(2.3169005) q[3];
sx q[3];
rz(-2.2970691) q[3];
sx q[3];
rz(2.9741014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.326062) q[0];
sx q[0];
rz(-2.675246) q[0];
sx q[0];
rz(-1.9933568) q[0];
rz(-2.8094021) q[1];
sx q[1];
rz(-0.59993184) q[1];
sx q[1];
rz(-1.0848328) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4702295) q[0];
sx q[0];
rz(-1.4309023) q[0];
sx q[0];
rz(-1.8663919) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33948989) q[2];
sx q[2];
rz(-0.84622806) q[2];
sx q[2];
rz(-0.17865114) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9124604) q[1];
sx q[1];
rz(-2.3082151) q[1];
sx q[1];
rz(1.499849) q[1];
rz(-pi) q[2];
rz(-2.6779019) q[3];
sx q[3];
rz(-1.1775908) q[3];
sx q[3];
rz(0.56305199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3499902) q[2];
sx q[2];
rz(-0.76452667) q[2];
sx q[2];
rz(3.1349658) q[2];
rz(0.22348063) q[3];
sx q[3];
rz(-1.0379182) q[3];
sx q[3];
rz(2.3124783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.70581907) q[0];
sx q[0];
rz(-2.3888102) q[0];
sx q[0];
rz(1.9594132) q[0];
rz(-1.301282) q[1];
sx q[1];
rz(-2.2186406) q[1];
sx q[1];
rz(-2.9822947) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7081122) q[0];
sx q[0];
rz(-1.5335521) q[0];
sx q[0];
rz(2.8539835) q[0];
rz(1.2550687) q[2];
sx q[2];
rz(-1.7509394) q[2];
sx q[2];
rz(2.2314928) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9827798) q[1];
sx q[1];
rz(-2.4430664) q[1];
sx q[1];
rz(0.79847269) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2377272) q[3];
sx q[3];
rz(-1.9124036) q[3];
sx q[3];
rz(1.2056392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84676877) q[2];
sx q[2];
rz(-0.20192768) q[2];
sx q[2];
rz(1.798604) q[2];
rz(-0.35074562) q[3];
sx q[3];
rz(-2.0028159) q[3];
sx q[3];
rz(-0.32575592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89966929) q[0];
sx q[0];
rz(-0.82829183) q[0];
sx q[0];
rz(-0.1668461) q[0];
rz(2.9150561) q[1];
sx q[1];
rz(-1.3275361) q[1];
sx q[1];
rz(-2.66364) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95889752) q[0];
sx q[0];
rz(-2.1271461) q[0];
sx q[0];
rz(-0.40370792) q[0];
rz(-pi) q[1];
rz(2.7331393) q[2];
sx q[2];
rz(-2.3801) q[2];
sx q[2];
rz(3.048427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8357827) q[1];
sx q[1];
rz(-2.1205582) q[1];
sx q[1];
rz(-0.98819403) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43059723) q[3];
sx q[3];
rz(-1.1599419) q[3];
sx q[3];
rz(1.2263067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53092521) q[2];
sx q[2];
rz(-1.3767367) q[2];
sx q[2];
rz(1.8492071) q[2];
rz(0.90062201) q[3];
sx q[3];
rz(-0.77018046) q[3];
sx q[3];
rz(-1.5463411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47386277) q[0];
sx q[0];
rz(-0.97309363) q[0];
sx q[0];
rz(-1.5245755) q[0];
rz(-1.4432817) q[1];
sx q[1];
rz(-0.89569211) q[1];
sx q[1];
rz(0.5009833) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058280073) q[0];
sx q[0];
rz(-0.73375637) q[0];
sx q[0];
rz(-1.0007798) q[0];
rz(-1.8800182) q[2];
sx q[2];
rz(-2.1078601) q[2];
sx q[2];
rz(-1.352965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9364719) q[1];
sx q[1];
rz(-2.3354482) q[1];
sx q[1];
rz(0.84772005) q[1];
x q[2];
rz(2.768211) q[3];
sx q[3];
rz(-2.8115559) q[3];
sx q[3];
rz(-0.16597834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3524126) q[2];
sx q[2];
rz(-0.12426201) q[2];
sx q[2];
rz(2.6782356) q[2];
rz(-3.0739259) q[3];
sx q[3];
rz(-1.8384408) q[3];
sx q[3];
rz(-2.9786003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25245923) q[0];
sx q[0];
rz(-1.8661789) q[0];
sx q[0];
rz(-0.5109936) q[0];
rz(-2.4976318) q[1];
sx q[1];
rz(-1.124758) q[1];
sx q[1];
rz(0.28800979) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32569186) q[0];
sx q[0];
rz(-1.7625232) q[0];
sx q[0];
rz(-1.5317388) q[0];
rz(-pi) q[1];
rz(0.63960774) q[2];
sx q[2];
rz(-0.98312964) q[2];
sx q[2];
rz(-1.4195051) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0102699) q[1];
sx q[1];
rz(-1.8449218) q[1];
sx q[1];
rz(2.5585973) q[1];
rz(-pi) q[2];
rz(-2.7631825) q[3];
sx q[3];
rz(-1.6749951) q[3];
sx q[3];
rz(0.7507594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1945232) q[2];
sx q[2];
rz(-1.9788519) q[2];
sx q[2];
rz(2.9266749) q[2];
rz(1.8042709) q[3];
sx q[3];
rz(-0.53842068) q[3];
sx q[3];
rz(0.18068331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8682078) q[0];
sx q[0];
rz(-2.2053563) q[0];
sx q[0];
rz(-1.0816164) q[0];
rz(-2.1514905) q[1];
sx q[1];
rz(-1.5847881) q[1];
sx q[1];
rz(-0.33499151) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2401857) q[0];
sx q[0];
rz(-1.6044582) q[0];
sx q[0];
rz(-1.2478527) q[0];
x q[1];
rz(0.0063517687) q[2];
sx q[2];
rz(-1.9926903) q[2];
sx q[2];
rz(2.5308501) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8674842) q[1];
sx q[1];
rz(-1.47164) q[1];
sx q[1];
rz(1.4851557) q[1];
rz(-pi) q[2];
rz(-0.0084235245) q[3];
sx q[3];
rz(-1.971774) q[3];
sx q[3];
rz(-3.0110735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0830903) q[2];
sx q[2];
rz(-2.6164656) q[2];
sx q[2];
rz(2.7527909) q[2];
rz(-0.36924103) q[3];
sx q[3];
rz(-0.25596127) q[3];
sx q[3];
rz(0.84454876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.9487172) q[0];
sx q[0];
rz(-2.188864) q[0];
sx q[0];
rz(2.4396851) q[0];
rz(-1.5319872) q[1];
sx q[1];
rz(-2.46789) q[1];
sx q[1];
rz(-0.38175499) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8030539) q[0];
sx q[0];
rz(-0.10333867) q[0];
sx q[0];
rz(0.65729143) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8605531) q[2];
sx q[2];
rz(-2.7157985) q[2];
sx q[2];
rz(-0.38044924) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6541958) q[1];
sx q[1];
rz(-1.8113448) q[1];
sx q[1];
rz(-1.3122561) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3424861) q[3];
sx q[3];
rz(-1.0630084) q[3];
sx q[3];
rz(-2.4158709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6751042) q[2];
sx q[2];
rz(-1.0305104) q[2];
sx q[2];
rz(-2.4934736) q[2];
rz(2.134792) q[3];
sx q[3];
rz(-1.9961793) q[3];
sx q[3];
rz(0.072331585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5299912) q[0];
sx q[0];
rz(-1.4016822) q[0];
sx q[0];
rz(0.66521426) q[0];
rz(-2.8499659) q[1];
sx q[1];
rz(-1.3529774) q[1];
sx q[1];
rz(-1.1381961) q[1];
rz(1.9251346) q[2];
sx q[2];
rz(-1.0783429) q[2];
sx q[2];
rz(1.1781296) q[2];
rz(-2.9838647) q[3];
sx q[3];
rz(-1.1554416) q[3];
sx q[3];
rz(-0.7994061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
