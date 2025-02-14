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
rz(3.0313015) q[0];
sx q[0];
rz(-1.8561441) q[0];
sx q[0];
rz(-2.8228446) q[0];
rz(4.3238001) q[1];
sx q[1];
rz(4.8053513) q[1];
sx q[1];
rz(4.3045192) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1308474) q[0];
sx q[0];
rz(-1.3358572) q[0];
sx q[0];
rz(-1.2287336) q[0];
rz(-3.0672795) q[2];
sx q[2];
rz(-2.9469159) q[2];
sx q[2];
rz(-3.1271324) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9388814) q[1];
sx q[1];
rz(-1.8716646) q[1];
sx q[1];
rz(1.8905413) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.031008677) q[3];
sx q[3];
rz(-0.57511273) q[3];
sx q[3];
rz(-1.9048579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.17243871) q[2];
sx q[2];
rz(-2.6814851) q[2];
sx q[2];
rz(-3.0058506) q[2];
rz(0.33484778) q[3];
sx q[3];
rz(-1.2232774) q[3];
sx q[3];
rz(2.7443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927476) q[0];
sx q[0];
rz(-1.0053585) q[0];
sx q[0];
rz(0.94648615) q[0];
rz(1.1874366) q[1];
sx q[1];
rz(-0.58381909) q[1];
sx q[1];
rz(0.28396398) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89881247) q[0];
sx q[0];
rz(-1.7526585) q[0];
sx q[0];
rz(2.9480431) q[0];
rz(-1.2912441) q[2];
sx q[2];
rz(-2.4233305) q[2];
sx q[2];
rz(-2.5606683) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.88017948) q[1];
sx q[1];
rz(-1.1481667) q[1];
sx q[1];
rz(-1.4201866) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3881298) q[3];
sx q[3];
rz(-1.8947269) q[3];
sx q[3];
rz(2.7517954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9339319) q[2];
sx q[2];
rz(-1.8457103) q[2];
sx q[2];
rz(2.3409823) q[2];
rz(1.8396395) q[3];
sx q[3];
rz(-1.1336361) q[3];
sx q[3];
rz(-3.083526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.12450739) q[0];
sx q[0];
rz(-0.43211102) q[0];
sx q[0];
rz(-2.3486163) q[0];
rz(0.68603459) q[1];
sx q[1];
rz(-1.425309) q[1];
sx q[1];
rz(-2.3763903) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3065465) q[0];
sx q[0];
rz(-1.6289428) q[0];
sx q[0];
rz(1.7821728) q[0];
rz(-pi) q[1];
rz(-0.051570895) q[2];
sx q[2];
rz(-1.5386536) q[2];
sx q[2];
rz(2.3598537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6531892) q[1];
sx q[1];
rz(-2.326366) q[1];
sx q[1];
rz(1.1865739) q[1];
rz(-pi) q[2];
rz(-2.2191677) q[3];
sx q[3];
rz(-2.2545345) q[3];
sx q[3];
rz(0.91372638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5746295) q[2];
sx q[2];
rz(-2.4506863) q[2];
sx q[2];
rz(-1.2464657) q[2];
rz(-1.1860819) q[3];
sx q[3];
rz(-1.3776774) q[3];
sx q[3];
rz(-0.20868364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74801159) q[0];
sx q[0];
rz(-1.3205386) q[0];
sx q[0];
rz(3.0082974) q[0];
rz(-0.26946274) q[1];
sx q[1];
rz(-0.91454426) q[1];
sx q[1];
rz(0.46636137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9892557) q[0];
sx q[0];
rz(-1.6902802) q[0];
sx q[0];
rz(-2.5463922) q[0];
rz(-2.1985717) q[2];
sx q[2];
rz(-1.4157214) q[2];
sx q[2];
rz(-0.53395203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0119446) q[1];
sx q[1];
rz(-0.62800558) q[1];
sx q[1];
rz(0.25917128) q[1];
x q[2];
rz(3.0250038) q[3];
sx q[3];
rz(-0.51672626) q[3];
sx q[3];
rz(-0.6445714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93126297) q[2];
sx q[2];
rz(-0.42079058) q[2];
sx q[2];
rz(-0.24410625) q[2];
rz(-1.3611475) q[3];
sx q[3];
rz(-2.1972392) q[3];
sx q[3];
rz(1.2066427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72531438) q[0];
sx q[0];
rz(-0.84753528) q[0];
sx q[0];
rz(-0.086294802) q[0];
rz(-1.3823973) q[1];
sx q[1];
rz(-1.0603797) q[1];
sx q[1];
rz(1.28654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3495958) q[0];
sx q[0];
rz(-1.9583869) q[0];
sx q[0];
rz(-2.9583065) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77877829) q[2];
sx q[2];
rz(-1.7302824) q[2];
sx q[2];
rz(2.3159112) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4205192) q[1];
sx q[1];
rz(-2.2700078) q[1];
sx q[1];
rz(0.10743227) q[1];
rz(-pi) q[2];
rz(-2.1819918) q[3];
sx q[3];
rz(-0.7905851) q[3];
sx q[3];
rz(0.039856002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7354108) q[2];
sx q[2];
rz(-1.7291131) q[2];
sx q[2];
rz(1.1253051) q[2];
rz(-2.3451037) q[3];
sx q[3];
rz(-0.73553604) q[3];
sx q[3];
rz(-2.9430732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35561246) q[0];
sx q[0];
rz(-2.8481843) q[0];
sx q[0];
rz(2.8134213) q[0];
rz(-2.7525821) q[1];
sx q[1];
rz(-1.5704472) q[1];
sx q[1];
rz(-1.6860115) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9302926) q[0];
sx q[0];
rz(-0.22575483) q[0];
sx q[0];
rz(0.36856099) q[0];
rz(3.0148618) q[2];
sx q[2];
rz(-1.7046844) q[2];
sx q[2];
rz(-2.9478879) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38992369) q[1];
sx q[1];
rz(-1.3397386) q[1];
sx q[1];
rz(-0.91737813) q[1];
rz(-pi) q[2];
rz(-0.041451575) q[3];
sx q[3];
rz(-1.2123479) q[3];
sx q[3];
rz(2.7849017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97633156) q[2];
sx q[2];
rz(-0.92504048) q[2];
sx q[2];
rz(-0.46249214) q[2];
rz(1.0235323) q[3];
sx q[3];
rz(-1.0547124) q[3];
sx q[3];
rz(-0.83824497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.8031215) q[0];
sx q[0];
rz(-1.4305038) q[0];
sx q[0];
rz(-0.17369239) q[0];
rz(-2.9202785) q[1];
sx q[1];
rz(-2.4559655) q[1];
sx q[1];
rz(1.7783222) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9598778) q[0];
sx q[0];
rz(-1.3893439) q[0];
sx q[0];
rz(-1.9812917) q[0];
x q[1];
rz(1.4837711) q[2];
sx q[2];
rz(-1.489893) q[2];
sx q[2];
rz(1.1822342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89243556) q[1];
sx q[1];
rz(-1.1022864) q[1];
sx q[1];
rz(2.4484642) q[1];
rz(-pi) q[2];
rz(2.2962988) q[3];
sx q[3];
rz(-2.7850683) q[3];
sx q[3];
rz(0.58024065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8814016) q[2];
sx q[2];
rz(-1.1611725) q[2];
sx q[2];
rz(-1.4455522) q[2];
rz(2.3675303) q[3];
sx q[3];
rz(-2.4249707) q[3];
sx q[3];
rz(-1.6707481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5792907) q[0];
sx q[0];
rz(-2.2023872) q[0];
sx q[0];
rz(-2.8939409) q[0];
rz(-2.472645) q[1];
sx q[1];
rz(-1.9525783) q[1];
sx q[1];
rz(-2.155969) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4307109) q[0];
sx q[0];
rz(-2.1257002) q[0];
sx q[0];
rz(1.4355833) q[0];
rz(0.56370391) q[2];
sx q[2];
rz(-1.8393751) q[2];
sx q[2];
rz(2.809066) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0528763) q[1];
sx q[1];
rz(-2.6357438) q[1];
sx q[1];
rz(-0.90152503) q[1];
x q[2];
rz(-1.5721442) q[3];
sx q[3];
rz(-2.165757) q[3];
sx q[3];
rz(-1.4524937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.140427) q[2];
sx q[2];
rz(-0.75158921) q[2];
sx q[2];
rz(1.2933732) q[2];
rz(-1.2398531) q[3];
sx q[3];
rz(-0.69609061) q[3];
sx q[3];
rz(2.4333439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2769315) q[0];
sx q[0];
rz(-1.4286574) q[0];
sx q[0];
rz(0.96555936) q[0];
rz(2.7652265) q[1];
sx q[1];
rz(-1.8195567) q[1];
sx q[1];
rz(-1.2300864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35621413) q[0];
sx q[0];
rz(-1.5206771) q[0];
sx q[0];
rz(0.26595195) q[0];
rz(2.4480341) q[2];
sx q[2];
rz(-2.6499814) q[2];
sx q[2];
rz(-2.7375321) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.849762) q[1];
sx q[1];
rz(-2.1755784) q[1];
sx q[1];
rz(-1.6563708) q[1];
rz(-1.6269095) q[3];
sx q[3];
rz(-0.75836042) q[3];
sx q[3];
rz(-2.2874934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2824715) q[2];
sx q[2];
rz(-0.66994795) q[2];
sx q[2];
rz(3.0100789) q[2];
rz(0.48111835) q[3];
sx q[3];
rz(-2.2115464) q[3];
sx q[3];
rz(0.49829811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.1199101) q[0];
sx q[0];
rz(-0.57761884) q[0];
sx q[0];
rz(2.6525894) q[0];
rz(-1.8148212) q[1];
sx q[1];
rz(-1.2811456) q[1];
sx q[1];
rz(0.16924032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45825935) q[0];
sx q[0];
rz(-1.1348327) q[0];
sx q[0];
rz(0.13157121) q[0];
x q[1];
rz(-0.18074482) q[2];
sx q[2];
rz(-1.482748) q[2];
sx q[2];
rz(3.0918025) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20470141) q[1];
sx q[1];
rz(-1.9866148) q[1];
sx q[1];
rz(-0.73483006) q[1];
x q[2];
rz(-2.333913) q[3];
sx q[3];
rz(-2.7480304) q[3];
sx q[3];
rz(-0.53787947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5532316) q[2];
sx q[2];
rz(-1.0855805) q[2];
sx q[2];
rz(2.9300743) q[2];
rz(-1.5254947) q[3];
sx q[3];
rz(-2.1414521) q[3];
sx q[3];
rz(-2.8741527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78032988) q[0];
sx q[0];
rz(-1.4365256) q[0];
sx q[0];
rz(-2.8509675) q[0];
rz(-2.3296539) q[1];
sx q[1];
rz(-0.62807905) q[1];
sx q[1];
rz(-1.5608578) q[1];
rz(-2.7079034) q[2];
sx q[2];
rz(-1.5885175) q[2];
sx q[2];
rz(1.3383404) q[2];
rz(2.3216861) q[3];
sx q[3];
rz(-1.0884566) q[3];
sx q[3];
rz(1.4940445) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
