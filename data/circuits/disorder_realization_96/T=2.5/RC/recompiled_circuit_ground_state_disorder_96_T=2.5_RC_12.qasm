OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.628217) q[0];
sx q[0];
rz(-1.4613419) q[0];
sx q[0];
rz(-2.2267377) q[0];
rz(2.272361) q[1];
sx q[1];
rz(-0.72328049) q[1];
sx q[1];
rz(-1.7887315) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9147089) q[0];
sx q[0];
rz(-2.0830724) q[0];
sx q[0];
rz(0.73988503) q[0];
x q[1];
rz(0.74341821) q[2];
sx q[2];
rz(-1.8525436) q[2];
sx q[2];
rz(2.8670058) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1754419) q[1];
sx q[1];
rz(-1.3805318) q[1];
sx q[1];
rz(-2.2504066) q[1];
rz(0.21463359) q[3];
sx q[3];
rz(-2.1691536) q[3];
sx q[3];
rz(-0.90902381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.73504084) q[2];
sx q[2];
rz(-0.81321365) q[2];
sx q[2];
rz(-2.1347894) q[2];
rz(-1.462228) q[3];
sx q[3];
rz(-1.4454602) q[3];
sx q[3];
rz(1.6026976) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2891069) q[0];
sx q[0];
rz(-2.2302581) q[0];
sx q[0];
rz(3.015633) q[0];
rz(1.1360629) q[1];
sx q[1];
rz(-1.845153) q[1];
sx q[1];
rz(-2.0819285) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9618767) q[0];
sx q[0];
rz(-1.5715412) q[0];
sx q[0];
rz(-3.1396833) q[0];
rz(-pi) q[1];
rz(-2.4723889) q[2];
sx q[2];
rz(-1.0604217) q[2];
sx q[2];
rz(2.3056895) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4821533) q[1];
sx q[1];
rz(-1.1350613) q[1];
sx q[1];
rz(-2.3547291) q[1];
x q[2];
rz(-2.7823388) q[3];
sx q[3];
rz(-1.6174966) q[3];
sx q[3];
rz(0.31622313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58553592) q[2];
sx q[2];
rz(-1.5289331) q[2];
sx q[2];
rz(-2.5118828) q[2];
rz(-1.2398047) q[3];
sx q[3];
rz(-0.74435484) q[3];
sx q[3];
rz(2.0201717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12104163) q[0];
sx q[0];
rz(-2.910055) q[0];
sx q[0];
rz(1.9112021) q[0];
rz(2.4379099) q[1];
sx q[1];
rz(-2.2129702) q[1];
sx q[1];
rz(1.5812662) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9890922) q[0];
sx q[0];
rz(-1.1772424) q[0];
sx q[0];
rz(3.1253184) q[0];
x q[1];
rz(2.9283728) q[2];
sx q[2];
rz(-0.74951321) q[2];
sx q[2];
rz(2.9507278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0496975) q[1];
sx q[1];
rz(-1.9822243) q[1];
sx q[1];
rz(2.0446214) q[1];
x q[2];
rz(2.7625691) q[3];
sx q[3];
rz(-2.3190917) q[3];
sx q[3];
rz(-1.8387972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8065717) q[2];
sx q[2];
rz(-2.0609914) q[2];
sx q[2];
rz(3.0873155) q[2];
rz(2.055376) q[3];
sx q[3];
rz(-2.351206) q[3];
sx q[3];
rz(2.6495892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739968) q[0];
sx q[0];
rz(-3.1258899) q[0];
sx q[0];
rz(2.1570461) q[0];
rz(-1.4060075) q[1];
sx q[1];
rz(-1.5407298) q[1];
sx q[1];
rz(-0.23449177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0011976) q[0];
sx q[0];
rz(-0.89522982) q[0];
sx q[0];
rz(2.7864149) q[0];
rz(-0.245204) q[2];
sx q[2];
rz(-1.2916995) q[2];
sx q[2];
rz(-0.14268219) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6813184) q[1];
sx q[1];
rz(-1.7806306) q[1];
sx q[1];
rz(-0.55822452) q[1];
x q[2];
rz(1.6164833) q[3];
sx q[3];
rz(-0.57668466) q[3];
sx q[3];
rz(-0.2660584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24388127) q[2];
sx q[2];
rz(-2.477024) q[2];
sx q[2];
rz(0.74771869) q[2];
rz(-2.054935) q[3];
sx q[3];
rz(-2.1472411) q[3];
sx q[3];
rz(-2.7896816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9413167) q[0];
sx q[0];
rz(-0.87012297) q[0];
sx q[0];
rz(-1.548832) q[0];
rz(0.13762711) q[1];
sx q[1];
rz(-2.177114) q[1];
sx q[1];
rz(-0.67759883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7279997) q[0];
sx q[0];
rz(-1.5847854) q[0];
sx q[0];
rz(1.5634693) q[0];
rz(2.9477411) q[2];
sx q[2];
rz(-1.7192013) q[2];
sx q[2];
rz(-2.1231164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92050925) q[1];
sx q[1];
rz(-1.5649321) q[1];
sx q[1];
rz(1.4339892) q[1];
x q[2];
rz(0.26664385) q[3];
sx q[3];
rz(-2.3867749) q[3];
sx q[3];
rz(2.9637869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.64348334) q[2];
sx q[2];
rz(-0.52916873) q[2];
sx q[2];
rz(0.37437487) q[2];
rz(-0.74786413) q[3];
sx q[3];
rz(-1.6467843) q[3];
sx q[3];
rz(1.098746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13136524) q[0];
sx q[0];
rz(-1.3656728) q[0];
sx q[0];
rz(-2.9170872) q[0];
rz(-2.7911216) q[1];
sx q[1];
rz(-1.1843362) q[1];
sx q[1];
rz(-1.9559466) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23100863) q[0];
sx q[0];
rz(-1.5155927) q[0];
sx q[0];
rz(-1.7189222) q[0];
rz(1.6669438) q[2];
sx q[2];
rz(-2.0324869) q[2];
sx q[2];
rz(0.27335121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9486432) q[1];
sx q[1];
rz(-1.99896) q[1];
sx q[1];
rz(-0.74775006) q[1];
x q[2];
rz(-0.77187431) q[3];
sx q[3];
rz(-2.3079685) q[3];
sx q[3];
rz(-1.9825359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6422358) q[2];
sx q[2];
rz(-0.88358742) q[2];
sx q[2];
rz(2.9767766) q[2];
rz(0.76588255) q[3];
sx q[3];
rz(-2.3131148) q[3];
sx q[3];
rz(-0.59144056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0484475) q[0];
sx q[0];
rz(-3.1146545) q[0];
sx q[0];
rz(0.41937605) q[0];
rz(1.5834437) q[1];
sx q[1];
rz(-0.42834586) q[1];
sx q[1];
rz(-1.8289808) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93363491) q[0];
sx q[0];
rz(-0.61481732) q[0];
sx q[0];
rz(-1.2864248) q[0];
rz(-pi) q[1];
rz(-2.3840226) q[2];
sx q[2];
rz(-2.295864) q[2];
sx q[2];
rz(-0.49882327) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.928682) q[1];
sx q[1];
rz(-2.8281459) q[1];
sx q[1];
rz(2.2261966) q[1];
rz(-0.26793865) q[3];
sx q[3];
rz(-1.6959785) q[3];
sx q[3];
rz(-0.94733688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7340362) q[2];
sx q[2];
rz(-2.4702256) q[2];
sx q[2];
rz(-2.7835795) q[2];
rz(-0.19896209) q[3];
sx q[3];
rz(-0.81148654) q[3];
sx q[3];
rz(0.097804047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1426549) q[0];
sx q[0];
rz(-1.0336579) q[0];
sx q[0];
rz(0.79310966) q[0];
rz(2.2456887) q[1];
sx q[1];
rz(-1.221012) q[1];
sx q[1];
rz(0.47826794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5332141) q[0];
sx q[0];
rz(-1.694931) q[0];
sx q[0];
rz(0.45093765) q[0];
rz(0.68217607) q[2];
sx q[2];
rz(-2.6016005) q[2];
sx q[2];
rz(0.021395187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7859808) q[1];
sx q[1];
rz(-1.1063687) q[1];
sx q[1];
rz(-2.8465438) q[1];
rz(-pi) q[2];
rz(-1.7630799) q[3];
sx q[3];
rz(-2.6182976) q[3];
sx q[3];
rz(1.7442601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.521495) q[2];
sx q[2];
rz(-1.5767117) q[2];
sx q[2];
rz(-3.133797) q[2];
rz(-1.1826285) q[3];
sx q[3];
rz(-1.3820442) q[3];
sx q[3];
rz(-1.8462605) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9919392) q[0];
sx q[0];
rz(-2.6872334) q[0];
sx q[0];
rz(0.40865189) q[0];
rz(1.0952134) q[1];
sx q[1];
rz(-1.6588255) q[1];
sx q[1];
rz(0.85831395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3214071) q[0];
sx q[0];
rz(-1.6371861) q[0];
sx q[0];
rz(1.7184577) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3038396) q[2];
sx q[2];
rz(-0.46707312) q[2];
sx q[2];
rz(2.2427223) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.44659251) q[1];
sx q[1];
rz(-1.4607251) q[1];
sx q[1];
rz(-1.7677444) q[1];
x q[2];
rz(2.8204042) q[3];
sx q[3];
rz(-2.5192755) q[3];
sx q[3];
rz(-2.8701607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7909164) q[2];
sx q[2];
rz(-0.75296593) q[2];
sx q[2];
rz(-0.6616627) q[2];
rz(2.3412797) q[3];
sx q[3];
rz(-1.9789663) q[3];
sx q[3];
rz(0.77320981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.8372339) q[0];
sx q[0];
rz(-0.14166129) q[0];
sx q[0];
rz(0.71453553) q[0];
rz(-2.6658658) q[1];
sx q[1];
rz(-2.5560502) q[1];
sx q[1];
rz(2.661396) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6640319) q[0];
sx q[0];
rz(-0.98958221) q[0];
sx q[0];
rz(-2.1607735) q[0];
rz(-pi) q[1];
rz(0.33145655) q[2];
sx q[2];
rz(-0.94539795) q[2];
sx q[2];
rz(-1.7632222) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.0968761) q[1];
sx q[1];
rz(-0.72511092) q[1];
sx q[1];
rz(2.7027351) q[1];
rz(-pi) q[2];
rz(-2.0065745) q[3];
sx q[3];
rz(-1.3131071) q[3];
sx q[3];
rz(-0.99761744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6314038) q[2];
sx q[2];
rz(-1.1978585) q[2];
sx q[2];
rz(-2.8655444) q[2];
rz(0.43802842) q[3];
sx q[3];
rz(-1.6249526) q[3];
sx q[3];
rz(2.5146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.940687) q[0];
sx q[0];
rz(-1.9872266) q[0];
sx q[0];
rz(-3.0845597) q[0];
rz(3.1055462) q[1];
sx q[1];
rz(-1.6135975) q[1];
sx q[1];
rz(1.5560908) q[1];
rz(1.0771607) q[2];
sx q[2];
rz(-1.8133327) q[2];
sx q[2];
rz(1.2015154) q[2];
rz(2.6556426) q[3];
sx q[3];
rz(-0.28578352) q[3];
sx q[3];
rz(-1.8430229) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
