OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(0.36261121) q[0];
rz(-2.2244722) q[1];
sx q[1];
rz(-2.6511104) q[1];
sx q[1];
rz(-2.7999556) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4877019) q[0];
sx q[0];
rz(-0.39429769) q[0];
sx q[0];
rz(2.8173692) q[0];
rz(-0.50589675) q[2];
sx q[2];
rz(-2.4873173) q[2];
sx q[2];
rz(1.8088532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29400723) q[1];
sx q[1];
rz(-0.60677401) q[1];
sx q[1];
rz(1.0311544) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0658675) q[3];
sx q[3];
rz(-2.6945811) q[3];
sx q[3];
rz(1.3017553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59387702) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(-0.95648742) q[2];
rz(-0.18125136) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(-1.6199002) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-3.0759838) q[0];
sx q[0];
rz(1.99615) q[0];
rz(-1.068813) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(-2.4157445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23913357) q[0];
sx q[0];
rz(-2.4665138) q[0];
sx q[0];
rz(-2.549987) q[0];
x q[1];
rz(2.5475694) q[2];
sx q[2];
rz(-1.960264) q[2];
sx q[2];
rz(-2.4246852) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.010163807) q[1];
sx q[1];
rz(-1.5457075) q[1];
sx q[1];
rz(1.8797092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29048357) q[3];
sx q[3];
rz(-0.64290128) q[3];
sx q[3];
rz(-3.052352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(-2.144311) q[2];
rz(-1.7287792) q[3];
sx q[3];
rz(-1.0935067) q[3];
sx q[3];
rz(-2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41855758) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(-1.8380802) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(-2.2361141) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14074621) q[0];
sx q[0];
rz(-1.329639) q[0];
sx q[0];
rz(1.3748037) q[0];
rz(-2.3891874) q[2];
sx q[2];
rz(-1.89626) q[2];
sx q[2];
rz(2.4706555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29618759) q[1];
sx q[1];
rz(-2.2803109) q[1];
sx q[1];
rz(-0.097464949) q[1];
rz(-pi) q[2];
rz(3.0985673) q[3];
sx q[3];
rz(-0.85775162) q[3];
sx q[3];
rz(-0.42792861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3811615) q[2];
sx q[2];
rz(-2.9734549) q[2];
sx q[2];
rz(-2.1543489) q[2];
rz(-3.0958214) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(-0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0320597) q[0];
sx q[0];
rz(-2.4592472) q[0];
sx q[0];
rz(3.0155244) q[0];
rz(0.18445045) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(1.1674081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0569939) q[0];
sx q[0];
rz(-1.7718678) q[0];
sx q[0];
rz(-2.3720471) q[0];
rz(-2.6378176) q[2];
sx q[2];
rz(-1.7283895) q[2];
sx q[2];
rz(0.7047082) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5645204) q[1];
sx q[1];
rz(-0.92622354) q[1];
sx q[1];
rz(1.9034027) q[1];
rz(-pi) q[2];
rz(-2.4950124) q[3];
sx q[3];
rz(-1.0160035) q[3];
sx q[3];
rz(-0.085689714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.8644774) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(-2.5100822) q[2];
rz(2.946092) q[3];
sx q[3];
rz(-1.86444) q[3];
sx q[3];
rz(-2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(0.50267977) q[0];
rz(1.1133105) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(-0.46708333) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0665849) q[0];
sx q[0];
rz(-1.9420615) q[0];
sx q[0];
rz(2.9270372) q[0];
x q[1];
rz(-2.5555243) q[2];
sx q[2];
rz(-1.35891) q[2];
sx q[2];
rz(-0.52473247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.621532) q[1];
sx q[1];
rz(-0.48081765) q[1];
sx q[1];
rz(2.0562999) q[1];
rz(-pi) q[2];
rz(1.3777556) q[3];
sx q[3];
rz(-0.69460624) q[3];
sx q[3];
rz(0.65929268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(2.6494027) q[2];
rz(2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(-1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.870134) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(-3.0555994) q[0];
rz(-1.4280691) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(-1.496398) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2917568) q[0];
sx q[0];
rz(-2.09016) q[0];
sx q[0];
rz(-2.915328) q[0];
rz(-pi) q[1];
rz(-2.9005269) q[2];
sx q[2];
rz(-1.0450372) q[2];
sx q[2];
rz(-2.966553) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8346356) q[1];
sx q[1];
rz(-2.7863057) q[1];
sx q[1];
rz(1.6389695) q[1];
x q[2];
rz(3.1268678) q[3];
sx q[3];
rz(-2.428599) q[3];
sx q[3];
rz(-1.7409489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4795586) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(-0.39624828) q[2];
rz(-2.5475492) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(-1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632161) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(-2.0492045) q[0];
rz(-1.3573525) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(0.011627442) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25225885) q[0];
sx q[0];
rz(-0.47124915) q[0];
sx q[0];
rz(2.528119) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.02545698) q[2];
sx q[2];
rz(-1.5989132) q[2];
sx q[2];
rz(-0.97719976) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1639451) q[1];
sx q[1];
rz(-1.5944591) q[1];
sx q[1];
rz(-1.0246828) q[1];
x q[2];
rz(2.743268) q[3];
sx q[3];
rz(-0.71361226) q[3];
sx q[3];
rz(2.6834965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.45550436) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(-2.8536076) q[2];
rz(-1.1503495) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(-0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.31036723) q[0];
sx q[0];
rz(-0.53124017) q[0];
sx q[0];
rz(1.2121375) q[0];
rz(0.37777004) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(1.4935965) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15570116) q[0];
sx q[0];
rz(-0.37030664) q[0];
sx q[0];
rz(-2.5168688) q[0];
x q[1];
rz(1.9940894) q[2];
sx q[2];
rz(-0.14483843) q[2];
sx q[2];
rz(-2.7763979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6495325) q[1];
sx q[1];
rz(-0.82693716) q[1];
sx q[1];
rz(0.56708132) q[1];
x q[2];
rz(1.7842403) q[3];
sx q[3];
rz(-0.5118013) q[3];
sx q[3];
rz(1.3444927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7887855) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(-1.8898233) q[2];
rz(-2.2552323) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(-3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66824526) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(-2.9633203) q[0];
rz(0.072862236) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(-1.1791621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3423791) q[0];
sx q[0];
rz(-1.2782492) q[0];
sx q[0];
rz(-0.53036687) q[0];
rz(0.91782848) q[2];
sx q[2];
rz(-2.4036916) q[2];
sx q[2];
rz(2.3248621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.804783) q[1];
sx q[1];
rz(-2.0840624) q[1];
sx q[1];
rz(-1.2745162) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2004553) q[3];
sx q[3];
rz(-2.5459873) q[3];
sx q[3];
rz(-2.5411118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69958413) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(-2.1098095) q[2];
rz(-1.7992841) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.853302) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(1.8348947) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(1.0888938) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13578829) q[0];
sx q[0];
rz(-3.0514768) q[0];
sx q[0];
rz(-2.3401005) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9577272) q[2];
sx q[2];
rz(-0.92391787) q[2];
sx q[2];
rz(1.136214) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2838193) q[1];
sx q[1];
rz(-2.3771264) q[1];
sx q[1];
rz(-2.43999) q[1];
x q[2];
rz(1.9128996) q[3];
sx q[3];
rz(-1.7497239) q[3];
sx q[3];
rz(0.10781328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6178199) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(-0.069996746) q[2];
rz(0.87674117) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-0.35818067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(-1.5169253) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(-2.0786053) q[2];
sx q[2];
rz(-2.4175736) q[2];
sx q[2];
rz(-2.8833817) q[2];
rz(-1.4528965) q[3];
sx q[3];
rz(-2.8040734) q[3];
sx q[3];
rz(0.87344195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
