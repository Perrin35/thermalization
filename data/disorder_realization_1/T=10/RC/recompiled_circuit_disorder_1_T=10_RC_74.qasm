OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(5.1685652) q[0];
sx q[0];
rz(9.4246372) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(-2.1773832) q[1];
sx q[1];
rz(1.1934086) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084273987) q[0];
sx q[0];
rz(-2.8017375) q[0];
sx q[0];
rz(2.1291332) q[0];
rz(-pi) q[1];
rz(2.675406) q[2];
sx q[2];
rz(-2.5417915) q[2];
sx q[2];
rz(2.8592062) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5692917) q[1];
sx q[1];
rz(-0.83218677) q[1];
sx q[1];
rz(2.4898847) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7883349) q[3];
sx q[3];
rz(-1.4635524) q[3];
sx q[3];
rz(-1.4835964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6821735) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(1.9127282) q[2];
rz(-1.4131644) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035778) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(2.1287825) q[0];
rz(-0.027659841) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-2.0181296) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2259953) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(2.0321839) q[0];
rz(-pi) q[1];
rz(-0.063617184) q[2];
sx q[2];
rz(-2.3607495) q[2];
sx q[2];
rz(-0.4175248) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23501913) q[1];
sx q[1];
rz(-1.0636914) q[1];
sx q[1];
rz(-2.1706235) q[1];
rz(-pi) q[2];
rz(-0.68617679) q[3];
sx q[3];
rz(-2.3585329) q[3];
sx q[3];
rz(-0.99883294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(2.222555) q[2];
rz(2.4675026) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(-1.8664237) q[0];
rz(0.69349849) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(1.1330053) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61921652) q[0];
sx q[0];
rz(-1.4285061) q[0];
sx q[0];
rz(-3.1382568) q[0];
rz(-pi) q[1];
rz(2.6373367) q[2];
sx q[2];
rz(-0.85668737) q[2];
sx q[2];
rz(0.91932636) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7743467) q[1];
sx q[1];
rz(-0.62421747) q[1];
sx q[1];
rz(1.063785) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62882256) q[3];
sx q[3];
rz(-2.0938794) q[3];
sx q[3];
rz(-0.76997013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2514078) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(3.1022762) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(-1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2599729) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(-1.1608634) q[0];
rz(-2.2456031) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(3.0060351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82942688) q[0];
sx q[0];
rz(-2.3608748) q[0];
sx q[0];
rz(0.45129816) q[0];
x q[1];
rz(-0.074684871) q[2];
sx q[2];
rz(-1.1876145) q[2];
sx q[2];
rz(2.243724) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5038824) q[1];
sx q[1];
rz(-1.4117068) q[1];
sx q[1];
rz(1.3824944) q[1];
rz(-pi) q[2];
rz(-1.8454148) q[3];
sx q[3];
rz(-1.4187958) q[3];
sx q[3];
rz(-2.8677468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23665145) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(-0.87990749) q[2];
rz(-0.044163477) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(0.28863171) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0376461) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(-3.0918616) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(-2.0577046) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.097682) q[0];
sx q[0];
rz(-1.327276) q[0];
sx q[0];
rz(2.9527412) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2976228) q[2];
sx q[2];
rz(-1.8130842) q[2];
sx q[2];
rz(-2.4730686) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0005972) q[1];
sx q[1];
rz(-0.53783572) q[1];
sx q[1];
rz(-1.7829893) q[1];
x q[2];
rz(-0.51387042) q[3];
sx q[3];
rz(-1.5268945) q[3];
sx q[3];
rz(-1.9572789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9087387) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(-0.24442913) q[2];
rz(2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2844834) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(0.094141468) q[0];
rz(-2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(2.24618) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46195128) q[0];
sx q[0];
rz(-1.2286751) q[0];
sx q[0];
rz(1.6109986) q[0];
rz(-0.32161153) q[2];
sx q[2];
rz(-2.4827637) q[2];
sx q[2];
rz(2.4668601) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15570607) q[1];
sx q[1];
rz(-2.967318) q[1];
sx q[1];
rz(1.8272912) q[1];
rz(1.6454837) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(-0.96935779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.133693) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(1.9512272) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(-0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068709277) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(-2.6224526) q[0];
rz(0.58147645) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(-1.8849467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089766895) q[0];
sx q[0];
rz(-1.3194114) q[0];
sx q[0];
rz(0.041640394) q[0];
x q[1];
rz(2.9314552) q[2];
sx q[2];
rz(-1.1147095) q[2];
sx q[2];
rz(1.2880585) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1756806) q[1];
sx q[1];
rz(-1.6351846) q[1];
sx q[1];
rz(2.8596911) q[1];
x q[2];
rz(-0.7092181) q[3];
sx q[3];
rz(-1.3243444) q[3];
sx q[3];
rz(0.49530464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1533623) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(-1.3640277) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(-1.6114657) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0664739) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(2.8334154) q[0];
rz(-3.0691052) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(0.38696188) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21803074) q[0];
sx q[0];
rz(-0.98919981) q[0];
sx q[0];
rz(-2.8704355) q[0];
rz(-pi) q[1];
rz(2.6110821) q[2];
sx q[2];
rz(-0.94420099) q[2];
sx q[2];
rz(0.26041398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-0.67589251) q[1];
sx q[1];
rz(-3.0216316) q[1];
rz(0.39051315) q[3];
sx q[3];
rz(-2.6115341) q[3];
sx q[3];
rz(-0.9583677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62362921) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(-2.7015838) q[2];
rz(2.4258339) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(2.0619152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(-1.0986885) q[0];
rz(2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(0.049302014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5987199) q[0];
sx q[0];
rz(-2.1763986) q[0];
sx q[0];
rz(2.4332895) q[0];
rz(-pi) q[1];
rz(2.2531613) q[2];
sx q[2];
rz(-1.1738136) q[2];
sx q[2];
rz(-1.9156485) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2652492) q[1];
sx q[1];
rz(-2.0298404) q[1];
sx q[1];
rz(0.52258073) q[1];
rz(-2.4008972) q[3];
sx q[3];
rz(-1.9441838) q[3];
sx q[3];
rz(1.545056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.07842841) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(-2.4196529) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9938875) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(-1.0797427) q[0];
rz(-1.059277) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(-1.4019029) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3831543) q[0];
sx q[0];
rz(-0.24807319) q[0];
sx q[0];
rz(-2.1092578) q[0];
rz(-0.18677588) q[2];
sx q[2];
rz(-1.5439856) q[2];
sx q[2];
rz(-1.1550127) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.460443) q[1];
sx q[1];
rz(-0.76421684) q[1];
sx q[1];
rz(2.0185673) q[1];
rz(-pi) q[2];
rz(-1.7781156) q[3];
sx q[3];
rz(-0.94982409) q[3];
sx q[3];
rz(-0.45351115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(1.6213017) q[2];
rz(2.5907717) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(-2.4441392) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.993492) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(-0.91611721) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(-2.4089087) q[2];
sx q[2];
rz(-2.1683243) q[2];
sx q[2];
rz(-1.8666946) q[2];
rz(1.5394474) q[3];
sx q[3];
rz(-0.51732224) q[3];
sx q[3];
rz(2.8696685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];