OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(0.00014076509) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(1.948184) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1228468) q[0];
sx q[0];
rz(-1.7483286) q[0];
sx q[0];
rz(-1.2794504) q[0];
rz(-0.46618669) q[2];
sx q[2];
rz(-0.59980118) q[2];
sx q[2];
rz(0.28238645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27543435) q[1];
sx q[1];
rz(-2.1992116) q[1];
sx q[1];
rz(-2.1584312) q[1];
rz(-pi) q[2];
rz(0.10981202) q[3];
sx q[3];
rz(-1.3545274) q[3];
sx q[3];
rz(-3.0307378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6821735) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(-1.2288644) q[2];
rz(-1.4131644) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035778) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(-1.0128101) q[0];
rz(0.027659841) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(-2.0181296) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2259953) q[0];
sx q[0];
rz(-1.5132656) q[0];
sx q[0];
rz(2.0321839) q[0];
rz(3.0779755) q[2];
sx q[2];
rz(-0.78084313) q[2];
sx q[2];
rz(0.4175248) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9532721) q[1];
sx q[1];
rz(-2.3768432) q[1];
sx q[1];
rz(-2.3482167) q[1];
x q[2];
rz(0.65621891) q[3];
sx q[3];
rz(-1.1074293) q[3];
sx q[3];
rz(-3.0955293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(-2.222555) q[2];
rz(-0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(1.526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(-2.0085874) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5223761) q[0];
sx q[0];
rz(-1.7130865) q[0];
sx q[0];
rz(-3.1382568) q[0];
rz(-pi) q[1];
rz(-2.0793545) q[2];
sx q[2];
rz(-0.8478176) q[2];
sx q[2];
rz(1.522097) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6269835) q[1];
sx q[1];
rz(-1.2830462) q[1];
sx q[1];
rz(2.1327553) q[1];
rz(-2.1902309) q[3];
sx q[3];
rz(-2.1054483) q[3];
sx q[3];
rz(-1.9922647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8901849) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.8481002) q[2];
rz(-0.039316468) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
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
rz(-1.4410102) q[1];
sx q[1];
rz(0.13555759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121658) q[0];
sx q[0];
rz(-2.3608748) q[0];
sx q[0];
rz(2.6902945) q[0];
x q[1];
rz(-1.9549471) q[2];
sx q[2];
rz(-1.640056) q[2];
sx q[2];
rz(-2.4406976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96326522) q[1];
sx q[1];
rz(-1.7566924) q[1];
sx q[1];
rz(-2.9796897) q[1];
rz(-1.2961779) q[3];
sx q[3];
rz(-1.7227968) q[3];
sx q[3];
rz(-2.8677468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(-2.2616852) q[2];
rz(-0.044163477) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(-2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.1039466) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(0.049731072) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(-2.0577046) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3736553) q[0];
sx q[0];
rz(-0.3070139) q[0];
sx q[0];
rz(2.2178749) q[0];
rz(-pi) q[1];
rz(-2.3124144) q[2];
sx q[2];
rz(-0.36311705) q[2];
sx q[2];
rz(-1.6104289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.52884) q[1];
sx q[1];
rz(-1.6788947) q[1];
sx q[1];
rz(1.0428863) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5203939) q[3];
sx q[3];
rz(-1.0574697) q[3];
sx q[3];
rz(0.36171519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9087387) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(0.24442913) q[2];
rz(-2.7092253) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(0.50306815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(1.8571092) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(3.0474512) q[0];
rz(2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(-2.24618) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46195128) q[0];
sx q[0];
rz(-1.2286751) q[0];
sx q[0];
rz(-1.530594) q[0];
rz(-pi) q[1];
rz(0.32161153) q[2];
sx q[2];
rz(-2.4827637) q[2];
sx q[2];
rz(-2.4668601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9858866) q[1];
sx q[1];
rz(-0.17427467) q[1];
sx q[1];
rz(-1.8272912) q[1];
x q[2];
rz(-1.496109) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(-0.96935779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.133693) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(0.80319476) q[2];
rz(1.9512272) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(-2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.068709277) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(-0.58147645) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(1.8849467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0518258) q[0];
sx q[0];
rz(-1.8221812) q[0];
sx q[0];
rz(3.0999523) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9314552) q[2];
sx q[2];
rz(-2.0268831) q[2];
sx q[2];
rz(1.8535341) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9659121) q[1];
sx q[1];
rz(-1.5064081) q[1];
sx q[1];
rz(-2.8596911) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4323746) q[3];
sx q[3];
rz(-1.3243444) q[3];
sx q[3];
rz(-2.646288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1533623) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(-1.3640277) q[2];
rz(-2.2310232) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(-1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0751188) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(2.8334154) q[0];
rz(0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(0.38696188) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9235619) q[0];
sx q[0];
rz(-0.98919981) q[0];
sx q[0];
rz(0.27115718) q[0];
rz(2.2690291) q[2];
sx q[2];
rz(-1.9930895) q[2];
sx q[2];
rz(-2.1625105) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8291694) q[1];
sx q[1];
rz(-2.2409391) q[1];
sx q[1];
rz(-1.6664684) q[1];
x q[2];
rz(0.39051315) q[3];
sx q[3];
rz(-2.6115341) q[3];
sx q[3];
rz(-0.9583677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(2.7015838) q[2];
rz(2.4258339) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(2.0619152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-3.0922906) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54287275) q[0];
sx q[0];
rz(-2.1763986) q[0];
sx q[0];
rz(-2.4332895) q[0];
x q[1];
rz(2.1575035) q[2];
sx q[2];
rz(-0.77312914) q[2];
sx q[2];
rz(-3.0423321) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8763435) q[1];
sx q[1];
rz(-2.0298404) q[1];
sx q[1];
rz(-0.52258073) q[1];
x q[2];
rz(-2.4008972) q[3];
sx q[3];
rz(-1.9441838) q[3];
sx q[3];
rz(1.545056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.07842841) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(2.4196529) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(-0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9938875) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(-2.06185) q[0];
rz(-2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(-1.4019029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7584383) q[0];
sx q[0];
rz(-0.24807319) q[0];
sx q[0];
rz(-1.0323348) q[0];
x q[1];
rz(2.9981668) q[2];
sx q[2];
rz(-2.9529245) q[2];
sx q[2];
rz(2.8667237) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2682174) q[1];
sx q[1];
rz(-2.2443319) q[1];
sx q[1];
rz(0.39336494) q[1];
rz(-pi) q[2];
rz(-1.363477) q[3];
sx q[3];
rz(-0.94982409) q[3];
sx q[3];
rz(-2.6880815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.520291) q[2];
rz(-2.5907717) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(-0.6974535) q[3];
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
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.993492) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(-2.2254754) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(0.73268391) q[2];
sx q[2];
rz(-2.1683243) q[2];
sx q[2];
rz(-1.8666946) q[2];
rz(3.1237596) q[3];
sx q[3];
rz(-2.087839) q[3];
sx q[3];
rz(2.9057333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
