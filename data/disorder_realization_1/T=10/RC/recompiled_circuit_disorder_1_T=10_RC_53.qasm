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
rz(-1.1146201) q[0];
sx q[0];
rz(-0.00014076509) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(-2.1773832) q[1];
sx q[1];
rz(1.1934086) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0187459) q[0];
sx q[0];
rz(-1.393264) q[0];
sx q[0];
rz(-1.8621423) q[0];
rz(-pi) q[1];
rz(-0.46618669) q[2];
sx q[2];
rz(-0.59980118) q[2];
sx q[2];
rz(-2.8592062) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4729893) q[1];
sx q[1];
rz(-2.0358634) q[1];
sx q[1];
rz(-0.71778645) q[1];
x q[2];
rz(-0.10981202) q[3];
sx q[3];
rz(-1.7870652) q[3];
sx q[3];
rz(0.11085489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45941916) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(-1.9127282) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035778) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(1.0128101) q[0];
rz(-0.027659841) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(1.123463) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37338049) q[0];
sx q[0];
rz(-1.1102311) q[0];
sx q[0];
rz(0.064231355) q[0];
rz(-1.6337109) q[2];
sx q[2];
rz(-2.3496369) q[2];
sx q[2];
rz(-2.6346249) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1883205) q[1];
sx q[1];
rz(-2.3768432) q[1];
sx q[1];
rz(-2.3482167) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[1];
rz(0.79364395) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(0.91903764) q[2];
rz(-0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27750257) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.2751689) q[0];
rz(-0.69349849) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(2.0085874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95110675) q[0];
sx q[0];
rz(-1.5740984) q[0];
sx q[0];
rz(1.4285054) q[0];
x q[1];
rz(0.50425597) q[2];
sx q[2];
rz(-2.2849053) q[2];
sx q[2];
rz(0.91932636) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9085711) q[1];
sx q[1];
rz(-2.1070707) q[1];
sx q[1];
rz(-0.33645333) q[1];
x q[2];
rz(-0.62882256) q[3];
sx q[3];
rz(-2.0938794) q[3];
sx q[3];
rz(-2.3716225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(-1.8481002) q[2];
rz(-3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816198) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(-1.9807293) q[0];
rz(-0.89598957) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(-3.0060351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82942688) q[0];
sx q[0];
rz(-0.78071785) q[0];
sx q[0];
rz(-0.45129816) q[0];
x q[1];
rz(-1.3877669) q[2];
sx q[2];
rz(-2.7515538) q[2];
sx q[2];
rz(0.70034617) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5038824) q[1];
sx q[1];
rz(-1.4117068) q[1];
sx q[1];
rz(-1.7590982) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0566063) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(1.3514951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23665145) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(-0.87990749) q[2];
rz(-0.044163477) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(-0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0376461) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(1.0132382) q[0];
rz(3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(-2.0577046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5686544) q[0];
sx q[0];
rz(-1.7540115) q[0];
sx q[0];
rz(-1.818548) q[0];
rz(2.8903557) q[2];
sx q[2];
rz(-1.8357956) q[2];
sx q[2];
rz(-0.83515177) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.14099546) q[1];
sx q[1];
rz(-2.6037569) q[1];
sx q[1];
rz(1.3586033) q[1];
rz(-pi) q[2];
rz(3.0524592) q[3];
sx q[3];
rz(-0.51557487) q[3];
sx q[3];
rz(2.6775132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9087387) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(2.8971635) q[2];
rz(-0.43236732) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(-0.50306815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(3.0474512) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(-0.89541268) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6796414) q[0];
sx q[0];
rz(-1.9129176) q[0];
sx q[0];
rz(-1.530594) q[0];
rz(-0.63352872) q[2];
sx q[2];
rz(-1.3760566) q[2];
sx q[2];
rz(-0.63846904) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9858866) q[1];
sx q[1];
rz(-0.17427467) q[1];
sx q[1];
rz(1.8272912) q[1];
rz(-pi) q[2];
rz(1.6454837) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(-0.96935779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.133693) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(-1.1903654) q[3];
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
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068709277) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(-0.51914006) q[0];
rz(2.5601162) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(1.2566459) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089766895) q[0];
sx q[0];
rz(-1.8221812) q[0];
sx q[0];
rz(-3.0999523) q[0];
rz(-pi) q[1];
rz(-1.16876) q[2];
sx q[2];
rz(-2.6425344) q[2];
sx q[2];
rz(-1.7390342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9659121) q[1];
sx q[1];
rz(-1.6351846) q[1];
sx q[1];
rz(-0.28190159) q[1];
rz(2.7729633) q[3];
sx q[3];
rz(-0.7437403) q[3];
sx q[3];
rz(-1.3524692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98823035) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(1.3640277) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(-1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0751188) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(0.30817729) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-1.0132353) q[1];
sx q[1];
rz(-2.7546308) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2012126) q[0];
sx q[0];
rz(-1.3450755) q[0];
sx q[0];
rz(2.1696027) q[0];
rz(-pi) q[1];
rz(-0.87256356) q[2];
sx q[2];
rz(-1.1485032) q[2];
sx q[2];
rz(2.1625105) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-0.67589251) q[1];
sx q[1];
rz(0.11996108) q[1];
x q[2];
rz(-2.7510795) q[3];
sx q[3];
rz(-2.6115341) q[3];
sx q[3];
rz(2.183225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62362921) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(0.44000885) q[2];
rz(2.4258339) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-2.0619152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14116645) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(-1.0986885) q[0];
rz(-0.72775841) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(-0.049302014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6599777) q[0];
sx q[0];
rz(-2.1351295) q[0];
sx q[0];
rz(0.83156395) q[0];
rz(2.2531613) q[2];
sx q[2];
rz(-1.1738136) q[2];
sx q[2];
rz(-1.9156485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2652492) q[1];
sx q[1];
rz(-1.1117522) q[1];
sx q[1];
rz(0.52258073) q[1];
x q[2];
rz(-2.6155869) q[3];
sx q[3];
rz(-2.3283539) q[3];
sx q[3];
rz(2.7362652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0631642) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(2.4196529) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(0.25434428) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14770517) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(1.0797427) q[0];
rz(1.059277) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(-1.4019029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3831543) q[0];
sx q[0];
rz(-2.8935195) q[0];
sx q[0];
rz(1.0323348) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18677588) q[2];
sx q[2];
rz(-1.597607) q[2];
sx q[2];
rz(-1.9865799) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68114963) q[1];
sx q[1];
rz(-2.3773758) q[1];
sx q[1];
rz(1.1230254) q[1];
rz(-pi) q[2];
rz(-1.7781156) q[3];
sx q[3];
rz(-2.1917686) q[3];
sx q[3];
rz(0.45351115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.520291) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(-2.4441392) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14810066) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(-0.91611721) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(0.79402906) q[2];
sx q[2];
rz(-2.2326438) q[2];
sx q[2];
rz(2.2868962) q[2];
rz(-0.017833088) q[3];
sx q[3];
rz(-2.087839) q[3];
sx q[3];
rz(2.9057333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
