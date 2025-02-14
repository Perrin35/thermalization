OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7944827) q[0];
sx q[0];
rz(-1.4248983) q[0];
sx q[0];
rz(0.27118924) q[0];
rz(3.5085161) q[1];
sx q[1];
rz(2.38382) q[1];
sx q[1];
rz(7.741306) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0163907) q[0];
sx q[0];
rz(-0.0076961829) q[0];
sx q[0];
rz(1.1622692) q[0];
rz(-1.0881937) q[2];
sx q[2];
rz(-1.4652243) q[2];
sx q[2];
rz(-2.2029049) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6071506) q[1];
sx q[1];
rz(-2.1297703) q[1];
sx q[1];
rz(-2.2453151) q[1];
rz(-1.5524158) q[3];
sx q[3];
rz(-1.2871398) q[3];
sx q[3];
rz(-0.90160927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1272507) q[2];
sx q[2];
rz(-2.47561) q[2];
sx q[2];
rz(1.2464397) q[2];
rz(1.0154826) q[3];
sx q[3];
rz(-2.6452439) q[3];
sx q[3];
rz(2.523876) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5366197) q[0];
sx q[0];
rz(-3.0687357) q[0];
sx q[0];
rz(2.4541722) q[0];
rz(-2.5303326) q[1];
sx q[1];
rz(-1.5115073) q[1];
sx q[1];
rz(-2.2329109) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4712418) q[0];
sx q[0];
rz(-2.2178239) q[0];
sx q[0];
rz(0.11423807) q[0];
x q[1];
rz(-2.1852399) q[2];
sx q[2];
rz(-1.4423352) q[2];
sx q[2];
rz(-0.10901448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79230358) q[1];
sx q[1];
rz(-0.97487236) q[1];
sx q[1];
rz(2.2265347) q[1];
x q[2];
rz(-0.55697316) q[3];
sx q[3];
rz(-1.4629629) q[3];
sx q[3];
rz(-1.7086017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0574135) q[2];
sx q[2];
rz(-0.97687352) q[2];
sx q[2];
rz(-1.1581536) q[2];
rz(-0.17364764) q[3];
sx q[3];
rz(-1.4519139) q[3];
sx q[3];
rz(2.4524073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77370206) q[0];
sx q[0];
rz(-2.9359718) q[0];
sx q[0];
rz(-0.50262991) q[0];
rz(1.3357119) q[1];
sx q[1];
rz(-0.10919658) q[1];
sx q[1];
rz(1.4044382) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.805946) q[0];
sx q[0];
rz(-0.61347658) q[0];
sx q[0];
rz(2.4015266) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3777755) q[2];
sx q[2];
rz(-1.2035596) q[2];
sx q[2];
rz(0.45805107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81922779) q[1];
sx q[1];
rz(-2.8016) q[1];
sx q[1];
rz(2.875706) q[1];
rz(-2.5504774) q[3];
sx q[3];
rz(-1.1151259) q[3];
sx q[3];
rz(-2.299813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3001083) q[2];
sx q[2];
rz(-0.69977641) q[2];
sx q[2];
rz(-0.67467275) q[2];
rz(0.35215968) q[3];
sx q[3];
rz(-1.1511185) q[3];
sx q[3];
rz(-1.5153511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34053764) q[0];
sx q[0];
rz(-0.041943701) q[0];
sx q[0];
rz(1.1658143) q[0];
rz(2.5862528) q[1];
sx q[1];
rz(-2.1644939) q[1];
sx q[1];
rz(1.4110483) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3379786) q[0];
sx q[0];
rz(-0.25000152) q[0];
sx q[0];
rz(-3.0743272) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9481692) q[2];
sx q[2];
rz(-2.0900332) q[2];
sx q[2];
rz(-0.80610454) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1852604) q[1];
sx q[1];
rz(-1.0721551) q[1];
sx q[1];
rz(-0.39875984) q[1];
rz(-2.3244546) q[3];
sx q[3];
rz(-1.7445095) q[3];
sx q[3];
rz(-0.44718633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0393151) q[2];
sx q[2];
rz(-0.5216051) q[2];
sx q[2];
rz(-2.6585141) q[2];
rz(-0.77110243) q[3];
sx q[3];
rz(-1.5604138) q[3];
sx q[3];
rz(1.7980827) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1191206) q[0];
sx q[0];
rz(-2.6799057) q[0];
sx q[0];
rz(-0.22892496) q[0];
rz(-1.4643033) q[1];
sx q[1];
rz(-0.56956446) q[1];
sx q[1];
rz(-0.48008188) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10413607) q[0];
sx q[0];
rz(-2.7711282) q[0];
sx q[0];
rz(-2.3568826) q[0];
rz(-1.8304906) q[2];
sx q[2];
rz(-1.1420446) q[2];
sx q[2];
rz(1.4285806) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0220974) q[1];
sx q[1];
rz(-1.9429617) q[1];
sx q[1];
rz(-1.6254025) q[1];
x q[2];
rz(-0.18365209) q[3];
sx q[3];
rz(-2.1878161) q[3];
sx q[3];
rz(-0.69818316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.071712581) q[2];
sx q[2];
rz(-2.0166848) q[2];
sx q[2];
rz(1.9055535) q[2];
rz(-1.49336) q[3];
sx q[3];
rz(-0.83087102) q[3];
sx q[3];
rz(1.4504455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1278095) q[0];
sx q[0];
rz(-0.74087983) q[0];
sx q[0];
rz(2.2770449) q[0];
rz(0.8283444) q[1];
sx q[1];
rz(-1.336785) q[1];
sx q[1];
rz(2.4116662) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975271) q[0];
sx q[0];
rz(-1.1008917) q[0];
sx q[0];
rz(0.92961981) q[0];
rz(-pi) q[1];
rz(1.3719158) q[2];
sx q[2];
rz(-2.2305397) q[2];
sx q[2];
rz(-2.4937862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.266842) q[1];
sx q[1];
rz(-1.0462985) q[1];
sx q[1];
rz(-0.28171087) q[1];
rz(-pi) q[2];
rz(0.72336332) q[3];
sx q[3];
rz(-1.5302684) q[3];
sx q[3];
rz(-1.2637637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.571903) q[2];
sx q[2];
rz(-2.9986585) q[2];
sx q[2];
rz(1.8611056) q[2];
rz(1.5661543) q[3];
sx q[3];
rz(-2.1011293) q[3];
sx q[3];
rz(0.86281002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60523072) q[0];
sx q[0];
rz(-2.1385758) q[0];
sx q[0];
rz(-2.5501116) q[0];
rz(-0.65762562) q[1];
sx q[1];
rz(-1.769519) q[1];
sx q[1];
rz(-0.11810158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0869136) q[0];
sx q[0];
rz(-0.42092338) q[0];
sx q[0];
rz(-2.551159) q[0];
rz(-pi) q[1];
rz(-3.0131091) q[2];
sx q[2];
rz(-2.621935) q[2];
sx q[2];
rz(2.3574587) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1527295) q[1];
sx q[1];
rz(-1.0514469) q[1];
sx q[1];
rz(1.4593655) q[1];
rz(-pi) q[2];
rz(0.92175092) q[3];
sx q[3];
rz(-2.8650682) q[3];
sx q[3];
rz(-0.10504237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9401231) q[2];
sx q[2];
rz(-0.38429364) q[2];
sx q[2];
rz(-1.2598134) q[2];
rz(-1.2093557) q[3];
sx q[3];
rz(-1.7989379) q[3];
sx q[3];
rz(2.2945837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2131293) q[0];
sx q[0];
rz(-2.6452112) q[0];
sx q[0];
rz(-1.5877566) q[0];
rz(-1.3262879) q[1];
sx q[1];
rz(-2.1554048) q[1];
sx q[1];
rz(2.3415668) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2934809) q[0];
sx q[0];
rz(-1.1901374) q[0];
sx q[0];
rz(-1.8649376) q[0];
x q[1];
rz(-2.0532908) q[2];
sx q[2];
rz(-2.2721599) q[2];
sx q[2];
rz(0.33368233) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86944888) q[1];
sx q[1];
rz(-0.71681685) q[1];
sx q[1];
rz(0.79811704) q[1];
rz(-pi) q[2];
rz(-2.127384) q[3];
sx q[3];
rz(-1.2383023) q[3];
sx q[3];
rz(1.0300762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.942261) q[2];
sx q[2];
rz(-2.3776725) q[2];
sx q[2];
rz(-1.2660816) q[2];
rz(2.074504) q[3];
sx q[3];
rz(-1.4184003) q[3];
sx q[3];
rz(3.0838695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1570194) q[0];
sx q[0];
rz(-0.83681256) q[0];
sx q[0];
rz(0.19503221) q[0];
rz(-2.2795279) q[1];
sx q[1];
rz(-1.74086) q[1];
sx q[1];
rz(-0.21924266) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040960366) q[0];
sx q[0];
rz(-1.0182683) q[0];
sx q[0];
rz(-0.29968963) q[0];
x q[1];
rz(1.931483) q[2];
sx q[2];
rz(-2.6668913) q[2];
sx q[2];
rz(2.0550967) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91867204) q[1];
sx q[1];
rz(-1.6038451) q[1];
sx q[1];
rz(-0.53302427) q[1];
rz(-1.0107993) q[3];
sx q[3];
rz(-0.83314291) q[3];
sx q[3];
rz(-3.0411947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7361043) q[2];
sx q[2];
rz(-2.2443503) q[2];
sx q[2];
rz(-1.415095) q[2];
rz(1.3380346) q[3];
sx q[3];
rz(-1.3058563) q[3];
sx q[3];
rz(-3.0726748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59704429) q[0];
sx q[0];
rz(-2.2291849) q[0];
sx q[0];
rz(1.9507116) q[0];
rz(1.4471588) q[1];
sx q[1];
rz(-1.2971327) q[1];
sx q[1];
rz(1.4318633) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3668684) q[0];
sx q[0];
rz(-0.69126832) q[0];
sx q[0];
rz(0.33095215) q[0];
rz(-pi) q[1];
rz(-1.9831311) q[2];
sx q[2];
rz(-1.025505) q[2];
sx q[2];
rz(-2.3067428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9976342) q[1];
sx q[1];
rz(-1.9737509) q[1];
sx q[1];
rz(-0.87586276) q[1];
rz(3.1095554) q[3];
sx q[3];
rz(-0.83870974) q[3];
sx q[3];
rz(-0.17994954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4900964) q[2];
sx q[2];
rz(-0.33612529) q[2];
sx q[2];
rz(-0.87745848) q[2];
rz(1.6089926) q[3];
sx q[3];
rz(-1.420615) q[3];
sx q[3];
rz(-2.1144313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.132591) q[0];
sx q[0];
rz(-1.7283716) q[0];
sx q[0];
rz(-0.51474095) q[0];
rz(0.71350907) q[1];
sx q[1];
rz(-2.3206354) q[1];
sx q[1];
rz(-1.6101507) q[1];
rz(0.3729214) q[2];
sx q[2];
rz(-1.3013617) q[2];
sx q[2];
rz(-0.19993275) q[2];
rz(1.7969098) q[3];
sx q[3];
rz(-1.3158847) q[3];
sx q[3];
rz(-2.55008) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
