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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1228468) q[0];
sx q[0];
rz(-1.7483286) q[0];
sx q[0];
rz(1.8621423) q[0];
x q[1];
rz(0.5483746) q[2];
sx q[2];
rz(-1.3142685) q[2];
sx q[2];
rz(0.89474364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27543435) q[1];
sx q[1];
rz(-2.1992116) q[1];
sx q[1];
rz(-0.98316146) q[1];
rz(-pi) q[2];
rz(-0.10981202) q[3];
sx q[3];
rz(-1.3545274) q[3];
sx q[3];
rz(3.0307378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.45941916) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(1.9127282) q[2];
rz(-1.7284283) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(-1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(2.1287825) q[0];
rz(0.027659841) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(2.0181296) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22966188) q[0];
sx q[0];
rz(-2.6768885) q[0];
sx q[0];
rz(-1.4421411) q[0];
x q[1];
rz(-1.5078817) q[2];
sx q[2];
rz(-2.3496369) q[2];
sx q[2];
rz(-0.5069678) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9532721) q[1];
sx q[1];
rz(-2.3768432) q[1];
sx q[1];
rz(-0.79337593) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68617679) q[3];
sx q[3];
rz(-0.78305972) q[3];
sx q[3];
rz(-2.1427597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(-0.91903764) q[2];
rz(-2.4675026) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(-1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27750257) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(1.1330053) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4988574) q[0];
sx q[0];
rz(-2.9992636) q[0];
sx q[0];
rz(1.5475153) q[0];
x q[1];
rz(-0.79046952) q[2];
sx q[2];
rz(-1.1970453) q[2];
sx q[2];
rz(2.8369396) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7743467) q[1];
sx q[1];
rz(-2.5173752) q[1];
sx q[1];
rz(-2.0778076) q[1];
rz(-pi) q[2];
rz(0.95136178) q[3];
sx q[3];
rz(-2.1054483) q[3];
sx q[3];
rz(-1.9922647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(3.1022762) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(1.8815276) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2599729) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(-1.9807293) q[0];
rz(2.2456031) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(3.0060351) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82942688) q[0];
sx q[0];
rz(-2.3608748) q[0];
sx q[0];
rz(2.6902945) q[0];
x q[1];
rz(1.7538257) q[2];
sx q[2];
rz(-0.39003885) q[2];
sx q[2];
rz(-0.70034617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1783274) q[1];
sx q[1];
rz(-1.3849003) q[1];
sx q[1];
rz(-2.9796897) q[1];
rz(-pi) q[2];
rz(-1.0566063) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(-1.7900975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9049412) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039466) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(-1.0132382) q[0];
rz(0.049731072) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(-1.0838881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3736553) q[0];
sx q[0];
rz(-2.8345788) q[0];
sx q[0];
rz(-2.2178749) q[0];
rz(-pi) q[1];
rz(2.3124144) q[2];
sx q[2];
rz(-2.7784756) q[2];
sx q[2];
rz(1.5311637) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0005972) q[1];
sx q[1];
rz(-0.53783572) q[1];
sx q[1];
rz(1.3586033) q[1];
x q[2];
rz(-1.6211987) q[3];
sx q[3];
rz(-1.0574697) q[3];
sx q[3];
rz(0.36171519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9087387) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(-0.24442913) q[2];
rz(0.43236732) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(-2.6385245) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(-0.094141468) q[0];
rz(-2.969818) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(-2.24618) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6796414) q[0];
sx q[0];
rz(-1.2286751) q[0];
sx q[0];
rz(1.6109986) q[0];
rz(-pi) q[1];
rz(-0.32161153) q[2];
sx q[2];
rz(-0.65882896) q[2];
sx q[2];
rz(0.67473251) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9858866) q[1];
sx q[1];
rz(-0.17427467) q[1];
sx q[1];
rz(1.3143015) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6454837) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(-2.1722349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.133693) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(2.3383979) q[2];
rz(1.1903654) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(0.51914006) q[0];
rz(-0.58147645) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(-1.8849467) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.885868) q[0];
sx q[0];
rz(-0.2547383) q[0];
sx q[0];
rz(1.4101009) q[0];
x q[1];
rz(1.1058544) q[2];
sx q[2];
rz(-1.7591811) q[2];
sx q[2];
rz(-0.18907324) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.17649594) q[1];
sx q[1];
rz(-2.8526222) q[1];
sx q[1];
rz(2.9138336) q[1];
x q[2];
rz(1.2506966) q[3];
sx q[3];
rz(-2.2543636) q[3];
sx q[3];
rz(-0.86910955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1533623) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(-1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0664739) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(0.30817729) q[0];
rz(-3.0691052) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(0.38696188) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68650866) q[0];
sx q[0];
rz(-2.5065656) q[0];
sx q[0];
rz(-1.9576661) q[0];
rz(0.87256356) q[2];
sx q[2];
rz(-1.1485032) q[2];
sx q[2];
rz(0.97908212) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8291694) q[1];
sx q[1];
rz(-2.2409391) q[1];
sx q[1];
rz(1.4751242) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.79027) q[3];
sx q[3];
rz(-1.0843127) q[3];
sx q[3];
rz(0.51318491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.62362921) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(2.7015838) q[2];
rz(-2.4258339) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.14116645) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(0.72775841) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(3.0922906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44137529) q[0];
sx q[0];
rz(-2.2451631) q[0];
sx q[0];
rz(-0.81654878) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88843139) q[2];
sx q[2];
rz(-1.9677791) q[2];
sx q[2];
rz(-1.2259442) q[2];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
x q[2];
rz(1.0827761) q[3];
sx q[3];
rz(-2.2501695) q[3];
sx q[3];
rz(-2.8454012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0631642) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(2.4196529) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(0.25434428) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(-1.059277) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(1.4019029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7584383) q[0];
sx q[0];
rz(-2.8935195) q[0];
sx q[0];
rz(2.1092578) q[0];
rz(0.18677588) q[2];
sx q[2];
rz(-1.5439856) q[2];
sx q[2];
rz(1.1550127) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.460443) q[1];
sx q[1];
rz(-0.76421684) q[1];
sx q[1];
rz(-1.1230254) q[1];
rz(-pi) q[2];
rz(-0.28016443) q[3];
sx q[3];
rz(-0.6503085) q[3];
sx q[3];
rz(0.10661099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4828651) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(-1.520291) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(0.6974535) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-2.3475636) q[2];
sx q[2];
rz(-2.2326438) q[2];
sx q[2];
rz(2.2868962) q[2];
rz(-2.0879073) q[3];
sx q[3];
rz(-1.5552945) q[3];
sx q[3];
rz(-1.8154715) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];