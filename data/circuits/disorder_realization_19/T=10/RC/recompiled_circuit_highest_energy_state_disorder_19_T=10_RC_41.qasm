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
rz(0.28111464) q[0];
sx q[0];
rz(-1.3507564) q[0];
sx q[0];
rz(-0.95066345) q[0];
rz(3.364346) q[1];
sx q[1];
rz(2.8877701) q[1];
sx q[1];
rz(8.756898) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0216103) q[0];
sx q[0];
rz(-1.9304781) q[0];
sx q[0];
rz(1.5403454) q[0];
x q[1];
rz(-3.1257258) q[2];
sx q[2];
rz(-1.7363225) q[2];
sx q[2];
rz(0.5918146) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3006134) q[1];
sx q[1];
rz(-1.6924227) q[1];
sx q[1];
rz(-0.089861265) q[1];
x q[2];
rz(-1.5667772) q[3];
sx q[3];
rz(-1.7815456) q[3];
sx q[3];
rz(0.8523418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5369947) q[2];
sx q[2];
rz(-2.1945685) q[2];
sx q[2];
rz(0.14744082) q[2];
rz(1.3747181) q[3];
sx q[3];
rz(-1.7539897) q[3];
sx q[3];
rz(-1.714777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0680256) q[0];
sx q[0];
rz(-1.0140714) q[0];
sx q[0];
rz(-0.2036988) q[0];
rz(3.0286466) q[1];
sx q[1];
rz(-2.4599383) q[1];
sx q[1];
rz(3.122094) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72806304) q[0];
sx q[0];
rz(-2.3747281) q[0];
sx q[0];
rz(-1.9770115) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2157529) q[2];
sx q[2];
rz(-0.40242919) q[2];
sx q[2];
rz(0.42237626) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9894101) q[1];
sx q[1];
rz(-1.0355897) q[1];
sx q[1];
rz(1.2104976) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56657378) q[3];
sx q[3];
rz(-2.1687897) q[3];
sx q[3];
rz(-1.2815042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.514275) q[2];
sx q[2];
rz(-1.5655727) q[2];
sx q[2];
rz(-0.19228284) q[2];
rz(-0.2187885) q[3];
sx q[3];
rz(-1.3008806) q[3];
sx q[3];
rz(-1.3505664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4718276) q[0];
sx q[0];
rz(-2.8596523) q[0];
sx q[0];
rz(-0.48926735) q[0];
rz(-1.6911814) q[1];
sx q[1];
rz(-1.9905118) q[1];
sx q[1];
rz(-2.6077008) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.691026) q[0];
sx q[0];
rz(-1.2493142) q[0];
sx q[0];
rz(1.5984999) q[0];
rz(-pi) q[1];
rz(3.0025173) q[2];
sx q[2];
rz(-1.5267045) q[2];
sx q[2];
rz(1.9382375) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3855672) q[1];
sx q[1];
rz(-1.7122762) q[1];
sx q[1];
rz(-1.1668851) q[1];
rz(-pi) q[2];
rz(2.4160014) q[3];
sx q[3];
rz(-1.707336) q[3];
sx q[3];
rz(1.2927593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78407946) q[2];
sx q[2];
rz(-2.1046905) q[2];
sx q[2];
rz(-0.64884031) q[2];
rz(0.31323788) q[3];
sx q[3];
rz(-2.1514838) q[3];
sx q[3];
rz(1.311897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.4637777) q[0];
sx q[0];
rz(-2.2721993) q[0];
sx q[0];
rz(-0.75611269) q[0];
rz(-0.4568049) q[1];
sx q[1];
rz(-2.017338) q[1];
sx q[1];
rz(0.66064984) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27460262) q[0];
sx q[0];
rz(-1.9714234) q[0];
sx q[0];
rz(1.3763877) q[0];
rz(1.0370937) q[2];
sx q[2];
rz(-0.42345475) q[2];
sx q[2];
rz(-1.6771984) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.403549) q[1];
sx q[1];
rz(-1.3894086) q[1];
sx q[1];
rz(0.84457835) q[1];
x q[2];
rz(1.1180498) q[3];
sx q[3];
rz(-1.2635316) q[3];
sx q[3];
rz(-1.6045517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75055355) q[2];
sx q[2];
rz(-1.4918574) q[2];
sx q[2];
rz(2.1569815) q[2];
rz(1.5084958) q[3];
sx q[3];
rz(-2.0506141) q[3];
sx q[3];
rz(0.35198894) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19531798) q[0];
sx q[0];
rz(-1.0349422) q[0];
sx q[0];
rz(0.6629194) q[0];
rz(-1.7341057) q[1];
sx q[1];
rz(-0.88277849) q[1];
sx q[1];
rz(2.1972806) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1809826) q[0];
sx q[0];
rz(-0.91413218) q[0];
sx q[0];
rz(1.3707536) q[0];
x q[1];
rz(2.3572982) q[2];
sx q[2];
rz(-0.81740618) q[2];
sx q[2];
rz(-2.9767175) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.014537) q[1];
sx q[1];
rz(-2.1528798) q[1];
sx q[1];
rz(2.0392557) q[1];
x q[2];
rz(-1.3317033) q[3];
sx q[3];
rz(-1.0464365) q[3];
sx q[3];
rz(-0.6839377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0974836) q[2];
sx q[2];
rz(-0.58610836) q[2];
sx q[2];
rz(-0.07494542) q[2];
rz(-0.99503851) q[3];
sx q[3];
rz(-1.1040265) q[3];
sx q[3];
rz(1.9869355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5465882) q[0];
sx q[0];
rz(-2.0363448) q[0];
sx q[0];
rz(-0.30584359) q[0];
rz(-0.88286895) q[1];
sx q[1];
rz(-0.63778937) q[1];
sx q[1];
rz(-0.84552228) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3061241) q[0];
sx q[0];
rz(-1.7996664) q[0];
sx q[0];
rz(2.2717975) q[0];
rz(-1.8550625) q[2];
sx q[2];
rz(-1.5696811) q[2];
sx q[2];
rz(-1.0363611) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9915139) q[1];
sx q[1];
rz(-2.2631915) q[1];
sx q[1];
rz(0.69025363) q[1];
rz(0.85725089) q[3];
sx q[3];
rz(-1.7117867) q[3];
sx q[3];
rz(0.14234358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4800097) q[2];
sx q[2];
rz(-2.2955387) q[2];
sx q[2];
rz(-1.0711077) q[2];
rz(2.8969911) q[3];
sx q[3];
rz(-0.94541234) q[3];
sx q[3];
rz(1.5379813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5474608) q[0];
sx q[0];
rz(-0.56532156) q[0];
sx q[0];
rz(0.93455899) q[0];
rz(2.4681828) q[1];
sx q[1];
rz(-0.62973657) q[1];
sx q[1];
rz(-0.096079439) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5058511) q[0];
sx q[0];
rz(-1.4560501) q[0];
sx q[0];
rz(3.1155354) q[0];
rz(1.7504018) q[2];
sx q[2];
rz(-1.7471004) q[2];
sx q[2];
rz(1.4063032) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9597577) q[1];
sx q[1];
rz(-1.3980306) q[1];
sx q[1];
rz(-0.29675014) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0453048) q[3];
sx q[3];
rz(-0.91595338) q[3];
sx q[3];
rz(1.8055467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38976321) q[2];
sx q[2];
rz(-1.7182173) q[2];
sx q[2];
rz(2.5982889) q[2];
rz(-0.027103847) q[3];
sx q[3];
rz(-0.47309858) q[3];
sx q[3];
rz(2.134034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.84697023) q[0];
sx q[0];
rz(-2.4610418) q[0];
sx q[0];
rz(-0.74817014) q[0];
rz(-1.6698042) q[1];
sx q[1];
rz(-1.4133778) q[1];
sx q[1];
rz(-0.69563037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3601937) q[0];
sx q[0];
rz(-1.3142141) q[0];
sx q[0];
rz(0.11959038) q[0];
rz(0.28098051) q[2];
sx q[2];
rz(-2.4644445) q[2];
sx q[2];
rz(0.80981648) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7456949) q[1];
sx q[1];
rz(-1.1198178) q[1];
sx q[1];
rz(-1.1956426) q[1];
rz(-pi) q[2];
rz(-3.031894) q[3];
sx q[3];
rz(-1.6792751) q[3];
sx q[3];
rz(1.2622076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2601605) q[2];
sx q[2];
rz(-0.83923927) q[2];
sx q[2];
rz(0.053675573) q[2];
rz(1.7565049) q[3];
sx q[3];
rz(-1.3421007) q[3];
sx q[3];
rz(-2.9746941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78228918) q[0];
sx q[0];
rz(-1.8817236) q[0];
sx q[0];
rz(2.2480929) q[0];
rz(2.9529849) q[1];
sx q[1];
rz(-0.74917561) q[1];
sx q[1];
rz(-2.9395054) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30464393) q[0];
sx q[0];
rz(-2.3173387) q[0];
sx q[0];
rz(0.040764256) q[0];
rz(-pi) q[1];
x q[1];
rz(0.049742266) q[2];
sx q[2];
rz(-1.047985) q[2];
sx q[2];
rz(-0.78429121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1043865) q[1];
sx q[1];
rz(-1.361711) q[1];
sx q[1];
rz(-1.0250183) q[1];
rz(-pi) q[2];
rz(-1.8340183) q[3];
sx q[3];
rz(-1.211245) q[3];
sx q[3];
rz(1.5375801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9933652) q[2];
sx q[2];
rz(-1.8348285) q[2];
sx q[2];
rz(-2.4723049) q[2];
rz(0.78684849) q[3];
sx q[3];
rz(-1.9653178) q[3];
sx q[3];
rz(-2.440786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65449077) q[0];
sx q[0];
rz(-2.632532) q[0];
sx q[0];
rz(-1.2855726) q[0];
rz(-2.9945943) q[1];
sx q[1];
rz(-1.9963943) q[1];
sx q[1];
rz(2.1910892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72785801) q[0];
sx q[0];
rz(-1.04299) q[0];
sx q[0];
rz(1.0038654) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6754402) q[2];
sx q[2];
rz(-0.54000914) q[2];
sx q[2];
rz(-1.4252942) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.353668) q[1];
sx q[1];
rz(-1.0312925) q[1];
sx q[1];
rz(-2.2108393) q[1];
rz(-pi) q[2];
rz(1.688429) q[3];
sx q[3];
rz(-1.7087987) q[3];
sx q[3];
rz(1.7774297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6805083) q[2];
sx q[2];
rz(-2.4775439) q[2];
sx q[2];
rz(-0.39592478) q[2];
rz(-0.33216533) q[3];
sx q[3];
rz(-1.9931404) q[3];
sx q[3];
rz(0.90126669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0300776) q[0];
sx q[0];
rz(-2.5319396) q[0];
sx q[0];
rz(1.4421705) q[0];
rz(3.1051927) q[1];
sx q[1];
rz(-1.2926688) q[1];
sx q[1];
rz(-1.7784437) q[1];
rz(-1.2716952) q[2];
sx q[2];
rz(-0.80812412) q[2];
sx q[2];
rz(-3.0199188) q[2];
rz(0.66988173) q[3];
sx q[3];
rz(-1.5276147) q[3];
sx q[3];
rz(2.0338175) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
