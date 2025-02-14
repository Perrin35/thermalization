OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7425793) q[0];
sx q[0];
rz(-1.9785545) q[0];
sx q[0];
rz(1.1385588) q[0];
rz(2.5055655) q[1];
sx q[1];
rz(-2.6316402) q[1];
sx q[1];
rz(-0.62503254) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5020385) q[0];
sx q[0];
rz(-2.9939277) q[0];
sx q[0];
rz(0.30252151) q[0];
rz(1.6051859) q[2];
sx q[2];
rz(-0.54044881) q[2];
sx q[2];
rz(2.2480223) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33746943) q[1];
sx q[1];
rz(-1.9826681) q[1];
sx q[1];
rz(-2.1571674) q[1];
rz(-2.9431683) q[3];
sx q[3];
rz(-2.3188667) q[3];
sx q[3];
rz(-2.2446475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2193489) q[2];
sx q[2];
rz(-0.85180989) q[2];
sx q[2];
rz(-2.7395978) q[2];
rz(2.3123907) q[3];
sx q[3];
rz(-1.821937) q[3];
sx q[3];
rz(1.7792262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7336693) q[0];
sx q[0];
rz(-0.55183721) q[0];
sx q[0];
rz(1.1056939) q[0];
rz(2.3528174) q[1];
sx q[1];
rz(-2.5194247) q[1];
sx q[1];
rz(0.50599352) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0536018) q[0];
sx q[0];
rz(-1.3357541) q[0];
sx q[0];
rz(1.4925692) q[0];
x q[1];
rz(0.77357341) q[2];
sx q[2];
rz(-2.1974652) q[2];
sx q[2];
rz(3.1322014) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11720235) q[1];
sx q[1];
rz(-1.3474696) q[1];
sx q[1];
rz(-1.9322371) q[1];
x q[2];
rz(-2.6384505) q[3];
sx q[3];
rz(-0.50853679) q[3];
sx q[3];
rz(-3.0340305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.904423) q[2];
sx q[2];
rz(-2.0776896) q[2];
sx q[2];
rz(-2.2696631) q[2];
rz(-1.0540086) q[3];
sx q[3];
rz(-2.0619312) q[3];
sx q[3];
rz(-1.4409298) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089040861) q[0];
sx q[0];
rz(-0.65610039) q[0];
sx q[0];
rz(-1.9144527) q[0];
rz(2.6350806) q[1];
sx q[1];
rz(-0.7917234) q[1];
sx q[1];
rz(0.92481771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8313356) q[0];
sx q[0];
rz(-1.7187722) q[0];
sx q[0];
rz(3.0382115) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6396937) q[2];
sx q[2];
rz(-1.1155978) q[2];
sx q[2];
rz(-2.8719939) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23942415) q[1];
sx q[1];
rz(-1.268506) q[1];
sx q[1];
rz(2.2151674) q[1];
rz(0.4731725) q[3];
sx q[3];
rz(-2.4235776) q[3];
sx q[3];
rz(1.1188467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0623124) q[2];
sx q[2];
rz(-1.541905) q[2];
sx q[2];
rz(2.6463032) q[2];
rz(1.7700206) q[3];
sx q[3];
rz(-2.3973231) q[3];
sx q[3];
rz(-1.8194958) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0168734) q[0];
sx q[0];
rz(-1.9833516) q[0];
sx q[0];
rz(-2.25368) q[0];
rz(1.3814231) q[1];
sx q[1];
rz(-2.1235178) q[1];
sx q[1];
rz(-2.6240614) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.237898) q[0];
sx q[0];
rz(-0.26142739) q[0];
sx q[0];
rz(-1.625007) q[0];
rz(0.7184658) q[2];
sx q[2];
rz(-1.9794165) q[2];
sx q[2];
rz(1.1312616) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0486517) q[1];
sx q[1];
rz(-3.042964) q[1];
sx q[1];
rz(1.4709378) q[1];
x q[2];
rz(-1.9512307) q[3];
sx q[3];
rz(-1.1305792) q[3];
sx q[3];
rz(0.39716431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.069328221) q[2];
sx q[2];
rz(-1.1257409) q[2];
sx q[2];
rz(-2.1605055) q[2];
rz(2.6632994) q[3];
sx q[3];
rz(-2.7893453) q[3];
sx q[3];
rz(2.6134764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.9454055) q[0];
sx q[0];
rz(-1.4237175) q[0];
sx q[0];
rz(-3.1086573) q[0];
rz(1.5252652) q[1];
sx q[1];
rz(-2.567629) q[1];
sx q[1];
rz(0.93564916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55452915) q[0];
sx q[0];
rz(-2.0407704) q[0];
sx q[0];
rz(1.5558467) q[0];
rz(-pi) q[1];
rz(-0.79765908) q[2];
sx q[2];
rz(-1.9211384) q[2];
sx q[2];
rz(-1.8961186) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.074240597) q[1];
sx q[1];
rz(-2.3462593) q[1];
sx q[1];
rz(-1.4834845) q[1];
rz(-2.557002) q[3];
sx q[3];
rz(-1.7742312) q[3];
sx q[3];
rz(3.0082321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4719438) q[2];
sx q[2];
rz(-2.7574597) q[2];
sx q[2];
rz(-1.6432537) q[2];
rz(0.60254997) q[3];
sx q[3];
rz(-2.3157401) q[3];
sx q[3];
rz(-1.92441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54685408) q[0];
sx q[0];
rz(-0.29404077) q[0];
sx q[0];
rz(0.038851693) q[0];
rz(2.7815869) q[1];
sx q[1];
rz(-1.729676) q[1];
sx q[1];
rz(-2.9927599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94939828) q[0];
sx q[0];
rz(-1.4440876) q[0];
sx q[0];
rz(1.697798) q[0];
x q[1];
rz(-2.932933) q[2];
sx q[2];
rz(-2.5620915) q[2];
sx q[2];
rz(-0.3267056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.56797261) q[1];
sx q[1];
rz(-2.180033) q[1];
sx q[1];
rz(-0.9602169) q[1];
rz(-0.58601309) q[3];
sx q[3];
rz(-1.8253606) q[3];
sx q[3];
rz(-0.76018047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0845324) q[2];
sx q[2];
rz(-2.3134505) q[2];
sx q[2];
rz(0.35489902) q[2];
rz(2.0830294) q[3];
sx q[3];
rz(-2.1954229) q[3];
sx q[3];
rz(-0.53275776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31838378) q[0];
sx q[0];
rz(-1.4223149) q[0];
sx q[0];
rz(-0.80867714) q[0];
rz(-0.093756229) q[1];
sx q[1];
rz(-2.3010727) q[1];
sx q[1];
rz(2.7309928) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6855934) q[0];
sx q[0];
rz(-2.5846722) q[0];
sx q[0];
rz(0.81057517) q[0];
x q[1];
rz(-2.8062079) q[2];
sx q[2];
rz(-0.68143564) q[2];
sx q[2];
rz(-0.49668703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.73170529) q[1];
sx q[1];
rz(-2.6123939) q[1];
sx q[1];
rz(-0.7284109) q[1];
rz(-pi) q[2];
rz(1.6924573) q[3];
sx q[3];
rz(-0.9100737) q[3];
sx q[3];
rz(-2.4796951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9984596) q[2];
sx q[2];
rz(-0.3930347) q[2];
sx q[2];
rz(3.105799) q[2];
rz(1.834747) q[3];
sx q[3];
rz(-1.1467609) q[3];
sx q[3];
rz(0.80865639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3150113) q[0];
sx q[0];
rz(-2.1666574) q[0];
sx q[0];
rz(-0.63631979) q[0];
rz(2.4782205) q[1];
sx q[1];
rz(-1.1095108) q[1];
sx q[1];
rz(2.5136307) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5945744) q[0];
sx q[0];
rz(-1.5421405) q[0];
sx q[0];
rz(0.076731971) q[0];
rz(0.3971252) q[2];
sx q[2];
rz(-1.4897122) q[2];
sx q[2];
rz(-0.18111595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0481215) q[1];
sx q[1];
rz(-1.9855762) q[1];
sx q[1];
rz(0.61474796) q[1];
rz(-pi) q[2];
rz(-0.18744882) q[3];
sx q[3];
rz(-1.5365684) q[3];
sx q[3];
rz(-1.6899861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5593354) q[2];
sx q[2];
rz(-0.11056837) q[2];
sx q[2];
rz(-1.7250693) q[2];
rz(-2.0420117) q[3];
sx q[3];
rz(-2.1018335) q[3];
sx q[3];
rz(-2.288868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.002554) q[0];
sx q[0];
rz(-2.250493) q[0];
sx q[0];
rz(-0.6148327) q[0];
rz(-0.50827208) q[1];
sx q[1];
rz(-0.95242396) q[1];
sx q[1];
rz(-2.8654548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8042223) q[0];
sx q[0];
rz(-1.7003248) q[0];
sx q[0];
rz(1.999755) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1008087) q[2];
sx q[2];
rz(-1.4333087) q[2];
sx q[2];
rz(-1.9078209) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1663637) q[1];
sx q[1];
rz(-2.1666636) q[1];
sx q[1];
rz(2.1286551) q[1];
rz(-0.34174796) q[3];
sx q[3];
rz(-0.93999388) q[3];
sx q[3];
rz(-0.83021008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.189956) q[2];
sx q[2];
rz(-2.1026244) q[2];
sx q[2];
rz(0.098380066) q[2];
rz(1.9010057) q[3];
sx q[3];
rz(-1.0616579) q[3];
sx q[3];
rz(3.0978751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8561309) q[0];
sx q[0];
rz(-2.0933445) q[0];
sx q[0];
rz(-2.7833126) q[0];
rz(-2.1367392) q[1];
sx q[1];
rz(-0.88468164) q[1];
sx q[1];
rz(0.34214941) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.005363883) q[0];
sx q[0];
rz(-2.507302) q[0];
sx q[0];
rz(2.7691288) q[0];
rz(2.1555485) q[2];
sx q[2];
rz(-1.2504203) q[2];
sx q[2];
rz(-1.4517871) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4163936) q[1];
sx q[1];
rz(-2.1402845) q[1];
sx q[1];
rz(-1.6983022) q[1];
x q[2];
rz(1.9072745) q[3];
sx q[3];
rz(-1.4732756) q[3];
sx q[3];
rz(2.7420126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7973914) q[2];
sx q[2];
rz(-2.6579393) q[2];
sx q[2];
rz(0.3717711) q[2];
rz(0.19112912) q[3];
sx q[3];
rz(-1.976795) q[3];
sx q[3];
rz(0.42310664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38901781) q[0];
sx q[0];
rz(-0.49674635) q[0];
sx q[0];
rz(2.4865271) q[0];
rz(-2.7816506) q[1];
sx q[1];
rz(-1.5725726) q[1];
sx q[1];
rz(-1.6307065) q[1];
rz(-1.7689479) q[2];
sx q[2];
rz(-2.8336278) q[2];
sx q[2];
rz(-2.8358493) q[2];
rz(3.1265518) q[3];
sx q[3];
rz(-2.3287366) q[3];
sx q[3];
rz(-0.47470075) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
