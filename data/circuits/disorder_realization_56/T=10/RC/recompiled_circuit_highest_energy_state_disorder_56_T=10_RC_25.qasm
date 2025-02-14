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
rz(2.7409878) q[0];
sx q[0];
rz(-0.40979835) q[0];
sx q[0];
rz(-2.0706489) q[0];
rz(-2.5228956) q[1];
sx q[1];
rz(-0.56801152) q[1];
sx q[1];
rz(-2.2979589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0742357) q[0];
sx q[0];
rz(-1.5412742) q[0];
sx q[0];
rz(-1.3062472) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29542342) q[2];
sx q[2];
rz(-1.43338) q[2];
sx q[2];
rz(-0.58770056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2045793) q[1];
sx q[1];
rz(-2.3790303) q[1];
sx q[1];
rz(2.049905) q[1];
rz(-pi) q[2];
rz(-1.0759962) q[3];
sx q[3];
rz(-1.1462948) q[3];
sx q[3];
rz(-2.692846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3598651) q[2];
sx q[2];
rz(-0.83911506) q[2];
sx q[2];
rz(-0.016999379) q[2];
rz(1.4012236) q[3];
sx q[3];
rz(-2.5798116) q[3];
sx q[3];
rz(1.9248272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725937) q[0];
sx q[0];
rz(-1.8091135) q[0];
sx q[0];
rz(0.78729415) q[0];
rz(0.17838082) q[1];
sx q[1];
rz(-1.9352103) q[1];
sx q[1];
rz(1.9040727) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4263692) q[0];
sx q[0];
rz(-0.58642095) q[0];
sx q[0];
rz(0.1118025) q[0];
rz(-pi) q[1];
rz(2.8493153) q[2];
sx q[2];
rz(-0.9773796) q[2];
sx q[2];
rz(3.127272) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.552265) q[1];
sx q[1];
rz(-2.5922212) q[1];
sx q[1];
rz(2.2087847) q[1];
rz(-2.7010553) q[3];
sx q[3];
rz(-1.8228662) q[3];
sx q[3];
rz(-2.1217176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42147192) q[2];
sx q[2];
rz(-1.6933491) q[2];
sx q[2];
rz(0.063701542) q[2];
rz(1.2740159) q[3];
sx q[3];
rz(-2.2985022) q[3];
sx q[3];
rz(1.564285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0434697) q[0];
sx q[0];
rz(-1.948057) q[0];
sx q[0];
rz(0.7435588) q[0];
rz(-2.6877563) q[1];
sx q[1];
rz(-1.791714) q[1];
sx q[1];
rz(-2.7226864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1322238) q[0];
sx q[0];
rz(-1.1282578) q[0];
sx q[0];
rz(2.924463) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9100283) q[2];
sx q[2];
rz(-0.76398173) q[2];
sx q[2];
rz(2.6126249) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6082904) q[1];
sx q[1];
rz(-2.4455297) q[1];
sx q[1];
rz(-1.1185557) q[1];
rz(-pi) q[2];
rz(2.2499372) q[3];
sx q[3];
rz(-1.4640995) q[3];
sx q[3];
rz(-1.5556435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0539187) q[2];
sx q[2];
rz(-0.81834617) q[2];
sx q[2];
rz(-1.9107001) q[2];
rz(2.6026717) q[3];
sx q[3];
rz(-1.6659707) q[3];
sx q[3];
rz(1.6787136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26381275) q[0];
sx q[0];
rz(-2.3872264) q[0];
sx q[0];
rz(-2.5394649) q[0];
rz(1.4471794) q[1];
sx q[1];
rz(-0.68232957) q[1];
sx q[1];
rz(0.86300659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2857638) q[0];
sx q[0];
rz(-1.4572954) q[0];
sx q[0];
rz(0.72300006) q[0];
x q[1];
rz(-0.90410467) q[2];
sx q[2];
rz(-0.46510425) q[2];
sx q[2];
rz(0.99921528) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.132024) q[1];
sx q[1];
rz(-1.3088262) q[1];
sx q[1];
rz(0.97673194) q[1];
x q[2];
rz(0.3163655) q[3];
sx q[3];
rz(-0.85815198) q[3];
sx q[3];
rz(-2.0848839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.065923125) q[2];
sx q[2];
rz(-1.1752335) q[2];
sx q[2];
rz(-2.2389331) q[2];
rz(-2.3265808) q[3];
sx q[3];
rz(-2.5383526) q[3];
sx q[3];
rz(-2.7284315) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0937423) q[0];
sx q[0];
rz(-0.047690064) q[0];
sx q[0];
rz(-0.93511859) q[0];
rz(-1.6085666) q[1];
sx q[1];
rz(-2.8159499) q[1];
sx q[1];
rz(-2.8757222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2926798) q[0];
sx q[0];
rz(-0.32912185) q[0];
sx q[0];
rz(1.2364682) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8498672) q[2];
sx q[2];
rz(-2.3014567) q[2];
sx q[2];
rz(-1.0026635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2468202) q[1];
sx q[1];
rz(-2.7547464) q[1];
sx q[1];
rz(0.43718613) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10281296) q[3];
sx q[3];
rz(-1.307337) q[3];
sx q[3];
rz(2.3578532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68044668) q[2];
sx q[2];
rz(-1.5004044) q[2];
sx q[2];
rz(-2.6908596) q[2];
rz(-1.9510673) q[3];
sx q[3];
rz(-0.7014941) q[3];
sx q[3];
rz(-2.7019971) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615237) q[0];
sx q[0];
rz(-2.838205) q[0];
sx q[0];
rz(-0.044483749) q[0];
rz(-0.26875177) q[1];
sx q[1];
rz(-0.93695295) q[1];
sx q[1];
rz(-0.43089795) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4558894) q[0];
sx q[0];
rz(-0.95818943) q[0];
sx q[0];
rz(-1.2489178) q[0];
rz(-pi) q[1];
rz(0.90164945) q[2];
sx q[2];
rz(-1.2001654) q[2];
sx q[2];
rz(1.2610051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2984276) q[1];
sx q[1];
rz(-2.3379605) q[1];
sx q[1];
rz(2.1433565) q[1];
rz(2.481302) q[3];
sx q[3];
rz(-1.5662282) q[3];
sx q[3];
rz(2.3167603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1162794) q[2];
sx q[2];
rz(-0.37386027) q[2];
sx q[2];
rz(-2.5171793) q[2];
rz(1.1011018) q[3];
sx q[3];
rz(-0.81239429) q[3];
sx q[3];
rz(-2.1541434) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91200149) q[0];
sx q[0];
rz(-1.0593375) q[0];
sx q[0];
rz(-0.68369317) q[0];
rz(1.4423485) q[1];
sx q[1];
rz(-2.3577299) q[1];
sx q[1];
rz(2.9396465) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1037052) q[0];
sx q[0];
rz(-1.5695086) q[0];
sx q[0];
rz(1.2426069) q[0];
x q[1];
rz(1.6228637) q[2];
sx q[2];
rz(-0.74801842) q[2];
sx q[2];
rz(-1.6287273) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8117042) q[1];
sx q[1];
rz(-1.0502019) q[1];
sx q[1];
rz(1.4938526) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45777623) q[3];
sx q[3];
rz(-1.6094064) q[3];
sx q[3];
rz(0.82070551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.086229) q[2];
sx q[2];
rz(-1.7217041) q[2];
sx q[2];
rz(1.8602547) q[2];
rz(-2.4751439) q[3];
sx q[3];
rz(-2.0446348) q[3];
sx q[3];
rz(-0.54653978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9796889) q[0];
sx q[0];
rz(-2.5894916) q[0];
sx q[0];
rz(2.6276278) q[0];
rz(-2.1886096) q[1];
sx q[1];
rz(-2.516808) q[1];
sx q[1];
rz(-2.0228588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6020929) q[0];
sx q[0];
rz(-1.148466) q[0];
sx q[0];
rz(-1.4813478) q[0];
x q[1];
rz(-2.0055112) q[2];
sx q[2];
rz(-0.82821199) q[2];
sx q[2];
rz(-1.0880053) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7756164) q[1];
sx q[1];
rz(-1.1630285) q[1];
sx q[1];
rz(-2.9080176) q[1];
rz(-1.9268613) q[3];
sx q[3];
rz(-2.0907932) q[3];
sx q[3];
rz(2.4442152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3861367) q[2];
sx q[2];
rz(-0.35217199) q[2];
sx q[2];
rz(2.238671) q[2];
rz(-1.4593982) q[3];
sx q[3];
rz(-1.7995588) q[3];
sx q[3];
rz(-2.3193147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4626386) q[0];
sx q[0];
rz(-0.32005388) q[0];
sx q[0];
rz(1.1902887) q[0];
rz(-2.8207488) q[1];
sx q[1];
rz(-1.4535934) q[1];
sx q[1];
rz(-1.920059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71305481) q[0];
sx q[0];
rz(-2.894141) q[0];
sx q[0];
rz(-0.58899458) q[0];
x q[1];
rz(1.2668816) q[2];
sx q[2];
rz(-0.7579782) q[2];
sx q[2];
rz(2.9044819) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9170984) q[1];
sx q[1];
rz(-0.82991822) q[1];
sx q[1];
rz(1.8486449) q[1];
x q[2];
rz(1.3876347) q[3];
sx q[3];
rz(-1.7622928) q[3];
sx q[3];
rz(0.93246704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2373576) q[2];
sx q[2];
rz(-2.9866437) q[2];
sx q[2];
rz(-2.8083189) q[2];
rz(-1.0171558) q[3];
sx q[3];
rz(-2.1867496) q[3];
sx q[3];
rz(2.8216383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79621133) q[0];
sx q[0];
rz(-2.6295202) q[0];
sx q[0];
rz(2.2228125) q[0];
rz(-2.9504919) q[1];
sx q[1];
rz(-0.42593503) q[1];
sx q[1];
rz(-2.3474615) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0427754) q[0];
sx q[0];
rz(-2.2198027) q[0];
sx q[0];
rz(-2.2521583) q[0];
rz(-2.7808519) q[2];
sx q[2];
rz(-1.3731925) q[2];
sx q[2];
rz(0.2195356) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0865143) q[1];
sx q[1];
rz(-2.3740413) q[1];
sx q[1];
rz(-1.1805736) q[1];
rz(-pi) q[2];
rz(1.3598241) q[3];
sx q[3];
rz(-1.8856241) q[3];
sx q[3];
rz(-1.4538159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0309151) q[2];
sx q[2];
rz(-0.6157178) q[2];
sx q[2];
rz(-0.75882971) q[2];
rz(-0.87016726) q[3];
sx q[3];
rz(-2.3535959) q[3];
sx q[3];
rz(-1.0677392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0017241521) q[0];
sx q[0];
rz(-1.4586466) q[0];
sx q[0];
rz(1.9139105) q[0];
rz(-2.7036746) q[1];
sx q[1];
rz(-1.4358078) q[1];
sx q[1];
rz(-1.5461071) q[1];
rz(-2.6545637) q[2];
sx q[2];
rz(-1.9430046) q[2];
sx q[2];
rz(-0.76443048) q[2];
rz(2.6043456) q[3];
sx q[3];
rz(-2.2840847) q[3];
sx q[3];
rz(1.7084697) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
