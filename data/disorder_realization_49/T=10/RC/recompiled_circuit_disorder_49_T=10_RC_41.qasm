OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5306659) q[0];
sx q[0];
rz(4.0806169) q[0];
sx q[0];
rz(9.4299849) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(4.2978573) q[1];
sx q[1];
rz(8.2351091) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.386728) q[0];
sx q[0];
rz(-1.0611738) q[0];
sx q[0];
rz(1.1763563) q[0];
rz(-0.66733811) q[2];
sx q[2];
rz(-0.22775209) q[2];
sx q[2];
rz(1.3078794) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3559239) q[1];
sx q[1];
rz(-2.7969116) q[1];
sx q[1];
rz(2.0211401) q[1];
rz(-pi) q[2];
rz(-2.1990843) q[3];
sx q[3];
rz(-1.2357467) q[3];
sx q[3];
rz(3.0818919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.71620119) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(2.544196) q[2];
rz(1.776009) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(-1.8610154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4269203) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(-2.8080217) q[0];
rz(2.0479653) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(-0.11322583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0229867) q[0];
sx q[0];
rz(-1.1456053) q[0];
sx q[0];
rz(-3.0990764) q[0];
x q[1];
rz(3.0683238) q[2];
sx q[2];
rz(-1.5904625) q[2];
sx q[2];
rz(-2.427223) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3791703) q[1];
sx q[1];
rz(-1.1602853) q[1];
sx q[1];
rz(-2.007184) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9085957) q[3];
sx q[3];
rz(-2.5099953) q[3];
sx q[3];
rz(0.88463569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(2.9120581) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(-3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-0.46626058) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(-2.0942028) q[1];
sx q[1];
rz(-1.2068345) q[1];
sx q[1];
rz(0.12761322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022284431) q[0];
sx q[0];
rz(-0.89953178) q[0];
sx q[0];
rz(2.9191201) q[0];
rz(-pi) q[1];
rz(-2.3348546) q[2];
sx q[2];
rz(-0.85959896) q[2];
sx q[2];
rz(-0.15653175) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0652005) q[1];
sx q[1];
rz(-0.25307357) q[1];
sx q[1];
rz(1.2364926) q[1];
rz(-0.25281275) q[3];
sx q[3];
rz(-0.95889927) q[3];
sx q[3];
rz(-1.52724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3601274) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(-0.34960738) q[3];
sx q[3];
rz(-2.4571556) q[3];
sx q[3];
rz(-1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67868245) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(2.983685) q[0];
rz(-2.7754916) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(-0.37240949) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2686553) q[0];
sx q[0];
rz(-0.71632179) q[0];
sx q[0];
rz(-1.3852081) q[0];
rz(-1.0834951) q[2];
sx q[2];
rz(-2.5022025) q[2];
sx q[2];
rz(-1.0207748) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.403703) q[1];
sx q[1];
rz(-0.40778128) q[1];
sx q[1];
rz(1.1853192) q[1];
x q[2];
rz(2.5468772) q[3];
sx q[3];
rz(-2.3200431) q[3];
sx q[3];
rz(2.6298414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41775122) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(0.27238971) q[2];
rz(2.2327936) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(-2.6866384) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(-1.1313261) q[0];
rz(2.5436026) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(-2.4647443) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35447435) q[0];
sx q[0];
rz(-1.1312383) q[0];
sx q[0];
rz(1.4835102) q[0];
rz(-pi) q[1];
rz(-1.3555688) q[2];
sx q[2];
rz(-2.2903246) q[2];
sx q[2];
rz(2.1446251) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9280793) q[1];
sx q[1];
rz(-0.25611862) q[1];
sx q[1];
rz(3.03979) q[1];
x q[2];
rz(2.7480514) q[3];
sx q[3];
rz(-1.593309) q[3];
sx q[3];
rz(-2.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1375492) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(-1.1452902) q[2];
rz(0.14906135) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(-2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1260219) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(1.4591249) q[0];
rz(-0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(-2.2084592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60676735) q[0];
sx q[0];
rz(-1.2762428) q[0];
sx q[0];
rz(1.0876925) q[0];
rz(-pi) q[1];
rz(-1.2007347) q[2];
sx q[2];
rz(-2.1413295) q[2];
sx q[2];
rz(-0.42966118) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3999403) q[1];
sx q[1];
rz(-1.5685023) q[1];
sx q[1];
rz(1.2029349) q[1];
rz(-pi) q[2];
rz(1.798614) q[3];
sx q[3];
rz(-0.67043793) q[3];
sx q[3];
rz(0.05712856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(-2.8549426) q[2];
rz(2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(-2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048112415) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(1.6749143) q[0];
rz(-1.02007) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(2.1405623) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46465835) q[0];
sx q[0];
rz(-1.6081728) q[0];
sx q[0];
rz(1.8514762) q[0];
rz(1.2008576) q[2];
sx q[2];
rz(-1.8897893) q[2];
sx q[2];
rz(1.621304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.24504532) q[1];
sx q[1];
rz(-2.7198615) q[1];
sx q[1];
rz(-1.1827724) q[1];
x q[2];
rz(-1.3592968) q[3];
sx q[3];
rz(-2.2106417) q[3];
sx q[3];
rz(1.1246455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.004185685) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(-3.0380847) q[2];
rz(0.59182709) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(-0.83759585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(-2.5296339) q[0];
rz(1.3423963) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-0.38988316) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1471257) q[0];
sx q[0];
rz(-1.2110706) q[0];
sx q[0];
rz(0.61867449) q[0];
rz(-1.9837603) q[2];
sx q[2];
rz(-1.0997084) q[2];
sx q[2];
rz(-2.4177891) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.08323616) q[1];
sx q[1];
rz(-2.5857158) q[1];
sx q[1];
rz(0.91169375) q[1];
x q[2];
rz(-1.893703) q[3];
sx q[3];
rz(-0.70526037) q[3];
sx q[3];
rz(0.44912072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8395681) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(-0.71425444) q[2];
rz(-0.8977302) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(0.57202655) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(-2.3642448) q[0];
rz(0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(-0.27639595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9696635) q[0];
sx q[0];
rz(-1.164325) q[0];
sx q[0];
rz(-1.7448234) q[0];
rz(-pi) q[1];
rz(-2.4403964) q[2];
sx q[2];
rz(-0.29507911) q[2];
sx q[2];
rz(1.3311177) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.285187) q[1];
sx q[1];
rz(-1.7289843) q[1];
sx q[1];
rz(-1.7941979) q[1];
rz(-pi) q[2];
rz(-0.3474465) q[3];
sx q[3];
rz(-2.5037519) q[3];
sx q[3];
rz(-0.70536246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2892264) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(3.0140871) q[2];
rz(-0.036711983) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7413095) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-2.7888443) q[0];
rz(2.5648975) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(-2.1113077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093826483) q[0];
sx q[0];
rz(-2.2802135) q[0];
sx q[0];
rz(-1.8295656) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2376386) q[2];
sx q[2];
rz(-2.2292456) q[2];
sx q[2];
rz(-1.2308987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87957571) q[1];
sx q[1];
rz(-0.55393065) q[1];
sx q[1];
rz(0.63417706) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2386205) q[3];
sx q[3];
rz(-2.1579451) q[3];
sx q[3];
rz(1.4004933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-0.34118787) q[2];
rz(-0.7406922) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(-3.0497131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13070233) q[0];
sx q[0];
rz(-1.9938835) q[0];
sx q[0];
rz(2.4947517) q[0];
rz(-2.4717992) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(0.43351602) q[2];
sx q[2];
rz(-1.3939861) q[2];
sx q[2];
rz(1.3409333) q[2];
rz(0.68941152) q[3];
sx q[3];
rz(-2.0333615) q[3];
sx q[3];
rz(0.61483308) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
