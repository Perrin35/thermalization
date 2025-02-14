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
rz(1.5647178) q[0];
sx q[0];
rz(-1.325542) q[0];
sx q[0];
rz(2.5817885) q[0];
rz(1.3979823) q[1];
sx q[1];
rz(-1.8140732) q[1];
sx q[1];
rz(-1.7550069) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5965393) q[0];
sx q[0];
rz(-0.60205108) q[0];
sx q[0];
rz(-1.1132147) q[0];
x q[1];
rz(-1.6835911) q[2];
sx q[2];
rz(-1.2769645) q[2];
sx q[2];
rz(-0.71982671) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.69405729) q[1];
sx q[1];
rz(-0.22429286) q[1];
sx q[1];
rz(-0.38557012) q[1];
rz(-pi) q[2];
rz(-1.8614426) q[3];
sx q[3];
rz(-1.564102) q[3];
sx q[3];
rz(0.98095241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6036512) q[2];
sx q[2];
rz(-0.62749046) q[2];
sx q[2];
rz(-2.6998935) q[2];
rz(-0.80080992) q[3];
sx q[3];
rz(-1.8947314) q[3];
sx q[3];
rz(-0.39113623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9073198) q[0];
sx q[0];
rz(-1.7664302) q[0];
sx q[0];
rz(-0.32670879) q[0];
rz(-1.3455343) q[1];
sx q[1];
rz(-1.5252472) q[1];
sx q[1];
rz(-0.84652841) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8252944) q[0];
sx q[0];
rz(-0.40028737) q[0];
sx q[0];
rz(1.2918703) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2791762) q[2];
sx q[2];
rz(-0.90718944) q[2];
sx q[2];
rz(0.070726591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0051028) q[1];
sx q[1];
rz(-1.041689) q[1];
sx q[1];
rz(-0.012261475) q[1];
x q[2];
rz(1.2760824) q[3];
sx q[3];
rz(-2.3008153) q[3];
sx q[3];
rz(-1.4286014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9091866) q[2];
sx q[2];
rz(-2.0441983) q[2];
sx q[2];
rz(0.48420134) q[2];
rz(-1.8353315) q[3];
sx q[3];
rz(-2.5808915) q[3];
sx q[3];
rz(1.9561214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66121286) q[0];
sx q[0];
rz(-0.73037195) q[0];
sx q[0];
rz(-2.6015306) q[0];
rz(-1.9909319) q[1];
sx q[1];
rz(-1.2723424) q[1];
sx q[1];
rz(2.5475492) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4580989) q[0];
sx q[0];
rz(-0.20266315) q[0];
sx q[0];
rz(0.99838241) q[0];
rz(-0.39806367) q[2];
sx q[2];
rz(-0.69483611) q[2];
sx q[2];
rz(0.24206012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6642528) q[1];
sx q[1];
rz(-2.3788733) q[1];
sx q[1];
rz(1.4422756) q[1];
x q[2];
rz(1.698231) q[3];
sx q[3];
rz(-2.6563247) q[3];
sx q[3];
rz(-1.3290153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8192886) q[2];
sx q[2];
rz(-1.015859) q[2];
sx q[2];
rz(2.4508396) q[2];
rz(-2.6639719) q[3];
sx q[3];
rz(-0.4126637) q[3];
sx q[3];
rz(-0.42948183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8884884) q[0];
sx q[0];
rz(-0.34585837) q[0];
sx q[0];
rz(2.389287) q[0];
rz(2.2350156) q[1];
sx q[1];
rz(-0.38175672) q[1];
sx q[1];
rz(0.22901542) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9045712) q[0];
sx q[0];
rz(-1.8638049) q[0];
sx q[0];
rz(-1.3364393) q[0];
x q[1];
rz(-2.7732753) q[2];
sx q[2];
rz(-2.5715156) q[2];
sx q[2];
rz(-1.3454448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.70826958) q[1];
sx q[1];
rz(-0.3618985) q[1];
sx q[1];
rz(-0.72558484) q[1];
rz(-pi) q[2];
rz(0.58521478) q[3];
sx q[3];
rz(-1.2622331) q[3];
sx q[3];
rz(1.2206961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0166246) q[2];
sx q[2];
rz(-0.11398537) q[2];
sx q[2];
rz(2.1086741) q[2];
rz(2.0756663) q[3];
sx q[3];
rz(-1.4706127) q[3];
sx q[3];
rz(1.8377931) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2579987) q[0];
sx q[0];
rz(-3.0553387) q[0];
sx q[0];
rz(-1.3294543) q[0];
rz(-0.52779245) q[1];
sx q[1];
rz(-1.8254435) q[1];
sx q[1];
rz(-2.7425308) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2891727) q[0];
sx q[0];
rz(-1.4932844) q[0];
sx q[0];
rz(1.9231251) q[0];
rz(-pi) q[1];
x q[1];
rz(2.40995) q[2];
sx q[2];
rz(-1.2564617) q[2];
sx q[2];
rz(0.86282496) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0001155) q[1];
sx q[1];
rz(-0.85073227) q[1];
sx q[1];
rz(-1.1517555) q[1];
rz(-1.3436704) q[3];
sx q[3];
rz(-2.4339614) q[3];
sx q[3];
rz(-0.16498241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.993678) q[2];
sx q[2];
rz(-2.3892011) q[2];
sx q[2];
rz(-2.1368775) q[2];
rz(0.29804722) q[3];
sx q[3];
rz(-1.3532676) q[3];
sx q[3];
rz(2.0731488) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3598969) q[0];
sx q[0];
rz(-1.5157461) q[0];
sx q[0];
rz(1.9500649) q[0];
rz(-0.53939348) q[1];
sx q[1];
rz(-2.0627779) q[1];
sx q[1];
rz(-0.76621145) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3193885) q[0];
sx q[0];
rz(-2.4766716) q[0];
sx q[0];
rz(0.98916905) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3503803) q[2];
sx q[2];
rz(-2.3177766) q[2];
sx q[2];
rz(-2.8132766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6886544) q[1];
sx q[1];
rz(-2.3804745) q[1];
sx q[1];
rz(-1.2378986) q[1];
rz(-pi) q[2];
rz(1.8155072) q[3];
sx q[3];
rz(-2.7608216) q[3];
sx q[3];
rz(2.0920193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.811565) q[2];
sx q[2];
rz(-1.4469688) q[2];
sx q[2];
rz(-0.88328973) q[2];
rz(-0.53089321) q[3];
sx q[3];
rz(-2.5575432) q[3];
sx q[3];
rz(2.365621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74060488) q[0];
sx q[0];
rz(-2.8477836) q[0];
sx q[0];
rz(2.4921932) q[0];
rz(-0.3736639) q[1];
sx q[1];
rz(-2.3039736) q[1];
sx q[1];
rz(1.1423473) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23655791) q[0];
sx q[0];
rz(-1.8982186) q[0];
sx q[0];
rz(-0.16752807) q[0];
rz(0.20988237) q[2];
sx q[2];
rz(-1.58605) q[2];
sx q[2];
rz(1.588087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2616407) q[1];
sx q[1];
rz(-1.2625719) q[1];
sx q[1];
rz(0.76517268) q[1];
x q[2];
rz(-3.0965831) q[3];
sx q[3];
rz(-1.112072) q[3];
sx q[3];
rz(-2.5436202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.00804) q[2];
sx q[2];
rz(-0.69591659) q[2];
sx q[2];
rz(2.7827061) q[2];
rz(-2.533203) q[3];
sx q[3];
rz(-1.5865934) q[3];
sx q[3];
rz(-0.86148328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1818039) q[0];
sx q[0];
rz(-0.53878468) q[0];
sx q[0];
rz(2.1754919) q[0];
rz(0.44202647) q[1];
sx q[1];
rz(-1.288488) q[1];
sx q[1];
rz(-0.40837049) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7743384) q[0];
sx q[0];
rz(-1.5975195) q[0];
sx q[0];
rz(-1.7033808) q[0];
rz(-pi) q[1];
rz(1.6138861) q[2];
sx q[2];
rz(-1.4231963) q[2];
sx q[2];
rz(-1.4606295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8345771) q[1];
sx q[1];
rz(-0.6041827) q[1];
sx q[1];
rz(-1.8240364) q[1];
rz(1.3134052) q[3];
sx q[3];
rz(-0.85581778) q[3];
sx q[3];
rz(-1.7325213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24892204) q[2];
sx q[2];
rz(-0.20888027) q[2];
sx q[2];
rz(-1.9425617) q[2];
rz(-1.8438953) q[3];
sx q[3];
rz(-2.0632931) q[3];
sx q[3];
rz(-3.1395932) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93160981) q[0];
sx q[0];
rz(-1.4298121) q[0];
sx q[0];
rz(-2.0344875) q[0];
rz(-0.18868407) q[1];
sx q[1];
rz(-0.37982267) q[1];
sx q[1];
rz(1.3245378) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077515006) q[0];
sx q[0];
rz(-2.0821838) q[0];
sx q[0];
rz(1.2426503) q[0];
rz(0.021090551) q[2];
sx q[2];
rz(-0.97375768) q[2];
sx q[2];
rz(1.8151547) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75211084) q[1];
sx q[1];
rz(-0.5403203) q[1];
sx q[1];
rz(1.575927) q[1];
rz(-pi) q[2];
rz(-2.0641493) q[3];
sx q[3];
rz(-2.2098118) q[3];
sx q[3];
rz(0.30917707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.79335) q[2];
sx q[2];
rz(-2.5620354) q[2];
sx q[2];
rz(-2.9938475) q[2];
rz(0.82792264) q[3];
sx q[3];
rz(-1.7235651) q[3];
sx q[3];
rz(2.2016321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4350236) q[0];
sx q[0];
rz(-0.35025418) q[0];
sx q[0];
rz(-1.5237923) q[0];
rz(-1.891547) q[1];
sx q[1];
rz(-2.3077272) q[1];
sx q[1];
rz(0.42125431) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6235879) q[0];
sx q[0];
rz(-2.5941018) q[0];
sx q[0];
rz(-2.3836552) q[0];
rz(0.6217072) q[2];
sx q[2];
rz(-1.6424456) q[2];
sx q[2];
rz(3.0376584) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6727461) q[1];
sx q[1];
rz(-1.4839659) q[1];
sx q[1];
rz(1.561291) q[1];
x q[2];
rz(1.3280844) q[3];
sx q[3];
rz(-1.1303217) q[3];
sx q[3];
rz(0.8410359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0201575) q[2];
sx q[2];
rz(-0.73548663) q[2];
sx q[2];
rz(0.45041931) q[2];
rz(1.2823229) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(2.569017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12238518) q[0];
sx q[0];
rz(-1.6078124) q[0];
sx q[0];
rz(-1.5851198) q[0];
rz(2.4284594) q[1];
sx q[1];
rz(-1.5819555) q[1];
sx q[1];
rz(-3.1176007) q[1];
rz(0.85999388) q[2];
sx q[2];
rz(-2.4888629) q[2];
sx q[2];
rz(0.98301907) q[2];
rz(1.4982669) q[3];
sx q[3];
rz(-1.4909496) q[3];
sx q[3];
rz(-0.95076233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
