OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61092678) q[0];
sx q[0];
rz(-0.93902421) q[0];
sx q[0];
rz(3.1363857) q[0];
rz(-1.3448673) q[1];
sx q[1];
rz(-1.1562647) q[1];
sx q[1];
rz(1.1896689) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0931041) q[0];
sx q[0];
rz(-2.5079873) q[0];
sx q[0];
rz(2.5392169) q[0];
rz(-pi) q[1];
rz(2.9615133) q[2];
sx q[2];
rz(-1.7109949) q[2];
sx q[2];
rz(2.7498498) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7835044) q[1];
sx q[1];
rz(-1.4231829) q[1];
sx q[1];
rz(-1.8834157) q[1];
x q[2];
rz(-2.1990843) q[3];
sx q[3];
rz(-1.9058459) q[3];
sx q[3];
rz(0.059700746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4253915) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(2.544196) q[2];
rz(-1.776009) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4269203) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(-0.33357099) q[0];
rz(2.0479653) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(-0.11322583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2213759) q[0];
sx q[0];
rz(-2.7144103) q[0];
sx q[0];
rz(1.4772052) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.073268854) q[2];
sx q[2];
rz(-1.5511302) q[2];
sx q[2];
rz(2.427223) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37564056) q[1];
sx q[1];
rz(-1.1728219) q[1];
sx q[1];
rz(2.6938733) q[1];
x q[2];
rz(2.9037335) q[3];
sx q[3];
rz(-0.97994643) q[3];
sx q[3];
rz(1.8464586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(-0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46626058) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(-2.0942028) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(-0.12761322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.688296) q[0];
sx q[0];
rz(-1.3971546) q[0];
sx q[0];
rz(-2.2542473) q[0];
rz(-2.4650376) q[2];
sx q[2];
rz(-2.1495719) q[2];
sx q[2];
rz(2.3253331) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.7319273) q[1];
sx q[1];
rz(-1.3320142) q[1];
sx q[1];
rz(0.084652918) q[1];
rz(-pi) q[2];
rz(-2.8887799) q[3];
sx q[3];
rz(-0.95889927) q[3];
sx q[3];
rz(-1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3601274) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-2.8313417) q[2];
rz(-2.7919853) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(-1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.4629102) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(-0.15790766) q[0];
rz(-0.36610106) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(-2.7691832) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87293738) q[0];
sx q[0];
rz(-0.71632179) q[0];
sx q[0];
rz(1.3852081) q[0];
rz(-pi) q[1];
rz(-2.0580975) q[2];
sx q[2];
rz(-2.5022025) q[2];
sx q[2];
rz(-2.1208178) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6652674) q[1];
sx q[1];
rz(-1.421126) q[1];
sx q[1];
rz(1.1900351) q[1];
rz(-pi) q[2];
rz(-2.5468772) q[3];
sx q[3];
rz(-2.3200431) q[3];
sx q[3];
rz(0.51175129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7238414) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(2.8692029) q[2];
rz(-2.2327936) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(0.4549543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.882778) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(-2.0102665) q[0];
rz(2.5436026) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(2.4647443) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35447435) q[0];
sx q[0];
rz(-2.0103543) q[0];
sx q[0];
rz(1.6580824) q[0];
rz(-pi) q[1];
rz(0.73111515) q[2];
sx q[2];
rz(-1.4094681) q[2];
sx q[2];
rz(0.71691712) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9280793) q[1];
sx q[1];
rz(-0.25611862) q[1];
sx q[1];
rz(-3.03979) q[1];
x q[2];
rz(-0.39354126) q[3];
sx q[3];
rz(-1.5482836) q[3];
sx q[3];
rz(2.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1375492) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(-1.9963025) q[2];
rz(2.9925313) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(1.0173652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260219) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(-1.4591249) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(2.2084592) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60676735) q[0];
sx q[0];
rz(-1.8653499) q[0];
sx q[0];
rz(2.0539001) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.538751) q[2];
sx q[2];
rz(-1.8800929) q[2];
sx q[2];
rz(1.7939292) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3999403) q[1];
sx q[1];
rz(-1.5730904) q[1];
sx q[1];
rz(1.2029349) q[1];
rz(2.9643781) q[3];
sx q[3];
rz(-0.92068499) q[3];
sx q[3];
rz(2.9110416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3559945) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(0.28665001) q[2];
rz(2.8921228) q[3];
sx q[3];
rz(-1.9727861) q[3];
sx q[3];
rz(-0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048112415) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(-1.4666784) q[0];
rz(-2.1215227) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(2.1405623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46465835) q[0];
sx q[0];
rz(-1.6081728) q[0];
sx q[0];
rz(1.8514762) q[0];
rz(-pi) q[1];
rz(2.8011488) q[2];
sx q[2];
rz(-1.2203487) q[2];
sx q[2];
rz(-0.0705138) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.66623464) q[1];
sx q[1];
rz(-1.1821786) q[1];
sx q[1];
rz(-2.9734441) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27490297) q[3];
sx q[3];
rz(-0.66920815) q[3];
sx q[3];
rz(1.4698524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.137407) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(-0.10350791) q[2];
rz(-2.5497656) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(-2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25933927) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(-0.6119588) q[0];
rz(-1.3423963) q[1];
sx q[1];
rz(-1.9806769) q[1];
sx q[1];
rz(2.7517095) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1471257) q[0];
sx q[0];
rz(-1.2110706) q[0];
sx q[0];
rz(2.5229182) q[0];
rz(-1.9837603) q[2];
sx q[2];
rz(-1.0997084) q[2];
sx q[2];
rz(-2.4177891) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65615678) q[1];
sx q[1];
rz(-2.0011144) q[1];
sx q[1];
rz(0.3635316) q[1];
rz(-pi) q[2];
x q[2];
rz(2.249986) q[3];
sx q[3];
rz(-1.3636175) q[3];
sx q[3];
rz(0.87219119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8395681) q[2];
sx q[2];
rz(-0.53047696) q[2];
sx q[2];
rz(2.4273382) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(2.5695661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(-0.77734787) q[0];
rz(0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(-0.27639595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9696635) q[0];
sx q[0];
rz(-1.9772677) q[0];
sx q[0];
rz(-1.7448234) q[0];
rz(0.22819613) q[2];
sx q[2];
rz(-1.3820634) q[2];
sx q[2];
rz(-0.91918321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2497219) q[1];
sx q[1];
rz(-2.8686214) q[1];
sx q[1];
rz(-2.194838) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3235839) q[3];
sx q[3];
rz(-2.1650378) q[3];
sx q[3];
rz(1.1288527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(3.0140871) q[2];
rz(0.036711983) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(-1.9071074) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7413095) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(0.35274831) q[0];
rz(2.5648975) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(2.1113077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.661783) q[0];
sx q[0];
rz(-0.74736887) q[0];
sx q[0];
rz(-2.8519147) q[0];
x q[1];
rz(1.2376386) q[2];
sx q[2];
rz(-2.2292456) q[2];
sx q[2];
rz(1.2308987) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5926338) q[1];
sx q[1];
rz(-1.1332129) q[1];
sx q[1];
rz(-1.9220819) q[1];
rz(-0.61334765) q[3];
sx q[3];
rz(-1.8457335) q[3];
sx q[3];
rz(0.018523286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(-0.7406922) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(-3.0497131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0108903) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(-2.4717992) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(0.40217051) q[2];
sx q[2];
rz(-0.46605863) q[2];
sx q[2];
rz(-3.0083187) q[2];
rz(-2.4767247) q[3];
sx q[3];
rz(-2.3330199) q[3];
sx q[3];
rz(-1.4521269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
