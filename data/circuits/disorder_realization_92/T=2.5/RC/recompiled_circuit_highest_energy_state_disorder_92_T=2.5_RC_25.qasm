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
rz(0.25843698) q[0];
sx q[0];
rz(-0.28764495) q[0];
sx q[0];
rz(2.0339461) q[0];
rz(-1.8183174) q[1];
sx q[1];
rz(-0.32185093) q[1];
sx q[1];
rz(0.62825656) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2520066) q[0];
sx q[0];
rz(-1.4528251) q[0];
sx q[0];
rz(-1.749176) q[0];
rz(-pi) q[1];
rz(-1.1789178) q[2];
sx q[2];
rz(-0.88681839) q[2];
sx q[2];
rz(-2.1563403) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3648277) q[1];
sx q[1];
rz(-1.7994295) q[1];
sx q[1];
rz(-1.5336114) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1048344) q[3];
sx q[3];
rz(-1.9566571) q[3];
sx q[3];
rz(-0.88632583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2884752) q[2];
sx q[2];
rz(-1.4318117) q[2];
sx q[2];
rz(2.8004069) q[2];
rz(0.8902542) q[3];
sx q[3];
rz(-1.2672) q[3];
sx q[3];
rz(0.53513479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61966908) q[0];
sx q[0];
rz(-0.6837908) q[0];
sx q[0];
rz(0.71440119) q[0];
rz(-1.5419386) q[1];
sx q[1];
rz(-0.49595141) q[1];
sx q[1];
rz(-0.094706789) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9908087) q[0];
sx q[0];
rz(-1.2512733) q[0];
sx q[0];
rz(3.0421769) q[0];
rz(0.33659597) q[2];
sx q[2];
rz(-0.99592745) q[2];
sx q[2];
rz(-2.6624555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75443711) q[1];
sx q[1];
rz(-1.913694) q[1];
sx q[1];
rz(-1.3198402) q[1];
rz(-pi) q[2];
rz(0.50928874) q[3];
sx q[3];
rz(-1.7096161) q[3];
sx q[3];
rz(1.2535439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16227214) q[2];
sx q[2];
rz(-1.8956192) q[2];
sx q[2];
rz(-0.47404131) q[2];
rz(-0.76215172) q[3];
sx q[3];
rz(-0.63665974) q[3];
sx q[3];
rz(-2.6704085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2786461) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(0.85292029) q[0];
rz(0.22311738) q[1];
sx q[1];
rz(-0.47839636) q[1];
sx q[1];
rz(-2.6257637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2775622) q[0];
sx q[0];
rz(-1.5169965) q[0];
sx q[0];
rz(3.0911616) q[0];
rz(-pi) q[1];
rz(2.2081096) q[2];
sx q[2];
rz(-1.3636944) q[2];
sx q[2];
rz(-2.0014971) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3577376) q[1];
sx q[1];
rz(-0.92354316) q[1];
sx q[1];
rz(2.3308862) q[1];
rz(-pi) q[2];
rz(-0.34308691) q[3];
sx q[3];
rz(-0.90215397) q[3];
sx q[3];
rz(0.17681387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3333266) q[2];
sx q[2];
rz(-2.1648679) q[2];
sx q[2];
rz(-0.094956368) q[2];
rz(-0.44148463) q[3];
sx q[3];
rz(-1.35651) q[3];
sx q[3];
rz(-1.9879742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3917291) q[0];
sx q[0];
rz(-0.57498217) q[0];
sx q[0];
rz(2.0854501) q[0];
rz(2.1978343) q[1];
sx q[1];
rz(-2.8902003) q[1];
sx q[1];
rz(-1.4518849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3817134) q[0];
sx q[0];
rz(-2.1739514) q[0];
sx q[0];
rz(-2.2908014) q[0];
rz(-pi) q[1];
rz(-1.9610434) q[2];
sx q[2];
rz(-2.0023291) q[2];
sx q[2];
rz(3.0600037) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6744212) q[1];
sx q[1];
rz(-1.8113664) q[1];
sx q[1];
rz(1.5340541) q[1];
rz(-pi) q[2];
rz(-3.0309276) q[3];
sx q[3];
rz(-1.5783211) q[3];
sx q[3];
rz(1.2570109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4132495) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(0.70017868) q[2];
rz(2.2705196) q[3];
sx q[3];
rz(-2.0756523) q[3];
sx q[3];
rz(-1.7843436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.0102608) q[0];
sx q[0];
rz(-2.6483674) q[0];
sx q[0];
rz(0.48466551) q[0];
rz(-0.87801814) q[1];
sx q[1];
rz(-1.9819825) q[1];
sx q[1];
rz(2.7427618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0244549) q[0];
sx q[0];
rz(-1.0255735) q[0];
sx q[0];
rz(1.412789) q[0];
x q[1];
rz(-2.4381934) q[2];
sx q[2];
rz(-1.3079155) q[2];
sx q[2];
rz(0.92353283) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8139161) q[1];
sx q[1];
rz(-2.3688201) q[1];
sx q[1];
rz(-1.0303968) q[1];
rz(-pi) q[2];
rz(2.5957803) q[3];
sx q[3];
rz(-2.1936676) q[3];
sx q[3];
rz(0.92178492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47348076) q[2];
sx q[2];
rz(-1.5405737) q[2];
sx q[2];
rz(0.87781805) q[2];
rz(2.7569438) q[3];
sx q[3];
rz(-2.2967702) q[3];
sx q[3];
rz(-1.9627242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8574852) q[0];
sx q[0];
rz(-0.50898886) q[0];
sx q[0];
rz(0.2704764) q[0];
rz(-2.8349304) q[1];
sx q[1];
rz(-2.6276734) q[1];
sx q[1];
rz(0.57672966) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78806308) q[0];
sx q[0];
rz(-1.3856141) q[0];
sx q[0];
rz(0.089140541) q[0];
rz(1.1398976) q[2];
sx q[2];
rz(-2.4043519) q[2];
sx q[2];
rz(-2.7708022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9432183) q[1];
sx q[1];
rz(-2.4369708) q[1];
sx q[1];
rz(2.2072029) q[1];
rz(0.0402952) q[3];
sx q[3];
rz(-0.3350122) q[3];
sx q[3];
rz(-0.96714902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3660672) q[2];
sx q[2];
rz(-1.0696573) q[2];
sx q[2];
rz(2.1712187) q[2];
rz(-0.35106418) q[3];
sx q[3];
rz(-2.909436) q[3];
sx q[3];
rz(2.8620913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.905726) q[0];
sx q[0];
rz(-0.38385639) q[0];
sx q[0];
rz(-0.50091499) q[0];
rz(0.72055912) q[1];
sx q[1];
rz(-2.2961398) q[1];
sx q[1];
rz(-2.2038584) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36762992) q[0];
sx q[0];
rz(-1.367462) q[0];
sx q[0];
rz(0.042026566) q[0];
rz(2.4831122) q[2];
sx q[2];
rz(-2.5665433) q[2];
sx q[2];
rz(2.5031896) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5590458) q[1];
sx q[1];
rz(-1.2101047) q[1];
sx q[1];
rz(-2.2466895) q[1];
rz(-pi) q[2];
rz(-0.012091919) q[3];
sx q[3];
rz(-2.4547057) q[3];
sx q[3];
rz(1.505132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0096036) q[2];
sx q[2];
rz(-2.9961573) q[2];
sx q[2];
rz(-2.2930938) q[2];
rz(2.2961473) q[3];
sx q[3];
rz(-0.65631056) q[3];
sx q[3];
rz(1.2257303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7312412) q[0];
sx q[0];
rz(-0.9587962) q[0];
sx q[0];
rz(-0.89286667) q[0];
rz(-3.0800173) q[1];
sx q[1];
rz(-0.67191809) q[1];
sx q[1];
rz(0.16199131) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6380986) q[0];
sx q[0];
rz(-1.5553405) q[0];
sx q[0];
rz(1.4686122) q[0];
rz(2.4396877) q[2];
sx q[2];
rz(-1.6584975) q[2];
sx q[2];
rz(3.1021995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2714241) q[1];
sx q[1];
rz(-1.7080016) q[1];
sx q[1];
rz(0.31183621) q[1];
rz(0.093858899) q[3];
sx q[3];
rz(-1.4673309) q[3];
sx q[3];
rz(-2.2458592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.050345) q[2];
sx q[2];
rz(-0.42751905) q[2];
sx q[2];
rz(0.72489911) q[2];
rz(0.27724087) q[3];
sx q[3];
rz(-0.96479779) q[3];
sx q[3];
rz(3.0550756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8898833) q[0];
sx q[0];
rz(-0.6530264) q[0];
sx q[0];
rz(-3.1157893) q[0];
rz(-1.39894) q[1];
sx q[1];
rz(-1.502123) q[1];
sx q[1];
rz(-3.0339411) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0749507) q[0];
sx q[0];
rz(-0.58569509) q[0];
sx q[0];
rz(-1.8138712) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1321105) q[2];
sx q[2];
rz(-1.138238) q[2];
sx q[2];
rz(-2.3985942) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9840309) q[1];
sx q[1];
rz(-1.5347297) q[1];
sx q[1];
rz(-0.47337516) q[1];
rz(-pi) q[2];
rz(3.1237431) q[3];
sx q[3];
rz(-1.4571725) q[3];
sx q[3];
rz(-2.2872432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4699576) q[2];
sx q[2];
rz(-1.2592955) q[2];
sx q[2];
rz(-2.0476445) q[2];
rz(0.57540244) q[3];
sx q[3];
rz(-2.3016774) q[3];
sx q[3];
rz(-2.0247139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6082918) q[0];
sx q[0];
rz(-2.5550483) q[0];
sx q[0];
rz(0.71304148) q[0];
rz(-2.9167922) q[1];
sx q[1];
rz(-0.59200042) q[1];
sx q[1];
rz(-2.0483268) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21205344) q[0];
sx q[0];
rz(-2.0418379) q[0];
sx q[0];
rz(-0.77917288) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20168882) q[2];
sx q[2];
rz(-1.301577) q[2];
sx q[2];
rz(1.9112916) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0436127) q[1];
sx q[1];
rz(-0.47046767) q[1];
sx q[1];
rz(-2.8941675) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1574) q[3];
sx q[3];
rz(-1.6781312) q[3];
sx q[3];
rz(-1.3232376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0010058086) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(-0.33983964) q[2];
rz(-3.0991683) q[3];
sx q[3];
rz(-1.8570615) q[3];
sx q[3];
rz(2.5137918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18702678) q[0];
sx q[0];
rz(-1.621959) q[0];
sx q[0];
rz(1.8931615) q[0];
rz(-2.3595702) q[1];
sx q[1];
rz(-2.2047058) q[1];
sx q[1];
rz(2.4155736) q[1];
rz(0.37219477) q[2];
sx q[2];
rz(-1.8831461) q[2];
sx q[2];
rz(2.1103721) q[2];
rz(-1.1449161) q[3];
sx q[3];
rz(-1.3419587) q[3];
sx q[3];
rz(1.6423196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
