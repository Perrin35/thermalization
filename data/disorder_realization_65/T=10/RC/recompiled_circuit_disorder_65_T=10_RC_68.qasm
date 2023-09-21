OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(-1.4927827) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(-1.0746497) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91905347) q[0];
sx q[0];
rz(-3.1226282) q[0];
sx q[0];
rz(2.8595964) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9825748) q[2];
sx q[2];
rz(-2.3454925) q[2];
sx q[2];
rz(-0.12479347) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5584471) q[1];
sx q[1];
rz(-1.4902643) q[1];
sx q[1];
rz(3.0196378) q[1];
rz(2.0243458) q[3];
sx q[3];
rz(-1.082886) q[3];
sx q[3];
rz(-0.32353668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(1.1734022) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56869498) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(-1.6404023) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06022913) q[0];
sx q[0];
rz(-1.4481359) q[0];
sx q[0];
rz(1.5457904) q[0];
x q[1];
rz(-2.8295849) q[2];
sx q[2];
rz(-1.1216251) q[2];
sx q[2];
rz(-2.1829407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3850296) q[1];
sx q[1];
rz(-2.2379413) q[1];
sx q[1];
rz(-1.7900311) q[1];
rz(-1.1797446) q[3];
sx q[3];
rz(-1.565233) q[3];
sx q[3];
rz(0.55913505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(1.6274874) q[2];
rz(1.9747915) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(-0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0062155) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(2.0626383) q[0];
rz(-1.9619933) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(1.5037781) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62961783) q[0];
sx q[0];
rz(-0.10911988) q[0];
sx q[0];
rz(2.2114121) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2321646) q[2];
sx q[2];
rz(-0.69790188) q[2];
sx q[2];
rz(2.1991889) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40239247) q[1];
sx q[1];
rz(-1.270073) q[1];
sx q[1];
rz(-1.6154358) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4012666) q[3];
sx q[3];
rz(-0.66344122) q[3];
sx q[3];
rz(2.0142114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7109795) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(2.3977996) q[2];
rz(-0.25742325) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(0.58810365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-1.051735) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.1490885) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3455032) q[0];
sx q[0];
rz(-0.26608276) q[0];
sx q[0];
rz(0.43568133) q[0];
rz(0.43196584) q[2];
sx q[2];
rz(-1.0216121) q[2];
sx q[2];
rz(0.84749046) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9149575) q[1];
sx q[1];
rz(-1.2775704) q[1];
sx q[1];
rz(-2.3073767) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37255128) q[3];
sx q[3];
rz(-2.3636892) q[3];
sx q[3];
rz(-2.3140964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(2.4728298) q[2];
rz(1.8956005) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(-0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4887061) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(-2.0892129) q[0];
rz(-1.2106238) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(-2.2573684) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8540871) q[0];
sx q[0];
rz(-0.26212087) q[0];
sx q[0];
rz(0.8192807) q[0];
x q[1];
rz(-0.40712936) q[2];
sx q[2];
rz(-1.1659157) q[2];
sx q[2];
rz(1.7083573) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65518889) q[1];
sx q[1];
rz(-2.9826394) q[1];
sx q[1];
rz(-2.7093191) q[1];
x q[2];
rz(-1.0377117) q[3];
sx q[3];
rz(-1.248705) q[3];
sx q[3];
rz(2.0562293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9050682) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(-2.0489676) q[2];
rz(2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(-2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4051751) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.3177692) q[0];
rz(2.4353943) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(2.172519) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5613865) q[0];
sx q[0];
rz(-1.611334) q[0];
sx q[0];
rz(-1.6798621) q[0];
rz(-pi) q[1];
rz(-0.078209608) q[2];
sx q[2];
rz(-2.4259896) q[2];
sx q[2];
rz(0.90925928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.59086696) q[1];
sx q[1];
rz(-0.64462763) q[1];
sx q[1];
rz(-2.4025737) q[1];
rz(-pi) q[2];
rz(-3*pi/11) q[3];
sx q[3];
rz(-2.487605) q[3];
sx q[3];
rz(-1.4531144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51450729) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(0.43760854) q[2];
rz(3.0833516) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(-1.7818041) q[0];
rz(1.4490022) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(-1.70599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5482225) q[0];
sx q[0];
rz(-2.3832294) q[0];
sx q[0];
rz(-2.3286657) q[0];
x q[1];
rz(-2.3195763) q[2];
sx q[2];
rz(-0.026958131) q[2];
sx q[2];
rz(2.6798623) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.83511931) q[1];
sx q[1];
rz(-1.8715579) q[1];
sx q[1];
rz(-1.1919828) q[1];
x q[2];
rz(-0.42922677) q[3];
sx q[3];
rz(-2.0373404) q[3];
sx q[3];
rz(0.45688094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4009565) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(0.17399542) q[2];
rz(-2.1334355) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.290264) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(-0.1329578) q[0];
rz(-0.48775396) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(-1.3652323) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5196461) q[0];
sx q[0];
rz(-1.6462109) q[0];
sx q[0];
rz(1.7072862) q[0];
x q[1];
rz(1.7403142) q[2];
sx q[2];
rz(-3.0119544) q[2];
sx q[2];
rz(-2.4469751) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9284968) q[1];
sx q[1];
rz(-2.1017385) q[1];
sx q[1];
rz(-1.6316158) q[1];
rz(2.7216464) q[3];
sx q[3];
rz(-1.7405602) q[3];
sx q[3];
rz(1.1743634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0662971) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-2.7189642) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(3.1388824) q[0];
sx q[0];
rz(-0.30160987) q[0];
sx q[0];
rz(-0.31059206) q[0];
rz(-2.8839135) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(-0.92393595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5101178) q[0];
sx q[0];
rz(-2.3259813) q[0];
sx q[0];
rz(0.60879137) q[0];
x q[1];
rz(-3.1191101) q[2];
sx q[2];
rz(-0.48869952) q[2];
sx q[2];
rz(-1.894001) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9388705) q[1];
sx q[1];
rz(-2.466423) q[1];
sx q[1];
rz(0.90941456) q[1];
x q[2];
rz(-2.1903673) q[3];
sx q[3];
rz(-2.1853672) q[3];
sx q[3];
rz(2.8545477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.7473934) q[2];
rz(1.1692858) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7681463) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(1.753153) q[0];
rz(-3.1346079) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(3.1034234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6258433) q[0];
sx q[0];
rz(-2.1793097) q[0];
sx q[0];
rz(3.1090816) q[0];
x q[1];
rz(-1.3456887) q[2];
sx q[2];
rz(-2.7499866) q[2];
sx q[2];
rz(3.0181146) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0109509) q[1];
sx q[1];
rz(-2.7656733) q[1];
sx q[1];
rz(-1.6846659) q[1];
rz(0.60572259) q[3];
sx q[3];
rz(-1.0014135) q[3];
sx q[3];
rz(-3.0063418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.96597153) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(-1.1981298) q[2];
rz(-2.0576058) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(0.23588022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(1.9772298) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(1.6434796) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(-1.5126363) q[2];
sx q[2];
rz(-1.460961) q[2];
sx q[2];
rz(0.030208781) q[2];
rz(-0.86919541) q[3];
sx q[3];
rz(-1.3251318) q[3];
sx q[3];
rz(-1.3983923) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
