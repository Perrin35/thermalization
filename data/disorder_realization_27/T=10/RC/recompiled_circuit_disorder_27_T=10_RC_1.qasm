OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.1383837) q[0];
sx q[0];
rz(-2.9870343) q[0];
sx q[0];
rz(2.4490693) q[0];
rz(5.0737557) q[1];
sx q[1];
rz(4.3901246) q[1];
sx q[1];
rz(7.6683383) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801075) q[0];
sx q[0];
rz(-2.2930817) q[0];
sx q[0];
rz(2.0018342) q[0];
rz(0.55180438) q[2];
sx q[2];
rz(-1.8066415) q[2];
sx q[2];
rz(-0.15197309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17978046) q[1];
sx q[1];
rz(-1.4084067) q[1];
sx q[1];
rz(0.23602545) q[1];
rz(-2.5296506) q[3];
sx q[3];
rz(-0.7512593) q[3];
sx q[3];
rz(-0.9179759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2549071) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(-0.20516667) q[2];
rz(2.3702879) q[3];
sx q[3];
rz(-2.3588534) q[3];
sx q[3];
rz(-1.1024968) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7339864) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(2.6876887) q[0];
rz(-2.1167963) q[1];
sx q[1];
rz(-0.4075993) q[1];
sx q[1];
rz(1.9143547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42392143) q[0];
sx q[0];
rz(-2.9917891) q[0];
sx q[0];
rz(-2.1013837) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71364673) q[2];
sx q[2];
rz(-2.4948641) q[2];
sx q[2];
rz(-1.5734067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8933567) q[1];
sx q[1];
rz(-1.3125988) q[1];
sx q[1];
rz(0.30171079) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.482588) q[3];
sx q[3];
rz(-0.25203029) q[3];
sx q[3];
rz(2.8641831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.041302117) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(2.5741637) q[2];
rz(2.7764017) q[3];
sx q[3];
rz(-1.4130211) q[3];
sx q[3];
rz(-2.173483) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48297468) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(0.89865249) q[0];
rz(-2.1458416) q[1];
sx q[1];
rz(-1.5834705) q[1];
sx q[1];
rz(2.8083037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22526564) q[0];
sx q[0];
rz(-0.81575459) q[0];
sx q[0];
rz(2.4832721) q[0];
x q[1];
rz(-0.097924175) q[2];
sx q[2];
rz(-1.3946748) q[2];
sx q[2];
rz(-2.2224561) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2644314) q[1];
sx q[1];
rz(-2.1794771) q[1];
sx q[1];
rz(2.9560637) q[1];
rz(-pi) q[2];
rz(-2.7379052) q[3];
sx q[3];
rz(-1.0567769) q[3];
sx q[3];
rz(2.6727563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68625346) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(-1.9906445) q[2];
rz(2.3006556) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(-1.1545198) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6999917) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(-2.4568795) q[0];
rz(2.1060064) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(-1.205014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8367774) q[0];
sx q[0];
rz(-0.69201058) q[0];
sx q[0];
rz(1.6230323) q[0];
x q[1];
rz(1.3022468) q[2];
sx q[2];
rz(-2.4797202) q[2];
sx q[2];
rz(-2.4115987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0282064) q[1];
sx q[1];
rz(-0.75062597) q[1];
sx q[1];
rz(-1.9764465) q[1];
rz(-pi) q[2];
rz(0.86592309) q[3];
sx q[3];
rz(-1.62543) q[3];
sx q[3];
rz(1.5254525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8923607) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(-0.37115804) q[2];
rz(-1.7403729) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(2.0223117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054759653) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(0.13312419) q[0];
rz(-0.99331028) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(-0.55508074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0017437) q[0];
sx q[0];
rz(-1.886133) q[0];
sx q[0];
rz(0.01339162) q[0];
rz(-pi) q[1];
x q[1];
rz(3.073408) q[2];
sx q[2];
rz(-1.1567332) q[2];
sx q[2];
rz(2.806864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2808387) q[1];
sx q[1];
rz(-0.44563952) q[1];
sx q[1];
rz(1.526236) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5186148) q[3];
sx q[3];
rz(-1.4900041) q[3];
sx q[3];
rz(0.83524708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.83539) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(3.0026657) q[2];
rz(2.1991918) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(-2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54365629) q[0];
sx q[0];
rz(-2.5718226) q[0];
sx q[0];
rz(-2.561835) q[0];
rz(3.014091) q[1];
sx q[1];
rz(-1.189905) q[1];
sx q[1];
rz(1.6019843) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.126639) q[0];
sx q[0];
rz(-0.74288988) q[0];
sx q[0];
rz(-1.2680306) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6164262) q[2];
sx q[2];
rz(-1.4403733) q[2];
sx q[2];
rz(-1.6861196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.786799) q[1];
sx q[1];
rz(-1.2878294) q[1];
sx q[1];
rz(2.016504) q[1];
x q[2];
rz(1.2195915) q[3];
sx q[3];
rz(-0.71577365) q[3];
sx q[3];
rz(0.90482611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55398983) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(-0.26947752) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(3.0814734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6948029) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(0.68429464) q[0];
rz(3.0220095) q[1];
sx q[1];
rz(-1.8493098) q[1];
sx q[1];
rz(0.51876846) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2001901) q[0];
sx q[0];
rz(-1.7046283) q[0];
sx q[0];
rz(2.3076513) q[0];
x q[1];
rz(-0.22612818) q[2];
sx q[2];
rz(-2.3580708) q[2];
sx q[2];
rz(-0.70623159) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8956633) q[1];
sx q[1];
rz(-2.3024125) q[1];
sx q[1];
rz(-2.500446) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.026168907) q[3];
sx q[3];
rz(-1.1702288) q[3];
sx q[3];
rz(2.8699584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8873022) q[2];
sx q[2];
rz(-1.3672978) q[2];
sx q[2];
rz(-0.56345144) q[2];
rz(3.0900132) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(-2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1812487) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(1.3990336) q[0];
rz(2.3545806) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(-2.3972437) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9104011) q[0];
sx q[0];
rz(-1.4230799) q[0];
sx q[0];
rz(-1.5157248) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9888399) q[2];
sx q[2];
rz(-2.0803917) q[2];
sx q[2];
rz(0.60915142) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62521711) q[1];
sx q[1];
rz(-2.5853734) q[1];
sx q[1];
rz(-2.7079627) q[1];
rz(-pi) q[2];
rz(1.5472502) q[3];
sx q[3];
rz(-1.1163201) q[3];
sx q[3];
rz(2.8794895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7156334) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(-2.5320833) q[2];
rz(-0.65731796) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9534) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(-1.6375861) q[0];
rz(-1.9001182) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(0.77493587) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3459754) q[0];
sx q[0];
rz(-0.91751639) q[0];
sx q[0];
rz(0.4972636) q[0];
rz(-pi) q[1];
rz(1.8911792) q[2];
sx q[2];
rz(-2.6777326) q[2];
sx q[2];
rz(1.3055717) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7271125) q[1];
sx q[1];
rz(-2.6467122) q[1];
sx q[1];
rz(0.40463573) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14258607) q[3];
sx q[3];
rz(-1.3538059) q[3];
sx q[3];
rz(0.91689527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.187414) q[2];
sx q[2];
rz(-0.22970197) q[2];
sx q[2];
rz(2.9837218) q[2];
rz(-1.9291417) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(-1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(2.9272595) q[0];
rz(2.4841323) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(-2.0956031) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56453088) q[0];
sx q[0];
rz(-0.37527592) q[0];
sx q[0];
rz(2.8481086) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59098737) q[2];
sx q[2];
rz(-1.5467484) q[2];
sx q[2];
rz(1.1514593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6194832) q[1];
sx q[1];
rz(-1.2088747) q[1];
sx q[1];
rz(-0.98720179) q[1];
rz(-0.40542116) q[3];
sx q[3];
rz(-1.5952205) q[3];
sx q[3];
rz(-0.072211857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41632286) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(-2.0521169) q[2];
rz(1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.2789223) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(1.5325585) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(-2.5074742) q[2];
sx q[2];
rz(-0.49370439) q[2];
sx q[2];
rz(1.6903071) q[2];
rz(-1.6339176) q[3];
sx q[3];
rz(-1.0216542) q[3];
sx q[3];
rz(1.7246104) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];