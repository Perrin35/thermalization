OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274974) q[0];
sx q[0];
rz(-0.56199718) q[0];
sx q[0];
rz(-2.9105817) q[0];
rz(-2.8958939) q[1];
sx q[1];
rz(-2.6872771) q[1];
sx q[1];
rz(1.8543724) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2184233) q[0];
sx q[0];
rz(-1.2997264) q[0];
sx q[0];
rz(-2.0725033) q[0];
rz(-pi) q[1];
rz(1.8051992) q[2];
sx q[2];
rz(-1.8097098) q[2];
sx q[2];
rz(-2.6221681) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3188149) q[1];
sx q[1];
rz(-0.19365573) q[1];
sx q[1];
rz(2.6050287) q[1];
rz(1.1679959) q[3];
sx q[3];
rz(-2.9514545) q[3];
sx q[3];
rz(-2.8791752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40453688) q[2];
sx q[2];
rz(-0.70713592) q[2];
sx q[2];
rz(-2.7810968) q[2];
rz(2.0170085) q[3];
sx q[3];
rz(-1.0485317) q[3];
sx q[3];
rz(1.8726965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8812113) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(-0.65688175) q[0];
rz(2.8086713) q[1];
sx q[1];
rz(-1.0924783) q[1];
sx q[1];
rz(0.93516707) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012197709) q[0];
sx q[0];
rz(-0.10888616) q[0];
sx q[0];
rz(0.28277855) q[0];
rz(-pi) q[1];
rz(-0.25721154) q[2];
sx q[2];
rz(-1.4818483) q[2];
sx q[2];
rz(-2.6572029) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6015687) q[1];
sx q[1];
rz(-1.717713) q[1];
sx q[1];
rz(1.5308892) q[1];
rz(1.7788497) q[3];
sx q[3];
rz(-1.2223772) q[3];
sx q[3];
rz(-1.5294242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3071345) q[2];
sx q[2];
rz(-1.5156526) q[2];
sx q[2];
rz(-1.2163986) q[2];
rz(2.6136716) q[3];
sx q[3];
rz(-2.1025751) q[3];
sx q[3];
rz(1.2549887) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6127748) q[0];
sx q[0];
rz(-2.3154494) q[0];
sx q[0];
rz(0.58498996) q[0];
rz(-0.12475573) q[1];
sx q[1];
rz(-0.59586066) q[1];
sx q[1];
rz(1.4488719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3901236) q[0];
sx q[0];
rz(-1.5737783) q[0];
sx q[0];
rz(3.1389159) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0694088) q[2];
sx q[2];
rz(-1.5572786) q[2];
sx q[2];
rz(1.1987606) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8068778) q[1];
sx q[1];
rz(-0.92744614) q[1];
sx q[1];
rz(2.592464) q[1];
x q[2];
rz(-2.4873729) q[3];
sx q[3];
rz(-0.4163792) q[3];
sx q[3];
rz(0.49921303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2584194) q[2];
sx q[2];
rz(-1.0018188) q[2];
sx q[2];
rz(1.6955356) q[2];
rz(-0.7575194) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(-2.3022046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3676753) q[0];
sx q[0];
rz(-0.92905074) q[0];
sx q[0];
rz(-2.0539334) q[0];
rz(0.37519535) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(2.1813724) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.023237) q[0];
sx q[0];
rz(-2.873422) q[0];
sx q[0];
rz(1.7212746) q[0];
rz(-1.4713418) q[2];
sx q[2];
rz(-1.6005472) q[2];
sx q[2];
rz(0.1386252) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.474388) q[1];
sx q[1];
rz(-1.1423813) q[1];
sx q[1];
rz(-2.0640255) q[1];
rz(2.8720565) q[3];
sx q[3];
rz(-1.663066) q[3];
sx q[3];
rz(0.59449457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20597657) q[2];
sx q[2];
rz(-1.2674067) q[2];
sx q[2];
rz(1.6440294) q[2];
rz(-2.2802672) q[3];
sx q[3];
rz(-2.438811) q[3];
sx q[3];
rz(2.6845045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.977026) q[0];
sx q[0];
rz(-0.67104665) q[0];
sx q[0];
rz(2.5373051) q[0];
rz(1.1445649) q[1];
sx q[1];
rz(-1.6555758) q[1];
sx q[1];
rz(-2.7191275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6910962) q[0];
sx q[0];
rz(-1.7393149) q[0];
sx q[0];
rz(0.082534747) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79549148) q[2];
sx q[2];
rz(-2.4078712) q[2];
sx q[2];
rz(2.7249641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8026655) q[1];
sx q[1];
rz(-2.4484854) q[1];
sx q[1];
rz(-3.0285804) q[1];
rz(-pi) q[2];
rz(-1.2286387) q[3];
sx q[3];
rz(-1.2099724) q[3];
sx q[3];
rz(-2.7865041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3182688) q[2];
sx q[2];
rz(-2.4906929) q[2];
sx q[2];
rz(-0.65625119) q[2];
rz(-1.6107791) q[3];
sx q[3];
rz(-1.8717513) q[3];
sx q[3];
rz(2.5608565) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73606473) q[0];
sx q[0];
rz(-2.1298213) q[0];
sx q[0];
rz(0.62498012) q[0];
rz(-0.75195733) q[1];
sx q[1];
rz(-1.0686921) q[1];
sx q[1];
rz(1.6065074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84941712) q[0];
sx q[0];
rz(-1.8596974) q[0];
sx q[0];
rz(-1.6221415) q[0];
rz(2.4703896) q[2];
sx q[2];
rz(-2.3103788) q[2];
sx q[2];
rz(2.9321456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8934763) q[1];
sx q[1];
rz(-1.3424338) q[1];
sx q[1];
rz(2.5121157) q[1];
rz(-2.5967136) q[3];
sx q[3];
rz(-1.0419462) q[3];
sx q[3];
rz(-1.8008055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1890761) q[2];
sx q[2];
rz(-1.7837046) q[2];
sx q[2];
rz(2.2124186) q[2];
rz(2.3164228) q[3];
sx q[3];
rz(-1.5114096) q[3];
sx q[3];
rz(-1.1024124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5850942) q[0];
sx q[0];
rz(-2.3346021) q[0];
sx q[0];
rz(1.7671385) q[0];
rz(0.51042405) q[1];
sx q[1];
rz(-0.7904895) q[1];
sx q[1];
rz(-0.62320954) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0334367) q[0];
sx q[0];
rz(-1.2386453) q[0];
sx q[0];
rz(2.9812814) q[0];
rz(2.8002732) q[2];
sx q[2];
rz(-2.6067197) q[2];
sx q[2];
rz(-2.24077) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0556866) q[1];
sx q[1];
rz(-1.6821563) q[1];
sx q[1];
rz(-1.7176877) q[1];
rz(2.9233452) q[3];
sx q[3];
rz(-0.75957662) q[3];
sx q[3];
rz(1.709136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0906543) q[2];
sx q[2];
rz(-0.81523681) q[2];
sx q[2];
rz(1.1673002) q[2];
rz(3.1189611) q[3];
sx q[3];
rz(-2.6204717) q[3];
sx q[3];
rz(-1.1289271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24564329) q[0];
sx q[0];
rz(-2.0063945) q[0];
sx q[0];
rz(1.0173215) q[0];
rz(-1.7474878) q[1];
sx q[1];
rz(-0.21109763) q[1];
sx q[1];
rz(2.5140433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.548259) q[0];
sx q[0];
rz(-2.152266) q[0];
sx q[0];
rz(-0.47971804) q[0];
x q[1];
rz(-1.6553788) q[2];
sx q[2];
rz(-0.89327565) q[2];
sx q[2];
rz(1.5280452) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4269451) q[1];
sx q[1];
rz(-1.8083739) q[1];
sx q[1];
rz(1.3690884) q[1];
x q[2];
rz(0.47431176) q[3];
sx q[3];
rz(-1.8758869) q[3];
sx q[3];
rz(-2.4449789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5158186) q[2];
sx q[2];
rz(-2.5067582) q[2];
sx q[2];
rz(-0.41075692) q[2];
rz(3.0913894) q[3];
sx q[3];
rz(-1.8886731) q[3];
sx q[3];
rz(-3.0661809) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5589767) q[0];
sx q[0];
rz(-2.2447383) q[0];
sx q[0];
rz(-0.26891747) q[0];
rz(0.31002632) q[1];
sx q[1];
rz(-2.0545484) q[1];
sx q[1];
rz(3.0115829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01154218) q[0];
sx q[0];
rz(-2.7461012) q[0];
sx q[0];
rz(-0.96590913) q[0];
rz(-0.35138826) q[2];
sx q[2];
rz(-1.2341502) q[2];
sx q[2];
rz(-1.2996246) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.42327729) q[1];
sx q[1];
rz(-1.0128824) q[1];
sx q[1];
rz(1.3807135) q[1];
x q[2];
rz(-0.9762398) q[3];
sx q[3];
rz(-2.3701982) q[3];
sx q[3];
rz(0.4679799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2076063) q[2];
sx q[2];
rz(-2.3748368) q[2];
sx q[2];
rz(3.0450191) q[2];
rz(1.6377595) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(0.79969978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12829517) q[0];
sx q[0];
rz(-0.18823637) q[0];
sx q[0];
rz(-2.6085594) q[0];
rz(-0.036272613) q[1];
sx q[1];
rz(-0.78250042) q[1];
sx q[1];
rz(1.689555) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2672469) q[0];
sx q[0];
rz(-1.9479935) q[0];
sx q[0];
rz(-0.95893559) q[0];
x q[1];
rz(0.31970892) q[2];
sx q[2];
rz(-2.8065348) q[2];
sx q[2];
rz(-0.31949319) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25305562) q[1];
sx q[1];
rz(-1.5408684) q[1];
sx q[1];
rz(-0.35945838) q[1];
rz(-0.17590268) q[3];
sx q[3];
rz(-2.6697192) q[3];
sx q[3];
rz(-2.7681818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4286246) q[2];
sx q[2];
rz(-2.4098101) q[2];
sx q[2];
rz(-3.1089605) q[2];
rz(-0.84515682) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(1.0802065) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7964771) q[0];
sx q[0];
rz(-0.65237541) q[0];
sx q[0];
rz(2.8271578) q[0];
rz(0.79700094) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(-1.8664411) q[2];
sx q[2];
rz(-1.0545803) q[2];
sx q[2];
rz(2.7068356) q[2];
rz(3.0211021) q[3];
sx q[3];
rz(-0.81462607) q[3];
sx q[3];
rz(-0.56584384) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
