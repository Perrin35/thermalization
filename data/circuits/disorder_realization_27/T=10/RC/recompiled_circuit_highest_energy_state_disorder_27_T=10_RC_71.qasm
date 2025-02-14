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
rz(2.8994695) q[0];
sx q[0];
rz(-2.4162633) q[0];
sx q[0];
rz(2.2348833) q[0];
rz(-0.46696219) q[1];
sx q[1];
rz(4.0443647) q[1];
sx q[1];
rz(9.180896) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8236602) q[0];
sx q[0];
rz(-0.77046227) q[0];
sx q[0];
rz(2.5567358) q[0];
x q[1];
rz(2.2026625) q[2];
sx q[2];
rz(-0.9900695) q[2];
sx q[2];
rz(1.1685358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6441124) q[1];
sx q[1];
rz(-1.6429158) q[1];
sx q[1];
rz(-2.9583272) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21190819) q[3];
sx q[3];
rz(-0.82900199) q[3];
sx q[3];
rz(0.96860368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4563518) q[2];
sx q[2];
rz(-2.6754003) q[2];
sx q[2];
rz(1.0249798) q[2];
rz(-3.0023365) q[3];
sx q[3];
rz(-1.9459629) q[3];
sx q[3];
rz(2.9774184) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76992947) q[0];
sx q[0];
rz(-2.0491845) q[0];
sx q[0];
rz(1.9364233) q[0];
rz(-1.6734164) q[1];
sx q[1];
rz(-1.877715) q[1];
sx q[1];
rz(1.42043) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44180621) q[0];
sx q[0];
rz(-1.4107293) q[0];
sx q[0];
rz(-2.3193441) q[0];
x q[1];
rz(1.3697789) q[2];
sx q[2];
rz(-0.2699142) q[2];
sx q[2];
rz(0.52598276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5788889) q[1];
sx q[1];
rz(-1.3486224) q[1];
sx q[1];
rz(1.9166757) q[1];
rz(-2.6254856) q[3];
sx q[3];
rz(-0.81023216) q[3];
sx q[3];
rz(-2.6339171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9517407) q[2];
sx q[2];
rz(-0.18988374) q[2];
sx q[2];
rz(-0.80113775) q[2];
rz(-2.7028911) q[3];
sx q[3];
rz(-1.2480241) q[3];
sx q[3];
rz(-2.0285105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37757117) q[0];
sx q[0];
rz(-0.29733297) q[0];
sx q[0];
rz(-1.1391033) q[0];
rz(-0.26690075) q[1];
sx q[1];
rz(-1.8903939) q[1];
sx q[1];
rz(-2.7762754) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1068126) q[0];
sx q[0];
rz(-2.0831554) q[0];
sx q[0];
rz(1.6842635) q[0];
x q[1];
rz(0.92411218) q[2];
sx q[2];
rz(-1.1050997) q[2];
sx q[2];
rz(2.6135824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68054861) q[1];
sx q[1];
rz(-0.37255424) q[1];
sx q[1];
rz(-1.7303321) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2868096) q[3];
sx q[3];
rz(-2.0826591) q[3];
sx q[3];
rz(1.3123684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2306564) q[2];
sx q[2];
rz(-0.74030423) q[2];
sx q[2];
rz(0.12348565) q[2];
rz(-2.2423045) q[3];
sx q[3];
rz(-1.6920009) q[3];
sx q[3];
rz(1.9062769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3673636) q[0];
sx q[0];
rz(-1.4643359) q[0];
sx q[0];
rz(2.2520219) q[0];
rz(2.3956237) q[1];
sx q[1];
rz(-1.4624701) q[1];
sx q[1];
rz(1.7568582) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2413413) q[0];
sx q[0];
rz(-1.5129733) q[0];
sx q[0];
rz(3.1067294) q[0];
x q[1];
rz(0.64176527) q[2];
sx q[2];
rz(-1.3744397) q[2];
sx q[2];
rz(1.8066483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7515604) q[1];
sx q[1];
rz(-2.1728325) q[1];
sx q[1];
rz(1.3436418) q[1];
rz(-1.6024578) q[3];
sx q[3];
rz(-1.0370812) q[3];
sx q[3];
rz(-0.16560743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7437462) q[2];
sx q[2];
rz(-0.83796871) q[2];
sx q[2];
rz(-2.1261334) q[2];
rz(3.1189392) q[3];
sx q[3];
rz(-1.3709143) q[3];
sx q[3];
rz(-0.13815752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7669693) q[0];
sx q[0];
rz(-1.2767108) q[0];
sx q[0];
rz(-0.20481566) q[0];
rz(-2.9264033) q[1];
sx q[1];
rz(-0.22061017) q[1];
sx q[1];
rz(2.2664216) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29680934) q[0];
sx q[0];
rz(-1.267485) q[0];
sx q[0];
rz(2.9861197) q[0];
rz(-0.1540252) q[2];
sx q[2];
rz(-0.8416033) q[2];
sx q[2];
rz(-0.83524365) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.2339911) q[1];
sx q[1];
rz(-2.7669766) q[1];
sx q[1];
rz(2.5020155) q[1];
rz(0.93980333) q[3];
sx q[3];
rz(-1.4550337) q[3];
sx q[3];
rz(0.94652688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.508226) q[2];
sx q[2];
rz(-2.2534011) q[2];
sx q[2];
rz(1.4166895) q[2];
rz(-2.21777) q[3];
sx q[3];
rz(-0.97771907) q[3];
sx q[3];
rz(-1.1317322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0210339) q[0];
sx q[0];
rz(-0.73796213) q[0];
sx q[0];
rz(2.293204) q[0];
rz(1.5658762) q[1];
sx q[1];
rz(-1.8137685) q[1];
sx q[1];
rz(-2.1956992) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83505982) q[0];
sx q[0];
rz(-0.9047821) q[0];
sx q[0];
rz(1.205797) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9376603) q[2];
sx q[2];
rz(-0.82101594) q[2];
sx q[2];
rz(-1.7199788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0578445) q[1];
sx q[1];
rz(-1.1473155) q[1];
sx q[1];
rz(0.64238703) q[1];
x q[2];
rz(-0.35870798) q[3];
sx q[3];
rz(-0.51878319) q[3];
sx q[3];
rz(-2.0174873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9653603) q[2];
sx q[2];
rz(-1.9051899) q[2];
sx q[2];
rz(-1.4005093) q[2];
rz(-1.6598353) q[3];
sx q[3];
rz(-1.7912495) q[3];
sx q[3];
rz(-2.9422133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7401212) q[0];
sx q[0];
rz(-2.6417702) q[0];
sx q[0];
rz(-1.2116145) q[0];
rz(0.47677332) q[1];
sx q[1];
rz(-1.8649273) q[1];
sx q[1];
rz(-1.635294) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20902987) q[0];
sx q[0];
rz(-1.2878006) q[0];
sx q[0];
rz(-2.1542633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2477307) q[2];
sx q[2];
rz(-1.3328716) q[2];
sx q[2];
rz(2.0428773) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7249178) q[1];
sx q[1];
rz(-2.43201) q[1];
sx q[1];
rz(0.1898983) q[1];
rz(-pi) q[2];
rz(2.6166418) q[3];
sx q[3];
rz(-0.54407507) q[3];
sx q[3];
rz(1.8067738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3393042) q[2];
sx q[2];
rz(-1.6243287) q[2];
sx q[2];
rz(-2.9586672) q[2];
rz(-2.4798992) q[3];
sx q[3];
rz(-1.2122093) q[3];
sx q[3];
rz(-0.70484149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45109192) q[0];
sx q[0];
rz(-1.6766312) q[0];
sx q[0];
rz(-2.9943384) q[0];
rz(-2.9946949) q[1];
sx q[1];
rz(-0.92120996) q[1];
sx q[1];
rz(1.9619092) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3203293) q[0];
sx q[0];
rz(-2.593145) q[0];
sx q[0];
rz(-2.4571277) q[0];
x q[1];
rz(-2.0496164) q[2];
sx q[2];
rz(-1.1767595) q[2];
sx q[2];
rz(1.0543329) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.98683607) q[1];
sx q[1];
rz(-2.1636226) q[1];
sx q[1];
rz(1.703771) q[1];
x q[2];
rz(1.4612064) q[3];
sx q[3];
rz(-2.3218751) q[3];
sx q[3];
rz(2.3070371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2008721) q[2];
sx q[2];
rz(-0.93738896) q[2];
sx q[2];
rz(2.3737523) q[2];
rz(0.52669865) q[3];
sx q[3];
rz(-1.0517164) q[3];
sx q[3];
rz(2.5832978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4621157) q[0];
sx q[0];
rz(-0.18590346) q[0];
sx q[0];
rz(-2.0515077) q[0];
rz(1.3500301) q[1];
sx q[1];
rz(-1.4005902) q[1];
sx q[1];
rz(-2.7166691) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37575484) q[0];
sx q[0];
rz(-1.100044) q[0];
sx q[0];
rz(-0.063675675) q[0];
x q[1];
rz(-0.64984729) q[2];
sx q[2];
rz(-1.6260626) q[2];
sx q[2];
rz(-1.856232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7560727) q[1];
sx q[1];
rz(-1.2282097) q[1];
sx q[1];
rz(-2.6581061) q[1];
x q[2];
rz(2.0941689) q[3];
sx q[3];
rz(-2.2096388) q[3];
sx q[3];
rz(-2.1372319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9099884) q[2];
sx q[2];
rz(-2.6308172) q[2];
sx q[2];
rz(-2.6166022) q[2];
rz(1.7867583) q[3];
sx q[3];
rz(-2.2460008) q[3];
sx q[3];
rz(1.3625328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43596426) q[0];
sx q[0];
rz(-0.66052496) q[0];
sx q[0];
rz(-0.10667644) q[0];
rz(-2.0872929) q[1];
sx q[1];
rz(-1.7091457) q[1];
sx q[1];
rz(1.7342825) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5085691) q[0];
sx q[0];
rz(-2.3481524) q[0];
sx q[0];
rz(2.2587848) q[0];
x q[1];
rz(-0.58511795) q[2];
sx q[2];
rz(-2.3986534) q[2];
sx q[2];
rz(1.9415426) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.095276621) q[1];
sx q[1];
rz(-1.7937978) q[1];
sx q[1];
rz(-0.15671244) q[1];
x q[2];
rz(2.1945215) q[3];
sx q[3];
rz(-0.92309381) q[3];
sx q[3];
rz(3.090911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4878896) q[2];
sx q[2];
rz(-2.6370878) q[2];
sx q[2];
rz(-0.27296909) q[2];
rz(-2.2702787) q[3];
sx q[3];
rz(-1.4600236) q[3];
sx q[3];
rz(-0.13450809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7731666) q[0];
sx q[0];
rz(-1.1933403) q[0];
sx q[0];
rz(-2.2884952) q[0];
rz(0.38289616) q[1];
sx q[1];
rz(-1.2877512) q[1];
sx q[1];
rz(3.126694) q[1];
rz(1.1287002) q[2];
sx q[2];
rz(-2.0876035) q[2];
sx q[2];
rz(0.58135396) q[2];
rz(-1.4703582) q[3];
sx q[3];
rz(-2.3681869) q[3];
sx q[3];
rz(-2.8988732) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
