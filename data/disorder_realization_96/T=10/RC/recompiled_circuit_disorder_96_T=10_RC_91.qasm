OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(-2.3242216) q[0];
rz(0.983239) q[1];
sx q[1];
rz(-0.53951889) q[1];
sx q[1];
rz(1.9411545) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49657208) q[0];
sx q[0];
rz(-1.8121769) q[0];
sx q[0];
rz(2.7215331) q[0];
rz(-pi) q[1];
rz(-1.0317694) q[2];
sx q[2];
rz(-1.6845778) q[2];
sx q[2];
rz(0.10345085) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5719205) q[1];
sx q[1];
rz(-1.8351646) q[1];
sx q[1];
rz(-2.4155248) q[1];
rz(-pi) q[2];
rz(-2.82253) q[3];
sx q[3];
rz(-0.82108077) q[3];
sx q[3];
rz(-1.0533489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0212705) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(-1.5585287) q[2];
rz(-2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(1.9460829) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(-2.6057459) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065952259) q[0];
sx q[0];
rz(-1.6951121) q[0];
sx q[0];
rz(0.9888222) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70948647) q[2];
sx q[2];
rz(-1.3777133) q[2];
sx q[2];
rz(-0.62765861) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10837308) q[1];
sx q[1];
rz(-0.48711005) q[1];
sx q[1];
rz(0.56652041) q[1];
rz(-1.4278533) q[3];
sx q[3];
rz(-1.5385475) q[3];
sx q[3];
rz(1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(0.066453233) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41985837) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(-2.9911175) q[0];
rz(0.45723215) q[1];
sx q[1];
rz(-2.229264) q[1];
sx q[1];
rz(3.1157852) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31608554) q[0];
sx q[0];
rz(-1.8206017) q[0];
sx q[0];
rz(-3.120963) q[0];
rz(-pi) q[1];
rz(0.089695887) q[2];
sx q[2];
rz(-1.8811474) q[2];
sx q[2];
rz(0.87583625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4986213) q[1];
sx q[1];
rz(-1.9830623) q[1];
sx q[1];
rz(0.69795124) q[1];
x q[2];
rz(3.0098626) q[3];
sx q[3];
rz(-1.9506491) q[3];
sx q[3];
rz(2.6704138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0187443) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(2.9600926) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(-1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.240775) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(2.2606842) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-2.6054629) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8857408) q[0];
sx q[0];
rz(-0.8937853) q[0];
sx q[0];
rz(-2.6141502) q[0];
rz(-pi) q[1];
rz(0.23668004) q[2];
sx q[2];
rz(-1.6263905) q[2];
sx q[2];
rz(2.6829164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.22524658) q[1];
sx q[1];
rz(-2.1919474) q[1];
sx q[1];
rz(-2.0398554) q[1];
x q[2];
rz(-0.36066182) q[3];
sx q[3];
rz(-0.70913991) q[3];
sx q[3];
rz(-0.053645596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46999103) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(0.70703834) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9064643) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(2.0902324) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(3.0984745) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9163141) q[0];
sx q[0];
rz(-1.7597223) q[0];
sx q[0];
rz(-2.180045) q[0];
rz(-pi) q[1];
rz(2.0448858) q[2];
sx q[2];
rz(-0.57069639) q[2];
sx q[2];
rz(-1.2793465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4823618) q[1];
sx q[1];
rz(-1.2995969) q[1];
sx q[1];
rz(1.2960474) q[1];
x q[2];
rz(1.7205914) q[3];
sx q[3];
rz(-0.78687039) q[3];
sx q[3];
rz(-2.8926135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0126426) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(0.35219231) q[2];
rz(-0.59018618) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(-2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5181638) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(-1.1556926) q[0];
rz(0.75025264) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(2.0828784) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67260325) q[0];
sx q[0];
rz(-2.6575343) q[0];
sx q[0];
rz(-2.6812535) q[0];
rz(-pi) q[1];
rz(1.2054339) q[2];
sx q[2];
rz(-2.0632671) q[2];
sx q[2];
rz(1.1577215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.42029542) q[1];
sx q[1];
rz(-1.8719721) q[1];
sx q[1];
rz(-2.9073614) q[1];
rz(-pi) q[2];
rz(0.70763564) q[3];
sx q[3];
rz(-1.71873) q[3];
sx q[3];
rz(-1.9572452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47026149) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(1.9865215) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.4412122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041615151) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(1.3776243) q[1];
sx q[1];
rz(-0.310985) q[1];
sx q[1];
rz(0.84164936) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0253042) q[0];
sx q[0];
rz(-1.3447273) q[0];
sx q[0];
rz(2.6095005) q[0];
x q[1];
rz(-2.5146671) q[2];
sx q[2];
rz(-1.4142087) q[2];
sx q[2];
rz(-2.2121034) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0824273) q[1];
sx q[1];
rz(-1.6470243) q[1];
sx q[1];
rz(-0.04870292) q[1];
x q[2];
rz(-0.21861403) q[3];
sx q[3];
rz(-0.84158763) q[3];
sx q[3];
rz(-1.5900172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.001361751) q[2];
sx q[2];
rz(-1.7732239) q[2];
sx q[2];
rz(-3.0916396) q[2];
rz(-2.4800381) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426303) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(1.7393973) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(2.7239674) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8861038) q[0];
sx q[0];
rz(-2.4574453) q[0];
sx q[0];
rz(2.4151001) q[0];
rz(-pi) q[1];
rz(1.8629486) q[2];
sx q[2];
rz(-2.7610965) q[2];
sx q[2];
rz(-0.7064864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0954674) q[1];
sx q[1];
rz(-1.0091262) q[1];
sx q[1];
rz(-2.5820877) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1975708) q[3];
sx q[3];
rz(-1.2959359) q[3];
sx q[3];
rz(1.1653792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79545704) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(1.0428838) q[2];
rz(0.67388326) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3089356) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(-0.62966627) q[0];
rz(2.5667403) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(0.94690698) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1614721) q[0];
sx q[0];
rz(-1.8959909) q[0];
sx q[0];
rz(1.8878216) q[0];
rz(-pi) q[1];
rz(-0.097054585) q[2];
sx q[2];
rz(-1.5311054) q[2];
sx q[2];
rz(-1.5511712) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7864101) q[1];
sx q[1];
rz(-1.2372412) q[1];
sx q[1];
rz(-1.4817609) q[1];
rz(-pi) q[2];
rz(-1.1287273) q[3];
sx q[3];
rz(-1.3978492) q[3];
sx q[3];
rz(2.5069619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5809014) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.8927195) q[2];
rz(-2.4272264) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-0.26556382) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666075) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(-0.65504909) q[0];
rz(0.89637268) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(2.4972829) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33675942) q[0];
sx q[0];
rz(-1.5157489) q[0];
sx q[0];
rz(0.17382646) q[0];
x q[1];
rz(2.9115453) q[2];
sx q[2];
rz(-2.2736079) q[2];
sx q[2];
rz(-0.41781296) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5035142) q[1];
sx q[1];
rz(-2.6942309) q[1];
sx q[1];
rz(-0.7034941) q[1];
rz(-pi) q[2];
x q[2];
rz(2.042291) q[3];
sx q[3];
rz(-1.1700556) q[3];
sx q[3];
rz(2.5221962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3506938) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(0.59990668) q[2];
rz(-2.24263) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8469289) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(2.9121493) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(-3.0145666) q[2];
sx q[2];
rz(-0.95022485) q[2];
sx q[2];
rz(-2.7765204) q[2];
rz(-2.7502144) q[3];
sx q[3];
rz(-0.93422514) q[3];
sx q[3];
rz(1.422062) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
