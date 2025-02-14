OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2990155) q[0];
sx q[0];
rz(3.314078) q[0];
sx q[0];
rz(9.4077851) q[0];
rz(3.6530082) q[1];
sx q[1];
rz(3.6379171) q[1];
sx q[1];
rz(12.080893) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66390178) q[0];
sx q[0];
rz(-1.6833651) q[0];
sx q[0];
rz(1.6825874) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83037776) q[2];
sx q[2];
rz(-0.86566209) q[2];
sx q[2];
rz(-0.0041088897) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79006195) q[1];
sx q[1];
rz(-2.2021501) q[1];
sx q[1];
rz(0.93598334) q[1];
x q[2];
rz(-3.0270013) q[3];
sx q[3];
rz(-1.262731) q[3];
sx q[3];
rz(1.9730572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39712507) q[2];
sx q[2];
rz(-3.0391389) q[2];
sx q[2];
rz(-0.79258072) q[2];
rz(0.98627311) q[3];
sx q[3];
rz(-1.4873742) q[3];
sx q[3];
rz(3.0084394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88696402) q[0];
sx q[0];
rz(-1.6034842) q[0];
sx q[0];
rz(-0.1057374) q[0];
rz(-2.3836783) q[1];
sx q[1];
rz(-0.71280232) q[1];
sx q[1];
rz(0.47592083) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3041046) q[0];
sx q[0];
rz(-1.6679576) q[0];
sx q[0];
rz(-2.9062953) q[0];
rz(0.49293002) q[2];
sx q[2];
rz(-1.4441737) q[2];
sx q[2];
rz(0.38304087) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8393499) q[1];
sx q[1];
rz(-1.005233) q[1];
sx q[1];
rz(-2.9981705) q[1];
x q[2];
rz(1.1844238) q[3];
sx q[3];
rz(-0.68102057) q[3];
sx q[3];
rz(-1.0134361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9953352) q[2];
sx q[2];
rz(-2.6671851) q[2];
sx q[2];
rz(1.8246626) q[2];
rz(2.3138192) q[3];
sx q[3];
rz(-2.0483978) q[3];
sx q[3];
rz(2.304346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71883172) q[0];
sx q[0];
rz(-1.3439002) q[0];
sx q[0];
rz(-1.9047009) q[0];
rz(-1.3924567) q[1];
sx q[1];
rz(-1.4207276) q[1];
sx q[1];
rz(0.69033355) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77237511) q[0];
sx q[0];
rz(-0.75389426) q[0];
sx q[0];
rz(-2.445141) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0341145) q[2];
sx q[2];
rz(-1.5672188) q[2];
sx q[2];
rz(1.6301375) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0597045) q[1];
sx q[1];
rz(-1.4041384) q[1];
sx q[1];
rz(2.9465066) q[1];
x q[2];
rz(-2.4776633) q[3];
sx q[3];
rz(-0.85278836) q[3];
sx q[3];
rz(-1.4527066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6231125) q[2];
sx q[2];
rz(-1.5223576) q[2];
sx q[2];
rz(3.0677262) q[2];
rz(1.9833924) q[3];
sx q[3];
rz(-1.9246212) q[3];
sx q[3];
rz(-1.4069675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2288007) q[0];
sx q[0];
rz(-3.085629) q[0];
sx q[0];
rz(-2.9486616) q[0];
rz(1.2436766) q[1];
sx q[1];
rz(-1.7453777) q[1];
sx q[1];
rz(2.9085433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10639272) q[0];
sx q[0];
rz(-1.7174722) q[0];
sx q[0];
rz(1.4752409) q[0];
rz(1.6150171) q[2];
sx q[2];
rz(-0.70927519) q[2];
sx q[2];
rz(-2.0486633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31598791) q[1];
sx q[1];
rz(-1.5145166) q[1];
sx q[1];
rz(0.31096267) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95121164) q[3];
sx q[3];
rz(-1.9635634) q[3];
sx q[3];
rz(-2.5529566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9243246) q[2];
sx q[2];
rz(-1.7622207) q[2];
sx q[2];
rz(0.064083727) q[2];
rz(0.25121769) q[3];
sx q[3];
rz(-0.46291864) q[3];
sx q[3];
rz(-0.29138756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083704405) q[0];
sx q[0];
rz(-1.8479481) q[0];
sx q[0];
rz(-1.2842913) q[0];
rz(-3.1341556) q[1];
sx q[1];
rz(-1.1859272) q[1];
sx q[1];
rz(-2.1844905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83643267) q[0];
sx q[0];
rz(-2.0070939) q[0];
sx q[0];
rz(0.35168217) q[0];
x q[1];
rz(2.1726274) q[2];
sx q[2];
rz(-2.0973679) q[2];
sx q[2];
rz(-0.32770448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6308208) q[1];
sx q[1];
rz(-1.8987406) q[1];
sx q[1];
rz(1.4664709) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4864205) q[3];
sx q[3];
rz(-1.7196894) q[3];
sx q[3];
rz(2.6096491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1288422) q[2];
sx q[2];
rz(-0.70657554) q[2];
sx q[2];
rz(-2.1112704) q[2];
rz(0.71850592) q[3];
sx q[3];
rz(-2.0593144) q[3];
sx q[3];
rz(-0.70657402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71501032) q[0];
sx q[0];
rz(-0.74786818) q[0];
sx q[0];
rz(-2.7744875) q[0];
rz(-2.7857065) q[1];
sx q[1];
rz(-1.117319) q[1];
sx q[1];
rz(-1.5918559) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7283495) q[0];
sx q[0];
rz(-1.6093045) q[0];
sx q[0];
rz(-2.4000077) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8086967) q[2];
sx q[2];
rz(-0.44778433) q[2];
sx q[2];
rz(-2.9749123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.15799668) q[1];
sx q[1];
rz(-1.3648197) q[1];
sx q[1];
rz(-0.41974824) q[1];
rz(-pi) q[2];
rz(0.14629062) q[3];
sx q[3];
rz(-0.47463575) q[3];
sx q[3];
rz(0.65505799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1648569) q[2];
sx q[2];
rz(-1.848449) q[2];
sx q[2];
rz(-0.0034275835) q[2];
rz(-0.88611832) q[3];
sx q[3];
rz(-1.0930073) q[3];
sx q[3];
rz(1.6974983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.43585983) q[0];
sx q[0];
rz(-1.6678565) q[0];
sx q[0];
rz(2.4263897) q[0];
rz(-0.12229478) q[1];
sx q[1];
rz(-2.1221752) q[1];
sx q[1];
rz(2.9939647) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1081165) q[0];
sx q[0];
rz(-1.8980935) q[0];
sx q[0];
rz(-1.9294338) q[0];
rz(-1.1806025) q[2];
sx q[2];
rz(-1.1653333) q[2];
sx q[2];
rz(-0.46145876) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1654037) q[1];
sx q[1];
rz(-0.87748412) q[1];
sx q[1];
rz(2.90392) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0711898) q[3];
sx q[3];
rz(-1.8569267) q[3];
sx q[3];
rz(-2.4918258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82835308) q[2];
sx q[2];
rz(-0.30892631) q[2];
sx q[2];
rz(-3.0301136) q[2];
rz(-0.76505032) q[3];
sx q[3];
rz(-1.7181516) q[3];
sx q[3];
rz(-3.1077207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9015123) q[0];
sx q[0];
rz(-0.96918786) q[0];
sx q[0];
rz(-1.0028268) q[0];
rz(0.78549939) q[1];
sx q[1];
rz(-1.6167043) q[1];
sx q[1];
rz(1.2202107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1357546) q[0];
sx q[0];
rz(-1.5501115) q[0];
sx q[0];
rz(0.51405859) q[0];
x q[1];
rz(2.1741232) q[2];
sx q[2];
rz(-0.10987205) q[2];
sx q[2];
rz(-2.6990698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0682098) q[1];
sx q[1];
rz(-1.0676976) q[1];
sx q[1];
rz(-0.36243172) q[1];
rz(-pi) q[2];
rz(-2.972159) q[3];
sx q[3];
rz(-2.1994576) q[3];
sx q[3];
rz(2.7055912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.23585524) q[2];
sx q[2];
rz(-1.7451655) q[2];
sx q[2];
rz(3.1077969) q[2];
rz(0.77364051) q[3];
sx q[3];
rz(-1.6126817) q[3];
sx q[3];
rz(2.0788367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77699023) q[0];
sx q[0];
rz(-2.5211625) q[0];
sx q[0];
rz(-2.799209) q[0];
rz(1.7204334) q[1];
sx q[1];
rz(-2.3183289) q[1];
sx q[1];
rz(-1.8470496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081704126) q[0];
sx q[0];
rz(-1.2098639) q[0];
sx q[0];
rz(1.7287219) q[0];
rz(0.85656389) q[2];
sx q[2];
rz(-2.1571549) q[2];
sx q[2];
rz(-0.34242737) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0659938) q[1];
sx q[1];
rz(-2.3477049) q[1];
sx q[1];
rz(2.9692279) q[1];
rz(-pi) q[2];
rz(2.8043069) q[3];
sx q[3];
rz(-1.7333859) q[3];
sx q[3];
rz(0.43460571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2932059) q[2];
sx q[2];
rz(-1.7561971) q[2];
sx q[2];
rz(2.6050341) q[2];
rz(-0.41283354) q[3];
sx q[3];
rz(-2.6091913) q[3];
sx q[3];
rz(-2.0280973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31529108) q[0];
sx q[0];
rz(-1.8719801) q[0];
sx q[0];
rz(1.4971365) q[0];
rz(0.81028384) q[1];
sx q[1];
rz(-1.5850001) q[1];
sx q[1];
rz(2.2946045) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6893456) q[0];
sx q[0];
rz(-0.23437491) q[0];
sx q[0];
rz(2.330392) q[0];
rz(-0.92310793) q[2];
sx q[2];
rz(-2.0294242) q[2];
sx q[2];
rz(3.0748526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3851953) q[1];
sx q[1];
rz(-2.1902731) q[1];
sx q[1];
rz(-1.2379012) q[1];
x q[2];
rz(1.296792) q[3];
sx q[3];
rz(-1.4782259) q[3];
sx q[3];
rz(2.0333729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0392796) q[2];
sx q[2];
rz(-1.3156834) q[2];
sx q[2];
rz(0.63423356) q[2];
rz(0.63747326) q[3];
sx q[3];
rz(-1.1280779) q[3];
sx q[3];
rz(-2.9243961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1524326) q[0];
sx q[0];
rz(-1.1897054) q[0];
sx q[0];
rz(-1.1783896) q[0];
rz(-2.4329026) q[1];
sx q[1];
rz(-1.3460174) q[1];
sx q[1];
rz(-1.8585471) q[1];
rz(-2.3914349) q[2];
sx q[2];
rz(-0.49396074) q[2];
sx q[2];
rz(0.33057292) q[2];
rz(-2.5051334) q[3];
sx q[3];
rz(-1.2598273) q[3];
sx q[3];
rz(-0.72478404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
