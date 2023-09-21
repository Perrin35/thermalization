OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(-2.6053083) q[0];
sx q[0];
rz(0.94776881) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(2.1138432) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42067895) q[0];
sx q[0];
rz(-1.5479711) q[0];
sx q[0];
rz(-2.7657397) q[0];
rz(-pi) q[1];
rz(1.0814632) q[2];
sx q[2];
rz(-1.3748193) q[2];
sx q[2];
rz(0.34055647) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.23119588) q[1];
sx q[1];
rz(-2.3379571) q[1];
sx q[1];
rz(2.5956144) q[1];
x q[2];
rz(2.0184228) q[3];
sx q[3];
rz(-2.1493704) q[3];
sx q[3];
rz(0.44958006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(1.0502846) q[2];
rz(-1.1132647) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072409078) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(-0.29775277) q[0];
rz(-2.521926) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(-1.108095) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0034804) q[0];
sx q[0];
rz(-1.7797884) q[0];
sx q[0];
rz(2.1955527) q[0];
x q[1];
rz(-0.96252243) q[2];
sx q[2];
rz(-1.9307185) q[2];
sx q[2];
rz(-0.72812176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3023492) q[1];
sx q[1];
rz(-1.1855372) q[1];
sx q[1];
rz(1.9531996) q[1];
x q[2];
rz(-2.8262029) q[3];
sx q[3];
rz(-0.72761256) q[3];
sx q[3];
rz(1.417158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6790598) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(-2.1662946) q[2];
rz(-0.9179999) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(-0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96238962) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(-2.5986824) q[0];
rz(-0.88223282) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(-0.96484819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436413) q[0];
sx q[0];
rz(-2.7086341) q[0];
sx q[0];
rz(-2.4233682) q[0];
rz(-pi) q[1];
rz(1.6134878) q[2];
sx q[2];
rz(-1.2741158) q[2];
sx q[2];
rz(1.3670849) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.112903) q[1];
sx q[1];
rz(-1.6351055) q[1];
sx q[1];
rz(1.6921922) q[1];
x q[2];
rz(1.987625) q[3];
sx q[3];
rz(-1.7524476) q[3];
sx q[3];
rz(0.75205284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0470011) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(2.0084521) q[2];
rz(2.9099693) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8621181) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(-1.3765155) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(-2.8994697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0644558) q[0];
sx q[0];
rz(-2.4479439) q[0];
sx q[0];
rz(-2.8855188) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7119615) q[2];
sx q[2];
rz(-1.5589082) q[2];
sx q[2];
rz(-2.5846543) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0603234) q[1];
sx q[1];
rz(-2.6594866) q[1];
sx q[1];
rz(-0.37270765) q[1];
x q[2];
rz(0.30947134) q[3];
sx q[3];
rz(-2.7052393) q[3];
sx q[3];
rz(2.5040124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1233998) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(3.0467395) q[2];
rz(-1.799396) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78385329) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(-2.7807996) q[0];
rz(1.3882673) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(-2.0070019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9103694) q[0];
sx q[0];
rz(-1.0783505) q[0];
sx q[0];
rz(-0.72427303) q[0];
rz(2.7822532) q[2];
sx q[2];
rz(-1.9760855) q[2];
sx q[2];
rz(2.5422424) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0134296) q[1];
sx q[1];
rz(-1.1268106) q[1];
sx q[1];
rz(1.3982989) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9633425) q[3];
sx q[3];
rz(-0.44294391) q[3];
sx q[3];
rz(-0.49803842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1363042) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(-0.57265442) q[2];
rz(-0.92875656) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3271493) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(1.2493398) q[0];
rz(-1.2231187) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(-1.1522326) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.457068) q[0];
sx q[0];
rz(-0.10604924) q[0];
sx q[0];
rz(-1.4507136) q[0];
x q[1];
rz(1.9828898) q[2];
sx q[2];
rz(-2.0372143) q[2];
sx q[2];
rz(-2.7127624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8730948) q[1];
sx q[1];
rz(-2.9151306) q[1];
sx q[1];
rz(0.48392673) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4310962) q[3];
sx q[3];
rz(-2.3264255) q[3];
sx q[3];
rz(-1.1991771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6727009) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(-2.1506298) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-2.794054) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836442) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(0.062285475) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(-2.7468162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1226574) q[0];
sx q[0];
rz(-1.1766953) q[0];
sx q[0];
rz(1.8271441) q[0];
rz(0.6997509) q[2];
sx q[2];
rz(-2.664898) q[2];
sx q[2];
rz(1.2303908) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3697606) q[1];
sx q[1];
rz(-0.50060111) q[1];
sx q[1];
rz(2.6066149) q[1];
rz(-pi) q[2];
rz(-1.6077605) q[3];
sx q[3];
rz(-1.4174403) q[3];
sx q[3];
rz(-1.6085898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64849598) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(0.27302343) q[2];
rz(-1.3027044) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(0.31204143) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(1.8564818) q[0];
rz(-1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1788951) q[0];
sx q[0];
rz(-0.51998752) q[0];
sx q[0];
rz(-0.87332256) q[0];
rz(-pi) q[1];
rz(-2.4457744) q[2];
sx q[2];
rz(-1.6098445) q[2];
sx q[2];
rz(-2.2829636) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6337874) q[1];
sx q[1];
rz(-1.2101189) q[1];
sx q[1];
rz(0.37640576) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31524661) q[3];
sx q[3];
rz(-1.1479706) q[3];
sx q[3];
rz(1.3413615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1061219) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(2.1179874) q[2];
rz(2.9566531) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(-2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(1.6212844) q[0];
rz(0.3301436) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(0.83713371) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9159106) q[0];
sx q[0];
rz(-1.7726232) q[0];
sx q[0];
rz(-0.68737824) q[0];
x q[1];
rz(0.57640055) q[2];
sx q[2];
rz(-1.8216368) q[2];
sx q[2];
rz(-0.62945156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9311034) q[1];
sx q[1];
rz(-1.6753907) q[1];
sx q[1];
rz(1.7235669) q[1];
x q[2];
rz(-1.0749531) q[3];
sx q[3];
rz(-1.8284203) q[3];
sx q[3];
rz(0.50592929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60124406) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(-0.99572292) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(-1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(1.0572222) q[0];
rz(3.0341042) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(2.1616139) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996795) q[0];
sx q[0];
rz(-1.4451888) q[0];
sx q[0];
rz(3.0903387) q[0];
rz(-1.0716295) q[2];
sx q[2];
rz(-2.0344779) q[2];
sx q[2];
rz(-2.0641363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6411297) q[1];
sx q[1];
rz(-0.49912057) q[1];
sx q[1];
rz(-1.0541381) q[1];
rz(-pi) q[2];
rz(1.8944593) q[3];
sx q[3];
rz(-1.666288) q[3];
sx q[3];
rz(0.0088866339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8367299) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(-1.0277964) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(0.28361472) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86823157) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(0.17765799) q[2];
sx q[2];
rz(-2.1226317) q[2];
sx q[2];
rz(2.5771099) q[2];
rz(0.50073033) q[3];
sx q[3];
rz(-1.2320319) q[3];
sx q[3];
rz(-2.4945955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
