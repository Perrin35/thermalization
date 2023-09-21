OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(-3.0769899) q[0];
sx q[0];
rz(3.1199772) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(1.8145476) q[1];
sx q[1];
rz(10.756455) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0927147) q[0];
sx q[0];
rz(-1.9034916) q[0];
sx q[0];
rz(-0.96941745) q[0];
rz(-pi) q[1];
rz(0.10279074) q[2];
sx q[2];
rz(-1.79709) q[2];
sx q[2];
rz(0.13470995) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.761258) q[1];
sx q[1];
rz(-0.77612703) q[1];
sx q[1];
rz(0.18591979) q[1];
x q[2];
rz(0.33820037) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(-0.81830922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33064476) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(1.0243246) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.8169096) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(-2.0781793) q[0];
rz(0.88513199) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(-3.1399472) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5981818) q[0];
sx q[0];
rz(-0.046292154) q[0];
sx q[0];
rz(-1.339519) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8284945) q[2];
sx q[2];
rz(-0.94444599) q[2];
sx q[2];
rz(-1.7110273) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8811223) q[1];
sx q[1];
rz(-0.80658856) q[1];
sx q[1];
rz(1.9116198) q[1];
rz(-pi) q[2];
rz(1.315829) q[3];
sx q[3];
rz(-1.4062738) q[3];
sx q[3];
rz(2.4428575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1559747) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(2.4334811) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.920632) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(0.8272585) q[0];
rz(-3.1365085) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(-1.089383) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82747805) q[0];
sx q[0];
rz(-0.70686045) q[0];
sx q[0];
rz(0.88721888) q[0];
rz(-0.59805869) q[2];
sx q[2];
rz(-1.9285678) q[2];
sx q[2];
rz(-1.9269112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.032420302) q[1];
sx q[1];
rz(-1.4854327) q[1];
sx q[1];
rz(1.1847772) q[1];
rz(-pi) q[2];
rz(-2.5855519) q[3];
sx q[3];
rz(-0.87198139) q[3];
sx q[3];
rz(1.3211105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0777145) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(-2.2568978) q[2];
rz(1.1832773) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5723715) q[0];
sx q[0];
rz(-0.5643934) q[0];
sx q[0];
rz(-0.72682056) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(0.35983905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7774178) q[0];
sx q[0];
rz(-0.67876498) q[0];
sx q[0];
rz(-1.2749519) q[0];
rz(-2.4904576) q[2];
sx q[2];
rz(-1.6614117) q[2];
sx q[2];
rz(3.0504984) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6105146) q[1];
sx q[1];
rz(-0.58137608) q[1];
sx q[1];
rz(0.73552144) q[1];
rz(-1.2218277) q[3];
sx q[3];
rz(-1.8604606) q[3];
sx q[3];
rz(-1.7581913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.56746733) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(2.3763669) q[2];
rz(-0.75677538) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(0.24100196) q[3];
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
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1446447) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(-1.416052) q[0];
rz(-0.36711806) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(-1.0353154) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0400378) q[0];
sx q[0];
rz(-0.77417513) q[0];
sx q[0];
rz(2.220962) q[0];
x q[1];
rz(1.5315227) q[2];
sx q[2];
rz(-2.0619259) q[2];
sx q[2];
rz(0.66205762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1086515) q[1];
sx q[1];
rz(-1.6254394) q[1];
sx q[1];
rz(0.095468949) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8661795) q[3];
sx q[3];
rz(-1.7731035) q[3];
sx q[3];
rz(-1.8154669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7845903) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(-2.8055577) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.4679573) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(-1.1451716) q[0];
rz(1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-0.2072269) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6672872) q[0];
sx q[0];
rz(-1.4315839) q[0];
sx q[0];
rz(-0.62367237) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5099498) q[2];
sx q[2];
rz(-0.70438671) q[2];
sx q[2];
rz(-2.7001675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6991899) q[1];
sx q[1];
rz(-1.8132205) q[1];
sx q[1];
rz(1.8276023) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3979785) q[3];
sx q[3];
rz(-2.380905) q[3];
sx q[3];
rz(-0.93483227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(0.72171372) q[2];
rz(-1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24467829) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(0.73202837) q[0];
rz(-0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(0.2917372) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7043982) q[0];
sx q[0];
rz(-1.5196374) q[0];
sx q[0];
rz(-0.41914661) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8740011) q[2];
sx q[2];
rz(-1.6723987) q[2];
sx q[2];
rz(0.50819699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8818672) q[1];
sx q[1];
rz(-1.423466) q[1];
sx q[1];
rz(2.8357361) q[1];
rz(0.50645701) q[3];
sx q[3];
rz(-2.5589057) q[3];
sx q[3];
rz(1.2014233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(2.3042802) q[2];
rz(1.1710179) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(-0.062019197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0983122) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(-3.035416) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(0.95867872) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230969) q[0];
sx q[0];
rz(-2.3765058) q[0];
sx q[0];
rz(2.3575248) q[0];
rz(-pi) q[1];
rz(1.0326951) q[2];
sx q[2];
rz(-1.3287462) q[2];
sx q[2];
rz(-2.6546216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.465185) q[1];
sx q[1];
rz(-2.0945815) q[1];
sx q[1];
rz(0.1233867) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.557425) q[3];
sx q[3];
rz(-1.6730047) q[3];
sx q[3];
rz(-2.971873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0083996) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(0.91910249) q[2];
rz(1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(-1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.33525) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(0.43689716) q[0];
rz(-2.4412952) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(-0.46554309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8252194) q[0];
sx q[0];
rz(-1.6557949) q[0];
sx q[0];
rz(-0.70789106) q[0];
rz(-pi) q[1];
rz(3.0131857) q[2];
sx q[2];
rz(-1.1616716) q[2];
sx q[2];
rz(2.1985334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1381582) q[1];
sx q[1];
rz(-2.5596566) q[1];
sx q[1];
rz(2.483063) q[1];
rz(-1.9302619) q[3];
sx q[3];
rz(-1.3350447) q[3];
sx q[3];
rz(2.204493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7312701) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.643606) q[2];
rz(-2.9368029) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64514226) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(-1.1219332) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(-0.5823935) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56349194) q[0];
sx q[0];
rz(-2.2390215) q[0];
sx q[0];
rz(-2.8816954) q[0];
x q[1];
rz(2.0612129) q[2];
sx q[2];
rz(-1.4368125) q[2];
sx q[2];
rz(-1.4659363) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0092587) q[1];
sx q[1];
rz(-1.7528105) q[1];
sx q[1];
rz(0.11972129) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12606975) q[3];
sx q[3];
rz(-2.9226916) q[3];
sx q[3];
rz(-0.88760469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1311243) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(-2.3790512) q[2];
rz(-1.4108346) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0653771) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(-1.8021884) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(-0.28152485) q[2];
sx q[2];
rz(-2.3073961) q[2];
sx q[2];
rz(-2.6632593) q[2];
rz(1.0675666) q[3];
sx q[3];
rz(-1.1451086) q[3];
sx q[3];
rz(1.2721636) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];