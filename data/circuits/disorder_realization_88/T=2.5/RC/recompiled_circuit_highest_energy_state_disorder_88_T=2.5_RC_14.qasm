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
rz(3.1495659) q[0];
sx q[0];
rz(6.719448) q[0];
sx q[0];
rz(8.414581) q[0];
rz(0.17207347) q[1];
sx q[1];
rz(7.0893256) q[1];
sx q[1];
rz(10.448699) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9925594) q[0];
sx q[0];
rz(-1.5639389) q[0];
sx q[0];
rz(-2.0586781) q[0];
rz(-pi) q[1];
rz(0.3918565) q[2];
sx q[2];
rz(-1.7855682) q[2];
sx q[2];
rz(1.2607167) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1271474) q[1];
sx q[1];
rz(-1.5902441) q[1];
sx q[1];
rz(-1.7313185) q[1];
rz(-1.6590622) q[3];
sx q[3];
rz(-1.3219943) q[3];
sx q[3];
rz(0.036347957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1930834) q[2];
sx q[2];
rz(-2.812959) q[2];
sx q[2];
rz(0.23552093) q[2];
rz(-2.113302) q[3];
sx q[3];
rz(-1.8833269) q[3];
sx q[3];
rz(1.9839015) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.704772) q[0];
sx q[0];
rz(-1.7828159) q[0];
sx q[0];
rz(0.92192465) q[0];
rz(-3.0251265) q[1];
sx q[1];
rz(-1.4215697) q[1];
sx q[1];
rz(0.88752735) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5320547) q[0];
sx q[0];
rz(-2.9132151) q[0];
sx q[0];
rz(-2.1837932) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9824184) q[2];
sx q[2];
rz(-0.35458699) q[2];
sx q[2];
rz(1.2693506) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4749467) q[1];
sx q[1];
rz(-1.6837059) q[1];
sx q[1];
rz(-2.5700047) q[1];
rz(-1.4027897) q[3];
sx q[3];
rz(-0.6989494) q[3];
sx q[3];
rz(-2.2022171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.55887115) q[2];
sx q[2];
rz(-2.116394) q[2];
sx q[2];
rz(2.9026046) q[2];
rz(-2.41921) q[3];
sx q[3];
rz(-0.84010774) q[3];
sx q[3];
rz(1.1649789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31579414) q[0];
sx q[0];
rz(-0.9676942) q[0];
sx q[0];
rz(-0.78829366) q[0];
rz(-1.3901419) q[1];
sx q[1];
rz(-2.7056521) q[1];
sx q[1];
rz(-1.4424666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9787131) q[0];
sx q[0];
rz(-1.7121234) q[0];
sx q[0];
rz(-0.070959758) q[0];
x q[1];
rz(-2.8462538) q[2];
sx q[2];
rz(-1.9304233) q[2];
sx q[2];
rz(-2.4725898) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.07787598) q[1];
sx q[1];
rz(-0.87802475) q[1];
sx q[1];
rz(-0.31097842) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9628311) q[3];
sx q[3];
rz(-1.074665) q[3];
sx q[3];
rz(-2.2876491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4060789) q[2];
sx q[2];
rz(-1.9903851) q[2];
sx q[2];
rz(0.25685143) q[2];
rz(0.86366051) q[3];
sx q[3];
rz(-1.0830027) q[3];
sx q[3];
rz(0.5736205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0470806) q[0];
sx q[0];
rz(-3.1086476) q[0];
sx q[0];
rz(1.5091913) q[0];
rz(0.31653658) q[1];
sx q[1];
rz(-1.5455952) q[1];
sx q[1];
rz(1.0155771) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94293649) q[0];
sx q[0];
rz(-1.1186677) q[0];
sx q[0];
rz(2.6200635) q[0];
rz(-2.8695753) q[2];
sx q[2];
rz(-0.81230703) q[2];
sx q[2];
rz(-0.35292948) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8997422) q[1];
sx q[1];
rz(-2.6961245) q[1];
sx q[1];
rz(2.1136606) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7523132) q[3];
sx q[3];
rz(-2.3861775) q[3];
sx q[3];
rz(0.28873539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29869002) q[2];
sx q[2];
rz(-2.0899253) q[2];
sx q[2];
rz(1.451937) q[2];
rz(1.9076094) q[3];
sx q[3];
rz(-1.0943639) q[3];
sx q[3];
rz(2.2598677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7538309) q[0];
sx q[0];
rz(-1.3085288) q[0];
sx q[0];
rz(1.0863139) q[0];
rz(-1.7408675) q[1];
sx q[1];
rz(-2.3298405) q[1];
sx q[1];
rz(3.1382255) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683748) q[0];
sx q[0];
rz(-1.5649319) q[0];
sx q[0];
rz(-0.59451367) q[0];
rz(2.1589627) q[2];
sx q[2];
rz(-1.8525043) q[2];
sx q[2];
rz(-2.1717193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6033404) q[1];
sx q[1];
rz(-1.0639186) q[1];
sx q[1];
rz(1.6350621) q[1];
rz(-3.1037381) q[3];
sx q[3];
rz(-2.541171) q[3];
sx q[3];
rz(-0.63943938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.1668261) q[2];
sx q[2];
rz(-0.70009118) q[2];
sx q[2];
rz(-0.33494803) q[2];
rz(2.8454928) q[3];
sx q[3];
rz(-0.94303232) q[3];
sx q[3];
rz(-0.12224841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5122546) q[0];
sx q[0];
rz(-2.6423995) q[0];
sx q[0];
rz(2.7235624) q[0];
rz(-0.66608518) q[1];
sx q[1];
rz(-1.8981551) q[1];
sx q[1];
rz(0.24806771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9739415) q[0];
sx q[0];
rz(-0.60567666) q[0];
sx q[0];
rz(0.96346897) q[0];
x q[1];
rz(-0.054912189) q[2];
sx q[2];
rz(-1.1217199) q[2];
sx q[2];
rz(-0.95740805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0937178) q[1];
sx q[1];
rz(-1.1823726) q[1];
sx q[1];
rz(-2.1234496) q[1];
rz(0.25603981) q[3];
sx q[3];
rz(-2.1845792) q[3];
sx q[3];
rz(-2.7308794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4416113) q[2];
sx q[2];
rz(-1.2629843) q[2];
sx q[2];
rz(-1.1629533) q[2];
rz(0.59398389) q[3];
sx q[3];
rz(-1.6006399) q[3];
sx q[3];
rz(-0.99015132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.537259) q[0];
sx q[0];
rz(-2.8626677) q[0];
sx q[0];
rz(-0.064924031) q[0];
rz(-2.7760432) q[1];
sx q[1];
rz(-1.431501) q[1];
sx q[1];
rz(-0.72986594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94862969) q[0];
sx q[0];
rz(-2.2092487) q[0];
sx q[0];
rz(0.13704152) q[0];
rz(-pi) q[1];
rz(-1.7634298) q[2];
sx q[2];
rz(-2.3245077) q[2];
sx q[2];
rz(-1.1258923) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0893475) q[1];
sx q[1];
rz(-2.5762853) q[1];
sx q[1];
rz(-1.3924527) q[1];
rz(-0.6859996) q[3];
sx q[3];
rz(-1.8988052) q[3];
sx q[3];
rz(-0.48638108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87021747) q[2];
sx q[2];
rz(-1.7729746) q[2];
sx q[2];
rz(-1.082487) q[2];
rz(1.6535951) q[3];
sx q[3];
rz(-0.36181417) q[3];
sx q[3];
rz(-0.35058072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7545886) q[0];
sx q[0];
rz(-0.8828187) q[0];
sx q[0];
rz(-0.3311232) q[0];
rz(1.4777199) q[1];
sx q[1];
rz(-2.176087) q[1];
sx q[1];
rz(0.32041916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5950111) q[0];
sx q[0];
rz(-3.0394331) q[0];
sx q[0];
rz(0.25545303) q[0];
x q[1];
rz(-2.3829303) q[2];
sx q[2];
rz(-0.21100036) q[2];
sx q[2];
rz(1.768844) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2843208) q[1];
sx q[1];
rz(-1.1445938) q[1];
sx q[1];
rz(-1.5878452) q[1];
rz(1.6410758) q[3];
sx q[3];
rz(-2.456074) q[3];
sx q[3];
rz(0.34834114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42935818) q[2];
sx q[2];
rz(-1.4171436) q[2];
sx q[2];
rz(2.7799535) q[2];
rz(-2.2278191) q[3];
sx q[3];
rz(-2.8445966) q[3];
sx q[3];
rz(-0.44338068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8906032) q[0];
sx q[0];
rz(-0.78762233) q[0];
sx q[0];
rz(0.11601624) q[0];
rz(1.8138255) q[1];
sx q[1];
rz(-1.1632183) q[1];
sx q[1];
rz(-1.9414925) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982298) q[0];
sx q[0];
rz(-1.589601) q[0];
sx q[0];
rz(0.19143611) q[0];
rz(1.7723254) q[2];
sx q[2];
rz(-2.4459165) q[2];
sx q[2];
rz(1.8303378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1397651) q[1];
sx q[1];
rz(-1.5561625) q[1];
sx q[1];
rz(0.56165265) q[1];
x q[2];
rz(-2.9755244) q[3];
sx q[3];
rz(-1.8072053) q[3];
sx q[3];
rz(-2.0061559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3001331) q[2];
sx q[2];
rz(-1.7068784) q[2];
sx q[2];
rz(0.33162281) q[2];
rz(-0.2837818) q[3];
sx q[3];
rz(-1.1016568) q[3];
sx q[3];
rz(-2.7200429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3109741) q[0];
sx q[0];
rz(-1.1116488) q[0];
sx q[0];
rz(-2.9750138) q[0];
rz(-0.84959787) q[1];
sx q[1];
rz(-2.6764937) q[1];
sx q[1];
rz(3.0679682) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62944996) q[0];
sx q[0];
rz(-0.76829043) q[0];
sx q[0];
rz(-2.7714202) q[0];
rz(-pi) q[1];
rz(-2.6643986) q[2];
sx q[2];
rz(-1.0721803) q[2];
sx q[2];
rz(2.3248364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60241781) q[1];
sx q[1];
rz(-0.6197558) q[1];
sx q[1];
rz(2.5386993) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62617733) q[3];
sx q[3];
rz(-1.9708084) q[3];
sx q[3];
rz(-1.7020066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.9061884) q[2];
sx q[2];
rz(-0.81934682) q[2];
sx q[2];
rz(0.6582312) q[2];
rz(1.8002347) q[3];
sx q[3];
rz(-2.1737289) q[3];
sx q[3];
rz(-0.85097504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8258719) q[0];
sx q[0];
rz(-0.82875874) q[0];
sx q[0];
rz(2.2525633) q[0];
rz(-0.52147621) q[1];
sx q[1];
rz(-1.5670525) q[1];
sx q[1];
rz(-1.570931) q[1];
rz(-0.53311382) q[2];
sx q[2];
rz(-0.86400142) q[2];
sx q[2];
rz(0.17178847) q[2];
rz(0.23052774) q[3];
sx q[3];
rz(-2.5067825) q[3];
sx q[3];
rz(1.4092177) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
