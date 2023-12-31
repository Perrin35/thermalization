OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(-2.196329) q[0];
sx q[0];
rz(0.52559108) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(2.2367509) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5072767) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(1.0213724) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3741671) q[2];
sx q[2];
rz(-1.7927205) q[2];
sx q[2];
rz(2.1264399) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20476725) q[1];
sx q[1];
rz(-1.1751886) q[1];
sx q[1];
rz(2.8863465) q[1];
rz(-pi) q[2];
rz(2.1078029) q[3];
sx q[3];
rz(-2.783943) q[3];
sx q[3];
rz(1.8358177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20377542) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(-3.0453483) q[2];
rz(2.105666) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7006943) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(-1.2922497) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(-2.4786425) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5040341) q[0];
sx q[0];
rz(-0.090888977) q[0];
sx q[0];
rz(-0.81764098) q[0];
rz(0.40132482) q[2];
sx q[2];
rz(-0.95762816) q[2];
sx q[2];
rz(-0.80712986) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.029164974) q[1];
sx q[1];
rz(-1.6532835) q[1];
sx q[1];
rz(2.6998181) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7327659) q[3];
sx q[3];
rz(-2.0521945) q[3];
sx q[3];
rz(0.045452047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.048916653) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-0.32901397) q[2];
rz(-2.4760903) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(1.3177419) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5927785) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(-2.9192525) q[0];
rz(-2.1242583) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(2.6229048) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1919353) q[0];
sx q[0];
rz(-1.9266832) q[0];
sx q[0];
rz(-2.3200672) q[0];
x q[1];
rz(2.8754183) q[2];
sx q[2];
rz(-1.6632348) q[2];
sx q[2];
rz(0.12649525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31836244) q[1];
sx q[1];
rz(-1.3981817) q[1];
sx q[1];
rz(-2.361972) q[1];
rz(0.054611562) q[3];
sx q[3];
rz(-1.5715989) q[3];
sx q[3];
rz(-0.69566788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2805933) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(3.0730491) q[2];
rz(-2.5391501) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(-3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.65072) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(0.17424507) q[0];
rz(-0.53025591) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(2.6285016) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67447353) q[0];
sx q[0];
rz(-1.4917443) q[0];
sx q[0];
rz(-0.62069686) q[0];
x q[1];
rz(-0.94167534) q[2];
sx q[2];
rz(-2.2721014) q[2];
sx q[2];
rz(2.6710682) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9070223) q[1];
sx q[1];
rz(-1.4325732) q[1];
sx q[1];
rz(2.6837818) q[1];
rz(-1.7771878) q[3];
sx q[3];
rz(-1.5406113) q[3];
sx q[3];
rz(2.9914732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6667368) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(2.9117057) q[2];
rz(0.41904467) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(-2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693562) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(0.67681926) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(0.61606032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4426684) q[0];
sx q[0];
rz(-1.6232423) q[0];
sx q[0];
rz(-0.035874493) q[0];
x q[1];
rz(-2.7763543) q[2];
sx q[2];
rz(-1.6825946) q[2];
sx q[2];
rz(2.0888084) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38356009) q[1];
sx q[1];
rz(-2.5639503) q[1];
sx q[1];
rz(0.010096117) q[1];
rz(-pi) q[2];
rz(-0.30619196) q[3];
sx q[3];
rz(-0.91727835) q[3];
sx q[3];
rz(-0.19838504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4074576) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(0.25536728) q[2];
rz(-1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(-0.42246321) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(-2.6690924) q[0];
rz(-0.52608144) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(-0.88476673) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17381829) q[0];
sx q[0];
rz(-0.38823715) q[0];
sx q[0];
rz(1.3740747) q[0];
x q[1];
rz(3.0753291) q[2];
sx q[2];
rz(-0.21050669) q[2];
sx q[2];
rz(-1.2103684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2192167) q[1];
sx q[1];
rz(-1.9799385) q[1];
sx q[1];
rz(2.8273696) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1690508) q[3];
sx q[3];
rz(-2.5707977) q[3];
sx q[3];
rz(-1.4626383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(-0.3113783) q[2];
rz(-1.7729676) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(-0.061070651) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(-0.73289245) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31156763) q[0];
sx q[0];
rz(-2.4550779) q[0];
sx q[0];
rz(2.4870883) q[0];
x q[1];
rz(-2.1204719) q[2];
sx q[2];
rz(-0.78133821) q[2];
sx q[2];
rz(2.6598425) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9961173) q[1];
sx q[1];
rz(-0.62578326) q[1];
sx q[1];
rz(1.5596418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.93768) q[3];
sx q[3];
rz(-0.11493472) q[3];
sx q[3];
rz(-2.4550408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(0.51458365) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(-2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.085389) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(2.432166) q[0];
rz(1.5052694) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(2.8628796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7473135) q[0];
sx q[0];
rz(-2.8907667) q[0];
sx q[0];
rz(-0.5745116) q[0];
rz(-0.5586241) q[2];
sx q[2];
rz(-1.7785636) q[2];
sx q[2];
rz(-2.7937906) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0875138) q[1];
sx q[1];
rz(-1.2607818) q[1];
sx q[1];
rz(-0.22884303) q[1];
rz(-pi) q[2];
rz(1.3474943) q[3];
sx q[3];
rz(-1.4947065) q[3];
sx q[3];
rz(3.1267816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30148208) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(-2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1466325) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(-0.96889281) q[0];
rz(0.47337198) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.999058) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.105285) q[0];
sx q[0];
rz(-1.4799911) q[0];
sx q[0];
rz(1.275195) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1387024) q[2];
sx q[2];
rz(-0.70334607) q[2];
sx q[2];
rz(-0.10290111) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4559608) q[1];
sx q[1];
rz(-2.5870442) q[1];
sx q[1];
rz(-2.901652) q[1];
rz(-pi) q[2];
rz(-0.90060602) q[3];
sx q[3];
rz(-1.4108676) q[3];
sx q[3];
rz(0.98774324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(2.1208105) q[2];
rz(2.8178689) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(2.1616515) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(-0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(1.2385626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0093875) q[0];
sx q[0];
rz(-1.3552249) q[0];
sx q[0];
rz(1.042004) q[0];
rz(-pi) q[1];
x q[1];
rz(3.086834) q[2];
sx q[2];
rz(-0.75373947) q[2];
sx q[2];
rz(0.71376212) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0760127) q[1];
sx q[1];
rz(-2.2231243) q[1];
sx q[1];
rz(-2.9933762) q[1];
rz(-1.402461) q[3];
sx q[3];
rz(-0.59909648) q[3];
sx q[3];
rz(-1.7408016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4548268) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(3.0977541) q[2];
rz(1.2010126) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(-2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4713521) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(-3.1162221) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(-1.17169) q[2];
sx q[2];
rz(-1.2882971) q[2];
sx q[2];
rz(-0.65297877) q[2];
rz(0.79904859) q[3];
sx q[3];
rz(-1.6934762) q[3];
sx q[3];
rz(1.3140524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
