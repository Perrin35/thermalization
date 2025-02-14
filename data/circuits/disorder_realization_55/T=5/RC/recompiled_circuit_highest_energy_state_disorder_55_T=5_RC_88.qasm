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
rz(-3.0814085) q[0];
sx q[0];
rz(-1.0589851) q[0];
sx q[0];
rz(1.0101779) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(2.8454236) q[1];
sx q[1];
rz(12.256395) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84278569) q[0];
sx q[0];
rz(-1.9901384) q[0];
sx q[0];
rz(-1.692665) q[0];
rz(-pi) q[1];
rz(-0.10711889) q[2];
sx q[2];
rz(-0.20640443) q[2];
sx q[2];
rz(-2.602488) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3821564) q[1];
sx q[1];
rz(-1.2431989) q[1];
sx q[1];
rz(0.20702481) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83155379) q[3];
sx q[3];
rz(-1.420701) q[3];
sx q[3];
rz(-2.7287366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4787204) q[2];
sx q[2];
rz(-0.24820776) q[2];
sx q[2];
rz(0.09566801) q[2];
rz(-0.11224789) q[3];
sx q[3];
rz(-2.2662558) q[3];
sx q[3];
rz(1.4314502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2410759) q[0];
sx q[0];
rz(-0.88923419) q[0];
sx q[0];
rz(-2.526793) q[0];
rz(-0.88848937) q[1];
sx q[1];
rz(-1.5526086) q[1];
sx q[1];
rz(-2.6420171) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.963012) q[0];
sx q[0];
rz(-2.7546429) q[0];
sx q[0];
rz(2.1301756) q[0];
rz(0.26878727) q[2];
sx q[2];
rz(-2.6013881) q[2];
sx q[2];
rz(2.7160783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37812472) q[1];
sx q[1];
rz(-2.2164882) q[1];
sx q[1];
rz(-0.21685361) q[1];
x q[2];
rz(-1.0811483) q[3];
sx q[3];
rz(-0.64125618) q[3];
sx q[3];
rz(-1.7288417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8932314) q[2];
sx q[2];
rz(-1.2373368) q[2];
sx q[2];
rz(-0.020817967) q[2];
rz(1.9013532) q[3];
sx q[3];
rz(-1.4695243) q[3];
sx q[3];
rz(-1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6388539) q[0];
sx q[0];
rz(-0.26832142) q[0];
sx q[0];
rz(-2.8079206) q[0];
rz(-0.73792136) q[1];
sx q[1];
rz(-1.2260022) q[1];
sx q[1];
rz(-3.0121682) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5851916) q[0];
sx q[0];
rz(-1.7298152) q[0];
sx q[0];
rz(2.2521453) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1652075) q[2];
sx q[2];
rz(-0.98862851) q[2];
sx q[2];
rz(0.79603031) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5065352) q[1];
sx q[1];
rz(-1.660907) q[1];
sx q[1];
rz(-1.4485301) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78227104) q[3];
sx q[3];
rz(-1.2201628) q[3];
sx q[3];
rz(-0.80814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1239329) q[2];
sx q[2];
rz(-1.5256226) q[2];
sx q[2];
rz(-1.7714436) q[2];
rz(-2.395199) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(-3.0588176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133404) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(0.56513894) q[0];
rz(0.42916974) q[1];
sx q[1];
rz(-0.64650911) q[1];
sx q[1];
rz(0.82829222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9681184) q[0];
sx q[0];
rz(-0.83195639) q[0];
sx q[0];
rz(-1.678874) q[0];
x q[1];
rz(-1.0617375) q[2];
sx q[2];
rz(-2.5974496) q[2];
sx q[2];
rz(2.6119815) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.1431378) q[1];
sx q[1];
rz(-1.3052485) q[1];
sx q[1];
rz(-0.50190718) q[1];
rz(-pi) q[2];
rz(-0.24803646) q[3];
sx q[3];
rz(-2.0741077) q[3];
sx q[3];
rz(2.4865884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4565178) q[2];
sx q[2];
rz(-1.061331) q[2];
sx q[2];
rz(1.431541) q[2];
rz(-0.93959129) q[3];
sx q[3];
rz(-0.51274931) q[3];
sx q[3];
rz(3.140894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13570304) q[0];
sx q[0];
rz(-2.8764184) q[0];
sx q[0];
rz(1.3294719) q[0];
rz(1.1520518) q[1];
sx q[1];
rz(-2.0538797) q[1];
sx q[1];
rz(-3.1413445) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6085998) q[0];
sx q[0];
rz(-1.3260576) q[0];
sx q[0];
rz(0.18969638) q[0];
x q[1];
rz(1.0365965) q[2];
sx q[2];
rz(-1.6459668) q[2];
sx q[2];
rz(-2.6435564) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2622288) q[1];
sx q[1];
rz(-1.2254224) q[1];
sx q[1];
rz(2.097258) q[1];
rz(-3.0609691) q[3];
sx q[3];
rz(-2.4084457) q[3];
sx q[3];
rz(-0.44185639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1382711) q[2];
sx q[2];
rz(-1.2500117) q[2];
sx q[2];
rz(-3.1357583) q[2];
rz(-0.88998574) q[3];
sx q[3];
rz(-0.68166387) q[3];
sx q[3];
rz(-2.0148923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936845) q[0];
sx q[0];
rz(-1.6981145) q[0];
sx q[0];
rz(-1.6538612) q[0];
rz(0.030390175) q[1];
sx q[1];
rz(-2.0149714) q[1];
sx q[1];
rz(-2.1991275) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68719351) q[0];
sx q[0];
rz(-1.3257349) q[0];
sx q[0];
rz(1.7646199) q[0];
x q[1];
rz(-2.8428188) q[2];
sx q[2];
rz(-1.3232097) q[2];
sx q[2];
rz(0.58591671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8854495) q[1];
sx q[1];
rz(-2.0383308) q[1];
sx q[1];
rz(0.3216775) q[1];
rz(1.0951772) q[3];
sx q[3];
rz(-1.8034569) q[3];
sx q[3];
rz(1.3568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5421062) q[2];
sx q[2];
rz(-2.3607871) q[2];
sx q[2];
rz(-1.3055118) q[2];
rz(-2.5944338) q[3];
sx q[3];
rz(-1.2055509) q[3];
sx q[3];
rz(2.3356596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18043537) q[0];
sx q[0];
rz(-1.0338217) q[0];
sx q[0];
rz(0.10051522) q[0];
rz(-2.6590977) q[1];
sx q[1];
rz(-1.1588187) q[1];
sx q[1];
rz(0.85711342) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2861569) q[0];
sx q[0];
rz(-1.248994) q[0];
sx q[0];
rz(0.67097539) q[0];
rz(-pi) q[1];
rz(1.521342) q[2];
sx q[2];
rz(-0.92453814) q[2];
sx q[2];
rz(0.44909278) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.65322058) q[1];
sx q[1];
rz(-1.778942) q[1];
sx q[1];
rz(0.4076165) q[1];
rz(-1.1825652) q[3];
sx q[3];
rz(-2.534158) q[3];
sx q[3];
rz(-1.5608112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40019217) q[2];
sx q[2];
rz(-0.97101784) q[2];
sx q[2];
rz(-0.3024438) q[2];
rz(2.1619469) q[3];
sx q[3];
rz(-2.1392348) q[3];
sx q[3];
rz(1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.361146) q[0];
sx q[0];
rz(-3.0028711) q[0];
sx q[0];
rz(0.030990344) q[0];
rz(2.567645) q[1];
sx q[1];
rz(-1.6684883) q[1];
sx q[1];
rz(-1.2773638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5073858) q[0];
sx q[0];
rz(-1.9045826) q[0];
sx q[0];
rz(-0.35315634) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8988016) q[2];
sx q[2];
rz(-1.0848248) q[2];
sx q[2];
rz(-0.15835855) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0884235) q[1];
sx q[1];
rz(-0.49234566) q[1];
sx q[1];
rz(1.7264992) q[1];
rz(1.2856917) q[3];
sx q[3];
rz(-1.630097) q[3];
sx q[3];
rz(1.7948732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.10451) q[2];
sx q[2];
rz(-3.0057378) q[2];
sx q[2];
rz(-2.6541397) q[2];
rz(-0.69665748) q[3];
sx q[3];
rz(-0.928855) q[3];
sx q[3];
rz(-0.063974403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071851991) q[0];
sx q[0];
rz(-0.47436473) q[0];
sx q[0];
rz(2.790614) q[0];
rz(-0.17414302) q[1];
sx q[1];
rz(-1.5085647) q[1];
sx q[1];
rz(1.7399656) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.642386) q[0];
sx q[0];
rz(-2.8642352) q[0];
sx q[0];
rz(-0.46210285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0171595) q[2];
sx q[2];
rz(-2.2702262) q[2];
sx q[2];
rz(-2.1289189) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1593139) q[1];
sx q[1];
rz(-1.754292) q[1];
sx q[1];
rz(0.029551701) q[1];
rz(-0.97934874) q[3];
sx q[3];
rz(-0.92100793) q[3];
sx q[3];
rz(-1.5929008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.031124) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(2.5049211) q[2];
rz(1.1104256) q[3];
sx q[3];
rz(-1.5444376) q[3];
sx q[3];
rz(-0.5184263) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7745895) q[0];
sx q[0];
rz(-0.29357266) q[0];
sx q[0];
rz(2.0830182) q[0];
rz(-3.1029347) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(-1.0640594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23180873) q[0];
sx q[0];
rz(-1.5663106) q[0];
sx q[0];
rz(1.5741328) q[0];
x q[1];
rz(2.3346155) q[2];
sx q[2];
rz(-0.39013559) q[2];
sx q[2];
rz(-1.5240508) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1326931) q[1];
sx q[1];
rz(-1.510396) q[1];
sx q[1];
rz(-1.5302883) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6820716) q[3];
sx q[3];
rz(-1.6020613) q[3];
sx q[3];
rz(2.7783436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.71558636) q[2];
sx q[2];
rz(-0.49804372) q[2];
sx q[2];
rz(3.0675724) q[2];
rz(2.7600539) q[3];
sx q[3];
rz(-1.8266725) q[3];
sx q[3];
rz(-2.7874302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1114125) q[0];
sx q[0];
rz(-1.3127865) q[0];
sx q[0];
rz(2.4319613) q[0];
rz(0.61182712) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(1.9520957) q[2];
sx q[2];
rz(-2.0443889) q[2];
sx q[2];
rz(-2.2257795) q[2];
rz(2.0430123) q[3];
sx q[3];
rz(-2.2618812) q[3];
sx q[3];
rz(-1.19899) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
