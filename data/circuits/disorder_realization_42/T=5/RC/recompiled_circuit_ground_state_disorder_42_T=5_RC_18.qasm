OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.2494994) q[0];
sx q[0];
rz(-1.9568994) q[0];
sx q[0];
rz(1.1353593) q[0];
rz(2.7197977) q[1];
sx q[1];
rz(-0.78039688) q[1];
sx q[1];
rz(-2.5193522) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61257768) q[0];
sx q[0];
rz(-1.4473995) q[0];
sx q[0];
rz(-1.4234659) q[0];
x q[1];
rz(-1.1235302) q[2];
sx q[2];
rz(-1.6058927) q[2];
sx q[2];
rz(0.54151565) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5075818) q[1];
sx q[1];
rz(-2.1651499) q[1];
sx q[1];
rz(-1.5848914) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22316605) q[3];
sx q[3];
rz(-1.0112178) q[3];
sx q[3];
rz(2.1055438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8047831) q[2];
sx q[2];
rz(-1.4562891) q[2];
sx q[2];
rz(-3.0554904) q[2];
rz(2.6297249) q[3];
sx q[3];
rz(-2.0972926) q[3];
sx q[3];
rz(2.4477203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43857273) q[0];
sx q[0];
rz(-2.1874671) q[0];
sx q[0];
rz(-0.91617209) q[0];
rz(0.29197261) q[1];
sx q[1];
rz(-2.5501854) q[1];
sx q[1];
rz(-2.3672262) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60304927) q[0];
sx q[0];
rz(-1.2197681) q[0];
sx q[0];
rz(-2.8013264) q[0];
rz(-pi) q[1];
rz(1.2254481) q[2];
sx q[2];
rz(-1.576445) q[2];
sx q[2];
rz(-0.7330837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9271597) q[1];
sx q[1];
rz(-1.6937052) q[1];
sx q[1];
rz(-1.6729808) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3182358) q[3];
sx q[3];
rz(-1.9867479) q[3];
sx q[3];
rz(0.49440629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.795221) q[2];
sx q[2];
rz(-2.5796311) q[2];
sx q[2];
rz(0.58491659) q[2];
rz(-1.4784721) q[3];
sx q[3];
rz(-1.1774747) q[3];
sx q[3];
rz(-1.4896721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65603489) q[0];
sx q[0];
rz(-2.1182883) q[0];
sx q[0];
rz(-2.5265332) q[0];
rz(0.3282322) q[1];
sx q[1];
rz(-1.0098207) q[1];
sx q[1];
rz(-1.0438017) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7481578) q[0];
sx q[0];
rz(-1.5183477) q[0];
sx q[0];
rz(1.5657052) q[0];
rz(-pi) q[1];
rz(-3.1310668) q[2];
sx q[2];
rz(-2.7132456) q[2];
sx q[2];
rz(2.7017252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7134628) q[1];
sx q[1];
rz(-2.2183811) q[1];
sx q[1];
rz(-0.50010276) q[1];
x q[2];
rz(0.7316771) q[3];
sx q[3];
rz(-2.6291375) q[3];
sx q[3];
rz(-2.3872914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7528167) q[2];
sx q[2];
rz(-0.39377585) q[2];
sx q[2];
rz(-2.8326995) q[2];
rz(-2.654352) q[3];
sx q[3];
rz(-1.4332708) q[3];
sx q[3];
rz(2.2545831) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3798645) q[0];
sx q[0];
rz(-2.3966615) q[0];
sx q[0];
rz(-2.6640025) q[0];
rz(1.8567122) q[1];
sx q[1];
rz(-1.711859) q[1];
sx q[1];
rz(-0.46599785) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4728238) q[0];
sx q[0];
rz(-1.4942989) q[0];
sx q[0];
rz(-0.55474218) q[0];
rz(-2.3498769) q[2];
sx q[2];
rz(-0.60388619) q[2];
sx q[2];
rz(1.1271267) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5090967) q[1];
sx q[1];
rz(-1.5358616) q[1];
sx q[1];
rz(2.3216802) q[1];
rz(-0.70798812) q[3];
sx q[3];
rz(-2.3205415) q[3];
sx q[3];
rz(2.193303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3232702) q[2];
sx q[2];
rz(-0.44760901) q[2];
sx q[2];
rz(1.4480048) q[2];
rz(-1.4130392) q[3];
sx q[3];
rz(-1.6084684) q[3];
sx q[3];
rz(-1.1477227) q[3];
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
rz(2.2959761) q[0];
sx q[0];
rz(-2.4281261) q[0];
sx q[0];
rz(0.071685858) q[0];
rz(1.6638311) q[1];
sx q[1];
rz(-1.7774589) q[1];
sx q[1];
rz(0.62612265) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6090995) q[0];
sx q[0];
rz(-2.3029444) q[0];
sx q[0];
rz(-3.0200543) q[0];
x q[1];
rz(-3.0813498) q[2];
sx q[2];
rz(-0.92298302) q[2];
sx q[2];
rz(-2.197165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71264919) q[1];
sx q[1];
rz(-0.82368851) q[1];
sx q[1];
rz(0.23703833) q[1];
rz(-0.463571) q[3];
sx q[3];
rz(-0.39059535) q[3];
sx q[3];
rz(2.5724346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46586299) q[2];
sx q[2];
rz(-0.53368038) q[2];
sx q[2];
rz(-1.6443058) q[2];
rz(-0.40694445) q[3];
sx q[3];
rz(-0.99138433) q[3];
sx q[3];
rz(-2.2814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96754718) q[0];
sx q[0];
rz(-1.2977192) q[0];
sx q[0];
rz(-0.26200727) q[0];
rz(1.2213446) q[1];
sx q[1];
rz(-1.5120993) q[1];
sx q[1];
rz(1.9409723) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32258275) q[0];
sx q[0];
rz(-0.8958502) q[0];
sx q[0];
rz(0.92827601) q[0];
rz(-pi) q[1];
rz(-0.85882218) q[2];
sx q[2];
rz(-1.1505652) q[2];
sx q[2];
rz(-0.65606299) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.15156721) q[1];
sx q[1];
rz(-0.56781102) q[1];
sx q[1];
rz(2.0358596) q[1];
rz(-pi) q[2];
rz(-0.035484826) q[3];
sx q[3];
rz(-1.662613) q[3];
sx q[3];
rz(-0.47210708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1205341) q[2];
sx q[2];
rz(-2.6594682) q[2];
sx q[2];
rz(-2.6684707) q[2];
rz(1.037723) q[3];
sx q[3];
rz(-2.7389052) q[3];
sx q[3];
rz(-2.0033526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4859908) q[0];
sx q[0];
rz(-0.36360535) q[0];
sx q[0];
rz(2.3714491) q[0];
rz(0.4153525) q[1];
sx q[1];
rz(-1.628592) q[1];
sx q[1];
rz(-1.279668) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9933536) q[0];
sx q[0];
rz(-0.73303427) q[0];
sx q[0];
rz(1.5336393) q[0];
x q[1];
rz(0.12721956) q[2];
sx q[2];
rz(-1.1423649) q[2];
sx q[2];
rz(-2.3495393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7918167) q[1];
sx q[1];
rz(-2.8703051) q[1];
sx q[1];
rz(-2.1955745) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.018536699) q[3];
sx q[3];
rz(-2.662475) q[3];
sx q[3];
rz(1.1441413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8659849) q[2];
sx q[2];
rz(-1.393968) q[2];
sx q[2];
rz(2.7371791) q[2];
rz(1.154493) q[3];
sx q[3];
rz(-1.5279852) q[3];
sx q[3];
rz(1.9563458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.481021) q[0];
sx q[0];
rz(-0.72089973) q[0];
sx q[0];
rz(2.639556) q[0];
rz(2.1227116) q[1];
sx q[1];
rz(-1.6682245) q[1];
sx q[1];
rz(0.92207164) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6630423) q[0];
sx q[0];
rz(-1.6660569) q[0];
sx q[0];
rz(-2.8956198) q[0];
rz(-pi) q[1];
rz(0.62412213) q[2];
sx q[2];
rz(-1.7723473) q[2];
sx q[2];
rz(-2.1132642) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0623297) q[1];
sx q[1];
rz(-3.0714066) q[1];
sx q[1];
rz(0.90026469) q[1];
rz(-pi) q[2];
x q[2];
rz(3.139134) q[3];
sx q[3];
rz(-1.4430678) q[3];
sx q[3];
rz(0.14593455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33793882) q[2];
sx q[2];
rz(-1.8622082) q[2];
sx q[2];
rz(-3.0214018) q[2];
rz(0.0029314824) q[3];
sx q[3];
rz(-0.35792297) q[3];
sx q[3];
rz(0.005793747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1015162) q[0];
sx q[0];
rz(-0.67779556) q[0];
sx q[0];
rz(-1.6325604) q[0];
rz(0.82935968) q[1];
sx q[1];
rz(-2.6396535) q[1];
sx q[1];
rz(-1.1306184) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7735159) q[0];
sx q[0];
rz(-2.6720071) q[0];
sx q[0];
rz(-2.3115524) q[0];
x q[1];
rz(2.4254029) q[2];
sx q[2];
rz(-1.1713059) q[2];
sx q[2];
rz(-1.3376118) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3663331) q[1];
sx q[1];
rz(-2.383552) q[1];
sx q[1];
rz(-2.232928) q[1];
rz(-pi) q[2];
rz(1.9958903) q[3];
sx q[3];
rz(-0.55820528) q[3];
sx q[3];
rz(2.3987215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.067554615) q[2];
sx q[2];
rz(-0.87104565) q[2];
sx q[2];
rz(1.3646431) q[2];
rz(-3.126295) q[3];
sx q[3];
rz(-2.3099895) q[3];
sx q[3];
rz(1.8196222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010654) q[0];
sx q[0];
rz(-2.4897713) q[0];
sx q[0];
rz(-0.1524674) q[0];
rz(-1.7561779) q[1];
sx q[1];
rz(-1.0453753) q[1];
sx q[1];
rz(-0.32858953) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5009896) q[0];
sx q[0];
rz(-0.016399212) q[0];
sx q[0];
rz(-2.5522405) q[0];
rz(-pi) q[1];
rz(1.3297021) q[2];
sx q[2];
rz(-2.3179991) q[2];
sx q[2];
rz(1.6833978) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0180792) q[1];
sx q[1];
rz(-1.5410551) q[1];
sx q[1];
rz(1.8090387) q[1];
rz(-2.685052) q[3];
sx q[3];
rz(-2.4437063) q[3];
sx q[3];
rz(1.0029576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80455989) q[2];
sx q[2];
rz(-1.1351265) q[2];
sx q[2];
rz(1.0767153) q[2];
rz(0.95450258) q[3];
sx q[3];
rz(-1.4408305) q[3];
sx q[3];
rz(0.29630989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3267219) q[0];
sx q[0];
rz(-0.9878511) q[0];
sx q[0];
rz(-1.7350154) q[0];
rz(-1.3181435) q[1];
sx q[1];
rz(-1.6271918) q[1];
sx q[1];
rz(0.79961332) q[1];
rz(-0.043371928) q[2];
sx q[2];
rz(-1.5159173) q[2];
sx q[2];
rz(2.2295502) q[2];
rz(-1.7794505) q[3];
sx q[3];
rz(-0.91051523) q[3];
sx q[3];
rz(-0.12403535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
