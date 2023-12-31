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
rz(-0.90484172) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68569505) q[0];
sx q[0];
rz(-2.064961) q[0];
sx q[0];
rz(0.49522884) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3741671) q[2];
sx q[2];
rz(-1.3488722) q[2];
sx q[2];
rz(2.1264399) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2658087) q[1];
sx q[1];
rz(-1.335653) q[1];
sx q[1];
rz(1.9782515) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1078029) q[3];
sx q[3];
rz(-2.783943) q[3];
sx q[3];
rz(1.8358177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.20377542) q[2];
sx q[2];
rz(-1.4062466) q[2];
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
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7006943) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(-1.2922497) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(2.4786425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6843296) q[0];
sx q[0];
rz(-1.6329137) q[0];
sx q[0];
rz(1.6371884) q[0];
rz(2.2234369) q[2];
sx q[2];
rz(-1.2456206) q[2];
sx q[2];
rz(2.6174389) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7141124) q[1];
sx q[1];
rz(-0.44891) q[1];
sx q[1];
rz(0.19101363) q[1];
x q[2];
rz(0.29942056) q[3];
sx q[3];
rz(-2.6357108) q[3];
sx q[3];
rz(2.8477856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.092676) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(2.8125787) q[2];
rz(-2.4760903) q[3];
sx q[3];
rz(-2.9232959) q[3];
sx q[3];
rz(1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5488141) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(-1.0173343) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(2.6229048) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9496574) q[0];
sx q[0];
rz(-1.2149095) q[0];
sx q[0];
rz(-2.3200672) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2661744) q[2];
sx q[2];
rz(-1.6632348) q[2];
sx q[2];
rz(3.0150974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0842523) q[1];
sx q[1];
rz(-2.3358314) q[1];
sx q[1];
rz(-1.8112103) q[1];
x q[2];
rz(3.1268901) q[3];
sx q[3];
rz(-3.0869752) q[3];
sx q[3];
rz(-0.86044776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8609994) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(3.0730491) q[2];
rz(-0.60244256) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(-3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.65072) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(2.9673476) q[0];
rz(0.53025591) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(-0.51309103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.355277) q[0];
sx q[0];
rz(-2.516541) q[0];
sx q[0];
rz(-0.13537188) q[0];
rz(0.60855234) q[2];
sx q[2];
rz(-0.90494472) q[2];
sx q[2];
rz(-2.7666639) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2345703) q[1];
sx q[1];
rz(-1.4325732) q[1];
sx q[1];
rz(-2.6837818) q[1];
rz(-pi) q[2];
rz(-0.030839132) q[3];
sx q[3];
rz(-1.7770924) q[3];
sx q[3];
rz(-1.7145969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47485581) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(-0.22988698) q[2];
rz(-2.722548) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(-0.45927799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6722365) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(2.4647734) q[0];
rz(2.6485486) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(-0.61606032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4426684) q[0];
sx q[0];
rz(-1.5183503) q[0];
sx q[0];
rz(0.035874493) q[0];
rz(-0.36523833) q[2];
sx q[2];
rz(-1.6825946) q[2];
sx q[2];
rz(-2.0888084) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7580326) q[1];
sx q[1];
rz(-2.5639503) q[1];
sx q[1];
rz(3.1314965) q[1];
x q[2];
rz(2.8354007) q[3];
sx q[3];
rz(-2.2243143) q[3];
sx q[3];
rz(-2.9432076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(0.25536728) q[2];
rz(-1.6051965) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(-2.6690924) q[0];
rz(-0.52608144) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(2.2568259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17381829) q[0];
sx q[0];
rz(-2.7533555) q[0];
sx q[0];
rz(-1.3740747) q[0];
x q[1];
rz(-0.21005819) q[2];
sx q[2];
rz(-1.5569599) q[2];
sx q[2];
rz(-2.7163598) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2192167) q[1];
sx q[1];
rz(-1.9799385) q[1];
sx q[1];
rz(-0.31422305) q[1];
x q[2];
rz(-2.8956036) q[3];
sx q[3];
rz(-1.0503328) q[3];
sx q[3];
rz(-0.99508475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(-1.7729676) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41806528) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(3.1014077) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(2.4087002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7948579) q[0];
sx q[0];
rz(-1.96694) q[0];
sx q[0];
rz(0.57647716) q[0];
rz(-pi) q[1];
rz(1.0211208) q[2];
sx q[2];
rz(-2.3602544) q[2];
sx q[2];
rz(0.48175016) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9961173) q[1];
sx q[1];
rz(-0.62578326) q[1];
sx q[1];
rz(-1.5819509) q[1];
rz(-pi) q[2];
rz(0.041386889) q[3];
sx q[3];
rz(-1.4635411) q[3];
sx q[3];
rz(1.0556575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0682893) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(0.51458365) q[2];
rz(-1.9040646) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(-2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.085389) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(-2.432166) q[0];
rz(1.5052694) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(-0.27871305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3942791) q[0];
sx q[0];
rz(-0.25082591) q[0];
sx q[0];
rz(2.5670811) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7630373) q[2];
sx q[2];
rz(-2.5494529) q[2];
sx q[2];
rz(-0.90422599) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41234327) q[1];
sx q[1];
rz(-1.7885498) q[1];
sx q[1];
rz(1.8885683) q[1];
rz(-pi) q[2];
x q[2];
rz(0.078019402) q[3];
sx q[3];
rz(-1.7934414) q[3];
sx q[3];
rz(-1.5732461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8401106) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(-2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99496019) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(0.96889281) q[0];
rz(2.6682207) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(-1.1425346) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0363077) q[0];
sx q[0];
rz(-1.6616016) q[0];
sx q[0];
rz(1.275195) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5683453) q[2];
sx q[2];
rz(-0.8674538) q[2];
sx q[2];
rz(3.0424812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7361703) q[1];
sx q[1];
rz(-2.1076964) q[1];
sx q[1];
rz(1.4246529) q[1];
rz(-pi) q[2];
rz(1.8248796) q[3];
sx q[3];
rz(-0.68613201) q[3];
sx q[3];
rz(-0.38476598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.896495) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(-1.0207821) q[2];
rz(2.8178689) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(3.1304205) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616515) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(-2.4401869) q[0];
rz(2.2258863) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(1.2385626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21047132) q[0];
sx q[0];
rz(-0.56715542) q[0];
sx q[0];
rz(1.9803067) q[0];
rz(-1.5194703) q[2];
sx q[2];
rz(-0.81846279) q[2];
sx q[2];
rz(-2.3527956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.065579942) q[1];
sx q[1];
rz(-2.2231243) q[1];
sx q[1];
rz(-0.14821649) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11390399) q[3];
sx q[3];
rz(-2.1602727) q[3];
sx q[3];
rz(-1.9437499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68676585) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-0.043838538) q[2];
rz(1.2010126) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4713521) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(3.1162221) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(-2.8364137) q[2];
sx q[2];
rz(-1.9532433) q[2];
sx q[2];
rz(0.80079186) q[2];
rz(1.3958037) q[3];
sx q[3];
rz(-2.3621515) q[3];
sx q[3];
rz(3.0099517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
