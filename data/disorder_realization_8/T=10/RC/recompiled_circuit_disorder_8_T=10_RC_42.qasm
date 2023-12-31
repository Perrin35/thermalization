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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5358937) q[0];
sx q[0];
rz(-0.68471691) q[0];
sx q[0];
rz(2.2936355) q[0];
x q[1];
rz(0.71360795) q[2];
sx q[2];
rz(-0.2954233) q[2];
sx q[2];
rz(1.7507391) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9368254) q[1];
sx q[1];
rz(-1.966404) q[1];
sx q[1];
rz(-0.25524615) q[1];
rz(-2.9526887) q[3];
sx q[3];
rz(-1.2652664) q[3];
sx q[3];
rz(-0.73959914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(-0.096244372) q[2];
rz(-1.0359267) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(-2.9878785) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44089833) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(-2.3764215) q[0];
rz(1.2922497) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(0.66295019) q[1];
rz(-pi/2) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-2.2234369) q[2];
sx q[2];
rz(-1.2456206) q[2];
sx q[2];
rz(-2.6174389) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1124277) q[1];
sx q[1];
rz(-1.4883092) q[1];
sx q[1];
rz(-0.44177456) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4088267) q[3];
sx q[3];
rz(-1.0893981) q[3];
sx q[3];
rz(-0.045452047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.048916653) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(2.4760903) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5927785) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(1.0173343) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(2.6229048) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9496574) q[0];
sx q[0];
rz(-1.2149095) q[0];
sx q[0];
rz(0.8215254) q[0];
rz(-pi) q[1];
rz(1.475004) q[2];
sx q[2];
rz(-1.3057858) q[2];
sx q[2];
rz(1.67213) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4244713) q[1];
sx q[1];
rz(-2.3470504) q[1];
sx q[1];
rz(-2.8984927) q[1];
x q[2];
rz(0.014702602) q[3];
sx q[3];
rz(-0.054617453) q[3];
sx q[3];
rz(2.2811449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2805933) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(2.5391501) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-1.4825876) q[1];
sx q[1];
rz(-2.6285016) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4671191) q[0];
sx q[0];
rz(-1.4917443) q[0];
sx q[0];
rz(-2.5208958) q[0];
rz(-pi) q[1];
rz(-2.5330403) q[2];
sx q[2];
rz(-0.90494472) q[2];
sx q[2];
rz(-2.7666639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.40401134) q[1];
sx q[1];
rz(-1.1176795) q[1];
sx q[1];
rz(1.7246507) q[1];
x q[2];
rz(3.1107535) q[3];
sx q[3];
rz(-1.3645002) q[3];
sx q[3];
rz(-1.4269958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6667368) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(2.722548) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(-2.6823147) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6722365) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(0.67681926) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(-0.61606032) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4426684) q[0];
sx q[0];
rz(-1.5183503) q[0];
sx q[0];
rz(3.1057182) q[0];
x q[1];
rz(2.8370503) q[2];
sx q[2];
rz(-0.38123044) q[2];
sx q[2];
rz(0.234137) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38356009) q[1];
sx q[1];
rz(-0.57764232) q[1];
sx q[1];
rz(-3.1314965) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1957937) q[3];
sx q[3];
rz(-2.4295394) q[3];
sx q[3];
rz(-0.28111162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(2.8862254) q[2];
rz(1.5363961) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42246321) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(2.6690924) q[0];
rz(-2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(-2.2568259) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038267604) q[0];
sx q[0];
rz(-1.190435) q[0];
sx q[0];
rz(-3.0618219) q[0];
rz(-pi) q[1];
rz(-2.9315345) q[2];
sx q[2];
rz(-1.5846328) q[2];
sx q[2];
rz(-2.7163598) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5342418) q[1];
sx q[1];
rz(-0.51041767) q[1];
sx q[1];
rz(0.95153248) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24598908) q[3];
sx q[3];
rz(-2.0912598) q[3];
sx q[3];
rz(2.1465079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3970967) q[2];
sx q[2];
rz(-0.8049736) q[2];
sx q[2];
rz(0.3113783) q[2];
rz(1.3686251) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.41806528) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(-3.080522) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(-0.73289245) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7948579) q[0];
sx q[0];
rz(-1.96694) q[0];
sx q[0];
rz(-2.5651155) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0211208) q[2];
sx q[2];
rz(-0.78133821) q[2];
sx q[2];
rz(0.48175016) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1317132) q[1];
sx q[1];
rz(-2.1965346) q[1];
sx q[1];
rz(0.0080607944) q[1];
x q[2];
rz(-3.1002058) q[3];
sx q[3];
rz(-1.4635411) q[3];
sx q[3];
rz(-2.0859352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0682893) q[2];
sx q[2];
rz(-1.1969593) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.085389) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(-0.7094267) q[0];
rz(1.6363232) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(-2.8628796) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73662169) q[0];
sx q[0];
rz(-1.706089) q[0];
sx q[0];
rz(-2.9297329) q[0];
rz(-pi) q[1];
rz(-0.3785554) q[2];
sx q[2];
rz(-0.59213973) q[2];
sx q[2];
rz(0.90422599) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7396001) q[1];
sx q[1];
rz(-2.7584689) q[1];
sx q[1];
rz(0.95462228) q[1];
x q[2];
rz(1.9023499) q[3];
sx q[3];
rz(-2.9058876) q[3];
sx q[3];
rz(-1.9086259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8401106) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(0.85514832) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(-2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99496019) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(0.96889281) q[0];
rz(-2.6682207) q[1];
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
rz(2.5794967) q[0];
sx q[0];
rz(-1.2764494) q[0];
sx q[0];
rz(3.0466945) q[0];
x q[1];
rz(-2.4382486) q[2];
sx q[2];
rz(-1.568927) q[2];
sx q[2];
rz(-1.4700996) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0513623) q[1];
sx q[1];
rz(-1.4453332) q[1];
sx q[1];
rz(-2.5999703) q[1];
rz(-pi) q[2];
rz(1.8248796) q[3];
sx q[3];
rz(-2.4554606) q[3];
sx q[3];
rz(0.38476598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.896495) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(-1.0207821) q[2];
rz(-2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(-0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1616515) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(-0.7014057) q[0];
rz(2.2258863) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(1.9030301) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0093875) q[0];
sx q[0];
rz(-1.7863677) q[0];
sx q[0];
rz(1.042004) q[0];
x q[1];
rz(3.086834) q[2];
sx q[2];
rz(-2.3878532) q[2];
sx q[2];
rz(2.4278305) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4148256) q[1];
sx q[1];
rz(-1.4531724) q[1];
sx q[1];
rz(2.2284501) q[1];
rz(-pi) q[2];
rz(0.978312) q[3];
sx q[3];
rz(-1.6654135) q[3];
sx q[3];
rz(2.832151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4548268) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(-1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(2.1380077) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6702406) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(-3.1162221) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(-2.8364137) q[2];
sx q[2];
rz(-1.9532433) q[2];
sx q[2];
rz(0.80079186) q[2];
rz(1.745789) q[3];
sx q[3];
rz(-0.77944118) q[3];
sx q[3];
rz(-0.13164095) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
