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
rz(1.4690118) q[0];
sx q[0];
rz(-0.91977683) q[0];
sx q[0];
rz(0.65499175) q[0];
rz(-1.4870149) q[1];
sx q[1];
rz(-1.4227285) q[1];
sx q[1];
rz(1.3230327) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9441863) q[0];
sx q[0];
rz(-1.8776769) q[0];
sx q[0];
rz(-2.5635864) q[0];
rz(-pi) q[1];
rz(1.8048067) q[2];
sx q[2];
rz(-1.8806226) q[2];
sx q[2];
rz(1.7676304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2046584) q[1];
sx q[1];
rz(-1.2950674) q[1];
sx q[1];
rz(0.10358056) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4134656) q[3];
sx q[3];
rz(-1.1525197) q[3];
sx q[3];
rz(-0.74467105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.677864) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(-2.4832671) q[2];
rz(-2.5074734) q[3];
sx q[3];
rz(-1.442797) q[3];
sx q[3];
rz(1.3707976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.051801) q[0];
sx q[0];
rz(-2.8103516) q[0];
sx q[0];
rz(-0.077089699) q[0];
rz(0.58473051) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(-1.1999406) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532419) q[0];
sx q[0];
rz(-2.6380499) q[0];
sx q[0];
rz(-1.153323) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7464094) q[2];
sx q[2];
rz(-2.9746911) q[2];
sx q[2];
rz(1.4457955) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6226095) q[1];
sx q[1];
rz(-0.72777343) q[1];
sx q[1];
rz(-2.1534797) q[1];
x q[2];
rz(-1.6317815) q[3];
sx q[3];
rz(-2.1290696) q[3];
sx q[3];
rz(1.3577485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8889019) q[2];
sx q[2];
rz(-2.1609047) q[2];
sx q[2];
rz(-0.67548951) q[2];
rz(0.0096970079) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(3.1233136) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460687) q[0];
sx q[0];
rz(-1.6697474) q[0];
sx q[0];
rz(2.38696) q[0];
rz(-3.1309639) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(2.1307814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36463144) q[0];
sx q[0];
rz(-1.2662132) q[0];
sx q[0];
rz(-2.0833894) q[0];
x q[1];
rz(-1.5358244) q[2];
sx q[2];
rz(-0.681923) q[2];
sx q[2];
rz(2.4833895) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6200808) q[1];
sx q[1];
rz(-2.3151868) q[1];
sx q[1];
rz(-1.6537731) q[1];
rz(0.63186462) q[3];
sx q[3];
rz(-2.2711828) q[3];
sx q[3];
rz(0.66955843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2406771) q[2];
sx q[2];
rz(-0.42833504) q[2];
sx q[2];
rz(0.50673103) q[2];
rz(-1.2725376) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(-0.27819628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.5946567) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(-2.0261672) q[0];
rz(-1.8755272) q[1];
sx q[1];
rz(-1.5645809) q[1];
sx q[1];
rz(-2.26684) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0698358) q[0];
sx q[0];
rz(-1.6594145) q[0];
sx q[0];
rz(1.3184692) q[0];
x q[1];
rz(3.1227448) q[2];
sx q[2];
rz(-2.224378) q[2];
sx q[2];
rz(-3.1342521) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.529244) q[1];
sx q[1];
rz(-1.9772286) q[1];
sx q[1];
rz(1.8608577) q[1];
x q[2];
rz(2.277209) q[3];
sx q[3];
rz(-1.4339851) q[3];
sx q[3];
rz(2.9446326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0227585) q[2];
sx q[2];
rz(-0.47670445) q[2];
sx q[2];
rz(-1.215722) q[2];
rz(1.3974961) q[3];
sx q[3];
rz(-1.5236866) q[3];
sx q[3];
rz(-1.7106748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8656411) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(1.4707461) q[0];
rz(-1.7631081) q[1];
sx q[1];
rz(-0.78796092) q[1];
sx q[1];
rz(1.7313622) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61527644) q[0];
sx q[0];
rz(-1.363376) q[0];
sx q[0];
rz(2.0735969) q[0];
rz(0.29844923) q[2];
sx q[2];
rz(-2.4458234) q[2];
sx q[2];
rz(-1.6844105) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.901625) q[1];
sx q[1];
rz(-2.8291467) q[1];
sx q[1];
rz(-3.1363673) q[1];
rz(-pi) q[2];
rz(2.7878109) q[3];
sx q[3];
rz(-1.1615586) q[3];
sx q[3];
rz(-1.080846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1269647) q[2];
sx q[2];
rz(-0.85863272) q[2];
sx q[2];
rz(1.8178168) q[2];
rz(1.7823559) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(2.7544379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1112261) q[0];
sx q[0];
rz(-1.68196) q[0];
sx q[0];
rz(3.1193745) q[0];
rz(0.70023316) q[1];
sx q[1];
rz(-1.2513688) q[1];
sx q[1];
rz(0.95692316) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6895804) q[0];
sx q[0];
rz(-2.4274024) q[0];
sx q[0];
rz(0.22269188) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9238052) q[2];
sx q[2];
rz(-1.3931257) q[2];
sx q[2];
rz(1.2051932) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14084223) q[1];
sx q[1];
rz(-0.4943119) q[1];
sx q[1];
rz(0.96338455) q[1];
rz(-2.2150008) q[3];
sx q[3];
rz(-0.53172382) q[3];
sx q[3];
rz(-2.7972972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2172829) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(1.699532) q[2];
rz(-1.1066655) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(-1.1520011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42665136) q[0];
sx q[0];
rz(-1.8639257) q[0];
sx q[0];
rz(-0.46698025) q[0];
rz(-1.3106208) q[1];
sx q[1];
rz(-0.95342031) q[1];
sx q[1];
rz(2.7511645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6990427) q[0];
sx q[0];
rz(-2.9973898) q[0];
sx q[0];
rz(0.22871916) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2288741) q[2];
sx q[2];
rz(-2.83005) q[2];
sx q[2];
rz(2.3780413) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5891287) q[1];
sx q[1];
rz(-2.2099582) q[1];
sx q[1];
rz(2.5659849) q[1];
rz(0.95818635) q[3];
sx q[3];
rz(-1.691992) q[3];
sx q[3];
rz(2.6874128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9499669) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(2.1245655) q[2];
rz(2.2022066) q[3];
sx q[3];
rz(-2.2685969) q[3];
sx q[3];
rz(-2.2123607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7153213) q[0];
sx q[0];
rz(-0.66522288) q[0];
sx q[0];
rz(-1.9287047) q[0];
rz(0.60316482) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(-2.718198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2379266) q[0];
sx q[0];
rz(-2.4366424) q[0];
sx q[0];
rz(1.6909107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1020057) q[2];
sx q[2];
rz(-2.9728242) q[2];
sx q[2];
rz(0.075254879) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4881542) q[1];
sx q[1];
rz(-0.27380015) q[1];
sx q[1];
rz(0.68922148) q[1];
rz(-pi) q[2];
rz(0.16984197) q[3];
sx q[3];
rz(-1.5780996) q[3];
sx q[3];
rz(0.23582349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6508871) q[2];
sx q[2];
rz(-0.31279534) q[2];
sx q[2];
rz(-0.97839626) q[2];
rz(2.4466416) q[3];
sx q[3];
rz(-1.2000822) q[3];
sx q[3];
rz(1.2945226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0438743) q[0];
sx q[0];
rz(-2.9053423) q[0];
sx q[0];
rz(0.043721113) q[0];
rz(1.1737191) q[1];
sx q[1];
rz(-2.3269188) q[1];
sx q[1];
rz(-0.79183212) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39003644) q[0];
sx q[0];
rz(-1.1564213) q[0];
sx q[0];
rz(-0.66585559) q[0];
rz(-pi) q[1];
rz(2.1570683) q[2];
sx q[2];
rz(-0.35536623) q[2];
sx q[2];
rz(-2.9815069) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83306317) q[1];
sx q[1];
rz(-1.2385784) q[1];
sx q[1];
rz(2.019472) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1926226) q[3];
sx q[3];
rz(-2.0190051) q[3];
sx q[3];
rz(2.1202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9110079) q[2];
sx q[2];
rz(-2.8318996) q[2];
sx q[2];
rz(-1.2370375) q[2];
rz(-0.48219314) q[3];
sx q[3];
rz(-2.2174406) q[3];
sx q[3];
rz(2.6528416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.024254) q[0];
sx q[0];
rz(-0.92702213) q[0];
sx q[0];
rz(-1.3379958) q[0];
rz(0.85211873) q[1];
sx q[1];
rz(-1.9133277) q[1];
sx q[1];
rz(-1.9911912) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.710708) q[0];
sx q[0];
rz(-1.0477433) q[0];
sx q[0];
rz(2.5780748) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.030376) q[2];
sx q[2];
rz(-1.885998) q[2];
sx q[2];
rz(1.0090131) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89429606) q[1];
sx q[1];
rz(-1.9622231) q[1];
sx q[1];
rz(-0.15116926) q[1];
x q[2];
rz(-0.97148599) q[3];
sx q[3];
rz(-1.8953634) q[3];
sx q[3];
rz(2.6443986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.368025) q[2];
sx q[2];
rz(-2.1358392) q[2];
sx q[2];
rz(-0.31039882) q[2];
rz(0.6271022) q[3];
sx q[3];
rz(-0.98395625) q[3];
sx q[3];
rz(-1.1712801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8695759) q[0];
sx q[0];
rz(-1.5675114) q[0];
sx q[0];
rz(-1.5935224) q[0];
rz(-1.9238453) q[1];
sx q[1];
rz(-1.6883313) q[1];
sx q[1];
rz(2.2399166) q[1];
rz(-0.2097585) q[2];
sx q[2];
rz(-1.0804313) q[2];
sx q[2];
rz(-2.6624138) q[2];
rz(0.32200702) q[3];
sx q[3];
rz(-1.5125572) q[3];
sx q[3];
rz(-1.5199979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
