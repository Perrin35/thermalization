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
rz(-1.6725809) q[0];
sx q[0];
rz(-2.2218158) q[0];
sx q[0];
rz(-0.65499175) q[0];
rz(-1.4870149) q[1];
sx q[1];
rz(-1.4227285) q[1];
sx q[1];
rz(1.3230327) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1974063) q[0];
sx q[0];
rz(-1.8776769) q[0];
sx q[0];
rz(-0.57800622) q[0];
rz(2.5147314) q[2];
sx q[2];
rz(-0.38598362) q[2];
sx q[2];
rz(0.71039334) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.536024) q[1];
sx q[1];
rz(-1.6704511) q[1];
sx q[1];
rz(-1.8479363) q[1];
rz(-pi) q[2];
rz(-1.0336797) q[3];
sx q[3];
rz(-0.91712807) q[3];
sx q[3];
rz(-1.9680227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46372867) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(0.65832552) q[2];
rz(-2.5074734) q[3];
sx q[3];
rz(-1.442797) q[3];
sx q[3];
rz(-1.770795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.051801) q[0];
sx q[0];
rz(-2.8103516) q[0];
sx q[0];
rz(0.077089699) q[0];
rz(-0.58473051) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(-1.9416521) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3886007) q[0];
sx q[0];
rz(-1.3738828) q[0];
sx q[0];
rz(1.1042751) q[0];
rz(-pi) q[1];
rz(-1.7464094) q[2];
sx q[2];
rz(-2.9746911) q[2];
sx q[2];
rz(1.4457955) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6226095) q[1];
sx q[1];
rz(-2.4138192) q[1];
sx q[1];
rz(-0.98811291) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6317815) q[3];
sx q[3];
rz(-2.1290696) q[3];
sx q[3];
rz(1.3577485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8889019) q[2];
sx q[2];
rz(-2.1609047) q[2];
sx q[2];
rz(-0.67548951) q[2];
rz(-3.1318956) q[3];
sx q[3];
rz(-0.24295013) q[3];
sx q[3];
rz(-3.1233136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6460687) q[0];
sx q[0];
rz(-1.4718453) q[0];
sx q[0];
rz(2.38696) q[0];
rz(3.1309639) q[1];
sx q[1];
rz(-1.7517215) q[1];
sx q[1];
rz(-1.0108112) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1026238) q[0];
sx q[0];
rz(-2.0576697) q[0];
sx q[0];
rz(-0.34619934) q[0];
rz(3.1132142) q[2];
sx q[2];
rz(-2.252223) q[2];
sx q[2];
rz(0.703237) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7422223) q[1];
sx q[1];
rz(-2.3934747) q[1];
sx q[1];
rz(-3.0518603) q[1];
rz(-pi) q[2];
rz(-2.1819886) q[3];
sx q[3];
rz(-0.90590796) q[3];
sx q[3];
rz(-1.5184107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2406771) q[2];
sx q[2];
rz(-0.42833504) q[2];
sx q[2];
rz(0.50673103) q[2];
rz(1.2725376) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(-2.8633964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5946567) q[0];
sx q[0];
rz(-1.9173859) q[0];
sx q[0];
rz(-2.0261672) q[0];
rz(1.2660654) q[1];
sx q[1];
rz(-1.5770117) q[1];
sx q[1];
rz(-0.87475264) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6653671) q[0];
sx q[0];
rz(-1.3194808) q[0];
sx q[0];
rz(0.091500207) q[0];
x q[1];
rz(-1.546193) q[2];
sx q[2];
rz(-2.4877791) q[2];
sx q[2];
rz(-0.023651274) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1072448) q[1];
sx q[1];
rz(-2.6470091) q[1];
sx q[1];
rz(0.58652189) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7797919) q[3];
sx q[3];
rz(-2.4242988) q[3];
sx q[3];
rz(1.532326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11883417) q[2];
sx q[2];
rz(-0.47670445) q[2];
sx q[2];
rz(1.9258707) q[2];
rz(-1.7440965) q[3];
sx q[3];
rz(-1.5236866) q[3];
sx q[3];
rz(1.4309179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8656411) q[0];
sx q[0];
rz(-3.0396099) q[0];
sx q[0];
rz(1.6708466) q[0];
rz(1.3784846) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(1.4102304) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3138235) q[0];
sx q[0];
rz(-2.6010989) q[0];
sx q[0];
rz(1.1590411) q[0];
x q[1];
rz(-1.8115787) q[2];
sx q[2];
rz(-0.91139873) q[2];
sx q[2];
rz(-1.8383775) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.901625) q[1];
sx q[1];
rz(-0.31244597) q[1];
sx q[1];
rz(0.0052253763) q[1];
rz(-pi) q[2];
rz(2.2447884) q[3];
sx q[3];
rz(-0.53433472) q[3];
sx q[3];
rz(-0.33269445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1269647) q[2];
sx q[2];
rz(-2.2829599) q[2];
sx q[2];
rz(-1.3237759) q[2];
rz(-1.3592367) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112261) q[0];
sx q[0];
rz(-1.4596326) q[0];
sx q[0];
rz(3.1193745) q[0];
rz(-0.70023316) q[1];
sx q[1];
rz(-1.2513688) q[1];
sx q[1];
rz(-0.95692316) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6895804) q[0];
sx q[0];
rz(-2.4274024) q[0];
sx q[0];
rz(-0.22269188) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3889203) q[2];
sx q[2];
rz(-1.3564912) q[2];
sx q[2];
rz(0.3265115) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14084223) q[1];
sx q[1];
rz(-0.4943119) q[1];
sx q[1];
rz(0.96338455) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33958667) q[3];
sx q[3];
rz(-1.1533779) q[3];
sx q[3];
rz(-1.0610895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2172829) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(1.4420606) q[2];
rz(-2.0349272) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(1.1520011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7149413) q[0];
sx q[0];
rz(-1.8639257) q[0];
sx q[0];
rz(0.46698025) q[0];
rz(-1.8309719) q[1];
sx q[1];
rz(-2.1881723) q[1];
sx q[1];
rz(-0.39042815) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6990427) q[0];
sx q[0];
rz(-2.9973898) q[0];
sx q[0];
rz(-0.22871916) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0340334) q[2];
sx q[2];
rz(-1.8637519) q[2];
sx q[2];
rz(0.40574542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5891287) q[1];
sx q[1];
rz(-2.2099582) q[1];
sx q[1];
rz(-2.5659849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14777811) q[3];
sx q[3];
rz(-0.96333233) q[3];
sx q[3];
rz(1.0318508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.19162576) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(2.1245655) q[2];
rz(-0.93938604) q[3];
sx q[3];
rz(-0.87299577) q[3];
sx q[3];
rz(-0.92923195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4262714) q[0];
sx q[0];
rz(-2.4763698) q[0];
sx q[0];
rz(1.2128879) q[0];
rz(-0.60316482) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(2.718198) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7170537) q[0];
sx q[0];
rz(-1.6485212) q[0];
sx q[0];
rz(-2.2721798) q[0];
x q[1];
rz(-1.4199443) q[2];
sx q[2];
rz(-1.4948339) q[2];
sx q[2];
rz(1.1830038) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.054476995) q[1];
sx q[1];
rz(-1.7810139) q[1];
sx q[1];
rz(1.3940548) q[1];
x q[2];
rz(1.5782062) q[3];
sx q[3];
rz(-1.4009589) q[3];
sx q[3];
rz(1.8078723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.49070552) q[2];
sx q[2];
rz(-2.8287973) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097718358) q[0];
sx q[0];
rz(-0.23625034) q[0];
sx q[0];
rz(-3.0978715) q[0];
rz(1.9678736) q[1];
sx q[1];
rz(-0.81467384) q[1];
sx q[1];
rz(-0.79183212) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39003644) q[0];
sx q[0];
rz(-1.1564213) q[0];
sx q[0];
rz(-0.66585559) q[0];
rz(-pi) q[1];
rz(-0.20251198) q[2];
sx q[2];
rz(-1.8648476) q[2];
sx q[2];
rz(-2.685315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89347311) q[1];
sx q[1];
rz(-1.1482825) q[1];
sx q[1];
rz(0.36568195) q[1];
x q[2];
rz(-2.6641162) q[3];
sx q[3];
rz(-1.9100185) q[3];
sx q[3];
rz(2.4216258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9110079) q[2];
sx q[2];
rz(-0.30969301) q[2];
sx q[2];
rz(-1.9045551) q[2];
rz(-2.6593995) q[3];
sx q[3];
rz(-2.2174406) q[3];
sx q[3];
rz(-2.6528416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.024254) q[0];
sx q[0];
rz(-2.2145705) q[0];
sx q[0];
rz(1.3379958) q[0];
rz(0.85211873) q[1];
sx q[1];
rz(-1.9133277) q[1];
sx q[1];
rz(-1.9911912) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.710708) q[0];
sx q[0];
rz(-2.0938494) q[0];
sx q[0];
rz(2.5780748) q[0];
rz(1.8878292) q[2];
sx q[2];
rz(-1.46508) q[2];
sx q[2];
rz(2.6144165) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4070421) q[1];
sx q[1];
rz(-1.4311387) q[1];
sx q[1];
rz(-1.1753083) q[1];
rz(-2.1701067) q[3];
sx q[3];
rz(-1.8953634) q[3];
sx q[3];
rz(0.49719405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77356768) q[2];
sx q[2];
rz(-1.0057534) q[2];
sx q[2];
rz(-2.8311938) q[2];
rz(-0.6271022) q[3];
sx q[3];
rz(-2.1576364) q[3];
sx q[3];
rz(1.9703126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8695759) q[0];
sx q[0];
rz(-1.5675114) q[0];
sx q[0];
rz(-1.5935224) q[0];
rz(1.2177474) q[1];
sx q[1];
rz(-1.6883313) q[1];
sx q[1];
rz(2.2399166) q[1];
rz(-2.9318342) q[2];
sx q[2];
rz(-2.0611613) q[2];
sx q[2];
rz(0.47917889) q[2];
rz(-1.6321833) q[3];
sx q[3];
rz(-1.2493549) q[3];
sx q[3];
rz(-3.1102104) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
