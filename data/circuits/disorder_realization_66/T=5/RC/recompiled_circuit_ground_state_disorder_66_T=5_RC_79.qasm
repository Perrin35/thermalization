OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.42216766) q[0];
sx q[0];
rz(7.0081975) q[0];
sx q[0];
rz(10.260697) q[0];
rz(0.82756502) q[1];
sx q[1];
rz(-2.8711072) q[1];
sx q[1];
rz(-0.40721133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4370928) q[0];
sx q[0];
rz(-1.8665534) q[0];
sx q[0];
rz(2.9192135) q[0];
rz(-pi) q[1];
rz(2.902342) q[2];
sx q[2];
rz(-1.2206589) q[2];
sx q[2];
rz(0.020933271) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5558124) q[1];
sx q[1];
rz(-1.4364873) q[1];
sx q[1];
rz(2.7463169) q[1];
rz(-pi) q[2];
rz(-0.9933957) q[3];
sx q[3];
rz(-0.90809408) q[3];
sx q[3];
rz(-0.71383324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.779458) q[2];
sx q[2];
rz(-0.3388277) q[2];
sx q[2];
rz(-0.96620488) q[2];
rz(-0.90159121) q[3];
sx q[3];
rz(-0.85351557) q[3];
sx q[3];
rz(1.0454267) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7411984) q[0];
sx q[0];
rz(-2.5303831) q[0];
sx q[0];
rz(0.10261593) q[0];
rz(-1.1688894) q[1];
sx q[1];
rz(-2.1714307) q[1];
sx q[1];
rz(2.6279367) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9540323) q[0];
sx q[0];
rz(-0.49278773) q[0];
sx q[0];
rz(1.378118) q[0];
rz(2.6046126) q[2];
sx q[2];
rz(-1.6936079) q[2];
sx q[2];
rz(1.7857743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8674492) q[1];
sx q[1];
rz(-1.2701057) q[1];
sx q[1];
rz(-0.97598981) q[1];
x q[2];
rz(1.1369785) q[3];
sx q[3];
rz(-1.2946715) q[3];
sx q[3];
rz(1.2009337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1932842) q[2];
sx q[2];
rz(-2.9479492) q[2];
sx q[2];
rz(-0.2085169) q[2];
rz(2.4559313) q[3];
sx q[3];
rz(-2.0833569) q[3];
sx q[3];
rz(-2.9764777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.044428069) q[0];
sx q[0];
rz(-1.8655638) q[0];
sx q[0];
rz(-1.9097419) q[0];
rz(-2.7992898) q[1];
sx q[1];
rz(-1.1016223) q[1];
sx q[1];
rz(1.5493578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4195815) q[0];
sx q[0];
rz(-1.3886459) q[0];
sx q[0];
rz(-2.8246882) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6913271) q[2];
sx q[2];
rz(-1.4107586) q[2];
sx q[2];
rz(-1.1628422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6173319) q[1];
sx q[1];
rz(-2.2843505) q[1];
sx q[1];
rz(-0.18194845) q[1];
rz(-pi) q[2];
rz(0.83264787) q[3];
sx q[3];
rz(-0.73048985) q[3];
sx q[3];
rz(2.5363142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0478829) q[2];
sx q[2];
rz(-1.1992531) q[2];
sx q[2];
rz(2.1236146) q[2];
rz(-1.4114981) q[3];
sx q[3];
rz(-0.74677765) q[3];
sx q[3];
rz(1.6104893) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8133076) q[0];
sx q[0];
rz(-1.3628553) q[0];
sx q[0];
rz(-0.12411975) q[0];
rz(2.2165551) q[1];
sx q[1];
rz(-1.3003474) q[1];
sx q[1];
rz(-2.4031924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3401854) q[0];
sx q[0];
rz(-2.3871195) q[0];
sx q[0];
rz(1.6615926) q[0];
rz(-2.5321711) q[2];
sx q[2];
rz(-1.2487234) q[2];
sx q[2];
rz(0.24064482) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.55220375) q[1];
sx q[1];
rz(-0.60221803) q[1];
sx q[1];
rz(2.7715384) q[1];
rz(-1.4703512) q[3];
sx q[3];
rz(-2.1444706) q[3];
sx q[3];
rz(0.84000444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.063252123) q[2];
sx q[2];
rz(-1.1880778) q[2];
sx q[2];
rz(1.1032907) q[2];
rz(2.8830146) q[3];
sx q[3];
rz(-2.0662112) q[3];
sx q[3];
rz(-2.8598089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2461808) q[0];
sx q[0];
rz(-0.87765944) q[0];
sx q[0];
rz(1.8038764) q[0];
rz(2.5216263) q[1];
sx q[1];
rz(-1.679531) q[1];
sx q[1];
rz(1.6372797) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65101609) q[0];
sx q[0];
rz(-1.2913112) q[0];
sx q[0];
rz(-2.9453251) q[0];
rz(0.049311801) q[2];
sx q[2];
rz(-2.4386536) q[2];
sx q[2];
rz(2.9300323) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8057509) q[1];
sx q[1];
rz(-2.1349135) q[1];
sx q[1];
rz(1.1104904) q[1];
rz(-2.3035731) q[3];
sx q[3];
rz(-1.4921643) q[3];
sx q[3];
rz(1.3669612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6164246) q[2];
sx q[2];
rz(-1.9889571) q[2];
sx q[2];
rz(-1.1241414) q[2];
rz(0.21640402) q[3];
sx q[3];
rz(-0.83306044) q[3];
sx q[3];
rz(2.0719349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.44806099) q[0];
sx q[0];
rz(-0.8140642) q[0];
sx q[0];
rz(2.6779209) q[0];
rz(-3.0060153) q[1];
sx q[1];
rz(-0.99628535) q[1];
sx q[1];
rz(-0.3840951) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3883108) q[0];
sx q[0];
rz(-2.658591) q[0];
sx q[0];
rz(0.24899434) q[0];
rz(-pi) q[1];
rz(-2.9303418) q[2];
sx q[2];
rz(-2.3659035) q[2];
sx q[2];
rz(0.83310177) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1138224) q[1];
sx q[1];
rz(-0.40317391) q[1];
sx q[1];
rz(-2.0084963) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0657173) q[3];
sx q[3];
rz(-1.0238227) q[3];
sx q[3];
rz(1.8219655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42528459) q[2];
sx q[2];
rz(-0.3766489) q[2];
sx q[2];
rz(-0.46871218) q[2];
rz(2.5675755) q[3];
sx q[3];
rz(-1.6711957) q[3];
sx q[3];
rz(-2.4911657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.9531517) q[0];
sx q[0];
rz(-1.3583536) q[0];
sx q[0];
rz(0.69001946) q[0];
rz(-2.4192339) q[1];
sx q[1];
rz(-2.0004309) q[1];
sx q[1];
rz(2.3648327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.748257) q[0];
sx q[0];
rz(-1.3384755) q[0];
sx q[0];
rz(1.4938225) q[0];
rz(-pi) q[1];
rz(0.52771126) q[2];
sx q[2];
rz(-0.65831682) q[2];
sx q[2];
rz(1.1084313) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.986514) q[1];
sx q[1];
rz(-2.6342794) q[1];
sx q[1];
rz(-2.581937) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5473821) q[3];
sx q[3];
rz(-0.6802313) q[3];
sx q[3];
rz(1.3772723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7338099) q[2];
sx q[2];
rz(-1.5668543) q[2];
sx q[2];
rz(-3.0877647) q[2];
rz(2.5194061) q[3];
sx q[3];
rz(-0.48431188) q[3];
sx q[3];
rz(-3.0622862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6513885) q[0];
sx q[0];
rz(-1.9909998) q[0];
sx q[0];
rz(-1.401061) q[0];
rz(2.1531064) q[1];
sx q[1];
rz(-1.8732312) q[1];
sx q[1];
rz(-0.40774694) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4238472) q[0];
sx q[0];
rz(-1.9739986) q[0];
sx q[0];
rz(-1.3333596) q[0];
rz(-pi) q[1];
rz(2.2456535) q[2];
sx q[2];
rz(-1.6701588) q[2];
sx q[2];
rz(2.4473913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11800471) q[1];
sx q[1];
rz(-1.1049926) q[1];
sx q[1];
rz(-2.3789469) q[1];
x q[2];
rz(-1.7090104) q[3];
sx q[3];
rz(-1.8402037) q[3];
sx q[3];
rz(0.22002735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72655788) q[2];
sx q[2];
rz(-2.0099535) q[2];
sx q[2];
rz(0.69990194) q[2];
rz(-1.0837726) q[3];
sx q[3];
rz(-2.2743069) q[3];
sx q[3];
rz(-2.9054902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20339762) q[0];
sx q[0];
rz(-1.6827826) q[0];
sx q[0];
rz(0.38988018) q[0];
rz(-1.1979206) q[1];
sx q[1];
rz(-0.094466297) q[1];
sx q[1];
rz(0.87497154) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41738865) q[0];
sx q[0];
rz(-1.5932788) q[0];
sx q[0];
rz(2.1940986) q[0];
rz(-pi) q[1];
rz(2.2570043) q[2];
sx q[2];
rz(-1.8914127) q[2];
sx q[2];
rz(0.55701643) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5510721) q[1];
sx q[1];
rz(-1.1391907) q[1];
sx q[1];
rz(1.9533532) q[1];
rz(-pi) q[2];
rz(-2.1375108) q[3];
sx q[3];
rz(-1.9329786) q[3];
sx q[3];
rz(-3.0423861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8209057) q[2];
sx q[2];
rz(-2.4497538) q[2];
sx q[2];
rz(0.64986491) q[2];
rz(-1.1167022) q[3];
sx q[3];
rz(-1.2973123) q[3];
sx q[3];
rz(2.9445649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7277302) q[0];
sx q[0];
rz(-0.05302269) q[0];
sx q[0];
rz(2.9470288) q[0];
rz(2.3692756) q[1];
sx q[1];
rz(-1.1338502) q[1];
sx q[1];
rz(-0.17759855) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1188162) q[0];
sx q[0];
rz(-1.5714065) q[0];
sx q[0];
rz(1.5691891) q[0];
rz(-2.3760711) q[2];
sx q[2];
rz(-2.3280848) q[2];
sx q[2];
rz(3.0610656) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0760484) q[1];
sx q[1];
rz(-2.0418344) q[1];
sx q[1];
rz(-0.7008497) q[1];
rz(-pi) q[2];
rz(0.94915139) q[3];
sx q[3];
rz(-0.92499925) q[3];
sx q[3];
rz(2.8634255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.78581587) q[2];
sx q[2];
rz(-0.59712258) q[2];
sx q[2];
rz(0.43006483) q[2];
rz(2.0307342) q[3];
sx q[3];
rz(-1.6836124) q[3];
sx q[3];
rz(-2.7250169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149902) q[0];
sx q[0];
rz(-1.5530598) q[0];
sx q[0];
rz(-0.1834827) q[0];
rz(1.968374) q[1];
sx q[1];
rz(-1.219974) q[1];
sx q[1];
rz(-3.0164607) q[1];
rz(2.0963893) q[2];
sx q[2];
rz(-1.3202616) q[2];
sx q[2];
rz(3.1332982) q[2];
rz(-2.2852702) q[3];
sx q[3];
rz(-1.1416937) q[3];
sx q[3];
rz(1.0037224) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
