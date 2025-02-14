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
rz(-2.2313843) q[0];
sx q[0];
rz(-2.3910523) q[0];
sx q[0];
rz(2.5810177) q[0];
rz(1.9380467) q[1];
sx q[1];
rz(-2.7653341) q[1];
sx q[1];
rz(2.9484152) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5995418) q[0];
sx q[0];
rz(-2.1452869) q[0];
sx q[0];
rz(-2.1035139) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2760787) q[2];
sx q[2];
rz(-1.4774472) q[2];
sx q[2];
rz(0.35356086) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5467599) q[1];
sx q[1];
rz(-1.7103737) q[1];
sx q[1];
rz(0.91800331) q[1];
rz(-0.47115342) q[3];
sx q[3];
rz(-2.7756925) q[3];
sx q[3];
rz(0.19972502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2304077) q[2];
sx q[2];
rz(-0.62778968) q[2];
sx q[2];
rz(-2.5837303) q[2];
rz(1.2281536) q[3];
sx q[3];
rz(-1.4598673) q[3];
sx q[3];
rz(-0.094999878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.2896344) q[0];
sx q[0];
rz(-1.8486706) q[0];
sx q[0];
rz(1.5767545) q[0];
rz(1.9948888) q[1];
sx q[1];
rz(-1.4253989) q[1];
sx q[1];
rz(-1.0316521) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080650157) q[0];
sx q[0];
rz(-2.553078) q[0];
sx q[0];
rz(2.0015765) q[0];
rz(2.7286058) q[2];
sx q[2];
rz(-1.3527762) q[2];
sx q[2];
rz(0.08139164) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.90742753) q[1];
sx q[1];
rz(-1.2303367) q[1];
sx q[1];
rz(0.84788604) q[1];
rz(1.7748717) q[3];
sx q[3];
rz(-2.0672247) q[3];
sx q[3];
rz(-1.6599409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8138294) q[2];
sx q[2];
rz(-0.29080614) q[2];
sx q[2];
rz(1.4168868) q[2];
rz(-2.318577) q[3];
sx q[3];
rz(-1.6133285) q[3];
sx q[3];
rz(0.64457646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3642035) q[0];
sx q[0];
rz(-1.1193898) q[0];
sx q[0];
rz(0.35225824) q[0];
rz(2.7525355) q[1];
sx q[1];
rz(-1.9720826) q[1];
sx q[1];
rz(-2.9152117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.987623) q[0];
sx q[0];
rz(-2.8806318) q[0];
sx q[0];
rz(-1.8960192) q[0];
rz(-2.1445072) q[2];
sx q[2];
rz(-1.8663532) q[2];
sx q[2];
rz(1.3546582) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55646642) q[1];
sx q[1];
rz(-2.0279998) q[1];
sx q[1];
rz(-0.5208421) q[1];
x q[2];
rz(-2.705615) q[3];
sx q[3];
rz(-1.5020026) q[3];
sx q[3];
rz(1.1747557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9809197) q[2];
sx q[2];
rz(-2.1143165) q[2];
sx q[2];
rz(-0.36160198) q[2];
rz(-2.8412039) q[3];
sx q[3];
rz(-1.8301679) q[3];
sx q[3];
rz(0.96295199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.921973) q[0];
sx q[0];
rz(-2.0491056) q[0];
sx q[0];
rz(0.12119448) q[0];
rz(1.1990064) q[1];
sx q[1];
rz(-1.0795178) q[1];
sx q[1];
rz(1.8669063) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.629144) q[0];
sx q[0];
rz(-1.7928018) q[0];
sx q[0];
rz(-1.52284) q[0];
x q[1];
rz(1.5869467) q[2];
sx q[2];
rz(-1.776219) q[2];
sx q[2];
rz(1.1753163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.53654799) q[1];
sx q[1];
rz(-1.0031317) q[1];
sx q[1];
rz(-1.8325388) q[1];
x q[2];
rz(-1.937247) q[3];
sx q[3];
rz(-2.7525558) q[3];
sx q[3];
rz(-1.4938482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.049383) q[2];
sx q[2];
rz(-1.694724) q[2];
sx q[2];
rz(-1.3775187) q[2];
rz(-1.1285909) q[3];
sx q[3];
rz(-0.94213525) q[3];
sx q[3];
rz(-2.3142864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2500797) q[0];
sx q[0];
rz(-2.2703607) q[0];
sx q[0];
rz(-2.0569892) q[0];
rz(2.4705823) q[1];
sx q[1];
rz(-1.5120466) q[1];
sx q[1];
rz(3.0674518) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4808124) q[0];
sx q[0];
rz(-2.5044821) q[0];
sx q[0];
rz(0.47218944) q[0];
rz(-pi) q[1];
rz(-2.3099075) q[2];
sx q[2];
rz(-2.7234969) q[2];
sx q[2];
rz(1.8036133) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.74846327) q[1];
sx q[1];
rz(-0.65915758) q[1];
sx q[1];
rz(2.0214969) q[1];
rz(-pi) q[2];
rz(-1.2957786) q[3];
sx q[3];
rz(-1.4396808) q[3];
sx q[3];
rz(1.0594378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.025617754) q[2];
sx q[2];
rz(-2.3927549) q[2];
sx q[2];
rz(-0.96308127) q[2];
rz(1.7668308) q[3];
sx q[3];
rz(-1.129496) q[3];
sx q[3];
rz(-0.53205427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2087723) q[0];
sx q[0];
rz(-0.29173478) q[0];
sx q[0];
rz(-1.8836841) q[0];
rz(2.3014297) q[1];
sx q[1];
rz(-1.1779307) q[1];
sx q[1];
rz(-2.449583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0449311) q[0];
sx q[0];
rz(-0.4293483) q[0];
sx q[0];
rz(-2.2675687) q[0];
rz(-pi) q[1];
rz(0.76592036) q[2];
sx q[2];
rz(-1.7007075) q[2];
sx q[2];
rz(-0.15405642) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4547448) q[1];
sx q[1];
rz(-2.1917289) q[1];
sx q[1];
rz(-2.5140155) q[1];
rz(-pi) q[2];
rz(-2.1282084) q[3];
sx q[3];
rz(-1.0104826) q[3];
sx q[3];
rz(-0.91065613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.17795263) q[2];
sx q[2];
rz(-1.0963564) q[2];
sx q[2];
rz(-0.5160416) q[2];
rz(1.7959203) q[3];
sx q[3];
rz(-0.44107744) q[3];
sx q[3];
rz(2.1399982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8733785) q[0];
sx q[0];
rz(-0.28609797) q[0];
sx q[0];
rz(1.0299261) q[0];
rz(0.65451199) q[1];
sx q[1];
rz(-1.8067358) q[1];
sx q[1];
rz(-0.17394224) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3315863) q[0];
sx q[0];
rz(-2.2751006) q[0];
sx q[0];
rz(-2.293701) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1047178) q[2];
sx q[2];
rz(-0.838111) q[2];
sx q[2];
rz(-1.378841) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6831931) q[1];
sx q[1];
rz(-1.8432861) q[1];
sx q[1];
rz(-2.2098455) q[1];
rz(2.6264013) q[3];
sx q[3];
rz(-1.8015007) q[3];
sx q[3];
rz(-2.9700792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8339463) q[2];
sx q[2];
rz(-0.53457326) q[2];
sx q[2];
rz(-1.3391116) q[2];
rz(-2.3329959) q[3];
sx q[3];
rz(-1.9924889) q[3];
sx q[3];
rz(-1.8554525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.9209552) q[0];
sx q[0];
rz(-1.0458825) q[0];
sx q[0];
rz(-1.9761696) q[0];
rz(-2.9479345) q[1];
sx q[1];
rz(-0.8332738) q[1];
sx q[1];
rz(-1.1509034) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5315463) q[0];
sx q[0];
rz(-0.70554107) q[0];
sx q[0];
rz(-0.39682589) q[0];
x q[1];
rz(1.3940349) q[2];
sx q[2];
rz(-2.3675248) q[2];
sx q[2];
rz(1.7947527) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.47025679) q[1];
sx q[1];
rz(-1.0733979) q[1];
sx q[1];
rz(2.6543544) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99566136) q[3];
sx q[3];
rz(-1.0049977) q[3];
sx q[3];
rz(-3.107323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9264441) q[2];
sx q[2];
rz(-2.1200659) q[2];
sx q[2];
rz(1.8683757) q[2];
rz(-0.5145973) q[3];
sx q[3];
rz(-0.319258) q[3];
sx q[3];
rz(0.74015051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99106818) q[0];
sx q[0];
rz(-0.06572289) q[0];
sx q[0];
rz(0.42732987) q[0];
rz(1.7006251) q[1];
sx q[1];
rz(-1.4375552) q[1];
sx q[1];
rz(0.89868054) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220427) q[0];
sx q[0];
rz(-1.2614095) q[0];
sx q[0];
rz(-1.1197107) q[0];
x q[1];
rz(-0.9814453) q[2];
sx q[2];
rz(-1.9716096) q[2];
sx q[2];
rz(-0.0055731853) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2262804) q[1];
sx q[1];
rz(-1.8050005) q[1];
sx q[1];
rz(2.9620902) q[1];
rz(0.54014787) q[3];
sx q[3];
rz(-1.1662048) q[3];
sx q[3];
rz(0.57309421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23540363) q[2];
sx q[2];
rz(-2.0960505) q[2];
sx q[2];
rz(1.40353) q[2];
rz(-0.67052001) q[3];
sx q[3];
rz(-1.5961921) q[3];
sx q[3];
rz(-1.3129354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32187605) q[0];
sx q[0];
rz(-0.96915594) q[0];
sx q[0];
rz(0.72625351) q[0];
rz(0.67604524) q[1];
sx q[1];
rz(-1.7773726) q[1];
sx q[1];
rz(2.3647251) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99771927) q[0];
sx q[0];
rz(-1.5251499) q[0];
sx q[0];
rz(2.2111228) q[0];
rz(-1.5258342) q[2];
sx q[2];
rz(-1.7167712) q[2];
sx q[2];
rz(3.1342497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2404691) q[1];
sx q[1];
rz(-1.8778342) q[1];
sx q[1];
rz(1.9528051) q[1];
x q[2];
rz(-1.0105206) q[3];
sx q[3];
rz(-2.0891857) q[3];
sx q[3];
rz(-0.7553525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0016510222) q[2];
sx q[2];
rz(-0.17639128) q[2];
sx q[2];
rz(1.4886935) q[2];
rz(-2.0148924) q[3];
sx q[3];
rz(-1.9727581) q[3];
sx q[3];
rz(-2.6039629) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61459944) q[0];
sx q[0];
rz(-2.499883) q[0];
sx q[0];
rz(-1.5747621) q[0];
rz(-0.78631403) q[1];
sx q[1];
rz(-0.17335261) q[1];
sx q[1];
rz(-0.39549624) q[1];
rz(-0.6720856) q[2];
sx q[2];
rz(-1.0955878) q[2];
sx q[2];
rz(0.53447117) q[2];
rz(-2.366131) q[3];
sx q[3];
rz(-1.3962593) q[3];
sx q[3];
rz(3.1229231) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
