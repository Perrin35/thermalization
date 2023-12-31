OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(1.8068846) q[0];
rz(2.788738) q[1];
sx q[1];
rz(3.3021441) q[1];
sx q[1];
rz(8.4488206) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0358409) q[0];
sx q[0];
rz(-2.5076206) q[0];
sx q[0];
rz(0.41497725) q[0];
x q[1];
rz(-0.33306723) q[2];
sx q[2];
rz(-2.1085848) q[2];
sx q[2];
rz(-1.2531467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0426892) q[1];
sx q[1];
rz(-0.96682036) q[1];
sx q[1];
rz(-0.16278111) q[1];
rz(-pi) q[2];
rz(0.39335143) q[3];
sx q[3];
rz(-0.8439807) q[3];
sx q[3];
rz(0.8918744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(-2.3577918) q[2];
rz(-0.33256724) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(-1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587104) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(2.9630307) q[0];
rz(1.3372955) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(0.006342412) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3865249) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(-0.13271876) q[0];
rz(-2.2750521) q[2];
sx q[2];
rz(-0.73359493) q[2];
sx q[2];
rz(-1.7689592) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9176863) q[1];
sx q[1];
rz(-0.31653857) q[1];
sx q[1];
rz(-2.9427337) q[1];
rz(0.56224058) q[3];
sx q[3];
rz(-1.7842818) q[3];
sx q[3];
rz(-1.436304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1938842) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(0.87810278) q[2];
rz(-0.86205035) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5851615) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(-3.1306144) q[0];
rz(-2.7745461) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(0.095741622) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.30004) q[0];
sx q[0];
rz(-1.7471572) q[0];
sx q[0];
rz(-0.15533133) q[0];
rz(-pi) q[1];
rz(-0.7272561) q[2];
sx q[2];
rz(-0.63823344) q[2];
sx q[2];
rz(-2.9122695) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0248191) q[1];
sx q[1];
rz(-1.0305335) q[1];
sx q[1];
rz(2.5953672) q[1];
rz(-1.6060353) q[3];
sx q[3];
rz(-2.1217151) q[3];
sx q[3];
rz(-2.2743724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0694971) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(-0.83479184) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(-2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(-2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(-1.0338763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82499408) q[0];
sx q[0];
rz(-0.98485095) q[0];
sx q[0];
rz(0.98866776) q[0];
rz(-pi) q[1];
rz(-0.99673523) q[2];
sx q[2];
rz(-1.1166995) q[2];
sx q[2];
rz(3.0388289) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66400601) q[1];
sx q[1];
rz(-2.9395182) q[1];
sx q[1];
rz(-1.2408153) q[1];
rz(-1.8978118) q[3];
sx q[3];
rz(-1.3033086) q[3];
sx q[3];
rz(0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(1.7948077) q[2];
rz(-2.7205617) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(-0.69798654) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088257) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-0.41473266) q[0];
rz(1.746009) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(0.57410747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87724553) q[0];
sx q[0];
rz(-0.97752042) q[0];
sx q[0];
rz(-0.1546774) q[0];
rz(2.4262869) q[2];
sx q[2];
rz(-1.5878521) q[2];
sx q[2];
rz(0.41358435) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.638006) q[1];
sx q[1];
rz(-1.2212911) q[1];
sx q[1];
rz(-1.867678) q[1];
rz(-pi) q[2];
rz(-0.19458171) q[3];
sx q[3];
rz(-0.91472018) q[3];
sx q[3];
rz(1.7723099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.627702) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742652) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(0.69818991) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(0.73227698) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9376611) q[0];
sx q[0];
rz(-1.4121778) q[0];
sx q[0];
rz(0.21477867) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76262577) q[2];
sx q[2];
rz(-1.0566933) q[2];
sx q[2];
rz(0.56157535) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9975035) q[1];
sx q[1];
rz(-0.92458506) q[1];
sx q[1];
rz(1.8540107) q[1];
rz(0.332571) q[3];
sx q[3];
rz(-0.59623527) q[3];
sx q[3];
rz(-2.3308844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(-0.30900624) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(-1.1674315) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56918615) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(2.9851595) q[0];
rz(-0.45174831) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-2.3715473) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2419287) q[0];
sx q[0];
rz(-0.97424346) q[0];
sx q[0];
rz(1.3616189) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1427878) q[2];
sx q[2];
rz(-0.6711798) q[2];
sx q[2];
rz(-2.3064248) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3502096) q[1];
sx q[1];
rz(-2.0767127) q[1];
sx q[1];
rz(-2.4398068) q[1];
x q[2];
rz(-0.50290147) q[3];
sx q[3];
rz(-1.8637878) q[3];
sx q[3];
rz(2.6754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(2.1264123) q[2];
rz(-1.7049568) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62676936) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(0.24169895) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(2.7752005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1358571) q[0];
sx q[0];
rz(-0.85548399) q[0];
sx q[0];
rz(2.4321796) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.512152) q[2];
sx q[2];
rz(-1.9295613) q[2];
sx q[2];
rz(-1.4251054) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7494292) q[1];
sx q[1];
rz(-0.87595075) q[1];
sx q[1];
rz(-1.7872582) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0189692) q[3];
sx q[3];
rz(-1.5295715) q[3];
sx q[3];
rz(-3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1239132) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-0.29385847) q[2];
rz(-3.1270694) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(-2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.83207399) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(0.88395399) q[0];
rz(2.4813095) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(0.79137897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26830772) q[0];
sx q[0];
rz(-0.76457667) q[0];
sx q[0];
rz(-1.1029878) q[0];
rz(-pi) q[1];
rz(2.0558526) q[2];
sx q[2];
rz(-1.9888478) q[2];
sx q[2];
rz(0.56318356) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89625724) q[1];
sx q[1];
rz(-2.5108813) q[1];
sx q[1];
rz(2.3920822) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16718849) q[3];
sx q[3];
rz(-1.0382004) q[3];
sx q[3];
rz(-2.5933468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(1.0212612) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(-2.5436201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026697712) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-0.60780203) q[0];
rz(2.9027477) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(-1.7609319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62771195) q[0];
sx q[0];
rz(-2.4988334) q[0];
sx q[0];
rz(1.7895262) q[0];
rz(-2.0062749) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(1.5362816) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8163562) q[1];
sx q[1];
rz(-1.163985) q[1];
sx q[1];
rz(2.7951294) q[1];
rz(-pi) q[2];
rz(-1.4019743) q[3];
sx q[3];
rz(-1.7572548) q[3];
sx q[3];
rz(2.8846915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(2.805368) q[2];
rz(2.0119038) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1375785) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(0.54429383) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(1.9359679) q[2];
sx q[2];
rz(-2.361459) q[2];
sx q[2];
rz(-0.38103719) q[2];
rz(1.929677) q[3];
sx q[3];
rz(-2.5419895) q[3];
sx q[3];
rz(0.06751577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
