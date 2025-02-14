OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3684664) q[0];
sx q[0];
rz(-2.3946895) q[0];
sx q[0];
rz(-0.84063831) q[0];
rz(-3.0199938) q[1];
sx q[1];
rz(-1.8687948) q[1];
sx q[1];
rz(2.8425541) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5344319) q[0];
sx q[0];
rz(-2.2272155) q[0];
sx q[0];
rz(0.54404152) q[0];
rz(-pi) q[1];
rz(1.4242572) q[2];
sx q[2];
rz(-2.6508109) q[2];
sx q[2];
rz(-2.7637568) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12451367) q[1];
sx q[1];
rz(-0.28055596) q[1];
sx q[1];
rz(2.6002167) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2848008) q[3];
sx q[3];
rz(-2.6702849) q[3];
sx q[3];
rz(-0.23391868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3806939) q[2];
sx q[2];
rz(-2.0913405) q[2];
sx q[2];
rz(2.8139581) q[2];
rz(1.7662175) q[3];
sx q[3];
rz(-1.7211434) q[3];
sx q[3];
rz(2.6473141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37680092) q[0];
sx q[0];
rz(-1.6015653) q[0];
sx q[0];
rz(2.2862527) q[0];
rz(3.1335462) q[1];
sx q[1];
rz(-1.2866373) q[1];
sx q[1];
rz(2.6904552) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6159003) q[0];
sx q[0];
rz(-0.85291686) q[0];
sx q[0];
rz(-0.20882102) q[0];
rz(-pi) q[1];
rz(-1.3786267) q[2];
sx q[2];
rz(-2.333667) q[2];
sx q[2];
rz(2.5469123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0135368) q[1];
sx q[1];
rz(-1.3971796) q[1];
sx q[1];
rz(-2.6980969) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69591095) q[3];
sx q[3];
rz(-2.5491121) q[3];
sx q[3];
rz(2.185948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.96192876) q[2];
sx q[2];
rz(-2.0018061) q[2];
sx q[2];
rz(3.0120604) q[2];
rz(-0.1772964) q[3];
sx q[3];
rz(-0.57759053) q[3];
sx q[3];
rz(3.0887443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.391908) q[0];
sx q[0];
rz(-0.41202298) q[0];
sx q[0];
rz(-2.5823197) q[0];
rz(-0.094712146) q[1];
sx q[1];
rz(-1.5385224) q[1];
sx q[1];
rz(-0.47725484) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4995296) q[0];
sx q[0];
rz(-1.9669878) q[0];
sx q[0];
rz(2.8642333) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37896403) q[2];
sx q[2];
rz(-0.38658374) q[2];
sx q[2];
rz(1.6340337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2169358) q[1];
sx q[1];
rz(-0.81243304) q[1];
sx q[1];
rz(-1.2595909) q[1];
rz(-0.84945143) q[3];
sx q[3];
rz(-1.1616544) q[3];
sx q[3];
rz(0.22709286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39530784) q[2];
sx q[2];
rz(-1.2664653) q[2];
sx q[2];
rz(-2.2115808) q[2];
rz(-1.8837455) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(3.0959082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78087085) q[0];
sx q[0];
rz(-1.5732795) q[0];
sx q[0];
rz(-1.038653) q[0];
rz(-2.5301798) q[1];
sx q[1];
rz(-0.88714209) q[1];
sx q[1];
rz(-0.92322737) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34133615) q[0];
sx q[0];
rz(-1.0025327) q[0];
sx q[0];
rz(-0.41109127) q[0];
rz(-3.1298248) q[2];
sx q[2];
rz(-1.9885049) q[2];
sx q[2];
rz(-2.2209446) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.23455305) q[1];
sx q[1];
rz(-2.4133293) q[1];
sx q[1];
rz(1.6978463) q[1];
rz(2.8037854) q[3];
sx q[3];
rz(-0.67854133) q[3];
sx q[3];
rz(-1.9935009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5103147) q[2];
sx q[2];
rz(-2.6738561) q[2];
sx q[2];
rz(0.54747096) q[2];
rz(3.0771717) q[3];
sx q[3];
rz(-1.3037325) q[3];
sx q[3];
rz(-2.2678383) q[3];
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
rz(pi/2) q[3];
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
rz(2.9861458) q[0];
sx q[0];
rz(-0.90279818) q[0];
sx q[0];
rz(-1.6814394) q[0];
rz(-1.0001146) q[1];
sx q[1];
rz(-2.2747048) q[1];
sx q[1];
rz(3.0217081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016235624) q[0];
sx q[0];
rz(-0.88847697) q[0];
sx q[0];
rz(-2.8518139) q[0];
x q[1];
rz(0.85544678) q[2];
sx q[2];
rz(-0.79022932) q[2];
sx q[2];
rz(1.2145192) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8441808) q[1];
sx q[1];
rz(-1.4919123) q[1];
sx q[1];
rz(-0.21551883) q[1];
x q[2];
rz(-1.0644887) q[3];
sx q[3];
rz(-2.1112006) q[3];
sx q[3];
rz(0.04549724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.039006058) q[2];
sx q[2];
rz(-0.90709364) q[2];
sx q[2];
rz(-0.26068035) q[2];
rz(0.87812224) q[3];
sx q[3];
rz(-1.9120646) q[3];
sx q[3];
rz(-0.78316435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53073019) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(-2.1395785) q[0];
rz(2.7643381) q[1];
sx q[1];
rz(-2.1899624) q[1];
sx q[1];
rz(0.9446876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024600895) q[0];
sx q[0];
rz(-0.62378609) q[0];
sx q[0];
rz(0.37895112) q[0];
rz(-pi) q[1];
rz(2.5414921) q[2];
sx q[2];
rz(-0.99008152) q[2];
sx q[2];
rz(2.3490259) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3516253) q[1];
sx q[1];
rz(-1.2500487) q[1];
sx q[1];
rz(-1.7299537) q[1];
rz(-pi) q[2];
rz(-0.11774534) q[3];
sx q[3];
rz(-1.0317993) q[3];
sx q[3];
rz(2.4460411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79823309) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(3.0976683) q[2];
rz(-1.6262866) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(1.4656434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8632904) q[0];
sx q[0];
rz(-2.589812) q[0];
sx q[0];
rz(-0.58468753) q[0];
rz(2.0901285) q[1];
sx q[1];
rz(-2.3221071) q[1];
sx q[1];
rz(-2.6712766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6795652) q[0];
sx q[0];
rz(-0.60067486) q[0];
sx q[0];
rz(2.8810049) q[0];
rz(2.3114071) q[2];
sx q[2];
rz(-2.6174394) q[2];
sx q[2];
rz(-2.1955127) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9661588) q[1];
sx q[1];
rz(-2.325238) q[1];
sx q[1];
rz(-2.6360126) q[1];
rz(3.028585) q[3];
sx q[3];
rz(-1.5586434) q[3];
sx q[3];
rz(-0.069610217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.26479244) q[2];
sx q[2];
rz(-2.6046643) q[2];
sx q[2];
rz(2.8301767) q[2];
rz(0.16820678) q[3];
sx q[3];
rz(-1.5663389) q[3];
sx q[3];
rz(0.28347191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6300221) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(-2.2054963) q[0];
rz(-0.49631897) q[1];
sx q[1];
rz(-0.66718188) q[1];
sx q[1];
rz(2.671303) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054809991) q[0];
sx q[0];
rz(-0.58712372) q[0];
sx q[0];
rz(2.5820288) q[0];
x q[1];
rz(-1.5580721) q[2];
sx q[2];
rz(-0.1977405) q[2];
sx q[2];
rz(-1.5997353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6185604) q[1];
sx q[1];
rz(-1.2292687) q[1];
sx q[1];
rz(-3.0681472) q[1];
rz(-2.5816282) q[3];
sx q[3];
rz(-0.32295152) q[3];
sx q[3];
rz(-1.1895869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63142598) q[2];
sx q[2];
rz(-1.7743856) q[2];
sx q[2];
rz(-1.4754971) q[2];
rz(0.79536074) q[3];
sx q[3];
rz(-2.9795591) q[3];
sx q[3];
rz(1.2696666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0805761) q[0];
sx q[0];
rz(-2.5503655) q[0];
sx q[0];
rz(0.06037816) q[0];
rz(-2.9810442) q[1];
sx q[1];
rz(-1.6146654) q[1];
sx q[1];
rz(2.1626332) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85764393) q[0];
sx q[0];
rz(-1.4547087) q[0];
sx q[0];
rz(2.5124971) q[0];
rz(-pi) q[1];
rz(1.9039745) q[2];
sx q[2];
rz(-0.71014437) q[2];
sx q[2];
rz(1.6802481) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.95293249) q[1];
sx q[1];
rz(-0.97129909) q[1];
sx q[1];
rz(-0.99296928) q[1];
rz(-pi) q[2];
rz(-1.6230725) q[3];
sx q[3];
rz(-2.5617122) q[3];
sx q[3];
rz(1.9573905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0742566) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(3.0687029) q[2];
rz(2.5439751) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(0.0028751956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.446796) q[0];
sx q[0];
rz(-2.1294761) q[0];
sx q[0];
rz(0.52892518) q[0];
rz(2.9534598) q[1];
sx q[1];
rz(-2.4362322) q[1];
sx q[1];
rz(0.58427748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0463294) q[0];
sx q[0];
rz(-1.3130616) q[0];
sx q[0];
rz(-1.3168174) q[0];
x q[1];
rz(-2.7033349) q[2];
sx q[2];
rz(-1.2152248) q[2];
sx q[2];
rz(-2.8388765) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1288695) q[1];
sx q[1];
rz(-1.1616316) q[1];
sx q[1];
rz(1.6175458) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8680471) q[3];
sx q[3];
rz(-0.97018948) q[3];
sx q[3];
rz(-0.41714868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75954413) q[2];
sx q[2];
rz(-2.1376762) q[2];
sx q[2];
rz(3.0697401) q[2];
rz(-1.0233277) q[3];
sx q[3];
rz(-1.5691248) q[3];
sx q[3];
rz(0.48880997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95579424) q[0];
sx q[0];
rz(-2.4837942) q[0];
sx q[0];
rz(2.876045) q[0];
rz(0.67509782) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(0.51495348) q[2];
sx q[2];
rz(-1.5025768) q[2];
sx q[2];
rz(1.2160355) q[2];
rz(-0.16719462) q[3];
sx q[3];
rz(-1.4124558) q[3];
sx q[3];
rz(2.4168933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
