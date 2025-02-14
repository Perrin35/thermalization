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
rz(-2.5833997) q[0];
sx q[0];
rz(-0.71891958) q[0];
sx q[0];
rz(2.4655226) q[0];
rz(0.27322912) q[1];
sx q[1];
rz(-2.1844808) q[1];
sx q[1];
rz(0.91125429) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34150037) q[0];
sx q[0];
rz(-1.5507292) q[0];
sx q[0];
rz(1.5930727) q[0];
rz(-pi) q[1];
rz(-0.010013117) q[2];
sx q[2];
rz(-0.62558936) q[2];
sx q[2];
rz(-2.6715476) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75685749) q[1];
sx q[1];
rz(-0.57418203) q[1];
sx q[1];
rz(-3.0822445) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1333255) q[3];
sx q[3];
rz(-0.76078868) q[3];
sx q[3];
rz(-0.68851346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91997826) q[2];
sx q[2];
rz(-2.6875434) q[2];
sx q[2];
rz(2.3294219) q[2];
rz(0.075856097) q[3];
sx q[3];
rz(-0.64801884) q[3];
sx q[3];
rz(-0.59504741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(2.6512063) q[0];
sx q[0];
rz(-2.6631329) q[0];
sx q[0];
rz(-1.2745717) q[0];
rz(-1.6365341) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(-1.9777745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5620016) q[0];
sx q[0];
rz(-1.1222071) q[0];
sx q[0];
rz(-0.513784) q[0];
rz(-0.69804783) q[2];
sx q[2];
rz(-0.62605941) q[2];
sx q[2];
rz(-1.2241838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5806942) q[1];
sx q[1];
rz(-0.95798641) q[1];
sx q[1];
rz(1.8263019) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9979565) q[3];
sx q[3];
rz(-1.2795951) q[3];
sx q[3];
rz(-1.8607651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85094467) q[2];
sx q[2];
rz(-3.1182351) q[2];
sx q[2];
rz(1.5811496) q[2];
rz(-0.23049878) q[3];
sx q[3];
rz(-2.0932308) q[3];
sx q[3];
rz(-0.44052625) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5480963) q[0];
sx q[0];
rz(-2.8289712) q[0];
sx q[0];
rz(-0.12259677) q[0];
rz(-0.7971881) q[1];
sx q[1];
rz(-1.0120579) q[1];
sx q[1];
rz(-3.0534993) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3159972) q[0];
sx q[0];
rz(-1.8261709) q[0];
sx q[0];
rz(-2.3488464) q[0];
rz(1.2253907) q[2];
sx q[2];
rz(-0.24125762) q[2];
sx q[2];
rz(0.119482) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8338203) q[1];
sx q[1];
rz(-1.4188179) q[1];
sx q[1];
rz(-3.003158) q[1];
rz(-pi) q[2];
rz(0.53998472) q[3];
sx q[3];
rz(-1.8740046) q[3];
sx q[3];
rz(-1.4410415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88283551) q[2];
sx q[2];
rz(-1.5828524) q[2];
sx q[2];
rz(1.8176414) q[2];
rz(-0.49805182) q[3];
sx q[3];
rz(-2.1784454) q[3];
sx q[3];
rz(-0.74670416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28904706) q[0];
sx q[0];
rz(-2.9883224) q[0];
sx q[0];
rz(2.9943941) q[0];
rz(-2.4920801) q[1];
sx q[1];
rz(-1.6308866) q[1];
sx q[1];
rz(-1.6726327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2788852) q[0];
sx q[0];
rz(-1.577958) q[0];
sx q[0];
rz(1.7938406) q[0];
rz(-pi) q[1];
rz(1.1403667) q[2];
sx q[2];
rz(-1.7375611) q[2];
sx q[2];
rz(2.8362136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8007954) q[1];
sx q[1];
rz(-0.95424517) q[1];
sx q[1];
rz(0.430213) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58261915) q[3];
sx q[3];
rz(-1.1938059) q[3];
sx q[3];
rz(-2.2347799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2088251) q[2];
sx q[2];
rz(-0.494445) q[2];
sx q[2];
rz(-2.8998568) q[2];
rz(-2.406481) q[3];
sx q[3];
rz(-1.7416411) q[3];
sx q[3];
rz(0.66489768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6786137) q[0];
sx q[0];
rz(-2.0942056) q[0];
sx q[0];
rz(-3.0188634) q[0];
rz(-0.8853451) q[1];
sx q[1];
rz(-2.131999) q[1];
sx q[1];
rz(-0.56040323) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2764328) q[0];
sx q[0];
rz(-1.7483646) q[0];
sx q[0];
rz(-2.0899169) q[0];
x q[1];
rz(2.8444195) q[2];
sx q[2];
rz(-1.3577596) q[2];
sx q[2];
rz(1.4122672) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97241114) q[1];
sx q[1];
rz(-1.3126846) q[1];
sx q[1];
rz(-1.0960137) q[1];
rz(-pi) q[2];
rz(1.089483) q[3];
sx q[3];
rz(-1.169724) q[3];
sx q[3];
rz(2.8477251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.6178599) q[2];
sx q[2];
rz(-0.79404074) q[2];
sx q[2];
rz(2.9368371) q[2];
rz(-0.93920416) q[3];
sx q[3];
rz(-1.771628) q[3];
sx q[3];
rz(-0.088976629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8646249) q[0];
sx q[0];
rz(-0.11700103) q[0];
sx q[0];
rz(-0.65698874) q[0];
rz(-2.9981546) q[1];
sx q[1];
rz(-1.5564432) q[1];
sx q[1];
rz(2.4030446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6929106) q[0];
sx q[0];
rz(-1.686272) q[0];
sx q[0];
rz(-0.67489745) q[0];
rz(-pi) q[1];
rz(1.9759278) q[2];
sx q[2];
rz(-2.6141593) q[2];
sx q[2];
rz(-2.1857174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.025102928) q[1];
sx q[1];
rz(-1.2560351) q[1];
sx q[1];
rz(-0.88309755) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70921405) q[3];
sx q[3];
rz(-1.1284168) q[3];
sx q[3];
rz(-2.1565587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2584381) q[2];
sx q[2];
rz(-2.482735) q[2];
sx q[2];
rz(0.0030041791) q[2];
rz(0.7891807) q[3];
sx q[3];
rz(-2.7572258) q[3];
sx q[3];
rz(1.9743617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060519144) q[0];
sx q[0];
rz(-2.7923212) q[0];
sx q[0];
rz(-0.14807598) q[0];
rz(-0.5332467) q[1];
sx q[1];
rz(-0.24594578) q[1];
sx q[1];
rz(1.5124793) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63321823) q[0];
sx q[0];
rz(-1.1092767) q[0];
sx q[0];
rz(1.2715696) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81432366) q[2];
sx q[2];
rz(-2.0481128) q[2];
sx q[2];
rz(-2.5392591) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5606192) q[1];
sx q[1];
rz(-1.6915913) q[1];
sx q[1];
rz(-1.2289775) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9307327) q[3];
sx q[3];
rz(-2.2818344) q[3];
sx q[3];
rz(-2.8703727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.44275722) q[2];
sx q[2];
rz(-1.4407225) q[2];
sx q[2];
rz(0.093078144) q[2];
rz(0.062151521) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(1.9183581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1208948) q[0];
sx q[0];
rz(-0.59497213) q[0];
sx q[0];
rz(-2.4769532) q[0];
rz(-1.7918034) q[1];
sx q[1];
rz(-2.2638958) q[1];
sx q[1];
rz(2.452449) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6038325) q[0];
sx q[0];
rz(-1.4546772) q[0];
sx q[0];
rz(2.8962747) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6479826) q[2];
sx q[2];
rz(-2.2008332) q[2];
sx q[2];
rz(-2.4929327) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7245111) q[1];
sx q[1];
rz(-1.8838091) q[1];
sx q[1];
rz(-1.4177807) q[1];
x q[2];
rz(-0.1647968) q[3];
sx q[3];
rz(-1.6305123) q[3];
sx q[3];
rz(-2.7415558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.050921116) q[2];
sx q[2];
rz(-2.4854269) q[2];
sx q[2];
rz(-1.7655168) q[2];
rz(1.9164267) q[3];
sx q[3];
rz(-2.4296032) q[3];
sx q[3];
rz(3.0998949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10251481) q[0];
sx q[0];
rz(-3.097105) q[0];
sx q[0];
rz(-2.5503889) q[0];
rz(2.8996331) q[1];
sx q[1];
rz(-1.0958902) q[1];
sx q[1];
rz(0.42344365) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11277448) q[0];
sx q[0];
rz(-1.8736381) q[0];
sx q[0];
rz(-1.8870728) q[0];
x q[1];
rz(-2.9482163) q[2];
sx q[2];
rz(-0.49623734) q[2];
sx q[2];
rz(-1.9691182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0656506) q[1];
sx q[1];
rz(-1.0417176) q[1];
sx q[1];
rz(-0.34543646) q[1];
x q[2];
rz(0.54355346) q[3];
sx q[3];
rz(-2.4494236) q[3];
sx q[3];
rz(2.8080733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.55687904) q[2];
sx q[2];
rz(-2.5355279) q[2];
sx q[2];
rz(1.9752183) q[2];
rz(1.1276468) q[3];
sx q[3];
rz(-0.96846628) q[3];
sx q[3];
rz(-0.52496547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0135076) q[0];
sx q[0];
rz(-2.7078244) q[0];
sx q[0];
rz(2.4684913) q[0];
rz(-0.85064864) q[1];
sx q[1];
rz(-1.3530082) q[1];
sx q[1];
rz(2.8180715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6768247) q[0];
sx q[0];
rz(-2.6861405) q[0];
sx q[0];
rz(-2.1663329) q[0];
rz(-0.27574972) q[2];
sx q[2];
rz(-2.033503) q[2];
sx q[2];
rz(-1.6088865) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21491815) q[1];
sx q[1];
rz(-1.2387382) q[1];
sx q[1];
rz(-2.3658793) q[1];
rz(-pi) q[2];
rz(-0.58384044) q[3];
sx q[3];
rz(-1.5773609) q[3];
sx q[3];
rz(-3.1232426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9234151) q[2];
sx q[2];
rz(-0.18109334) q[2];
sx q[2];
rz(3.0695445) q[2];
rz(1.0142903) q[3];
sx q[3];
rz(-0.91599661) q[3];
sx q[3];
rz(-0.54916507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88076787) q[0];
sx q[0];
rz(-1.4243955) q[0];
sx q[0];
rz(-2.0212174) q[0];
rz(-0.33521677) q[1];
sx q[1];
rz(-2.4991279) q[1];
sx q[1];
rz(2.0229708) q[1];
rz(1.8374683) q[2];
sx q[2];
rz(-0.49378569) q[2];
sx q[2];
rz(2.7616382) q[2];
rz(-2.9670197) q[3];
sx q[3];
rz(-0.1452419) q[3];
sx q[3];
rz(-2.8066842) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
