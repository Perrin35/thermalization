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
rz(0.76437104) q[0];
sx q[0];
rz(1.8005014) q[0];
sx q[0];
rz(10.330248) q[0];
rz(1.5068997) q[1];
sx q[1];
rz(-2.1807179) q[1];
sx q[1];
rz(-1.6406055) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5643217) q[0];
sx q[0];
rz(-2.4080896) q[0];
sx q[0];
rz(1.2726239) q[0];
x q[1];
rz(-0.17619422) q[2];
sx q[2];
rz(-1.1374047) q[2];
sx q[2];
rz(2.7769763) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6480744) q[1];
sx q[1];
rz(-1.5941125) q[1];
sx q[1];
rz(3.0991395) q[1];
rz(-pi) q[2];
rz(1.0505643) q[3];
sx q[3];
rz(-2.5135698) q[3];
sx q[3];
rz(1.3768547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91244873) q[2];
sx q[2];
rz(-1.7488166) q[2];
sx q[2];
rz(0.32192117) q[2];
rz(2.1866482) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(-1.1436852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2442653) q[0];
sx q[0];
rz(-1.0538415) q[0];
sx q[0];
rz(-0.51189297) q[0];
rz(-1.5882675) q[1];
sx q[1];
rz(-2.6542108) q[1];
sx q[1];
rz(1.1221251) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47488362) q[0];
sx q[0];
rz(-2.8132952) q[0];
sx q[0];
rz(-0.31995456) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.000359) q[2];
sx q[2];
rz(-1.6074174) q[2];
sx q[2];
rz(0.0096732339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1074688) q[1];
sx q[1];
rz(-1.9975047) q[1];
sx q[1];
rz(0.36860768) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0997222) q[3];
sx q[3];
rz(-2.2212914) q[3];
sx q[3];
rz(-2.4471403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3731709) q[2];
sx q[2];
rz(-1.827652) q[2];
sx q[2];
rz(0.46305099) q[2];
rz(0.57524663) q[3];
sx q[3];
rz(-1.3653711) q[3];
sx q[3];
rz(-0.099253207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4082044) q[0];
sx q[0];
rz(-2.2938804) q[0];
sx q[0];
rz(0.43011618) q[0];
rz(0.46562132) q[1];
sx q[1];
rz(-0.72048134) q[1];
sx q[1];
rz(-2.5045085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89403215) q[0];
sx q[0];
rz(-1.5567509) q[0];
sx q[0];
rz(-1.5298903) q[0];
rz(-pi) q[1];
rz(1.5970206) q[2];
sx q[2];
rz(-2.3560212) q[2];
sx q[2];
rz(0.77411133) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8608537) q[1];
sx q[1];
rz(-0.85608427) q[1];
sx q[1];
rz(1.9016985) q[1];
x q[2];
rz(-2.5175321) q[3];
sx q[3];
rz(-2.5284323) q[3];
sx q[3];
rz(-1.1414736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35885262) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(0.72506881) q[2];
rz(-1.3310165) q[3];
sx q[3];
rz(-2.1264117) q[3];
sx q[3];
rz(-0.63841188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26938874) q[0];
sx q[0];
rz(-1.6870455) q[0];
sx q[0];
rz(-1.4440906) q[0];
rz(-3.0929502) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(-2.8526502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.671593) q[0];
sx q[0];
rz(-2.6942188) q[0];
sx q[0];
rz(0.98487052) q[0];
x q[1];
rz(-1.7554531) q[2];
sx q[2];
rz(-2.9884905) q[2];
sx q[2];
rz(0.36400041) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1565981) q[1];
sx q[1];
rz(-2.9356476) q[1];
sx q[1];
rz(1.0099645) q[1];
rz(-pi) q[2];
x q[2];
rz(0.069372481) q[3];
sx q[3];
rz(-0.66646229) q[3];
sx q[3];
rz(-1.762076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85663969) q[2];
sx q[2];
rz(-0.53264561) q[2];
sx q[2];
rz(0.3802158) q[2];
rz(2.5034261) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(1.1285454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.3087092) q[0];
sx q[0];
rz(-0.55279151) q[0];
sx q[0];
rz(1.0943476) q[0];
rz(2.265918) q[1];
sx q[1];
rz(-0.62962571) q[1];
sx q[1];
rz(1.099115) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35661429) q[0];
sx q[0];
rz(-1.4334502) q[0];
sx q[0];
rz(0.38909586) q[0];
x q[1];
rz(0.34277041) q[2];
sx q[2];
rz(-1.2452599) q[2];
sx q[2];
rz(0.24791644) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0620516) q[1];
sx q[1];
rz(-1.4126443) q[1];
sx q[1];
rz(-2.1248716) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3640251) q[3];
sx q[3];
rz(-0.66477699) q[3];
sx q[3];
rz(-0.95486508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2431474) q[2];
sx q[2];
rz(-1.3546319) q[2];
sx q[2];
rz(0.97664991) q[2];
rz(-2.6099033) q[3];
sx q[3];
rz(-2.6952126) q[3];
sx q[3];
rz(0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831182) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(-1.2716768) q[0];
rz(-0.42287982) q[1];
sx q[1];
rz(-1.5300749) q[1];
sx q[1];
rz(-1.5333102) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5137702) q[0];
sx q[0];
rz(-1.5585527) q[0];
sx q[0];
rz(-1.5041385) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0094290084) q[2];
sx q[2];
rz(-1.6688188) q[2];
sx q[2];
rz(2.470051) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9257838) q[1];
sx q[1];
rz(-1.5152351) q[1];
sx q[1];
rz(-2.8381398) q[1];
rz(-pi) q[2];
x q[2];
rz(1.815237) q[3];
sx q[3];
rz(-0.77271739) q[3];
sx q[3];
rz(-0.48549197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5292458) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(0.38522729) q[2];
rz(3.0569844) q[3];
sx q[3];
rz(-0.43360964) q[3];
sx q[3];
rz(-0.87578526) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92152921) q[0];
sx q[0];
rz(-1.055287) q[0];
sx q[0];
rz(0.39749843) q[0];
rz(-2.6037604) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(0.0040815512) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794681) q[0];
sx q[0];
rz(-1.4433658) q[0];
sx q[0];
rz(-2.3412933) q[0];
rz(2.053431) q[2];
sx q[2];
rz(-2.3950999) q[2];
sx q[2];
rz(1.5160402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.34985182) q[1];
sx q[1];
rz(-2.428973) q[1];
sx q[1];
rz(-0.60421677) q[1];
rz(-0.8936196) q[3];
sx q[3];
rz(-1.0056408) q[3];
sx q[3];
rz(2.2806185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6323382) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(2.9285367) q[2];
rz(-0.80900711) q[3];
sx q[3];
rz(-0.041497858) q[3];
sx q[3];
rz(1.8607148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1236561) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(-1.9785471) q[0];
rz(-0.2991547) q[1];
sx q[1];
rz(-1.9367633) q[1];
sx q[1];
rz(-0.97602731) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7892864) q[0];
sx q[0];
rz(-1.3255766) q[0];
sx q[0];
rz(-1.6847436) q[0];
x q[1];
rz(-1.7485745) q[2];
sx q[2];
rz(-1.4206829) q[2];
sx q[2];
rz(-2.3689388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73607549) q[1];
sx q[1];
rz(-1.8092522) q[1];
sx q[1];
rz(-1.2082491) q[1];
rz(-pi) q[2];
rz(-0.42404948) q[3];
sx q[3];
rz(-1.3433787) q[3];
sx q[3];
rz(0.018674803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30190793) q[2];
sx q[2];
rz(-2.8161616) q[2];
sx q[2];
rz(0.1304661) q[2];
rz(1.8179551) q[3];
sx q[3];
rz(-1.2490844) q[3];
sx q[3];
rz(-0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2015304) q[0];
sx q[0];
rz(-2.764743) q[0];
sx q[0];
rz(-1.4991722) q[0];
rz(2.6014853) q[1];
sx q[1];
rz(-2.3441548) q[1];
sx q[1];
rz(-1.3444208) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7758862) q[0];
sx q[0];
rz(-1.2822064) q[0];
sx q[0];
rz(-1.9802753) q[0];
rz(-0.70602472) q[2];
sx q[2];
rz(-0.75504061) q[2];
sx q[2];
rz(-1.4817099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0186179) q[1];
sx q[1];
rz(-0.6280762) q[1];
sx q[1];
rz(1.1372034) q[1];
rz(0.88560652) q[3];
sx q[3];
rz(-0.84784283) q[3];
sx q[3];
rz(-0.9566488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.53238791) q[2];
sx q[2];
rz(-0.33668533) q[2];
sx q[2];
rz(-1.4538291) q[2];
rz(0.42803556) q[3];
sx q[3];
rz(-1.5513709) q[3];
sx q[3];
rz(1.9868896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353272) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(-1.6499299) q[0];
rz(-0.55117575) q[1];
sx q[1];
rz(-1.0124413) q[1];
sx q[1];
rz(-2.5490882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1051286) q[0];
sx q[0];
rz(-1.8499287) q[0];
sx q[0];
rz(1.3238086) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76982381) q[2];
sx q[2];
rz(-2.4425659) q[2];
sx q[2];
rz(-1.5567284) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3687835) q[1];
sx q[1];
rz(-1.9516489) q[1];
sx q[1];
rz(-1.8541149) q[1];
rz(-pi) q[2];
x q[2];
rz(0.09162147) q[3];
sx q[3];
rz(-2.6406248) q[3];
sx q[3];
rz(-0.17551455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42980117) q[2];
sx q[2];
rz(-1.0415123) q[2];
sx q[2];
rz(3.0832624) q[2];
rz(-0.37060261) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(0.0037732865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2881099) q[0];
sx q[0];
rz(-1.3567038) q[0];
sx q[0];
rz(0.10467341) q[0];
rz(1.3660322) q[1];
sx q[1];
rz(-2.3401101) q[1];
sx q[1];
rz(0.95536864) q[1];
rz(0.38717196) q[2];
sx q[2];
rz(-2.5324814) q[2];
sx q[2];
rz(-1.0003288) q[2];
rz(-0.52363734) q[3];
sx q[3];
rz(-2.0515531) q[3];
sx q[3];
rz(2.6691347) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
