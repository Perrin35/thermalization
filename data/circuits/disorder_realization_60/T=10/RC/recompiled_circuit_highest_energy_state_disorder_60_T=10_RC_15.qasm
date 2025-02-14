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
rz(4.6484923) q[1];
sx q[1];
rz(5.3223106) q[1];
sx q[1];
rz(7.9237908) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.577271) q[0];
sx q[0];
rz(-2.4080896) q[0];
sx q[0];
rz(-1.2726239) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2086966) q[2];
sx q[2];
rz(-2.675867) q[2];
sx q[2];
rz(-0.036368792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6480744) q[1];
sx q[1];
rz(-1.5474802) q[1];
sx q[1];
rz(-3.0991395) q[1];
rz(-pi) q[2];
rz(0.34637536) q[3];
sx q[3];
rz(-1.0358255) q[3];
sx q[3];
rz(-0.76081027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2291439) q[2];
sx q[2];
rz(-1.7488166) q[2];
sx q[2];
rz(-2.8196715) q[2];
rz(0.95494444) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(-1.9979075) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2442653) q[0];
sx q[0];
rz(-1.0538415) q[0];
sx q[0];
rz(2.6296997) q[0];
rz(-1.5533252) q[1];
sx q[1];
rz(-2.6542108) q[1];
sx q[1];
rz(2.0194676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3998386) q[0];
sx q[0];
rz(-1.6723833) q[0];
sx q[0];
rz(-2.8288657) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4830515) q[2];
sx q[2];
rz(-2.7105687) q[2];
sx q[2];
rz(-1.4813678) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8367887) q[1];
sx q[1];
rz(-1.9049909) q[1];
sx q[1];
rz(2.0242974) q[1];
rz(-pi) q[2];
rz(2.0418704) q[3];
sx q[3];
rz(-2.2212914) q[3];
sx q[3];
rz(0.6944523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3731709) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(-2.6785417) q[2];
rz(-2.566346) q[3];
sx q[3];
rz(-1.3653711) q[3];
sx q[3];
rz(3.0423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4082044) q[0];
sx q[0];
rz(-0.84771228) q[0];
sx q[0];
rz(0.43011618) q[0];
rz(-2.6759713) q[1];
sx q[1];
rz(-0.72048134) q[1];
sx q[1];
rz(0.63708416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67733902) q[0];
sx q[0];
rz(-1.6116983) q[0];
sx q[0];
rz(-3.1275355) q[0];
x q[1];
rz(3.1153684) q[2];
sx q[2];
rz(-2.3560239) q[2];
sx q[2];
rz(0.81119591) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7635258) q[1];
sx q[1];
rz(-0.77516205) q[1];
sx q[1];
rz(0.35825348) q[1];
rz(-pi) q[2];
rz(1.1807084) q[3];
sx q[3];
rz(-1.0849139) q[3];
sx q[3];
rz(0.41956832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.78274) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(-2.4165238) q[2];
rz(1.3310165) q[3];
sx q[3];
rz(-1.015181) q[3];
sx q[3];
rz(2.5031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.26938874) q[0];
sx q[0];
rz(-1.4545472) q[0];
sx q[0];
rz(-1.697502) q[0];
rz(3.0929502) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(2.8526502) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64000696) q[0];
sx q[0];
rz(-1.8123535) q[0];
sx q[0];
rz(1.1904741) q[0];
rz(1.4202576) q[2];
sx q[2];
rz(-1.5427914) q[2];
sx q[2];
rz(1.3893407) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98499456) q[1];
sx q[1];
rz(-2.9356476) q[1];
sx q[1];
rz(2.1316281) q[1];
rz(-pi) q[2];
rz(2.4763002) q[3];
sx q[3];
rz(-1.5279309) q[3];
sx q[3];
rz(-2.8957518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.85663969) q[2];
sx q[2];
rz(-2.608947) q[2];
sx q[2];
rz(0.3802158) q[2];
rz(-2.5034261) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(2.0130472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83288348) q[0];
sx q[0];
rz(-2.5888011) q[0];
sx q[0];
rz(1.0943476) q[0];
rz(0.87567466) q[1];
sx q[1];
rz(-2.5119669) q[1];
sx q[1];
rz(-2.0424776) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1581061) q[0];
sx q[0];
rz(-1.1855584) q[0];
sx q[0];
rz(-1.7190821) q[0];
rz(-pi) q[1];
rz(-2.7988222) q[2];
sx q[2];
rz(-1.8963328) q[2];
sx q[2];
rz(2.8936762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5531986) q[1];
sx q[1];
rz(-1.0244245) q[1];
sx q[1];
rz(0.1853892) q[1];
x q[2];
rz(0.15954475) q[3];
sx q[3];
rz(-0.92260781) q[3];
sx q[3];
rz(1.2153347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2431474) q[2];
sx q[2];
rz(-1.3546319) q[2];
sx q[2];
rz(0.97664991) q[2];
rz(2.6099033) q[3];
sx q[3];
rz(-0.44638005) q[3];
sx q[3];
rz(-2.4272052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0584745) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(1.8699159) q[0];
rz(-0.42287982) q[1];
sx q[1];
rz(-1.6115178) q[1];
sx q[1];
rz(-1.6082825) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9031808) q[0];
sx q[0];
rz(-3.0738214) q[0];
sx q[0];
rz(-1.3890024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6663867) q[2];
sx q[2];
rz(-0.098473452) q[2];
sx q[2];
rz(-0.57548875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3376006) q[1];
sx q[1];
rz(-1.873766) q[1];
sx q[1];
rz(-1.6290118) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9098784) q[3];
sx q[3];
rz(-2.3149256) q[3];
sx q[3];
rz(0.8207013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5292458) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(-0.38522729) q[2];
rz(-3.0569844) q[3];
sx q[3];
rz(-2.707983) q[3];
sx q[3];
rz(-0.87578526) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2200634) q[0];
sx q[0];
rz(-2.0863057) q[0];
sx q[0];
rz(2.7440942) q[0];
rz(0.53783224) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(-3.1375111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57848141) q[0];
sx q[0];
rz(-0.77881506) q[0];
sx q[0];
rz(1.7527197) q[0];
rz(-0.40553826) q[2];
sx q[2];
rz(-0.92541646) q[2];
sx q[2];
rz(-2.1359512) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4392802) q[1];
sx q[1];
rz(-1.951362) q[1];
sx q[1];
rz(-2.5234532) q[1];
rz(0.68304045) q[3];
sx q[3];
rz(-2.1284102) q[3];
sx q[3];
rz(-2.838359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5092545) q[2];
sx q[2];
rz(-1.961144) q[2];
sx q[2];
rz(2.9285367) q[2];
rz(-0.80900711) q[3];
sx q[3];
rz(-3.1000948) q[3];
sx q[3];
rz(1.2808778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1236561) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(-1.1630455) q[0];
rz(-2.842438) q[1];
sx q[1];
rz(-1.9367633) q[1];
sx q[1];
rz(0.97602731) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3487745) q[0];
sx q[0];
rz(-2.8716757) q[0];
sx q[0];
rz(2.7151373) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86324228) q[2];
sx q[2];
rz(-0.23216557) q[2];
sx q[2];
rz(-1.6492998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27777003) q[1];
sx q[1];
rz(-0.43102095) q[1];
sx q[1];
rz(0.96993877) q[1];
rz(-pi) q[2];
rz(-1.3221413) q[3];
sx q[3];
rz(-1.9832522) q[3];
sx q[3];
rz(1.4506884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8396847) q[2];
sx q[2];
rz(-2.8161616) q[2];
sx q[2];
rz(-3.0111266) q[2];
rz(-1.8179551) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(-0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94006222) q[0];
sx q[0];
rz(-2.764743) q[0];
sx q[0];
rz(-1.6424204) q[0];
rz(2.6014853) q[1];
sx q[1];
rz(-2.3441548) q[1];
sx q[1];
rz(-1.3444208) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3657065) q[0];
sx q[0];
rz(-1.2822064) q[0];
sx q[0];
rz(1.9802753) q[0];
x q[1];
rz(-2.5201468) q[2];
sx q[2];
rz(-2.0315731) q[2];
sx q[2];
rz(-2.4968392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.498913) q[1];
sx q[1];
rz(-1.0084001) q[1];
sx q[1];
rz(-0.29611961) q[1];
x q[2];
rz(0.85050438) q[3];
sx q[3];
rz(-2.0652186) q[3];
sx q[3];
rz(-1.1098343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6092047) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(-1.6877635) q[2];
rz(2.7135571) q[3];
sx q[3];
rz(-1.5513709) q[3];
sx q[3];
rz(1.154703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353272) q[0];
sx q[0];
rz(-2.689671) q[0];
sx q[0];
rz(-1.4916627) q[0];
rz(2.5904169) q[1];
sx q[1];
rz(-1.0124413) q[1];
sx q[1];
rz(-2.5490882) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77757072) q[0];
sx q[0];
rz(-2.7710272) q[0];
sx q[0];
rz(2.4353566) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3717688) q[2];
sx q[2];
rz(-0.69902674) q[2];
sx q[2];
rz(-1.5567284) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7043276) q[1];
sx q[1];
rz(-2.6710837) q[1];
sx q[1];
rz(0.60948845) q[1];
x q[2];
rz(1.6208526) q[3];
sx q[3];
rz(-2.069469) q[3];
sx q[3];
rz(0.071144516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42980117) q[2];
sx q[2];
rz(-1.0415123) q[2];
sx q[2];
rz(0.05833021) q[2];
rz(0.37060261) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(-0.0037732865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.8534828) q[0];
sx q[0];
rz(-1.7848889) q[0];
sx q[0];
rz(-3.0369192) q[0];
rz(-1.3660322) q[1];
sx q[1];
rz(-0.80148253) q[1];
sx q[1];
rz(-2.186224) q[1];
rz(1.8283394) q[2];
sx q[2];
rz(-1.0124442) q[2];
sx q[2];
rz(-0.53895216) q[2];
rz(1.0287063) q[3];
sx q[3];
rz(-1.1114612) q[3];
sx q[3];
rz(-2.3041861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
