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
rz(1.9873729) q[0];
sx q[0];
rz(-1.6183102) q[0];
sx q[0];
rz(-0.14259882) q[0];
rz(-2.5481186) q[1];
sx q[1];
rz(-0.65294099) q[1];
sx q[1];
rz(-2.0963734) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3680688) q[0];
sx q[0];
rz(-2.0288101) q[0];
sx q[0];
rz(-2.5450443) q[0];
x q[1];
rz(-0.036816557) q[2];
sx q[2];
rz(-0.93339257) q[2];
sx q[2];
rz(0.76257818) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7683623) q[1];
sx q[1];
rz(-1.6715896) q[1];
sx q[1];
rz(-1.412315) q[1];
rz(-3.0229085) q[3];
sx q[3];
rz(-2.136345) q[3];
sx q[3];
rz(-2.8142336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1656701) q[2];
sx q[2];
rz(-1.0079577) q[2];
sx q[2];
rz(-0.21318501) q[2];
rz(-2.2255157) q[3];
sx q[3];
rz(-0.35917425) q[3];
sx q[3];
rz(-0.38755125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.85584545) q[0];
sx q[0];
rz(-0.88749945) q[0];
sx q[0];
rz(-2.6708653) q[0];
rz(1.1031021) q[1];
sx q[1];
rz(-0.50869894) q[1];
sx q[1];
rz(0.57977605) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.267201) q[0];
sx q[0];
rz(-1.6306595) q[0];
sx q[0];
rz(-2.6169928) q[0];
rz(-pi) q[1];
rz(-2.6227399) q[2];
sx q[2];
rz(-2.5343072) q[2];
sx q[2];
rz(-2.3407206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6439542) q[1];
sx q[1];
rz(-2.4259287) q[1];
sx q[1];
rz(0.96244754) q[1];
rz(0.96478077) q[3];
sx q[3];
rz(-0.44108118) q[3];
sx q[3];
rz(-1.4853322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9819928) q[2];
sx q[2];
rz(-2.2856568) q[2];
sx q[2];
rz(-0.095005438) q[2];
rz(-2.5659918) q[3];
sx q[3];
rz(-0.78229457) q[3];
sx q[3];
rz(-1.4466205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59548241) q[0];
sx q[0];
rz(-2.4330916) q[0];
sx q[0];
rz(2.8693759) q[0];
rz(-3.0369924) q[1];
sx q[1];
rz(-2.101254) q[1];
sx q[1];
rz(-1.4629755) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42727106) q[0];
sx q[0];
rz(-2.2939689) q[0];
sx q[0];
rz(0.066825213) q[0];
x q[1];
rz(-0.06426364) q[2];
sx q[2];
rz(-2.3552364) q[2];
sx q[2];
rz(-2.3051895) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6879287) q[1];
sx q[1];
rz(-1.8519823) q[1];
sx q[1];
rz(0.67800501) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9862254) q[3];
sx q[3];
rz(-0.37305635) q[3];
sx q[3];
rz(-3.1257926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0649123) q[2];
sx q[2];
rz(-1.7408337) q[2];
sx q[2];
rz(0.40985516) q[2];
rz(3.0900433) q[3];
sx q[3];
rz(-2.1487273) q[3];
sx q[3];
rz(2.4079017) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9728397) q[0];
sx q[0];
rz(-0.77744716) q[0];
sx q[0];
rz(0.69946104) q[0];
rz(0.93060023) q[1];
sx q[1];
rz(-1.365463) q[1];
sx q[1];
rz(-1.0994937) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.612815) q[0];
sx q[0];
rz(-1.8474425) q[0];
sx q[0];
rz(2.1601281) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9540953) q[2];
sx q[2];
rz(-2.5591203) q[2];
sx q[2];
rz(-0.7815643) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9605254) q[1];
sx q[1];
rz(-1.2511504) q[1];
sx q[1];
rz(1.245371) q[1];
x q[2];
rz(0.64208018) q[3];
sx q[3];
rz(-1.2386981) q[3];
sx q[3];
rz(-1.4043216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3804271) q[2];
sx q[2];
rz(-1.0720422) q[2];
sx q[2];
rz(-0.94044828) q[2];
rz(-0.56194168) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(0.57653069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(3.1171653) q[0];
sx q[0];
rz(-0.086240135) q[0];
sx q[0];
rz(-1.0166136) q[0];
rz(-0.92574614) q[1];
sx q[1];
rz(-0.75030202) q[1];
sx q[1];
rz(3.064916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31875089) q[0];
sx q[0];
rz(-2.7116576) q[0];
sx q[0];
rz(-1.1640401) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.935604) q[2];
sx q[2];
rz(-1.811337) q[2];
sx q[2];
rz(-2.5702916) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3073342) q[1];
sx q[1];
rz(-1.3594207) q[1];
sx q[1];
rz(2.1366875) q[1];
rz(-pi) q[2];
x q[2];
rz(1.881095) q[3];
sx q[3];
rz(-2.3633133) q[3];
sx q[3];
rz(-2.8097514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3735247) q[2];
sx q[2];
rz(-0.7879476) q[2];
sx q[2];
rz(3.0184271) q[2];
rz(-2.4672616) q[3];
sx q[3];
rz(-2.6682523) q[3];
sx q[3];
rz(2.3136852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8318361) q[0];
sx q[0];
rz(-0.18853822) q[0];
sx q[0];
rz(-2.6605666) q[0];
rz(-2.9006529) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(-0.13490881) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8277934) q[0];
sx q[0];
rz(-2.302357) q[0];
sx q[0];
rz(-2.2025432) q[0];
x q[1];
rz(-0.53510869) q[2];
sx q[2];
rz(-1.4819179) q[2];
sx q[2];
rz(0.8366226) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72249216) q[1];
sx q[1];
rz(-2.3416391) q[1];
sx q[1];
rz(1.1833619) q[1];
x q[2];
rz(2.9929912) q[3];
sx q[3];
rz(-0.45431787) q[3];
sx q[3];
rz(-3.076072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2065108) q[2];
sx q[2];
rz(-0.93588459) q[2];
sx q[2];
rz(1.3235462) q[2];
rz(-0.27213085) q[3];
sx q[3];
rz(-2.2836253) q[3];
sx q[3];
rz(0.14565295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018933522) q[0];
sx q[0];
rz(-2.3441362) q[0];
sx q[0];
rz(-2.4830699) q[0];
rz(-1.5497442) q[1];
sx q[1];
rz(-1.0127944) q[1];
sx q[1];
rz(2.6351567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.63686) q[0];
sx q[0];
rz(-1.9397282) q[0];
sx q[0];
rz(1.641666) q[0];
rz(-pi) q[1];
rz(-1.5120878) q[2];
sx q[2];
rz(-0.48248267) q[2];
sx q[2];
rz(0.62472945) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5665019) q[1];
sx q[1];
rz(-0.12037863) q[1];
sx q[1];
rz(-0.68415227) q[1];
rz(1.5293838) q[3];
sx q[3];
rz(-0.99185252) q[3];
sx q[3];
rz(1.675663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3353614) q[2];
sx q[2];
rz(-0.74593097) q[2];
sx q[2];
rz(-1.8719505) q[2];
rz(-1.7757724) q[3];
sx q[3];
rz(-0.013805496) q[3];
sx q[3];
rz(-2.5891916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.7954471) q[0];
sx q[0];
rz(-0.2427635) q[0];
sx q[0];
rz(0.41918293) q[0];
rz(0.090713352) q[1];
sx q[1];
rz(-0.9181298) q[1];
sx q[1];
rz(1.027164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61495435) q[0];
sx q[0];
rz(-1.1246343) q[0];
sx q[0];
rz(0.63132186) q[0];
rz(-pi) q[1];
rz(-2.807051) q[2];
sx q[2];
rz(-0.024166232) q[2];
sx q[2];
rz(-2.6466359) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13090573) q[1];
sx q[1];
rz(-1.7831149) q[1];
sx q[1];
rz(-0.10879083) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6375362) q[3];
sx q[3];
rz(-0.68664521) q[3];
sx q[3];
rz(0.29717577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1586228) q[2];
sx q[2];
rz(-1.0047487) q[2];
sx q[2];
rz(0.99009222) q[2];
rz(-2.8236735) q[3];
sx q[3];
rz(-0.81219321) q[3];
sx q[3];
rz(0.038343553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7541499) q[0];
sx q[0];
rz(-2.4328572) q[0];
sx q[0];
rz(-2.5842174) q[0];
rz(0.36224657) q[1];
sx q[1];
rz(-1.4366415) q[1];
sx q[1];
rz(-0.16709669) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97826577) q[0];
sx q[0];
rz(-1.3045132) q[0];
sx q[0];
rz(-0.31727088) q[0];
x q[1];
rz(1.3873151) q[2];
sx q[2];
rz(-0.87786061) q[2];
sx q[2];
rz(-2.3132581) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5569589) q[1];
sx q[1];
rz(-1.8773729) q[1];
sx q[1];
rz(0.68309288) q[1];
rz(-pi) q[2];
rz(2.2542893) q[3];
sx q[3];
rz(-1.512882) q[3];
sx q[3];
rz(1.7660727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1780213) q[2];
sx q[2];
rz(-2.9489297) q[2];
sx q[2];
rz(2.7131405) q[2];
rz(-2.4109449) q[3];
sx q[3];
rz(-1.080039) q[3];
sx q[3];
rz(2.54125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927602) q[0];
sx q[0];
rz(-1.4880143) q[0];
sx q[0];
rz(-1.1195419) q[0];
rz(-0.66850942) q[1];
sx q[1];
rz(-2.5709277) q[1];
sx q[1];
rz(2.8817435) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90454532) q[0];
sx q[0];
rz(-2.0656842) q[0];
sx q[0];
rz(-0.2965692) q[0];
x q[1];
rz(-2.2462559) q[2];
sx q[2];
rz(-0.5731715) q[2];
sx q[2];
rz(-1.8942539) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6653241) q[1];
sx q[1];
rz(-0.55108738) q[1];
sx q[1];
rz(0.70161201) q[1];
rz(-pi) q[2];
rz(-2.9838461) q[3];
sx q[3];
rz(-2.254527) q[3];
sx q[3];
rz(-0.58303787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4708289) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(2.0551576) q[2];
rz(0.49232617) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(0.72211784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5928741) q[0];
sx q[0];
rz(-1.6328136) q[0];
sx q[0];
rz(1.6528224) q[0];
rz(-1.8863574) q[1];
sx q[1];
rz(-1.8240758) q[1];
sx q[1];
rz(0.12771894) q[1];
rz(1.3484501) q[2];
sx q[2];
rz(-2.7334474) q[2];
sx q[2];
rz(1.853142) q[2];
rz(-1.4896354) q[3];
sx q[3];
rz(-2.4774144) q[3];
sx q[3];
rz(-1.7437205) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
