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
rz(-0.51802975) q[0];
sx q[0];
rz(1.1413483) q[0];
sx q[0];
rz(9.4649014) q[0];
rz(0.24815458) q[1];
sx q[1];
rz(5.2207898) q[1];
sx q[1];
rz(9.3407486) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87977982) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(-0.69095822) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7438091) q[2];
sx q[2];
rz(-1.7353463) q[2];
sx q[2];
rz(1.0918504) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1511615) q[1];
sx q[1];
rz(-0.50499454) q[1];
sx q[1];
rz(2.7846365) q[1];
rz(-pi) q[2];
rz(2.5474882) q[3];
sx q[3];
rz(-0.66451529) q[3];
sx q[3];
rz(2.6474109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6930801) q[2];
sx q[2];
rz(-1.3658407) q[2];
sx q[2];
rz(0.25126323) q[2];
rz(1.1723899) q[3];
sx q[3];
rz(-1.1029714) q[3];
sx q[3];
rz(-3.0917061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6867111) q[0];
sx q[0];
rz(-2.659681) q[0];
sx q[0];
rz(1.7676109) q[0];
rz(2.8107367) q[1];
sx q[1];
rz(-1.7287946) q[1];
sx q[1];
rz(-2.8673598) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5170256) q[0];
sx q[0];
rz(-2.8170878) q[0];
sx q[0];
rz(-2.9175678) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5636257) q[2];
sx q[2];
rz(-1.0375377) q[2];
sx q[2];
rz(-0.28927848) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3220904) q[1];
sx q[1];
rz(-1.7871457) q[1];
sx q[1];
rz(2.809193) q[1];
rz(1.1521454) q[3];
sx q[3];
rz(-1.6466221) q[3];
sx q[3];
rz(2.0328498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8809044) q[2];
sx q[2];
rz(-1.3284677) q[2];
sx q[2];
rz(-1.0880967) q[2];
rz(-1.043383) q[3];
sx q[3];
rz(-0.16845307) q[3];
sx q[3];
rz(-0.25922957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.008721) q[0];
sx q[0];
rz(-0.50478029) q[0];
sx q[0];
rz(-2.0399427) q[0];
rz(0.77611008) q[1];
sx q[1];
rz(-1.848369) q[1];
sx q[1];
rz(2.4993842) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74321824) q[0];
sx q[0];
rz(-1.9954329) q[0];
sx q[0];
rz(-1.0845127) q[0];
rz(-pi) q[1];
rz(1.3234336) q[2];
sx q[2];
rz(-2.2752834) q[2];
sx q[2];
rz(3.041628) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57802618) q[1];
sx q[1];
rz(-1.1456923) q[1];
sx q[1];
rz(0.11386392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65861838) q[3];
sx q[3];
rz(-2.082654) q[3];
sx q[3];
rz(0.030125253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8316101) q[2];
sx q[2];
rz(-1.9881366) q[2];
sx q[2];
rz(-2.2115121) q[2];
rz(-2.3275404) q[3];
sx q[3];
rz(-1.9508773) q[3];
sx q[3];
rz(1.2306151) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.880068) q[0];
sx q[0];
rz(-1.1898758) q[0];
sx q[0];
rz(0.50450182) q[0];
rz(-0.91744676) q[1];
sx q[1];
rz(-2.4023299) q[1];
sx q[1];
rz(2.9333072) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.298131) q[0];
sx q[0];
rz(-1.4220474) q[0];
sx q[0];
rz(-2.896033) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4217626) q[2];
sx q[2];
rz(-1.3583169) q[2];
sx q[2];
rz(1.8393387) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54305824) q[1];
sx q[1];
rz(-2.0428791) q[1];
sx q[1];
rz(-0.92453875) q[1];
rz(-pi) q[2];
rz(0.55996079) q[3];
sx q[3];
rz(-1.3508288) q[3];
sx q[3];
rz(0.44693963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8248799) q[2];
sx q[2];
rz(-2.0471639) q[2];
sx q[2];
rz(1.1701976) q[2];
rz(-2.1757388) q[3];
sx q[3];
rz(-1.0713157) q[3];
sx q[3];
rz(-2.9962208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.589094) q[0];
sx q[0];
rz(-0.67411244) q[0];
sx q[0];
rz(-2.476995) q[0];
rz(0.26847863) q[1];
sx q[1];
rz(-1.45739) q[1];
sx q[1];
rz(-2.1427515) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12832175) q[0];
sx q[0];
rz(-1.4185462) q[0];
sx q[0];
rz(0.98441846) q[0];
x q[1];
rz(-0.13370698) q[2];
sx q[2];
rz(-1.1503476) q[2];
sx q[2];
rz(-2.9131817) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4722574) q[1];
sx q[1];
rz(-0.20755033) q[1];
sx q[1];
rz(2.186354) q[1];
rz(-3.1206297) q[3];
sx q[3];
rz(-2.2624894) q[3];
sx q[3];
rz(1.616459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1642248) q[2];
sx q[2];
rz(-2.2013142) q[2];
sx q[2];
rz(0.11717907) q[2];
rz(-0.059727877) q[3];
sx q[3];
rz(-1.2195339) q[3];
sx q[3];
rz(-0.81502771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77434671) q[0];
sx q[0];
rz(-0.4137488) q[0];
sx q[0];
rz(1.2600979) q[0];
rz(0.11481181) q[1];
sx q[1];
rz(-1.3453307) q[1];
sx q[1];
rz(0.095349163) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0336825) q[0];
sx q[0];
rz(-2.9386407) q[0];
sx q[0];
rz(-3.0686707) q[0];
rz(-pi) q[1];
rz(-1.0667858) q[2];
sx q[2];
rz(-1.7169478) q[2];
sx q[2];
rz(2.604217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5231875) q[1];
sx q[1];
rz(-1.1506936) q[1];
sx q[1];
rz(2.0610496) q[1];
rz(1.5524149) q[3];
sx q[3];
rz(-0.52747969) q[3];
sx q[3];
rz(-2.7888586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0872588) q[2];
sx q[2];
rz(-2.2237873) q[2];
sx q[2];
rz(-0.66162649) q[2];
rz(-2.3719487) q[3];
sx q[3];
rz(-1.6911643) q[3];
sx q[3];
rz(2.2746287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0041644) q[0];
sx q[0];
rz(-0.34808174) q[0];
sx q[0];
rz(-2.3837756) q[0];
rz(-3.1163395) q[1];
sx q[1];
rz(-2.7460637) q[1];
sx q[1];
rz(1.9662366) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1036183) q[0];
sx q[0];
rz(-2.0424065) q[0];
sx q[0];
rz(-1.1704117) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98837672) q[2];
sx q[2];
rz(-2.3727131) q[2];
sx q[2];
rz(-1.066756) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3010294) q[1];
sx q[1];
rz(-1.9604654) q[1];
sx q[1];
rz(0.79857608) q[1];
x q[2];
rz(-0.42899363) q[3];
sx q[3];
rz(-1.5499664) q[3];
sx q[3];
rz(-0.57383895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5777099) q[2];
sx q[2];
rz(-1.6174199) q[2];
sx q[2];
rz(-1.9263402) q[2];
rz(0.75872129) q[3];
sx q[3];
rz(-1.7251549) q[3];
sx q[3];
rz(1.5569713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4828846) q[0];
sx q[0];
rz(-1.8526798) q[0];
sx q[0];
rz(-0.38920745) q[0];
rz(3.0534741) q[1];
sx q[1];
rz(-1.7033109) q[1];
sx q[1];
rz(1.5315936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54231468) q[0];
sx q[0];
rz(-1.3345584) q[0];
sx q[0];
rz(-1.3482984) q[0];
rz(-0.56014748) q[2];
sx q[2];
rz(-0.74290007) q[2];
sx q[2];
rz(-0.35336885) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1080146) q[1];
sx q[1];
rz(-0.44541767) q[1];
sx q[1];
rz(-3.0207795) q[1];
rz(2.2810488) q[3];
sx q[3];
rz(-1.2757974) q[3];
sx q[3];
rz(-0.33911946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4072998) q[2];
sx q[2];
rz(-0.6137085) q[2];
sx q[2];
rz(1.8180004) q[2];
rz(1.7981516) q[3];
sx q[3];
rz(-0.7898134) q[3];
sx q[3];
rz(2.8048803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6650498) q[0];
sx q[0];
rz(-2.2344868) q[0];
sx q[0];
rz(-2.6764349) q[0];
rz(3.0910885) q[1];
sx q[1];
rz(-0.89021325) q[1];
sx q[1];
rz(0.49096289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4024593) q[0];
sx q[0];
rz(-1.2417034) q[0];
sx q[0];
rz(-1.9080481) q[0];
rz(-pi) q[1];
rz(-0.30740909) q[2];
sx q[2];
rz(-1.4239862) q[2];
sx q[2];
rz(1.6830144) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95066324) q[1];
sx q[1];
rz(-1.0196389) q[1];
sx q[1];
rz(0.74857708) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6816986) q[3];
sx q[3];
rz(-1.2909799) q[3];
sx q[3];
rz(-0.96074897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.65840536) q[2];
sx q[2];
rz(-0.5961954) q[2];
sx q[2];
rz(-2.2779706) q[2];
rz(-2.9577799) q[3];
sx q[3];
rz(-2.176216) q[3];
sx q[3];
rz(0.98627728) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7154536) q[0];
sx q[0];
rz(-1.5197536) q[0];
sx q[0];
rz(-2.5157628) q[0];
rz(-0.65678701) q[1];
sx q[1];
rz(-0.96581179) q[1];
sx q[1];
rz(-1.9331369) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10851181) q[0];
sx q[0];
rz(-1.7038561) q[0];
sx q[0];
rz(-0.10255524) q[0];
x q[1];
rz(2.7371128) q[2];
sx q[2];
rz(-2.4411429) q[2];
sx q[2];
rz(3.1062448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6555935) q[1];
sx q[1];
rz(-1.5734377) q[1];
sx q[1];
rz(3.0551617) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5700611) q[3];
sx q[3];
rz(-1.4213955) q[3];
sx q[3];
rz(-0.3954486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5829696) q[2];
sx q[2];
rz(-1.8331336) q[2];
sx q[2];
rz(0.51897955) q[2];
rz(-0.44273043) q[3];
sx q[3];
rz(-1.818592) q[3];
sx q[3];
rz(1.7980661) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9291572) q[0];
sx q[0];
rz(-2.1514308) q[0];
sx q[0];
rz(1.9707752) q[0];
rz(0.15405542) q[1];
sx q[1];
rz(-2.2046721) q[1];
sx q[1];
rz(2.2278723) q[1];
rz(2.6250203) q[2];
sx q[2];
rz(-1.5348829) q[2];
sx q[2];
rz(1.665967) q[2];
rz(2.7967706) q[3];
sx q[3];
rz(-2.214684) q[3];
sx q[3];
rz(-0.82292258) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
