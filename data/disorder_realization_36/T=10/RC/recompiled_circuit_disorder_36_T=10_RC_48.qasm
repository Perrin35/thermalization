OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(-2.6262992) q[0];
rz(1.3442858) q[1];
sx q[1];
rz(-2.9872515) q[1];
sx q[1];
rz(0.57758346) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89864697) q[0];
sx q[0];
rz(-2.4624914) q[0];
sx q[0];
rz(2.8773017) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0628113) q[2];
sx q[2];
rz(-2.0487818) q[2];
sx q[2];
rz(2.9302772) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9211728) q[1];
sx q[1];
rz(-0.49282679) q[1];
sx q[1];
rz(-2.1652031) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76831423) q[3];
sx q[3];
rz(-2.3294805) q[3];
sx q[3];
rz(2.1384359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.704533) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(0.68721592) q[2];
rz(1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9706443) q[0];
sx q[0];
rz(-1.0785372) q[0];
sx q[0];
rz(-1.8815536) q[0];
rz(-1.0062224) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(-2.2959183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082367912) q[0];
sx q[0];
rz(-1.1384283) q[0];
sx q[0];
rz(-0.57605497) q[0];
rz(-pi) q[1];
rz(0.13375608) q[2];
sx q[2];
rz(-1.4417366) q[2];
sx q[2];
rz(1.918902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.37296346) q[1];
sx q[1];
rz(-1.0781204) q[1];
sx q[1];
rz(-1.3586587) q[1];
rz(-2.2250697) q[3];
sx q[3];
rz(-2.407981) q[3];
sx q[3];
rz(-2.3678126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78850293) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(0.4804002) q[2];
rz(-1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(-3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-1.0864331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7074791) q[0];
sx q[0];
rz(-2.0322324) q[0];
sx q[0];
rz(-3.1397318) q[0];
x q[1];
rz(-0.92347446) q[2];
sx q[2];
rz(-1.662865) q[2];
sx q[2];
rz(-1.1084686) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.31049) q[1];
sx q[1];
rz(-1.5655787) q[1];
sx q[1];
rz(2.8527841) q[1];
x q[2];
rz(-2.9680786) q[3];
sx q[3];
rz(-1.1713235) q[3];
sx q[3];
rz(0.80843335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0125668) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(-1.770299) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(2.3390521) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(-0.048359811) q[0];
rz(2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(1.4455459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018054124) q[0];
sx q[0];
rz(-0.92045438) q[0];
sx q[0];
rz(1.7975848) q[0];
rz(-pi) q[1];
rz(0.21638685) q[2];
sx q[2];
rz(-1.4828223) q[2];
sx q[2];
rz(2.3073334) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9963035) q[1];
sx q[1];
rz(-2.9128296) q[1];
sx q[1];
rz(2.8537675) q[1];
rz(-pi) q[2];
rz(2.0914145) q[3];
sx q[3];
rz(-1.8225267) q[3];
sx q[3];
rz(0.43581918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.066102862) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(-2.1172822) q[2];
rz(-1.5284437) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(-2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4399453) q[0];
sx q[0];
rz(-2.3174536) q[0];
sx q[0];
rz(1.2874999) q[0];
rz(-0.31907407) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(0.85420001) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720649) q[0];
sx q[0];
rz(-2.2846691) q[0];
sx q[0];
rz(3.1307427) q[0];
x q[1];
rz(-1.008026) q[2];
sx q[2];
rz(-2.3587583) q[2];
sx q[2];
rz(0.40751878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2367868) q[1];
sx q[1];
rz(-0.85610897) q[1];
sx q[1];
rz(2.3805815) q[1];
x q[2];
rz(-1.0075931) q[3];
sx q[3];
rz(-1.0922722) q[3];
sx q[3];
rz(-1.0790881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.95191082) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(-2.5642776) q[2];
rz(2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(3.1381056) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(3.0694718) q[0];
rz(2.0347319) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(3.0153826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0059347) q[0];
sx q[0];
rz(-1.5392443) q[0];
sx q[0];
rz(0.21261442) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3665479) q[2];
sx q[2];
rz(-1.8473052) q[2];
sx q[2];
rz(-2.0691878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5050161) q[1];
sx q[1];
rz(-0.86537433) q[1];
sx q[1];
rz(-0.24800639) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9494709) q[3];
sx q[3];
rz(-1.2479094) q[3];
sx q[3];
rz(-2.9411112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(-0.827968) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5489952) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(-0.098408498) q[0];
rz(-1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-0.55955204) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37724272) q[0];
sx q[0];
rz(-2.5284405) q[0];
sx q[0];
rz(-2.9186547) q[0];
rz(-pi) q[1];
rz(1.6227116) q[2];
sx q[2];
rz(-1.3097109) q[2];
sx q[2];
rz(1.3494929) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0911078) q[1];
sx q[1];
rz(-2.9851966) q[1];
sx q[1];
rz(-0.97538235) q[1];
x q[2];
rz(-0.26569326) q[3];
sx q[3];
rz(-0.64722792) q[3];
sx q[3];
rz(-2.2006187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6825535) q[2];
sx q[2];
rz(-1.8457396) q[2];
sx q[2];
rz(2.7977978) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664117) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(0.73076105) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(2.2699845) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9974737) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(0.075450443) q[0];
rz(-2.7360104) q[2];
sx q[2];
rz(-1.7103346) q[2];
sx q[2];
rz(1.1428733) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.65054446) q[1];
sx q[1];
rz(-2.3666413) q[1];
sx q[1];
rz(-2.6130555) q[1];
x q[2];
rz(-1.7730764) q[3];
sx q[3];
rz(-2.5296122) q[3];
sx q[3];
rz(0.66093854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72620755) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(-2.0020206) q[2];
rz(1.4987882) q[3];
sx q[3];
rz(-2.747624) q[3];
sx q[3];
rz(-0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20294872) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(-1.2217481) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(-1.5244012) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.411392) q[0];
sx q[0];
rz(-3.054266) q[0];
sx q[0];
rz(-1.2325531) q[0];
rz(1.379307) q[2];
sx q[2];
rz(-1.101149) q[2];
sx q[2];
rz(0.42275235) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2251687) q[1];
sx q[1];
rz(-1.7290219) q[1];
sx q[1];
rz(1.2374864) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41096656) q[3];
sx q[3];
rz(-0.68115679) q[3];
sx q[3];
rz(2.4364803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1200072) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(-0.35153708) q[2];
rz(-2.0848138) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(-0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64365023) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(1.836401) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(2.8881853) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5288552) q[0];
sx q[0];
rz(-2.7326072) q[0];
sx q[0];
rz(-1.3091875) q[0];
x q[1];
rz(1.1240187) q[2];
sx q[2];
rz(-2.7687216) q[2];
sx q[2];
rz(-0.98137059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7518172) q[1];
sx q[1];
rz(-0.42418617) q[1];
sx q[1];
rz(0.61702375) q[1];
rz(-pi) q[2];
rz(0.6110544) q[3];
sx q[3];
rz(-1.9285893) q[3];
sx q[3];
rz(-2.4886481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5499251) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(-2.2425966) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(-0.7111711) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(-2.3989427) q[2];
sx q[2];
rz(-2.7682318) q[2];
sx q[2];
rz(-2.9329185) q[2];
rz(0.75541227) q[3];
sx q[3];
rz(-1.0676386) q[3];
sx q[3];
rz(-0.044015351) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
