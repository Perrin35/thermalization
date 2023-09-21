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
rz(0.51529348) q[0];
rz(4.4858785) q[1];
sx q[1];
rz(2.9872515) q[1];
sx q[1];
rz(6.8607688) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89864697) q[0];
sx q[0];
rz(-0.67910128) q[0];
sx q[0];
rz(0.26429096) q[0];
rz(-0.53519997) q[2];
sx q[2];
rz(-1.1241962) q[2];
sx q[2];
rz(1.5314147) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2540993) q[1];
sx q[1];
rz(-1.8389529) q[1];
sx q[1];
rz(1.9894132) q[1];
x q[2];
rz(2.4926315) q[3];
sx q[3];
rz(-2.0994086) q[3];
sx q[3];
rz(-1.1543857) q[3];
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
rz(2.1263188) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(3.0190873) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(1.2600391) q[0];
rz(-2.1353703) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(0.84567436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082367912) q[0];
sx q[0];
rz(-1.1384283) q[0];
sx q[0];
rz(2.5655377) q[0];
x q[1];
rz(3.0078366) q[2];
sx q[2];
rz(-1.6998561) q[2];
sx q[2];
rz(1.918902) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37296346) q[1];
sx q[1];
rz(-1.0781204) q[1];
sx q[1];
rz(1.3586587) q[1];
x q[2];
rz(0.94988471) q[3];
sx q[3];
rz(-1.1511027) q[3];
sx q[3];
rz(1.3148395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78850293) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(-2.6611924) q[2];
rz(-1.7885615) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(-1.9539179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95124328) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(2.3685266) q[0];
rz(-3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(2.0551596) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7032996) q[0];
sx q[0];
rz(-0.46143954) q[0];
sx q[0];
rz(-1.5745387) q[0];
rz(-0.11523192) q[2];
sx q[2];
rz(-0.92667246) q[2];
sx q[2];
rz(2.6098721) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4003488) q[1];
sx q[1];
rz(-1.2819918) q[1];
sx q[1];
rz(1.5762394) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1658737) q[3];
sx q[3];
rz(-1.4110663) q[3];
sx q[3];
rz(2.3111642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0125668) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(-1.770299) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(0.80254054) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(0.048359811) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(1.6960467) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915582) q[0];
sx q[0];
rz(-1.3908747) q[0];
sx q[0];
rz(2.4787089) q[0];
rz(-pi) q[1];
rz(-1.4807329) q[2];
sx q[2];
rz(-1.7863331) q[2];
sx q[2];
rz(-2.4243674) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1452892) q[1];
sx q[1];
rz(-2.9128296) q[1];
sx q[1];
rz(-2.8537675) q[1];
rz(-pi) q[2];
rz(1.0501782) q[3];
sx q[3];
rz(-1.319066) q[3];
sx q[3];
rz(0.43581918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.066102862) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(-2.1172822) q[2];
rz(1.6131489) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(-2.9083692) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.8540927) q[0];
rz(-0.31907407) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(-0.85420001) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.655495) q[0];
sx q[0];
rz(-2.4276519) q[0];
sx q[0];
rz(1.5833202) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6536077) q[2];
sx q[2];
rz(-2.2099566) q[2];
sx q[2];
rz(1.1346863) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2367868) q[1];
sx q[1];
rz(-2.2854837) q[1];
sx q[1];
rz(2.3805815) q[1];
x q[2];
rz(2.5913127) q[3];
sx q[3];
rz(-1.0770505) q[3];
sx q[3];
rz(0.20875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.95191082) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(2.5642776) q[2];
rz(-0.50950766) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(2.234941) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0034870738) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(0.072120897) q[0];
rz(-2.0347319) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(3.0153826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5835411) q[0];
sx q[0];
rz(-1.3582894) q[0];
sx q[0];
rz(-1.538518) q[0];
rz(-2.7563165) q[2];
sx q[2];
rz(-0.81309536) q[2];
sx q[2];
rz(0.77020459) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26423745) q[1];
sx q[1];
rz(-2.4009631) q[1];
sx q[1];
rz(-1.8514368) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0522271) q[3];
sx q[3];
rz(-0.37399451) q[3];
sx q[3];
rz(0.75043375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90298992) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(2.3664756) q[2];
rz(2.3136247) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(-1.1221788) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59259748) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(3.0431842) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-2.5820406) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0347621) q[0];
sx q[0];
rz(-2.1666399) q[0];
sx q[0];
rz(-1.4164657) q[0];
x q[1];
rz(0.19182972) q[2];
sx q[2];
rz(-2.8755113) q[2];
sx q[2];
rz(-1.9907469) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0911078) q[1];
sx q[1];
rz(-2.9851966) q[1];
sx q[1];
rz(-0.97538235) q[1];
x q[2];
rz(0.63013245) q[3];
sx q[3];
rz(-1.7297941) q[3];
sx q[3];
rz(-0.4160479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45903912) q[2];
sx q[2];
rz(-1.8457396) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664117) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-2.4108316) q[0];
rz(0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(2.2699845) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14411892) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(0.075450443) q[0];
rz(1.7224738) q[2];
sx q[2];
rz(-1.9722087) q[2];
sx q[2];
rz(0.36827189) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4910482) q[1];
sx q[1];
rz(-2.3666413) q[1];
sx q[1];
rz(2.6130555) q[1];
rz(0.96846795) q[3];
sx q[3];
rz(-1.6864711) q[3];
sx q[3];
rz(-1.0761716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.72620755) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(-2.0020206) q[2];
rz(-1.6428044) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.9386439) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(1.2217481) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(1.6171914) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4964543) q[0];
sx q[0];
rz(-1.5997412) q[0];
sx q[0];
rz(-1.653198) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7828214) q[2];
sx q[2];
rz(-2.6371187) q[2];
sx q[2];
rz(2.3141253) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0601378) q[1];
sx q[1];
rz(-0.36768915) q[1];
sx q[1];
rz(-1.1170438) q[1];
rz(-pi) q[2];
rz(1.8839621) q[3];
sx q[3];
rz(-2.1861665) q[3];
sx q[3];
rz(-2.9477011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1200072) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(-2.7900556) q[2];
rz(-1.0567788) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(-2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4979424) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.836401) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(-0.25340733) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5288552) q[0];
sx q[0];
rz(-0.40898541) q[0];
sx q[0];
rz(1.3091875) q[0];
rz(2.9741653) q[2];
sx q[2];
rz(-1.9055467) q[2];
sx q[2];
rz(1.4565005) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7518172) q[1];
sx q[1];
rz(-0.42418617) q[1];
sx q[1];
rz(0.61702375) q[1];
rz(-pi) q[2];
rz(-2.5640423) q[3];
sx q[3];
rz(-0.69637075) q[3];
sx q[3];
rz(-1.7601354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59166756) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-0.5029451) q[2];
rz(-0.89899603) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(-1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9713365) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(2.4304216) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(2.3989427) q[2];
sx q[2];
rz(-0.37336083) q[2];
sx q[2];
rz(0.20867418) q[2];
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