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
rz(0.69005203) q[0];
sx q[0];
rz(-0.85059387) q[0];
sx q[0];
rz(-2.9414862) q[0];
rz(-0.19835681) q[1];
sx q[1];
rz(5.7601647) q[1];
sx q[1];
rz(10.756607) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9560741) q[0];
sx q[0];
rz(-1.8022493) q[0];
sx q[0];
rz(2.0218019) q[0];
rz(2.0609786) q[2];
sx q[2];
rz(-1.0287544) q[2];
sx q[2];
rz(1.3203743) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.15439776) q[1];
sx q[1];
rz(-1.5227277) q[1];
sx q[1];
rz(2.4632719) q[1];
rz(1.0947919) q[3];
sx q[3];
rz(-1.6443569) q[3];
sx q[3];
rz(1.5293763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0693822) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(-1.3145831) q[2];
rz(-2.2383111) q[3];
sx q[3];
rz(-1.3638146) q[3];
sx q[3];
rz(-2.7956853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8697934) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(-0.85579175) q[0];
rz(2.3675512) q[1];
sx q[1];
rz(-1.9776521) q[1];
sx q[1];
rz(2.3537297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5408527) q[0];
sx q[0];
rz(-2.8898281) q[0];
sx q[0];
rz(-1.4815848) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.79736) q[2];
sx q[2];
rz(-1.6107127) q[2];
sx q[2];
rz(2.0502063) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98709244) q[1];
sx q[1];
rz(-0.56114158) q[1];
sx q[1];
rz(-0.89548703) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66707261) q[3];
sx q[3];
rz(-1.7366341) q[3];
sx q[3];
rz(1.9465684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.065980109) q[2];
sx q[2];
rz(-1.2300666) q[2];
sx q[2];
rz(2.7454929) q[2];
rz(-1.5710477) q[3];
sx q[3];
rz(-1.5038265) q[3];
sx q[3];
rz(-1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1876672) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(-3.0809825) q[0];
rz(-2.4568779) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(0.33904591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.402917) q[0];
sx q[0];
rz(-2.5581723) q[0];
sx q[0];
rz(0.03761272) q[0];
x q[1];
rz(-1.4475559) q[2];
sx q[2];
rz(-1.7584821) q[2];
sx q[2];
rz(-2.0579684) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10968929) q[1];
sx q[1];
rz(-1.8432143) q[1];
sx q[1];
rz(-1.2485571) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77683461) q[3];
sx q[3];
rz(-1.562444) q[3];
sx q[3];
rz(1.9296822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6591349) q[2];
sx q[2];
rz(-1.8637916) q[2];
sx q[2];
rz(-2.6711312) q[2];
rz(-0.86236924) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7793133) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(2.733574) q[0];
rz(1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(1.8203576) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119638) q[0];
sx q[0];
rz(-0.59998673) q[0];
sx q[0];
rz(2.5866051) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2259813) q[2];
sx q[2];
rz(-1.639099) q[2];
sx q[2];
rz(-0.11924041) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27038747) q[1];
sx q[1];
rz(-0.14158881) q[1];
sx q[1];
rz(2.6006727) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3636171) q[3];
sx q[3];
rz(-2.5730592) q[3];
sx q[3];
rz(1.3028627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9185751) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(-2.4187386) q[2];
rz(-0.84960788) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.9549114) q[0];
sx q[0];
rz(-2.6246922) q[0];
sx q[0];
rz(2.5816259) q[0];
rz(0.2050744) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(-0.29108873) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2084853) q[0];
sx q[0];
rz(-1.7400842) q[0];
sx q[0];
rz(-2.8693958) q[0];
rz(0.061441378) q[2];
sx q[2];
rz(-2.3413918) q[2];
sx q[2];
rz(-0.35129181) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0138058) q[1];
sx q[1];
rz(-1.769279) q[1];
sx q[1];
rz(-0.86680331) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7311312) q[3];
sx q[3];
rz(-0.82538) q[3];
sx q[3];
rz(-1.3813409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9186972) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(2.843294) q[2];
rz(-1.9722624) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30763141) q[0];
sx q[0];
rz(-1.6809604) q[0];
sx q[0];
rz(-0.6024012) q[0];
rz(0.3715474) q[1];
sx q[1];
rz(-1.5448152) q[1];
sx q[1];
rz(2.1098302) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4435972) q[0];
sx q[0];
rz(-1.5549608) q[0];
sx q[0];
rz(1.7450733) q[0];
rz(-pi) q[1];
rz(-0.4571896) q[2];
sx q[2];
rz(-0.44538272) q[2];
sx q[2];
rz(3.074844) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6556485) q[1];
sx q[1];
rz(-1.7595425) q[1];
sx q[1];
rz(0.48635482) q[1];
rz(-pi) q[2];
rz(2.2797645) q[3];
sx q[3];
rz(-2.3579881) q[3];
sx q[3];
rz(2.6605822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.80456698) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(1.6229013) q[2];
rz(-2.2931781) q[3];
sx q[3];
rz(-0.87526667) q[3];
sx q[3];
rz(3.1019822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1161716) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(0.99789944) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(1.7707228) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28579636) q[0];
sx q[0];
rz(-1.5573475) q[0];
sx q[0];
rz(0.6310668) q[0];
rz(-pi) q[1];
rz(0.15372686) q[2];
sx q[2];
rz(-1.9228553) q[2];
sx q[2];
rz(-1.597126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6295885) q[1];
sx q[1];
rz(-1.3301714) q[1];
sx q[1];
rz(-2.5179067) q[1];
rz(-pi) q[2];
rz(0.73565817) q[3];
sx q[3];
rz(-2.1950985) q[3];
sx q[3];
rz(-1.5746547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85713282) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(0.029646309) q[2];
rz(-0.61521411) q[3];
sx q[3];
rz(-1.3968202) q[3];
sx q[3];
rz(-2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0520332) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(0.7830559) q[0];
rz(2.0217333) q[1];
sx q[1];
rz(-0.66711396) q[1];
sx q[1];
rz(-1.3979744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2180041) q[0];
sx q[0];
rz(-0.62380416) q[0];
sx q[0];
rz(0.82360928) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2863761) q[2];
sx q[2];
rz(-1.3462974) q[2];
sx q[2];
rz(2.7530991) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82247558) q[1];
sx q[1];
rz(-2.2386595) q[1];
sx q[1];
rz(-3.130359) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3572631) q[3];
sx q[3];
rz(-1.4186067) q[3];
sx q[3];
rz(-2.1209716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8457501) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(1.3207377) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9418697) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(2.1871908) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(-1.0844213) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53756489) q[0];
sx q[0];
rz(-1.5286235) q[0];
sx q[0];
rz(2.2794363) q[0];
rz(1.8516638) q[2];
sx q[2];
rz(-1.3251588) q[2];
sx q[2];
rz(0.40058595) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22218765) q[1];
sx q[1];
rz(-2.0424941) q[1];
sx q[1];
rz(-3.0604048) q[1];
rz(1.3865269) q[3];
sx q[3];
rz(-0.83067719) q[3];
sx q[3];
rz(-2.4746462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.083100975) q[2];
sx q[2];
rz(-0.19187555) q[2];
sx q[2];
rz(2.7044738) q[2];
rz(1.7982177) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(0.04537151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042353543) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(1.4105256) q[0];
rz(0.22077416) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(-2.208362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8659023) q[0];
sx q[0];
rz(-2.3161016) q[0];
sx q[0];
rz(1.5522214) q[0];
rz(-pi) q[1];
rz(0.61874091) q[2];
sx q[2];
rz(-1.8685942) q[2];
sx q[2];
rz(-2.1190475) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.90551841) q[1];
sx q[1];
rz(-0.80414861) q[1];
sx q[1];
rz(0.18412904) q[1];
x q[2];
rz(2.203412) q[3];
sx q[3];
rz(-1.4804192) q[3];
sx q[3];
rz(2.602102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.431388) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(-0.63146511) q[2];
rz(-2.4225875) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(-2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63856335) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(-0.44166625) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(-1.1900525) q[2];
sx q[2];
rz(-1.5413918) q[2];
sx q[2];
rz(-0.96985758) q[2];
rz(-2.9819103) q[3];
sx q[3];
rz(-0.57542141) q[3];
sx q[3];
rz(0.66948359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
