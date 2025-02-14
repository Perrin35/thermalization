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
rz(0.12953144) q[0];
sx q[0];
rz(-0.13106267) q[0];
sx q[0];
rz(2.5370606) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(2.8837535) q[1];
sx q[1];
rz(16.937994) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2211232) q[0];
sx q[0];
rz(-1.2386432) q[0];
sx q[0];
rz(-1.3060547) q[0];
rz(-pi) q[1];
rz(-2.8240439) q[2];
sx q[2];
rz(-2.5032836) q[2];
sx q[2];
rz(2.2090863) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9046272) q[1];
sx q[1];
rz(-0.89348999) q[1];
sx q[1];
rz(0.64579247) q[1];
rz(-0.12768605) q[3];
sx q[3];
rz(-2.0632072) q[3];
sx q[3];
rz(-3.1000053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0543542) q[2];
sx q[2];
rz(-1.1776935) q[2];
sx q[2];
rz(0.13620201) q[2];
rz(2.6916091) q[3];
sx q[3];
rz(-2.4583702) q[3];
sx q[3];
rz(2.3581678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1012652) q[0];
sx q[0];
rz(-2.3430921) q[0];
sx q[0];
rz(0.26327565) q[0];
rz(2.1030262) q[1];
sx q[1];
rz(-0.34663215) q[1];
sx q[1];
rz(2.3668049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8717125) q[0];
sx q[0];
rz(-0.96202606) q[0];
sx q[0];
rz(-1.3193498) q[0];
rz(2.4500174) q[2];
sx q[2];
rz(-1.0333158) q[2];
sx q[2];
rz(2.7730758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6331433) q[1];
sx q[1];
rz(-1.6078498) q[1];
sx q[1];
rz(-1.1913185) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70330422) q[3];
sx q[3];
rz(-1.9580132) q[3];
sx q[3];
rz(-0.73782286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7210377) q[2];
sx q[2];
rz(-1.6246395) q[2];
sx q[2];
rz(2.5431385) q[2];
rz(1.362644) q[3];
sx q[3];
rz(-2.2612031) q[3];
sx q[3];
rz(2.8708598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8552928) q[0];
sx q[0];
rz(-2.6787651) q[0];
sx q[0];
rz(-1.5077952) q[0];
rz(0.15448054) q[1];
sx q[1];
rz(-1.721761) q[1];
sx q[1];
rz(-0.78937626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75784439) q[0];
sx q[0];
rz(-1.0417291) q[0];
sx q[0];
rz(-1.4704513) q[0];
rz(1.71307) q[2];
sx q[2];
rz(-1.8492438) q[2];
sx q[2];
rz(1.5253893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7676413) q[1];
sx q[1];
rz(-2.1726296) q[1];
sx q[1];
rz(2.7695719) q[1];
rz(-pi) q[2];
rz(-2.6487938) q[3];
sx q[3];
rz(-2.2792553) q[3];
sx q[3];
rz(0.81564834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35884759) q[2];
sx q[2];
rz(-2.2300356) q[2];
sx q[2];
rz(1.4770799) q[2];
rz(-1.9278256) q[3];
sx q[3];
rz(-1.8002847) q[3];
sx q[3];
rz(1.414813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7724991) q[0];
sx q[0];
rz(-1.0558244) q[0];
sx q[0];
rz(2.5308727) q[0];
rz(-0.67131132) q[1];
sx q[1];
rz(-1.6339615) q[1];
sx q[1];
rz(-1.3549365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63263921) q[0];
sx q[0];
rz(-2.0193978) q[0];
sx q[0];
rz(2.6030356) q[0];
x q[1];
rz(-1.1995537) q[2];
sx q[2];
rz(-1.4535731) q[2];
sx q[2];
rz(0.10088149) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0555041) q[1];
sx q[1];
rz(-1.633606) q[1];
sx q[1];
rz(-1.8163866) q[1];
rz(3.0464744) q[3];
sx q[3];
rz(-0.98604938) q[3];
sx q[3];
rz(-2.1383536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9352202) q[2];
sx q[2];
rz(-2.361203) q[2];
sx q[2];
rz(-2.8727403) q[2];
rz(2.3937461) q[3];
sx q[3];
rz(-0.81955376) q[3];
sx q[3];
rz(-1.4486754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95626962) q[0];
sx q[0];
rz(-1.3893501) q[0];
sx q[0];
rz(0.67778936) q[0];
rz(-0.14450821) q[1];
sx q[1];
rz(-2.3542207) q[1];
sx q[1];
rz(1.6544624) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22600284) q[0];
sx q[0];
rz(-1.5446413) q[0];
sx q[0];
rz(1.4348861) q[0];
rz(1.9727835) q[2];
sx q[2];
rz(-1.5176519) q[2];
sx q[2];
rz(1.0099908) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2320391) q[1];
sx q[1];
rz(-1.1097481) q[1];
sx q[1];
rz(2.1559245) q[1];
rz(-pi) q[2];
x q[2];
rz(1.946536) q[3];
sx q[3];
rz(-1.5259966) q[3];
sx q[3];
rz(0.57023772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3378478) q[2];
sx q[2];
rz(-0.2636815) q[2];
sx q[2];
rz(2.7860876) q[2];
rz(1.0454987) q[3];
sx q[3];
rz(-1.239536) q[3];
sx q[3];
rz(0.56226468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97583714) q[0];
sx q[0];
rz(-2.5582357) q[0];
sx q[0];
rz(-3.052886) q[0];
rz(-1.6070131) q[1];
sx q[1];
rz(-1.3101703) q[1];
sx q[1];
rz(-0.075693695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82268084) q[0];
sx q[0];
rz(-1.3782129) q[0];
sx q[0];
rz(-2.4680424) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83916243) q[2];
sx q[2];
rz(-1.5328836) q[2];
sx q[2];
rz(2.1597852) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5548426) q[1];
sx q[1];
rz(-1.7126709) q[1];
sx q[1];
rz(2.1290522) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47786062) q[3];
sx q[3];
rz(-0.63242542) q[3];
sx q[3];
rz(2.6101108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3169516) q[2];
sx q[2];
rz(-1.3836766) q[2];
sx q[2];
rz(2.1023882) q[2];
rz(2.9292987) q[3];
sx q[3];
rz(-1.1024691) q[3];
sx q[3];
rz(1.4958517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742663) q[0];
sx q[0];
rz(-0.75263158) q[0];
sx q[0];
rz(2.0972032) q[0];
rz(0.77230612) q[1];
sx q[1];
rz(-2.4817395) q[1];
sx q[1];
rz(-0.73371249) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6873467) q[0];
sx q[0];
rz(-0.92763072) q[0];
sx q[0];
rz(0.9439133) q[0];
rz(1.9416503) q[2];
sx q[2];
rz(-2.7463253) q[2];
sx q[2];
rz(-2.7140009) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4419879) q[1];
sx q[1];
rz(-1.8612222) q[1];
sx q[1];
rz(-1.6922349) q[1];
rz(0.83603199) q[3];
sx q[3];
rz(-1.6510909) q[3];
sx q[3];
rz(-1.3463595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5250728) q[2];
sx q[2];
rz(-0.51023444) q[2];
sx q[2];
rz(2.5329242) q[2];
rz(-2.0742553) q[3];
sx q[3];
rz(-2.267024) q[3];
sx q[3];
rz(2.5909891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6390425) q[0];
sx q[0];
rz(-2.783343) q[0];
sx q[0];
rz(-1.4962366) q[0];
rz(1.7773588) q[1];
sx q[1];
rz(-0.61449209) q[1];
sx q[1];
rz(-0.38633698) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57418121) q[0];
sx q[0];
rz(-1.8339388) q[0];
sx q[0];
rz(0.69460034) q[0];
rz(-pi) q[1];
rz(-1.3851829) q[2];
sx q[2];
rz(-1.5020554) q[2];
sx q[2];
rz(2.5752714) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5959634) q[1];
sx q[1];
rz(-1.9013604) q[1];
sx q[1];
rz(1.4413553) q[1];
x q[2];
rz(-0.058606996) q[3];
sx q[3];
rz(-0.59133619) q[3];
sx q[3];
rz(-0.014733068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6284457) q[2];
sx q[2];
rz(-1.7690423) q[2];
sx q[2];
rz(2.9885542) q[2];
rz(-0.89093527) q[3];
sx q[3];
rz(-0.95518437) q[3];
sx q[3];
rz(-2.4002767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1560169) q[0];
sx q[0];
rz(-0.22933904) q[0];
sx q[0];
rz(0.34737059) q[0];
rz(0.037847606) q[1];
sx q[1];
rz(-1.1696576) q[1];
sx q[1];
rz(-0.52245021) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.396401) q[0];
sx q[0];
rz(-2.7184882) q[0];
sx q[0];
rz(1.3518831) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1407479) q[2];
sx q[2];
rz(-2.7760996) q[2];
sx q[2];
rz(-1.5216684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47679893) q[1];
sx q[1];
rz(-1.069624) q[1];
sx q[1];
rz(-1.0887926) q[1];
rz(-pi) q[2];
rz(-1.0112004) q[3];
sx q[3];
rz(-1.921145) q[3];
sx q[3];
rz(2.2009785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5106421) q[2];
sx q[2];
rz(-2.6783671) q[2];
sx q[2];
rz(-1.9412712) q[2];
rz(0.55780324) q[3];
sx q[3];
rz(-1.6226945) q[3];
sx q[3];
rz(2.7571078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8484304) q[0];
sx q[0];
rz(-0.98386216) q[0];
sx q[0];
rz(1.4137319) q[0];
rz(0.14104715) q[1];
sx q[1];
rz(-1.3755362) q[1];
sx q[1];
rz(0.65473762) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065718) q[0];
sx q[0];
rz(-1.362934) q[0];
sx q[0];
rz(1.3285358) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59487307) q[2];
sx q[2];
rz(-2.0030177) q[2];
sx q[2];
rz(3.1405666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.480994) q[1];
sx q[1];
rz(-2.8428322) q[1];
sx q[1];
rz(2.1887652) q[1];
x q[2];
rz(0.68924381) q[3];
sx q[3];
rz(-2.3955477) q[3];
sx q[3];
rz(0.41710284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7901788) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(-0.72511017) q[2];
rz(0.890598) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(-2.5543673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17726041) q[0];
sx q[0];
rz(-1.6163419) q[0];
sx q[0];
rz(1.1905715) q[0];
rz(2.019885) q[1];
sx q[1];
rz(-1.5678761) q[1];
sx q[1];
rz(-0.97490464) q[1];
rz(1.2091985) q[2];
sx q[2];
rz(-2.2022916) q[2];
sx q[2];
rz(-1.794338) q[2];
rz(-0.18425758) q[3];
sx q[3];
rz(-2.4202716) q[3];
sx q[3];
rz(0.088868401) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
