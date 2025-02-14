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
rz(-1.1542198) q[0];
sx q[0];
rz(-1.5232824) q[0];
sx q[0];
rz(-2.9989938) q[0];
rz(0.59347403) q[1];
sx q[1];
rz(-2.4886517) q[1];
sx q[1];
rz(-1.0452193) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0526177) q[0];
sx q[0];
rz(-1.0427022) q[0];
sx q[0];
rz(2.108197) q[0];
rz(-pi) q[1];
rz(-1.6204616) q[2];
sx q[2];
rz(-0.6383183) q[2];
sx q[2];
rz(2.3172003) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3821311) q[1];
sx q[1];
rz(-2.9540017) q[1];
sx q[1];
rz(1.0008903) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7552283) q[3];
sx q[3];
rz(-2.5650555) q[3];
sx q[3];
rz(-0.54631305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97592252) q[2];
sx q[2];
rz(-2.133635) q[2];
sx q[2];
rz(-0.21318501) q[2];
rz(-0.91607696) q[3];
sx q[3];
rz(-0.35917425) q[3];
sx q[3];
rz(-2.7540414) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85584545) q[0];
sx q[0];
rz(-0.88749945) q[0];
sx q[0];
rz(0.47072738) q[0];
rz(-1.1031021) q[1];
sx q[1];
rz(-0.50869894) q[1];
sx q[1];
rz(-0.57977605) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3382028) q[0];
sx q[0];
rz(-2.0943599) q[0];
sx q[0];
rz(-1.5016599) q[0];
rz(-pi) q[1];
rz(-1.2389555) q[2];
sx q[2];
rz(-1.0523426) q[2];
sx q[2];
rz(1.4084852) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3892711) q[1];
sx q[1];
rz(-2.139341) q[1];
sx q[1];
rz(-0.46115498) q[1];
rz(-pi) q[2];
rz(-0.26269368) q[3];
sx q[3];
rz(-1.212271) q[3];
sx q[3];
rz(-2.3101447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9819928) q[2];
sx q[2];
rz(-0.85593587) q[2];
sx q[2];
rz(3.0465872) q[2];
rz(2.5659918) q[3];
sx q[3];
rz(-0.78229457) q[3];
sx q[3];
rz(1.4466205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5461102) q[0];
sx q[0];
rz(-2.4330916) q[0];
sx q[0];
rz(-0.27221671) q[0];
rz(-0.1046003) q[1];
sx q[1];
rz(-2.101254) q[1];
sx q[1];
rz(-1.6786172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7143216) q[0];
sx q[0];
rz(-0.84762379) q[0];
sx q[0];
rz(0.066825213) q[0];
rz(-pi) q[1];
rz(-0.06426364) q[2];
sx q[2];
rz(-0.78635629) q[2];
sx q[2];
rz(2.3051895) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.21468563) q[1];
sx q[1];
rz(-0.72539293) q[1];
sx q[1];
rz(-2.7100485) q[1];
rz(-pi) q[2];
rz(1.5103064) q[3];
sx q[3];
rz(-1.9391427) q[3];
sx q[3];
rz(0.18243901) q[3];
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
rz(-0.05154933) q[3];
sx q[3];
rz(-2.1487273) q[3];
sx q[3];
rz(2.4079017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.16875295) q[0];
sx q[0];
rz(-0.77744716) q[0];
sx q[0];
rz(2.4421316) q[0];
rz(2.2109924) q[1];
sx q[1];
rz(-1.7761296) q[1];
sx q[1];
rz(2.0420989) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65419009) q[0];
sx q[0];
rz(-0.64399566) q[0];
sx q[0];
rz(-1.0985159) q[0];
rz(1.9540953) q[2];
sx q[2];
rz(-2.5591203) q[2];
sx q[2];
rz(-0.7815643) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9605254) q[1];
sx q[1];
rz(-1.2511504) q[1];
sx q[1];
rz(1.8962217) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6190998) q[3];
sx q[3];
rz(-2.429649) q[3];
sx q[3];
rz(-2.5639736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3804271) q[2];
sx q[2];
rz(-1.0720422) q[2];
sx q[2];
rz(-0.94044828) q[2];
rz(-0.56194168) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(-2.565062) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024427323) q[0];
sx q[0];
rz(-0.086240135) q[0];
sx q[0];
rz(-1.0166136) q[0];
rz(-0.92574614) q[1];
sx q[1];
rz(-0.75030202) q[1];
sx q[1];
rz(-0.076676682) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6252687) q[0];
sx q[0];
rz(-1.7364572) q[0];
sx q[0];
rz(1.9693841) q[0];
x q[1];
rz(-1.8163278) q[2];
sx q[2];
rz(-1.7707728) q[2];
sx q[2];
rz(-0.94975805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3073342) q[1];
sx q[1];
rz(-1.3594207) q[1];
sx q[1];
rz(-1.0049051) q[1];
rz(-2.8491945) q[3];
sx q[3];
rz(-2.3030223) q[3];
sx q[3];
rz(2.3866619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3735247) q[2];
sx q[2];
rz(-2.3536451) q[2];
sx q[2];
rz(-3.0184271) q[2];
rz(-0.67433107) q[3];
sx q[3];
rz(-0.47334039) q[3];
sx q[3];
rz(2.3136852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3097565) q[0];
sx q[0];
rz(-0.18853822) q[0];
sx q[0];
rz(-0.48102608) q[0];
rz(2.9006529) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(0.13490881) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8277934) q[0];
sx q[0];
rz(-2.302357) q[0];
sx q[0];
rz(0.93904943) q[0];
rz(-pi) q[1];
rz(-2.606484) q[2];
sx q[2];
rz(-1.6596748) q[2];
sx q[2];
rz(0.8366226) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9489158) q[1];
sx q[1];
rz(-2.2971616) q[1];
sx q[1];
rz(-2.7706258) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9929912) q[3];
sx q[3];
rz(-2.6872748) q[3];
sx q[3];
rz(3.076072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9350819) q[2];
sx q[2];
rz(-0.93588459) q[2];
sx q[2];
rz(1.8180465) q[2];
rz(0.27213085) q[3];
sx q[3];
rz(-0.85796732) q[3];
sx q[3];
rz(-2.9959397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1226591) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(0.65852273) q[0];
rz(-1.5918484) q[1];
sx q[1];
rz(-1.0127944) q[1];
sx q[1];
rz(0.50643593) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0404694) q[0];
sx q[0];
rz(-1.5047025) q[0];
sx q[0];
rz(2.7718146) q[0];
x q[1];
rz(-2.0525706) q[2];
sx q[2];
rz(-1.543569) q[2];
sx q[2];
rz(0.99808642) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87878972) q[1];
sx q[1];
rz(-1.6639941) q[1];
sx q[1];
rz(-1.4944939) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.562256) q[3];
sx q[3];
rz(-1.5361353) q[3];
sx q[3];
rz(-3.0140585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3353614) q[2];
sx q[2];
rz(-2.3956617) q[2];
sx q[2];
rz(1.8719505) q[2];
rz(-1.7757724) q[3];
sx q[3];
rz(-3.1277872) q[3];
sx q[3];
rz(2.5891916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34614554) q[0];
sx q[0];
rz(-2.8988291) q[0];
sx q[0];
rz(-2.7224097) q[0];
rz(3.0508793) q[1];
sx q[1];
rz(-2.2234629) q[1];
sx q[1];
rz(-2.1144287) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8801422) q[0];
sx q[0];
rz(-2.1322726) q[0];
sx q[0];
rz(2.1057157) q[0];
rz(-pi) q[1];
rz(1.5787324) q[2];
sx q[2];
rz(-1.5936226) q[2];
sx q[2];
rz(-2.9812682) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.724714) q[1];
sx q[1];
rz(-1.6771349) q[1];
sx q[1];
rz(1.7843397) q[1];
rz(-2.2563491) q[3];
sx q[3];
rz(-1.5285057) q[3];
sx q[3];
rz(-1.2219747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98296982) q[2];
sx q[2];
rz(-2.136844) q[2];
sx q[2];
rz(0.99009222) q[2];
rz(0.31791911) q[3];
sx q[3];
rz(-2.3293994) q[3];
sx q[3];
rz(3.1032491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38744277) q[0];
sx q[0];
rz(-0.70873547) q[0];
sx q[0];
rz(-0.55737525) q[0];
rz(0.36224657) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(0.16709669) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.635256) q[0];
sx q[0];
rz(-1.8765159) q[0];
sx q[0];
rz(1.291226) q[0];
rz(-0.70126798) q[2];
sx q[2];
rz(-1.7116364) q[2];
sx q[2];
rz(-0.62447157) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22699478) q[1];
sx q[1];
rz(-2.216504) q[1];
sx q[1];
rz(-1.9582932) q[1];
x q[2];
rz(-0.07463712) q[3];
sx q[3];
rz(-0.88866808) q[3];
sx q[3];
rz(-2.9934237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96357137) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(2.7131405) q[2];
rz(0.7306478) q[3];
sx q[3];
rz(-2.0615536) q[3];
sx q[3];
rz(0.60034269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927602) q[0];
sx q[0];
rz(-1.4880143) q[0];
sx q[0];
rz(1.1195419) q[0];
rz(2.4730832) q[1];
sx q[1];
rz(-2.5709277) q[1];
sx q[1];
rz(2.8817435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52212673) q[0];
sx q[0];
rz(-1.3106951) q[0];
sx q[0];
rz(-1.0570265) q[0];
rz(-pi) q[1];
rz(2.2462559) q[2];
sx q[2];
rz(-0.5731715) q[2];
sx q[2];
rz(-1.2473388) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7185134) q[1];
sx q[1];
rz(-1.9155518) q[1];
sx q[1];
rz(-0.43889795) q[1];
rz(-pi) q[2];
rz(-0.15774653) q[3];
sx q[3];
rz(-0.88706568) q[3];
sx q[3];
rz(2.5585548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6707637) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(-2.0551576) q[2];
rz(0.49232617) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(-2.4194748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54871854) q[0];
sx q[0];
rz(-1.508779) q[0];
sx q[0];
rz(-1.4887703) q[0];
rz(-1.8863574) q[1];
sx q[1];
rz(-1.8240758) q[1];
sx q[1];
rz(0.12771894) q[1];
rz(-1.1716539) q[2];
sx q[2];
rz(-1.4831586) q[2];
sx q[2];
rz(-3.0638564) q[2];
rz(-0.90821785) q[3];
sx q[3];
rz(-1.6207909) q[3];
sx q[3];
rz(-0.1089628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
