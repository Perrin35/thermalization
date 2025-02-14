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
rz(0.59347403) q[1];
sx q[1];
rz(3.7945336) q[1];
sx q[1];
rz(8.3795587) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3616989) q[0];
sx q[0];
rz(-2.4068659) q[0];
sx q[0];
rz(2.421342) q[0];
x q[1];
rz(0.93306834) q[2];
sx q[2];
rz(-1.5412113) q[2];
sx q[2];
rz(-2.3552908) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3821311) q[1];
sx q[1];
rz(-2.9540017) q[1];
sx q[1];
rz(2.1407024) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7552283) q[3];
sx q[3];
rz(-0.57653713) q[3];
sx q[3];
rz(2.5952796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97592252) q[2];
sx q[2];
rz(-1.0079577) q[2];
sx q[2];
rz(0.21318501) q[2];
rz(-0.91607696) q[3];
sx q[3];
rz(-2.7824184) q[3];
sx q[3];
rz(-0.38755125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85584545) q[0];
sx q[0];
rz(-0.88749945) q[0];
sx q[0];
rz(2.6708653) q[0];
rz(1.1031021) q[1];
sx q[1];
rz(-0.50869894) q[1];
sx q[1];
rz(0.57977605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.267201) q[0];
sx q[0];
rz(-1.6306595) q[0];
sx q[0];
rz(2.6169928) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2389555) q[2];
sx q[2];
rz(-1.0523426) q[2];
sx q[2];
rz(-1.4084852) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3892711) q[1];
sx q[1];
rz(-2.139341) q[1];
sx q[1];
rz(-2.6804377) q[1];
x q[2];
rz(-1.9409431) q[3];
sx q[3];
rz(-1.3251628) q[3];
sx q[3];
rz(0.64526886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9819928) q[2];
sx q[2];
rz(-0.85593587) q[2];
sx q[2];
rz(-3.0465872) q[2];
rz(0.57560086) q[3];
sx q[3];
rz(-0.78229457) q[3];
sx q[3];
rz(-1.4466205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59548241) q[0];
sx q[0];
rz(-2.4330916) q[0];
sx q[0];
rz(-2.8693759) q[0];
rz(0.1046003) q[1];
sx q[1];
rz(-1.0403386) q[1];
sx q[1];
rz(-1.6786172) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8151096) q[0];
sx q[0];
rz(-0.72569752) q[0];
sx q[0];
rz(1.6463019) q[0];
rz(-pi) q[1];
rz(3.077329) q[2];
sx q[2];
rz(-2.3552364) q[2];
sx q[2];
rz(0.83640316) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21468563) q[1];
sx q[1];
rz(-0.72539293) q[1];
sx q[1];
rz(2.7100485) q[1];
rz(-pi) q[2];
rz(-2.7726309) q[3];
sx q[3];
rz(-1.6272244) q[3];
sx q[3];
rz(-1.4101613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0766803) q[2];
sx q[2];
rz(-1.7408337) q[2];
sx q[2];
rz(-0.40985516) q[2];
rz(-0.05154933) q[3];
sx q[3];
rz(-2.1487273) q[3];
sx q[3];
rz(2.4079017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.16875295) q[0];
sx q[0];
rz(-0.77744716) q[0];
sx q[0];
rz(-2.4421316) q[0];
rz(-2.2109924) q[1];
sx q[1];
rz(-1.7761296) q[1];
sx q[1];
rz(-2.0420989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52877766) q[0];
sx q[0];
rz(-1.2941501) q[0];
sx q[0];
rz(-0.98146455) q[0];
rz(-0.24153562) q[2];
sx q[2];
rz(-1.0354098) q[2];
sx q[2];
rz(-2.8098742) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9605254) q[1];
sx q[1];
rz(-1.2511504) q[1];
sx q[1];
rz(1.8962217) q[1];
rz(-pi) q[2];
rz(2.4995125) q[3];
sx q[3];
rz(-1.9028946) q[3];
sx q[3];
rz(1.737271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3804271) q[2];
sx q[2];
rz(-1.0720422) q[2];
sx q[2];
rz(2.2011444) q[2];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
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
rz(2.2158465) q[1];
sx q[1];
rz(-0.75030202) q[1];
sx q[1];
rz(-0.076676682) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0177855) q[0];
sx q[0];
rz(-1.9636256) q[0];
sx q[0];
rz(-0.17946243) q[0];
x q[1];
rz(-0.87575298) q[2];
sx q[2];
rz(-0.31538439) q[2];
sx q[2];
rz(-1.8502667) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5591084) q[1];
sx q[1];
rz(-2.541572) q[1];
sx q[1];
rz(-1.9514854) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2604976) q[3];
sx q[3];
rz(-0.7782794) q[3];
sx q[3];
rz(2.8097514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3735247) q[2];
sx q[2];
rz(-2.3536451) q[2];
sx q[2];
rz(3.0184271) q[2];
rz(2.4672616) q[3];
sx q[3];
rz(-0.47334039) q[3];
sx q[3];
rz(-0.82790747) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3097565) q[0];
sx q[0];
rz(-2.9530544) q[0];
sx q[0];
rz(2.6605666) q[0];
rz(2.9006529) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(-3.0066838) q[1];
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
rz(2.606484) q[2];
sx q[2];
rz(-1.6596748) q[2];
sx q[2];
rz(-0.8366226) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9489158) q[1];
sx q[1];
rz(-2.2971616) q[1];
sx q[1];
rz(2.7706258) q[1];
rz(-pi) q[2];
rz(1.4986131) q[3];
sx q[3];
rz(-2.0197387) q[3];
sx q[3];
rz(0.23060966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9350819) q[2];
sx q[2];
rz(-2.2057081) q[2];
sx q[2];
rz(-1.3235462) q[2];
rz(-0.27213085) q[3];
sx q[3];
rz(-0.85796732) q[3];
sx q[3];
rz(-0.14565295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1226591) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(2.4830699) q[0];
rz(1.5918484) q[1];
sx q[1];
rz(-2.1287983) q[1];
sx q[1];
rz(0.50643593) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50473266) q[0];
sx q[0];
rz(-1.2018645) q[0];
sx q[0];
rz(1.4999267) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1108702) q[2];
sx q[2];
rz(-2.0523768) q[2];
sx q[2];
rz(0.55847634) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5665019) q[1];
sx q[1];
rz(-3.021214) q[1];
sx q[1];
rz(0.68415227) q[1];
x q[2];
rz(1.6122088) q[3];
sx q[3];
rz(-0.99185252) q[3];
sx q[3];
rz(1.4659297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8062313) q[2];
sx q[2];
rz(-0.74593097) q[2];
sx q[2];
rz(1.2696421) q[2];
rz(-1.3658203) q[3];
sx q[3];
rz(-3.1277872) q[3];
sx q[3];
rz(-2.5891916) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34614554) q[0];
sx q[0];
rz(-2.8988291) q[0];
sx q[0];
rz(2.7224097) q[0];
rz(0.090713352) q[1];
sx q[1];
rz(-2.2234629) q[1];
sx q[1];
rz(-1.027164) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8801422) q[0];
sx q[0];
rz(-1.0093201) q[0];
sx q[0];
rz(-2.1057157) q[0];
rz(-1.5628603) q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60909437) q[1];
sx q[1];
rz(-2.9033992) q[1];
sx q[1];
rz(1.1042117) q[1];
x q[2];
rz(1.5040565) q[3];
sx q[3];
rz(-0.68664521) q[3];
sx q[3];
rz(0.29717577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98296982) q[2];
sx q[2];
rz(-1.0047487) q[2];
sx q[2];
rz(0.99009222) q[2];
rz(-2.8236735) q[3];
sx q[3];
rz(-2.3293994) q[3];
sx q[3];
rz(-0.038343553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7541499) q[0];
sx q[0];
rz(-2.4328572) q[0];
sx q[0];
rz(-0.55737525) q[0];
rz(0.36224657) q[1];
sx q[1];
rz(-1.4366415) q[1];
sx q[1];
rz(-0.16709669) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97826577) q[0];
sx q[0];
rz(-1.8370795) q[0];
sx q[0];
rz(-0.31727088) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4403247) q[2];
sx q[2];
rz(-1.7116364) q[2];
sx q[2];
rz(2.5171211) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5846338) q[1];
sx q[1];
rz(-1.8773729) q[1];
sx q[1];
rz(-2.4584998) q[1];
x q[2];
rz(3.0669555) q[3];
sx q[3];
rz(-2.2529246) q[3];
sx q[3];
rz(2.9934237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.96357137) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(-2.7131405) q[2];
rz(-0.7306478) q[3];
sx q[3];
rz(-2.0615536) q[3];
sx q[3];
rz(-0.60034269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.44883248) q[0];
sx q[0];
rz(-1.6535783) q[0];
sx q[0];
rz(2.0220508) q[0];
rz(-2.4730832) q[1];
sx q[1];
rz(-2.5709277) q[1];
sx q[1];
rz(-2.8817435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4762806) q[0];
sx q[0];
rz(-2.5710433) q[0];
sx q[0];
rz(-2.0670939) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.75801) q[2];
sx q[2];
rz(-2.0077939) q[2];
sx q[2];
rz(-2.008977) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4230792) q[1];
sx q[1];
rz(-1.2260409) q[1];
sx q[1];
rz(-2.7026947) q[1];
x q[2];
rz(-2.2606528) q[3];
sx q[3];
rz(-1.6928813) q[3];
sx q[3];
rz(1.087904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6707637) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(2.0551576) q[2];
rz(-2.6492665) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(0.72211784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54871854) q[0];
sx q[0];
rz(-1.508779) q[0];
sx q[0];
rz(-1.4887703) q[0];
rz(1.2552352) q[1];
sx q[1];
rz(-1.8240758) q[1];
sx q[1];
rz(0.12771894) q[1];
rz(-1.7931425) q[2];
sx q[2];
rz(-2.7334474) q[2];
sx q[2];
rz(1.853142) q[2];
rz(-2.2333748) q[3];
sx q[3];
rz(-1.5208018) q[3];
sx q[3];
rz(3.0326299) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
