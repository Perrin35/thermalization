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
rz(0.81646252) q[0];
sx q[0];
rz(3.2433885) q[0];
sx q[0];
rz(9.9643702) q[0];
rz(0.49630961) q[1];
sx q[1];
rz(-0.30975431) q[1];
sx q[1];
rz(0.53909477) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2236299) q[0];
sx q[0];
rz(-1.5453891) q[0];
sx q[0];
rz(0.44141234) q[0];
rz(-pi) q[1];
rz(0.38549785) q[2];
sx q[2];
rz(-2.756167) q[2];
sx q[2];
rz(2.5138234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2869563) q[1];
sx q[1];
rz(-2.7883734) q[1];
sx q[1];
rz(-1.1562721) q[1];
rz(2.21978) q[3];
sx q[3];
rz(-1.505449) q[3];
sx q[3];
rz(1.9416888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3794136) q[2];
sx q[2];
rz(-2.2654686) q[2];
sx q[2];
rz(-1.1890821) q[2];
rz(1.9573697) q[3];
sx q[3];
rz(-0.90457478) q[3];
sx q[3];
rz(-1.5949465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14520833) q[0];
sx q[0];
rz(-1.5897911) q[0];
sx q[0];
rz(1.8183964) q[0];
rz(-2.6595751) q[1];
sx q[1];
rz(-2.2299485) q[1];
sx q[1];
rz(0.97420305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7071814) q[0];
sx q[0];
rz(-1.392258) q[0];
sx q[0];
rz(-2.2212127) q[0];
rz(-pi) q[1];
rz(-2.2902238) q[2];
sx q[2];
rz(-0.15002827) q[2];
sx q[2];
rz(-1.139037) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99813733) q[1];
sx q[1];
rz(-1.0364729) q[1];
sx q[1];
rz(-1.8443395) q[1];
rz(-0.5663381) q[3];
sx q[3];
rz(-1.7741331) q[3];
sx q[3];
rz(1.7254064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1366068) q[2];
sx q[2];
rz(-2.5544781) q[2];
sx q[2];
rz(-0.74964398) q[2];
rz(2.7096115) q[3];
sx q[3];
rz(-1.1155198) q[3];
sx q[3];
rz(1.3031134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62464803) q[0];
sx q[0];
rz(-0.2722781) q[0];
sx q[0];
rz(2.5352449) q[0];
rz(-1.3336522) q[1];
sx q[1];
rz(-1.2140112) q[1];
sx q[1];
rz(-2.8820754) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1071651) q[0];
sx q[0];
rz(-1.3239064) q[0];
sx q[0];
rz(-0.16040032) q[0];
x q[1];
rz(-2.9664842) q[2];
sx q[2];
rz(-1.3710183) q[2];
sx q[2];
rz(1.4005043) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3644581) q[1];
sx q[1];
rz(-1.196047) q[1];
sx q[1];
rz(-2.0261637) q[1];
rz(0.43478888) q[3];
sx q[3];
rz(-1.0514976) q[3];
sx q[3];
rz(3.0360707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8848662) q[2];
sx q[2];
rz(-1.3639516) q[2];
sx q[2];
rz(-1.5941031) q[2];
rz(-0.78440845) q[3];
sx q[3];
rz(-1.6268077) q[3];
sx q[3];
rz(-1.5659531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49482685) q[0];
sx q[0];
rz(-1.0972728) q[0];
sx q[0];
rz(-2.8821017) q[0];
rz(-1.9208113) q[1];
sx q[1];
rz(-0.44993284) q[1];
sx q[1];
rz(-1.4422013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3078559) q[0];
sx q[0];
rz(-2.4015732) q[0];
sx q[0];
rz(0.83777512) q[0];
rz(-0.97564189) q[2];
sx q[2];
rz(-1.6791428) q[2];
sx q[2];
rz(2.8841022) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.759023) q[1];
sx q[1];
rz(-1.3241395) q[1];
sx q[1];
rz(0.36212977) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4846701) q[3];
sx q[3];
rz(-0.68400506) q[3];
sx q[3];
rz(-2.4620584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1059025) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(-2.7316459) q[2];
rz(-1.9299054) q[3];
sx q[3];
rz(-2.4544921) q[3];
sx q[3];
rz(1.80779) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164417) q[0];
sx q[0];
rz(-0.71617675) q[0];
sx q[0];
rz(0.27405611) q[0];
rz(0.43824276) q[1];
sx q[1];
rz(-1.2934877) q[1];
sx q[1];
rz(2.4058707) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11392191) q[0];
sx q[0];
rz(-0.7589853) q[0];
sx q[0];
rz(3.0149197) q[0];
x q[1];
rz(-2.9922036) q[2];
sx q[2];
rz(-1.427703) q[2];
sx q[2];
rz(-0.18994513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.84668881) q[1];
sx q[1];
rz(-1.3219993) q[1];
sx q[1];
rz(1.2685246) q[1];
x q[2];
rz(-0.51208074) q[3];
sx q[3];
rz(-0.98618531) q[3];
sx q[3];
rz(2.0090112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.938544) q[2];
sx q[2];
rz(-1.4578578) q[2];
sx q[2];
rz(2.5277444) q[2];
rz(2.7326873) q[3];
sx q[3];
rz(-0.96858612) q[3];
sx q[3];
rz(2.3992505) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3530389) q[0];
sx q[0];
rz(-0.63816324) q[0];
sx q[0];
rz(-1.2370538) q[0];
rz(-1.4452112) q[1];
sx q[1];
rz(-0.45672363) q[1];
sx q[1];
rz(2.9663185) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0171368) q[0];
sx q[0];
rz(-1.5457898) q[0];
sx q[0];
rz(-1.3978005) q[0];
rz(-pi) q[1];
rz(-2.4779123) q[2];
sx q[2];
rz(-0.77518565) q[2];
sx q[2];
rz(2.7221687) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3289017) q[1];
sx q[1];
rz(-2.0022503) q[1];
sx q[1];
rz(0.56568362) q[1];
x q[2];
rz(-2.3168269) q[3];
sx q[3];
rz(-0.76465339) q[3];
sx q[3];
rz(0.95461707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.96482977) q[2];
sx q[2];
rz(-1.333933) q[2];
sx q[2];
rz(0.97243398) q[2];
rz(-1.4771627) q[3];
sx q[3];
rz(-1.9921314) q[3];
sx q[3];
rz(-0.76178637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2870188) q[0];
sx q[0];
rz(-3.119097) q[0];
sx q[0];
rz(2.7884685) q[0];
rz(-1.0278206) q[1];
sx q[1];
rz(-1.2007583) q[1];
sx q[1];
rz(-1.0110528) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7582304) q[0];
sx q[0];
rz(-2.7840021) q[0];
sx q[0];
rz(0.041186853) q[0];
x q[1];
rz(-2.2765144) q[2];
sx q[2];
rz(-2.0524244) q[2];
sx q[2];
rz(-2.8291246) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8450635) q[1];
sx q[1];
rz(-1.4140714) q[1];
sx q[1];
rz(1.2178376) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8500438) q[3];
sx q[3];
rz(-1.8729775) q[3];
sx q[3];
rz(-2.9460689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.31356835) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(1.774452) q[2];
rz(1.2683055) q[3];
sx q[3];
rz(-1.6627848) q[3];
sx q[3];
rz(-0.84793276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0276412) q[0];
sx q[0];
rz(-2.5340762) q[0];
sx q[0];
rz(-2.0116346) q[0];
rz(-2.9226411) q[1];
sx q[1];
rz(-1.4066701) q[1];
sx q[1];
rz(2.2311282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51800358) q[0];
sx q[0];
rz(-1.2388889) q[0];
sx q[0];
rz(0.62923543) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6633419) q[2];
sx q[2];
rz(-1.4946117) q[2];
sx q[2];
rz(0.26631276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76669756) q[1];
sx q[1];
rz(-1.3037221) q[1];
sx q[1];
rz(1.0407991) q[1];
rz(-pi) q[2];
rz(-1.1651785) q[3];
sx q[3];
rz(-2.2542103) q[3];
sx q[3];
rz(-2.3644517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8255446) q[2];
sx q[2];
rz(-0.80222183) q[2];
sx q[2];
rz(0.82553378) q[2];
rz(-0.46142203) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(-0.44574827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-2.5017186) q[0];
sx q[0];
rz(-1.3383144) q[0];
sx q[0];
rz(-2.3097532) q[0];
rz(-0.72744751) q[1];
sx q[1];
rz(-1.7214382) q[1];
sx q[1];
rz(0.096253455) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1548658) q[0];
sx q[0];
rz(-2.4270202) q[0];
sx q[0];
rz(-2.3113219) q[0];
x q[1];
rz(-0.23230884) q[2];
sx q[2];
rz(-0.26089982) q[2];
sx q[2];
rz(-1.0123073) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2029496) q[1];
sx q[1];
rz(-1.661282) q[1];
sx q[1];
rz(-1.7825148) q[1];
x q[2];
rz(1.3745802) q[3];
sx q[3];
rz(-1.726578) q[3];
sx q[3];
rz(2.2957612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7353797) q[2];
sx q[2];
rz(-1.7071743) q[2];
sx q[2];
rz(-1.5926788) q[2];
rz(2.5088572) q[3];
sx q[3];
rz(-1.2318719) q[3];
sx q[3];
rz(-2.8922141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.7525472) q[0];
sx q[0];
rz(-2.1010375) q[0];
sx q[0];
rz(2.7885875) q[0];
rz(-2.8189335) q[1];
sx q[1];
rz(-2.8162075) q[1];
sx q[1];
rz(0.37193146) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7854561) q[0];
sx q[0];
rz(-1.2704986) q[0];
sx q[0];
rz(-0.2072643) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3916799) q[2];
sx q[2];
rz(-1.2460016) q[2];
sx q[2];
rz(1.3472404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5185753) q[1];
sx q[1];
rz(-1.2378344) q[1];
sx q[1];
rz(-1.740231) q[1];
rz(-pi) q[2];
rz(0.84623611) q[3];
sx q[3];
rz(-1.5414943) q[3];
sx q[3];
rz(1.3361564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9670664) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(-1.8444427) q[2];
rz(-0.8477115) q[3];
sx q[3];
rz(-1.1573557) q[3];
sx q[3];
rz(-2.1367836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2321155) q[0];
sx q[0];
rz(-1.7625325) q[0];
sx q[0];
rz(0.43176227) q[0];
rz(1.7535946) q[1];
sx q[1];
rz(-1.2139865) q[1];
sx q[1];
rz(-0.075275631) q[1];
rz(-0.70099945) q[2];
sx q[2];
rz(-0.95078118) q[2];
sx q[2];
rz(-2.4581428) q[2];
rz(1.9907436) q[3];
sx q[3];
rz(-0.81880488) q[3];
sx q[3];
rz(2.8863751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
