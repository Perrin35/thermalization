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
rz(-0.4001652) q[0];
sx q[0];
rz(-0.53181177) q[0];
sx q[0];
rz(2.3398633) q[0];
rz(3.7705295) q[1];
sx q[1];
rz(4.9744199) q[1];
sx q[1];
rz(9.0355102) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2118195) q[0];
sx q[0];
rz(-0.7302098) q[0];
sx q[0];
rz(-0.41359652) q[0];
x q[1];
rz(-1.5250144) q[2];
sx q[2];
rz(-1.3897224) q[2];
sx q[2];
rz(1.2200955) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1133742) q[1];
sx q[1];
rz(-2.1686834) q[1];
sx q[1];
rz(-2.1917927) q[1];
rz(-pi) q[2];
rz(2.2981529) q[3];
sx q[3];
rz(-1.5411845) q[3];
sx q[3];
rz(-0.51472291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4929216) q[2];
sx q[2];
rz(-0.45009437) q[2];
sx q[2];
rz(-2.8802803) q[2];
rz(2.785545) q[3];
sx q[3];
rz(-1.1084403) q[3];
sx q[3];
rz(2.0480919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68591958) q[0];
sx q[0];
rz(-0.57722592) q[0];
sx q[0];
rz(-0.24396954) q[0];
rz(-2.7120554) q[1];
sx q[1];
rz(-1.4725087) q[1];
sx q[1];
rz(0.74091774) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1572621) q[0];
sx q[0];
rz(-1.2411012) q[0];
sx q[0];
rz(1.1801721) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.039224) q[2];
sx q[2];
rz(-0.38313619) q[2];
sx q[2];
rz(-1.3956816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0216995) q[1];
sx q[1];
rz(-2.5346273) q[1];
sx q[1];
rz(-0.40268437) q[1];
x q[2];
rz(-2.4262364) q[3];
sx q[3];
rz(-1.9543867) q[3];
sx q[3];
rz(-2.5959248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.7168768) q[2];
sx q[2];
rz(-1.5224785) q[2];
sx q[2];
rz(-0.87872046) q[2];
rz(0.97484136) q[3];
sx q[3];
rz(-1.8941555) q[3];
sx q[3];
rz(1.461302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9413302) q[0];
sx q[0];
rz(-0.55568475) q[0];
sx q[0];
rz(-0.19361198) q[0];
rz(-0.10784736) q[1];
sx q[1];
rz(-0.71433181) q[1];
sx q[1];
rz(2.329619) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0064751) q[0];
sx q[0];
rz(-2.0251946) q[0];
sx q[0];
rz(-3.019089) q[0];
rz(-pi) q[1];
rz(2.3247955) q[2];
sx q[2];
rz(-1.1650165) q[2];
sx q[2];
rz(-2.7614373) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35235876) q[1];
sx q[1];
rz(-1.339777) q[1];
sx q[1];
rz(-1.1670626) q[1];
x q[2];
rz(-1.1319086) q[3];
sx q[3];
rz(-1.4445856) q[3];
sx q[3];
rz(-2.6708093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0632443) q[2];
sx q[2];
rz(-1.2816659) q[2];
sx q[2];
rz(-3.1380624) q[2];
rz(2.6053536) q[3];
sx q[3];
rz(-0.25501525) q[3];
sx q[3];
rz(2.6030538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2980767) q[0];
sx q[0];
rz(-1.902782) q[0];
sx q[0];
rz(2.3115944) q[0];
rz(-1.7371381) q[1];
sx q[1];
rz(-2.2300215) q[1];
sx q[1];
rz(-2.8703168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41554444) q[0];
sx q[0];
rz(-1.7689287) q[0];
sx q[0];
rz(1.3130441) q[0];
x q[1];
rz(-0.098255946) q[2];
sx q[2];
rz(-0.66731221) q[2];
sx q[2];
rz(1.4405533) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30078706) q[1];
sx q[1];
rz(-1.5983762) q[1];
sx q[1];
rz(-0.3572398) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3448571) q[3];
sx q[3];
rz(-1.4599783) q[3];
sx q[3];
rz(-3.022242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5156877) q[2];
sx q[2];
rz(-2.7582464) q[2];
sx q[2];
rz(2.3670727) q[2];
rz(2.1382616) q[3];
sx q[3];
rz(-0.39350513) q[3];
sx q[3];
rz(2.5956608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.14458732) q[0];
sx q[0];
rz(-2.044675) q[0];
sx q[0];
rz(-1.0928094) q[0];
rz(1.3358759) q[1];
sx q[1];
rz(-2.3922258) q[1];
sx q[1];
rz(-2.5978079) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0896801) q[0];
sx q[0];
rz(-2.8767715) q[0];
sx q[0];
rz(-1.3604683) q[0];
rz(-pi) q[1];
rz(1.6751965) q[2];
sx q[2];
rz(-2.246665) q[2];
sx q[2];
rz(0.75985786) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1256225) q[1];
sx q[1];
rz(-1.984375) q[1];
sx q[1];
rz(0.57012635) q[1];
rz(-pi) q[2];
rz(-0.95044391) q[3];
sx q[3];
rz(-1.5717829) q[3];
sx q[3];
rz(-2.7867449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.233923) q[2];
sx q[2];
rz(-2.449514) q[2];
sx q[2];
rz(-0.27883369) q[2];
rz(-2.8558266) q[3];
sx q[3];
rz(-0.91364342) q[3];
sx q[3];
rz(-2.7260776) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7061507) q[0];
sx q[0];
rz(-2.4984062) q[0];
sx q[0];
rz(-2.1162794) q[0];
rz(-3.0030491) q[1];
sx q[1];
rz(-1.5899315) q[1];
sx q[1];
rz(2.9655546) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4885725) q[0];
sx q[0];
rz(-0.58642846) q[0];
sx q[0];
rz(2.6676548) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7819809) q[2];
sx q[2];
rz(-1.7533025) q[2];
sx q[2];
rz(0.75342853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0149918) q[1];
sx q[1];
rz(-0.53829191) q[1];
sx q[1];
rz(-2.018257) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17060664) q[3];
sx q[3];
rz(-1.1511369) q[3];
sx q[3];
rz(-2.4158286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2328918) q[2];
sx q[2];
rz(-0.42753926) q[2];
sx q[2];
rz(2.3517081) q[2];
rz(-0.32957736) q[3];
sx q[3];
rz(-1.3720082) q[3];
sx q[3];
rz(2.4901539) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71094197) q[0];
sx q[0];
rz(-2.9281404) q[0];
sx q[0];
rz(-2.1479204) q[0];
rz(2.1395394) q[1];
sx q[1];
rz(-2.1498945) q[1];
sx q[1];
rz(1.4659945) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4167673) q[0];
sx q[0];
rz(-0.91499889) q[0];
sx q[0];
rz(1.7955128) q[0];
x q[1];
rz(-1.1118717) q[2];
sx q[2];
rz(-1.6739769) q[2];
sx q[2];
rz(2.9261719) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0198145) q[1];
sx q[1];
rz(-1.7516859) q[1];
sx q[1];
rz(-1.719285) q[1];
rz(2.0577578) q[3];
sx q[3];
rz(-2.3407702) q[3];
sx q[3];
rz(1.6826676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34622908) q[2];
sx q[2];
rz(-2.2987404) q[2];
sx q[2];
rz(-0.25974926) q[2];
rz(0.78884697) q[3];
sx q[3];
rz(-2.2987821) q[3];
sx q[3];
rz(-1.747725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2271093) q[0];
sx q[0];
rz(-1.8552584) q[0];
sx q[0];
rz(2.1504543) q[0];
rz(-1.8046851) q[1];
sx q[1];
rz(-0.82695812) q[1];
sx q[1];
rz(1.7327259) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1240272) q[0];
sx q[0];
rz(-0.47241898) q[0];
sx q[0];
rz(2.8851465) q[0];
rz(-pi) q[1];
rz(1.7193871) q[2];
sx q[2];
rz(-2.598306) q[2];
sx q[2];
rz(-3.0843184) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.046890586) q[1];
sx q[1];
rz(-2.6841503) q[1];
sx q[1];
rz(-2.4777458) q[1];
x q[2];
rz(-2.9314133) q[3];
sx q[3];
rz(-2.4553799) q[3];
sx q[3];
rz(2.8061638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1651429) q[2];
sx q[2];
rz(-1.4850478) q[2];
sx q[2];
rz(0.58839798) q[2];
rz(-2.8277561) q[3];
sx q[3];
rz(-2.2321759) q[3];
sx q[3];
rz(0.61217827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9070076) q[0];
sx q[0];
rz(-1.122965) q[0];
sx q[0];
rz(-0.054397415) q[0];
rz(0.76039487) q[1];
sx q[1];
rz(-0.94579983) q[1];
sx q[1];
rz(1.0992345) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71610036) q[0];
sx q[0];
rz(-1.7370701) q[0];
sx q[0];
rz(0.99207741) q[0];
rz(-pi) q[1];
rz(2.3407158) q[2];
sx q[2];
rz(-0.95716864) q[2];
sx q[2];
rz(2.9517945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1832378) q[1];
sx q[1];
rz(-2.3545165) q[1];
sx q[1];
rz(1.8182204) q[1];
rz(1.4599007) q[3];
sx q[3];
rz(-1.520021) q[3];
sx q[3];
rz(-0.61855701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.86193639) q[2];
sx q[2];
rz(-1.5402237) q[2];
sx q[2];
rz(0.81735617) q[2];
rz(-0.014184626) q[3];
sx q[3];
rz(-0.395702) q[3];
sx q[3];
rz(1.0891677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31187439) q[0];
sx q[0];
rz(-1.051396) q[0];
sx q[0];
rz(-2.9201087) q[0];
rz(0.43137506) q[1];
sx q[1];
rz(-0.90310496) q[1];
sx q[1];
rz(1.0198786) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021602623) q[0];
sx q[0];
rz(-0.051030397) q[0];
sx q[0];
rz(2.7330815) q[0];
x q[1];
rz(0.47616565) q[2];
sx q[2];
rz(-0.22606255) q[2];
sx q[2];
rz(3.0788596) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4012667) q[1];
sx q[1];
rz(-1.5366842) q[1];
sx q[1];
rz(-1.9032065) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76935919) q[3];
sx q[3];
rz(-0.1884547) q[3];
sx q[3];
rz(-2.3615866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43883103) q[2];
sx q[2];
rz(-1.9352545) q[2];
sx q[2];
rz(-1.2552525) q[2];
rz(-0.025764763) q[3];
sx q[3];
rz(-1.8119101) q[3];
sx q[3];
rz(2.6222031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230566) q[0];
sx q[0];
rz(-2.3522455) q[0];
sx q[0];
rz(1.1242207) q[0];
rz(-2.8070246) q[1];
sx q[1];
rz(-0.48929132) q[1];
sx q[1];
rz(-1.0390859) q[1];
rz(0.4998906) q[2];
sx q[2];
rz(-2.6251948) q[2];
sx q[2];
rz(0.37175409) q[2];
rz(2.2065577) q[3];
sx q[3];
rz(-0.72590704) q[3];
sx q[3];
rz(1.4480729) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
