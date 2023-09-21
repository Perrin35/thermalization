OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(7.1516501) q[0];
sx q[0];
rz(9.2317543) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2642759) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(0.067406128) q[0];
x q[1];
rz(-2.1630152) q[2];
sx q[2];
rz(-2.7263612) q[2];
sx q[2];
rz(2.4821971) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5176757) q[1];
sx q[1];
rz(-0.78849925) q[1];
sx q[1];
rz(-2.4967525) q[1];
rz(2.7636823) q[3];
sx q[3];
rz(-1.0381191) q[3];
sx q[3];
rz(-2.7045254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1258939) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(1.0788318) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9803479) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(0.20794491) q[0];
rz(2.5646599) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(1.6764486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0365021) q[0];
sx q[0];
rz(-2.4173792) q[0];
sx q[0];
rz(3.1378531) q[0];
rz(-pi) q[1];
rz(-2.0287839) q[2];
sx q[2];
rz(-2.6070242) q[2];
sx q[2];
rz(2.2528258) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1190471) q[1];
sx q[1];
rz(-2.3733768) q[1];
sx q[1];
rz(-0.70336282) q[1];
x q[2];
rz(1.2259237) q[3];
sx q[3];
rz(-0.85540918) q[3];
sx q[3];
rz(-1.2777002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1229822) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(-0.1097651) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(-2.6254568) q[0];
rz(2.5667045) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(-2.3410472) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.752906) q[0];
sx q[0];
rz(-1.1645002) q[0];
sx q[0];
rz(-0.15588649) q[0];
x q[1];
rz(-2.3163047) q[2];
sx q[2];
rz(-2.2824259) q[2];
sx q[2];
rz(0.7691783) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2785981) q[1];
sx q[1];
rz(-1.1578214) q[1];
sx q[1];
rz(1.5763271) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.839746) q[3];
sx q[3];
rz(-0.89248025) q[3];
sx q[3];
rz(2.9523926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.67733726) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.4250071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99252218) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(2.676679) q[0];
rz(-0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(2.0565313) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5269055) q[0];
sx q[0];
rz(-2.5589716) q[0];
sx q[0];
rz(-0.42838642) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2041353) q[2];
sx q[2];
rz(-1.8030093) q[2];
sx q[2];
rz(2.1249078) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3738721) q[1];
sx q[1];
rz(-1.4704736) q[1];
sx q[1];
rz(-2.9994681) q[1];
rz(-0.61120175) q[3];
sx q[3];
rz(-1.9107606) q[3];
sx q[3];
rz(1.9577648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(-0.68112779) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(-0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(-2.1599105) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22589332) q[0];
sx q[0];
rz(-0.75076538) q[0];
sx q[0];
rz(-2.4322926) q[0];
x q[1];
rz(0.9972516) q[2];
sx q[2];
rz(-2.855636) q[2];
sx q[2];
rz(2.2821102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9202068) q[1];
sx q[1];
rz(-2.033794) q[1];
sx q[1];
rz(1.4057926) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97335191) q[3];
sx q[3];
rz(-1.6320328) q[3];
sx q[3];
rz(2.6254568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(-0.48941082) q[2];
rz(-2.126746) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(-1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(0.16528027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344791) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(-1.6001742) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4539102) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(0.24535594) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89989923) q[1];
sx q[1];
rz(-0.69679931) q[1];
sx q[1];
rz(0.65710575) q[1];
x q[2];
rz(1.1214439) q[3];
sx q[3];
rz(-1.3580139) q[3];
sx q[3];
rz(-1.3495812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(1.4432663) q[2];
rz(-0.36744395) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7845602) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(2.7923287) q[0];
rz(-0.7473942) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-2.4051037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38024494) q[0];
sx q[0];
rz(-1.5878829) q[0];
sx q[0];
rz(1.6197617) q[0];
rz(-0.24918208) q[2];
sx q[2];
rz(-1.2417925) q[2];
sx q[2];
rz(0.56275425) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4941102) q[1];
sx q[1];
rz(-0.78299114) q[1];
sx q[1];
rz(-0.25581911) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69865366) q[3];
sx q[3];
rz(-2.5210288) q[3];
sx q[3];
rz(2.8750949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(0.53156701) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(-1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
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
rz(3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(2.334306) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(0.92179006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1944273) q[0];
sx q[0];
rz(-2.2568586) q[0];
sx q[0];
rz(-1.7454733) q[0];
rz(2.3085262) q[2];
sx q[2];
rz(-1.3563915) q[2];
sx q[2];
rz(-3.0467141) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.31729749) q[1];
sx q[1];
rz(-1.3166787) q[1];
sx q[1];
rz(-1.7829249) q[1];
rz(-pi) q[2];
rz(1.7796302) q[3];
sx q[3];
rz(-2.4224671) q[3];
sx q[3];
rz(1.5987087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(-1.3575859) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(-1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8885324) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(-1.4755479) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(-2.4618861) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4363791) q[0];
sx q[0];
rz(-2.0853015) q[0];
sx q[0];
rz(-1.8041457) q[0];
rz(1.7858511) q[2];
sx q[2];
rz(-0.70677033) q[2];
sx q[2];
rz(3.0422473) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6913773) q[1];
sx q[1];
rz(-0.9141578) q[1];
sx q[1];
rz(-2.0037829) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4623975) q[3];
sx q[3];
rz(-1.067357) q[3];
sx q[3];
rz(0.60590832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(-2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0891721) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(2.4269379) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(1.9649327) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0929716) q[0];
sx q[0];
rz(-1.4074568) q[0];
sx q[0];
rz(-3.0924348) q[0];
rz(-pi) q[1];
rz(2.2628225) q[2];
sx q[2];
rz(-0.76239097) q[2];
sx q[2];
rz(-0.15904418) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1706108) q[1];
sx q[1];
rz(-0.5628399) q[1];
sx q[1];
rz(-1.7112205) q[1];
rz(2.9858399) q[3];
sx q[3];
rz(-0.51625801) q[3];
sx q[3];
rz(1.932365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(-1.127355) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42416278) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-1.2926119) q[2];
sx q[2];
rz(-0.85360151) q[2];
sx q[2];
rz(-3.1384946) q[2];
rz(-3.1191961) q[3];
sx q[3];
rz(-0.35084421) q[3];
sx q[3];
rz(2.4435333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];