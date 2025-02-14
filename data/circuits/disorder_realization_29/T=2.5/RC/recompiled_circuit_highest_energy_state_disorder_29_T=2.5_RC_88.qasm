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
rz(-1.8972549) q[0];
sx q[0];
rz(-2.3490348) q[0];
sx q[0];
rz(-2.8886524) q[0];
rz(-0.77415544) q[1];
sx q[1];
rz(-0.3781265) q[1];
sx q[1];
rz(0.3869431) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16432504) q[0];
sx q[0];
rz(-1.8935793) q[0];
sx q[0];
rz(2.8708145) q[0];
rz(-0.57780452) q[2];
sx q[2];
rz(-1.6434323) q[2];
sx q[2];
rz(-0.60608038) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7751037) q[1];
sx q[1];
rz(-1.6548043) q[1];
sx q[1];
rz(1.6506399) q[1];
rz(-pi) q[2];
rz(-0.67220848) q[3];
sx q[3];
rz(-0.70939964) q[3];
sx q[3];
rz(-1.0649452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0953377) q[2];
sx q[2];
rz(-2.3825808) q[2];
sx q[2];
rz(1.7389899) q[2];
rz(-1.4143573) q[3];
sx q[3];
rz(-2.4284913) q[3];
sx q[3];
rz(2.6426017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4614748) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(0.71037355) q[0];
rz(-2.12517) q[1];
sx q[1];
rz(-1.8417532) q[1];
sx q[1];
rz(2.3764835) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0046526) q[0];
sx q[0];
rz(-1.0400794) q[0];
sx q[0];
rz(-1.9416481) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4508453) q[2];
sx q[2];
rz(-2.1155069) q[2];
sx q[2];
rz(-2.4360457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5580235) q[1];
sx q[1];
rz(-0.92881948) q[1];
sx q[1];
rz(2.7318673) q[1];
x q[2];
rz(0.14987544) q[3];
sx q[3];
rz(-1.2399763) q[3];
sx q[3];
rz(-0.56043032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0244828) q[2];
sx q[2];
rz(-0.64579248) q[2];
sx q[2];
rz(-2.4369241) q[2];
rz(-0.17413983) q[3];
sx q[3];
rz(-0.87248674) q[3];
sx q[3];
rz(-0.38159889) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0524549) q[0];
sx q[0];
rz(-1.8280886) q[0];
sx q[0];
rz(-2.254159) q[0];
rz(1.8916091) q[1];
sx q[1];
rz(-1.2108112) q[1];
sx q[1];
rz(1.7807622) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0075390752) q[0];
sx q[0];
rz(-2.0749395) q[0];
sx q[0];
rz(2.0043122) q[0];
x q[1];
rz(1.4329696) q[2];
sx q[2];
rz(-1.4915474) q[2];
sx q[2];
rz(2.2887231) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4579297) q[1];
sx q[1];
rz(-2.3291846) q[1];
sx q[1];
rz(0.13813604) q[1];
x q[2];
rz(1.674747) q[3];
sx q[3];
rz(-0.39109215) q[3];
sx q[3];
rz(1.1268738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36797324) q[2];
sx q[2];
rz(-2.7757288) q[2];
sx q[2];
rz(-1.492929) q[2];
rz(2.6922373) q[3];
sx q[3];
rz(-1.8466693) q[3];
sx q[3];
rz(1.2202643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.087695) q[0];
sx q[0];
rz(-2.5515285) q[0];
sx q[0];
rz(1.7328523) q[0];
rz(0.98211163) q[1];
sx q[1];
rz(-1.5487919) q[1];
sx q[1];
rz(1.7353479) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2236693) q[0];
sx q[0];
rz(-1.6817998) q[0];
sx q[0];
rz(1.7864321) q[0];
rz(-pi) q[1];
rz(-2.2135229) q[2];
sx q[2];
rz(-0.065970369) q[2];
sx q[2];
rz(2.3118491) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0413734) q[1];
sx q[1];
rz(-2.5587359) q[1];
sx q[1];
rz(1.1122056) q[1];
rz(2.8729183) q[3];
sx q[3];
rz(-0.93963913) q[3];
sx q[3];
rz(0.97418601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9881607) q[2];
sx q[2];
rz(-2.1250696) q[2];
sx q[2];
rz(2.9856258) q[2];
rz(2.4912452) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(-0.11421886) q[3];
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
rz(-pi/2) q[3];
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
rz(-2.0680189) q[0];
sx q[0];
rz(-1.8719712) q[0];
sx q[0];
rz(-0.20075783) q[0];
rz(1.9643895) q[1];
sx q[1];
rz(-2.2889844) q[1];
sx q[1];
rz(1.9815365) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33452144) q[0];
sx q[0];
rz(-1.40309) q[0];
sx q[0];
rz(-1.3434065) q[0];
x q[1];
rz(-3.026408) q[2];
sx q[2];
rz(-1.6403926) q[2];
sx q[2];
rz(-1.7382966) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7457218) q[1];
sx q[1];
rz(-1.9211287) q[1];
sx q[1];
rz(-1.2474023) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1006159) q[3];
sx q[3];
rz(-1.3966832) q[3];
sx q[3];
rz(-1.8903738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7572299) q[2];
sx q[2];
rz(-2.195916) q[2];
sx q[2];
rz(0.45822701) q[2];
rz(0.630817) q[3];
sx q[3];
rz(-1.9061371) q[3];
sx q[3];
rz(0.49921504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-0.45796564) q[0];
sx q[0];
rz(-0.85245913) q[0];
sx q[0];
rz(3.1347347) q[0];
rz(2.9099756) q[1];
sx q[1];
rz(-1.7624785) q[1];
sx q[1];
rz(1.0622271) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10929617) q[0];
sx q[0];
rz(-1.8826906) q[0];
sx q[0];
rz(1.4410785) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7825323) q[2];
sx q[2];
rz(-1.0811624) q[2];
sx q[2];
rz(-0.28070606) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5353229) q[1];
sx q[1];
rz(-0.71782604) q[1];
sx q[1];
rz(-2.7503783) q[1];
x q[2];
rz(-1.3788578) q[3];
sx q[3];
rz(-0.71832359) q[3];
sx q[3];
rz(-3.0305221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0007533) q[2];
sx q[2];
rz(-1.8893628) q[2];
sx q[2];
rz(1.8737277) q[2];
rz(-0.86152348) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(-1.7027732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.93290257) q[0];
sx q[0];
rz(-2.8080495) q[0];
sx q[0];
rz(-0.21555899) q[0];
rz(-2.6453099) q[1];
sx q[1];
rz(-1.1880621) q[1];
sx q[1];
rz(-0.15942474) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2941344) q[0];
sx q[0];
rz(-2.205664) q[0];
sx q[0];
rz(-0.0032739689) q[0];
x q[1];
rz(-0.47093289) q[2];
sx q[2];
rz(-1.7438981) q[2];
sx q[2];
rz(1.7488232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7466399) q[1];
sx q[1];
rz(-1.3356908) q[1];
sx q[1];
rz(3.0993153) q[1];
rz(1.0219021) q[3];
sx q[3];
rz(-1.1682142) q[3];
sx q[3];
rz(0.91820133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7364007) q[2];
sx q[2];
rz(-2.0173732) q[2];
sx q[2];
rz(-0.51652017) q[2];
rz(1.2107595) q[3];
sx q[3];
rz(-1.4529994) q[3];
sx q[3];
rz(1.3794544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4844168) q[0];
sx q[0];
rz(-3.0653937) q[0];
sx q[0];
rz(-0.54022378) q[0];
rz(-1.6172488) q[1];
sx q[1];
rz(-1.0018307) q[1];
sx q[1];
rz(1.2670955) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62477797) q[0];
sx q[0];
rz(-2.5886726) q[0];
sx q[0];
rz(0.42527683) q[0];
rz(0.77719633) q[2];
sx q[2];
rz(-1.0269477) q[2];
sx q[2];
rz(1.0842619) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9766751) q[1];
sx q[1];
rz(-1.7150214) q[1];
sx q[1];
rz(0.37101908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2121939) q[3];
sx q[3];
rz(-1.2817973) q[3];
sx q[3];
rz(2.3232587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2544864) q[2];
sx q[2];
rz(-1.9820513) q[2];
sx q[2];
rz(0.17733388) q[2];
rz(1.0698498) q[3];
sx q[3];
rz(-1.0515352) q[3];
sx q[3];
rz(-2.9593318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808481) q[0];
sx q[0];
rz(-1.9875263) q[0];
sx q[0];
rz(3.0408707) q[0];
rz(1.992647) q[1];
sx q[1];
rz(-1.1958242) q[1];
sx q[1];
rz(1.1740059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7851104) q[0];
sx q[0];
rz(-2.504341) q[0];
sx q[0];
rz(-2.8204945) q[0];
rz(-3.0413925) q[2];
sx q[2];
rz(-0.76876193) q[2];
sx q[2];
rz(0.59610808) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.18408882) q[1];
sx q[1];
rz(-0.042467707) q[1];
sx q[1];
rz(2.0464508) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6948918) q[3];
sx q[3];
rz(-2.4467203) q[3];
sx q[3];
rz(-1.2563039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9775057) q[2];
sx q[2];
rz(-1.4848494) q[2];
sx q[2];
rz(0.40317765) q[2];
rz(-2.4781503) q[3];
sx q[3];
rz(-2.5622029) q[3];
sx q[3];
rz(-3.1373533) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4106301) q[0];
sx q[0];
rz(-1.0799438) q[0];
sx q[0];
rz(-0.078911111) q[0];
rz(0.90947378) q[1];
sx q[1];
rz(-2.0957004) q[1];
sx q[1];
rz(-0.39631072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8412668) q[0];
sx q[0];
rz(-1.5551651) q[0];
sx q[0];
rz(-0.0090740694) q[0];
rz(0.79371039) q[2];
sx q[2];
rz(-1.5657445) q[2];
sx q[2];
rz(1.9144626) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.24416284) q[1];
sx q[1];
rz(-1.8051935) q[1];
sx q[1];
rz(0.48903709) q[1];
rz(-pi) q[2];
rz(2.1911133) q[3];
sx q[3];
rz(-1.1071604) q[3];
sx q[3];
rz(0.5666545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4497455) q[2];
sx q[2];
rz(-0.88374603) q[2];
sx q[2];
rz(0.79908243) q[2];
rz(-0.07829047) q[3];
sx q[3];
rz(-1.2543863) q[3];
sx q[3];
rz(-2.4962795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.66452022) q[0];
sx q[0];
rz(-1.7417396) q[0];
sx q[0];
rz(-2.59792) q[0];
rz(-1.1935344) q[1];
sx q[1];
rz(-1.2763034) q[1];
sx q[1];
rz(-2.6313849) q[1];
rz(0.1706201) q[2];
sx q[2];
rz(-2.0224051) q[2];
sx q[2];
rz(-0.97013459) q[2];
rz(2.067749) q[3];
sx q[3];
rz(-2.2324149) q[3];
sx q[3];
rz(0.80692337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
