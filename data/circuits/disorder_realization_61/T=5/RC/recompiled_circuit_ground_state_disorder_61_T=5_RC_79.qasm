OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0680024) q[0];
sx q[0];
rz(-1.468714) q[0];
sx q[0];
rz(0.52748632) q[0];
rz(-2.5975851) q[1];
sx q[1];
rz(-0.44096947) q[1];
sx q[1];
rz(0.0013594065) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0802976) q[0];
sx q[0];
rz(-1.5004145) q[0];
sx q[0];
rz(-0.0091307739) q[0];
rz(-2.8090879) q[2];
sx q[2];
rz(-2.5026179) q[2];
sx q[2];
rz(-0.19579355) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4365523) q[1];
sx q[1];
rz(-2.4715354) q[1];
sx q[1];
rz(2.4923508) q[1];
x q[2];
rz(1.6) q[3];
sx q[3];
rz(-1.3714223) q[3];
sx q[3];
rz(0.77413156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9881543) q[2];
sx q[2];
rz(-2.4997288) q[2];
sx q[2];
rz(-1.7786857) q[2];
rz(0.07987944) q[3];
sx q[3];
rz(-2.2814543) q[3];
sx q[3];
rz(3.1352654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7774571) q[0];
sx q[0];
rz(-1.2292925) q[0];
sx q[0];
rz(-0.51032132) q[0];
rz(1.3279042) q[1];
sx q[1];
rz(-2.4102305) q[1];
sx q[1];
rz(2.7046943) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6104159) q[0];
sx q[0];
rz(-2.057909) q[0];
sx q[0];
rz(0.25419828) q[0];
rz(-pi) q[1];
rz(-1.8067752) q[2];
sx q[2];
rz(-2.0046104) q[2];
sx q[2];
rz(3.0133875) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.092097923) q[1];
sx q[1];
rz(-1.5391401) q[1];
sx q[1];
rz(2.3718216) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9741967) q[3];
sx q[3];
rz(-1.6296436) q[3];
sx q[3];
rz(1.7047391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83866155) q[2];
sx q[2];
rz(-2.624056) q[2];
sx q[2];
rz(-2.9663864) q[2];
rz(0.13904275) q[3];
sx q[3];
rz(-1.6486721) q[3];
sx q[3];
rz(0.89865941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43056968) q[0];
sx q[0];
rz(-0.88453203) q[0];
sx q[0];
rz(2.7810466) q[0];
rz(1.4973466) q[1];
sx q[1];
rz(-1.2253573) q[1];
sx q[1];
rz(-2.5490733) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07537341) q[0];
sx q[0];
rz(-0.96642113) q[0];
sx q[0];
rz(0.79232596) q[0];
x q[1];
rz(-2.9354344) q[2];
sx q[2];
rz(-0.64616424) q[2];
sx q[2];
rz(0.29096067) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6806023) q[1];
sx q[1];
rz(-0.36961886) q[1];
sx q[1];
rz(1.2301983) q[1];
x q[2];
rz(2.4094916) q[3];
sx q[3];
rz(-0.59716254) q[3];
sx q[3];
rz(2.5954318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3561919) q[2];
sx q[2];
rz(-2.5742026) q[2];
sx q[2];
rz(-0.94245911) q[2];
rz(2.5475907) q[3];
sx q[3];
rz(-2.3022251) q[3];
sx q[3];
rz(0.12282898) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9451611) q[0];
sx q[0];
rz(-2.5158947) q[0];
sx q[0];
rz(-2.4932267) q[0];
rz(-2.2658589) q[1];
sx q[1];
rz(-1.6301194) q[1];
sx q[1];
rz(1.7902364) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60015357) q[0];
sx q[0];
rz(-2.3838288) q[0];
sx q[0];
rz(1.6275854) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8114795) q[2];
sx q[2];
rz(-1.7042024) q[2];
sx q[2];
rz(-1.2776983) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9571805) q[1];
sx q[1];
rz(-2.4676128) q[1];
sx q[1];
rz(1.5247702) q[1];
rz(-pi) q[2];
rz(-0.13552119) q[3];
sx q[3];
rz(-0.38619872) q[3];
sx q[3];
rz(3.112325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9883092) q[2];
sx q[2];
rz(-0.51258665) q[2];
sx q[2];
rz(1.2549866) q[2];
rz(-0.60162383) q[3];
sx q[3];
rz(-2.2646077) q[3];
sx q[3];
rz(1.4582483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.0643784) q[0];
sx q[0];
rz(-2.0266396) q[0];
sx q[0];
rz(3.0078122) q[0];
rz(-1.4859707) q[1];
sx q[1];
rz(-0.32951117) q[1];
sx q[1];
rz(2.9697184) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8012708) q[0];
sx q[0];
rz(-1.9376799) q[0];
sx q[0];
rz(-1.0563068) q[0];
rz(-pi) q[1];
rz(1.1823229) q[2];
sx q[2];
rz(-2.3459593) q[2];
sx q[2];
rz(2.1284136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1273074) q[1];
sx q[1];
rz(-1.2218857) q[1];
sx q[1];
rz(2.9800913) q[1];
rz(-pi) q[2];
rz(0.47953812) q[3];
sx q[3];
rz(-1.1894636) q[3];
sx q[3];
rz(1.4863734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6493426) q[2];
sx q[2];
rz(-0.52521962) q[2];
sx q[2];
rz(1.4225175) q[2];
rz(-2.9366142) q[3];
sx q[3];
rz(-2.180438) q[3];
sx q[3];
rz(-2.3375296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64254665) q[0];
sx q[0];
rz(-2.5048984) q[0];
sx q[0];
rz(0.194304) q[0];
rz(-2.4957472) q[1];
sx q[1];
rz(-1.9648809) q[1];
sx q[1];
rz(-2.7913854) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6381643) q[0];
sx q[0];
rz(-1.0409875) q[0];
sx q[0];
rz(1.4132771) q[0];
rz(-pi) q[1];
rz(-2.3268125) q[2];
sx q[2];
rz(-1.0836923) q[2];
sx q[2];
rz(1.8419982) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0547452) q[1];
sx q[1];
rz(-1.655181) q[1];
sx q[1];
rz(-2.6230472) q[1];
rz(-pi) q[2];
rz(-3.0549235) q[3];
sx q[3];
rz(-2.4695167) q[3];
sx q[3];
rz(1.8928526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.084006) q[2];
sx q[2];
rz(-0.39613327) q[2];
sx q[2];
rz(0.32220379) q[2];
rz(-2.3816439) q[3];
sx q[3];
rz(-2.6252169) q[3];
sx q[3];
rz(-3.0907536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8754804) q[0];
sx q[0];
rz(-0.66829824) q[0];
sx q[0];
rz(-0.71568263) q[0];
rz(-3.074805) q[1];
sx q[1];
rz(-0.5611493) q[1];
sx q[1];
rz(0.36398789) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92236906) q[0];
sx q[0];
rz(-2.4441193) q[0];
sx q[0];
rz(1.0154025) q[0];
rz(-pi) q[1];
rz(-2.5189713) q[2];
sx q[2];
rz(-1.004815) q[2];
sx q[2];
rz(-0.37812585) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1541031) q[1];
sx q[1];
rz(-1.5871154) q[1];
sx q[1];
rz(0.76413566) q[1];
x q[2];
rz(2.4956483) q[3];
sx q[3];
rz(-1.6872354) q[3];
sx q[3];
rz(-1.2723591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5169342) q[2];
sx q[2];
rz(-0.90183574) q[2];
sx q[2];
rz(-1.2348403) q[2];
rz(0.45461795) q[3];
sx q[3];
rz(-1.8574628) q[3];
sx q[3];
rz(-0.07479085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.5816077) q[0];
sx q[0];
rz(-3.0261664) q[0];
sx q[0];
rz(-2.5009632) q[0];
rz(-0.36264125) q[1];
sx q[1];
rz(-2.8022712) q[1];
sx q[1];
rz(1.6222662) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9576833) q[0];
sx q[0];
rz(-1.0445692) q[0];
sx q[0];
rz(-2.7503114) q[0];
rz(-pi) q[1];
rz(-2.4108886) q[2];
sx q[2];
rz(-1.7901229) q[2];
sx q[2];
rz(-2.4958269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.84953972) q[1];
sx q[1];
rz(-0.89703945) q[1];
sx q[1];
rz(-0.11815355) q[1];
rz(-1.7666807) q[3];
sx q[3];
rz(-2.6614411) q[3];
sx q[3];
rz(-2.5993915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92056876) q[2];
sx q[2];
rz(-1.8683044) q[2];
sx q[2];
rz(-1.9839239) q[2];
rz(-2.8122592) q[3];
sx q[3];
rz(-2.9521827) q[3];
sx q[3];
rz(0.91387373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2374903) q[0];
sx q[0];
rz(-0.79039031) q[0];
sx q[0];
rz(0.41123408) q[0];
rz(2.5143738) q[1];
sx q[1];
rz(-2.5721481) q[1];
sx q[1];
rz(-1.8597182) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86253765) q[0];
sx q[0];
rz(-1.9337961) q[0];
sx q[0];
rz(-1.3393066) q[0];
rz(-pi) q[1];
rz(2.5847264) q[2];
sx q[2];
rz(-2.1376196) q[2];
sx q[2];
rz(-1.8447529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3157278) q[1];
sx q[1];
rz(-2.4701354) q[1];
sx q[1];
rz(-1.3888098) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9427845) q[3];
sx q[3];
rz(-1.947177) q[3];
sx q[3];
rz(2.876296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7810479) q[2];
sx q[2];
rz(-1.5881528) q[2];
sx q[2];
rz(-2.1629199) q[2];
rz(0.22748889) q[3];
sx q[3];
rz(-2.3922908) q[3];
sx q[3];
rz(-2.6849875) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50642538) q[0];
sx q[0];
rz(-3.083482) q[0];
sx q[0];
rz(0.79570049) q[0];
rz(-0.52867633) q[1];
sx q[1];
rz(-2.1541336) q[1];
sx q[1];
rz(2.9472369) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7995509) q[0];
sx q[0];
rz(-0.38368762) q[0];
sx q[0];
rz(1.8162526) q[0];
rz(-2.7493619) q[2];
sx q[2];
rz(-1.298045) q[2];
sx q[2];
rz(-0.14209948) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72994191) q[1];
sx q[1];
rz(-2.6788708) q[1];
sx q[1];
rz(-2.8587841) q[1];
rz(-0.80598237) q[3];
sx q[3];
rz(-2.4641086) q[3];
sx q[3];
rz(-2.4777378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84870321) q[2];
sx q[2];
rz(-1.6699426) q[2];
sx q[2];
rz(-0.29937747) q[2];
rz(0.40144604) q[3];
sx q[3];
rz(-0.58599389) q[3];
sx q[3];
rz(2.6011023) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15783475) q[0];
sx q[0];
rz(-0.83136375) q[0];
sx q[0];
rz(-1.2663483) q[0];
rz(-0.10706317) q[1];
sx q[1];
rz(-2.3657847) q[1];
sx q[1];
rz(1.8484144) q[1];
rz(0.54218311) q[2];
sx q[2];
rz(-1.5246921) q[2];
sx q[2];
rz(0.17490457) q[2];
rz(-0.43184256) q[3];
sx q[3];
rz(-2.8819537) q[3];
sx q[3];
rz(1.0859539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
