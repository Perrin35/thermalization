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
rz(-2.1844644) q[0];
sx q[0];
rz(-1.6532093) q[0];
sx q[0];
rz(-1.1487577) q[0];
rz(0.59612885) q[1];
sx q[1];
rz(-0.40607536) q[1];
sx q[1];
rz(-1.4546855) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5243149) q[0];
sx q[0];
rz(-1.9348659) q[0];
sx q[0];
rz(0.84366701) q[0];
rz(1.4613011) q[2];
sx q[2];
rz(-0.81297648) q[2];
sx q[2];
rz(2.2535498) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0599938) q[1];
sx q[1];
rz(-2.3761056) q[1];
sx q[1];
rz(-0.19123921) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0166854) q[3];
sx q[3];
rz(-2.4797399) q[3];
sx q[3];
rz(3.007909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9601606) q[2];
sx q[2];
rz(-2.1655607) q[2];
sx q[2];
rz(-1.206548) q[2];
rz(-1.9801961) q[3];
sx q[3];
rz(-2.5738218) q[3];
sx q[3];
rz(1.5709706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96381617) q[0];
sx q[0];
rz(-2.0160567) q[0];
sx q[0];
rz(-0.97208446) q[0];
rz(-2.8731335) q[1];
sx q[1];
rz(-1.700289) q[1];
sx q[1];
rz(-0.45480248) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9869315) q[0];
sx q[0];
rz(-1.4442181) q[0];
sx q[0];
rz(-1.2722227) q[0];
rz(-pi) q[1];
rz(1.9591397) q[2];
sx q[2];
rz(-1.4253311) q[2];
sx q[2];
rz(-0.59284537) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6356335) q[1];
sx q[1];
rz(-1.1065496) q[1];
sx q[1];
rz(-1.6281158) q[1];
x q[2];
rz(0.56549325) q[3];
sx q[3];
rz(-1.0988937) q[3];
sx q[3];
rz(1.7902692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0316169) q[2];
sx q[2];
rz(-2.4133108) q[2];
sx q[2];
rz(0.99962437) q[2];
rz(0.085974606) q[3];
sx q[3];
rz(-1.5558259) q[3];
sx q[3];
rz(-0.011367817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40083945) q[0];
sx q[0];
rz(-1.4680306) q[0];
sx q[0];
rz(-2.9174347) q[0];
rz(-1.5466746) q[1];
sx q[1];
rz(-1.5432576) q[1];
sx q[1];
rz(-1.8348144) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49840333) q[0];
sx q[0];
rz(-2.1191594) q[0];
sx q[0];
rz(-1.4275309) q[0];
rz(0.73422696) q[2];
sx q[2];
rz(-0.9817183) q[2];
sx q[2];
rz(0.40225077) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.94470384) q[1];
sx q[1];
rz(-0.60736362) q[1];
sx q[1];
rz(0.783226) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66931386) q[3];
sx q[3];
rz(-1.770854) q[3];
sx q[3];
rz(1.1049096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.811502) q[2];
sx q[2];
rz(-1.6197438) q[2];
sx q[2];
rz(2.1211993) q[2];
rz(-1.8118106) q[3];
sx q[3];
rz(-2.0164169) q[3];
sx q[3];
rz(-1.6167195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40600768) q[0];
sx q[0];
rz(-1.0574874) q[0];
sx q[0];
rz(1.965858) q[0];
rz(2.8658087) q[1];
sx q[1];
rz(-0.83668721) q[1];
sx q[1];
rz(-0.63794678) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0887749) q[0];
sx q[0];
rz(-0.034398641) q[0];
sx q[0];
rz(0.31487353) q[0];
rz(0.3023382) q[2];
sx q[2];
rz(-2.4983642) q[2];
sx q[2];
rz(2.0143721) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.040311) q[1];
sx q[1];
rz(-2.4263546) q[1];
sx q[1];
rz(0.50846993) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8627871) q[3];
sx q[3];
rz(-2.8099217) q[3];
sx q[3];
rz(-1.3607894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5059169) q[2];
sx q[2];
rz(-1.5779147) q[2];
sx q[2];
rz(1.6944616) q[2];
rz(2.9315089) q[3];
sx q[3];
rz(-2.688372) q[3];
sx q[3];
rz(-2.213018) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90784812) q[0];
sx q[0];
rz(-2.9892428) q[0];
sx q[0];
rz(1.0688758) q[0];
rz(1.8416038) q[1];
sx q[1];
rz(-1.4062107) q[1];
sx q[1];
rz(1.9357505) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.962709) q[0];
sx q[0];
rz(-0.41085748) q[0];
sx q[0];
rz(1.7715447) q[0];
rz(-pi) q[1];
rz(2.2236125) q[2];
sx q[2];
rz(-1.5783572) q[2];
sx q[2];
rz(-0.078469097) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5609392) q[1];
sx q[1];
rz(-1.4672797) q[1];
sx q[1];
rz(0.85080457) q[1];
x q[2];
rz(0.48128328) q[3];
sx q[3];
rz(-1.7168772) q[3];
sx q[3];
rz(-1.2969601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24870366) q[2];
sx q[2];
rz(-2.7619669) q[2];
sx q[2];
rz(-3.0730263) q[2];
rz(0.93962234) q[3];
sx q[3];
rz(-1.7328123) q[3];
sx q[3];
rz(1.9127964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3786316) q[0];
sx q[0];
rz(-1.7447504) q[0];
sx q[0];
rz(2.6201541) q[0];
rz(1.6261082) q[1];
sx q[1];
rz(-1.3896959) q[1];
sx q[1];
rz(2.5159786) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8396225) q[0];
sx q[0];
rz(-1.5396984) q[0];
sx q[0];
rz(1.0307701) q[0];
rz(-0.37920239) q[2];
sx q[2];
rz(-1.450945) q[2];
sx q[2];
rz(0.7650991) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.82078104) q[1];
sx q[1];
rz(-2.6412475) q[1];
sx q[1];
rz(1.9920177) q[1];
x q[2];
rz(1.5099105) q[3];
sx q[3];
rz(-1.1495483) q[3];
sx q[3];
rz(2.6256068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1336512) q[2];
sx q[2];
rz(-0.59591728) q[2];
sx q[2];
rz(3.0360743) q[2];
rz(0.96406913) q[3];
sx q[3];
rz(-1.9897507) q[3];
sx q[3];
rz(-2.156637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0725726) q[0];
sx q[0];
rz(-2.9264086) q[0];
sx q[0];
rz(-1.7286812) q[0];
rz(2.7627796) q[1];
sx q[1];
rz(-1.8870995) q[1];
sx q[1];
rz(-2.5884195) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2755097) q[0];
sx q[0];
rz(-1.9635154) q[0];
sx q[0];
rz(0.87529542) q[0];
rz(-3.000053) q[2];
sx q[2];
rz(-1.6602548) q[2];
sx q[2];
rz(2.3613561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87482086) q[1];
sx q[1];
rz(-1.2874585) q[1];
sx q[1];
rz(1.5288407) q[1];
x q[2];
rz(-2.2284248) q[3];
sx q[3];
rz(-0.32507691) q[3];
sx q[3];
rz(-2.7108148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6806014) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(-0.41213948) q[2];
rz(-0.66017094) q[3];
sx q[3];
rz(-0.5414525) q[3];
sx q[3];
rz(2.6485543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2039345) q[0];
sx q[0];
rz(-1.0451319) q[0];
sx q[0];
rz(-2.9397021) q[0];
rz(1.0578602) q[1];
sx q[1];
rz(-2.7767534) q[1];
sx q[1];
rz(-2.8813664) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3500811) q[0];
sx q[0];
rz(-0.66257325) q[0];
sx q[0];
rz(-0.64278472) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.260602) q[2];
sx q[2];
rz(-0.63369232) q[2];
sx q[2];
rz(-1.6735759) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59108666) q[1];
sx q[1];
rz(-2.6896853) q[1];
sx q[1];
rz(-0.68303789) q[1];
rz(2.9256938) q[3];
sx q[3];
rz(-0.67189081) q[3];
sx q[3];
rz(3.1300621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1883833) q[2];
sx q[2];
rz(-1.7069495) q[2];
sx q[2];
rz(-2.5592213) q[2];
rz(2.6992056) q[3];
sx q[3];
rz(-1.9059537) q[3];
sx q[3];
rz(-2.3053187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22170947) q[0];
sx q[0];
rz(-1.6136805) q[0];
sx q[0];
rz(-1.7145994) q[0];
rz(-1.3555948) q[1];
sx q[1];
rz(-1.6537063) q[1];
sx q[1];
rz(1.4367163) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4886996) q[0];
sx q[0];
rz(-0.32036361) q[0];
sx q[0];
rz(1.5071177) q[0];
x q[1];
rz(1.8709917) q[2];
sx q[2];
rz(-1.1087024) q[2];
sx q[2];
rz(2.0661046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1940665) q[1];
sx q[1];
rz(-1.3296275) q[1];
sx q[1];
rz(0.014726676) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0441699) q[3];
sx q[3];
rz(-1.7071402) q[3];
sx q[3];
rz(0.52140331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9068678) q[2];
sx q[2];
rz(-1.2805254) q[2];
sx q[2];
rz(1.6024164) q[2];
rz(-0.60798821) q[3];
sx q[3];
rz(-1.4092849) q[3];
sx q[3];
rz(-0.68979818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.36668396) q[0];
sx q[0];
rz(-0.66119778) q[0];
sx q[0];
rz(0.82345024) q[0];
rz(-0.89556328) q[1];
sx q[1];
rz(-1.3500373) q[1];
sx q[1];
rz(2.7604738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51682794) q[0];
sx q[0];
rz(-0.73113686) q[0];
sx q[0];
rz(-2.4419489) q[0];
rz(2.7398749) q[2];
sx q[2];
rz(-1.416393) q[2];
sx q[2];
rz(-0.32333514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5504073) q[1];
sx q[1];
rz(-2.1169615) q[1];
sx q[1];
rz(0.3355432) q[1];
x q[2];
rz(-2.1249843) q[3];
sx q[3];
rz(-1.018152) q[3];
sx q[3];
rz(1.2070398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33619189) q[2];
sx q[2];
rz(-0.41173428) q[2];
sx q[2];
rz(0.98708785) q[2];
rz(-1.6030715) q[3];
sx q[3];
rz(-2.5542732) q[3];
sx q[3];
rz(0.2400329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22454746) q[0];
sx q[0];
rz(-1.2397091) q[0];
sx q[0];
rz(1.5201257) q[0];
rz(-2.4317901) q[1];
sx q[1];
rz(-2.5820422) q[1];
sx q[1];
rz(-1.6834264) q[1];
rz(2.8483656) q[2];
sx q[2];
rz(-0.37920375) q[2];
sx q[2];
rz(-1.3398021) q[2];
rz(-0.98593203) q[3];
sx q[3];
rz(-2.2227047) q[3];
sx q[3];
rz(-0.97490132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
