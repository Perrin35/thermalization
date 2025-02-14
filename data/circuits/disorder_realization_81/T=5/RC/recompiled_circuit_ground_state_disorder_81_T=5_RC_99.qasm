OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12299744) q[0];
sx q[0];
rz(-1.7641492) q[0];
sx q[0];
rz(-2.7121845) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(7.9908854) q[1];
sx q[1];
rz(7.8828852) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.854092) q[0];
sx q[0];
rz(-1.2852291) q[0];
sx q[0];
rz(0.2233611) q[0];
x q[1];
rz(1.6060271) q[2];
sx q[2];
rz(-1.4270947) q[2];
sx q[2];
rz(-1.3006388) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0555901) q[1];
sx q[1];
rz(-1.1430972) q[1];
sx q[1];
rz(1.6040638) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8640932) q[3];
sx q[3];
rz(-1.2240922) q[3];
sx q[3];
rz(-1.5145921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30449197) q[2];
sx q[2];
rz(-1.203275) q[2];
sx q[2];
rz(2.4486747) q[2];
rz(2.9414226) q[3];
sx q[3];
rz(-2.9514511) q[3];
sx q[3];
rz(-1.1339124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8189341) q[0];
sx q[0];
rz(-2.942473) q[0];
sx q[0];
rz(-1.06485) q[0];
rz(2.5191567) q[1];
sx q[1];
rz(-1.348) q[1];
sx q[1];
rz(-2.650824) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7486742) q[0];
sx q[0];
rz(-3.1160917) q[0];
sx q[0];
rz(2.7703417) q[0];
x q[1];
rz(-2.7193473) q[2];
sx q[2];
rz(-1.7755055) q[2];
sx q[2];
rz(-1.7038356) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.307617) q[1];
sx q[1];
rz(-0.81980356) q[1];
sx q[1];
rz(-1.4403254) q[1];
x q[2];
rz(1.7401668) q[3];
sx q[3];
rz(-1.6182096) q[3];
sx q[3];
rz(2.5518038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48851594) q[2];
sx q[2];
rz(-1.642903) q[2];
sx q[2];
rz(0.24031362) q[2];
rz(-2.6416685) q[3];
sx q[3];
rz(-0.58568716) q[3];
sx q[3];
rz(1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2568473) q[0];
sx q[0];
rz(-2.186543) q[0];
sx q[0];
rz(0.98168674) q[0];
rz(-0.49199545) q[1];
sx q[1];
rz(-2.056608) q[1];
sx q[1];
rz(1.5031776) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5971736) q[0];
sx q[0];
rz(-2.0453296) q[0];
sx q[0];
rz(-2.3719792) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.515834) q[2];
sx q[2];
rz(-1.171249) q[2];
sx q[2];
rz(1.1732701) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43736378) q[1];
sx q[1];
rz(-1.2384602) q[1];
sx q[1];
rz(1.9802914) q[1];
rz(1.668456) q[3];
sx q[3];
rz(-1.5916232) q[3];
sx q[3];
rz(-1.5725003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65695277) q[2];
sx q[2];
rz(-0.52565614) q[2];
sx q[2];
rz(-3.0204115) q[2];
rz(2.1335404) q[3];
sx q[3];
rz(-1.1153778) q[3];
sx q[3];
rz(1.0813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-2.0693531) q[0];
sx q[0];
rz(-0.00088748137) q[0];
sx q[0];
rz(2.6278611) q[0];
rz(-0.95798245) q[1];
sx q[1];
rz(-1.999141) q[1];
sx q[1];
rz(1.4428008) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7126677) q[0];
sx q[0];
rz(-2.6972572) q[0];
sx q[0];
rz(0.39074583) q[0];
x q[1];
rz(1.7920234) q[2];
sx q[2];
rz(-2.0906679) q[2];
sx q[2];
rz(1.2017565) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.629238) q[1];
sx q[1];
rz(-1.6529473) q[1];
sx q[1];
rz(2.1869529) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.47993) q[3];
sx q[3];
rz(-1.6473624) q[3];
sx q[3];
rz(-2.9947502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5865667) q[2];
sx q[2];
rz(-2.1777966) q[2];
sx q[2];
rz(1.8761934) q[2];
rz(-1.966656) q[3];
sx q[3];
rz(-2.411071) q[3];
sx q[3];
rz(-1.1635715) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2555399) q[0];
sx q[0];
rz(-0.20296725) q[0];
sx q[0];
rz(1.8170005) q[0];
rz(-0.96877226) q[1];
sx q[1];
rz(-1.4854393) q[1];
sx q[1];
rz(1.0383777) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7234112) q[0];
sx q[0];
rz(-1.1901642) q[0];
sx q[0];
rz(1.7443954) q[0];
rz(2.3164301) q[2];
sx q[2];
rz(-1.0010011) q[2];
sx q[2];
rz(-1.2078326) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6250302) q[1];
sx q[1];
rz(-1.5296401) q[1];
sx q[1];
rz(-2.7709318) q[1];
x q[2];
rz(1.0594924) q[3];
sx q[3];
rz(-0.97701525) q[3];
sx q[3];
rz(0.64835129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5021299) q[2];
sx q[2];
rz(-0.30439964) q[2];
sx q[2];
rz(0.89140618) q[2];
rz(-1.2616448) q[3];
sx q[3];
rz(-1.6311878) q[3];
sx q[3];
rz(2.4752899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0797743) q[0];
sx q[0];
rz(-0.49538716) q[0];
sx q[0];
rz(-0.99639446) q[0];
rz(2.2415316) q[1];
sx q[1];
rz(-2.3521017) q[1];
sx q[1];
rz(0.17734227) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19939199) q[0];
sx q[0];
rz(-0.15155242) q[0];
sx q[0];
rz(-1.5723537) q[0];
x q[1];
rz(-1.5992237) q[2];
sx q[2];
rz(-0.6387944) q[2];
sx q[2];
rz(2.6306689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4882433) q[1];
sx q[1];
rz(-1.0945787) q[1];
sx q[1];
rz(-2.819092) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0445339) q[3];
sx q[3];
rz(-1.0639816) q[3];
sx q[3];
rz(-0.34891303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8511054) q[2];
sx q[2];
rz(-1.1815716) q[2];
sx q[2];
rz(2.8372852) q[2];
rz(0.32637063) q[3];
sx q[3];
rz(-0.54309741) q[3];
sx q[3];
rz(0.91199005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.07311634) q[0];
sx q[0];
rz(-2.731972) q[0];
sx q[0];
rz(-2.2090744) q[0];
rz(-0.069123507) q[1];
sx q[1];
rz(-0.2176452) q[1];
sx q[1];
rz(2.1014012) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0737586) q[0];
sx q[0];
rz(-0.76601765) q[0];
sx q[0];
rz(-3.104408) q[0];
x q[1];
rz(-2.9670976) q[2];
sx q[2];
rz(-1.3180079) q[2];
sx q[2];
rz(-2.6215907) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55807796) q[1];
sx q[1];
rz(-0.78418523) q[1];
sx q[1];
rz(0.23434831) q[1];
rz(-pi) q[2];
rz(-0.64994855) q[3];
sx q[3];
rz(-2.3293265) q[3];
sx q[3];
rz(2.5358729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.13504623) q[2];
sx q[2];
rz(-0.73837215) q[2];
sx q[2];
rz(0.73299232) q[2];
rz(1.0829571) q[3];
sx q[3];
rz(-2.0854009) q[3];
sx q[3];
rz(-2.1347031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873544) q[0];
sx q[0];
rz(-3.0576958) q[0];
sx q[0];
rz(0.49302897) q[0];
rz(2.9250277) q[1];
sx q[1];
rz(-1.1410057) q[1];
sx q[1];
rz(-2.1176178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1676583) q[0];
sx q[0];
rz(-1.5827145) q[0];
sx q[0];
rz(2.5456136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6844953) q[2];
sx q[2];
rz(-1.9010882) q[2];
sx q[2];
rz(1.2680666) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1223381) q[1];
sx q[1];
rz(-2.755909) q[1];
sx q[1];
rz(-0.38094873) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29924691) q[3];
sx q[3];
rz(-0.95029921) q[3];
sx q[3];
rz(2.3503691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0852802) q[2];
sx q[2];
rz(-1.6322501) q[2];
sx q[2];
rz(-0.67576605) q[2];
rz(0.41302776) q[3];
sx q[3];
rz(-1.4512117) q[3];
sx q[3];
rz(0.54615027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6290879) q[0];
sx q[0];
rz(-2.4830723) q[0];
sx q[0];
rz(0.034544695) q[0];
rz(-3.0870364) q[1];
sx q[1];
rz(-1.1628393) q[1];
sx q[1];
rz(-2.3202855) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23653447) q[0];
sx q[0];
rz(-2.5590959) q[0];
sx q[0];
rz(0.18180099) q[0];
rz(-3.05282) q[2];
sx q[2];
rz(-1.5615511) q[2];
sx q[2];
rz(2.451596) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2945613) q[1];
sx q[1];
rz(-2.9518513) q[1];
sx q[1];
rz(-2.4824597) q[1];
x q[2];
rz(2.6590632) q[3];
sx q[3];
rz(-2.165613) q[3];
sx q[3];
rz(-1.9687259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7811232) q[2];
sx q[2];
rz(-1.0162153) q[2];
sx q[2];
rz(0.13295573) q[2];
rz(-0.41108701) q[3];
sx q[3];
rz(-1.6229595) q[3];
sx q[3];
rz(-1.0857922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046831176) q[0];
sx q[0];
rz(-0.60698858) q[0];
sx q[0];
rz(1.2588311) q[0];
rz(-1.6350485) q[1];
sx q[1];
rz(-2.4848487) q[1];
sx q[1];
rz(0.81319317) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55594873) q[0];
sx q[0];
rz(-0.085594479) q[0];
sx q[0];
rz(-2.7187347) q[0];
x q[1];
rz(-1.4237464) q[2];
sx q[2];
rz(-0.67710256) q[2];
sx q[2];
rz(-2.8092217) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.836536) q[1];
sx q[1];
rz(-0.99545762) q[1];
sx q[1];
rz(2.1161377) q[1];
rz(2.68937) q[3];
sx q[3];
rz(-1.6912231) q[3];
sx q[3];
rz(2.6129006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4930111) q[2];
sx q[2];
rz(-1.1375256) q[2];
sx q[2];
rz(-1.0515593) q[2];
rz(1.018853) q[3];
sx q[3];
rz(-0.84210432) q[3];
sx q[3];
rz(-2.4837608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1558341) q[0];
sx q[0];
rz(-2.4875165) q[0];
sx q[0];
rz(-2.2610337) q[0];
rz(-1.8219933) q[1];
sx q[1];
rz(-1.1450014) q[1];
sx q[1];
rz(-3.0897279) q[1];
rz(0.20959494) q[2];
sx q[2];
rz(-1.4996698) q[2];
sx q[2];
rz(-1.3347129) q[2];
rz(-2.8494358) q[3];
sx q[3];
rz(-2.9360129) q[3];
sx q[3];
rz(2.871411) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
