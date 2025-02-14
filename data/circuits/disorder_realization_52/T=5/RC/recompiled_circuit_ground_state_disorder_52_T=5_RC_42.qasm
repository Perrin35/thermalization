OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0527394) q[0];
sx q[0];
rz(-2.7272447) q[0];
sx q[0];
rz(-2.156884) q[0];
rz(2.2766278) q[1];
sx q[1];
rz(-2.2948761) q[1];
sx q[1];
rz(2.4196978) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5402031) q[0];
sx q[0];
rz(-2.5059675) q[0];
sx q[0];
rz(0.74936015) q[0];
x q[1];
rz(0.38227947) q[2];
sx q[2];
rz(-2.3685902) q[2];
sx q[2];
rz(-0.24958463) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5923323) q[1];
sx q[1];
rz(-1.499814) q[1];
sx q[1];
rz(0.029220079) q[1];
x q[2];
rz(-1.2087565) q[3];
sx q[3];
rz(-2.7305121) q[3];
sx q[3];
rz(-0.53888881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7559173) q[2];
sx q[2];
rz(-2.6188681) q[2];
sx q[2];
rz(-0.44504607) q[2];
rz(1.2612777) q[3];
sx q[3];
rz(-1.4538366) q[3];
sx q[3];
rz(-1.5104347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.6181013) q[0];
sx q[0];
rz(-1.0502676) q[0];
sx q[0];
rz(-0.65446788) q[0];
rz(2.0596313) q[1];
sx q[1];
rz(-1.8257717) q[1];
sx q[1];
rz(1.6771603) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0078115) q[0];
sx q[0];
rz(-2.0424583) q[0];
sx q[0];
rz(-0.056344294) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6763386) q[2];
sx q[2];
rz(-0.60655884) q[2];
sx q[2];
rz(-2.7594729) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38514458) q[1];
sx q[1];
rz(-0.83122453) q[1];
sx q[1];
rz(-1.6265197) q[1];
x q[2];
rz(2.5188279) q[3];
sx q[3];
rz(-0.55985427) q[3];
sx q[3];
rz(0.44677681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.64708853) q[2];
sx q[2];
rz(-0.096179811) q[2];
sx q[2];
rz(-0.75867009) q[2];
rz(-1.8132973) q[3];
sx q[3];
rz(-2.0601065) q[3];
sx q[3];
rz(-1.7648511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7741622) q[0];
sx q[0];
rz(-1.6663015) q[0];
sx q[0];
rz(3.1040763) q[0];
rz(-2.1027749) q[1];
sx q[1];
rz(-2.807834) q[1];
sx q[1];
rz(-0.85938251) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5512269) q[0];
sx q[0];
rz(-0.96053505) q[0];
sx q[0];
rz(1.8009737) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3808566) q[2];
sx q[2];
rz(-1.4773792) q[2];
sx q[2];
rz(0.38047516) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22166477) q[1];
sx q[1];
rz(-1.8916191) q[1];
sx q[1];
rz(-2.0179835) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66481164) q[3];
sx q[3];
rz(-0.71352173) q[3];
sx q[3];
rz(-0.87687311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77217707) q[2];
sx q[2];
rz(-2.206216) q[2];
sx q[2];
rz(-2.3773362) q[2];
rz(-0.94591004) q[3];
sx q[3];
rz(-1.2063682) q[3];
sx q[3];
rz(-0.8980155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1979444) q[0];
sx q[0];
rz(-0.89364377) q[0];
sx q[0];
rz(1.8290895) q[0];
rz(0.30458826) q[1];
sx q[1];
rz(-2.1423788) q[1];
sx q[1];
rz(0.33590683) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2755796) q[0];
sx q[0];
rz(-0.073750138) q[0];
sx q[0];
rz(-1.3763675) q[0];
x q[1];
rz(2.1763889) q[2];
sx q[2];
rz(-1.196567) q[2];
sx q[2];
rz(0.40520129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2780032) q[1];
sx q[1];
rz(-1.9801323) q[1];
sx q[1];
rz(0.54249842) q[1];
rz(1.6275938) q[3];
sx q[3];
rz(-1.3496163) q[3];
sx q[3];
rz(3.055453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37561068) q[2];
sx q[2];
rz(-0.42372647) q[2];
sx q[2];
rz(1.83164) q[2];
rz(1.2951819) q[3];
sx q[3];
rz(-2.3056307) q[3];
sx q[3];
rz(1.0999934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.497371) q[0];
sx q[0];
rz(-0.27232429) q[0];
sx q[0];
rz(2.3261133) q[0];
rz(0.71715912) q[1];
sx q[1];
rz(-1.6970789) q[1];
sx q[1];
rz(-1.1517634) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8659521) q[0];
sx q[0];
rz(-1.3530088) q[0];
sx q[0];
rz(2.6135246) q[0];
rz(1.3281335) q[2];
sx q[2];
rz(-1.8662226) q[2];
sx q[2];
rz(0.60305078) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.30241769) q[1];
sx q[1];
rz(-1.2726899) q[1];
sx q[1];
rz(0.83946622) q[1];
rz(2.8397039) q[3];
sx q[3];
rz(-1.535901) q[3];
sx q[3];
rz(-0.85012415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0723476) q[2];
sx q[2];
rz(-2.3657511) q[2];
sx q[2];
rz(-1.5838712) q[2];
rz(-2.1690058) q[3];
sx q[3];
rz(-1.9105304) q[3];
sx q[3];
rz(-1.4371654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4438181) q[0];
sx q[0];
rz(-0.2727209) q[0];
sx q[0];
rz(-2.4716603) q[0];
rz(2.2386235) q[1];
sx q[1];
rz(-2.4882856) q[1];
sx q[1];
rz(3.0526551) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1207711) q[0];
sx q[0];
rz(-2.1820314) q[0];
sx q[0];
rz(2.8957248) q[0];
rz(-pi) q[1];
rz(-1.3207664) q[2];
sx q[2];
rz(-2.2019115) q[2];
sx q[2];
rz(2.4970412) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8008314) q[1];
sx q[1];
rz(-0.25909875) q[1];
sx q[1];
rz(-1.8482261) q[1];
rz(-pi) q[2];
rz(-3.012978) q[3];
sx q[3];
rz(-2.0484784) q[3];
sx q[3];
rz(-1.1123808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4778121) q[2];
sx q[2];
rz(-1.0441484) q[2];
sx q[2];
rz(0.021154724) q[2];
rz(0.92709213) q[3];
sx q[3];
rz(-1.6095716) q[3];
sx q[3];
rz(0.65897861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0878736) q[0];
sx q[0];
rz(-1.3608195) q[0];
sx q[0];
rz(-1.1370283) q[0];
rz(-2.673705) q[1];
sx q[1];
rz(-1.4257905) q[1];
sx q[1];
rz(0.55380026) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1374986) q[0];
sx q[0];
rz(-1.2826254) q[0];
sx q[0];
rz(0.96364809) q[0];
rz(-pi) q[1];
rz(2.7871889) q[2];
sx q[2];
rz(-1.1605673) q[2];
sx q[2];
rz(-1.1910476) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3674161) q[1];
sx q[1];
rz(-1.2700081) q[1];
sx q[1];
rz(1.0298301) q[1];
rz(-pi) q[2];
rz(0.061896996) q[3];
sx q[3];
rz(-2.1445285) q[3];
sx q[3];
rz(-5*pi/6) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0849453) q[2];
sx q[2];
rz(-2.2979996) q[2];
sx q[2];
rz(-2.4455369) q[2];
rz(-2.7737235) q[3];
sx q[3];
rz(-2.0335679) q[3];
sx q[3];
rz(-2.8554816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2680161) q[0];
sx q[0];
rz(-1.753267) q[0];
sx q[0];
rz(0.8152813) q[0];
rz(0.45375991) q[1];
sx q[1];
rz(-0.90180698) q[1];
sx q[1];
rz(-2.544983) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94039791) q[0];
sx q[0];
rz(-1.0642645) q[0];
sx q[0];
rz(0.80123637) q[0];
rz(-pi) q[1];
rz(2.3109092) q[2];
sx q[2];
rz(-1.4998337) q[2];
sx q[2];
rz(1.1712779) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69747671) q[1];
sx q[1];
rz(-1.2713065) q[1];
sx q[1];
rz(-2.1397352) q[1];
x q[2];
rz(-0.96489789) q[3];
sx q[3];
rz(-2.4526074) q[3];
sx q[3];
rz(0.24619547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0348908) q[2];
sx q[2];
rz(-1.7009578) q[2];
sx q[2];
rz(-1.9218669) q[2];
rz(1.2446416) q[3];
sx q[3];
rz(-1.605426) q[3];
sx q[3];
rz(-0.19967782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.225746) q[0];
sx q[0];
rz(-0.93996489) q[0];
sx q[0];
rz(-1.655727) q[0];
rz(-1.111221) q[1];
sx q[1];
rz(-1.8467555) q[1];
sx q[1];
rz(1.750754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6538453) q[0];
sx q[0];
rz(-1.5731205) q[0];
sx q[0];
rz(0.33727686) q[0];
rz(-pi) q[1];
rz(-2.7971091) q[2];
sx q[2];
rz(-1.4622258) q[2];
sx q[2];
rz(0.44652938) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7906906) q[1];
sx q[1];
rz(-1.6054253) q[1];
sx q[1];
rz(-1.778855) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2038571) q[3];
sx q[3];
rz(-1.8979567) q[3];
sx q[3];
rz(0.1042455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.86604649) q[2];
sx q[2];
rz(-1.486843) q[2];
sx q[2];
rz(0.1869959) q[2];
rz(0.20728076) q[3];
sx q[3];
rz(-0.63396251) q[3];
sx q[3];
rz(2.0980339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2279376) q[0];
sx q[0];
rz(-2.3917103) q[0];
sx q[0];
rz(-3.1392745) q[0];
rz(-1.9193316) q[1];
sx q[1];
rz(-1.5879251) q[1];
sx q[1];
rz(1.7244171) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9681536) q[0];
sx q[0];
rz(-1.441297) q[0];
sx q[0];
rz(-0.20682516) q[0];
rz(0.50062407) q[2];
sx q[2];
rz(-2.5412895) q[2];
sx q[2];
rz(-0.17498091) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0884888) q[1];
sx q[1];
rz(-1.3657059) q[1];
sx q[1];
rz(0.049909485) q[1];
x q[2];
rz(-2.9692698) q[3];
sx q[3];
rz(-2.1502521) q[3];
sx q[3];
rz(-0.31238279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2554539) q[2];
sx q[2];
rz(-1.5172493) q[2];
sx q[2];
rz(-0.21190602) q[2];
rz(2.741277) q[3];
sx q[3];
rz(-0.90462697) q[3];
sx q[3];
rz(-1.3306085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048007456) q[0];
sx q[0];
rz(-1.891991) q[0];
sx q[0];
rz(1.4461507) q[0];
rz(-2.6493337) q[1];
sx q[1];
rz(-1.6620363) q[1];
sx q[1];
rz(-1.5580039) q[1];
rz(-2.9335446) q[2];
sx q[2];
rz(-1.3606417) q[2];
sx q[2];
rz(-0.66858895) q[2];
rz(-0.51385469) q[3];
sx q[3];
rz(-1.1198291) q[3];
sx q[3];
rz(-1.8201943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
