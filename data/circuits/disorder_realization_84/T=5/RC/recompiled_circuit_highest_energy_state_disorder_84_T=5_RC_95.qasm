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
rz(-0.33997384) q[0];
sx q[0];
rz(-0.61859328) q[0];
sx q[0];
rz(0.11237385) q[0];
rz(0.94315851) q[1];
sx q[1];
rz(4.637742) q[1];
sx q[1];
rz(14.129402) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0709597) q[0];
sx q[0];
rz(-0.68365133) q[0];
sx q[0];
rz(2.4029469) q[0];
rz(-0.19313397) q[2];
sx q[2];
rz(-2.4330957) q[2];
sx q[2];
rz(-0.091781864) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0940548) q[1];
sx q[1];
rz(-0.93701053) q[1];
sx q[1];
rz(-0.82828133) q[1];
rz(-pi) q[2];
rz(-2.2500185) q[3];
sx q[3];
rz(-1.2812616) q[3];
sx q[3];
rz(2.0016157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9872687) q[2];
sx q[2];
rz(-1.4877886) q[2];
sx q[2];
rz(-0.38779116) q[2];
rz(0.44529861) q[3];
sx q[3];
rz(-2.6058091) q[3];
sx q[3];
rz(-2.9545412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33291373) q[0];
sx q[0];
rz(-0.50030047) q[0];
sx q[0];
rz(2.152541) q[0];
rz(-1.5679081) q[1];
sx q[1];
rz(-2.3830919) q[1];
sx q[1];
rz(0.067795098) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2665015) q[0];
sx q[0];
rz(-1.9034732) q[0];
sx q[0];
rz(1.4473296) q[0];
rz(-pi) q[1];
rz(0.8348454) q[2];
sx q[2];
rz(-2.7335484) q[2];
sx q[2];
rz(2.7960475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9631182) q[1];
sx q[1];
rz(-0.97260288) q[1];
sx q[1];
rz(1.7439227) q[1];
rz(3.0127497) q[3];
sx q[3];
rz(-0.93183231) q[3];
sx q[3];
rz(2.2966677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0109791) q[2];
sx q[2];
rz(-2.4108672) q[2];
sx q[2];
rz(2.0558426) q[2];
rz(-2.7340414) q[3];
sx q[3];
rz(-1.4980039) q[3];
sx q[3];
rz(1.3688068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8959592) q[0];
sx q[0];
rz(-2.2838554) q[0];
sx q[0];
rz(1.4602383) q[0];
rz(-2.0270089) q[1];
sx q[1];
rz(-2.3268301) q[1];
sx q[1];
rz(2.1688555) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4125318) q[0];
sx q[0];
rz(-1.5648989) q[0];
sx q[0];
rz(0.062882857) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7454992) q[2];
sx q[2];
rz(-2.1427769) q[2];
sx q[2];
rz(0.33656578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3537843) q[1];
sx q[1];
rz(-2.5671509) q[1];
sx q[1];
rz(0.30409388) q[1];
rz(-pi) q[2];
rz(2.2422355) q[3];
sx q[3];
rz(-1.4345508) q[3];
sx q[3];
rz(1.5929102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14791791) q[2];
sx q[2];
rz(-2.6213054) q[2];
sx q[2];
rz(2.6479123) q[2];
rz(0.70519051) q[3];
sx q[3];
rz(-2.5934936) q[3];
sx q[3];
rz(-1.7459474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.482835) q[0];
sx q[0];
rz(-2.5107497) q[0];
sx q[0];
rz(0.054542907) q[0];
rz(1.4788117) q[1];
sx q[1];
rz(-0.80816591) q[1];
sx q[1];
rz(-0.85173839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8438827) q[0];
sx q[0];
rz(-0.57115924) q[0];
sx q[0];
rz(-0.63262748) q[0];
x q[1];
rz(2.7481467) q[2];
sx q[2];
rz(-1.473765) q[2];
sx q[2];
rz(1.513956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0151685) q[1];
sx q[1];
rz(-0.42747091) q[1];
sx q[1];
rz(0.36001701) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5879166) q[3];
sx q[3];
rz(-2.3683057) q[3];
sx q[3];
rz(2.0594085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30107522) q[2];
sx q[2];
rz(-0.12945759) q[2];
sx q[2];
rz(-1.2811309) q[2];
rz(1.2766131) q[3];
sx q[3];
rz(-1.8254779) q[3];
sx q[3];
rz(0.94094706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5285444) q[0];
sx q[0];
rz(-3.0597882) q[0];
sx q[0];
rz(1.645389) q[0];
rz(1.243535) q[1];
sx q[1];
rz(-1.4356109) q[1];
sx q[1];
rz(0.03874716) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17651788) q[0];
sx q[0];
rz(-1.825188) q[0];
sx q[0];
rz(-0.1381257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1705519) q[2];
sx q[2];
rz(-2.0482302) q[2];
sx q[2];
rz(3.1078452) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97315362) q[1];
sx q[1];
rz(-1.8171696) q[1];
sx q[1];
rz(-2.8484341) q[1];
rz(-pi) q[2];
rz(-0.30743044) q[3];
sx q[3];
rz(-1.0671158) q[3];
sx q[3];
rz(1.4952457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1425928) q[2];
sx q[2];
rz(-2.4685389) q[2];
sx q[2];
rz(-0.9781982) q[2];
rz(0.081196872) q[3];
sx q[3];
rz(-1.2673763) q[3];
sx q[3];
rz(0.80011884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4303495) q[0];
sx q[0];
rz(-0.62380236) q[0];
sx q[0];
rz(-1.2687564) q[0];
rz(-2.2637892) q[1];
sx q[1];
rz(-1.9453847) q[1];
sx q[1];
rz(-2.9315604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18845226) q[0];
sx q[0];
rz(-2.3658381) q[0];
sx q[0];
rz(-0.95073582) q[0];
x q[1];
rz(0.66234373) q[2];
sx q[2];
rz(-2.6731632) q[2];
sx q[2];
rz(1.1520907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.91466838) q[1];
sx q[1];
rz(-1.4836425) q[1];
sx q[1];
rz(0.29030771) q[1];
rz(-pi) q[2];
x q[2];
rz(2.00497) q[3];
sx q[3];
rz(-1.3011271) q[3];
sx q[3];
rz(2.8119905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7741125) q[2];
sx q[2];
rz(-1.3051278) q[2];
sx q[2];
rz(1.0779862) q[2];
rz(-0.20400253) q[3];
sx q[3];
rz(-1.9670468) q[3];
sx q[3];
rz(-1.3998122) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.792895) q[0];
sx q[0];
rz(-2.5070511) q[0];
sx q[0];
rz(0.10575159) q[0];
rz(1.3144685) q[1];
sx q[1];
rz(-1.0377089) q[1];
sx q[1];
rz(0.43397841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2600685) q[0];
sx q[0];
rz(-2.6556042) q[0];
sx q[0];
rz(1.6218779) q[0];
rz(-pi) q[1];
rz(0.28027541) q[2];
sx q[2];
rz(-1.2024101) q[2];
sx q[2];
rz(-1.0944674) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24819198) q[1];
sx q[1];
rz(-1.6576644) q[1];
sx q[1];
rz(-0.76201622) q[1];
rz(-pi) q[2];
rz(0.26386719) q[3];
sx q[3];
rz(-1.4949867) q[3];
sx q[3];
rz(0.77944642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5964261) q[2];
sx q[2];
rz(-0.20238987) q[2];
sx q[2];
rz(-1.8983967) q[2];
rz(3.0498114) q[3];
sx q[3];
rz(-1.4877078) q[3];
sx q[3];
rz(1.6091326) q[3];
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
rz(-1.5246326) q[0];
sx q[0];
rz(-1.4240823) q[0];
sx q[0];
rz(-2.4272163) q[0];
rz(-2.3081035) q[1];
sx q[1];
rz(-1.0316713) q[1];
sx q[1];
rz(-2.3918236) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3692255) q[0];
sx q[0];
rz(-2.0207267) q[0];
sx q[0];
rz(2.1693267) q[0];
x q[1];
rz(-0.60192666) q[2];
sx q[2];
rz(-1.2363803) q[2];
sx q[2];
rz(0.95171164) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52940375) q[1];
sx q[1];
rz(-0.85221186) q[1];
sx q[1];
rz(0.3297594) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.060018221) q[3];
sx q[3];
rz(-1.3207153) q[3];
sx q[3];
rz(1.3986971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82911503) q[2];
sx q[2];
rz(-1.1822327) q[2];
sx q[2];
rz(-1.0324837) q[2];
rz(-1.9898344) q[3];
sx q[3];
rz(-2.0572898) q[3];
sx q[3];
rz(-1.1758218) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51079232) q[0];
sx q[0];
rz(-1.3277338) q[0];
sx q[0];
rz(-3.0603141) q[0];
rz(-0.47668138) q[1];
sx q[1];
rz(-1.3414693) q[1];
sx q[1];
rz(3.124602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.995249) q[0];
sx q[0];
rz(-2.0011847) q[0];
sx q[0];
rz(2.0018863) q[0];
rz(-2.4912253) q[2];
sx q[2];
rz(-1.128049) q[2];
sx q[2];
rz(0.75185094) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.022804755) q[1];
sx q[1];
rz(-0.14789109) q[1];
sx q[1];
rz(-2.1639362) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6829302) q[3];
sx q[3];
rz(-1.6067351) q[3];
sx q[3];
rz(3.1089715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9019258) q[2];
sx q[2];
rz(-1.5319752) q[2];
sx q[2];
rz(-1.7473183) q[2];
rz(1.4782921) q[3];
sx q[3];
rz(-0.52634382) q[3];
sx q[3];
rz(-0.31005508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33206853) q[0];
sx q[0];
rz(-0.14241756) q[0];
sx q[0];
rz(-2.944067) q[0];
rz(-1.4203513) q[1];
sx q[1];
rz(-2.0887801) q[1];
sx q[1];
rz(-1.24409) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3408139) q[0];
sx q[0];
rz(-1.7127303) q[0];
sx q[0];
rz(-2.5020788) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9483673) q[2];
sx q[2];
rz(-0.8329637) q[2];
sx q[2];
rz(1.1242614) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4287891) q[1];
sx q[1];
rz(-1.8658651) q[1];
sx q[1];
rz(-0.45810385) q[1];
rz(-pi) q[2];
rz(-2.1498419) q[3];
sx q[3];
rz(-1.5958991) q[3];
sx q[3];
rz(-1.2601528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5856005) q[2];
sx q[2];
rz(-2.8828794) q[2];
sx q[2];
rz(0.82536215) q[2];
rz(-2.182492) q[3];
sx q[3];
rz(-0.93950713) q[3];
sx q[3];
rz(-1.3033029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30967228) q[0];
sx q[0];
rz(-2.8490424) q[0];
sx q[0];
rz(1.9953315) q[0];
rz(-1.4347026) q[1];
sx q[1];
rz(-2.2271894) q[1];
sx q[1];
rz(0.31417876) q[1];
rz(-2.564128) q[2];
sx q[2];
rz(-0.86311917) q[2];
sx q[2];
rz(-1.6481332) q[2];
rz(1.0657749) q[3];
sx q[3];
rz(-1.4151971) q[3];
sx q[3];
rz(-0.74756037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
