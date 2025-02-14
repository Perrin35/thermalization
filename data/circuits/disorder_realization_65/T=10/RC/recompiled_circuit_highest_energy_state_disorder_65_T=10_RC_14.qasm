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
rz(0.49864545) q[0];
sx q[0];
rz(-2.5957624) q[0];
sx q[0];
rz(-0.11830615) q[0];
rz(-2.2220597) q[1];
sx q[1];
rz(-1.3032721) q[1];
sx q[1];
rz(-1.9272756) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8144381) q[0];
sx q[0];
rz(-2.0025578) q[0];
sx q[0];
rz(0.23349725) q[0];
x q[1];
rz(2.5164175) q[2];
sx q[2];
rz(-2.0707651) q[2];
sx q[2];
rz(1.878405) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8439629) q[1];
sx q[1];
rz(-1.7073892) q[1];
sx q[1];
rz(-0.94874391) q[1];
x q[2];
rz(-1.8439775) q[3];
sx q[3];
rz(-1.5915378) q[3];
sx q[3];
rz(-1.0047877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.97293568) q[2];
sx q[2];
rz(-2.0106222) q[2];
sx q[2];
rz(-1.0966148) q[2];
rz(-0.24122572) q[3];
sx q[3];
rz(-1.7719519) q[3];
sx q[3];
rz(2.296804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0519401) q[0];
sx q[0];
rz(-2.0515428) q[0];
sx q[0];
rz(0.3845149) q[0];
rz(0.3081201) q[1];
sx q[1];
rz(-0.48079753) q[1];
sx q[1];
rz(-0.62517977) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97323167) q[0];
sx q[0];
rz(-0.28795469) q[0];
sx q[0];
rz(-3.0298427) q[0];
x q[1];
rz(2.9956498) q[2];
sx q[2];
rz(-2.4193442) q[2];
sx q[2];
rz(-0.099522245) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5601666) q[1];
sx q[1];
rz(-2.4871965) q[1];
sx q[1];
rz(-3.0017716) q[1];
x q[2];
rz(-0.49021198) q[3];
sx q[3];
rz(-2.1313088) q[3];
sx q[3];
rz(1.1309159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1326617) q[2];
sx q[2];
rz(-2.2098358) q[2];
sx q[2];
rz(2.8399732) q[2];
rz(0.53145069) q[3];
sx q[3];
rz(-0.88574946) q[3];
sx q[3];
rz(-2.0187995) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6279491) q[0];
sx q[0];
rz(-1.0114089) q[0];
sx q[0];
rz(-1.9568141) q[0];
rz(2.9966677) q[1];
sx q[1];
rz(-1.8791608) q[1];
sx q[1];
rz(-1.5941934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3959409) q[0];
sx q[0];
rz(-0.81123039) q[0];
sx q[0];
rz(0.70180362) q[0];
x q[1];
rz(-2.9582439) q[2];
sx q[2];
rz(-1.1969222) q[2];
sx q[2];
rz(3.0484859) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.545274) q[1];
sx q[1];
rz(-0.68877586) q[1];
sx q[1];
rz(-0.94249814) q[1];
x q[2];
rz(2.4800042) q[3];
sx q[3];
rz(-2.6433655) q[3];
sx q[3];
rz(0.029904043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3323815) q[2];
sx q[2];
rz(-0.96497649) q[2];
sx q[2];
rz(-1.3261718) q[2];
rz(0.91340804) q[3];
sx q[3];
rz(-1.923442) q[3];
sx q[3];
rz(-0.53328812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.42295414) q[0];
sx q[0];
rz(-0.47684968) q[0];
sx q[0];
rz(-1.3612716) q[0];
rz(0.71296972) q[1];
sx q[1];
rz(-2.0013335) q[1];
sx q[1];
rz(2.4580809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28302971) q[0];
sx q[0];
rz(-0.57335317) q[0];
sx q[0];
rz(1.5893776) q[0];
x q[1];
rz(-0.2095378) q[2];
sx q[2];
rz(-1.7895797) q[2];
sx q[2];
rz(-1.0971341) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55895581) q[1];
sx q[1];
rz(-0.92933944) q[1];
sx q[1];
rz(-2.2333686) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38280757) q[3];
sx q[3];
rz(-0.26255739) q[3];
sx q[3];
rz(-2.5298339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5061364) q[2];
sx q[2];
rz(-2.7244302) q[2];
sx q[2];
rz(-0.48640856) q[2];
rz(-2.6247315) q[3];
sx q[3];
rz(-1.9603399) q[3];
sx q[3];
rz(-2.1080871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085623398) q[0];
sx q[0];
rz(-1.6317246) q[0];
sx q[0];
rz(-1.7686718) q[0];
rz(-1.7347451) q[1];
sx q[1];
rz(-0.56766784) q[1];
sx q[1];
rz(-0.47703201) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041056496) q[0];
sx q[0];
rz(-1.9874128) q[0];
sx q[0];
rz(-2.9265397) q[0];
rz(-pi) q[1];
rz(-2.7482175) q[2];
sx q[2];
rz(-1.0842536) q[2];
sx q[2];
rz(1.7560619) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9368767) q[1];
sx q[1];
rz(-0.5936247) q[1];
sx q[1];
rz(2.1267664) q[1];
rz(-pi) q[2];
rz(-0.78272696) q[3];
sx q[3];
rz(-1.0102135) q[3];
sx q[3];
rz(-2.4334139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10739747) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(-2.6743215) q[2];
rz(0.16921903) q[3];
sx q[3];
rz(-1.4110112) q[3];
sx q[3];
rz(2.8362078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53448236) q[0];
sx q[0];
rz(-2.2619673) q[0];
sx q[0];
rz(-0.18950263) q[0];
rz(-0.27319187) q[1];
sx q[1];
rz(-1.0739645) q[1];
sx q[1];
rz(-2.3409519) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36291262) q[0];
sx q[0];
rz(-1.8040468) q[0];
sx q[0];
rz(-0.044384758) q[0];
rz(-pi) q[1];
rz(2.8932299) q[2];
sx q[2];
rz(-1.6585322) q[2];
sx q[2];
rz(1.1540292) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2786497) q[1];
sx q[1];
rz(-0.23506308) q[1];
sx q[1];
rz(-2.1549015) q[1];
rz(2.6012035) q[3];
sx q[3];
rz(-1.1307934) q[3];
sx q[3];
rz(-0.42127452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2022986) q[2];
sx q[2];
rz(-0.17387986) q[2];
sx q[2];
rz(0.536971) q[2];
rz(-1.8592853) q[3];
sx q[3];
rz(-1.3064281) q[3];
sx q[3];
rz(-1.6622701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.110431) q[0];
sx q[0];
rz(-0.81542504) q[0];
sx q[0];
rz(1.331331) q[0];
rz(1.7000465) q[1];
sx q[1];
rz(-1.6621637) q[1];
sx q[1];
rz(-1.9190681) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5957907) q[0];
sx q[0];
rz(-1.8078601) q[0];
sx q[0];
rz(1.14181) q[0];
rz(-pi) q[1];
rz(2.5897854) q[2];
sx q[2];
rz(-1.6344317) q[2];
sx q[2];
rz(0.76630521) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8795506) q[1];
sx q[1];
rz(-0.73644222) q[1];
sx q[1];
rz(-1.8950736) q[1];
rz(-pi) q[2];
rz(1.8172713) q[3];
sx q[3];
rz(-2.685084) q[3];
sx q[3];
rz(0.64303165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9407201) q[2];
sx q[2];
rz(-0.49786374) q[2];
sx q[2];
rz(2.8087924) q[2];
rz(-2.8424272) q[3];
sx q[3];
rz(-1.9009512) q[3];
sx q[3];
rz(3.1201194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0126295) q[0];
sx q[0];
rz(-0.7681995) q[0];
sx q[0];
rz(-2.8666038) q[0];
rz(0.081427447) q[1];
sx q[1];
rz(-2.3402201) q[1];
sx q[1];
rz(-1.3224695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3429936) q[0];
sx q[0];
rz(-1.8019774) q[0];
sx q[0];
rz(-3.0804599) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1312348) q[2];
sx q[2];
rz(-2.4050674) q[2];
sx q[2];
rz(-2.9721952) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6024633) q[1];
sx q[1];
rz(-0.88383288) q[1];
sx q[1];
rz(-1.0844356) q[1];
rz(2.2127719) q[3];
sx q[3];
rz(-0.95212338) q[3];
sx q[3];
rz(0.88255461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.72121173) q[2];
sx q[2];
rz(-3.0215441) q[2];
sx q[2];
rz(-1.7397286) q[2];
rz(1.2841355) q[3];
sx q[3];
rz(-1.9522342) q[3];
sx q[3];
rz(0.0349667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5597124) q[0];
sx q[0];
rz(-1.5927097) q[0];
sx q[0];
rz(-2.6522719) q[0];
rz(3.1210461) q[1];
sx q[1];
rz(-1.9053883) q[1];
sx q[1];
rz(-0.70125854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72557455) q[0];
sx q[0];
rz(-1.3301992) q[0];
sx q[0];
rz(-2.1558236) q[0];
rz(-pi) q[1];
rz(-0.047766165) q[2];
sx q[2];
rz(-0.86776185) q[2];
sx q[2];
rz(-3.0813129) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.012246) q[1];
sx q[1];
rz(-1.80629) q[1];
sx q[1];
rz(-0.78462623) q[1];
x q[2];
rz(2.4550426) q[3];
sx q[3];
rz(-1.1995763) q[3];
sx q[3];
rz(2.8438501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7265085) q[2];
sx q[2];
rz(-2.5231611) q[2];
sx q[2];
rz(2.6611967) q[2];
rz(-1.4091617) q[3];
sx q[3];
rz(-1.3536072) q[3];
sx q[3];
rz(-2.8044146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82520634) q[0];
sx q[0];
rz(-1.9012863) q[0];
sx q[0];
rz(2.7301042) q[0];
rz(1.2443789) q[1];
sx q[1];
rz(-1.139541) q[1];
sx q[1];
rz(2.8172475) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56670226) q[0];
sx q[0];
rz(-1.5181203) q[0];
sx q[0];
rz(1.943669) q[0];
x q[1];
rz(1.2191992) q[2];
sx q[2];
rz(-1.8931539) q[2];
sx q[2];
rz(-0.11478648) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0592546) q[1];
sx q[1];
rz(-1.6538701) q[1];
sx q[1];
rz(-2.1359813) q[1];
x q[2];
rz(-3.0342372) q[3];
sx q[3];
rz(-0.12190652) q[3];
sx q[3];
rz(2.0481244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1856498) q[2];
sx q[2];
rz(-1.1142542) q[2];
sx q[2];
rz(3.0246217) q[2];
rz(2.1864435) q[3];
sx q[3];
rz(-0.17894608) q[3];
sx q[3];
rz(-2.5999033) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2627926) q[0];
sx q[0];
rz(-2.265082) q[0];
sx q[0];
rz(-2.2646917) q[0];
rz(-2.0211438) q[1];
sx q[1];
rz(-0.81605492) q[1];
sx q[1];
rz(1.1647404) q[1];
rz(1.0638857) q[2];
sx q[2];
rz(-1.3851266) q[2];
sx q[2];
rz(1.2141488) q[2];
rz(1.3001623) q[3];
sx q[3];
rz(-1.5951372) q[3];
sx q[3];
rz(-1.5230877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
