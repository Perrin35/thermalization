OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1144855) q[0];
sx q[0];
rz(-0.067582421) q[0];
sx q[0];
rz(-0.58900589) q[0];
rz(2.0677805) q[1];
sx q[1];
rz(-1.685073) q[1];
sx q[1];
rz(0.20751247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4366174) q[0];
sx q[0];
rz(-0.31148404) q[0];
sx q[0];
rz(-2.2732449) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6512538) q[2];
sx q[2];
rz(-1.7535889) q[2];
sx q[2];
rz(-1.8535943) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0879651) q[1];
sx q[1];
rz(-3.1131272) q[1];
sx q[1];
rz(-1.6971807) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21353586) q[3];
sx q[3];
rz(-2.6693404) q[3];
sx q[3];
rz(-2.6708598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8091858) q[2];
sx q[2];
rz(-3.1250592) q[2];
sx q[2];
rz(0.22745505) q[2];
rz(-0.65666947) q[3];
sx q[3];
rz(-0.69461099) q[3];
sx q[3];
rz(-2.606126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6856132) q[0];
sx q[0];
rz(-0.032129012) q[0];
sx q[0];
rz(-0.44422126) q[0];
rz(1.1086858) q[1];
sx q[1];
rz(-1.8007092) q[1];
sx q[1];
rz(1.9734372) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5801808) q[0];
sx q[0];
rz(-1.7636443) q[0];
sx q[0];
rz(0.69985244) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45380317) q[2];
sx q[2];
rz(-1.6112932) q[2];
sx q[2];
rz(-1.7292505) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25623577) q[1];
sx q[1];
rz(-2.0200122) q[1];
sx q[1];
rz(-0.083811772) q[1];
x q[2];
rz(-0.91932591) q[3];
sx q[3];
rz(-0.72685234) q[3];
sx q[3];
rz(-0.60692274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68219677) q[2];
sx q[2];
rz(-1.0707868) q[2];
sx q[2];
rz(0.11967858) q[2];
rz(0.97673544) q[3];
sx q[3];
rz(-0.027438199) q[3];
sx q[3];
rz(-1.1915092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.479849) q[0];
sx q[0];
rz(-0.16525826) q[0];
sx q[0];
rz(1.4953493) q[0];
rz(2.7961075) q[1];
sx q[1];
rz(-2.5472239) q[1];
sx q[1];
rz(2.3631309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6133967) q[0];
sx q[0];
rz(-1.5776618) q[0];
sx q[0];
rz(0.055533218) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3216281) q[2];
sx q[2];
rz(-1.1021309) q[2];
sx q[2];
rz(-2.6033273) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0960041) q[1];
sx q[1];
rz(-2.0109777) q[1];
sx q[1];
rz(0.11543302) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4069294) q[3];
sx q[3];
rz(-2.0941995) q[3];
sx q[3];
rz(-1.911834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8567132) q[2];
sx q[2];
rz(-0.036980193) q[2];
sx q[2];
rz(-1.7516837) q[2];
rz(3.1189647) q[3];
sx q[3];
rz(-0.32630625) q[3];
sx q[3];
rz(0.41445318) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49700272) q[0];
sx q[0];
rz(-0.039660064) q[0];
sx q[0];
rz(2.6139551) q[0];
rz(2.7854994) q[1];
sx q[1];
rz(-1.6368607) q[1];
sx q[1];
rz(1.5632695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8601807) q[0];
sx q[0];
rz(-1.1426272) q[0];
sx q[0];
rz(2.5606945) q[0];
x q[1];
rz(1.6354792) q[2];
sx q[2];
rz(-0.89406653) q[2];
sx q[2];
rz(-1.6874718) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.17828748) q[1];
sx q[1];
rz(-1.7366341) q[1];
sx q[1];
rz(2.01009) q[1];
rz(-pi) q[2];
rz(1.7308233) q[3];
sx q[3];
rz(-2.6607294) q[3];
sx q[3];
rz(2.4208925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.57033527) q[2];
sx q[2];
rz(-3.1331077) q[2];
sx q[2];
rz(-3.0961032) q[2];
rz(-0.82413834) q[3];
sx q[3];
rz(-2.6233311) q[3];
sx q[3];
rz(-2.1477108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5387886) q[0];
sx q[0];
rz(-2.0863918) q[0];
sx q[0];
rz(1.7428727) q[0];
rz(2.973373) q[1];
sx q[1];
rz(-1.8054447) q[1];
sx q[1];
rz(2.9199563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0960064) q[0];
sx q[0];
rz(-2.0915338) q[0];
sx q[0];
rz(-2.8686499) q[0];
rz(-pi) q[1];
rz(1.7847117) q[2];
sx q[2];
rz(-1.6375721) q[2];
sx q[2];
rz(-0.58148384) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4532104) q[1];
sx q[1];
rz(-0.4391138) q[1];
sx q[1];
rz(-1.6279787) q[1];
rz(1.0030909) q[3];
sx q[3];
rz(-2.4074887) q[3];
sx q[3];
rz(1.0455893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98442709) q[2];
sx q[2];
rz(-3.114279) q[2];
sx q[2];
rz(-1.0013162) q[2];
rz(-2.2973513) q[3];
sx q[3];
rz(-0.068035754) q[3];
sx q[3];
rz(0.74725738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7897414) q[0];
sx q[0];
rz(-2.8838367) q[0];
sx q[0];
rz(-1.3699654) q[0];
rz(0.17579707) q[1];
sx q[1];
rz(-1.5134696) q[1];
sx q[1];
rz(-2.083185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87230667) q[0];
sx q[0];
rz(-0.71505419) q[0];
sx q[0];
rz(-1.1319517) q[0];
x q[1];
rz(-1.5878229) q[2];
sx q[2];
rz(-2.5527708) q[2];
sx q[2];
rz(2.7463809) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0479916) q[1];
sx q[1];
rz(-0.8719043) q[1];
sx q[1];
rz(-1.2456139) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.109531) q[3];
sx q[3];
rz(-0.70343164) q[3];
sx q[3];
rz(-0.7348133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3093962) q[2];
sx q[2];
rz(-3.1377073) q[2];
sx q[2];
rz(-2.2957809) q[2];
rz(-0.63515615) q[3];
sx q[3];
rz(-0.35717765) q[3];
sx q[3];
rz(0.14491189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7307067) q[0];
sx q[0];
rz(-0.17179739) q[0];
sx q[0];
rz(-2.8806277) q[0];
rz(-1.7102526) q[1];
sx q[1];
rz(-0.15161082) q[1];
sx q[1];
rz(1.3014911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1430014) q[0];
sx q[0];
rz(-2.3218985) q[0];
sx q[0];
rz(-2.2065333) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13589047) q[2];
sx q[2];
rz(-3.0294578) q[2];
sx q[2];
rz(2.2035315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83205081) q[1];
sx q[1];
rz(-2.0027035) q[1];
sx q[1];
rz(-1.5276272) q[1];
rz(-pi) q[2];
rz(1.5742794) q[3];
sx q[3];
rz(-0.26393587) q[3];
sx q[3];
rz(1.6613632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5362376) q[2];
sx q[2];
rz(-1.321512) q[2];
sx q[2];
rz(-3.1129254) q[2];
rz(3.0981433) q[3];
sx q[3];
rz(-0.236792) q[3];
sx q[3];
rz(-1.4437599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.068977721) q[0];
sx q[0];
rz(-0.34729877) q[0];
sx q[0];
rz(0.29796991) q[0];
rz(-1.4120253) q[1];
sx q[1];
rz(-2.7822918) q[1];
sx q[1];
rz(1.5922458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2320579) q[0];
sx q[0];
rz(-1.2611054) q[0];
sx q[0];
rz(-2.880611) q[0];
rz(-1.554073) q[2];
sx q[2];
rz(-1.564476) q[2];
sx q[2];
rz(0.14594742) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9142278) q[1];
sx q[1];
rz(-2.4691205) q[1];
sx q[1];
rz(2.6537031) q[1];
rz(0.91101908) q[3];
sx q[3];
rz(-0.35189128) q[3];
sx q[3];
rz(-1.4705197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90184244) q[2];
sx q[2];
rz(-0.032278927) q[2];
sx q[2];
rz(0.65995222) q[2];
rz(-2.3615071) q[3];
sx q[3];
rz(-1.8476906) q[3];
sx q[3];
rz(1.9735146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6338947) q[0];
sx q[0];
rz(-3.1051023) q[0];
sx q[0];
rz(-0.80398917) q[0];
rz(2.9102303) q[1];
sx q[1];
rz(-1.7295126) q[1];
sx q[1];
rz(3.0661327) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.696414) q[0];
sx q[0];
rz(-1.7761765) q[0];
sx q[0];
rz(2.7778047) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.6663975e-05) q[2];
sx q[2];
rz(-1.5701887) q[2];
sx q[2];
rz(3.031784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4659404) q[1];
sx q[1];
rz(-0.91788252) q[1];
sx q[1];
rz(-0.5799579) q[1];
rz(-pi) q[2];
rz(-1.4867501) q[3];
sx q[3];
rz(-2.7527697) q[3];
sx q[3];
rz(-3.0879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.017645322) q[2];
sx q[2];
rz(-2.6804774) q[2];
sx q[2];
rz(0.41575113) q[2];
rz(1.6425284) q[3];
sx q[3];
rz(-3.1406904) q[3];
sx q[3];
rz(-2.3568962) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7163664) q[0];
sx q[0];
rz(-3.0854736) q[0];
sx q[0];
rz(2.8477493) q[0];
rz(1.5360688) q[1];
sx q[1];
rz(-0.87296456) q[1];
sx q[1];
rz(-1.4651089) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11742453) q[0];
sx q[0];
rz(-0.10483042) q[0];
sx q[0];
rz(0.094314055) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60428166) q[2];
sx q[2];
rz(-1.2245032) q[2];
sx q[2];
rz(1.1187579) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.067071384) q[1];
sx q[1];
rz(-1.5980729) q[1];
sx q[1];
rz(2.4949353) q[1];
rz(2.5876643) q[3];
sx q[3];
rz(-0.0084453765) q[3];
sx q[3];
rz(-2.3609207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6426223) q[2];
sx q[2];
rz(-2.5140646) q[2];
sx q[2];
rz(1.9210531) q[2];
rz(2.7484861) q[3];
sx q[3];
rz(-0.016963907) q[3];
sx q[3];
rz(0.19256798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4569693) q[0];
sx q[0];
rz(-1.4263117) q[0];
sx q[0];
rz(-0.34761467) q[0];
rz(-2.5034703) q[1];
sx q[1];
rz(-0.77824021) q[1];
sx q[1];
rz(2.4660769) q[1];
rz(-0.38370321) q[2];
sx q[2];
rz(-0.20846882) q[2];
sx q[2];
rz(-0.61412703) q[2];
rz(-1.5777231) q[3];
sx q[3];
rz(-1.4580115) q[3];
sx q[3];
rz(1.9409822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
