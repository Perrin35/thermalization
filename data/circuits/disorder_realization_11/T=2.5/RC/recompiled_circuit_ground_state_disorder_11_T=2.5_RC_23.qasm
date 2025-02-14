OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0271072) q[0];
sx q[0];
rz(-3.0740102) q[0];
sx q[0];
rz(-2.5525868) q[0];
rz(2.0677805) q[1];
sx q[1];
rz(-1.685073) q[1];
sx q[1];
rz(-2.9340802) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70497528) q[0];
sx q[0];
rz(-2.8301086) q[0];
sx q[0];
rz(-2.2732449) q[0];
x q[1];
rz(-0.41012058) q[2];
sx q[2];
rz(-2.9420576) q[2];
sx q[2];
rz(-1.4360957) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0536276) q[1];
sx q[1];
rz(-3.1131272) q[1];
sx q[1];
rz(1.444412) q[1];
x q[2];
rz(-1.6786241) q[3];
sx q[3];
rz(-2.031481) q[3];
sx q[3];
rz(0.70957795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8091858) q[2];
sx q[2];
rz(-0.016533479) q[2];
sx q[2];
rz(-0.22745505) q[2];
rz(-0.65666947) q[3];
sx q[3];
rz(-2.4469817) q[3];
sx q[3];
rz(-0.53546661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.4559795) q[0];
sx q[0];
rz(-3.1094636) q[0];
sx q[0];
rz(-0.44422126) q[0];
rz(-1.1086858) q[1];
sx q[1];
rz(-1.8007092) q[1];
sx q[1];
rz(-1.9734372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16938528) q[0];
sx q[0];
rz(-0.88645259) q[0];
sx q[0];
rz(1.3208525) q[0];
rz(0.092165784) q[2];
sx q[2];
rz(-2.6861114) q[2];
sx q[2];
rz(-0.24126894) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2780952) q[1];
sx q[1];
rz(-1.6462763) q[1];
sx q[1];
rz(-1.1202034) q[1];
x q[2];
rz(-0.95529629) q[3];
sx q[3];
rz(-1.1560734) q[3];
sx q[3];
rz(-2.6956468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4593959) q[2];
sx q[2];
rz(-2.0708059) q[2];
sx q[2];
rz(-0.11967858) q[2];
rz(0.97673544) q[3];
sx q[3];
rz(-0.027438199) q[3];
sx q[3];
rz(-1.1915092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6617436) q[0];
sx q[0];
rz(-2.9763344) q[0];
sx q[0];
rz(-1.4953493) q[0];
rz(0.34548512) q[1];
sx q[1];
rz(-0.59436878) q[1];
sx q[1];
rz(-0.77846175) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.528196) q[0];
sx q[0];
rz(-1.5776618) q[0];
sx q[0];
rz(-0.055533218) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81996452) q[2];
sx q[2];
rz(-2.0394618) q[2];
sx q[2];
rz(-0.53826538) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9215007) q[1];
sx q[1];
rz(-2.6874872) q[1];
sx q[1];
rz(1.8106255) q[1];
rz(-pi) q[2];
rz(2.4308079) q[3];
sx q[3];
rz(-0.87276283) q[3];
sx q[3];
rz(2.9772515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8567132) q[2];
sx q[2];
rz(-0.036980193) q[2];
sx q[2];
rz(1.7516837) q[2];
rz(0.022627929) q[3];
sx q[3];
rz(-0.32630625) q[3];
sx q[3];
rz(2.7271395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6445899) q[0];
sx q[0];
rz(-0.039660064) q[0];
sx q[0];
rz(-0.52763754) q[0];
rz(2.7854994) q[1];
sx q[1];
rz(-1.6368607) q[1];
sx q[1];
rz(-1.5783232) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023285063) q[0];
sx q[0];
rz(-2.0935128) q[0];
sx q[0];
rz(1.0710595) q[0];
x q[1];
rz(3.061297) q[2];
sx q[2];
rz(-0.67932898) q[2];
sx q[2];
rz(1.3510564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4112579) q[1];
sx q[1];
rz(-0.46763845) q[1];
sx q[1];
rz(1.1958666) q[1];
rz(-pi) q[2];
rz(1.0951871) q[3];
sx q[3];
rz(-1.6445673) q[3];
sx q[3];
rz(2.1493585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5712574) q[2];
sx q[2];
rz(-3.1331077) q[2];
sx q[2];
rz(3.0961032) q[2];
rz(0.82413834) q[3];
sx q[3];
rz(-2.6233311) q[3];
sx q[3];
rz(2.1477108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5387886) q[0];
sx q[0];
rz(-2.0863918) q[0];
sx q[0];
rz(-1.3987199) q[0];
rz(0.16821965) q[1];
sx q[1];
rz(-1.8054447) q[1];
sx q[1];
rz(0.22163637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0960064) q[0];
sx q[0];
rz(-2.0915338) q[0];
sx q[0];
rz(2.8686499) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.068328502) q[2];
sx q[2];
rz(-1.7842275) q[2];
sx q[2];
rz(-2.1377856) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.065818345) q[1];
sx q[1];
rz(-1.5464968) q[1];
sx q[1];
rz(-2.009281) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6898674) q[3];
sx q[3];
rz(-2.171031) q[3];
sx q[3];
rz(2.8056895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98442709) q[2];
sx q[2];
rz(-0.027313622) q[2];
sx q[2];
rz(1.0013162) q[2];
rz(2.2973513) q[3];
sx q[3];
rz(-3.0735569) q[3];
sx q[3];
rz(0.74725738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7897414) q[0];
sx q[0];
rz(-2.8838367) q[0];
sx q[0];
rz(1.3699654) q[0];
rz(-2.9657956) q[1];
sx q[1];
rz(-1.5134696) q[1];
sx q[1];
rz(-2.083185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1025126) q[0];
sx q[0];
rz(-1.853117) q[0];
sx q[0];
rz(2.2369871) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5878229) q[2];
sx q[2];
rz(-0.58882182) q[2];
sx q[2];
rz(0.39521171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0936011) q[1];
sx q[1];
rz(-0.8719043) q[1];
sx q[1];
rz(1.2456139) q[1];
x q[2];
rz(-2.0320616) q[3];
sx q[3];
rz(-2.438161) q[3];
sx q[3];
rz(2.4067794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3093962) q[2];
sx q[2];
rz(-0.0038853566) q[2];
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
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(-1.8401015) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8230391) q[0];
sx q[0];
rz(-2.1995499) q[0];
sx q[0];
rz(0.56644337) q[0];
rz(-3.030483) q[2];
sx q[2];
rz(-1.5859563) q[2];
sx q[2];
rz(0.76778256) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4123685) q[1];
sx q[1];
rz(-0.43392402) q[1];
sx q[1];
rz(3.0482376) q[1];
x q[2];
rz(-1.306862) q[3];
sx q[3];
rz(-1.571705) q[3];
sx q[3];
rz(3.0543882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6053551) q[2];
sx q[2];
rz(-1.8200807) q[2];
sx q[2];
rz(3.1129254) q[2];
rz(0.043449314) q[3];
sx q[3];
rz(-2.9048007) q[3];
sx q[3];
rz(-1.4437599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068977721) q[0];
sx q[0];
rz(-2.7942939) q[0];
sx q[0];
rz(2.8436227) q[0];
rz(-1.4120253) q[1];
sx q[1];
rz(-2.7822918) q[1];
sx q[1];
rz(-1.5493468) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2320579) q[0];
sx q[0];
rz(-1.2611054) q[0];
sx q[0];
rz(-0.26098164) q[0];
rz(-pi) q[1];
rz(0.0063212379) q[2];
sx q[2];
rz(-1.5875193) q[2];
sx q[2];
rz(1.716638) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3181646) q[1];
sx q[1];
rz(-2.1534502) q[1];
sx q[1];
rz(1.9280487) q[1];
rz(-0.91101908) q[3];
sx q[3];
rz(-0.35189128) q[3];
sx q[3];
rz(1.4705197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2397502) q[2];
sx q[2];
rz(-3.1093137) q[2];
sx q[2];
rz(0.65995222) q[2];
rz(2.3615071) q[3];
sx q[3];
rz(-1.2939021) q[3];
sx q[3];
rz(1.9735146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50769794) q[0];
sx q[0];
rz(-3.1051023) q[0];
sx q[0];
rz(-2.3376035) q[0];
rz(-2.9102303) q[1];
sx q[1];
rz(-1.4120801) q[1];
sx q[1];
rz(-0.075459935) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4451786) q[0];
sx q[0];
rz(-1.7761765) q[0];
sx q[0];
rz(0.36378797) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.141556) q[2];
sx q[2];
rz(-1.571404) q[2];
sx q[2];
rz(3.031784) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8676843) q[1];
sx q[1];
rz(-1.1204506) q[1];
sx q[1];
rz(-0.83012786) q[1];
rz(-1.4867501) q[3];
sx q[3];
rz(-2.7527697) q[3];
sx q[3];
rz(-3.0879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.017645322) q[2];
sx q[2];
rz(-2.6804774) q[2];
sx q[2];
rz(0.41575113) q[2];
rz(1.6425284) q[3];
sx q[3];
rz(-3.1406904) q[3];
sx q[3];
rz(0.7846964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7163664) q[0];
sx q[0];
rz(-3.0854736) q[0];
sx q[0];
rz(-0.29384336) q[0];
rz(1.6055239) q[1];
sx q[1];
rz(-0.87296456) q[1];
sx q[1];
rz(-1.6764838) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595725) q[0];
sx q[0];
rz(-1.5806507) q[0];
sx q[0];
rz(-3.0372247) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9840327) q[2];
sx q[2];
rz(-1.0069478) q[2];
sx q[2];
rz(-0.68222943) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6017573) q[1];
sx q[1];
rz(-0.64714995) q[1];
sx q[1];
rz(3.0963417) q[1];
x q[2];
rz(1.5752389) q[3];
sx q[3];
rz(-1.5636139) q[3];
sx q[3];
rz(-1.8069763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6426223) q[2];
sx q[2];
rz(-2.5140646) q[2];
sx q[2];
rz(-1.9210531) q[2];
rz(2.7484861) q[3];
sx q[3];
rz(-0.016963907) q[3];
sx q[3];
rz(0.19256798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-2.7578894) q[2];
sx q[2];
rz(-2.9331238) q[2];
sx q[2];
rz(2.5274656) q[2];
rz(0.11278747) q[3];
sx q[3];
rz(-1.5639136) q[3];
sx q[3];
rz(-2.7706272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
