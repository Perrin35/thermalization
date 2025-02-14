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
rz(-7.3569975) q[1];
sx q[1];
rz(4.8266657) q[1];
sx q[1];
rz(6.0756728) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5975345) q[0];
sx q[0];
rz(-1.7701214) q[0];
sx q[0];
rz(-1.3298278) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18337266) q[2];
sx q[2];
rz(-1.6499106) q[2];
sx q[2];
rz(-0.26814207) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0536276) q[1];
sx q[1];
rz(-0.028465406) q[1];
sx q[1];
rz(1.444412) q[1];
rz(-pi) q[2];
rz(1.4629685) q[3];
sx q[3];
rz(-2.031481) q[3];
sx q[3];
rz(0.70957795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8091858) q[2];
sx q[2];
rz(-3.1250592) q[2];
sx q[2];
rz(2.9141376) q[2];
rz(2.4849232) q[3];
sx q[3];
rz(-0.69461099) q[3];
sx q[3];
rz(0.53546661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4559795) q[0];
sx q[0];
rz(-3.1094636) q[0];
sx q[0];
rz(0.44422126) q[0];
rz(2.0329068) q[1];
sx q[1];
rz(-1.3408835) q[1];
sx q[1];
rz(1.9734372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21441169) q[0];
sx q[0];
rz(-0.72158883) q[0];
sx q[0];
rz(-2.847228) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5257448) q[2];
sx q[2];
rz(-2.0241996) q[2];
sx q[2];
rz(-0.13870961) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6942648) q[1];
sx q[1];
rz(-0.45644293) q[1];
sx q[1];
rz(-1.7427299) q[1];
x q[2];
rz(0.95529629) q[3];
sx q[3];
rz(-1.1560734) q[3];
sx q[3];
rz(-0.44594582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4593959) q[2];
sx q[2];
rz(-2.0708059) q[2];
sx q[2];
rz(0.11967858) q[2];
rz(2.1648572) q[3];
sx q[3];
rz(-3.1141545) q[3];
sx q[3];
rz(1.9500835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6617436) q[0];
sx q[0];
rz(-2.9763344) q[0];
sx q[0];
rz(-1.6462434) q[0];
rz(-2.7961075) q[1];
sx q[1];
rz(-2.5472239) q[1];
sx q[1];
rz(-2.3631309) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0986106) q[0];
sx q[0];
rz(-1.5152644) q[0];
sx q[0];
rz(1.5639202) q[0];
rz(-pi) q[1];
rz(0.81996452) q[2];
sx q[2];
rz(-1.1021309) q[2];
sx q[2];
rz(0.53826538) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.22009199) q[1];
sx q[1];
rz(-0.45410546) q[1];
sx q[1];
rz(-1.3309671) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71078474) q[3];
sx q[3];
rz(-0.87276283) q[3];
sx q[3];
rz(-0.16434114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8567132) q[2];
sx q[2];
rz(-3.1046125) q[2];
sx q[2];
rz(1.3899089) q[2];
rz(0.022627929) q[3];
sx q[3];
rz(-0.32630625) q[3];
sx q[3];
rz(2.7271395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(2.6445899) q[0];
sx q[0];
rz(-3.1019326) q[0];
sx q[0];
rz(0.52763754) q[0];
rz(-2.7854994) q[1];
sx q[1];
rz(-1.6368607) q[1];
sx q[1];
rz(-1.5632695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1183076) q[0];
sx q[0];
rz(-1.0480799) q[0];
sx q[0];
rz(2.0705332) q[0];
rz(-pi) q[1];
rz(3.061297) q[2];
sx q[2];
rz(-0.67932898) q[2];
sx q[2];
rz(-1.7905362) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8265022) q[1];
sx q[1];
rz(-1.13794) q[1];
sx q[1];
rz(-2.9587246) q[1];
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
rz(0.99388188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5387886) q[0];
sx q[0];
rz(-2.0863918) q[0];
sx q[0];
rz(1.3987199) q[0];
rz(-0.16821965) q[1];
sx q[1];
rz(-1.3361479) q[1];
sx q[1];
rz(0.22163637) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7547627) q[0];
sx q[0];
rz(-1.3347751) q[0];
sx q[0];
rz(2.1079662) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7847117) q[2];
sx q[2];
rz(-1.6375721) q[2];
sx q[2];
rz(0.58148384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.065818345) q[1];
sx q[1];
rz(-1.5950959) q[1];
sx q[1];
rz(-2.009281) q[1];
rz(2.6898674) q[3];
sx q[3];
rz(-2.171031) q[3];
sx q[3];
rz(-0.33590318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1571656) q[2];
sx q[2];
rz(-0.027313622) q[2];
sx q[2];
rz(1.0013162) q[2];
rz(-0.84424132) q[3];
sx q[3];
rz(-0.068035754) q[3];
sx q[3];
rz(-0.74725738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35185128) q[0];
sx q[0];
rz(-2.8838367) q[0];
sx q[0];
rz(1.3699654) q[0];
rz(0.17579707) q[1];
sx q[1];
rz(-1.5134696) q[1];
sx q[1];
rz(1.0584077) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0390801) q[0];
sx q[0];
rz(-1.2884757) q[0];
sx q[0];
rz(-2.2369871) q[0];
x q[1];
rz(1.5878229) q[2];
sx q[2];
rz(-2.5527708) q[2];
sx q[2];
rz(0.39521171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0479916) q[1];
sx q[1];
rz(-2.2696884) q[1];
sx q[1];
rz(-1.2456139) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2203683) q[3];
sx q[3];
rz(-1.8628253) q[3];
sx q[3];
rz(-2.667922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83219641) q[2];
sx q[2];
rz(-3.1377073) q[2];
sx q[2];
rz(-0.84581172) q[2];
rz(2.5064365) q[3];
sx q[3];
rz(-0.35717765) q[3];
sx q[3];
rz(0.14491189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7307067) q[0];
sx q[0];
rz(-0.17179739) q[0];
sx q[0];
rz(2.8806277) q[0];
rz(1.7102526) q[1];
sx q[1];
rz(-2.9899818) q[1];
sx q[1];
rz(-1.8401015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9985913) q[0];
sx q[0];
rz(-0.81969417) q[0];
sx q[0];
rz(-2.2065333) q[0];
rz(-pi) q[1];
rz(-1.5555423) q[2];
sx q[2];
rz(-1.4596995) q[2];
sx q[2];
rz(2.3402702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72922414) q[1];
sx q[1];
rz(-0.43392402) q[1];
sx q[1];
rz(3.0482376) q[1];
rz(-pi) q[2];
rz(1.5673133) q[3];
sx q[3];
rz(-2.8776568) q[3];
sx q[3];
rz(1.6613632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5362376) q[2];
sx q[2];
rz(-1.8200807) q[2];
sx q[2];
rz(-0.028667299) q[2];
rz(3.0981433) q[3];
sx q[3];
rz(-0.236792) q[3];
sx q[3];
rz(1.6978327) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0726149) q[0];
sx q[0];
rz(-2.7942939) q[0];
sx q[0];
rz(0.29796991) q[0];
rz(1.4120253) q[1];
sx q[1];
rz(-0.35930082) q[1];
sx q[1];
rz(1.5922458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41995364) q[0];
sx q[0];
rz(-1.8190939) q[0];
sx q[0];
rz(1.8906276) q[0];
rz(-1.554073) q[2];
sx q[2];
rz(-1.5771167) q[2];
sx q[2];
rz(2.9956452) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1916532) q[1];
sx q[1];
rz(-1.8671163) q[1];
sx q[1];
rz(0.61298989) q[1];
rz(-1.8531591) q[3];
sx q[3];
rz(-1.783665) q[3];
sx q[3];
rz(2.4119056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2397502) q[2];
sx q[2];
rz(-0.032278927) q[2];
sx q[2];
rz(-2.4816404) q[2];
rz(-0.7800855) q[3];
sx q[3];
rz(-1.8476906) q[3];
sx q[3];
rz(-1.9735146) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6338947) q[0];
sx q[0];
rz(-3.1051023) q[0];
sx q[0];
rz(2.3376035) q[0];
rz(-0.2313624) q[1];
sx q[1];
rz(-1.4120801) q[1];
sx q[1];
rz(-3.0661327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20310878) q[0];
sx q[0];
rz(-1.214998) q[0];
sx q[0];
rz(-1.3514765) q[0];
rz(-pi) q[1];
rz(-1.5701887) q[2];
sx q[2];
rz(-1.5707597) q[2];
sx q[2];
rz(-1.680605) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.14743155) q[1];
sx q[1];
rz(-2.2975337) q[1];
sx q[1];
rz(0.94908157) q[1];
rz(-1.9583804) q[3];
sx q[3];
rz(-1.5389666) q[3];
sx q[3];
rz(-1.4393161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1239473) q[2];
sx q[2];
rz(-0.46111527) q[2];
sx q[2];
rz(-0.41575113) q[2];
rz(1.6425284) q[3];
sx q[3];
rz(-0.00090229546) q[3];
sx q[3];
rz(2.3568962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(2.8477493) q[0];
rz(1.5360688) q[1];
sx q[1];
rz(-0.87296456) q[1];
sx q[1];
rz(1.6764838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241681) q[0];
sx q[0];
rz(-3.0367622) q[0];
sx q[0];
rz(-0.094314055) q[0];
x q[1];
rz(1.15756) q[2];
sx q[2];
rz(-1.0069478) q[2];
sx q[2];
rz(2.4593632) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.067071384) q[1];
sx q[1];
rz(-1.5980729) q[1];
sx q[1];
rz(2.4949353) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5876643) q[3];
sx q[3];
rz(-3.1331473) q[3];
sx q[3];
rz(2.3609207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6426223) q[2];
sx q[2];
rz(-2.5140646) q[2];
sx q[2];
rz(-1.2205396) q[2];
rz(0.39310655) q[3];
sx q[3];
rz(-0.016963907) q[3];
sx q[3];
rz(2.9490247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6846234) q[0];
sx q[0];
rz(-1.715281) q[0];
sx q[0];
rz(2.793978) q[0];
rz(-2.5034703) q[1];
sx q[1];
rz(-0.77824021) q[1];
sx q[1];
rz(2.4660769) q[1];
rz(2.9478922) q[2];
sx q[2];
rz(-1.4932409) q[2];
sx q[2];
rz(1.3328339) q[2];
rz(-3.0288052) q[3];
sx q[3];
rz(-1.5639136) q[3];
sx q[3];
rz(-2.7706272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
