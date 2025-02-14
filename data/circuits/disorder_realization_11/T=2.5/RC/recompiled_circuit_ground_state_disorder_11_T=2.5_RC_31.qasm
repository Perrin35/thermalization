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
rz(3.0740102) q[0];
sx q[0];
rz(10.013784) q[0];
rz(2.0677805) q[1];
sx q[1];
rz(-1.685073) q[1];
sx q[1];
rz(0.20751247) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02188259) q[0];
sx q[0];
rz(-1.8069022) q[0];
sx q[0];
rz(0.2050928) q[0];
rz(-0.41012058) q[2];
sx q[2];
rz(-0.1995351) q[2];
sx q[2];
rz(-1.705497) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3564975) q[1];
sx q[1];
rz(-1.5743839) q[1];
sx q[1];
rz(-1.5990348) q[1];
rz(2.6785844) q[3];
sx q[3];
rz(-1.6673458) q[3];
sx q[3];
rz(0.90930401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33240685) q[2];
sx q[2];
rz(-3.1250592) q[2];
sx q[2];
rz(-0.22745505) q[2];
rz(0.65666947) q[3];
sx q[3];
rz(-0.69461099) q[3];
sx q[3];
rz(2.606126) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6856132) q[0];
sx q[0];
rz(-3.1094636) q[0];
sx q[0];
rz(-2.6973714) q[0];
rz(-2.0329068) q[1];
sx q[1];
rz(-1.8007092) q[1];
sx q[1];
rz(1.9734372) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9722074) q[0];
sx q[0];
rz(-2.2551401) q[0];
sx q[0];
rz(-1.3208525) q[0];
x q[1];
rz(1.6158478) q[2];
sx q[2];
rz(-1.117393) q[2];
sx q[2];
rz(3.002883) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6942648) q[1];
sx q[1];
rz(-2.6851497) q[1];
sx q[1];
rz(1.7427299) q[1];
x q[2];
rz(-0.91932591) q[3];
sx q[3];
rz(-2.4147403) q[3];
sx q[3];
rz(-2.5346699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4593959) q[2];
sx q[2];
rz(-1.0707868) q[2];
sx q[2];
rz(-0.11967858) q[2];
rz(2.1648572) q[3];
sx q[3];
rz(-3.1141545) q[3];
sx q[3];
rz(1.9500835) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.479849) q[0];
sx q[0];
rz(-0.16525826) q[0];
sx q[0];
rz(1.6462434) q[0];
rz(2.7961075) q[1];
sx q[1];
rz(-2.5472239) q[1];
sx q[1];
rz(-0.77846175) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6133967) q[0];
sx q[0];
rz(-1.5639308) q[0];
sx q[0];
rz(-3.0860594) q[0];
x q[1];
rz(2.3216281) q[2];
sx q[2];
rz(-2.0394618) q[2];
sx q[2];
rz(-2.6033273) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0960041) q[1];
sx q[1];
rz(-2.0109777) q[1];
sx q[1];
rz(0.11543302) q[1];
rz(-pi) q[2];
rz(-2.2317844) q[3];
sx q[3];
rz(-2.1902553) q[3];
sx q[3];
rz(0.76515686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8567132) q[2];
sx q[2];
rz(-3.1046125) q[2];
sx q[2];
rz(1.7516837) q[2];
rz(-3.1189647) q[3];
sx q[3];
rz(-0.32630625) q[3];
sx q[3];
rz(-0.41445318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6445899) q[0];
sx q[0];
rz(-3.1019326) q[0];
sx q[0];
rz(2.6139551) q[0];
rz(-0.35609326) q[1];
sx q[1];
rz(-1.5047319) q[1];
sx q[1];
rz(-1.5632695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85335818) q[0];
sx q[0];
rz(-2.4348867) q[0];
sx q[0];
rz(2.447829) q[0];
x q[1];
rz(-1.5061134) q[2];
sx q[2];
rz(-2.2475261) q[2];
sx q[2];
rz(-1.4541208) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4112579) q[1];
sx q[1];
rz(-2.6739542) q[1];
sx q[1];
rz(-1.1958666) q[1];
rz(-3.0586518) q[3];
sx q[3];
rz(-1.0965875) q[3];
sx q[3];
rz(-0.54061962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5712574) q[2];
sx q[2];
rz(-0.0084849914) q[2];
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
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5387886) q[0];
sx q[0];
rz(-2.0863918) q[0];
sx q[0];
rz(-1.7428727) q[0];
rz(-2.973373) q[1];
sx q[1];
rz(-1.8054447) q[1];
sx q[1];
rz(-2.9199563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5835041) q[0];
sx q[0];
rz(-2.5595491) q[0];
sx q[0];
rz(-2.01016) q[0];
rz(-pi) q[1];
rz(-1.356881) q[2];
sx q[2];
rz(-1.5040205) q[2];
sx q[2];
rz(-2.5601088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.065818345) q[1];
sx q[1];
rz(-1.5950959) q[1];
sx q[1];
rz(2.009281) q[1];
x q[2];
rz(1.0030909) q[3];
sx q[3];
rz(-2.4074887) q[3];
sx q[3];
rz(1.0455893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1571656) q[2];
sx q[2];
rz(-3.114279) q[2];
sx q[2];
rz(2.1402764) q[2];
rz(0.84424132) q[3];
sx q[3];
rz(-0.068035754) q[3];
sx q[3];
rz(0.74725738) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35185128) q[0];
sx q[0];
rz(-2.8838367) q[0];
sx q[0];
rz(-1.7716273) q[0];
rz(2.9657956) q[1];
sx q[1];
rz(-1.628123) q[1];
sx q[1];
rz(-2.083185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1025126) q[0];
sx q[0];
rz(-1.2884757) q[0];
sx q[0];
rz(-2.2369871) q[0];
x q[1];
rz(-2.1595512) q[2];
sx q[2];
rz(-1.5802523) q[2];
sx q[2];
rz(-1.9801677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5762944) q[1];
sx q[1];
rz(-2.382462) q[1];
sx q[1];
rz(0.36328333) q[1];
x q[2];
rz(-0.92122434) q[3];
sx q[3];
rz(-1.8628253) q[3];
sx q[3];
rz(-0.47367063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83219641) q[2];
sx q[2];
rz(-0.0038853566) q[2];
sx q[2];
rz(-2.2957809) q[2];
rz(-2.5064365) q[3];
sx q[3];
rz(-0.35717765) q[3];
sx q[3];
rz(2.9966808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7307067) q[0];
sx q[0];
rz(-2.9697953) q[0];
sx q[0];
rz(-2.8806277) q[0];
rz(-1.7102526) q[1];
sx q[1];
rz(-2.9899818) q[1];
sx q[1];
rz(-1.3014911) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8230391) q[0];
sx q[0];
rz(-2.1995499) q[0];
sx q[0];
rz(-0.56644337) q[0];
x q[1];
rz(-3.030483) q[2];
sx q[2];
rz(-1.5556364) q[2];
sx q[2];
rz(-0.76778256) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72922414) q[1];
sx q[1];
rz(-0.43392402) q[1];
sx q[1];
rz(-3.0482376) q[1];
x q[2];
rz(1.306862) q[3];
sx q[3];
rz(-1.571705) q[3];
sx q[3];
rz(0.087204427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6053551) q[2];
sx q[2];
rz(-1.8200807) q[2];
sx q[2];
rz(3.1129254) q[2];
rz(3.0981433) q[3];
sx q[3];
rz(-2.9048007) q[3];
sx q[3];
rz(1.4437599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068977721) q[0];
sx q[0];
rz(-2.7942939) q[0];
sx q[0];
rz(-2.8436227) q[0];
rz(1.4120253) q[1];
sx q[1];
rz(-2.7822918) q[1];
sx q[1];
rz(-1.5922458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2320579) q[0];
sx q[0];
rz(-1.8804872) q[0];
sx q[0];
rz(0.26098164) q[0];
rz(-pi) q[1];
rz(1.5875196) q[2];
sx q[2];
rz(-1.5771167) q[2];
sx q[2];
rz(-0.14594742) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1916532) q[1];
sx q[1];
rz(-1.2744764) q[1];
sx q[1];
rz(-2.5286028) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91101908) q[3];
sx q[3];
rz(-2.7897014) q[3];
sx q[3];
rz(-1.4705197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90184244) q[2];
sx q[2];
rz(-3.1093137) q[2];
sx q[2];
rz(-0.65995222) q[2];
rz(-2.3615071) q[3];
sx q[3];
rz(-1.8476906) q[3];
sx q[3];
rz(-1.1680781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50769794) q[0];
sx q[0];
rz(-0.036490353) q[0];
sx q[0];
rz(2.3376035) q[0];
rz(-0.2313624) q[1];
sx q[1];
rz(-1.7295126) q[1];
sx q[1];
rz(-0.075459935) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7754525) q[0];
sx q[0];
rz(-0.41549993) q[0];
sx q[0];
rz(0.52966161) q[0];
rz(1.5105336) q[2];
sx q[2];
rz(-0.00060877006) q[2];
sx q[2];
rz(-3.0920467) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4659404) q[1];
sx q[1];
rz(-0.91788252) q[1];
sx q[1];
rz(-0.5799579) q[1];
rz(-pi) q[2];
rz(1.1832123) q[3];
sx q[3];
rz(-1.5389666) q[3];
sx q[3];
rz(-1.4393161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.017645322) q[2];
sx q[2];
rz(-0.46111527) q[2];
sx q[2];
rz(0.41575113) q[2];
rz(1.6425284) q[3];
sx q[3];
rz(-0.00090229546) q[3];
sx q[3];
rz(2.3568962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7163664) q[0];
sx q[0];
rz(-0.056119053) q[0];
sx q[0];
rz(-2.8477493) q[0];
rz(1.6055239) q[1];
sx q[1];
rz(-2.2686281) q[1];
sx q[1];
rz(1.6764838) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241681) q[0];
sx q[0];
rz(-0.10483042) q[0];
sx q[0];
rz(-0.094314055) q[0];
rz(-pi) q[1];
rz(-0.60428166) q[2];
sx q[2];
rz(-1.9170894) q[2];
sx q[2];
rz(1.1187579) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6017573) q[1];
sx q[1];
rz(-2.4944427) q[1];
sx q[1];
rz(0.045250968) q[1];
rz(0.0071825414) q[3];
sx q[3];
rz(-1.5752388) q[3];
sx q[3];
rz(2.9053807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6426223) q[2];
sx q[2];
rz(-0.62752807) q[2];
sx q[2];
rz(1.2205396) q[2];
rz(-0.39310655) q[3];
sx q[3];
rz(-0.016963907) q[3];
sx q[3];
rz(-2.9490247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4569693) q[0];
sx q[0];
rz(-1.715281) q[0];
sx q[0];
rz(2.793978) q[0];
rz(2.5034703) q[1];
sx q[1];
rz(-2.3633524) q[1];
sx q[1];
rz(-0.6755158) q[1];
rz(-2.9478922) q[2];
sx q[2];
rz(-1.6483518) q[2];
sx q[2];
rz(-1.8087587) q[2];
rz(0.06107851) q[3];
sx q[3];
rz(-0.11299639) q[3];
sx q[3];
rz(-1.1391409) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
