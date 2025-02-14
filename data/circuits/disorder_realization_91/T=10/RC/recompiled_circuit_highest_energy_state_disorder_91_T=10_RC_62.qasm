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
rz(0.089601547) q[0];
sx q[0];
rz(-0.90553415) q[0];
sx q[0];
rz(-0.18984689) q[0];
rz(2.7449961) q[1];
sx q[1];
rz(-0.34062579) q[1];
sx q[1];
rz(2.5315898) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7977915) q[0];
sx q[0];
rz(-0.34489783) q[0];
sx q[0];
rz(2.3380322) q[0];
x q[1];
rz(1.1930614) q[2];
sx q[2];
rz(-0.47171041) q[2];
sx q[2];
rz(-0.38208252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0428517) q[1];
sx q[1];
rz(-1.0611649) q[1];
sx q[1];
rz(-2.3364728) q[1];
rz(-pi) q[2];
rz(-1.0759495) q[3];
sx q[3];
rz(-0.54022721) q[3];
sx q[3];
rz(-1.5264508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6923328) q[2];
sx q[2];
rz(-1.924694) q[2];
sx q[2];
rz(-0.24162351) q[2];
rz(0.062945098) q[3];
sx q[3];
rz(-1.2017622) q[3];
sx q[3];
rz(0.81301779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6908506) q[0];
sx q[0];
rz(-1.8618604) q[0];
sx q[0];
rz(-2.8513841) q[0];
rz(-0.33503512) q[1];
sx q[1];
rz(-0.93828833) q[1];
sx q[1];
rz(-2.8163574) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9262518) q[0];
sx q[0];
rz(-1.6115973) q[0];
sx q[0];
rz(-2.6759139) q[0];
rz(0.26518719) q[2];
sx q[2];
rz(-0.90861215) q[2];
sx q[2];
rz(-1.4410849) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8222386) q[1];
sx q[1];
rz(-1.7387973) q[1];
sx q[1];
rz(0.15138123) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65167221) q[3];
sx q[3];
rz(-1.0965986) q[3];
sx q[3];
rz(1.4693361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8364496) q[2];
sx q[2];
rz(-0.93905753) q[2];
sx q[2];
rz(-2.0024425) q[2];
rz(2.3818453) q[3];
sx q[3];
rz(-0.097948827) q[3];
sx q[3];
rz(-0.37794149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515795) q[0];
sx q[0];
rz(-2.4115998) q[0];
sx q[0];
rz(-2.3082025) q[0];
rz(3.1381798) q[1];
sx q[1];
rz(-0.5178057) q[1];
sx q[1];
rz(1.980967) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4118372) q[0];
sx q[0];
rz(-1.5416391) q[0];
sx q[0];
rz(1.2458318) q[0];
x q[1];
rz(2.4936475) q[2];
sx q[2];
rz(-1.3908236) q[2];
sx q[2];
rz(-0.56862105) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6033481) q[1];
sx q[1];
rz(-1.5711492) q[1];
sx q[1];
rz(1.805549) q[1];
x q[2];
rz(-2.0187601) q[3];
sx q[3];
rz(-2.7395757) q[3];
sx q[3];
rz(-2.1340919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87355906) q[2];
sx q[2];
rz(-1.7256836) q[2];
sx q[2];
rz(0.1669008) q[2];
rz(1.9514826) q[3];
sx q[3];
rz(-0.24470617) q[3];
sx q[3];
rz(1.083583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13389182) q[0];
sx q[0];
rz(-1.5950483) q[0];
sx q[0];
rz(-0.85064763) q[0];
rz(0.8029241) q[1];
sx q[1];
rz(-1.9839169) q[1];
sx q[1];
rz(2.0096774) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9659289) q[0];
sx q[0];
rz(-0.53232876) q[0];
sx q[0];
rz(1.2086297) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2766303) q[2];
sx q[2];
rz(-1.7554211) q[2];
sx q[2];
rz(-2.612584) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2537576) q[1];
sx q[1];
rz(-1.3030392) q[1];
sx q[1];
rz(-2.5303809) q[1];
rz(-pi) q[2];
rz(2.411667) q[3];
sx q[3];
rz(-1.2582964) q[3];
sx q[3];
rz(-0.70532986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5256646) q[2];
sx q[2];
rz(-1.5027081) q[2];
sx q[2];
rz(0.28044236) q[2];
rz(2.7403455) q[3];
sx q[3];
rz(-2.8420227) q[3];
sx q[3];
rz(1.7846599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2111557) q[0];
sx q[0];
rz(-1.7481952) q[0];
sx q[0];
rz(2.9537971) q[0];
rz(-0.64741778) q[1];
sx q[1];
rz(-2.1719666) q[1];
sx q[1];
rz(-2.5513249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85344106) q[0];
sx q[0];
rz(-2.7803505) q[0];
sx q[0];
rz(1.6608742) q[0];
rz(-1.9365385) q[2];
sx q[2];
rz(-2.3252333) q[2];
sx q[2];
rz(-2.2162645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94475585) q[1];
sx q[1];
rz(-1.4902672) q[1];
sx q[1];
rz(-1.5425986) q[1];
rz(-1.3594317) q[3];
sx q[3];
rz(-2.2244033) q[3];
sx q[3];
rz(-0.048649064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4831627) q[2];
sx q[2];
rz(-1.1700609) q[2];
sx q[2];
rz(2.9675193) q[2];
rz(2.6744794) q[3];
sx q[3];
rz(-2.4090448) q[3];
sx q[3];
rz(2.9113801) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1751404) q[0];
sx q[0];
rz(-1.7564961) q[0];
sx q[0];
rz(1.0873644) q[0];
rz(0.1611791) q[1];
sx q[1];
rz(-1.168707) q[1];
sx q[1];
rz(2.2081614) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1891834) q[0];
sx q[0];
rz(-1.8851452) q[0];
sx q[0];
rz(2.0670003) q[0];
rz(1.8699588) q[2];
sx q[2];
rz(-2.751902) q[2];
sx q[2];
rz(2.1914633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.63042917) q[1];
sx q[1];
rz(-1.1463209) q[1];
sx q[1];
rz(1.2066578) q[1];
rz(-0.94901086) q[3];
sx q[3];
rz(-1.0836667) q[3];
sx q[3];
rz(-2.0355527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33438385) q[2];
sx q[2];
rz(-1.4391359) q[2];
sx q[2];
rz(-1.2558233) q[2];
rz(3.1387591) q[3];
sx q[3];
rz(-0.56003672) q[3];
sx q[3];
rz(2.475256) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46347076) q[0];
sx q[0];
rz(-3.1246298) q[0];
sx q[0];
rz(2.8609138) q[0];
rz(-2.8374788) q[1];
sx q[1];
rz(-0.76144832) q[1];
sx q[1];
rz(1.6573409) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.723618) q[0];
sx q[0];
rz(-3.1153346) q[0];
sx q[0];
rz(2.2951207) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7353457) q[2];
sx q[2];
rz(-0.79716792) q[2];
sx q[2];
rz(0.88594243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9439895) q[1];
sx q[1];
rz(-1.0609975) q[1];
sx q[1];
rz(-2.5453525) q[1];
rz(3.0077107) q[3];
sx q[3];
rz(-1.1999793) q[3];
sx q[3];
rz(1.6755392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41875276) q[2];
sx q[2];
rz(-2.3374228) q[2];
sx q[2];
rz(2.1201102) q[2];
rz(-0.57182765) q[3];
sx q[3];
rz(-2.573206) q[3];
sx q[3];
rz(-0.9182601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21917139) q[0];
sx q[0];
rz(-0.89184856) q[0];
sx q[0];
rz(-0.047155596) q[0];
rz(-1.7258518) q[1];
sx q[1];
rz(-1.163131) q[1];
sx q[1];
rz(2.7498551) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3035256) q[0];
sx q[0];
rz(-0.45539344) q[0];
sx q[0];
rz(3.0440211) q[0];
rz(-pi) q[1];
rz(1.91003) q[2];
sx q[2];
rz(-1.7493141) q[2];
sx q[2];
rz(2.0925131) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.246576) q[1];
sx q[1];
rz(-1.5575174) q[1];
sx q[1];
rz(-3.1311656) q[1];
rz(-2.7487031) q[3];
sx q[3];
rz(-1.7709915) q[3];
sx q[3];
rz(1.4999215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91531104) q[2];
sx q[2];
rz(-2.2690124) q[2];
sx q[2];
rz(-0.20381168) q[2];
rz(-2.5053744) q[3];
sx q[3];
rz(-1.7812984) q[3];
sx q[3];
rz(-3.0740331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.5930475) q[0];
sx q[0];
rz(-0.021634463) q[0];
sx q[0];
rz(-2.0245323) q[0];
rz(-1.8695658) q[1];
sx q[1];
rz(-0.8465603) q[1];
sx q[1];
rz(0.55714947) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2761189) q[0];
sx q[0];
rz(-1.1412924) q[0];
sx q[0];
rz(0.19709023) q[0];
rz(-pi) q[1];
rz(-1.4049872) q[2];
sx q[2];
rz(-1.6631544) q[2];
sx q[2];
rz(-0.044008642) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3342569) q[1];
sx q[1];
rz(-0.42667056) q[1];
sx q[1];
rz(-0.92324275) q[1];
rz(-2.5999913) q[3];
sx q[3];
rz(-1.8159051) q[3];
sx q[3];
rz(-2.4211943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.40994) q[2];
sx q[2];
rz(-0.2356379) q[2];
sx q[2];
rz(-1.5059936) q[2];
rz(2.951494) q[3];
sx q[3];
rz(-2.5458769) q[3];
sx q[3];
rz(-3.0715517) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6079123) q[0];
sx q[0];
rz(-0.30617014) q[0];
sx q[0];
rz(-2.876907) q[0];
rz(-0.39868042) q[1];
sx q[1];
rz(-2.61187) q[1];
sx q[1];
rz(-2.9369489) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35025233) q[0];
sx q[0];
rz(-0.58316427) q[0];
sx q[0];
rz(1.2009317) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25252931) q[2];
sx q[2];
rz(-2.1871532) q[2];
sx q[2];
rz(-1.4184619) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61581332) q[1];
sx q[1];
rz(-0.74545292) q[1];
sx q[1];
rz(1.6604056) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7162303) q[3];
sx q[3];
rz(-1.0974636) q[3];
sx q[3];
rz(2.4504599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9201422) q[2];
sx q[2];
rz(-1.7998989) q[2];
sx q[2];
rz(-1.7500925) q[2];
rz(2.5642388) q[3];
sx q[3];
rz(-0.57830638) q[3];
sx q[3];
rz(-2.3277843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6089384) q[0];
sx q[0];
rz(-1.2464936) q[0];
sx q[0];
rz(-2.0410224) q[0];
rz(0.34250034) q[1];
sx q[1];
rz(-0.96440146) q[1];
sx q[1];
rz(-1.0920116) q[1];
rz(2.9663646) q[2];
sx q[2];
rz(-1.874085) q[2];
sx q[2];
rz(0.25599538) q[2];
rz(-0.92855056) q[3];
sx q[3];
rz(-0.83893574) q[3];
sx q[3];
rz(-0.75225603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
