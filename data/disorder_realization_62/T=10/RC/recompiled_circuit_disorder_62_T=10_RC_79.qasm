OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(-2.9641889) q[0];
sx q[0];
rz(2.0071964) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(-2.477975) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952607) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(2.6803826) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77779777) q[2];
sx q[2];
rz(-1.3133089) q[2];
sx q[2];
rz(2.4553026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9259778) q[1];
sx q[1];
rz(-1.8324592) q[1];
sx q[1];
rz(1.91933) q[1];
rz(-0.10856467) q[3];
sx q[3];
rz(-1.3285471) q[3];
sx q[3];
rz(3.0755175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.87876451) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(-0.051068548) q[2];
rz(-0.55705327) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(-1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5550845) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(-0.82582981) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(-1.9155496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2873043) q[0];
sx q[0];
rz(-1.1278296) q[0];
sx q[0];
rz(1.3501549) q[0];
x q[1];
rz(0.77483564) q[2];
sx q[2];
rz(-2.2171387) q[2];
sx q[2];
rz(-1.6006084) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5920168) q[1];
sx q[1];
rz(-1.5712275) q[1];
sx q[1];
rz(-1.3509343) q[1];
rz(-pi) q[2];
rz(0.0094718178) q[3];
sx q[3];
rz(-0.99820271) q[3];
sx q[3];
rz(0.026281683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(0.33102316) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(-1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(1.249041) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(-1.0294611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80699608) q[0];
sx q[0];
rz(-0.96141978) q[0];
sx q[0];
rz(2.8732804) q[0];
x q[1];
rz(1.2842032) q[2];
sx q[2];
rz(-2.2006052) q[2];
sx q[2];
rz(2.0558002) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0975115) q[1];
sx q[1];
rz(-1.71002) q[1];
sx q[1];
rz(2.984725) q[1];
x q[2];
rz(0.73996468) q[3];
sx q[3];
rz(-1.4810586) q[3];
sx q[3];
rz(-1.2147853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.133698) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(-0.50764817) q[2];
rz(-1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(-1.5456276) q[0];
rz(-2.0987299) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(1.5159336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0632616) q[0];
sx q[0];
rz(-0.69561361) q[0];
sx q[0];
rz(0.22380933) q[0];
rz(3.0558673) q[2];
sx q[2];
rz(-2.1282196) q[2];
sx q[2];
rz(-0.9312219) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.792946) q[1];
sx q[1];
rz(-2.209084) q[1];
sx q[1];
rz(-1.2299041) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5528615) q[3];
sx q[3];
rz(-1.3874467) q[3];
sx q[3];
rz(2.6453032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3778014) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(-0.14136782) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448733) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(1.5198583) q[0];
rz(0.63201085) q[1];
sx q[1];
rz(-0.72223392) q[1];
sx q[1];
rz(-2.246726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1631854) q[0];
sx q[0];
rz(-0.95046959) q[0];
sx q[0];
rz(-1.0582256) q[0];
rz(-pi) q[1];
rz(0.72563719) q[2];
sx q[2];
rz(-1.1302395) q[2];
sx q[2];
rz(-2.5051136) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4836854) q[1];
sx q[1];
rz(-0.92805082) q[1];
sx q[1];
rz(2.583858) q[1];
x q[2];
rz(-2.8551293) q[3];
sx q[3];
rz(-1.5781919) q[3];
sx q[3];
rz(1.8252107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.616509) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(1.1425225) q[2];
rz(-0.74674314) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557945) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.3775795) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(2.6766052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12596345) q[0];
sx q[0];
rz(-2.0032126) q[0];
sx q[0];
rz(0.36884357) q[0];
rz(-pi) q[1];
rz(-1.497252) q[2];
sx q[2];
rz(-1.2300756) q[2];
sx q[2];
rz(1.7433012) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2496693) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(-0.93530099) q[1];
rz(-pi) q[2];
rz(1.3607849) q[3];
sx q[3];
rz(-0.87267733) q[3];
sx q[3];
rz(2.7426646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49089367) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(0.4450376) q[2];
rz(0.93368357) q[3];
sx q[3];
rz(-1.7088339) q[3];
sx q[3];
rz(0.26708189) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12748195) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(-1.7215464) q[0];
rz(3.1177915) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(0.15596095) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0322745) q[0];
sx q[0];
rz(-0.30800691) q[0];
sx q[0];
rz(2.5449635) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6290226) q[2];
sx q[2];
rz(-1.9480431) q[2];
sx q[2];
rz(0.23114983) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8487726) q[1];
sx q[1];
rz(-0.18310586) q[1];
sx q[1];
rz(-1.2700901) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3994) q[3];
sx q[3];
rz(-0.67897292) q[3];
sx q[3];
rz(-1.2364482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.069313958) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496562) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(0.33777133) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(0.24857323) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08911207) q[0];
sx q[0];
rz(-2.5626474) q[0];
sx q[0];
rz(-2.6738033) q[0];
x q[1];
rz(-2.6612501) q[2];
sx q[2];
rz(-2.0311653) q[2];
sx q[2];
rz(0.72887052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0174745) q[1];
sx q[1];
rz(-1.0867456) q[1];
sx q[1];
rz(-2.2494621) q[1];
x q[2];
rz(1.400984) q[3];
sx q[3];
rz(-1.6465934) q[3];
sx q[3];
rz(-0.71789391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2934072) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(-0.48197204) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(-1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329353) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(-2.8588262) q[0];
rz(0.70156082) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(-1.823002) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5134207) q[0];
sx q[0];
rz(-1.4903729) q[0];
sx q[0];
rz(1.1286939) q[0];
rz(-pi) q[1];
rz(0.96005) q[2];
sx q[2];
rz(-2.1337482) q[2];
sx q[2];
rz(2.8658531) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5732167) q[1];
sx q[1];
rz(-1.425256) q[1];
sx q[1];
rz(2.2578866) q[1];
rz(2.7301844) q[3];
sx q[3];
rz(-0.79380006) q[3];
sx q[3];
rz(1.0994764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7302154) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-0.39917699) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(-1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76319641) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(-1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(-1.261196) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6529918) q[0];
sx q[0];
rz(-1.6761259) q[0];
sx q[0];
rz(-0.41098849) q[0];
rz(1.5485498) q[2];
sx q[2];
rz(-2.1689479) q[2];
sx q[2];
rz(-1.0215789) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0295769) q[1];
sx q[1];
rz(-2.638991) q[1];
sx q[1];
rz(-0.13346787) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4576549) q[3];
sx q[3];
rz(-1.6947019) q[3];
sx q[3];
rz(2.370196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(-1.2184881) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8284843) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(0.60824153) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(1.6115887) q[2];
sx q[2];
rz(-0.57208021) q[2];
sx q[2];
rz(-0.013442599) q[2];
rz(2.9363587) q[3];
sx q[3];
rz(-0.60094613) q[3];
sx q[3];
rz(-0.080106674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];