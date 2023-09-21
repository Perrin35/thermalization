OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4857793) q[0];
sx q[0];
rz(-1.9137148) q[0];
sx q[0];
rz(-1.2113843) q[0];
rz(-0.70648944) q[2];
sx q[2];
rz(-0.90105614) q[2];
sx q[2];
rz(1.1342088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5174487) q[1];
sx q[1];
rz(-1.2341208) q[1];
sx q[1];
rz(1.8271853) q[1];
rz(-2.0516112) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(-2.4570176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2661665) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(2.0092633) q[2];
rz(1.6752361) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(1.0124538) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448626) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-2.9557513) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(2.9247608) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29138716) q[0];
sx q[0];
rz(-0.7190401) q[0];
sx q[0];
rz(1.1262116) q[0];
rz(-pi) q[1];
rz(-1.0160604) q[2];
sx q[2];
rz(-1.9811355) q[2];
sx q[2];
rz(1.0085269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5237907) q[1];
sx q[1];
rz(-2.4317867) q[1];
sx q[1];
rz(1.8989423) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84083765) q[3];
sx q[3];
rz(-0.84078046) q[3];
sx q[3];
rz(0.18883146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.310114) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(-1.8537834) q[2];
rz(-2.3790322) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-2.8365703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6771616) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(2.537354) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(0.93260971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9661449) q[0];
sx q[0];
rz(-1.517059) q[0];
sx q[0];
rz(-1.8885814) q[0];
rz(2.2432125) q[2];
sx q[2];
rz(-0.2873688) q[2];
sx q[2];
rz(2.3786366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.697726) q[1];
sx q[1];
rz(-1.6892471) q[1];
sx q[1];
rz(-2.5812134) q[1];
rz(-pi) q[2];
rz(-0.44585769) q[3];
sx q[3];
rz(-2.0816396) q[3];
sx q[3];
rz(2.0542991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.147826) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(2.0489342) q[2];
rz(0.5422194) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(-2.6413667) q[0];
rz(-2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(1.4979699) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26101199) q[0];
sx q[0];
rz(-1.8646761) q[0];
sx q[0];
rz(1.0512933) q[0];
rz(-1.285032) q[2];
sx q[2];
rz(-2.9432202) q[2];
sx q[2];
rz(2.6464268) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3264309) q[1];
sx q[1];
rz(-1.3109428) q[1];
sx q[1];
rz(-1.3304779) q[1];
rz(-pi) q[2];
rz(-0.76969947) q[3];
sx q[3];
rz(-0.61004988) q[3];
sx q[3];
rz(-1.1806012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(0.70181075) q[2];
rz(-2.3102405) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9005301) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(2.3262614) q[0];
rz(1.6197846) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(1.048208) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3807555) q[0];
sx q[0];
rz(-1.8341944) q[0];
sx q[0];
rz(0.25816985) q[0];
rz(-pi) q[1];
rz(-1.5585209) q[2];
sx q[2];
rz(-0.91931146) q[2];
sx q[2];
rz(2.2001681) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2506927) q[1];
sx q[1];
rz(-0.39784583) q[1];
sx q[1];
rz(-0.65244168) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0190373) q[3];
sx q[3];
rz(-0.69283797) q[3];
sx q[3];
rz(-2.2321731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6158225) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(-2.0416416) q[2];
rz(0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(-2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0734171) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(0.90240479) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(3.0117603) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79486217) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(2.5559588) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1075222) q[2];
sx q[2];
rz(-2.495129) q[2];
sx q[2];
rz(-1.549364) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1631158) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(-0.58660581) q[1];
rz(-0.31452175) q[3];
sx q[3];
rz(-0.57146996) q[3];
sx q[3];
rz(0.86626569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(-0.20425805) q[2];
rz(-1.9355109) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234574) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.4468505) q[0];
rz(-1.2591259) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(-0.68626219) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2080363) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(-0.74525381) q[0];
x q[1];
rz(-2.1967728) q[2];
sx q[2];
rz(-2.2959024) q[2];
sx q[2];
rz(0.041989728) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.64993091) q[1];
sx q[1];
rz(-1.6451391) q[1];
sx q[1];
rz(1.7564303) q[1];
rz(-pi) q[2];
rz(0.42240123) q[3];
sx q[3];
rz(-1.8654612) q[3];
sx q[3];
rz(2.9871123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69616047) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(-0.0017722842) q[2];
rz(0.56162515) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6034265) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(-1.4461393) q[0];
rz(-2.360545) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(1.6400281) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7088889) q[0];
sx q[0];
rz(-0.49312691) q[0];
sx q[0];
rz(-1.6119484) q[0];
rz(-3.115032) q[2];
sx q[2];
rz(-1.5090669) q[2];
sx q[2];
rz(0.82690566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0280684) q[1];
sx q[1];
rz(-1.4816195) q[1];
sx q[1];
rz(0.32797565) q[1];
rz(-pi) q[2];
rz(2.1498508) q[3];
sx q[3];
rz(-1.0240882) q[3];
sx q[3];
rz(-2.9255097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(1.8224576) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8050352) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.2040899) q[0];
rz(2.7583292) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4597804) q[0];
sx q[0];
rz(-0.60428719) q[0];
sx q[0];
rz(-2.2629645) q[0];
rz(-pi) q[1];
rz(0.69182379) q[2];
sx q[2];
rz(-1.8080538) q[2];
sx q[2];
rz(2.0123864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.44605276) q[1];
sx q[1];
rz(-2.3791168) q[1];
sx q[1];
rz(-1.6474849) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6978108) q[3];
sx q[3];
rz(-0.1212596) q[3];
sx q[3];
rz(-1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(-1.8593672) q[2];
rz(1.4964237) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6431817) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(-2.9472651) q[0];
rz(-2.1037897) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(1.0338354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4170096) q[0];
sx q[0];
rz(-1.4405182) q[0];
sx q[0];
rz(-2.2307322) q[0];
x q[1];
rz(0.98722234) q[2];
sx q[2];
rz(-0.7910896) q[2];
sx q[2];
rz(-0.9466048) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0018113) q[1];
sx q[1];
rz(-1.0011295) q[1];
sx q[1];
rz(-3.0210178) q[1];
x q[2];
rz(0.71318993) q[3];
sx q[3];
rz(-1.7367559) q[3];
sx q[3];
rz(2.1470269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(0.6357843) q[2];
rz(0.27030269) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(-1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.4355961) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(0.031899115) q[2];
sx q[2];
rz(-0.96822856) q[2];
sx q[2];
rz(-0.4005489) q[2];
rz(-0.070449645) q[3];
sx q[3];
rz(-1.2415213) q[3];
sx q[3];
rz(0.51125676) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];