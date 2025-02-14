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
rz(1.7776547) q[0];
sx q[0];
rz(-0.41843709) q[0];
sx q[0];
rz(-1.3016181) q[0];
rz(0.16297451) q[1];
sx q[1];
rz(-1.4459223) q[1];
sx q[1];
rz(0.016782848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.421464) q[0];
sx q[0];
rz(-1.4961494) q[0];
sx q[0];
rz(1.2089085) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7503337) q[2];
sx q[2];
rz(-1.5870567) q[2];
sx q[2];
rz(3.1237512) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.0019579828) q[1];
sx q[1];
rz(-0.0070925709) q[1];
sx q[1];
rz(-2.3435739) q[1];
x q[2];
rz(-0.15233202) q[3];
sx q[3];
rz(-1.6248584) q[3];
sx q[3];
rz(-0.033933725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5590543) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(-1.5622697) q[2];
rz(-2.9301379) q[3];
sx q[3];
rz(-0.00051694218) q[3];
sx q[3];
rz(0.16099425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6120537) q[0];
sx q[0];
rz(-0.27612975) q[0];
sx q[0];
rz(1.8147234) q[0];
rz(0.57948411) q[1];
sx q[1];
rz(-3.1377628) q[1];
sx q[1];
rz(-0.63900596) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20637437) q[0];
sx q[0];
rz(-2.2274096) q[0];
sx q[0];
rz(-2.4103863) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0206356) q[2];
sx q[2];
rz(-1.5874344) q[2];
sx q[2];
rz(-1.5543092) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.745805) q[1];
sx q[1];
rz(-1.5823872) q[1];
sx q[1];
rz(0.017000217) q[1];
rz(-pi) q[2];
rz(-3.0538179) q[3];
sx q[3];
rz(-2.3592161) q[3];
sx q[3];
rz(-2.6915472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4352033) q[2];
sx q[2];
rz(-3.0053164) q[2];
sx q[2];
rz(-1.603568) q[2];
rz(1.5609353) q[3];
sx q[3];
rz(-0.014336421) q[3];
sx q[3];
rz(-3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2494217) q[0];
sx q[0];
rz(-2.6264661) q[0];
sx q[0];
rz(0.38145915) q[0];
rz(-2.4341266) q[1];
sx q[1];
rz(-0.019376945) q[1];
sx q[1];
rz(-1.1245419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8302001) q[0];
sx q[0];
rz(-1.8233577) q[0];
sx q[0];
rz(-2.8930032) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6857621) q[2];
sx q[2];
rz(-1.5966151) q[2];
sx q[2];
rz(-1.7242277) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8360236) q[1];
sx q[1];
rz(-1.5828805) q[1];
sx q[1];
rz(-3.0786242) q[1];
x q[2];
rz(2.400983) q[3];
sx q[3];
rz(-0.69878529) q[3];
sx q[3];
rz(1.6388338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6996998) q[2];
sx q[2];
rz(-3.1293588) q[2];
sx q[2];
rz(-3.0921248) q[2];
rz(-0.60702819) q[3];
sx q[3];
rz(-0.0012461239) q[3];
sx q[3];
rz(-1.2073257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(0.98168755) q[0];
sx q[0];
rz(-0.1683546) q[0];
sx q[0];
rz(-3.1244151) q[0];
rz(2.848564) q[1];
sx q[1];
rz(-2.3511062) q[1];
sx q[1];
rz(-1.5944098) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3348986) q[0];
sx q[0];
rz(-2.3133754) q[0];
sx q[0];
rz(-0.81636274) q[0];
rz(0.36357673) q[2];
sx q[2];
rz(-2.5695317) q[2];
sx q[2];
rz(-0.043302082) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80611011) q[1];
sx q[1];
rz(-1.4467738) q[1];
sx q[1];
rz(-1.5107442) q[1];
rz(1.5610351) q[3];
sx q[3];
rz(-0.92484821) q[3];
sx q[3];
rz(-2.498359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8881417) q[2];
sx q[2];
rz(-0.47051045) q[2];
sx q[2];
rz(-0.68224254) q[2];
rz(0.055179723) q[3];
sx q[3];
rz(-0.0076871593) q[3];
sx q[3];
rz(1.8365708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3799915) q[0];
sx q[0];
rz(-2.70607) q[0];
sx q[0];
rz(0.4250266) q[0];
rz(1.60166) q[1];
sx q[1];
rz(-0.48307499) q[1];
sx q[1];
rz(2.3262598) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5402712) q[0];
sx q[0];
rz(-3.0794332) q[0];
sx q[0];
rz(-0.17610733) q[0];
x q[1];
rz(-0.0032665423) q[2];
sx q[2];
rz(-1.5940394) q[2];
sx q[2];
rz(2.2540664) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22742352) q[1];
sx q[1];
rz(-2.9939751) q[1];
sx q[1];
rz(-0.62995894) q[1];
x q[2];
rz(-1.6922705) q[3];
sx q[3];
rz(-1.2290579) q[3];
sx q[3];
rz(-0.69167826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0067979) q[2];
sx q[2];
rz(-3.1289913) q[2];
sx q[2];
rz(1.4726144) q[2];
rz(2.2115479) q[3];
sx q[3];
rz(-3.127122) q[3];
sx q[3];
rz(0.84021935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3227661) q[0];
sx q[0];
rz(-0.023094026) q[0];
sx q[0];
rz(-1.7148788) q[0];
rz(-2.4115883) q[1];
sx q[1];
rz(-2.5529824) q[1];
sx q[1];
rz(2.0402562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0520605) q[0];
sx q[0];
rz(-0.87498481) q[0];
sx q[0];
rz(-2.7874283) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9726102) q[2];
sx q[2];
rz(-1.797953) q[2];
sx q[2];
rz(-0.81534602) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8423398) q[1];
sx q[1];
rz(-1.4758759) q[1];
sx q[1];
rz(-1.4466982) q[1];
rz(-pi) q[2];
rz(3.1093842) q[3];
sx q[3];
rz(-1.695249) q[3];
sx q[3];
rz(-2.2197753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4250028) q[2];
sx q[2];
rz(-0.060881946) q[2];
sx q[2];
rz(1.3072183) q[2];
rz(0.37846765) q[3];
sx q[3];
rz(-3.1186447) q[3];
sx q[3];
rz(0.65346658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3017479) q[0];
sx q[0];
rz(-1.2675588) q[0];
sx q[0];
rz(0.92754716) q[0];
rz(-1.7842267) q[1];
sx q[1];
rz(-0.83186847) q[1];
sx q[1];
rz(1.6102788) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5132234) q[0];
sx q[0];
rz(-1.5931061) q[0];
sx q[0];
rz(1.6091899) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3258097) q[2];
sx q[2];
rz(-1.9568482) q[2];
sx q[2];
rz(0.99967848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5802637) q[1];
sx q[1];
rz(-3.0124843) q[1];
sx q[1];
rz(1.5498954) q[1];
rz(-0.31503265) q[3];
sx q[3];
rz(-1.0803534) q[3];
sx q[3];
rz(-0.88446188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7808468) q[2];
sx q[2];
rz(-3.1372034) q[2];
sx q[2];
rz(1.8878262) q[2];
rz(-2.4474261) q[3];
sx q[3];
rz(-2.4024051) q[3];
sx q[3];
rz(-0.23301253) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6041782) q[0];
sx q[0];
rz(-1.0056714) q[0];
sx q[0];
rz(2.1042714) q[0];
rz(1.6088156) q[1];
sx q[1];
rz(-2.9210126) q[1];
sx q[1];
rz(1.6745837) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98437418) q[0];
sx q[0];
rz(-1.7658002) q[0];
sx q[0];
rz(0.068679811) q[0];
x q[1];
rz(-0.27811082) q[2];
sx q[2];
rz(-1.5928972) q[2];
sx q[2];
rz(1.8430361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3008219) q[1];
sx q[1];
rz(-1.5696708) q[1];
sx q[1];
rz(-0.00043934396) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19123938) q[3];
sx q[3];
rz(-2.6789224) q[3];
sx q[3];
rz(-1.3706051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1303225) q[2];
sx q[2];
rz(-0.20877561) q[2];
sx q[2];
rz(-3.0976963) q[2];
rz(-0.52055001) q[3];
sx q[3];
rz(-3.136941) q[3];
sx q[3];
rz(1.4588149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875824) q[0];
sx q[0];
rz(-3.1391322) q[0];
sx q[0];
rz(-1.3186697) q[0];
rz(-1.7240546) q[1];
sx q[1];
rz(-2.8520165) q[1];
sx q[1];
rz(-1.5444548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.330724) q[0];
sx q[0];
rz(-1.4557342) q[0];
sx q[0];
rz(-0.48145357) q[0];
x q[1];
rz(-0.6575281) q[2];
sx q[2];
rz(-1.7818461) q[2];
sx q[2];
rz(-1.5653277) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.456157) q[1];
sx q[1];
rz(-2.7815869) q[1];
sx q[1];
rz(-2.092157) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4539976) q[3];
sx q[3];
rz(-1.15117) q[3];
sx q[3];
rz(2.9586723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80457193) q[2];
sx q[2];
rz(-1.8420409) q[2];
sx q[2];
rz(2.9447832) q[2];
rz(-1.9491516) q[3];
sx q[3];
rz(-0.20659031) q[3];
sx q[3];
rz(0.20096745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30109677) q[0];
sx q[0];
rz(-1.3571285) q[0];
sx q[0];
rz(-1.1916196) q[0];
rz(-1.6169029) q[1];
sx q[1];
rz(-0.646851) q[1];
sx q[1];
rz(1.5764538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7708262) q[0];
sx q[0];
rz(-1.2773371) q[0];
sx q[0];
rz(-1.7459041) q[0];
rz(-pi) q[1];
rz(1.1249816) q[2];
sx q[2];
rz(-2.6458394) q[2];
sx q[2];
rz(3.1204566) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31101878) q[1];
sx q[1];
rz(-1.569848) q[1];
sx q[1];
rz(1.5712954) q[1];
rz(-pi) q[2];
rz(-3.0281046) q[3];
sx q[3];
rz(-1.4295414) q[3];
sx q[3];
rz(-0.065441386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9578751) q[2];
sx q[2];
rz(-2.5536733) q[2];
sx q[2];
rz(-1.4584165) q[2];
rz(0.030473907) q[3];
sx q[3];
rz(-3.1320429) q[3];
sx q[3];
rz(-2.9392346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4771582) q[0];
sx q[0];
rz(-1.3344593) q[0];
sx q[0];
rz(1.6819171) q[0];
rz(-1.5741813) q[1];
sx q[1];
rz(-1.3290783) q[1];
sx q[1];
rz(-3.0507416) q[1];
rz(-1.6363999) q[2];
sx q[2];
rz(-3.0658683) q[2];
sx q[2];
rz(0.22584596) q[2];
rz(2.1021765) q[3];
sx q[3];
rz(-2.3134091) q[3];
sx q[3];
rz(-0.29534657) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
