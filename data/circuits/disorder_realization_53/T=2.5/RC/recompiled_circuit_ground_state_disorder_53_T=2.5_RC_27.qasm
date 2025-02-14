OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98874918) q[0];
sx q[0];
rz(-0.56892836) q[0];
sx q[0];
rz(-0.95175728) q[0];
rz(0.8079575) q[1];
sx q[1];
rz(-3.0264049) q[1];
sx q[1];
rz(0.99339956) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0885926) q[0];
sx q[0];
rz(-1.3198084) q[0];
sx q[0];
rz(-0.95337501) q[0];
rz(-2.1615218) q[2];
sx q[2];
rz(-0.23071846) q[2];
sx q[2];
rz(1.7144817) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6019878) q[1];
sx q[1];
rz(-1.0430416) q[1];
sx q[1];
rz(-0.62202203) q[1];
x q[2];
rz(-0.66178721) q[3];
sx q[3];
rz(-1.4628304) q[3];
sx q[3];
rz(1.3707096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6318165) q[2];
sx q[2];
rz(-2.5413373) q[2];
sx q[2];
rz(0.14524761) q[2];
rz(-0.97233573) q[3];
sx q[3];
rz(-1.626868) q[3];
sx q[3];
rz(1.8632896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.4165322) q[0];
sx q[0];
rz(-0.57045492) q[0];
sx q[0];
rz(-0.78713) q[0];
rz(-2.9540673) q[1];
sx q[1];
rz(-1.7555321) q[1];
sx q[1];
rz(-3.131391) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30287179) q[0];
sx q[0];
rz(-0.10940675) q[0];
sx q[0];
rz(1.6463237) q[0];
rz(2.0554916) q[2];
sx q[2];
rz(-1.5124413) q[2];
sx q[2];
rz(1.397152) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.71828498) q[1];
sx q[1];
rz(-1.3415636) q[1];
sx q[1];
rz(-1.4450816) q[1];
x q[2];
rz(-0.83008978) q[3];
sx q[3];
rz(-1.4674868) q[3];
sx q[3];
rz(1.938397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0804704) q[2];
sx q[2];
rz(-3.0219813) q[2];
sx q[2];
rz(-1.1629026) q[2];
rz(1.2969147) q[3];
sx q[3];
rz(-1.4328522) q[3];
sx q[3];
rz(-3.1356623) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6892683) q[0];
sx q[0];
rz(-2.5732714) q[0];
sx q[0];
rz(-2.7962621) q[0];
rz(-0.055222424) q[1];
sx q[1];
rz(-0.60631141) q[1];
sx q[1];
rz(0.20923722) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.077186) q[0];
sx q[0];
rz(-1.4923054) q[0];
sx q[0];
rz(-2.2839943) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.018304869) q[2];
sx q[2];
rz(-1.6430815) q[2];
sx q[2];
rz(0.56906453) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1112695) q[1];
sx q[1];
rz(-1.7220338) q[1];
sx q[1];
rz(2.4200685) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8992527) q[3];
sx q[3];
rz(-2.3018357) q[3];
sx q[3];
rz(3.1104607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0914586) q[2];
sx q[2];
rz(-1.4860934) q[2];
sx q[2];
rz(2.6996108) q[2];
rz(2.6637391) q[3];
sx q[3];
rz(-1.7171532) q[3];
sx q[3];
rz(1.8817687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.0100937) q[0];
sx q[0];
rz(-2.5690014) q[0];
sx q[0];
rz(-1.9189438) q[0];
rz(-1.6652416) q[1];
sx q[1];
rz(-1.0226701) q[1];
sx q[1];
rz(2.8388035) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091508301) q[0];
sx q[0];
rz(-1.2909527) q[0];
sx q[0];
rz(-2.742663) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6499182) q[2];
sx q[2];
rz(-1.6058358) q[2];
sx q[2];
rz(-1.6043881) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.11377442) q[1];
sx q[1];
rz(-2.824265) q[1];
sx q[1];
rz(-1.7185663) q[1];
rz(-pi) q[2];
rz(0.62742558) q[3];
sx q[3];
rz(-1.0194821) q[3];
sx q[3];
rz(2.3105636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0786232) q[2];
sx q[2];
rz(-1.4457694) q[2];
sx q[2];
rz(-1.7930188) q[2];
rz(0.013269987) q[3];
sx q[3];
rz(-2.477406) q[3];
sx q[3];
rz(-1.9224904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.211798) q[0];
sx q[0];
rz(-3.0199265) q[0];
sx q[0];
rz(-1.9450872) q[0];
rz(-0.78367805) q[1];
sx q[1];
rz(-1.0107026) q[1];
sx q[1];
rz(2.3616621) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62155277) q[0];
sx q[0];
rz(-2.2889458) q[0];
sx q[0];
rz(1.2214127) q[0];
x q[1];
rz(-0.11855306) q[2];
sx q[2];
rz(-2.701512) q[2];
sx q[2];
rz(-1.5755744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1577252) q[1];
sx q[1];
rz(-1.9541426) q[1];
sx q[1];
rz(2.4216837) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50239222) q[3];
sx q[3];
rz(-2.5414782) q[3];
sx q[3];
rz(-0.49139088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1448867) q[2];
sx q[2];
rz(-0.11652623) q[2];
sx q[2];
rz(0.04280002) q[2];
rz(1.5329817) q[3];
sx q[3];
rz(-1.3442842) q[3];
sx q[3];
rz(2.2406421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3881989) q[0];
sx q[0];
rz(-2.2008984) q[0];
sx q[0];
rz(0.70145506) q[0];
rz(2.258621) q[1];
sx q[1];
rz(-1.3621623) q[1];
sx q[1];
rz(-2.0515474) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42692894) q[0];
sx q[0];
rz(-1.8986456) q[0];
sx q[0];
rz(-1.9886243) q[0];
rz(-pi) q[1];
rz(2.6323363) q[2];
sx q[2];
rz(-0.79051566) q[2];
sx q[2];
rz(2.5495286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5970814) q[1];
sx q[1];
rz(-1.7131355) q[1];
sx q[1];
rz(-0.049171731) q[1];
rz(-pi) q[2];
rz(0.16441508) q[3];
sx q[3];
rz(-1.3344171) q[3];
sx q[3];
rz(-2.6316093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8604454) q[2];
sx q[2];
rz(-0.88938418) q[2];
sx q[2];
rz(-1.3775728) q[2];
rz(-3.0456165) q[3];
sx q[3];
rz(-2.4317661) q[3];
sx q[3];
rz(2.6095552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9607361) q[0];
sx q[0];
rz(-2.2611698) q[0];
sx q[0];
rz(-1.9218504) q[0];
rz(-0.60918728) q[1];
sx q[1];
rz(-1.6222619) q[1];
sx q[1];
rz(2.6763197) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0739586) q[0];
sx q[0];
rz(-2.2797757) q[0];
sx q[0];
rz(-0.65517212) q[0];
x q[1];
rz(-2.005079) q[2];
sx q[2];
rz(-2.4854276) q[2];
sx q[2];
rz(2.9796571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45256685) q[1];
sx q[1];
rz(-1.6151607) q[1];
sx q[1];
rz(-2.0173128) q[1];
rz(-1.3749397) q[3];
sx q[3];
rz(-0.80482414) q[3];
sx q[3];
rz(0.63279654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51218) q[2];
sx q[2];
rz(-2.8359154) q[2];
sx q[2];
rz(2.7194887) q[2];
rz(3.1230538) q[3];
sx q[3];
rz(-1.5471349) q[3];
sx q[3];
rz(-0.6137994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2043532) q[0];
sx q[0];
rz(-0.48870191) q[0];
sx q[0];
rz(-2.3204284) q[0];
rz(-2.4429854) q[1];
sx q[1];
rz(-2.3804074) q[1];
sx q[1];
rz(-3.1331114) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7754448) q[0];
sx q[0];
rz(-0.73655948) q[0];
sx q[0];
rz(-2.6790274) q[0];
x q[1];
rz(2.7445265) q[2];
sx q[2];
rz(-2.0609359) q[2];
sx q[2];
rz(-1.9352143) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0331133) q[1];
sx q[1];
rz(-2.3244807) q[1];
sx q[1];
rz(1.2014649) q[1];
x q[2];
rz(-2.2147785) q[3];
sx q[3];
rz(-1.5379526) q[3];
sx q[3];
rz(-2.42086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2681793) q[2];
sx q[2];
rz(-0.11205967) q[2];
sx q[2];
rz(2.6123135) q[2];
rz(-1.4389634) q[3];
sx q[3];
rz(-1.8301423) q[3];
sx q[3];
rz(-2.8992991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5684123) q[0];
sx q[0];
rz(-0.72887623) q[0];
sx q[0];
rz(0.53701425) q[0];
rz(0.88045398) q[1];
sx q[1];
rz(-1.302779) q[1];
sx q[1];
rz(0.74768487) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76067257) q[0];
sx q[0];
rz(-2.2539646) q[0];
sx q[0];
rz(2.0051413) q[0];
rz(-pi) q[1];
rz(2.0751817) q[2];
sx q[2];
rz(-0.73200544) q[2];
sx q[2];
rz(2.5901834) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0084298531) q[1];
sx q[1];
rz(-1.8122261) q[1];
sx q[1];
rz(-1.2397093) q[1];
x q[2];
rz(-1.8219833) q[3];
sx q[3];
rz(-2.1681941) q[3];
sx q[3];
rz(-2.8772014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65413862) q[2];
sx q[2];
rz(-0.70709252) q[2];
sx q[2];
rz(-2.1716165) q[2];
rz(-1.6449432) q[3];
sx q[3];
rz(-1.0140489) q[3];
sx q[3];
rz(1.0223201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3269761) q[0];
sx q[0];
rz(-1.6642445) q[0];
sx q[0];
rz(-2.5388663) q[0];
rz(3.1372435) q[1];
sx q[1];
rz(-1.8252204) q[1];
sx q[1];
rz(-1.1109005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0051188) q[0];
sx q[0];
rz(-2.2181151) q[0];
sx q[0];
rz(2.3939483) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52934016) q[2];
sx q[2];
rz(-0.94369859) q[2];
sx q[2];
rz(-0.68466078) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7676684) q[1];
sx q[1];
rz(-1.6925188) q[1];
sx q[1];
rz(-0.59438225) q[1];
x q[2];
rz(-2.5488126) q[3];
sx q[3];
rz(-2.5912632) q[3];
sx q[3];
rz(-1.569848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6415619) q[2];
sx q[2];
rz(-1.2189453) q[2];
sx q[2];
rz(0.59263372) q[2];
rz(-0.57989132) q[3];
sx q[3];
rz(-2.4903844) q[3];
sx q[3];
rz(1.7013223) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92689571) q[0];
sx q[0];
rz(-2.0669282) q[0];
sx q[0];
rz(2.1672473) q[0];
rz(0.018085619) q[1];
sx q[1];
rz(-2.7985202) q[1];
sx q[1];
rz(0.090029686) q[1];
rz(2.8333673) q[2];
sx q[2];
rz(-1.8147974) q[2];
sx q[2];
rz(-2.360156) q[2];
rz(-2.0578961) q[3];
sx q[3];
rz(-2.2550411) q[3];
sx q[3];
rz(-2.5451345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
