OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8074789) q[0];
sx q[0];
rz(2.4456094) q[0];
sx q[0];
rz(13.453475) q[0];
rz(-0.71169418) q[1];
sx q[1];
rz(5.2206098) q[1];
sx q[1];
rz(9.0945464) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45107156) q[0];
sx q[0];
rz(-2.1981648) q[0];
sx q[0];
rz(2.5894707) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9357968) q[2];
sx q[2];
rz(-2.096039) q[2];
sx q[2];
rz(0.15196063) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.020259) q[1];
sx q[1];
rz(-1.4687045) q[1];
sx q[1];
rz(2.6481204) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0631271) q[3];
sx q[3];
rz(-2.0002504) q[3];
sx q[3];
rz(-1.0965018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7547138) q[2];
sx q[2];
rz(-2.7140996) q[2];
sx q[2];
rz(-1.0270366) q[2];
rz(-0.88296452) q[3];
sx q[3];
rz(-1.3436147) q[3];
sx q[3];
rz(0.55356717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7850194) q[0];
sx q[0];
rz(-3.0632601) q[0];
sx q[0];
rz(-0.26764348) q[0];
rz(1.74125) q[1];
sx q[1];
rz(-0.61336556) q[1];
sx q[1];
rz(0.5563446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9984335) q[0];
sx q[0];
rz(-2.0578009) q[0];
sx q[0];
rz(0.28227455) q[0];
rz(-pi) q[1];
rz(-1.2265967) q[2];
sx q[2];
rz(-2.2453893) q[2];
sx q[2];
rz(0.062002484) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2241321) q[1];
sx q[1];
rz(-1.9497834) q[1];
sx q[1];
rz(2.1716592) q[1];
rz(-0.90640599) q[3];
sx q[3];
rz(-1.0944546) q[3];
sx q[3];
rz(-1.9204467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9050425) q[2];
sx q[2];
rz(-1.6908129) q[2];
sx q[2];
rz(-1.8897918) q[2];
rz(1.8768138) q[3];
sx q[3];
rz(-2.4301961) q[3];
sx q[3];
rz(1.0714162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.65997893) q[0];
sx q[0];
rz(-0.55584207) q[0];
sx q[0];
rz(-2.3023093) q[0];
rz(-1.6370157) q[1];
sx q[1];
rz(-1.4312276) q[1];
sx q[1];
rz(1.7617216) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2502852) q[0];
sx q[0];
rz(-1.4472741) q[0];
sx q[0];
rz(3.059292) q[0];
rz(-pi) q[1];
rz(-1.9154869) q[2];
sx q[2];
rz(-0.57248964) q[2];
sx q[2];
rz(2.9782611) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4288001) q[1];
sx q[1];
rz(-1.8903435) q[1];
sx q[1];
rz(1.0733114) q[1];
rz(-pi) q[2];
rz(2.3776235) q[3];
sx q[3];
rz(-2.4141365) q[3];
sx q[3];
rz(-0.025345552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69924361) q[2];
sx q[2];
rz(-1.9244497) q[2];
sx q[2];
rz(-2.4231363) q[2];
rz(-2.5213304) q[3];
sx q[3];
rz(-2.2060427) q[3];
sx q[3];
rz(0.81234318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45850596) q[0];
sx q[0];
rz(-1.2309256) q[0];
sx q[0];
rz(-1.4076391) q[0];
rz(2.5900224) q[1];
sx q[1];
rz(-3.0307814) q[1];
sx q[1];
rz(1.3161906) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5410446) q[0];
sx q[0];
rz(-1.8529612) q[0];
sx q[0];
rz(2.1938384) q[0];
rz(1.5274996) q[2];
sx q[2];
rz(-0.27972066) q[2];
sx q[2];
rz(-0.86076517) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11034695) q[1];
sx q[1];
rz(-1.4390042) q[1];
sx q[1];
rz(-2.736997) q[1];
rz(-0.67263453) q[3];
sx q[3];
rz(-2.1676873) q[3];
sx q[3];
rz(1.2850645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5482285) q[2];
sx q[2];
rz(-2.7061988) q[2];
sx q[2];
rz(-2.1335874) q[2];
rz(-2.8625782) q[3];
sx q[3];
rz(-2.3271826) q[3];
sx q[3];
rz(0.45849714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27610436) q[0];
sx q[0];
rz(-1.8115598) q[0];
sx q[0];
rz(2.8651067) q[0];
rz(-2.8522988) q[1];
sx q[1];
rz(-1.9464867) q[1];
sx q[1];
rz(2.8791265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55928265) q[0];
sx q[0];
rz(-1.773343) q[0];
sx q[0];
rz(-1.3900533) q[0];
rz(-pi) q[1];
rz(0.97945531) q[2];
sx q[2];
rz(-1.5797414) q[2];
sx q[2];
rz(1.971955) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.561132) q[1];
sx q[1];
rz(-1.9101686) q[1];
sx q[1];
rz(-2.5575693) q[1];
rz(-0.57837242) q[3];
sx q[3];
rz(-1.6993853) q[3];
sx q[3];
rz(-1.8687539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3929954) q[2];
sx q[2];
rz(-0.59610072) q[2];
sx q[2];
rz(1.7945012) q[2];
rz(-0.10284452) q[3];
sx q[3];
rz(-1.2343854) q[3];
sx q[3];
rz(-1.6245406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(3.1384001) q[0];
sx q[0];
rz(-1.6380558) q[0];
sx q[0];
rz(0.060977161) q[0];
rz(-1.2334476) q[1];
sx q[1];
rz(-1.7465218) q[1];
sx q[1];
rz(1.5256418) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4075482) q[0];
sx q[0];
rz(-1.5015232) q[0];
sx q[0];
rz(-2.6114527) q[0];
rz(-pi) q[1];
rz(3.1093237) q[2];
sx q[2];
rz(-2.319558) q[2];
sx q[2];
rz(-2.8634594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4021637) q[1];
sx q[1];
rz(-1.8212178) q[1];
sx q[1];
rz(3.0274505) q[1];
x q[2];
rz(-1.7876225) q[3];
sx q[3];
rz(-2.4783756) q[3];
sx q[3];
rz(1.7182223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6700651) q[2];
sx q[2];
rz(-1.7933041) q[2];
sx q[2];
rz(1.8111551) q[2];
rz(2.0165675) q[3];
sx q[3];
rz(-0.87441134) q[3];
sx q[3];
rz(-0.037954656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1838945) q[0];
sx q[0];
rz(-1.6989468) q[0];
sx q[0];
rz(2.884602) q[0];
rz(-2.7067302) q[1];
sx q[1];
rz(-2.3915274) q[1];
sx q[1];
rz(0.90352568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436179) q[0];
sx q[0];
rz(-1.5666641) q[0];
sx q[0];
rz(-1.6630737) q[0];
x q[1];
rz(-0.14638325) q[2];
sx q[2];
rz(-2.3902262) q[2];
sx q[2];
rz(-2.9379972) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7419646) q[1];
sx q[1];
rz(-0.077721462) q[1];
sx q[1];
rz(-0.6570973) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89374505) q[3];
sx q[3];
rz(-0.66679685) q[3];
sx q[3];
rz(-2.9725985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8187108) q[2];
sx q[2];
rz(-1.633753) q[2];
sx q[2];
rz(-1.7203077) q[2];
rz(0.70982248) q[3];
sx q[3];
rz(-1.1387419) q[3];
sx q[3];
rz(2.6074954) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98599559) q[0];
sx q[0];
rz(-2.7105712) q[0];
sx q[0];
rz(1.0850329) q[0];
rz(0.64708465) q[1];
sx q[1];
rz(-1.9661463) q[1];
sx q[1];
rz(0.064402493) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0376799) q[0];
sx q[0];
rz(-1.896201) q[0];
sx q[0];
rz(-0.66489451) q[0];
rz(-0.66775457) q[2];
sx q[2];
rz(-2.1760094) q[2];
sx q[2];
rz(-2.3507694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.850122) q[1];
sx q[1];
rz(-0.89813342) q[1];
sx q[1];
rz(2.5419278) q[1];
x q[2];
rz(0.91592083) q[3];
sx q[3];
rz(-0.43817156) q[3];
sx q[3];
rz(-0.34153356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.25962466) q[2];
sx q[2];
rz(-0.3825863) q[2];
sx q[2];
rz(0.043370334) q[2];
rz(2.4671593) q[3];
sx q[3];
rz(-1.3573656) q[3];
sx q[3];
rz(1.1446713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6944273) q[0];
sx q[0];
rz(-2.5612216) q[0];
sx q[0];
rz(-0.32661435) q[0];
rz(1.0940374) q[1];
sx q[1];
rz(-2.2973165) q[1];
sx q[1];
rz(2.6577267) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992386) q[0];
sx q[0];
rz(-1.0629553) q[0];
sx q[0];
rz(-0.81224982) q[0];
x q[1];
rz(2.6695374) q[2];
sx q[2];
rz(-0.46923551) q[2];
sx q[2];
rz(-3.0735441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6623508) q[1];
sx q[1];
rz(-1.3562981) q[1];
sx q[1];
rz(2.3420326) q[1];
rz(2.4993091) q[3];
sx q[3];
rz(-1.7574199) q[3];
sx q[3];
rz(0.24333653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6924374) q[2];
sx q[2];
rz(-2.1026976) q[2];
sx q[2];
rz(-3.0628487) q[2];
rz(-0.76830831) q[3];
sx q[3];
rz(-1.2945622) q[3];
sx q[3];
rz(-2.9505742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.47062) q[0];
sx q[0];
rz(-2.5363531) q[0];
sx q[0];
rz(0.11669267) q[0];
rz(0.8600046) q[1];
sx q[1];
rz(-1.8000894) q[1];
sx q[1];
rz(-2.2681627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0104994) q[0];
sx q[0];
rz(-2.2981055) q[0];
sx q[0];
rz(-2.8167679) q[0];
rz(-0.53773138) q[2];
sx q[2];
rz(-2.5641547) q[2];
sx q[2];
rz(0.87960192) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22120093) q[1];
sx q[1];
rz(-1.6947692) q[1];
sx q[1];
rz(2.6392379) q[1];
rz(-1.6386581) q[3];
sx q[3];
rz(-0.96047365) q[3];
sx q[3];
rz(-0.44059702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42150911) q[2];
sx q[2];
rz(-0.24622157) q[2];
sx q[2];
rz(2.1923547) q[2];
rz(-1.2468437) q[3];
sx q[3];
rz(-1.6722164) q[3];
sx q[3];
rz(0.59103549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0134037) q[0];
sx q[0];
rz(-1.0008151) q[0];
sx q[0];
rz(-1.6304954) q[0];
rz(-2.7746101) q[1];
sx q[1];
rz(-1.3584247) q[1];
sx q[1];
rz(-0.57986837) q[1];
rz(2.8383474) q[2];
sx q[2];
rz(-1.1228391) q[2];
sx q[2];
rz(-0.11930677) q[2];
rz(-0.49574481) q[3];
sx q[3];
rz(-0.75400018) q[3];
sx q[3];
rz(-0.41771623) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
