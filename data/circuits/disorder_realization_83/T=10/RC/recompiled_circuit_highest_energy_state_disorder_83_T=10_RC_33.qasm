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
rz(2.7163765) q[0];
sx q[0];
rz(-1.334231) q[0];
sx q[0];
rz(0.38036007) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(2.9549197) q[1];
sx q[1];
rz(8.5001707) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43792576) q[0];
sx q[0];
rz(-0.80090085) q[0];
sx q[0];
rz(-2.6573703) q[0];
rz(-pi) q[1];
rz(2.9745462) q[2];
sx q[2];
rz(-1.0606597) q[2];
sx q[2];
rz(1.270592) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.315359) q[1];
sx q[1];
rz(-0.72972882) q[1];
sx q[1];
rz(1.5118096) q[1];
x q[2];
rz(-3.0713586) q[3];
sx q[3];
rz(-3.0133504) q[3];
sx q[3];
rz(1.2508891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1251462) q[2];
sx q[2];
rz(-0.95637286) q[2];
sx q[2];
rz(0.98801405) q[2];
rz(-2.9347349) q[3];
sx q[3];
rz(-1.7346953) q[3];
sx q[3];
rz(2.701581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4619231) q[0];
sx q[0];
rz(-1.4689057) q[0];
sx q[0];
rz(-2.8894506) q[0];
rz(3.120046) q[1];
sx q[1];
rz(-0.69308678) q[1];
sx q[1];
rz(-2.2861939) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6538174) q[0];
sx q[0];
rz(-0.027055351) q[0];
sx q[0];
rz(1.6271126) q[0];
x q[1];
rz(-2.5874675) q[2];
sx q[2];
rz(-1.9337855) q[2];
sx q[2];
rz(2.3100694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8586848) q[1];
sx q[1];
rz(-2.5577684) q[1];
sx q[1];
rz(1.1633384) q[1];
rz(-pi) q[2];
rz(-1.9048433) q[3];
sx q[3];
rz(-2.2260087) q[3];
sx q[3];
rz(-1.1545336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1713193) q[2];
sx q[2];
rz(-0.0042985175) q[2];
sx q[2];
rz(0.26101905) q[2];
rz(-0.039693443) q[3];
sx q[3];
rz(-1.3986162) q[3];
sx q[3];
rz(0.54350129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-2.4259341) q[0];
sx q[0];
rz(-1.6682699) q[0];
sx q[0];
rz(0.69450992) q[0];
rz(2.014324) q[1];
sx q[1];
rz(-0.85113168) q[1];
sx q[1];
rz(2.1515501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2106101) q[0];
sx q[0];
rz(-1.4890122) q[0];
sx q[0];
rz(2.7922996) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9431413) q[2];
sx q[2];
rz(-0.99764148) q[2];
sx q[2];
rz(2.1015374) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.5824008) q[1];
sx q[1];
rz(-1.2123931) q[1];
sx q[1];
rz(2.6646975) q[1];
rz(-2.18543) q[3];
sx q[3];
rz(-2.3376138) q[3];
sx q[3];
rz(-1.1756736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29828829) q[2];
sx q[2];
rz(-2.1683606) q[2];
sx q[2];
rz(0.49015552) q[2];
rz(2.7847024) q[3];
sx q[3];
rz(-2.7481952) q[3];
sx q[3];
rz(1.8753768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.056034293) q[0];
sx q[0];
rz(-1.1857251) q[0];
sx q[0];
rz(0.84247843) q[0];
rz(-0.8017686) q[1];
sx q[1];
rz(-0.27920488) q[1];
sx q[1];
rz(-2.3559949) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2548837) q[0];
sx q[0];
rz(-2.2629316) q[0];
sx q[0];
rz(-1.2413625) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7217721) q[2];
sx q[2];
rz(-2.2822126) q[2];
sx q[2];
rz(0.47455088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3093058) q[1];
sx q[1];
rz(-2.5046299) q[1];
sx q[1];
rz(2.6061406) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3219452) q[3];
sx q[3];
rz(-1.420212) q[3];
sx q[3];
rz(2.8898847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12513146) q[2];
sx q[2];
rz(-2.3931914) q[2];
sx q[2];
rz(1.008519) q[2];
rz(-0.85159167) q[3];
sx q[3];
rz(-0.75751704) q[3];
sx q[3];
rz(-2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9652902) q[0];
sx q[0];
rz(-2.1602614) q[0];
sx q[0];
rz(-1.0705795) q[0];
rz(0.77752441) q[1];
sx q[1];
rz(-2.2309525) q[1];
sx q[1];
rz(-2.1308965) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010653822) q[0];
sx q[0];
rz(-1.9294327) q[0];
sx q[0];
rz(2.6333195) q[0];
rz(-pi) q[1];
rz(1.7456876) q[2];
sx q[2];
rz(-1.5064459) q[2];
sx q[2];
rz(-1.627587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1221558) q[1];
sx q[1];
rz(-2.2194462) q[1];
sx q[1];
rz(-0.98585702) q[1];
rz(2.9845731) q[3];
sx q[3];
rz(-2.5738705) q[3];
sx q[3];
rz(-2.1322676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8396478) q[2];
sx q[2];
rz(-3.0366615) q[2];
sx q[2];
rz(-1.5555752) q[2];
rz(1.0157061) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(1.4122081) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509336) q[0];
sx q[0];
rz(-0.18336329) q[0];
sx q[0];
rz(2.3405128) q[0];
rz(3.0084897) q[1];
sx q[1];
rz(-1.8201273) q[1];
sx q[1];
rz(-2.7395693) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4426081) q[0];
sx q[0];
rz(-0.8359209) q[0];
sx q[0];
rz(-1.1135654) q[0];
rz(-pi) q[1];
rz(-3.1394464) q[2];
sx q[2];
rz(-1.8615926) q[2];
sx q[2];
rz(-2.4418497) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1014598) q[1];
sx q[1];
rz(-1.8454843) q[1];
sx q[1];
rz(2.0707612) q[1];
x q[2];
rz(0.65816718) q[3];
sx q[3];
rz(-1.6699381) q[3];
sx q[3];
rz(1.0189354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.63783995) q[2];
sx q[2];
rz(-1.736234) q[2];
sx q[2];
rz(2.3657738) q[2];
rz(-1.6860298) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(2.1337401) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044947226) q[0];
sx q[0];
rz(-2.2438887) q[0];
sx q[0];
rz(0.67895472) q[0];
rz(-1.8567765) q[1];
sx q[1];
rz(-0.7130475) q[1];
sx q[1];
rz(-1.7624034) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75684568) q[0];
sx q[0];
rz(-2.3492536) q[0];
sx q[0];
rz(-1.3261262) q[0];
rz(-pi) q[1];
rz(-3.0581362) q[2];
sx q[2];
rz(-1.8493358) q[2];
sx q[2];
rz(0.73501529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1759291) q[1];
sx q[1];
rz(-2.0116076) q[1];
sx q[1];
rz(-0.55056527) q[1];
rz(-pi) q[2];
rz(-1.3706743) q[3];
sx q[3];
rz(-1.516023) q[3];
sx q[3];
rz(-1.9211662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8046367) q[2];
sx q[2];
rz(-0.75504428) q[2];
sx q[2];
rz(1.9026559) q[2];
rz(-1.8094481) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(2.8968887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181344) q[0];
sx q[0];
rz(-1.1661538) q[0];
sx q[0];
rz(-0.13846692) q[0];
rz(-0.18374099) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(-0.26652452) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0154769) q[0];
sx q[0];
rz(-1.4086282) q[0];
sx q[0];
rz(1.2408942) q[0];
rz(1.1600003) q[2];
sx q[2];
rz(-0.60705429) q[2];
sx q[2];
rz(1.0714873) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6147242) q[1];
sx q[1];
rz(-0.42335948) q[1];
sx q[1];
rz(-1.4335267) q[1];
x q[2];
rz(3.1317461) q[3];
sx q[3];
rz(-0.72692633) q[3];
sx q[3];
rz(1.7444514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0900241) q[2];
sx q[2];
rz(-1.8847621) q[2];
sx q[2];
rz(-2.9070692) q[2];
rz(1.6656434) q[3];
sx q[3];
rz(-1.6022976) q[3];
sx q[3];
rz(-0.75638151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.804857) q[0];
sx q[0];
rz(-2.2754301) q[0];
sx q[0];
rz(2.3418703) q[0];
rz(2.5526478) q[1];
sx q[1];
rz(-2.2942693) q[1];
sx q[1];
rz(-0.2074997) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0711771) q[0];
sx q[0];
rz(-1.1780329) q[0];
sx q[0];
rz(1.7123187) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8325808) q[2];
sx q[2];
rz(-0.88243077) q[2];
sx q[2];
rz(-0.38028827) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7991528) q[1];
sx q[1];
rz(-0.9121597) q[1];
sx q[1];
rz(2.6299146) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5141684) q[3];
sx q[3];
rz(-2.671173) q[3];
sx q[3];
rz(-0.75477598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3346682) q[2];
sx q[2];
rz(-2.25756) q[2];
sx q[2];
rz(2.7743288) q[2];
rz(1.4372829) q[3];
sx q[3];
rz(-1.1742914) q[3];
sx q[3];
rz(-1.1737163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87840286) q[0];
sx q[0];
rz(-1.0028361) q[0];
sx q[0];
rz(-2.3314085) q[0];
rz(1.1148249) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(-0.55327639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033940559) q[0];
sx q[0];
rz(-0.93354152) q[0];
sx q[0];
rz(-0.069038387) q[0];
x q[1];
rz(2.3478339) q[2];
sx q[2];
rz(-2.0362034) q[2];
sx q[2];
rz(-2.2528354) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3685963) q[1];
sx q[1];
rz(-2.8008411) q[1];
sx q[1];
rz(-2.1590538) q[1];
rz(-pi) q[2];
rz(-2.3555967) q[3];
sx q[3];
rz(-1.3162426) q[3];
sx q[3];
rz(2.8543775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0193923) q[2];
sx q[2];
rz(-0.70851749) q[2];
sx q[2];
rz(-2.9019287) q[2];
rz(1.3278809) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(0.46752587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0620621) q[0];
sx q[0];
rz(-1.3539599) q[0];
sx q[0];
rz(0.39394105) q[0];
rz(0.65555864) q[1];
sx q[1];
rz(-1.1759023) q[1];
sx q[1];
rz(-2.1942153) q[1];
rz(1.2389567) q[2];
sx q[2];
rz(-1.17579) q[2];
sx q[2];
rz(-0.50926846) q[2];
rz(1.3169133) q[3];
sx q[3];
rz(-1.0051654) q[3];
sx q[3];
rz(-1.2011423) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
