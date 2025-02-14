OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.06935057) q[0];
sx q[0];
rz(7.3370186) q[0];
sx q[0];
rz(10.673527) q[0];
rz(1.9034003) q[1];
sx q[1];
rz(-0.44013953) q[1];
sx q[1];
rz(-1.8990489) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7799067) q[0];
sx q[0];
rz(-0.70915993) q[0];
sx q[0];
rz(-2.5283924) q[0];
x q[1];
rz(-3.0468175) q[2];
sx q[2];
rz(-0.80033007) q[2];
sx q[2];
rz(-2.2781585) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5091376) q[1];
sx q[1];
rz(-1.5702562) q[1];
sx q[1];
rz(2.8191791) q[1];
x q[2];
rz(2.6655198) q[3];
sx q[3];
rz(-0.66573788) q[3];
sx q[3];
rz(1.2749323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.093546346) q[2];
sx q[2];
rz(-2.8600433) q[2];
sx q[2];
rz(-1.9264889) q[2];
rz(1.9563227) q[3];
sx q[3];
rz(-0.90648854) q[3];
sx q[3];
rz(-1.5989446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854061) q[0];
sx q[0];
rz(-0.39919272) q[0];
sx q[0];
rz(-1.4584374) q[0];
rz(-2.9996808) q[1];
sx q[1];
rz(-1.2071573) q[1];
sx q[1];
rz(1.9292319) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.713341) q[0];
sx q[0];
rz(-1.9104092) q[0];
sx q[0];
rz(-0.88554145) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0040063695) q[2];
sx q[2];
rz(-2.5494583) q[2];
sx q[2];
rz(-0.59321813) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3691359) q[1];
sx q[1];
rz(-1.6650272) q[1];
sx q[1];
rz(0.89558954) q[1];
rz(1.2207915) q[3];
sx q[3];
rz(-2.1560966) q[3];
sx q[3];
rz(-3.0452951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.045018) q[2];
sx q[2];
rz(-0.18288945) q[2];
sx q[2];
rz(-0.96192399) q[2];
rz(-0.19790459) q[3];
sx q[3];
rz(-1.3062545) q[3];
sx q[3];
rz(-1.643868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69979954) q[0];
sx q[0];
rz(-0.72023359) q[0];
sx q[0];
rz(-1.066347) q[0];
rz(-2.9329246) q[1];
sx q[1];
rz(-0.51869789) q[1];
sx q[1];
rz(-2.4934703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3828165) q[0];
sx q[0];
rz(-2.1471239) q[0];
sx q[0];
rz(0.062950171) q[0];
rz(-2.6510973) q[2];
sx q[2];
rz(-1.6332262) q[2];
sx q[2];
rz(1.708781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7327506) q[1];
sx q[1];
rz(-0.84580814) q[1];
sx q[1];
rz(-2.6413061) q[1];
rz(-pi) q[2];
rz(-3.0410679) q[3];
sx q[3];
rz(-2.4320452) q[3];
sx q[3];
rz(1.2840934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2842747) q[2];
sx q[2];
rz(-1.4633353) q[2];
sx q[2];
rz(-0.54764444) q[2];
rz(2.137843) q[3];
sx q[3];
rz(-2.818483) q[3];
sx q[3];
rz(2.9799262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52466398) q[0];
sx q[0];
rz(-2.014761) q[0];
sx q[0];
rz(-2.774985) q[0];
rz(-1.8560575) q[1];
sx q[1];
rz(-1.6211685) q[1];
sx q[1];
rz(-2.3191648) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7290958) q[0];
sx q[0];
rz(-1.3733882) q[0];
sx q[0];
rz(1.6486932) q[0];
x q[1];
rz(-1.3011312) q[2];
sx q[2];
rz(-1.7326151) q[2];
sx q[2];
rz(1.9105034) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.3318204) q[1];
sx q[1];
rz(-2.0533516) q[1];
sx q[1];
rz(0.81565522) q[1];
x q[2];
rz(-2.8662258) q[3];
sx q[3];
rz(-1.5689225) q[3];
sx q[3];
rz(-3.0089859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.503868) q[2];
sx q[2];
rz(-0.34054264) q[2];
sx q[2];
rz(1.5558422) q[2];
rz(1.9147035) q[3];
sx q[3];
rz(-2.5033958) q[3];
sx q[3];
rz(1.3251806) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1894839) q[0];
sx q[0];
rz(-1.1186849) q[0];
sx q[0];
rz(-0.9504016) q[0];
rz(-0.19418007) q[1];
sx q[1];
rz(-1.8870528) q[1];
sx q[1];
rz(0.28087428) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062894363) q[0];
sx q[0];
rz(-0.67606976) q[0];
sx q[0];
rz(-1.1514787) q[0];
rz(2.9921586) q[2];
sx q[2];
rz(-1.5583056) q[2];
sx q[2];
rz(-0.51692671) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7627613) q[1];
sx q[1];
rz(-1.647659) q[1];
sx q[1];
rz(-0.084048653) q[1];
rz(-pi) q[2];
rz(-2.836197) q[3];
sx q[3];
rz(-1.0825233) q[3];
sx q[3];
rz(2.0528169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6516271) q[2];
sx q[2];
rz(-1.4474892) q[2];
sx q[2];
rz(0.70072407) q[2];
rz(-1.0939595) q[3];
sx q[3];
rz(-2.14812) q[3];
sx q[3];
rz(2.3195364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99264282) q[0];
sx q[0];
rz(-1.0473017) q[0];
sx q[0];
rz(-2.6348422) q[0];
rz(1.1243593) q[1];
sx q[1];
rz(-1.9729112) q[1];
sx q[1];
rz(-0.36875025) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0848727) q[0];
sx q[0];
rz(-2.4115926) q[0];
sx q[0];
rz(-1.2873136) q[0];
rz(-pi) q[1];
rz(-1.1738335) q[2];
sx q[2];
rz(-0.78437524) q[2];
sx q[2];
rz(2.337817) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1146093) q[1];
sx q[1];
rz(-2.8294246) q[1];
sx q[1];
rz(-1.8725558) q[1];
rz(0.82225497) q[3];
sx q[3];
rz(-2.5516627) q[3];
sx q[3];
rz(-1.1168036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2345978) q[2];
sx q[2];
rz(-0.58595053) q[2];
sx q[2];
rz(2.9242945) q[2];
rz(-1.6654738) q[3];
sx q[3];
rz(-0.81513351) q[3];
sx q[3];
rz(-0.34172094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9076964) q[0];
sx q[0];
rz(-0.32231575) q[0];
sx q[0];
rz(2.3436558) q[0];
rz(-1.7012874) q[1];
sx q[1];
rz(-1.0402352) q[1];
sx q[1];
rz(-0.61680102) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9413404) q[0];
sx q[0];
rz(-1.495122) q[0];
sx q[0];
rz(0.6348293) q[0];
rz(-pi) q[1];
rz(-2.6404713) q[2];
sx q[2];
rz(-0.62868147) q[2];
sx q[2];
rz(-1.4528265) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0500461) q[1];
sx q[1];
rz(-1.3669479) q[1];
sx q[1];
rz(-0.14728228) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7019667) q[3];
sx q[3];
rz(-2.1118374) q[3];
sx q[3];
rz(-2.5000985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0133609) q[2];
sx q[2];
rz(-1.1527454) q[2];
sx q[2];
rz(2.4398003) q[2];
rz(-2.2026786) q[3];
sx q[3];
rz(-0.89812583) q[3];
sx q[3];
rz(1.2452589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.668648) q[0];
sx q[0];
rz(-2.9845147) q[0];
sx q[0];
rz(-3.0182086) q[0];
rz(2.3911047) q[1];
sx q[1];
rz(-1.4742943) q[1];
sx q[1];
rz(1.0184681) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7582507) q[0];
sx q[0];
rz(-1.6755548) q[0];
sx q[0];
rz(-1.7927756) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46393259) q[2];
sx q[2];
rz(-1.1008769) q[2];
sx q[2];
rz(-1.4063094) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2606517) q[1];
sx q[1];
rz(-2.6545432) q[1];
sx q[1];
rz(-3.0609291) q[1];
rz(-0.94208053) q[3];
sx q[3];
rz(-1.6438369) q[3];
sx q[3];
rz(-1.4062509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5156775) q[2];
sx q[2];
rz(-1.6720142) q[2];
sx q[2];
rz(-2.6061626) q[2];
rz(1.911602) q[3];
sx q[3];
rz(-2.1401236) q[3];
sx q[3];
rz(0.011822239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5984421) q[0];
sx q[0];
rz(-2.4813528) q[0];
sx q[0];
rz(-1.862233) q[0];
rz(-1.8774425) q[1];
sx q[1];
rz(-1.8708355) q[1];
sx q[1];
rz(0.15170161) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063552) q[0];
sx q[0];
rz(-1.6280109) q[0];
sx q[0];
rz(-0.39852674) q[0];
rz(-2.9931917) q[2];
sx q[2];
rz(-3.0062468) q[2];
sx q[2];
rz(-2.2035408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38167414) q[1];
sx q[1];
rz(-1.2592788) q[1];
sx q[1];
rz(-0.98656922) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0284729) q[3];
sx q[3];
rz(-0.92307011) q[3];
sx q[3];
rz(0.59779378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.58574) q[2];
sx q[2];
rz(-0.93583217) q[2];
sx q[2];
rz(-1.9110511) q[2];
rz(-1.7963643) q[3];
sx q[3];
rz(-1.3627005) q[3];
sx q[3];
rz(1.2983373) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96974385) q[0];
sx q[0];
rz(-1.8010704) q[0];
sx q[0];
rz(-2.6961683) q[0];
rz(-2.716966) q[1];
sx q[1];
rz(-1.5698965) q[1];
sx q[1];
rz(2.5158023) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15881702) q[0];
sx q[0];
rz(-2.3390769) q[0];
sx q[0];
rz(0.83440749) q[0];
rz(1.2568057) q[2];
sx q[2];
rz(-2.6382448) q[2];
sx q[2];
rz(-1.7679917) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0808472) q[1];
sx q[1];
rz(-1.0506127) q[1];
sx q[1];
rz(0.29436269) q[1];
x q[2];
rz(2.8938953) q[3];
sx q[3];
rz(-0.15532914) q[3];
sx q[3];
rz(-2.4904136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9618535) q[2];
sx q[2];
rz(-2.1135795) q[2];
sx q[2];
rz(-1.6801768) q[2];
rz(0.05750582) q[3];
sx q[3];
rz(-1.688136) q[3];
sx q[3];
rz(0.97829515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4089324) q[0];
sx q[0];
rz(-1.4807777) q[0];
sx q[0];
rz(-2.5430191) q[0];
rz(-0.21845017) q[1];
sx q[1];
rz(-1.5427867) q[1];
sx q[1];
rz(-1.3839518) q[1];
rz(-1.4870395) q[2];
sx q[2];
rz(-2.1154006) q[2];
sx q[2];
rz(-1.5474609) q[2];
rz(-0.784008) q[3];
sx q[3];
rz(-1.6029458) q[3];
sx q[3];
rz(-2.8523469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
