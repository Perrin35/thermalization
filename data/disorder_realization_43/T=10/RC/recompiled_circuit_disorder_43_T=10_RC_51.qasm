OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2551978) q[0];
sx q[0];
rz(-1.891341) q[0];
sx q[0];
rz(1.3347081) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(0.97595739) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6061493) q[0];
sx q[0];
rz(-0.9978928) q[0];
sx q[0];
rz(-1.2826305) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1337778) q[2];
sx q[2];
rz(-1.8554167) q[2];
sx q[2];
rz(-0.1422589) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0426892) q[1];
sx q[1];
rz(-2.1747723) q[1];
sx q[1];
rz(0.16278111) q[1];
rz(-0.8043886) q[3];
sx q[3];
rz(-1.8612923) q[3];
sx q[3];
rz(-0.94798541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41539899) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(2.8090254) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(-1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48288229) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(-2.9630307) q[0];
rz(-1.8042971) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(3.1352502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9316677) q[0];
sx q[0];
rz(-1.7010265) q[0];
sx q[0];
rz(-1.3757214) q[0];
rz(0.96887529) q[2];
sx q[2];
rz(-1.1224147) q[2];
sx q[2];
rz(-2.7768163) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9176863) q[1];
sx q[1];
rz(-2.8250541) q[1];
sx q[1];
rz(-2.9427337) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8216324) q[3];
sx q[3];
rz(-2.1187966) q[3];
sx q[3];
rz(2.8743924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1938842) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(0.87810278) q[2];
rz(0.86205035) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(-2.7745461) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(0.095741622) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.028713) q[0];
sx q[0];
rz(-0.23447795) q[0];
sx q[0];
rz(-2.2857091) q[0];
rz(-pi) q[1];
rz(2.6355866) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(-2.4207052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0248191) q[1];
sx q[1];
rz(-2.1110592) q[1];
sx q[1];
rz(2.5953672) q[1];
rz(-pi) q[2];
rz(-1.6060353) q[3];
sx q[3];
rz(-2.1217151) q[3];
sx q[3];
rz(-2.2743724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0694971) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(2.3068008) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(-2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9469706) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(-2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(-1.0338763) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945275) q[0];
sx q[0];
rz(-2.3405502) q[0];
sx q[0];
rz(2.4497776) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6150377) q[2];
sx q[2];
rz(-2.0806081) q[2];
sx q[2];
rz(1.9499792) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.66400601) q[1];
sx q[1];
rz(-0.2020745) q[1];
sx q[1];
rz(-1.9007773) q[1];
rz(-1.8978118) q[3];
sx q[3];
rz(-1.838284) q[3];
sx q[3];
rz(-0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(-2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088257) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-0.41473266) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(-0.57410747) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2643471) q[0];
sx q[0];
rz(-0.97752042) q[0];
sx q[0];
rz(-2.9869153) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5482043) q[2];
sx q[2];
rz(-0.85561692) q[2];
sx q[2];
rz(-1.1720282) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1787604) q[1];
sx q[1];
rz(-1.292359) q[1];
sx q[1];
rz(0.36414418) q[1];
rz(1.3247213) q[3];
sx q[3];
rz(-0.68021357) q[3];
sx q[3];
rz(-2.084793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.85270143) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(-1.5138907) q[2];
rz(-0.04143516) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(-3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067327499) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(-2.4434027) q[0];
rz(0.16695887) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(-0.73227698) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66757827) q[0];
sx q[0];
rz(-1.3587553) q[0];
sx q[0];
rz(1.7330806) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4562624) q[2];
sx q[2];
rz(-2.2517859) q[2];
sx q[2];
rz(-2.6076917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9975035) q[1];
sx q[1];
rz(-0.92458506) q[1];
sx q[1];
rz(1.8540107) q[1];
x q[2];
rz(-0.5703339) q[3];
sx q[3];
rz(-1.7551646) q[3];
sx q[3];
rz(2.1031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7531062) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(1.9801271) q[2];
rz(-0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(-1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5724065) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(0.15643315) q[0];
rz(-2.6898443) q[1];
sx q[1];
rz(-0.86507559) q[1];
sx q[1];
rz(0.77004534) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5891588) q[0];
sx q[0];
rz(-1.3981515) q[0];
sx q[0];
rz(2.5347559) q[0];
rz(-0.98203512) q[2];
sx q[2];
rz(-1.2274449) q[2];
sx q[2];
rz(-2.8729168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7913831) q[1];
sx q[1];
rz(-2.0767127) q[1];
sx q[1];
rz(0.70178589) q[1];
x q[2];
rz(1.2392063) q[3];
sx q[3];
rz(-1.0911897) q[3];
sx q[3];
rz(-2.1945206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4975171) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(-1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62676936) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(2.8998937) q[0];
rz(-2.4027951) q[1];
sx q[1];
rz(-2.6530478) q[1];
sx q[1];
rz(0.36639211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9241656) q[0];
sx q[0];
rz(-2.1806742) q[0];
sx q[0];
rz(0.92745552) q[0];
rz(0.56717746) q[2];
sx q[2];
rz(-2.4293373) q[2];
sx q[2];
rz(2.8380307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0802676) q[1];
sx q[1];
rz(-0.72239164) q[1];
sx q[1];
rz(-2.8894043) q[1];
rz(-pi) q[2];
rz(1.6493158) q[3];
sx q[3];
rz(-0.55320569) q[3];
sx q[3];
rz(-1.5403403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-0.29385847) q[2];
rz(-0.014523225) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095187) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(2.2576387) q[0];
rz(0.66028315) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(0.79137897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34245472) q[0];
sx q[0];
rz(-2.236811) q[0];
sx q[0];
rz(-0.40823437) q[0];
x q[1];
rz(2.6762166) q[2];
sx q[2];
rz(-1.1306136) q[2];
sx q[2];
rz(-2.344775) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1014175) q[1];
sx q[1];
rz(-1.1244332) q[1];
sx q[1];
rz(-2.0324213) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16718849) q[3];
sx q[3];
rz(-1.0382004) q[3];
sx q[3];
rz(0.54824588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.066594921) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(-2.1203314) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1148949) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-2.5337906) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.7609319) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7847583) q[0];
sx q[0];
rz(-0.94576242) q[0];
sx q[0];
rz(-0.16108315) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5153377) q[2];
sx q[2];
rz(-0.43606731) q[2];
sx q[2];
rz(-0.01576327) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3874515) q[1];
sx q[1];
rz(-1.8879461) q[1];
sx q[1];
rz(1.1412568) q[1];
x q[2];
rz(-1.4019743) q[3];
sx q[3];
rz(-1.3843378) q[3];
sx q[3];
rz(0.25690119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3537139) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(2.805368) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(-1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375785) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(-0.54429383) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(-1.2056247) q[2];
sx q[2];
rz(-2.361459) q[2];
sx q[2];
rz(-0.38103719) q[2];
rz(-1.929677) q[3];
sx q[3];
rz(-0.59960312) q[3];
sx q[3];
rz(-3.0740769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
