OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(-2.2572416) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8774672) q[0];
sx q[0];
rz(-1.2116417) q[0];
sx q[0];
rz(2.9496664) q[0];
rz(-pi) q[1];
rz(1.5760954) q[2];
sx q[2];
rz(-2.1359348) q[2];
sx q[2];
rz(1.1331171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2734387) q[1];
sx q[1];
rz(-1.8167129) q[1];
sx q[1];
rz(1.6383365) q[1];
rz(-pi) q[2];
rz(2.7847071) q[3];
sx q[3];
rz(-0.89324739) q[3];
sx q[3];
rz(-0.68912904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.78645906) q[2];
sx q[2];
rz(-2.3278475) q[2];
sx q[2];
rz(-0.65594977) q[2];
rz(-1.9338699) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(0.99457994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8939963) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(0.43352747) q[0];
rz(-2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(0.0072335009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7318856) q[0];
sx q[0];
rz(-1.1792372) q[0];
sx q[0];
rz(-0.10086985) q[0];
rz(-pi) q[1];
rz(-2.8480808) q[2];
sx q[2];
rz(-1.8732757) q[2];
sx q[2];
rz(2.7345865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4870287) q[1];
sx q[1];
rz(-1.9332758) q[1];
sx q[1];
rz(-1.3351721) q[1];
rz(-pi) q[2];
rz(-0.98874436) q[3];
sx q[3];
rz(-0.75609222) q[3];
sx q[3];
rz(-0.83272979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(-0.69765222) q[2];
rz(-0.12156045) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(-2.8377623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6737297) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(-1.3695705) q[0];
rz(1.9000152) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(-0.26161584) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75152552) q[0];
sx q[0];
rz(-1.5887504) q[0];
sx q[0];
rz(0.0096339027) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0171266) q[2];
sx q[2];
rz(-0.83041588) q[2];
sx q[2];
rz(0.24701961) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2178104) q[1];
sx q[1];
rz(-2.3465956) q[1];
sx q[1];
rz(1.7407106) q[1];
rz(1.6022127) q[3];
sx q[3];
rz(-2.0835702) q[3];
sx q[3];
rz(-2.153742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4613142) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(2.938081) q[2];
rz(0.92173785) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1638284) q[0];
sx q[0];
rz(-1.5638567) q[0];
sx q[0];
rz(1.571636) q[0];
rz(1.0034026) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(-1.8932231) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0715863) q[0];
sx q[0];
rz(-1.6674111) q[0];
sx q[0];
rz(-3.0759401) q[0];
rz(1.7522881) q[2];
sx q[2];
rz(-1.776812) q[2];
sx q[2];
rz(0.69603053) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.26989386) q[1];
sx q[1];
rz(-1.4596241) q[1];
sx q[1];
rz(-2.486869) q[1];
rz(-1.3094041) q[3];
sx q[3];
rz(-2.2769615) q[3];
sx q[3];
rz(0.9725001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1127597) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(-0.83703414) q[2];
rz(1.2083496) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(0.77520448) q[0];
rz(-0.40183055) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(2.2391589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1989312) q[0];
sx q[0];
rz(-0.64093243) q[0];
sx q[0];
rz(3.065227) q[0];
rz(-pi) q[1];
rz(-0.47307737) q[2];
sx q[2];
rz(-0.50191754) q[2];
sx q[2];
rz(0.15854533) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.202995) q[1];
sx q[1];
rz(-1.0293759) q[1];
sx q[1];
rz(-0.90403892) q[1];
x q[2];
rz(3.1387572) q[3];
sx q[3];
rz(-1.2540725) q[3];
sx q[3];
rz(-1.0032652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19501413) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(-1.1966594) q[2];
rz(-1.6992016) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(1.4590013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3301795) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(0.50317558) q[0];
rz(1.6852089) q[1];
sx q[1];
rz(-1.0738942) q[1];
sx q[1];
rz(-2.9398289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1267032) q[0];
sx q[0];
rz(-0.12173437) q[0];
sx q[0];
rz(-0.99862167) q[0];
x q[1];
rz(-0.47735881) q[2];
sx q[2];
rz(-1.9743894) q[2];
sx q[2];
rz(-1.3178283) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0338248) q[1];
sx q[1];
rz(-2.0467842) q[1];
sx q[1];
rz(1.0955217) q[1];
rz(-pi) q[2];
rz(2.4089036) q[3];
sx q[3];
rz(-1.2142039) q[3];
sx q[3];
rz(-2.4458812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9266944) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(-0.26724896) q[2];
rz(0.8231419) q[3];
sx q[3];
rz(-3.0624793) q[3];
sx q[3];
rz(2.2657623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(1.293752) q[0];
sx q[0];
rz(-0.1593312) q[0];
sx q[0];
rz(3.0840432) q[0];
rz(-1.4808222) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(0.94271359) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1377624) q[0];
sx q[0];
rz(-0.98235213) q[0];
sx q[0];
rz(2.0501191) q[0];
rz(0.96037453) q[2];
sx q[2];
rz(-0.27563169) q[2];
sx q[2];
rz(-0.20197091) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47989935) q[1];
sx q[1];
rz(-2.6091895) q[1];
sx q[1];
rz(2.2366621) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0122635) q[3];
sx q[3];
rz(-2.5848476) q[3];
sx q[3];
rz(0.41551513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1192347) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(2.7589202) q[2];
rz(-2.102397) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(-0.97810811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7169749) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(2.8714645) q[0];
rz(-2.5121636) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(0.28392917) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1536381) q[0];
sx q[0];
rz(-1.1414014) q[0];
sx q[0];
rz(0.9149787) q[0];
x q[1];
rz(0.59992744) q[2];
sx q[2];
rz(-1.4495965) q[2];
sx q[2];
rz(-1.8184513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64360196) q[1];
sx q[1];
rz(-1.0712578) q[1];
sx q[1];
rz(-2.1010146) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1030032) q[3];
sx q[3];
rz(-2.4814197) q[3];
sx q[3];
rz(2.431228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5876864) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(1.930687) q[2];
rz(2.8816913) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(-0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07638409) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(2.912345) q[0];
rz(-2.8385838) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(1.680826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0871353) q[0];
sx q[0];
rz(-1.5883049) q[0];
sx q[0];
rz(-2.8580335) q[0];
x q[1];
rz(-1.7503579) q[2];
sx q[2];
rz(-1.944724) q[2];
sx q[2];
rz(1.1073081) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8252392) q[1];
sx q[1];
rz(-1.5484527) q[1];
sx q[1];
rz(-2.6166726) q[1];
x q[2];
rz(-1.9577515) q[3];
sx q[3];
rz(-2.317252) q[3];
sx q[3];
rz(-1.6798306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4298657) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(2.7097278) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(-0.7152043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(-2.8046872) q[0];
rz(2.9341872) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(-0.38063231) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1487409) q[0];
sx q[0];
rz(-1.7721575) q[0];
sx q[0];
rz(-3.0701748) q[0];
x q[1];
rz(0.74110465) q[2];
sx q[2];
rz(-1.1198824) q[2];
sx q[2];
rz(1.2620743) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14643529) q[1];
sx q[1];
rz(-0.53968118) q[1];
sx q[1];
rz(-1.9567009) q[1];
rz(-pi) q[2];
rz(2.5467039) q[3];
sx q[3];
rz(-2.4747362) q[3];
sx q[3];
rz(-0.24645933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(2.005119) q[2];
rz(0.040955695) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8052335) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(2.1144755) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(-1.3409875) q[2];
sx q[2];
rz(-1.8125712) q[2];
sx q[2];
rz(-0.3704091) q[2];
rz(-1.1313664) q[3];
sx q[3];
rz(-1.3888748) q[3];
sx q[3];
rz(0.5007762) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];