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
rz(-2.8357808) q[0];
sx q[0];
rz(-0.6337136) q[0];
sx q[0];
rz(-0.60854882) q[0];
rz(2.4701056) q[1];
sx q[1];
rz(-0.66317135) q[1];
sx q[1];
rz(1.6821678) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7072351) q[0];
sx q[0];
rz(-2.5706814) q[0];
sx q[0];
rz(1.1263322) q[0];
rz(-pi) q[1];
rz(-2.9179108) q[2];
sx q[2];
rz(-2.0111901) q[2];
sx q[2];
rz(-0.5612095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83460036) q[1];
sx q[1];
rz(-1.1728334) q[1];
sx q[1];
rz(0.88944056) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18078928) q[3];
sx q[3];
rz(-0.76813625) q[3];
sx q[3];
rz(3.1282792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3964316) q[2];
sx q[2];
rz(-0.67239422) q[2];
sx q[2];
rz(0.82610899) q[2];
rz(-2.7769026) q[3];
sx q[3];
rz(-1.6148022) q[3];
sx q[3];
rz(2.8472692) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9593338) q[0];
sx q[0];
rz(-1.2428281) q[0];
sx q[0];
rz(0.44573319) q[0];
rz(0.590473) q[1];
sx q[1];
rz(-2.1069374) q[1];
sx q[1];
rz(2.7139434) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3217309) q[0];
sx q[0];
rz(-2.68235) q[0];
sx q[0];
rz(1.4646572) q[0];
rz(-pi) q[1];
rz(-0.81645963) q[2];
sx q[2];
rz(-2.0000946) q[2];
sx q[2];
rz(0.72390899) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1688437) q[1];
sx q[1];
rz(-2.4764937) q[1];
sx q[1];
rz(1.14023) q[1];
rz(-pi) q[2];
rz(-0.0049875445) q[3];
sx q[3];
rz(-0.79646275) q[3];
sx q[3];
rz(2.4465268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6815971) q[2];
sx q[2];
rz(-2.5619016) q[2];
sx q[2];
rz(1.0350636) q[2];
rz(1.5929818) q[3];
sx q[3];
rz(-0.17954738) q[3];
sx q[3];
rz(-0.88919324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.0276133) q[0];
sx q[0];
rz(-0.005682156) q[0];
sx q[0];
rz(-0.11153829) q[0];
rz(1.5533718) q[1];
sx q[1];
rz(-0.40451834) q[1];
sx q[1];
rz(2.9389971) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83765471) q[0];
sx q[0];
rz(-1.1753774) q[0];
sx q[0];
rz(-2.9962792) q[0];
rz(-1.2738704) q[2];
sx q[2];
rz(-1.2938754) q[2];
sx q[2];
rz(-2.2001147) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21661585) q[1];
sx q[1];
rz(-1.0362715) q[1];
sx q[1];
rz(0.85304867) q[1];
x q[2];
rz(-2.5509994) q[3];
sx q[3];
rz(-1.2150798) q[3];
sx q[3];
rz(1.4416665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4488039) q[2];
sx q[2];
rz(-1.0705798) q[2];
sx q[2];
rz(-0.98133522) q[2];
rz(-2.8547309) q[3];
sx q[3];
rz(-0.57932866) q[3];
sx q[3];
rz(0.19744344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2364748) q[0];
sx q[0];
rz(-3.0906257) q[0];
sx q[0];
rz(-0.68848759) q[0];
rz(0.14432898) q[1];
sx q[1];
rz(-1.9537484) q[1];
sx q[1];
rz(-0.7483288) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7649224) q[0];
sx q[0];
rz(-0.85205305) q[0];
sx q[0];
rz(-1.5115949) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2950254) q[2];
sx q[2];
rz(-0.39381105) q[2];
sx q[2];
rz(-1.9005601) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0341558) q[1];
sx q[1];
rz(-1.7567051) q[1];
sx q[1];
rz(2.8380945) q[1];
rz(-pi) q[2];
rz(-2.2714628) q[3];
sx q[3];
rz(-1.3324454) q[3];
sx q[3];
rz(-2.7173661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.051108483) q[2];
sx q[2];
rz(-1.8277617) q[2];
sx q[2];
rz(-1.2468106) q[2];
rz(-0.11940739) q[3];
sx q[3];
rz(-1.9236919) q[3];
sx q[3];
rz(-0.30521211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7684105) q[0];
sx q[0];
rz(-0.38653448) q[0];
sx q[0];
rz(2.5266732) q[0];
rz(0.83456314) q[1];
sx q[1];
rz(-1.7465697) q[1];
sx q[1];
rz(0.83647299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8744321) q[0];
sx q[0];
rz(-1.5615441) q[0];
sx q[0];
rz(-1.6051588) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4725288) q[2];
sx q[2];
rz(-2.0537801) q[2];
sx q[2];
rz(1.3945182) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78042049) q[1];
sx q[1];
rz(-2.5675312) q[1];
sx q[1];
rz(1.7758382) q[1];
rz(-pi) q[2];
rz(2.3079268) q[3];
sx q[3];
rz(-2.1676237) q[3];
sx q[3];
rz(2.4868929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5323083) q[2];
sx q[2];
rz(-0.60695761) q[2];
sx q[2];
rz(-2.2212846) q[2];
rz(-2.3406384) q[3];
sx q[3];
rz(-2.6194173) q[3];
sx q[3];
rz(2.8029158) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9698708) q[0];
sx q[0];
rz(-1.0665749) q[0];
sx q[0];
rz(-2.875476) q[0];
rz(-1.4316106) q[1];
sx q[1];
rz(-2.0224729) q[1];
sx q[1];
rz(2.5439579) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1475567) q[0];
sx q[0];
rz(-1.7606252) q[0];
sx q[0];
rz(-2.5231383) q[0];
x q[1];
rz(0.061068717) q[2];
sx q[2];
rz(-2.960223) q[2];
sx q[2];
rz(-0.86038113) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1823765) q[1];
sx q[1];
rz(-2.2666449) q[1];
sx q[1];
rz(-1.6435502) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9868578) q[3];
sx q[3];
rz(-1.8879963) q[3];
sx q[3];
rz(-1.4829319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1282244) q[2];
sx q[2];
rz(-0.2672264) q[2];
sx q[2];
rz(-2.5426148) q[2];
rz(0.61601764) q[3];
sx q[3];
rz(-0.44121656) q[3];
sx q[3];
rz(1.0429355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1004341) q[0];
sx q[0];
rz(-2.2293595) q[0];
sx q[0];
rz(2.8225733) q[0];
rz(-3.0306385) q[1];
sx q[1];
rz(-1.3919421) q[1];
sx q[1];
rz(-2.9999733) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34401152) q[0];
sx q[0];
rz(-0.89056361) q[0];
sx q[0];
rz(-1.4466831) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0044195375) q[2];
sx q[2];
rz(-1.0194155) q[2];
sx q[2];
rz(-2.1194292) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1846022) q[1];
sx q[1];
rz(-2.5612368) q[1];
sx q[1];
rz(-2.7932211) q[1];
rz(-pi) q[2];
rz(-1.8291437) q[3];
sx q[3];
rz(-0.53505361) q[3];
sx q[3];
rz(2.759397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55613279) q[2];
sx q[2];
rz(-2.1161049) q[2];
sx q[2];
rz(0.11036135) q[2];
rz(2.6335671) q[3];
sx q[3];
rz(-0.042162687) q[3];
sx q[3];
rz(1.149811) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0069649) q[0];
sx q[0];
rz(-0.65161324) q[0];
sx q[0];
rz(1.9223258) q[0];
rz(-0.02267516) q[1];
sx q[1];
rz(-0.89212787) q[1];
sx q[1];
rz(0.55331826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5641877) q[0];
sx q[0];
rz(-0.66364849) q[0];
sx q[0];
rz(-1.7148561) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0108286) q[2];
sx q[2];
rz(-2.4891395) q[2];
sx q[2];
rz(-0.7459695) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5071514) q[1];
sx q[1];
rz(-0.63903016) q[1];
sx q[1];
rz(2.5386058) q[1];
rz(-1.2723921) q[3];
sx q[3];
rz(-1.8074028) q[3];
sx q[3];
rz(2.2158156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7360709) q[2];
sx q[2];
rz(-2.5260479) q[2];
sx q[2];
rz(1.7919398) q[2];
rz(3.0810629) q[3];
sx q[3];
rz(-2.2318201) q[3];
sx q[3];
rz(2.5133666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8796006) q[0];
sx q[0];
rz(-2.1475726) q[0];
sx q[0];
rz(0.0066198786) q[0];
rz(2.5026542) q[1];
sx q[1];
rz(-0.21955755) q[1];
sx q[1];
rz(0.42983291) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13794261) q[0];
sx q[0];
rz(-1.2009504) q[0];
sx q[0];
rz(-1.957264) q[0];
x q[1];
rz(-2.0833932) q[2];
sx q[2];
rz(-2.655003) q[2];
sx q[2];
rz(1.8712107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.51825) q[1];
sx q[1];
rz(-2.3584705) q[1];
sx q[1];
rz(3.1108969) q[1];
x q[2];
rz(-0.58542975) q[3];
sx q[3];
rz(-1.5133891) q[3];
sx q[3];
rz(-1.2186183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2868353) q[2];
sx q[2];
rz(-1.9660051) q[2];
sx q[2];
rz(3.1268934) q[2];
rz(1.1356575) q[3];
sx q[3];
rz(-2.0292) q[3];
sx q[3];
rz(1.8358102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056674615) q[0];
sx q[0];
rz(-2.882353) q[0];
sx q[0];
rz(-1.2925451) q[0];
rz(-2.9855285) q[1];
sx q[1];
rz(-0.81764692) q[1];
sx q[1];
rz(2.3110716) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1500524) q[0];
sx q[0];
rz(-0.73313289) q[0];
sx q[0];
rz(1.7197648) q[0];
x q[1];
rz(1.8761329) q[2];
sx q[2];
rz(-0.80080253) q[2];
sx q[2];
rz(-2.3126471) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3780508) q[1];
sx q[1];
rz(-1.6507859) q[1];
sx q[1];
rz(-1.7879809) q[1];
rz(-pi) q[2];
rz(1.3545348) q[3];
sx q[3];
rz(-1.2118846) q[3];
sx q[3];
rz(-0.87417904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0326651) q[2];
sx q[2];
rz(-0.69643164) q[2];
sx q[2];
rz(2.6984974) q[2];
rz(-2.9923934) q[3];
sx q[3];
rz(-0.34330338) q[3];
sx q[3];
rz(-0.33299115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14015848) q[0];
sx q[0];
rz(-2.5509111) q[0];
sx q[0];
rz(2.1668707) q[0];
rz(0.99826605) q[1];
sx q[1];
rz(-1.9010192) q[1];
sx q[1];
rz(-0.98767282) q[1];
rz(2.6218565) q[2];
sx q[2];
rz(-2.7512641) q[2];
sx q[2];
rz(0.47533585) q[2];
rz(2.0988437) q[3];
sx q[3];
rz(-0.73560148) q[3];
sx q[3];
rz(-0.92675496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
