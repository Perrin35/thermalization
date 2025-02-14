OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.26164126396179) q[0];
sx q[0];
rz(2.99614417751367) q[0];
sx q[0];
rz(9.15691891907855) q[0];
rz(-1.56758666038513) q[1];
sx q[1];
rz(6.11463800271089) q[1];
sx q[1];
rz(10.0028789400975) q[1];
cx q[1],q[0];
rz(3.01658391952515) q[0];
sx q[0];
rz(1.84682312806184) q[0];
sx q[0];
rz(10.1191239714543) q[0];
rz(-4.0139307975769) q[2];
sx q[2];
rz(4.81670430501039) q[2];
sx q[2];
rz(11.3084531783979) q[2];
cx q[2],q[1];
rz(-0.720880270004272) q[1];
sx q[1];
rz(1.13627496560151) q[1];
sx q[1];
rz(8.97352094053432) q[1];
rz(2.43791937828064) q[3];
sx q[3];
rz(1.44385853608186) q[3];
sx q[3];
rz(7.23978922366306) q[3];
cx q[3],q[2];
rz(-0.198700398206711) q[2];
sx q[2];
rz(3.55271285970742) q[2];
sx q[2];
rz(7.94873604773685) q[2];
rz(-0.5792515873909) q[3];
sx q[3];
rz(4.29605618317659) q[3];
sx q[3];
rz(10.0907020926396) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-3.99934792518616) q[0];
sx q[0];
rz(2.09601608117158) q[0];
sx q[0];
rz(12.6215114355008) q[0];
rz(2.22712302207947) q[1];
sx q[1];
rz(1.44174590905244) q[1];
sx q[1];
rz(8.20897636412784) q[1];
cx q[1],q[0];
rz(2.89126348495483) q[0];
sx q[0];
rz(3.07606524427468) q[0];
sx q[0];
rz(8.14130625724002) q[0];
rz(0.492320835590363) q[2];
sx q[2];
rz(5.18687930901582) q[2];
sx q[2];
rz(10.4077918887059) q[2];
cx q[2],q[1];
rz(-5.46912002563477) q[1];
sx q[1];
rz(1.17209974129731) q[1];
sx q[1];
rz(15.6879472494046) q[1];
rz(0.137414142489433) q[3];
sx q[3];
rz(4.26628163655336) q[3];
sx q[3];
rz(9.72143164872333) q[3];
cx q[3],q[2];
rz(1.93322825431824) q[2];
sx q[2];
rz(2.62319818337495) q[2];
sx q[2];
rz(8.80420104264423) q[2];
rz(1.68794512748718) q[3];
sx q[3];
rz(4.17332151730592) q[3];
sx q[3];
rz(8.34778163432285) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.8715912103653) q[0];
sx q[0];
rz(3.47296878893907) q[0];
sx q[0];
rz(10.9351639509122) q[0];
rz(-0.673914670944214) q[1];
sx q[1];
rz(4.98804381688172) q[1];
sx q[1];
rz(7.90855536460086) q[1];
cx q[1],q[0];
rz(3.22979521751404) q[0];
sx q[0];
rz(3.67276737292344) q[0];
sx q[0];
rz(6.23402879237338) q[0];
rz(0.497559785842896) q[2];
sx q[2];
rz(5.58157530625398) q[2];
sx q[2];
rz(9.32945363073751) q[2];
cx q[2],q[1];
rz(-0.618058145046234) q[1];
sx q[1];
rz(8.94064250786836) q[1];
sx q[1];
rz(11.4995071649472) q[1];
rz(-1.11567759513855) q[3];
sx q[3];
rz(3.74607482750947) q[3];
sx q[3];
rz(6.9678160905759) q[3];
cx q[3],q[2];
rz(3.52660202980042) q[2];
sx q[2];
rz(4.74147322972352) q[2];
sx q[2];
rz(11.9008462190549) q[2];
rz(-0.743379831314087) q[3];
sx q[3];
rz(4.26267448266084) q[3];
sx q[3];
rz(9.65229583381816) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(4.51254892349243) q[0];
sx q[0];
rz(4.19407847722108) q[0];
sx q[0];
rz(10.7250301599424) q[0];
rz(-0.119967319071293) q[1];
sx q[1];
rz(4.22213867505128) q[1];
sx q[1];
rz(8.37803289889499) q[1];
cx q[1],q[0];
rz(0.204897314310074) q[0];
sx q[0];
rz(6.29539600213105) q[0];
sx q[0];
rz(12.3320994138639) q[0];
rz(0.414183706045151) q[2];
sx q[2];
rz(2.3799171765619) q[2];
sx q[2];
rz(7.98792479037448) q[2];
cx q[2],q[1];
rz(1.50685977935791) q[1];
sx q[1];
rz(4.76748636563356) q[1];
sx q[1];
rz(9.91139051913425) q[1];
rz(-2.01882123947144) q[3];
sx q[3];
rz(5.25445452530915) q[3];
sx q[3];
rz(6.04152867793247) q[3];
cx q[3],q[2];
rz(-1.34694242477417) q[2];
sx q[2];
rz(4.20349160035188) q[2];
sx q[2];
rz(14.9364061117093) q[2];
rz(-1.05093419551849) q[3];
sx q[3];
rz(4.17509749730165) q[3];
sx q[3];
rz(8.23138782977268) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.134503841400146) q[0];
sx q[0];
rz(4.83518114884431) q[0];
sx q[0];
rz(8.61949328183337) q[0];
rz(-1.44364202022552) q[1];
sx q[1];
rz(2.3256358226114) q[1];
sx q[1];
rz(9.46124137043163) q[1];
cx q[1],q[0];
rz(2.54924607276917) q[0];
sx q[0];
rz(1.00394645531709) q[0];
sx q[0];
rz(7.79667375086948) q[0];
rz(2.57339525222778) q[2];
sx q[2];
rz(4.65468910534913) q[2];
sx q[2];
rz(7.23215363024875) q[2];
cx q[2],q[1];
rz(3.37373900413513) q[1];
sx q[1];
rz(4.65922358830506) q[1];
sx q[1];
rz(6.76055786608859) q[1];
rz(1.32124722003937) q[3];
sx q[3];
rz(6.87565437157685) q[3];
sx q[3];
rz(8.40998539923831) q[3];
cx q[3],q[2];
rz(-2.52638506889343) q[2];
sx q[2];
rz(4.26630905468995) q[2];
sx q[2];
rz(10.0014012217443) q[2];
rz(-1.54551005363464) q[3];
sx q[3];
rz(2.28710642655427) q[3];
sx q[3];
rz(8.84790012835666) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.110351897776127) q[0];
sx q[0];
rz(4.62915817101533) q[0];
sx q[0];
rz(10.0245737194936) q[0];
rz(2.33926296234131) q[1];
sx q[1];
rz(5.1914440711313) q[1];
sx q[1];
rz(6.93101522921726) q[1];
cx q[1],q[0];
rz(1.2596652507782) q[0];
sx q[0];
rz(0.818109424906321) q[0];
sx q[0];
rz(11.9415876626889) q[0];
rz(6.31605386734009) q[2];
sx q[2];
rz(0.210156591730662) q[2];
sx q[2];
rz(7.88009951113864) q[2];
cx q[2],q[1];
rz(-3.17761826515198) q[1];
sx q[1];
rz(8.75401035149629) q[1];
sx q[1];
rz(8.34098765849277) q[1];
rz(2.30238890647888) q[3];
sx q[3];
rz(5.43619075615937) q[3];
sx q[3];
rz(9.81150690316364) q[3];
cx q[3],q[2];
rz(-0.997377932071686) q[2];
sx q[2];
rz(3.51472941239411) q[2];
sx q[2];
rz(11.4691440820615) q[2];
rz(-2.14134740829468) q[3];
sx q[3];
rz(2.21792432864244) q[3];
sx q[3];
rz(11.2606590747754) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.03449869155884) q[0];
sx q[0];
rz(4.29602900345857) q[0];
sx q[0];
rz(12.7655634641568) q[0];
rz(1.84326028823853) q[1];
sx q[1];
rz(3.50270891387994) q[1];
sx q[1];
rz(10.0668200015943) q[1];
cx q[1],q[0];
rz(0.108109712600708) q[0];
sx q[0];
rz(2.38315215905244) q[0];
sx q[0];
rz(9.70152417420551) q[0];
rz(3.66942739486694) q[2];
sx q[2];
rz(5.0758000930124) q[2];
sx q[2];
rz(9.37236250414654) q[2];
cx q[2],q[1];
rz(0.0470280013978481) q[1];
sx q[1];
rz(0.215122612314769) q[1];
sx q[1];
rz(7.46925184725925) q[1];
rz(-0.45002955198288) q[3];
sx q[3];
rz(4.89511838753755) q[3];
sx q[3];
rz(9.16624791025325) q[3];
cx q[3],q[2];
rz(0.228849977254868) q[2];
sx q[2];
rz(4.20727637608583) q[2];
sx q[2];
rz(7.45749685763522) q[2];
rz(-0.922853887081146) q[3];
sx q[3];
rz(3.57920426328714) q[3];
sx q[3];
rz(9.25547370909854) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.66284239292145) q[0];
sx q[0];
rz(4.62324705918367) q[0];
sx q[0];
rz(10.42394307851) q[0];
rz(0.83813464641571) q[1];
sx q[1];
rz(2.23313877184922) q[1];
sx q[1];
rz(10.450322008125) q[1];
cx q[1],q[0];
rz(1.98601722717285) q[0];
sx q[0];
rz(3.97442838748033) q[0];
sx q[0];
rz(10.0823409914891) q[0];
rz(2.88372206687927) q[2];
sx q[2];
rz(4.96553030808503) q[2];
sx q[2];
rz(17.9115533590238) q[2];
cx q[2],q[1];
rz(1.22378122806549) q[1];
sx q[1];
rz(7.46125522454316) q[1];
sx q[1];
rz(6.81948325633212) q[1];
rz(2.28979229927063) q[3];
sx q[3];
rz(5.13490513165528) q[3];
sx q[3];
rz(9.05501473545238) q[3];
cx q[3],q[2];
rz(3.12135577201843) q[2];
sx q[2];
rz(1.78422239621217) q[2];
sx q[2];
rz(9.64197041689559) q[2];
rz(-2.00026655197144) q[3];
sx q[3];
rz(3.5165825506025) q[3];
sx q[3];
rz(8.50962898730441) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.71042573451996) q[0];
sx q[0];
rz(4.43581870396669) q[0];
sx q[0];
rz(9.44898744522735) q[0];
rz(-1.09120333194733) q[1];
sx q[1];
rz(1.11874357064302) q[1];
sx q[1];
rz(10.2388543248098) q[1];
cx q[1],q[0];
rz(-0.592107892036438) q[0];
sx q[0];
rz(3.85744163592393) q[0];
sx q[0];
rz(11.1264476537625) q[0];
rz(0.700687527656555) q[2];
sx q[2];
rz(3.89716974099214) q[2];
sx q[2];
rz(12.9016668558042) q[2];
cx q[2],q[1];
rz(-4.33545827865601) q[1];
sx q[1];
rz(4.71366146405274) q[1];
sx q[1];
rz(8.66588840483829) q[1];
rz(-0.357528805732727) q[3];
sx q[3];
rz(3.38487889071042) q[3];
sx q[3];
rz(11.099561905853) q[3];
cx q[3],q[2];
rz(-2.7476863861084) q[2];
sx q[2];
rz(4.00363621314103) q[2];
sx q[2];
rz(14.0964298009793) q[2];
rz(4.1927433013916) q[3];
sx q[3];
rz(8.0839001258188) q[3];
sx q[3];
rz(9.52625053971216) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.53973388671875) q[0];
sx q[0];
rz(4.1325468142801) q[0];
sx q[0];
rz(9.90193799733325) q[0];
rz(1.59021437168121) q[1];
sx q[1];
rz(4.80421462853486) q[1];
sx q[1];
rz(11.2593155860822) q[1];
cx q[1],q[0];
rz(2.40687489509583) q[0];
sx q[0];
rz(3.99715301592881) q[0];
sx q[0];
rz(11.0555875062863) q[0];
rz(2.38081955909729) q[2];
sx q[2];
rz(-0.562134353322438) q[2];
sx q[2];
rz(9.29105967878505) q[2];
cx q[2],q[1];
rz(3.30227065086365) q[1];
sx q[1];
rz(4.18749192555482) q[1];
sx q[1];
rz(7.15763232707187) q[1];
rz(1.9973646402359) q[3];
sx q[3];
rz(5.27678099473054) q[3];
sx q[3];
rz(13.6697568654935) q[3];
cx q[3],q[2];
rz(9.59083271026611) q[2];
sx q[2];
rz(7.29878011544282) q[2];
sx q[2];
rz(11.8139931917112) q[2];
rz(3.25423240661621) q[3];
sx q[3];
rz(6.10930863221223) q[3];
sx q[3];
rz(8.26203272341892) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.22396290302277) q[0];
sx q[0];
rz(3.61979061563546) q[0];
sx q[0];
rz(5.8702659368436) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.0351939201355) q[1];
sx q[1];
rz(0.57065740426118) q[1];
sx q[1];
rz(6.39096615313693) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(4.17224931716919) q[2];
sx q[2];
rz(6.66655317147309) q[2];
sx q[2];
rz(7.7646785736005) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.685580134391785) q[3];
sx q[3];
rz(5.35194841225679) q[3];
sx q[3];
rz(7.41265962123081) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
