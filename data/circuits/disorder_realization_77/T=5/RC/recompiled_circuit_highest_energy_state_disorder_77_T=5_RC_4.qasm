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
rz(2.0692985) q[0];
sx q[0];
rz(3.3252) q[0];
sx q[0];
rz(12.73458) q[0];
rz(1.8311365) q[1];
sx q[1];
rz(-1.1163196) q[1];
sx q[1];
rz(1.1180374) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2471825) q[0];
sx q[0];
rz(-1.4617111) q[0];
sx q[0];
rz(2.8139958) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9498212) q[2];
sx q[2];
rz(-0.91345913) q[2];
sx q[2];
rz(-2.71703) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.11169) q[1];
sx q[1];
rz(-2.3808834) q[1];
sx q[1];
rz(-3.1055443) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1088896) q[3];
sx q[3];
rz(-0.30131626) q[3];
sx q[3];
rz(2.1010984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2844124) q[2];
sx q[2];
rz(-0.99049157) q[2];
sx q[2];
rz(1.277479) q[2];
rz(2.5340581) q[3];
sx q[3];
rz(-1.7888864) q[3];
sx q[3];
rz(1.7551306) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69433895) q[0];
sx q[0];
rz(-2.4632833) q[0];
sx q[0];
rz(0.69764486) q[0];
rz(2.8088226) q[1];
sx q[1];
rz(-1.1079301) q[1];
sx q[1];
rz(-0.3322126) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2282917) q[0];
sx q[0];
rz(-1.6915404) q[0];
sx q[0];
rz(3.0101454) q[0];
x q[1];
rz(-1.2142014) q[2];
sx q[2];
rz(-2.5348672) q[2];
sx q[2];
rz(1.0068829) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7113641) q[1];
sx q[1];
rz(-1.1540742) q[1];
sx q[1];
rz(0.27208379) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4139514) q[3];
sx q[3];
rz(-1.8944358) q[3];
sx q[3];
rz(-2.0726085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4405926) q[2];
sx q[2];
rz(-1.168246) q[2];
sx q[2];
rz(-2.8458703) q[2];
rz(-3.131007) q[3];
sx q[3];
rz(-0.9897832) q[3];
sx q[3];
rz(2.4081965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6538786) q[0];
sx q[0];
rz(-2.6831477) q[0];
sx q[0];
rz(0.53471765) q[0];
rz(-1.1010822) q[1];
sx q[1];
rz(-1.4357166) q[1];
sx q[1];
rz(2.2284177) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74643) q[0];
sx q[0];
rz(-0.60747738) q[0];
sx q[0];
rz(2.1588401) q[0];
rz(1.0139364) q[2];
sx q[2];
rz(-1.7903613) q[2];
sx q[2];
rz(-2.1923421) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40620921) q[1];
sx q[1];
rz(-0.92572955) q[1];
sx q[1];
rz(-2.2570133) q[1];
rz(-1.4891154) q[3];
sx q[3];
rz(-0.73058999) q[3];
sx q[3];
rz(-3.0944892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5494626) q[2];
sx q[2];
rz(-0.89207804) q[2];
sx q[2];
rz(0.38258749) q[2];
rz(1.4866746) q[3];
sx q[3];
rz(-2.8463709) q[3];
sx q[3];
rz(2.2216643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.611068) q[0];
sx q[0];
rz(-1.684573) q[0];
sx q[0];
rz(-2.0953505) q[0];
rz(2.5776082) q[1];
sx q[1];
rz(-2.1395186) q[1];
sx q[1];
rz(-1.8246338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4119349) q[0];
sx q[0];
rz(-2.1163116) q[0];
sx q[0];
rz(-2.3778937) q[0];
rz(-pi) q[1];
rz(2.2235419) q[2];
sx q[2];
rz(-1.3155283) q[2];
sx q[2];
rz(-1.7328615) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35933381) q[1];
sx q[1];
rz(-0.31681199) q[1];
sx q[1];
rz(-0.13765361) q[1];
rz(2.7183771) q[3];
sx q[3];
rz(-0.82941705) q[3];
sx q[3];
rz(-2.2256782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4468533) q[2];
sx q[2];
rz(-1.4444192) q[2];
sx q[2];
rz(-2.1430338) q[2];
rz(-0.57991943) q[3];
sx q[3];
rz(-1.6603575) q[3];
sx q[3];
rz(-1.8370321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9272083) q[0];
sx q[0];
rz(-1.2283607) q[0];
sx q[0];
rz(2.4217915) q[0];
rz(-2.39847) q[1];
sx q[1];
rz(-1.9235976) q[1];
sx q[1];
rz(3.0845089) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6030977) q[0];
sx q[0];
rz(-2.7481226) q[0];
sx q[0];
rz(-2.7843977) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75760773) q[2];
sx q[2];
rz(-1.3146244) q[2];
sx q[2];
rz(-1.3776101) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.47810208) q[1];
sx q[1];
rz(-1.851972) q[1];
sx q[1];
rz(1.4040787) q[1];
rz(2.1716398) q[3];
sx q[3];
rz(-0.68884736) q[3];
sx q[3];
rz(0.057202489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4382978) q[2];
sx q[2];
rz(-1.9652099) q[2];
sx q[2];
rz(1.2698184) q[2];
rz(-0.54369175) q[3];
sx q[3];
rz(-2.4589977) q[3];
sx q[3];
rz(1.773275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0755587) q[0];
sx q[0];
rz(-2.7664001) q[0];
sx q[0];
rz(-2.913108) q[0];
rz(-2.1560419) q[1];
sx q[1];
rz(-1.5551714) q[1];
sx q[1];
rz(1.9353345) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0096480308) q[0];
sx q[0];
rz(-1.4033269) q[0];
sx q[0];
rz(-0.33958774) q[0];
rz(-pi) q[1];
rz(0.22322519) q[2];
sx q[2];
rz(-1.8656261) q[2];
sx q[2];
rz(-2.1708174) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9282189) q[1];
sx q[1];
rz(-0.62053939) q[1];
sx q[1];
rz(-1.7804407) q[1];
rz(-1.4643741) q[3];
sx q[3];
rz(-2.1118577) q[3];
sx q[3];
rz(2.953767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0050521684) q[2];
sx q[2];
rz(-1.2753992) q[2];
sx q[2];
rz(-1.8955815) q[2];
rz(0.15240845) q[3];
sx q[3];
rz(-2.992575) q[3];
sx q[3];
rz(-0.77308956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3054672) q[0];
sx q[0];
rz(-1.9274599) q[0];
sx q[0];
rz(0.25492302) q[0];
rz(2.1615324) q[1];
sx q[1];
rz(-1.3420811) q[1];
sx q[1];
rz(-0.23475383) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2571208) q[0];
sx q[0];
rz(-2.5354374) q[0];
sx q[0];
rz(-2.0091093) q[0];
x q[1];
rz(0.17189856) q[2];
sx q[2];
rz(-1.8020639) q[2];
sx q[2];
rz(2.8794546) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77490451) q[1];
sx q[1];
rz(-1.8340381) q[1];
sx q[1];
rz(2.7145492) q[1];
rz(-1.4244377) q[3];
sx q[3];
rz(-2.4597733) q[3];
sx q[3];
rz(-2.5949929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.047082575) q[2];
sx q[2];
rz(-0.42028704) q[2];
sx q[2];
rz(1.9455998) q[2];
rz(-2.1434873) q[3];
sx q[3];
rz(-1.866303) q[3];
sx q[3];
rz(0.78817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81414476) q[0];
sx q[0];
rz(-2.271551) q[0];
sx q[0];
rz(-0.63234627) q[0];
rz(-1.5517722) q[1];
sx q[1];
rz(-1.7951782) q[1];
sx q[1];
rz(1.8428892) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6734553) q[0];
sx q[0];
rz(-2.0103665) q[0];
sx q[0];
rz(2.8643231) q[0];
rz(1.5517637) q[2];
sx q[2];
rz(-0.90523273) q[2];
sx q[2];
rz(1.9766146) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2181311) q[1];
sx q[1];
rz(-2.0773481) q[1];
sx q[1];
rz(1.523237) q[1];
x q[2];
rz(0.15141204) q[3];
sx q[3];
rz(-1.5431343) q[3];
sx q[3];
rz(-0.20229483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87045264) q[2];
sx q[2];
rz(-2.2282034) q[2];
sx q[2];
rz(0.75922981) q[2];
rz(-2.7325654) q[3];
sx q[3];
rz(-1.7274011) q[3];
sx q[3];
rz(-0.11212382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2893534) q[0];
sx q[0];
rz(-1.2989346) q[0];
sx q[0];
rz(2.7759743) q[0];
rz(-3.1386555) q[1];
sx q[1];
rz(-1.6117088) q[1];
sx q[1];
rz(1.0681577) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097314358) q[0];
sx q[0];
rz(-1.8476163) q[0];
sx q[0];
rz(2.5702049) q[0];
rz(-pi) q[1];
rz(-0.83311501) q[2];
sx q[2];
rz(-2.1016869) q[2];
sx q[2];
rz(-2.7754727) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.059243971) q[1];
sx q[1];
rz(-1.5326534) q[1];
sx q[1];
rz(-1.9583942) q[1];
rz(0.13309388) q[3];
sx q[3];
rz(-2.9381917) q[3];
sx q[3];
rz(1.668556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.7630161) q[2];
sx q[2];
rz(-2.2282659) q[2];
sx q[2];
rz(0.98706377) q[2];
rz(2.3297564) q[3];
sx q[3];
rz(-2.2631009) q[3];
sx q[3];
rz(1.2442376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0768452) q[0];
sx q[0];
rz(-2.5223314) q[0];
sx q[0];
rz(1.6774119) q[0];
rz(-0.96425104) q[1];
sx q[1];
rz(-1.3136656) q[1];
sx q[1];
rz(-1.3826694) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2428873) q[0];
sx q[0];
rz(-0.90685493) q[0];
sx q[0];
rz(-2.6580145) q[0];
rz(-pi) q[1];
rz(0.37071642) q[2];
sx q[2];
rz(-1.3138736) q[2];
sx q[2];
rz(3.1213426) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6056435) q[1];
sx q[1];
rz(-1.1472196) q[1];
sx q[1];
rz(-2.1602195) q[1];
rz(-pi) q[2];
rz(0.85282214) q[3];
sx q[3];
rz(-2.3896165) q[3];
sx q[3];
rz(-0.2096006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1796403) q[2];
sx q[2];
rz(-1.3397168) q[2];
sx q[2];
rz(-1.4947653) q[2];
rz(-2.0252114) q[3];
sx q[3];
rz(-0.79018441) q[3];
sx q[3];
rz(-2.5377929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8981001) q[0];
sx q[0];
rz(-1.6742764) q[0];
sx q[0];
rz(1.8754638) q[0];
rz(0.19337868) q[1];
sx q[1];
rz(-2.0722176) q[1];
sx q[1];
rz(-1.8373012) q[1];
rz(-0.33756145) q[2];
sx q[2];
rz(-2.3783219) q[2];
sx q[2];
rz(-1.572883) q[2];
rz(-2.6287617) q[3];
sx q[3];
rz(-0.75779946) q[3];
sx q[3];
rz(2.1661314) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
