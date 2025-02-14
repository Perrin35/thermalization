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
rz(-1.0722941) q[0];
sx q[0];
rz(-0.18360734) q[0];
sx q[0];
rz(2.973383) q[0];
rz(-1.3104562) q[1];
sx q[1];
rz(-2.0252731) q[1];
sx q[1];
rz(2.0235553) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89441012) q[0];
sx q[0];
rz(-1.4617111) q[0];
sx q[0];
rz(-0.32759687) q[0];
rz(-pi) q[1];
rz(0.69324915) q[2];
sx q[2];
rz(-1.8680671) q[2];
sx q[2];
rz(1.7566441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.11169) q[1];
sx q[1];
rz(-2.3808834) q[1];
sx q[1];
rz(0.036048338) q[1];
rz(1.1088896) q[3];
sx q[3];
rz(-2.8402764) q[3];
sx q[3];
rz(-2.1010984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85718021) q[2];
sx q[2];
rz(-2.1511011) q[2];
sx q[2];
rz(1.8641137) q[2];
rz(-2.5340581) q[3];
sx q[3];
rz(-1.3527063) q[3];
sx q[3];
rz(1.7551306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4472537) q[0];
sx q[0];
rz(-0.67830938) q[0];
sx q[0];
rz(2.4439478) q[0];
rz(2.8088226) q[1];
sx q[1];
rz(-1.1079301) q[1];
sx q[1];
rz(2.8093801) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.913301) q[0];
sx q[0];
rz(-1.4500522) q[0];
sx q[0];
rz(3.0101454) q[0];
x q[1];
rz(0.23770419) q[2];
sx q[2];
rz(-1.0071041) q[2];
sx q[2];
rz(1.4326045) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0337341) q[1];
sx q[1];
rz(-0.49328557) q[1];
sx q[1];
rz(-2.1164338) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9929797) q[3];
sx q[3];
rz(-2.2530972) q[3];
sx q[3];
rz(2.9158031) q[3];
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
rz(0.2957224) q[2];
rz(-3.131007) q[3];
sx q[3];
rz(-2.1518095) q[3];
sx q[3];
rz(0.73339614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6538786) q[0];
sx q[0];
rz(-0.45844498) q[0];
sx q[0];
rz(2.606875) q[0];
rz(-2.0405105) q[1];
sx q[1];
rz(-1.4357166) q[1];
sx q[1];
rz(0.91317493) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0771784) q[0];
sx q[0];
rz(-2.065669) q[0];
sx q[0];
rz(0.36806182) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9703254) q[2];
sx q[2];
rz(-0.59430423) q[2];
sx q[2];
rz(2.8565796) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6221546) q[1];
sx q[1];
rz(-1.039912) q[1];
sx q[1];
rz(0.77150788) q[1];
rz(-pi) q[2];
rz(-1.4891154) q[3];
sx q[3];
rz(-2.4110027) q[3];
sx q[3];
rz(3.0944892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5494626) q[2];
sx q[2];
rz(-0.89207804) q[2];
sx q[2];
rz(2.7590052) q[2];
rz(-1.4866746) q[3];
sx q[3];
rz(-2.8463709) q[3];
sx q[3];
rz(0.91992831) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.611068) q[0];
sx q[0];
rz(-1.4570197) q[0];
sx q[0];
rz(-1.0462421) q[0];
rz(0.56398448) q[1];
sx q[1];
rz(-1.0020741) q[1];
sx q[1];
rz(-1.8246338) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3377258) q[0];
sx q[0];
rz(-0.90529862) q[0];
sx q[0];
rz(-2.4212877) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8242048) q[2];
sx q[2];
rz(-2.1989951) q[2];
sx q[2];
rz(-2.788822) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6374938) q[1];
sx q[1];
rz(-1.257084) q[1];
sx q[1];
rz(-1.5258386) q[1];
rz(-pi) q[2];
rz(-0.78329194) q[3];
sx q[3];
rz(-1.2630594) q[3];
sx q[3];
rz(2.7819992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6947393) q[2];
sx q[2];
rz(-1.6971735) q[2];
sx q[2];
rz(0.99855885) q[2];
rz(-0.57991943) q[3];
sx q[3];
rz(-1.6603575) q[3];
sx q[3];
rz(-1.8370321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21438433) q[0];
sx q[0];
rz(-1.913232) q[0];
sx q[0];
rz(0.71980113) q[0];
rz(-0.7431227) q[1];
sx q[1];
rz(-1.9235976) q[1];
sx q[1];
rz(0.057083759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29963076) q[0];
sx q[0];
rz(-1.4363382) q[0];
sx q[0];
rz(2.7706783) q[0];
x q[1];
rz(-0.75760773) q[2];
sx q[2];
rz(-1.3146244) q[2];
sx q[2];
rz(-1.3776101) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.47810208) q[1];
sx q[1];
rz(-1.2896207) q[1];
sx q[1];
rz(1.4040787) q[1];
x q[2];
rz(0.96995285) q[3];
sx q[3];
rz(-2.4527453) q[3];
sx q[3];
rz(0.057202489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4382978) q[2];
sx q[2];
rz(-1.1763828) q[2];
sx q[2];
rz(1.2698184) q[2];
rz(2.5979009) q[3];
sx q[3];
rz(-2.4589977) q[3];
sx q[3];
rz(1.773275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0755587) q[0];
sx q[0];
rz(-0.37519255) q[0];
sx q[0];
rz(2.913108) q[0];
rz(-0.98555073) q[1];
sx q[1];
rz(-1.5551714) q[1];
sx q[1];
rz(1.2062581) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6199667) q[0];
sx q[0];
rz(-1.2361467) q[0];
sx q[0];
rz(1.7482032) q[0];
rz(-pi) q[1];
rz(-1.2689077) q[2];
sx q[2];
rz(-1.3573555) q[2];
sx q[2];
rz(-0.53415307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4691733) q[1];
sx q[1];
rz(-2.1757728) q[1];
sx q[1];
rz(0.14766001) q[1];
rz(-pi) q[2];
rz(-1.4643741) q[3];
sx q[3];
rz(-2.1118577) q[3];
sx q[3];
rz(2.953767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3054672) q[0];
sx q[0];
rz(-1.2141328) q[0];
sx q[0];
rz(2.8866696) q[0];
rz(-0.98006025) q[1];
sx q[1];
rz(-1.7995116) q[1];
sx q[1];
rz(0.23475383) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8229651) q[0];
sx q[0];
rz(-1.8150094) q[0];
sx q[0];
rz(-2.1313214) q[0];
rz(-pi) q[1];
rz(-1.336194) q[2];
sx q[2];
rz(-1.7380746) q[2];
sx q[2];
rz(-1.3484311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8262485) q[1];
sx q[1];
rz(-2.6442206) q[1];
sx q[1];
rz(-2.5647463) q[1];
rz(-pi) q[2];
rz(0.8942306) q[3];
sx q[3];
rz(-1.4787592) q[3];
sx q[3];
rz(-0.91023723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.047082575) q[2];
sx q[2];
rz(-0.42028704) q[2];
sx q[2];
rz(-1.9455998) q[2];
rz(0.99810537) q[3];
sx q[3];
rz(-1.2752897) q[3];
sx q[3];
rz(-0.78817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3274479) q[0];
sx q[0];
rz(-0.87004167) q[0];
sx q[0];
rz(2.5092464) q[0];
rz(-1.5898204) q[1];
sx q[1];
rz(-1.3464144) q[1];
sx q[1];
rz(-1.2987035) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0576026) q[0];
sx q[0];
rz(-2.6267536) q[0];
sx q[0];
rz(2.0979416) q[0];
x q[1];
rz(-2.4759411) q[2];
sx q[2];
rz(-1.5558262) q[2];
sx q[2];
rz(-2.7240208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2181311) q[1];
sx q[1];
rz(-1.0642446) q[1];
sx q[1];
rz(-1.6183557) q[1];
x q[2];
rz(2.9901806) q[3];
sx q[3];
rz(-1.5431343) q[3];
sx q[3];
rz(0.20229483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.27114) q[2];
sx q[2];
rz(-2.2282034) q[2];
sx q[2];
rz(2.3823628) q[2];
rz(-0.40902725) q[3];
sx q[3];
rz(-1.4141915) q[3];
sx q[3];
rz(3.0294688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8522393) q[0];
sx q[0];
rz(-1.2989346) q[0];
sx q[0];
rz(0.36561832) q[0];
rz(0.002937142) q[1];
sx q[1];
rz(-1.5298839) q[1];
sx q[1];
rz(2.073435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4941752) q[0];
sx q[0];
rz(-1.023698) q[0];
sx q[0];
rz(-1.8965333) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83311501) q[2];
sx q[2];
rz(-1.0399057) q[2];
sx q[2];
rz(2.7754727) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.60469) q[1];
sx q[1];
rz(-0.3893756) q[1];
sx q[1];
rz(-1.4701718) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13309388) q[3];
sx q[3];
rz(-0.20340098) q[3];
sx q[3];
rz(-1.4730367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.7630161) q[2];
sx q[2];
rz(-2.2282659) q[2];
sx q[2];
rz(0.98706377) q[2];
rz(-0.8118363) q[3];
sx q[3];
rz(-0.87849179) q[3];
sx q[3];
rz(1.8973551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0768452) q[0];
sx q[0];
rz(-2.5223314) q[0];
sx q[0];
rz(1.6774119) q[0];
rz(2.1773416) q[1];
sx q[1];
rz(-1.827927) q[1];
sx q[1];
rz(-1.7589232) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98507568) q[0];
sx q[0];
rz(-1.1958953) q[0];
sx q[0];
rz(2.2945754) q[0];
rz(-pi) q[1];
rz(-1.8455454) q[2];
sx q[2];
rz(-1.2128069) q[2];
sx q[2];
rz(-1.6894947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9082038) q[1];
sx q[1];
rz(-2.1022132) q[1];
sx q[1];
rz(0.49698331) q[1];
rz(-pi) q[2];
rz(-0.85282214) q[3];
sx q[3];
rz(-2.3896165) q[3];
sx q[3];
rz(-2.9319921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9619523) q[2];
sx q[2];
rz(-1.3397168) q[2];
sx q[2];
rz(1.4947653) q[2];
rz(2.0252114) q[3];
sx q[3];
rz(-2.3514082) q[3];
sx q[3];
rz(-2.5377929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8981001) q[0];
sx q[0];
rz(-1.4673163) q[0];
sx q[0];
rz(-1.2661288) q[0];
rz(2.948214) q[1];
sx q[1];
rz(-1.0693751) q[1];
sx q[1];
rz(1.3042915) q[1];
rz(1.2639574) q[2];
sx q[2];
rz(-0.86021341) q[2];
sx q[2];
rz(-2.025069) q[2];
rz(-2.4520649) q[3];
sx q[3];
rz(-1.22682) q[3];
sx q[3];
rz(0.20709909) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
