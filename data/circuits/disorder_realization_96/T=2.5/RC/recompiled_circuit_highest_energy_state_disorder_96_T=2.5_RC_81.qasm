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
rz(1.9198298) q[0];
sx q[0];
rz(-1.3718995) q[0];
sx q[0];
rz(-2.0576117) q[0];
rz(-0.78217512) q[1];
sx q[1];
rz(-2.6922234) q[1];
sx q[1];
rz(0.78701293) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5840658) q[0];
sx q[0];
rz(-1.1827994) q[0];
sx q[0];
rz(-0.16325133) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.044620847) q[2];
sx q[2];
rz(-0.89552829) q[2];
sx q[2];
rz(-1.9683128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1842332) q[1];
sx q[1];
rz(-2.5788947) q[1];
sx q[1];
rz(-1.8398647) q[1];
x q[2];
rz(-1.161133) q[3];
sx q[3];
rz(-2.4142263) q[3];
sx q[3];
rz(-1.1687708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7996173) q[2];
sx q[2];
rz(-1.819333) q[2];
sx q[2];
rz(-0.70023099) q[2];
rz(-1.3618943) q[3];
sx q[3];
rz(-1.2075862) q[3];
sx q[3];
rz(-3.0119058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48548651) q[0];
sx q[0];
rz(-0.34657297) q[0];
sx q[0];
rz(-1.19278) q[0];
rz(0.15395173) q[1];
sx q[1];
rz(-0.62869453) q[1];
sx q[1];
rz(-1.8923632) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.10478) q[0];
sx q[0];
rz(-1.6333155) q[0];
sx q[0];
rz(0.17351242) q[0];
rz(-0.93358139) q[2];
sx q[2];
rz(-1.6145448) q[2];
sx q[2];
rz(1.9160567) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4255562) q[1];
sx q[1];
rz(-1.2662983) q[1];
sx q[1];
rz(1.6554805) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1129679) q[3];
sx q[3];
rz(-0.71958032) q[3];
sx q[3];
rz(2.6828669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93231702) q[2];
sx q[2];
rz(-1.9220158) q[2];
sx q[2];
rz(-2.2323214) q[2];
rz(-0.13898177) q[3];
sx q[3];
rz(-0.2551955) q[3];
sx q[3];
rz(2.2709258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.22055498) q[0];
sx q[0];
rz(-1.166936) q[0];
sx q[0];
rz(-1.9535109) q[0];
rz(2.2782169) q[1];
sx q[1];
rz(-1.8977576) q[1];
sx q[1];
rz(-2.1720502) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.506451) q[0];
sx q[0];
rz(-1.7291082) q[0];
sx q[0];
rz(0.16198762) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1968422) q[2];
sx q[2];
rz(-1.8965169) q[2];
sx q[2];
rz(2.0545127) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2565793) q[1];
sx q[1];
rz(-1.6672214) q[1];
sx q[1];
rz(0.23553358) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2754955) q[3];
sx q[3];
rz(-1.1148165) q[3];
sx q[3];
rz(-2.7957819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.36677507) q[2];
sx q[2];
rz(-2.884951) q[2];
sx q[2];
rz(0.65566629) q[2];
rz(1.6171148) q[3];
sx q[3];
rz(-1.7888125) q[3];
sx q[3];
rz(-2.2881983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8758288) q[0];
sx q[0];
rz(-1.1273071) q[0];
sx q[0];
rz(1.458459) q[0];
rz(-1.726285) q[1];
sx q[1];
rz(-1.7067319) q[1];
sx q[1];
rz(1.8728135) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9139713) q[0];
sx q[0];
rz(-0.73621553) q[0];
sx q[0];
rz(0.52991863) q[0];
rz(2.2591296) q[2];
sx q[2];
rz(-2.3849786) q[2];
sx q[2];
rz(0.34230912) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8377338) q[1];
sx q[1];
rz(-1.3875675) q[1];
sx q[1];
rz(-2.2732123) q[1];
rz(-pi) q[2];
rz(-1.9400408) q[3];
sx q[3];
rz(-0.24414586) q[3];
sx q[3];
rz(0.63318726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4844674) q[2];
sx q[2];
rz(-1.3544269) q[2];
sx q[2];
rz(-1.2628868) q[2];
rz(2.9386988) q[3];
sx q[3];
rz(-2.1094567) q[3];
sx q[3];
rz(-1.8111022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6577067) q[0];
sx q[0];
rz(-1.0840451) q[0];
sx q[0];
rz(-0.44802353) q[0];
rz(1.1605284) q[1];
sx q[1];
rz(-2.0612165) q[1];
sx q[1];
rz(-2.0652658) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3884013) q[0];
sx q[0];
rz(-1.7832169) q[0];
sx q[0];
rz(-2.302049) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14400173) q[2];
sx q[2];
rz(-0.33641923) q[2];
sx q[2];
rz(2.9478204) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68491918) q[1];
sx q[1];
rz(-3.0701048) q[1];
sx q[1];
rz(1.0409357) q[1];
rz(0.020298941) q[3];
sx q[3];
rz(-0.95385984) q[3];
sx q[3];
rz(0.20245173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6059604) q[2];
sx q[2];
rz(-1.863402) q[2];
sx q[2];
rz(-0.34206259) q[2];
rz(2.0096807) q[3];
sx q[3];
rz(-2.0703546) q[3];
sx q[3];
rz(-0.38558495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015942052) q[0];
sx q[0];
rz(-2.5378939) q[0];
sx q[0];
rz(-1.1018671) q[0];
rz(0.28533882) q[1];
sx q[1];
rz(-0.97831786) q[1];
sx q[1];
rz(-0.91011059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1130331) q[0];
sx q[0];
rz(-0.1814778) q[0];
sx q[0];
rz(-0.89361287) q[0];
x q[1];
rz(0.093488113) q[2];
sx q[2];
rz(-1.6851236) q[2];
sx q[2];
rz(-1.0766407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.460798) q[1];
sx q[1];
rz(-0.44728794) q[1];
sx q[1];
rz(1.7030925) q[1];
rz(-2.9858382) q[3];
sx q[3];
rz(-1.7475583) q[3];
sx q[3];
rz(-2.8946517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82112271) q[2];
sx q[2];
rz(-1.7197101) q[2];
sx q[2];
rz(2.087743) q[2];
rz(2.7779135) q[3];
sx q[3];
rz(-1.6561008) q[3];
sx q[3];
rz(-2.6065839) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66154552) q[0];
sx q[0];
rz(-2.1286025) q[0];
sx q[0];
rz(0.15740982) q[0];
rz(-0.62955457) q[1];
sx q[1];
rz(-1.5903571) q[1];
sx q[1];
rz(-3.0622283) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12211563) q[0];
sx q[0];
rz(-1.6412897) q[0];
sx q[0];
rz(1.6623783) q[0];
x q[1];
rz(0.39418295) q[2];
sx q[2];
rz(-1.7346343) q[2];
sx q[2];
rz(-0.25111408) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5461681) q[1];
sx q[1];
rz(-1.079266) q[1];
sx q[1];
rz(1.7173871) q[1];
rz(-pi) q[2];
rz(-1.6116051) q[3];
sx q[3];
rz(-2.5548078) q[3];
sx q[3];
rz(-0.87848488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6307512) q[2];
sx q[2];
rz(-1.8264822) q[2];
sx q[2];
rz(0.17244478) q[2];
rz(2.4991961) q[3];
sx q[3];
rz(-0.96697092) q[3];
sx q[3];
rz(2.7914458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.073535) q[0];
sx q[0];
rz(-1.8386766) q[0];
sx q[0];
rz(2.5750343) q[0];
rz(2.9134275) q[1];
sx q[1];
rz(-1.6888432) q[1];
sx q[1];
rz(-1.3120922) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24893783) q[0];
sx q[0];
rz(-0.83278197) q[0];
sx q[0];
rz(-0.33794649) q[0];
rz(-3.0242537) q[2];
sx q[2];
rz(-0.84716958) q[2];
sx q[2];
rz(0.1186419) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58103463) q[1];
sx q[1];
rz(-1.0315064) q[1];
sx q[1];
rz(-2.5453387) q[1];
rz(1.4890461) q[3];
sx q[3];
rz(-2.7513732) q[3];
sx q[3];
rz(-1.3345598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18099004) q[2];
sx q[2];
rz(-1.484551) q[2];
sx q[2];
rz(2.9385938) q[2];
rz(0.8391884) q[3];
sx q[3];
rz(-2.8036717) q[3];
sx q[3];
rz(0.85116974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39076552) q[0];
sx q[0];
rz(-0.44742328) q[0];
sx q[0];
rz(-2.9946193) q[0];
rz(-2.7087063) q[1];
sx q[1];
rz(-1.1914445) q[1];
sx q[1];
rz(-1.254522) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4449249) q[0];
sx q[0];
rz(-2.4913209) q[0];
sx q[0];
rz(0.5954605) q[0];
x q[1];
rz(-1.3225088) q[2];
sx q[2];
rz(-1.5380368) q[2];
sx q[2];
rz(2.0078307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.69392943) q[1];
sx q[1];
rz(-1.7433102) q[1];
sx q[1];
rz(-1.9206143) q[1];
rz(-pi) q[2];
rz(-0.16814516) q[3];
sx q[3];
rz(-1.5922308) q[3];
sx q[3];
rz(2.2529064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36048421) q[2];
sx q[2];
rz(-1.0177178) q[2];
sx q[2];
rz(2.600889) q[2];
rz(-2.8402719) q[3];
sx q[3];
rz(-2.7408528) q[3];
sx q[3];
rz(-1.8436684) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28894579) q[0];
sx q[0];
rz(-1.1971373) q[0];
sx q[0];
rz(-2.4135015) q[0];
rz(2.6241265) q[1];
sx q[1];
rz(-1.4897852) q[1];
sx q[1];
rz(0.99658406) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7383125) q[0];
sx q[0];
rz(-1.2265674) q[0];
sx q[0];
rz(-3.1136373) q[0];
rz(-0.15309454) q[2];
sx q[2];
rz(-0.587112) q[2];
sx q[2];
rz(-2.401641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7000632) q[1];
sx q[1];
rz(-0.25694381) q[1];
sx q[1];
rz(-1.5154626) q[1];
rz(-pi) q[2];
rz(0.18325652) q[3];
sx q[3];
rz(-1.6695287) q[3];
sx q[3];
rz(1.9016641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95801282) q[2];
sx q[2];
rz(-1.5637584) q[2];
sx q[2];
rz(3.0713522) q[2];
rz(3.1349414) q[3];
sx q[3];
rz(-2.9700322) q[3];
sx q[3];
rz(0.10449617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2086647) q[0];
sx q[0];
rz(-1.4069955) q[0];
sx q[0];
rz(1.7250742) q[0];
rz(-0.76864645) q[1];
sx q[1];
rz(-1.9192764) q[1];
sx q[1];
rz(2.8142014) q[1];
rz(-1.458039) q[2];
sx q[2];
rz(-1.5509477) q[2];
sx q[2];
rz(1.0890065) q[2];
rz(0.27003084) q[3];
sx q[3];
rz(-1.9237917) q[3];
sx q[3];
rz(-1.6648116) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
