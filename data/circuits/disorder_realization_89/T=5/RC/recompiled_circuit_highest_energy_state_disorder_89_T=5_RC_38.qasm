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
rz(1.8472449) q[0];
sx q[0];
rz(-1.6142774) q[0];
sx q[0];
rz(-0.18728988) q[0];
rz(-1.3297431) q[1];
sx q[1];
rz(-2.5951374) q[1];
sx q[1];
rz(1.3475013) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4418877) q[0];
sx q[0];
rz(-1.0093657) q[0];
sx q[0];
rz(0.092362837) q[0];
rz(1.0521837) q[2];
sx q[2];
rz(-1.901327) q[2];
sx q[2];
rz(1.7180819) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.042990265) q[1];
sx q[1];
rz(-0.8892376) q[1];
sx q[1];
rz(3.1206162) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7745733) q[3];
sx q[3];
rz(-1.8223267) q[3];
sx q[3];
rz(-1.2937781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61907855) q[2];
sx q[2];
rz(-1.9587025) q[2];
sx q[2];
rz(0.28156933) q[2];
rz(1.36739) q[3];
sx q[3];
rz(-1.2308246) q[3];
sx q[3];
rz(2.863133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53681579) q[0];
sx q[0];
rz(-2.8474658) q[0];
sx q[0];
rz(-1.7915223) q[0];
rz(3.0916832) q[1];
sx q[1];
rz(-2.8027746) q[1];
sx q[1];
rz(-1.8972338) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8690714) q[0];
sx q[0];
rz(-1.6289113) q[0];
sx q[0];
rz(-0.34942594) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18644615) q[2];
sx q[2];
rz(-1.0294339) q[2];
sx q[2];
rz(0.10554927) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.10319732) q[1];
sx q[1];
rz(-1.8382676) q[1];
sx q[1];
rz(2.9385376) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1265321) q[3];
sx q[3];
rz(-2.0685398) q[3];
sx q[3];
rz(2.5023482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0758948) q[2];
sx q[2];
rz(-1.1827129) q[2];
sx q[2];
rz(-0.069570216) q[2];
rz(-0.73404297) q[3];
sx q[3];
rz(-0.78248048) q[3];
sx q[3];
rz(-2.5231584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0760536) q[0];
sx q[0];
rz(-0.19936182) q[0];
sx q[0];
rz(2.8850436) q[0];
rz(1.8007295) q[1];
sx q[1];
rz(-2.1840265) q[1];
sx q[1];
rz(-1.0440913) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097819177) q[0];
sx q[0];
rz(-1.4025084) q[0];
sx q[0];
rz(1.8884482) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22473904) q[2];
sx q[2];
rz(-1.7794975) q[2];
sx q[2];
rz(-1.8930184) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0716182) q[1];
sx q[1];
rz(-2.721792) q[1];
sx q[1];
rz(2.384925) q[1];
rz(-0.14117853) q[3];
sx q[3];
rz(-1.0021082) q[3];
sx q[3];
rz(-1.6077667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23036817) q[2];
sx q[2];
rz(-1.0736829) q[2];
sx q[2];
rz(-2.4252841) q[2];
rz(-2.6835942) q[3];
sx q[3];
rz(-1.4665946) q[3];
sx q[3];
rz(2.9108293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3476747) q[0];
sx q[0];
rz(-1.0031928) q[0];
sx q[0];
rz(-0.66116655) q[0];
rz(-1.8874774) q[1];
sx q[1];
rz(-1.5690119) q[1];
sx q[1];
rz(1.5987781) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7044012) q[0];
sx q[0];
rz(-1.7462303) q[0];
sx q[0];
rz(-1.3153005) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6933367) q[2];
sx q[2];
rz(-1.0545306) q[2];
sx q[2];
rz(-0.22656952) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1261166) q[1];
sx q[1];
rz(-0.86864352) q[1];
sx q[1];
rz(-1.6534991) q[1];
x q[2];
rz(-1.0374674) q[3];
sx q[3];
rz(-1.2308443) q[3];
sx q[3];
rz(-1.252591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25518146) q[2];
sx q[2];
rz(-0.35021439) q[2];
sx q[2];
rz(2.0070455) q[2];
rz(0.55401951) q[3];
sx q[3];
rz(-1.3628682) q[3];
sx q[3];
rz(-0.58486432) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8447113) q[0];
sx q[0];
rz(-3.1394594) q[0];
sx q[0];
rz(-1.1312436) q[0];
rz(3.0834815) q[1];
sx q[1];
rz(-1.8986214) q[1];
sx q[1];
rz(-1.4505454) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7052225) q[0];
sx q[0];
rz(-0.76226018) q[0];
sx q[0];
rz(-0.99135474) q[0];
rz(-pi) q[1];
rz(-0.76116107) q[2];
sx q[2];
rz(-1.3532172) q[2];
sx q[2];
rz(-2.0812931) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6213773) q[1];
sx q[1];
rz(-1.2437553) q[1];
sx q[1];
rz(1.3205297) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41330321) q[3];
sx q[3];
rz(-2.163475) q[3];
sx q[3];
rz(-1.0091708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45945534) q[2];
sx q[2];
rz(-2.5294999) q[2];
sx q[2];
rz(1.3539782) q[2];
rz(-0.55245429) q[3];
sx q[3];
rz(-2.2508636) q[3];
sx q[3];
rz(-2.6433105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.246493) q[0];
sx q[0];
rz(-0.95588481) q[0];
sx q[0];
rz(1.0981052) q[0];
rz(-1.2753298) q[1];
sx q[1];
rz(-1.9018491) q[1];
sx q[1];
rz(2.1001508) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7324243) q[0];
sx q[0];
rz(-2.083626) q[0];
sx q[0];
rz(-0.70968117) q[0];
x q[1];
rz(2.8813754) q[2];
sx q[2];
rz(-0.12311664) q[2];
sx q[2];
rz(-0.80568635) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2944156) q[1];
sx q[1];
rz(-1.4694459) q[1];
sx q[1];
rz(-1.0560929) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8547229) q[3];
sx q[3];
rz(-2.2755945) q[3];
sx q[3];
rz(1.4554086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8424524) q[2];
sx q[2];
rz(-1.9879397) q[2];
sx q[2];
rz(0.057272043) q[2];
rz(2.9412269) q[3];
sx q[3];
rz(-0.55557957) q[3];
sx q[3];
rz(-0.19181767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5361236) q[0];
sx q[0];
rz(-2.496688) q[0];
sx q[0];
rz(0.38240018) q[0];
rz(2.1740225) q[1];
sx q[1];
rz(-0.41125527) q[1];
sx q[1];
rz(1.254902) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8604831) q[0];
sx q[0];
rz(-2.5279928) q[0];
sx q[0];
rz(3.0439418) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91302769) q[2];
sx q[2];
rz(-1.0464729) q[2];
sx q[2];
rz(1.1990449) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6204024) q[1];
sx q[1];
rz(-2.5772021) q[1];
sx q[1];
rz(-0.82427967) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77882336) q[3];
sx q[3];
rz(-1.7692788) q[3];
sx q[3];
rz(2.4183111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3249698) q[2];
sx q[2];
rz(-0.35494706) q[2];
sx q[2];
rz(2.2897282) q[2];
rz(-2.8999117) q[3];
sx q[3];
rz(-1.207374) q[3];
sx q[3];
rz(2.2309301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3333176) q[0];
sx q[0];
rz(-0.18263826) q[0];
sx q[0];
rz(-2.0727378) q[0];
rz(-1.3775657) q[1];
sx q[1];
rz(-2.865538) q[1];
sx q[1];
rz(2.4527803) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6816872) q[0];
sx q[0];
rz(-1.2121768) q[0];
sx q[0];
rz(1.9214517) q[0];
rz(-pi) q[1];
x q[1];
rz(2.942932) q[2];
sx q[2];
rz(-0.61214329) q[2];
sx q[2];
rz(-2.8191301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4885284) q[1];
sx q[1];
rz(-1.9540857) q[1];
sx q[1];
rz(-1.2252818) q[1];
x q[2];
rz(1.7433002) q[3];
sx q[3];
rz(-2.093708) q[3];
sx q[3];
rz(-2.3397818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0869202) q[2];
sx q[2];
rz(-2.4189147) q[2];
sx q[2];
rz(-1.5773391) q[2];
rz(-0.65586048) q[3];
sx q[3];
rz(-1.7003931) q[3];
sx q[3];
rz(-2.3724469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5250788) q[0];
sx q[0];
rz(-1.0322796) q[0];
sx q[0];
rz(2.9922564) q[0];
rz(-2.0556045) q[1];
sx q[1];
rz(-0.95567742) q[1];
sx q[1];
rz(0.48233262) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59583658) q[0];
sx q[0];
rz(-2.0424423) q[0];
sx q[0];
rz(-3.0455941) q[0];
rz(1.8033837) q[2];
sx q[2];
rz(-1.5637737) q[2];
sx q[2];
rz(-1.0237895) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9378523) q[1];
sx q[1];
rz(-2.3934869) q[1];
sx q[1];
rz(0.40740537) q[1];
x q[2];
rz(1.6685358) q[3];
sx q[3];
rz(-1.9888788) q[3];
sx q[3];
rz(-0.12471499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0444191) q[2];
sx q[2];
rz(-0.83690518) q[2];
sx q[2];
rz(0.4591628) q[2];
rz(3.0432213) q[3];
sx q[3];
rz(-1.2854853) q[3];
sx q[3];
rz(-1.2488825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2761053) q[0];
sx q[0];
rz(-1.8806172) q[0];
sx q[0];
rz(0.14044811) q[0];
rz(1.9239931) q[1];
sx q[1];
rz(-2.0001037) q[1];
sx q[1];
rz(1.5509031) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86027788) q[0];
sx q[0];
rz(-1.6522121) q[0];
sx q[0];
rz(-2.415262) q[0];
rz(0.17284278) q[2];
sx q[2];
rz(-1.3982492) q[2];
sx q[2];
rz(2.3177528) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.019245) q[1];
sx q[1];
rz(-0.68210318) q[1];
sx q[1];
rz(0.42241272) q[1];
rz(-0.89992739) q[3];
sx q[3];
rz(-1.560671) q[3];
sx q[3];
rz(2.4163818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6259674) q[2];
sx q[2];
rz(-2.1145861) q[2];
sx q[2];
rz(-1.94858) q[2];
rz(2.6595645) q[3];
sx q[3];
rz(-0.41976443) q[3];
sx q[3];
rz(-2.1217864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(2.1358418) q[0];
sx q[0];
rz(-2.0477722) q[0];
sx q[0];
rz(2.3070106) q[0];
rz(2.3781378) q[1];
sx q[1];
rz(-1.8370942) q[1];
sx q[1];
rz(-2.7948517) q[1];
rz(3.1196567) q[2];
sx q[2];
rz(-1.4540945) q[2];
sx q[2];
rz(-3.1380063) q[2];
rz(0.59377311) q[3];
sx q[3];
rz(-1.9304921) q[3];
sx q[3];
rz(-1.0747128) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
