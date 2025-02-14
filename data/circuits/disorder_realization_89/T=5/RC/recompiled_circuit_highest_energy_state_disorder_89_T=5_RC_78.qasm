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
rz(-1.2943478) q[0];
sx q[0];
rz(-1.5273153) q[0];
sx q[0];
rz(-2.9543028) q[0];
rz(-1.3297431) q[1];
sx q[1];
rz(-2.5951374) q[1];
sx q[1];
rz(1.3475013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8719604) q[0];
sx q[0];
rz(-0.56817164) q[0];
sx q[0];
rz(1.7164036) q[0];
rz(-0.37624575) q[2];
sx q[2];
rz(-1.0828138) q[2];
sx q[2];
rz(-0.03586344) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0986024) q[1];
sx q[1];
rz(-0.8892376) q[1];
sx q[1];
rz(-0.020976457) q[1];
x q[2];
rz(2.884976) q[3];
sx q[3];
rz(-1.768075) q[3];
sx q[3];
rz(2.9159604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5225141) q[2];
sx q[2];
rz(-1.1828902) q[2];
sx q[2];
rz(-2.8600233) q[2];
rz(1.7742026) q[3];
sx q[3];
rz(-1.2308246) q[3];
sx q[3];
rz(-2.863133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53681579) q[0];
sx q[0];
rz(-2.8474658) q[0];
sx q[0];
rz(1.3500704) q[0];
rz(-0.049909441) q[1];
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
rz(-1.2725213) q[0];
sx q[0];
rz(-1.6289113) q[0];
sx q[0];
rz(2.7921667) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8698385) q[2];
sx q[2];
rz(-0.56952945) q[2];
sx q[2];
rz(-0.45646748) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3765722) q[1];
sx q[1];
rz(-2.8072551) q[1];
sx q[1];
rz(-2.2052081) q[1];
rz(-pi) q[2];
rz(-0.66920993) q[3];
sx q[3];
rz(-2.4871181) q[3];
sx q[3];
rz(1.4231285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0758948) q[2];
sx q[2];
rz(-1.1827129) q[2];
sx q[2];
rz(0.069570216) q[2];
rz(0.73404297) q[3];
sx q[3];
rz(-2.3591122) q[3];
sx q[3];
rz(-2.5231584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.065539) q[0];
sx q[0];
rz(-2.9422308) q[0];
sx q[0];
rz(-0.25654909) q[0];
rz(1.3408631) q[1];
sx q[1];
rz(-2.1840265) q[1];
sx q[1];
rz(1.0440913) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0018104) q[0];
sx q[0];
rz(-0.35813791) q[0];
sx q[0];
rz(-2.0689808) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3568722) q[2];
sx q[2];
rz(-1.7905777) q[2];
sx q[2];
rz(-2.8666988) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0699745) q[1];
sx q[1];
rz(-2.721792) q[1];
sx q[1];
rz(-0.75666766) q[1];
x q[2];
rz(-1.7875015) q[3];
sx q[3];
rz(-2.5575221) q[3];
sx q[3];
rz(-1.3497373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23036817) q[2];
sx q[2];
rz(-2.0679097) q[2];
sx q[2];
rz(-0.71630859) q[2];
rz(0.45799842) q[3];
sx q[3];
rz(-1.4665946) q[3];
sx q[3];
rz(-0.23076335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3476747) q[0];
sx q[0];
rz(-2.1383998) q[0];
sx q[0];
rz(2.4804261) q[0];
rz(-1.8874774) q[1];
sx q[1];
rz(-1.5725807) q[1];
sx q[1];
rz(-1.5987781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5970478) q[0];
sx q[0];
rz(-2.8327541) q[0];
sx q[0];
rz(0.95914532) q[0];
x q[1];
rz(-1.0087246) q[2];
sx q[2];
rz(-1.1843388) q[2];
sx q[2];
rz(1.5773048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.143121) q[1];
sx q[1];
rz(-0.70618343) q[1];
sx q[1];
rz(-0.097340214) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1041252) q[3];
sx q[3];
rz(-1.9107483) q[3];
sx q[3];
rz(1.252591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8864112) q[2];
sx q[2];
rz(-0.35021439) q[2];
sx q[2];
rz(-2.0070455) q[2];
rz(-0.55401951) q[3];
sx q[3];
rz(-1.7787245) q[3];
sx q[3];
rz(-0.58486432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2968813) q[0];
sx q[0];
rz(-3.1394594) q[0];
sx q[0];
rz(1.1312436) q[0];
rz(3.0834815) q[1];
sx q[1];
rz(-1.8986214) q[1];
sx q[1];
rz(-1.4505454) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3076271) q[0];
sx q[0];
rz(-1.183032) q[0];
sx q[0];
rz(0.8967171) q[0];
x q[1];
rz(-1.8671473) q[2];
sx q[2];
rz(-0.83186281) q[2];
sx q[2];
rz(-2.8339164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.031353019) q[1];
sx q[1];
rz(-1.3340557) q[1];
sx q[1];
rz(-2.8048022) q[1];
x q[2];
rz(-0.93673076) q[3];
sx q[3];
rz(-1.2311683) q[3];
sx q[3];
rz(2.820211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45945534) q[2];
sx q[2];
rz(-0.61209279) q[2];
sx q[2];
rz(1.7876145) q[2];
rz(-2.5891384) q[3];
sx q[3];
rz(-0.8907291) q[3];
sx q[3];
rz(0.49828211) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.246493) q[0];
sx q[0];
rz(-2.1857078) q[0];
sx q[0];
rz(-2.0434875) q[0];
rz(-1.8662628) q[1];
sx q[1];
rz(-1.9018491) q[1];
sx q[1];
rz(-2.1001508) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5811036) q[0];
sx q[0];
rz(-0.96700089) q[0];
sx q[0];
rz(2.2093456) q[0];
rz(-1.5389693) q[2];
sx q[2];
rz(-1.6897481) q[2];
sx q[2];
rz(0.54357375) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66646092) q[1];
sx q[1];
rz(-1.0589927) q[1];
sx q[1];
rz(-3.0252837) q[1];
rz(-pi) q[2];
rz(-1.8919577) q[3];
sx q[3];
rz(-0.75157673) q[3];
sx q[3];
rz(-2.1134714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8424524) q[2];
sx q[2];
rz(-1.9879397) q[2];
sx q[2];
rz(0.057272043) q[2];
rz(-0.20036571) q[3];
sx q[3];
rz(-0.55557957) q[3];
sx q[3];
rz(-0.19181767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6054691) q[0];
sx q[0];
rz(-2.496688) q[0];
sx q[0];
rz(2.7591925) q[0];
rz(-2.1740225) q[1];
sx q[1];
rz(-0.41125527) q[1];
sx q[1];
rz(1.8866906) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8604831) q[0];
sx q[0];
rz(-2.5279928) q[0];
sx q[0];
rz(-0.09765082) q[0];
x q[1];
rz(2.228565) q[2];
sx q[2];
rz(-1.0464729) q[2];
sx q[2];
rz(-1.1990449) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.7130512) q[1];
sx q[1];
rz(-1.9425434) q[1];
sx q[1];
rz(-1.1357689) q[1];
x q[2];
rz(-0.27886919) q[3];
sx q[3];
rz(-0.79852102) q[3];
sx q[3];
rz(0.65031933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3249698) q[2];
sx q[2];
rz(-0.35494706) q[2];
sx q[2];
rz(-2.2897282) q[2];
rz(-2.8999117) q[3];
sx q[3];
rz(-1.207374) q[3];
sx q[3];
rz(2.2309301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8082751) q[0];
sx q[0];
rz(-2.9589544) q[0];
sx q[0];
rz(2.0727378) q[0];
rz(1.764027) q[1];
sx q[1];
rz(-0.27605468) q[1];
sx q[1];
rz(0.68881234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.903021) q[0];
sx q[0];
rz(-1.2433143) q[0];
sx q[0];
rz(-0.37974289) q[0];
rz(-pi) q[1];
rz(1.7084863) q[2];
sx q[2];
rz(-0.97239649) q[2];
sx q[2];
rz(0.56365594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.090052) q[1];
sx q[1];
rz(-1.2512815) q[1];
sx q[1];
rz(-2.7367112) q[1];
rz(-1.7433002) q[3];
sx q[3];
rz(-1.0478847) q[3];
sx q[3];
rz(0.80181087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0546725) q[2];
sx q[2];
rz(-0.72267795) q[2];
sx q[2];
rz(1.5642536) q[2];
rz(-0.65586048) q[3];
sx q[3];
rz(-1.4411996) q[3];
sx q[3];
rz(-0.76914579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61651388) q[0];
sx q[0];
rz(-2.109313) q[0];
sx q[0];
rz(0.14933625) q[0];
rz(2.0556045) q[1];
sx q[1];
rz(-0.95567742) q[1];
sx q[1];
rz(2.65926) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5457561) q[0];
sx q[0];
rz(-2.0424423) q[0];
sx q[0];
rz(-3.0455941) q[0];
x q[1];
rz(1.6012548) q[2];
sx q[2];
rz(-0.23269146) q[2];
sx q[2];
rz(0.57664492) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73585549) q[1];
sx q[1];
rz(-0.89621004) q[1];
sx q[1];
rz(-1.2184185) q[1];
x q[2];
rz(1.6685358) q[3];
sx q[3];
rz(-1.9888788) q[3];
sx q[3];
rz(3.0168777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.09717354) q[2];
sx q[2];
rz(-0.83690518) q[2];
sx q[2];
rz(-2.6824299) q[2];
rz(3.0432213) q[3];
sx q[3];
rz(-1.2854853) q[3];
sx q[3];
rz(-1.2488825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8654873) q[0];
sx q[0];
rz(-1.8806172) q[0];
sx q[0];
rz(3.0011445) q[0];
rz(1.2175995) q[1];
sx q[1];
rz(-1.1414889) q[1];
sx q[1];
rz(1.5509031) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5031932) q[0];
sx q[0];
rz(-0.84740438) q[0];
sx q[0];
rz(1.4620848) q[0];
rz(-pi) q[1];
rz(-1.7459007) q[2];
sx q[2];
rz(-1.7410472) q[2];
sx q[2];
rz(-0.71699063) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.925736) q[1];
sx q[1];
rz(-1.832215) q[1];
sx q[1];
rz(-2.5040313) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89992739) q[3];
sx q[3];
rz(-1.560671) q[3];
sx q[3];
rz(0.72521082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6259674) q[2];
sx q[2];
rz(-2.1145861) q[2];
sx q[2];
rz(-1.1930126) q[2];
rz(0.48202816) q[3];
sx q[3];
rz(-2.7218282) q[3];
sx q[3];
rz(1.0198063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1358418) q[0];
sx q[0];
rz(-2.0477722) q[0];
sx q[0];
rz(2.3070106) q[0];
rz(-2.3781378) q[1];
sx q[1];
rz(-1.3044985) q[1];
sx q[1];
rz(0.34674092) q[1];
rz(1.687526) q[2];
sx q[2];
rz(-1.5490096) q[2];
sx q[2];
rz(-1.5646554) q[2];
rz(0.59178036) q[3];
sx q[3];
rz(-2.4587678) q[3];
sx q[3];
rz(-3.126161) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
