OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(0.82013446) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(3.0537002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.983842) q[0];
sx q[0];
rz(-1.112126) q[0];
sx q[0];
rz(1.0755324) q[0];
rz(-1.0432265) q[2];
sx q[2];
rz(-1.1967778) q[2];
sx q[2];
rz(-0.56127115) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.80884113) q[1];
sx q[1];
rz(-1.0665227) q[1];
sx q[1];
rz(-1.5682975) q[1];
rz(-pi) q[2];
rz(-1.8294931) q[3];
sx q[3];
rz(-1.5651349) q[3];
sx q[3];
rz(-2.2744846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0007881) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(0.56935707) q[2];
rz(1.5287483) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(-1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5434718) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(0.55364451) q[0];
rz(-1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.233261) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2376018) q[0];
sx q[0];
rz(-1.0407789) q[0];
sx q[0];
rz(-0.7941829) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9821175) q[2];
sx q[2];
rz(-2.4607686) q[2];
sx q[2];
rz(-0.28123873) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1439647) q[1];
sx q[1];
rz(-1.7768754) q[1];
sx q[1];
rz(3.1254083) q[1];
rz(-1.8096829) q[3];
sx q[3];
rz(-2.2242862) q[3];
sx q[3];
rz(-1.7042004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0040940293) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(-3.1393576) q[2];
rz(-2.3114752) q[3];
sx q[3];
rz(-0.79289645) q[3];
sx q[3];
rz(1.8921651) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24901351) q[0];
sx q[0];
rz(-2.5254288) q[0];
sx q[0];
rz(0.85154831) q[0];
rz(2.3705204) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(-0.071203701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.486858) q[0];
sx q[0];
rz(-1.3329957) q[0];
sx q[0];
rz(2.9391857) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0729012) q[2];
sx q[2];
rz(-2.6325668) q[2];
sx q[2];
rz(0.44391649) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60656602) q[1];
sx q[1];
rz(-0.92381322) q[1];
sx q[1];
rz(-1.7391298) q[1];
rz(-0.86665385) q[3];
sx q[3];
rz(-2.1979077) q[3];
sx q[3];
rz(-1.8303527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8292134) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(2.6181347) q[2];
rz(0.33189014) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(-2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.7317384) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(-0.82558924) q[0];
rz(-1.9748953) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(-1.4473787) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71890812) q[0];
sx q[0];
rz(-1.3410765) q[0];
sx q[0];
rz(0.16750383) q[0];
x q[1];
rz(-1.6214192) q[2];
sx q[2];
rz(-0.17855893) q[2];
sx q[2];
rz(-0.52390097) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73562685) q[1];
sx q[1];
rz(-1.4117129) q[1];
sx q[1];
rz(-2.2707978) q[1];
rz(-pi) q[2];
rz(-0.36376603) q[3];
sx q[3];
rz(-1.7978661) q[3];
sx q[3];
rz(-1.0432537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46431413) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(-2.9966667) q[2];
rz(1.011301) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(-3.1269126) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1458364) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(-1.1902635) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(-1.9285944) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66515231) q[0];
sx q[0];
rz(-1.9422918) q[0];
sx q[0];
rz(1.2907589) q[0];
rz(-pi) q[1];
rz(0.91765399) q[2];
sx q[2];
rz(-2.0954164) q[2];
sx q[2];
rz(-1.7078924) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2781093) q[1];
sx q[1];
rz(-0.28243318) q[1];
sx q[1];
rz(-2.7062098) q[1];
rz(1.8323684) q[3];
sx q[3];
rz(-1.3993169) q[3];
sx q[3];
rz(1.8985626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.50903901) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(2.5115013) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(-0.8852638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1155788) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(-2.8552326) q[0];
rz(-0.59965602) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(0.58475959) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.023914) q[0];
sx q[0];
rz(-0.6891793) q[0];
sx q[0];
rz(2.6373449) q[0];
x q[1];
rz(1.5969959) q[2];
sx q[2];
rz(-0.95653406) q[2];
sx q[2];
rz(-0.40528497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6962291) q[1];
sx q[1];
rz(-2.3949957) q[1];
sx q[1];
rz(1.1397051) q[1];
rz(-2.8486409) q[3];
sx q[3];
rz(-0.17360273) q[3];
sx q[3];
rz(-1.1374744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9770603) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(1.5131081) q[2];
rz(2.1785054) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(-0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(-2.5928296) q[0];
rz(-0.051963003) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(-0.60639492) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0604027) q[0];
sx q[0];
rz(-0.58433825) q[0];
sx q[0];
rz(-1.1330182) q[0];
x q[1];
rz(0.37574558) q[2];
sx q[2];
rz(-0.94259113) q[2];
sx q[2];
rz(-1.1496161) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.079271) q[1];
sx q[1];
rz(-0.81264001) q[1];
sx q[1];
rz(0.15516849) q[1];
rz(2.7860836) q[3];
sx q[3];
rz(-1.9572557) q[3];
sx q[3];
rz(-0.10458065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2018044) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-2.6331804) q[2];
rz(-1.1172179) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-0.33928076) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(2.3876277) q[0];
rz(0.94999653) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(-2.9139013) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3823711) q[0];
sx q[0];
rz(-1.641944) q[0];
sx q[0];
rz(1.9842149) q[0];
rz(-pi) q[1];
rz(-2.4661602) q[2];
sx q[2];
rz(-2.5917412) q[2];
sx q[2];
rz(2.7477802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88340064) q[1];
sx q[1];
rz(-1.0419783) q[1];
sx q[1];
rz(-1.8316815) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55898198) q[3];
sx q[3];
rz(-1.6056656) q[3];
sx q[3];
rz(-2.1502286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.078538744) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(2.3030346) q[2];
rz(2.7770384) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(-0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.39695981) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(-0.80378419) q[0];
rz(2.0869758) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(-1.1484336) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.046612) q[0];
sx q[0];
rz(-1.5946832) q[0];
sx q[0];
rz(-1.6427965) q[0];
x q[1];
rz(-2.0954275) q[2];
sx q[2];
rz(-1.8881919) q[2];
sx q[2];
rz(-1.5898926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.941274) q[1];
sx q[1];
rz(-2.4793844) q[1];
sx q[1];
rz(0.83076417) q[1];
rz(-pi) q[2];
rz(-1.5356482) q[3];
sx q[3];
rz(-1.2190378) q[3];
sx q[3];
rz(0.37946057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0917197) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(2.3633374) q[2];
rz(-0.19566472) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8675999) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(0.25319779) q[0];
rz(-2.4018438) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(0.69828066) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2525576) q[0];
sx q[0];
rz(-1.6425942) q[0];
sx q[0];
rz(0.30307146) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1168078) q[2];
sx q[2];
rz(-1.4966655) q[2];
sx q[2];
rz(-0.63400808) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89430289) q[1];
sx q[1];
rz(-1.3548046) q[1];
sx q[1];
rz(2.0386019) q[1];
x q[2];
rz(1.2020338) q[3];
sx q[3];
rz(-0.59951111) q[3];
sx q[3];
rz(2.5332019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0239821) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(0.84021604) q[2];
rz(0.64030567) q[3];
sx q[3];
rz(-2.2655723) q[3];
sx q[3];
rz(-1.9742112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.2587851) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(-0.029126833) q[1];
sx q[1];
rz(-3.0283785) q[1];
sx q[1];
rz(1.821847) q[1];
rz(0.60224709) q[2];
sx q[2];
rz(-1.2802274) q[2];
sx q[2];
rz(1.0614492) q[2];
rz(-0.15016951) q[3];
sx q[3];
rz(-0.88610813) q[3];
sx q[3];
rz(-0.025972493) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
