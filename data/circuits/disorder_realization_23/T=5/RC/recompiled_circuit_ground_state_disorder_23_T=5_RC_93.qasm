OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57103676) q[0];
sx q[0];
rz(3.8138226) q[0];
sx q[0];
rz(11.091118) q[0];
rz(0.9530468) q[1];
sx q[1];
rz(3.3068125) q[1];
sx q[1];
rz(10.693476) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30820134) q[0];
sx q[0];
rz(-1.1723639) q[0];
sx q[0];
rz(1.3223668) q[0];
rz(-2.2242429) q[2];
sx q[2];
rz(-1.9731865) q[2];
sx q[2];
rz(1.0831329) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0787934) q[1];
sx q[1];
rz(-1.4401762) q[1];
sx q[1];
rz(2.1427633) q[1];
rz(2.8976909) q[3];
sx q[3];
rz(-0.77195215) q[3];
sx q[3];
rz(-0.011079196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1573726) q[2];
sx q[2];
rz(-2.7998388) q[2];
sx q[2];
rz(2.3483707) q[2];
rz(-2.4685517) q[3];
sx q[3];
rz(-2.0561736) q[3];
sx q[3];
rz(-1.8719155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9840045) q[0];
sx q[0];
rz(-2.3389811) q[0];
sx q[0];
rz(-0.60930914) q[0];
rz(-0.20978236) q[1];
sx q[1];
rz(-2.3502626) q[1];
sx q[1];
rz(0.9598859) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578589) q[0];
sx q[0];
rz(-1.4918033) q[0];
sx q[0];
rz(-1.5261488) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0591956) q[2];
sx q[2];
rz(-1.6325743) q[2];
sx q[2];
rz(0.90324963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9896796) q[1];
sx q[1];
rz(-1.6335347) q[1];
sx q[1];
rz(-2.3082118) q[1];
rz(-0.87374346) q[3];
sx q[3];
rz(-0.62752073) q[3];
sx q[3];
rz(0.88117679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3557055) q[2];
sx q[2];
rz(-2.3023119) q[2];
sx q[2];
rz(2.6500224) q[2];
rz(-0.0056886557) q[3];
sx q[3];
rz(-2.0523043) q[3];
sx q[3];
rz(-2.6277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095057644) q[0];
sx q[0];
rz(-2.6061366) q[0];
sx q[0];
rz(2.3989578) q[0];
rz(-0.62478089) q[1];
sx q[1];
rz(-2.2014328) q[1];
sx q[1];
rz(2.0754441) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10342205) q[0];
sx q[0];
rz(-0.26525149) q[0];
sx q[0];
rz(1.0031423) q[0];
rz(0.19045388) q[2];
sx q[2];
rz(-1.2553658) q[2];
sx q[2];
rz(-2.4078232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.448145) q[1];
sx q[1];
rz(-1.9733917) q[1];
sx q[1];
rz(-0.13685302) q[1];
rz(-pi) q[2];
rz(-2.0622018) q[3];
sx q[3];
rz(-1.3771476) q[3];
sx q[3];
rz(2.1313063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2523969) q[2];
sx q[2];
rz(-0.19813457) q[2];
sx q[2];
rz(-0.55257094) q[2];
rz(-2.6342454) q[3];
sx q[3];
rz(-1.0891958) q[3];
sx q[3];
rz(-0.56940091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4813389) q[0];
sx q[0];
rz(-0.57732552) q[0];
sx q[0];
rz(1.8950155) q[0];
rz(1.735911) q[1];
sx q[1];
rz(-0.77189267) q[1];
sx q[1];
rz(2.4235922) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7557186) q[0];
sx q[0];
rz(-1.0654952) q[0];
sx q[0];
rz(2.2875209) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23913582) q[2];
sx q[2];
rz(-2.1209361) q[2];
sx q[2];
rz(0.58708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79149063) q[1];
sx q[1];
rz(-1.4254505) q[1];
sx q[1];
rz(2.3816557) q[1];
x q[2];
rz(-0.15851373) q[3];
sx q[3];
rz(-1.4441381) q[3];
sx q[3];
rz(-1.4617006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2672853) q[2];
sx q[2];
rz(-0.75338805) q[2];
sx q[2];
rz(2.4370952) q[2];
rz(2.7590175) q[3];
sx q[3];
rz(-1.7035328) q[3];
sx q[3];
rz(-2.2883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48069561) q[0];
sx q[0];
rz(-0.040204164) q[0];
sx q[0];
rz(0.30886343) q[0];
rz(-0.94912306) q[1];
sx q[1];
rz(-0.32052761) q[1];
sx q[1];
rz(0.55108756) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8723485) q[0];
sx q[0];
rz(-1.9595032) q[0];
sx q[0];
rz(1.7354911) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91294841) q[2];
sx q[2];
rz(-1.5163053) q[2];
sx q[2];
rz(1.8860429) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6092564) q[1];
sx q[1];
rz(-0.21363959) q[1];
sx q[1];
rz(-2.3268301) q[1];
rz(-pi) q[2];
rz(-2.9948009) q[3];
sx q[3];
rz(-0.78189584) q[3];
sx q[3];
rz(1.1679389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71517313) q[2];
sx q[2];
rz(-0.74843633) q[2];
sx q[2];
rz(-0.61881649) q[2];
rz(-0.62301451) q[3];
sx q[3];
rz(-2.2996733) q[3];
sx q[3];
rz(2.714341) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1022559) q[0];
sx q[0];
rz(-2.4781041) q[0];
sx q[0];
rz(2.6797507) q[0];
rz(-2.6448008) q[1];
sx q[1];
rz(-1.3235612) q[1];
sx q[1];
rz(1.7740446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0557308) q[0];
sx q[0];
rz(-1.6501969) q[0];
sx q[0];
rz(-2.2532796) q[0];
x q[1];
rz(-1.8470079) q[2];
sx q[2];
rz(-0.51339692) q[2];
sx q[2];
rz(-0.16977507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.313844) q[1];
sx q[1];
rz(-1.5987458) q[1];
sx q[1];
rz(-1.8716783) q[1];
rz(-pi) q[2];
rz(2.9799906) q[3];
sx q[3];
rz(-2.186612) q[3];
sx q[3];
rz(-0.018425724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0355012) q[2];
sx q[2];
rz(-2.0095299) q[2];
sx q[2];
rz(-0.79037017) q[2];
rz(2.5370989) q[3];
sx q[3];
rz(-1.0249745) q[3];
sx q[3];
rz(-2.7214971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6919493) q[0];
sx q[0];
rz(-2.255891) q[0];
sx q[0];
rz(2.7272136) q[0];
rz(2.5212133) q[1];
sx q[1];
rz(-1.9199771) q[1];
sx q[1];
rz(-1.5588123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5904424) q[0];
sx q[0];
rz(-1.4865459) q[0];
sx q[0];
rz(0.1131375) q[0];
rz(-0.16838603) q[2];
sx q[2];
rz(-0.10659519) q[2];
sx q[2];
rz(1.835198) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.61395105) q[1];
sx q[1];
rz(-1.5802586) q[1];
sx q[1];
rz(-1.5130312) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2101733) q[3];
sx q[3];
rz(-2.7322331) q[3];
sx q[3];
rz(1.6334074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9862426) q[2];
sx q[2];
rz(-2.4537931) q[2];
sx q[2];
rz(2.9637994) q[2];
rz(-0.47884652) q[3];
sx q[3];
rz(-1.1136473) q[3];
sx q[3];
rz(-0.068304121) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9947522) q[0];
sx q[0];
rz(-0.9786334) q[0];
sx q[0];
rz(-3.0290208) q[0];
rz(-0.62537891) q[1];
sx q[1];
rz(-1.1275147) q[1];
sx q[1];
rz(-0.88923997) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5035985) q[0];
sx q[0];
rz(-1.6498008) q[0];
sx q[0];
rz(-0.043009506) q[0];
rz(2.8577789) q[2];
sx q[2];
rz(-0.2760016) q[2];
sx q[2];
rz(0.93940473) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9611284) q[1];
sx q[1];
rz(-0.58149177) q[1];
sx q[1];
rz(-1.1322458) q[1];
x q[2];
rz(-1.8440644) q[3];
sx q[3];
rz(-2.8972761) q[3];
sx q[3];
rz(-0.48297802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.858295) q[2];
sx q[2];
rz(-0.12792835) q[2];
sx q[2];
rz(-2.8683786) q[2];
rz(0.22935271) q[3];
sx q[3];
rz(-1.554824) q[3];
sx q[3];
rz(1.2685512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.6943618) q[0];
sx q[0];
rz(-2.2706967) q[0];
sx q[0];
rz(-0.91621512) q[0];
rz(-0.18237309) q[1];
sx q[1];
rz(-2.9353751) q[1];
sx q[1];
rz(1.1590385) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75814509) q[0];
sx q[0];
rz(-2.8383857) q[0];
sx q[0];
rz(-2.4571553) q[0];
x q[1];
rz(3.0460814) q[2];
sx q[2];
rz(-1.6217124) q[2];
sx q[2];
rz(-1.7820304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9578286) q[1];
sx q[1];
rz(-1.2570341) q[1];
sx q[1];
rz(-2.9741964) q[1];
x q[2];
rz(2.4156225) q[3];
sx q[3];
rz(-1.0683224) q[3];
sx q[3];
rz(-0.53086764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.75659043) q[2];
sx q[2];
rz(-0.80005163) q[2];
sx q[2];
rz(2.8263367) q[2];
rz(2.5473525) q[3];
sx q[3];
rz(-0.88163328) q[3];
sx q[3];
rz(-2.5074904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88157982) q[0];
sx q[0];
rz(-2.4729112) q[0];
sx q[0];
rz(2.9823533) q[0];
rz(-1.0881933) q[1];
sx q[1];
rz(-0.91251487) q[1];
sx q[1];
rz(-2.870627) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20734678) q[0];
sx q[0];
rz(-2.0306315) q[0];
sx q[0];
rz(2.8608972) q[0];
rz(-pi) q[1];
rz(2.4687211) q[2];
sx q[2];
rz(-1.4904163) q[2];
sx q[2];
rz(-0.56416914) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1098699) q[1];
sx q[1];
rz(-0.9389482) q[1];
sx q[1];
rz(-1.9919591) q[1];
rz(2.6600574) q[3];
sx q[3];
rz(-1.8706483) q[3];
sx q[3];
rz(-1.4389584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8012041) q[2];
sx q[2];
rz(-2.3694254) q[2];
sx q[2];
rz(-2.8188952) q[2];
rz(0.05803756) q[3];
sx q[3];
rz(-2.3400584) q[3];
sx q[3];
rz(0.29840741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.67644453) q[0];
sx q[0];
rz(-1.5100751) q[0];
sx q[0];
rz(-0.81612192) q[0];
rz(-0.25794087) q[1];
sx q[1];
rz(-1.1358658) q[1];
sx q[1];
rz(1.6075016) q[1];
rz(1.8564687) q[2];
sx q[2];
rz(-1.8485879) q[2];
sx q[2];
rz(-1.585203) q[2];
rz(-0.65597038) q[3];
sx q[3];
rz(-1.9269939) q[3];
sx q[3];
rz(-3.1300598) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
