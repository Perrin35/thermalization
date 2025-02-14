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
rz(1.4662161) q[0];
sx q[0];
rz(5.3386547) q[0];
sx q[0];
rz(9.6261779) q[0];
rz(4.2138777) q[1];
sx q[1];
rz(3.9044851) q[1];
sx q[1];
rz(4.8621096) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22682193) q[0];
sx q[0];
rz(-2.1807007) q[0];
sx q[0];
rz(0.066909153) q[0];
x q[1];
rz(0.51789639) q[2];
sx q[2];
rz(-0.96518789) q[2];
sx q[2];
rz(0.54010682) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.101946) q[1];
sx q[1];
rz(-1.1803487) q[1];
sx q[1];
rz(2.4389308) q[1];
rz(-pi) q[2];
rz(-0.33290036) q[3];
sx q[3];
rz(-1.1984636) q[3];
sx q[3];
rz(-1.2941192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1066771) q[2];
sx q[2];
rz(-2.3372529) q[2];
sx q[2];
rz(-2.8625281) q[2];
rz(1.2532824) q[3];
sx q[3];
rz(-2.8918355) q[3];
sx q[3];
rz(2.9899924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5308373) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(-0.24638677) q[0];
rz(-1.1950182) q[1];
sx q[1];
rz(-0.89391005) q[1];
sx q[1];
rz(-1.6349207) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6794066) q[0];
sx q[0];
rz(-2.6194318) q[0];
sx q[0];
rz(2.0798111) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8737239) q[2];
sx q[2];
rz(-1.5298577) q[2];
sx q[2];
rz(2.8589279) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7333201) q[1];
sx q[1];
rz(-2.773914) q[1];
sx q[1];
rz(1.7321083) q[1];
rz(-pi) q[2];
rz(2.0923268) q[3];
sx q[3];
rz(-1.9535011) q[3];
sx q[3];
rz(-2.692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7404777) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(-1.231989) q[2];
rz(1.1022107) q[3];
sx q[3];
rz(-0.004318459) q[3];
sx q[3];
rz(1.8050885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-0.44322893) q[0];
sx q[0];
rz(-0.6969499) q[0];
sx q[0];
rz(-2.2337636) q[0];
rz(-2.5746131) q[1];
sx q[1];
rz(-2.1588529) q[1];
sx q[1];
rz(3.0388015) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5823321) q[0];
sx q[0];
rz(-2.1301422) q[0];
sx q[0];
rz(-1.3224959) q[0];
x q[1];
rz(-3.0605761) q[2];
sx q[2];
rz(-1.6747432) q[2];
sx q[2];
rz(0.31652094) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5788865) q[1];
sx q[1];
rz(-1.9745312) q[1];
sx q[1];
rz(2.5078234) q[1];
rz(-2.8277367) q[3];
sx q[3];
rz(-2.102763) q[3];
sx q[3];
rz(0.44600081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8583782) q[2];
sx q[2];
rz(-2.3329222) q[2];
sx q[2];
rz(1.4374479) q[2];
rz(0.39250675) q[3];
sx q[3];
rz(-1.1350574) q[3];
sx q[3];
rz(-0.2984305) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0482386) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(1.9527973) q[0];
rz(-0.28611723) q[1];
sx q[1];
rz(-2.5240099) q[1];
sx q[1];
rz(-1.5123222) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85507369) q[0];
sx q[0];
rz(-0.60337702) q[0];
sx q[0];
rz(2.0749178) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3740923) q[2];
sx q[2];
rz(-1.3250809) q[2];
sx q[2];
rz(1.9227288) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.69763598) q[1];
sx q[1];
rz(-0.62503615) q[1];
sx q[1];
rz(1.2500863) q[1];
x q[2];
rz(-1.9387454) q[3];
sx q[3];
rz(-0.5371437) q[3];
sx q[3];
rz(0.14432553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35147038) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(2.4646087) q[2];
rz(-1.8949159) q[3];
sx q[3];
rz(-1.8691984) q[3];
sx q[3];
rz(1.5164794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93382728) q[0];
sx q[0];
rz(-2.8849869) q[0];
sx q[0];
rz(0.024913464) q[0];
rz(1.0554396) q[1];
sx q[1];
rz(-2.1565304) q[1];
sx q[1];
rz(-2.8581462) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6993857) q[0];
sx q[0];
rz(-1.8309347) q[0];
sx q[0];
rz(-0.26630638) q[0];
rz(-pi) q[1];
x q[1];
rz(3.010036) q[2];
sx q[2];
rz(-1.3765928) q[2];
sx q[2];
rz(2.921791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3230309) q[1];
sx q[1];
rz(-0.50172443) q[1];
sx q[1];
rz(1.5072359) q[1];
rz(0.89783606) q[3];
sx q[3];
rz(-1.8404318) q[3];
sx q[3];
rz(1.8561038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31263605) q[2];
sx q[2];
rz(-2.7241311) q[2];
sx q[2];
rz(-2.8157595) q[2];
rz(1.2827986) q[3];
sx q[3];
rz(-1.9841586) q[3];
sx q[3];
rz(1.5939943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40389898) q[0];
sx q[0];
rz(-1.4827381) q[0];
sx q[0];
rz(-2.7666336) q[0];
rz(-0.47863475) q[1];
sx q[1];
rz(-0.67239434) q[1];
sx q[1];
rz(0.37014827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3309114) q[0];
sx q[0];
rz(-3.1157012) q[0];
sx q[0];
rz(2.6736892) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13932087) q[2];
sx q[2];
rz(-1.3874386) q[2];
sx q[2];
rz(-0.56688165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.0035841442) q[1];
sx q[1];
rz(-0.89911997) q[1];
sx q[1];
rz(0.9673311) q[1];
rz(-0.42780723) q[3];
sx q[3];
rz(-1.5662346) q[3];
sx q[3];
rz(3.004194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39773539) q[2];
sx q[2];
rz(-0.19146679) q[2];
sx q[2];
rz(1.7605304) q[2];
rz(-2.0928275) q[3];
sx q[3];
rz(-1.3449113) q[3];
sx q[3];
rz(-0.63961187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8884856) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(-0.92415586) q[0];
rz(2.2249075) q[1];
sx q[1];
rz(-1.0662096) q[1];
sx q[1];
rz(-1.7402657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9144672) q[0];
sx q[0];
rz(-2.0642529) q[0];
sx q[0];
rz(1.5052133) q[0];
rz(-1.7426874) q[2];
sx q[2];
rz(-0.92280932) q[2];
sx q[2];
rz(-1.6164219) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1223334) q[1];
sx q[1];
rz(-2.386555) q[1];
sx q[1];
rz(0.23438677) q[1];
x q[2];
rz(1.0823233) q[3];
sx q[3];
rz(-1.7632177) q[3];
sx q[3];
rz(0.50979641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2283198) q[2];
sx q[2];
rz(-1.6797804) q[2];
sx q[2];
rz(1.9165967) q[2];
rz(-0.0082958881) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(0.86709658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5834354) q[0];
sx q[0];
rz(-2.3112516) q[0];
sx q[0];
rz(0.48026568) q[0];
rz(2.9971314) q[1];
sx q[1];
rz(-0.90765777) q[1];
sx q[1];
rz(-0.86437782) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68587559) q[0];
sx q[0];
rz(-0.17016958) q[0];
sx q[0];
rz(1.8491114) q[0];
rz(0.34274613) q[2];
sx q[2];
rz(-3.0411093) q[2];
sx q[2];
rz(0.1048987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0740114) q[1];
sx q[1];
rz(-1.0339104) q[1];
sx q[1];
rz(2.8606961) q[1];
x q[2];
rz(-0.79473991) q[3];
sx q[3];
rz(-0.38215548) q[3];
sx q[3];
rz(2.4433745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9353443) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(0.63507357) q[2];
rz(-2.5687929) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(0.87535453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0315345) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(-3.0173259) q[0];
rz(0.36525137) q[1];
sx q[1];
rz(-1.7144014) q[1];
sx q[1];
rz(1.431538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54473684) q[0];
sx q[0];
rz(-1.3787621) q[0];
sx q[0];
rz(1.7512683) q[0];
rz(-1.0042436) q[2];
sx q[2];
rz(-0.99663094) q[2];
sx q[2];
rz(3.1176709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6874976) q[1];
sx q[1];
rz(-0.75088596) q[1];
sx q[1];
rz(-2.6862957) q[1];
rz(-pi) q[2];
rz(0.67691524) q[3];
sx q[3];
rz(-0.64161086) q[3];
sx q[3];
rz(0.42719242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8753836) q[2];
sx q[2];
rz(-2.2201846) q[2];
sx q[2];
rz(-1.9688152) q[2];
rz(1.4069936) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(1.9269358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59239546) q[0];
sx q[0];
rz(-2.0068491) q[0];
sx q[0];
rz(-2.5469575) q[0];
rz(-2.1356964) q[1];
sx q[1];
rz(-1.2024095) q[1];
sx q[1];
rz(-2.2524021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2834463) q[0];
sx q[0];
rz(-2.915417) q[0];
sx q[0];
rz(1.7986138) q[0];
rz(-pi) q[1];
rz(-3.0709549) q[2];
sx q[2];
rz(-1.1114235) q[2];
sx q[2];
rz(2.2769711) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8831625) q[1];
sx q[1];
rz(-3.0297802) q[1];
sx q[1];
rz(-2.1360141) q[1];
rz(-pi) q[2];
x q[2];
rz(0.066915705) q[3];
sx q[3];
rz(-2.624238) q[3];
sx q[3];
rz(-0.87393119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95120007) q[2];
sx q[2];
rz(-1.9659646) q[2];
sx q[2];
rz(0.54538837) q[2];
rz(0.96327463) q[3];
sx q[3];
rz(-1.9656209) q[3];
sx q[3];
rz(1.4639328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49190285) q[0];
sx q[0];
rz(-2.2408673) q[0];
sx q[0];
rz(1.3229205) q[0];
rz(-0.84359618) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(2.0100935) q[2];
sx q[2];
rz(-1.7514624) q[2];
sx q[2];
rz(-3.024586) q[2];
rz(-0.58335494) q[3];
sx q[3];
rz(-1.1370549) q[3];
sx q[3];
rz(2.1255253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
