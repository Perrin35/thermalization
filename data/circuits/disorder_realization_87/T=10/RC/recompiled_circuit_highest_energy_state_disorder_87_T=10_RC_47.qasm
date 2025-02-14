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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0312289) q[0];
sx q[0];
rz(-2.5284934) q[0];
sx q[0];
rz(-1.6661864) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95008738) q[2];
sx q[2];
rz(-2.3664775) q[2];
sx q[2];
rz(2.895854) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78064976) q[1];
sx q[1];
rz(-2.2113178) q[1];
sx q[1];
rz(-2.065413) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33290036) q[3];
sx q[3];
rz(-1.9431291) q[3];
sx q[3];
rz(1.8474735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1066771) q[2];
sx q[2];
rz(-0.80433977) q[2];
sx q[2];
rz(0.2790645) q[2];
rz(-1.8883102) q[3];
sx q[3];
rz(-0.24975714) q[3];
sx q[3];
rz(-2.9899924) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6107553) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(0.24638677) q[0];
rz(-1.1950182) q[1];
sx q[1];
rz(-2.2476826) q[1];
sx q[1];
rz(1.6349207) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55913299) q[0];
sx q[0];
rz(-1.8163067) q[0];
sx q[0];
rz(-1.1051635) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2678688) q[2];
sx q[2];
rz(-1.5298577) q[2];
sx q[2];
rz(2.8589279) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8283702) q[1];
sx q[1];
rz(-1.6285609) q[1];
sx q[1];
rz(1.2074797) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0923268) q[3];
sx q[3];
rz(-1.9535011) q[3];
sx q[3];
rz(-2.692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40111497) q[2];
sx q[2];
rz(-2.1666708) q[2];
sx q[2];
rz(-1.231989) q[2];
rz(1.1022107) q[3];
sx q[3];
rz(-3.1372742) q[3];
sx q[3];
rz(1.3365041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6983637) q[0];
sx q[0];
rz(-2.4446428) q[0];
sx q[0];
rz(-0.90782905) q[0];
rz(0.56697956) q[1];
sx q[1];
rz(-0.98273977) q[1];
sx q[1];
rz(-3.0388015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9963275) q[0];
sx q[0];
rz(-1.3609556) q[0];
sx q[0];
rz(0.57344936) q[0];
rz(-0.91102131) q[2];
sx q[2];
rz(-0.13170031) q[2];
sx q[2];
rz(0.98051276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6425848) q[1];
sx q[1];
rz(-0.73621589) q[1];
sx q[1];
rz(0.6249439) q[1];
x q[2];
rz(-1.0877043) q[3];
sx q[3];
rz(-2.5317041) q[3];
sx q[3];
rz(2.1263378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8583782) q[2];
sx q[2];
rz(-0.80867043) q[2];
sx q[2];
rz(1.7041448) q[2];
rz(-0.39250675) q[3];
sx q[3];
rz(-2.0065353) q[3];
sx q[3];
rz(-0.2984305) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0482386) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(-1.1887953) q[0];
rz(-0.28611723) q[1];
sx q[1];
rz(-2.5240099) q[1];
sx q[1];
rz(-1.5123222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8522569) q[0];
sx q[0];
rz(-1.8484383) q[0];
sx q[0];
rz(-2.1136978) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3740923) q[2];
sx q[2];
rz(-1.3250809) q[2];
sx q[2];
rz(-1.9227288) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5315631) q[1];
sx q[1];
rz(-1.3852784) q[1];
sx q[1];
rz(2.1712135) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9387454) q[3];
sx q[3];
rz(-0.5371437) q[3];
sx q[3];
rz(2.9972671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7901223) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(-0.67698395) q[2];
rz(1.8949159) q[3];
sx q[3];
rz(-1.8691984) q[3];
sx q[3];
rz(-1.5164794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93382728) q[0];
sx q[0];
rz(-2.8849869) q[0];
sx q[0];
rz(3.1166792) q[0];
rz(-2.086153) q[1];
sx q[1];
rz(-2.1565304) q[1];
sx q[1];
rz(-2.8581462) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058539778) q[0];
sx q[0];
rz(-1.3136567) q[0];
sx q[0];
rz(1.8399946) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1589958) q[2];
sx q[2];
rz(-0.23410205) q[2];
sx q[2];
rz(-0.38123075) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.746096) q[1];
sx q[1];
rz(-2.0714134) q[1];
sx q[1];
rz(-3.1067645) q[1];
x q[2];
rz(-0.89783606) q[3];
sx q[3];
rz(-1.3011609) q[3];
sx q[3];
rz(-1.2854888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8289566) q[2];
sx q[2];
rz(-2.7241311) q[2];
sx q[2];
rz(-2.8157595) q[2];
rz(-1.8587941) q[3];
sx q[3];
rz(-1.9841586) q[3];
sx q[3];
rz(1.5939943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7376937) q[0];
sx q[0];
rz(-1.4827381) q[0];
sx q[0];
rz(-0.37495908) q[0];
rz(-0.47863475) q[1];
sx q[1];
rz(-2.4691983) q[1];
sx q[1];
rz(-0.37014827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86287303) q[0];
sx q[0];
rz(-1.5476883) q[0];
sx q[0];
rz(-1.5824759) q[0];
rz(-0.92801969) q[2];
sx q[2];
rz(-0.22980873) q[2];
sx q[2];
rz(1.2224876) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83997516) q[1];
sx q[1];
rz(-2.2711922) q[1];
sx q[1];
rz(0.61995929) q[1];
x q[2];
rz(-0.010995098) q[3];
sx q[3];
rz(-2.7137626) q[3];
sx q[3];
rz(1.7181991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7438573) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.3810623) q[2];
rz(-1.0487652) q[3];
sx q[3];
rz(-1.3449113) q[3];
sx q[3];
rz(-2.5019808) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25310707) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(-0.92415586) q[0];
rz(0.91668516) q[1];
sx q[1];
rz(-1.0662096) q[1];
sx q[1];
rz(-1.4013269) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2271254) q[0];
sx q[0];
rz(-2.0642529) q[0];
sx q[0];
rz(1.6363793) q[0];
rz(-pi) q[1];
rz(-0.22221128) q[2];
sx q[2];
rz(-2.4743818) q[2];
sx q[2];
rz(-1.3363584) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8055011) q[1];
sx q[1];
rz(-0.84118836) q[1];
sx q[1];
rz(-1.7859687) q[1];
rz(2.0592693) q[3];
sx q[3];
rz(-1.3783749) q[3];
sx q[3];
rz(-2.6317962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2283198) q[2];
sx q[2];
rz(-1.6797804) q[2];
sx q[2];
rz(1.224996) q[2];
rz(3.1332968) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(0.86709658) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55815721) q[0];
sx q[0];
rz(-2.3112516) q[0];
sx q[0];
rz(2.661327) q[0];
rz(0.14446124) q[1];
sx q[1];
rz(-0.90765777) q[1];
sx q[1];
rz(-2.2772148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7378993) q[0];
sx q[0];
rz(-1.4072352) q[0];
sx q[0];
rz(3.0944194) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7988465) q[2];
sx q[2];
rz(-0.10048332) q[2];
sx q[2];
rz(0.1048987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64329735) q[1];
sx q[1];
rz(-1.8113231) q[1];
sx q[1];
rz(-2.125419) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79473991) q[3];
sx q[3];
rz(-2.7594372) q[3];
sx q[3];
rz(0.69821815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9353443) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(-0.63507357) q[2];
rz(-0.57279974) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(2.2662381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0315345) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(3.0173259) q[0];
rz(0.36525137) q[1];
sx q[1];
rz(-1.7144014) q[1];
sx q[1];
rz(1.431538) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5968558) q[0];
sx q[0];
rz(-1.7628306) q[0];
sx q[0];
rz(1.7512683) q[0];
rz(1.0042436) q[2];
sx q[2];
rz(-2.1449617) q[2];
sx q[2];
rz(3.1176709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6811583) q[1];
sx q[1];
rz(-1.8755113) q[1];
sx q[1];
rz(2.4439993) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67691524) q[3];
sx q[3];
rz(-0.64161086) q[3];
sx q[3];
rz(-2.7144002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8753836) q[2];
sx q[2];
rz(-0.92140809) q[2];
sx q[2];
rz(1.1727775) q[2];
rz(1.4069936) q[3];
sx q[3];
rz(-1.3318136) q[3];
sx q[3];
rz(-1.9269358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239546) q[0];
sx q[0];
rz(-2.0068491) q[0];
sx q[0];
rz(2.5469575) q[0];
rz(1.0058962) q[1];
sx q[1];
rz(-1.2024095) q[1];
sx q[1];
rz(0.88919052) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93485661) q[0];
sx q[0];
rz(-1.5201269) q[0];
sx q[0];
rz(-1.3502747) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4290733) q[2];
sx q[2];
rz(-0.46438875) q[2];
sx q[2];
rz(-1.022867) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8915705) q[1];
sx q[1];
rz(-1.6305939) q[1];
sx q[1];
rz(-1.476261) q[1];
rz(-pi) q[2];
rz(3.0746769) q[3];
sx q[3];
rz(-2.624238) q[3];
sx q[3];
rz(-2.2676615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95120007) q[2];
sx q[2];
rz(-1.9659646) q[2];
sx q[2];
rz(0.54538837) q[2];
rz(2.178318) q[3];
sx q[3];
rz(-1.9656209) q[3];
sx q[3];
rz(1.6776599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6496898) q[0];
sx q[0];
rz(-0.90072537) q[0];
sx q[0];
rz(-1.8186722) q[0];
rz(2.2979965) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(-1.1314992) q[2];
sx q[2];
rz(-1.7514624) q[2];
sx q[2];
rz(-3.024586) q[2];
rz(-0.69915184) q[3];
sx q[3];
rz(-0.71157645) q[3];
sx q[3];
rz(-0.012307766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
