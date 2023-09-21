OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5972714) q[0];
sx q[0];
rz(-0.40364021) q[0];
sx q[0];
rz(-0.37024745) q[0];
rz(0.20180841) q[1];
sx q[1];
rz(-1.2528074) q[1];
sx q[1];
rz(1.7226146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0050632) q[0];
sx q[0];
rz(-0.29924527) q[0];
sx q[0];
rz(-0.61009272) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6909157) q[2];
sx q[2];
rz(-2.5833231) q[2];
sx q[2];
rz(2.6435341) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4048684) q[1];
sx q[1];
rz(-2.1202592) q[1];
sx q[1];
rz(-3.0800372) q[1];
x q[2];
rz(-2.8098104) q[3];
sx q[3];
rz(-0.79091573) q[3];
sx q[3];
rz(-0.055981759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3966763) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(-2.7837226) q[2];
rz(-2.9499124) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.1382004) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(-0.068280846) q[0];
rz(-2.1349019) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(2.4904747) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40819528) q[0];
sx q[0];
rz(-1.7719643) q[0];
sx q[0];
rz(0.64543076) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6158193) q[2];
sx q[2];
rz(-1.5915046) q[2];
sx q[2];
rz(0.5171585) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.843833) q[1];
sx q[1];
rz(-2.4503772) q[1];
sx q[1];
rz(2.3046231) q[1];
x q[2];
rz(0.39204709) q[3];
sx q[3];
rz(-2.6169176) q[3];
sx q[3];
rz(2.3384561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-0.15658997) q[2];
sx q[2];
rz(0.90332705) q[2];
rz(-0.7545169) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2925401) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(1.6429398) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-1.0292056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.878359) q[0];
sx q[0];
rz(-1.9089111) q[0];
sx q[0];
rz(2.6792206) q[0];
rz(1.9919954) q[2];
sx q[2];
rz(-1.8987978) q[2];
sx q[2];
rz(1.0331819) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.610565) q[1];
sx q[1];
rz(-0.37279168) q[1];
sx q[1];
rz(-0.60786604) q[1];
rz(-pi) q[2];
rz(1.3658931) q[3];
sx q[3];
rz(-0.5538867) q[3];
sx q[3];
rz(2.3212471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2091973) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.8445245) q[2];
rz(-2.5629937) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.2066752) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(-2.7329965) q[0];
rz(-1.4273377) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(-2.2600007) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20760575) q[0];
sx q[0];
rz(-2.2946977) q[0];
sx q[0];
rz(-1.8007604) q[0];
x q[1];
rz(-2.2927631) q[2];
sx q[2];
rz(-2.1049044) q[2];
sx q[2];
rz(1.7497077) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2243005) q[1];
sx q[1];
rz(-0.91887337) q[1];
sx q[1];
rz(0.93445458) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6374171) q[3];
sx q[3];
rz(-1.9873709) q[3];
sx q[3];
rz(-1.4814324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0980229) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(1.1126474) q[2];
rz(0.55551314) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1495789) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(-1.7768815) q[0];
rz(0.31750202) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(-0.11725765) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1299767) q[0];
sx q[0];
rz(-0.59966171) q[0];
sx q[0];
rz(2.9090803) q[0];
x q[1];
rz(1.3567032) q[2];
sx q[2];
rz(-0.97852409) q[2];
sx q[2];
rz(-1.585373) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14355625) q[1];
sx q[1];
rz(-1.7245502) q[1];
sx q[1];
rz(-2.9424332) q[1];
x q[2];
rz(-0.076905964) q[3];
sx q[3];
rz(-2.1566026) q[3];
sx q[3];
rz(-2.2211071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34510288) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(-1.3797181) q[2];
rz(1.951925) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(-0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18774408) q[0];
sx q[0];
rz(-1.6402316) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(1.1114936) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(3.0336753) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2544125) q[0];
sx q[0];
rz(-0.88258703) q[0];
sx q[0];
rz(-1.956091) q[0];
rz(-3.1015322) q[2];
sx q[2];
rz(-1.5594348) q[2];
sx q[2];
rz(0.60919112) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2006827) q[1];
sx q[1];
rz(-2.4001277) q[1];
sx q[1];
rz(2.7156668) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2599254) q[3];
sx q[3];
rz(-1.1881184) q[3];
sx q[3];
rz(1.4985639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(-2.0325913) q[2];
rz(-1.2935982) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368822) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(0.014904508) q[0];
rz(2.7203454) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(2.3419535) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22916238) q[0];
sx q[0];
rz(-1.8254571) q[0];
sx q[0];
rz(-2.4486662) q[0];
rz(-1.2866576) q[2];
sx q[2];
rz(-1.8737027) q[2];
sx q[2];
rz(-2.6877833) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9866242) q[1];
sx q[1];
rz(-1.5151549) q[1];
sx q[1];
rz(2.9961622) q[1];
rz(0.13749595) q[3];
sx q[3];
rz(-2.1318448) q[3];
sx q[3];
rz(1.4631127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(-2.2195393) q[2];
rz(-1.3098035) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9668982) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(-0.2640557) q[0];
rz(-1.7407725) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(-0.31731269) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50002977) q[0];
sx q[0];
rz(-0.39425685) q[0];
sx q[0];
rz(1.7612329) q[0];
rz(-1.1515456) q[2];
sx q[2];
rz(-2.7630685) q[2];
sx q[2];
rz(-2.7839157) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6416157) q[1];
sx q[1];
rz(-2.9950905) q[1];
sx q[1];
rz(0.67540692) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64538892) q[3];
sx q[3];
rz(-1.7040952) q[3];
sx q[3];
rz(-0.79750878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7060966) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(-1.6395817) q[2];
rz(2.8912985) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(1.1423473) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2450927) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(-1.7171575) q[0];
rz(-2.6622488) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(-3.0260578) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45142052) q[0];
sx q[0];
rz(-0.85309404) q[0];
sx q[0];
rz(-2.3510128) q[0];
x q[1];
rz(-2.8250474) q[2];
sx q[2];
rz(-2.4015421) q[2];
sx q[2];
rz(-0.22063247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9563453) q[1];
sx q[1];
rz(-2.6223409) q[1];
sx q[1];
rz(-0.96038702) q[1];
rz(-0.69620903) q[3];
sx q[3];
rz(-1.6632102) q[3];
sx q[3];
rz(-0.19314167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66403786) q[2];
sx q[2];
rz(-2.0319735) q[2];
sx q[2];
rz(1.4257365) q[2];
rz(1.7410295) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(-2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(0.26915959) q[0];
rz(-2.0458938) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(-1.7620618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1901967) q[0];
sx q[0];
rz(-1.7362036) q[0];
sx q[0];
rz(3.0347546) q[0];
x q[1];
rz(1.7494406) q[2];
sx q[2];
rz(-0.42334712) q[2];
sx q[2];
rz(1.64738) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60296842) q[1];
sx q[1];
rz(-1.5403962) q[1];
sx q[1];
rz(-1.7441185) q[1];
rz(-0.74897154) q[3];
sx q[3];
rz(-0.43991551) q[3];
sx q[3];
rz(-2.9797152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.049008869) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(-1.5157549) q[2];
rz(-1.9745291) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-2.1102171) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2172858) q[0];
sx q[0];
rz(-1.8368245) q[0];
sx q[0];
rz(0.40689847) q[0];
rz(-2.6869607) q[1];
sx q[1];
rz(-2.0352719) q[1];
sx q[1];
rz(-0.24771053) q[1];
rz(2.968593) q[2];
sx q[2];
rz(-2.0377918) q[2];
sx q[2];
rz(1.3755058) q[2];
rz(1.3116253) q[3];
sx q[3];
rz(-1.5117241) q[3];
sx q[3];
rz(-1.1321887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
