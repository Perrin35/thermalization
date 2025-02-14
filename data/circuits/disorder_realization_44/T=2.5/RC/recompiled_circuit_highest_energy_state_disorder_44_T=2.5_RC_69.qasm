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
rz(-0.62701464) q[0];
sx q[0];
rz(3.7778683) q[0];
sx q[0];
rz(10.502622) q[0];
rz(1.0765422) q[1];
sx q[1];
rz(-0.62907469) q[1];
sx q[1];
rz(-1.0318626) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916578) q[0];
sx q[0];
rz(-0.0041793267) q[0];
sx q[0];
rz(3.0463534) q[0];
x q[1];
rz(-2.7912895) q[2];
sx q[2];
rz(-1.7997348) q[2];
sx q[2];
rz(-0.13267429) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2950287) q[1];
sx q[1];
rz(-1.16638) q[1];
sx q[1];
rz(1.9454898) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4937819) q[3];
sx q[3];
rz(-1.9056013) q[3];
sx q[3];
rz(0.96559262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2304113) q[2];
sx q[2];
rz(-1.958467) q[2];
sx q[2];
rz(0.46768701) q[2];
rz(-0.45105252) q[3];
sx q[3];
rz(-0.39877287) q[3];
sx q[3];
rz(0.27802813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4462047) q[0];
sx q[0];
rz(-0.22286335) q[0];
sx q[0];
rz(0.15431246) q[0];
rz(0.34126869) q[1];
sx q[1];
rz(-2.7879265) q[1];
sx q[1];
rz(-2.7214859) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1695009) q[0];
sx q[0];
rz(-1.7636714) q[0];
sx q[0];
rz(2.4986096) q[0];
x q[1];
rz(-2.2809199) q[2];
sx q[2];
rz(-1.6063074) q[2];
sx q[2];
rz(1.6954317) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1359033) q[1];
sx q[1];
rz(-1.5694251) q[1];
sx q[1];
rz(2.6708009) q[1];
rz(-pi) q[2];
rz(1.4147725) q[3];
sx q[3];
rz(-2.3120263) q[3];
sx q[3];
rz(-1.6680731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61098617) q[2];
sx q[2];
rz(-0.89908081) q[2];
sx q[2];
rz(2.3685624) q[2];
rz(-2.2981339) q[3];
sx q[3];
rz(-1.5638899) q[3];
sx q[3];
rz(0.57190603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.968349) q[0];
sx q[0];
rz(-2.1143715) q[0];
sx q[0];
rz(-2.9395043) q[0];
rz(-2.659722) q[1];
sx q[1];
rz(-2.9402969) q[1];
sx q[1];
rz(0.906382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68708352) q[0];
sx q[0];
rz(-2.2453899) q[0];
sx q[0];
rz(-0.97954155) q[0];
x q[1];
rz(-0.34445539) q[2];
sx q[2];
rz(-2.0018501) q[2];
sx q[2];
rz(1.5004917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6828595) q[1];
sx q[1];
rz(-0.56536973) q[1];
sx q[1];
rz(2.8394147) q[1];
rz(1.866147) q[3];
sx q[3];
rz(-1.4100299) q[3];
sx q[3];
rz(-2.8637342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0704982) q[2];
sx q[2];
rz(-1.1134032) q[2];
sx q[2];
rz(-2.5615198) q[2];
rz(1.5813367) q[3];
sx q[3];
rz(-1.7507078) q[3];
sx q[3];
rz(-1.4238547) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4105014) q[0];
sx q[0];
rz(-0.061256496) q[0];
sx q[0];
rz(-3.0938003) q[0];
rz(-1.2740678) q[1];
sx q[1];
rz(-1.2893226) q[1];
sx q[1];
rz(0.045225708) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84789373) q[0];
sx q[0];
rz(-1.7195104) q[0];
sx q[0];
rz(-0.92647628) q[0];
x q[1];
rz(2.2061646) q[2];
sx q[2];
rz(-1.8385988) q[2];
sx q[2];
rz(1.8647461) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2312551) q[1];
sx q[1];
rz(-1.5738166) q[1];
sx q[1];
rz(-1.5763723) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2269656) q[3];
sx q[3];
rz(-2.4883929) q[3];
sx q[3];
rz(0.40756215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4101326) q[2];
sx q[2];
rz(-1.8803909) q[2];
sx q[2];
rz(0.88978466) q[2];
rz(2.5951923) q[3];
sx q[3];
rz(-1.0011287) q[3];
sx q[3];
rz(0.31118292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50362098) q[0];
sx q[0];
rz(-1.6872971) q[0];
sx q[0];
rz(-1.7244435) q[0];
rz(-2.0218938) q[1];
sx q[1];
rz(-1.3816625) q[1];
sx q[1];
rz(1.959257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1954297) q[0];
sx q[0];
rz(-0.39349213) q[0];
sx q[0];
rz(-0.18184967) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77967092) q[2];
sx q[2];
rz(-2.5354249) q[2];
sx q[2];
rz(0.066628284) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4880138) q[1];
sx q[1];
rz(-1.5216478) q[1];
sx q[1];
rz(2.4202511) q[1];
rz(-2.6695146) q[3];
sx q[3];
rz(-2.3340477) q[3];
sx q[3];
rz(-2.9387752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.917439) q[2];
sx q[2];
rz(-0.76442337) q[2];
sx q[2];
rz(2.7177641) q[2];
rz(2.0297) q[3];
sx q[3];
rz(-0.49911505) q[3];
sx q[3];
rz(2.4901938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3424585) q[0];
sx q[0];
rz(-0.0722216) q[0];
sx q[0];
rz(2.1891731) q[0];
rz(-2.3656288) q[1];
sx q[1];
rz(-2.2064078) q[1];
sx q[1];
rz(-1.5663358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99235202) q[0];
sx q[0];
rz(-1.5910205) q[0];
sx q[0];
rz(-3.1035191) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96409728) q[2];
sx q[2];
rz(-1.6025873) q[2];
sx q[2];
rz(-0.83039415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.414622) q[1];
sx q[1];
rz(-1.3515673) q[1];
sx q[1];
rz(2.5269233) q[1];
rz(-2.9569472) q[3];
sx q[3];
rz(-1.0870516) q[3];
sx q[3];
rz(0.32102206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39885193) q[2];
sx q[2];
rz(-0.80790085) q[2];
sx q[2];
rz(0.74637949) q[2];
rz(1.0910723) q[3];
sx q[3];
rz(-1.4336136) q[3];
sx q[3];
rz(0.55238849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93254507) q[0];
sx q[0];
rz(-3.0952251) q[0];
sx q[0];
rz(-0.5433425) q[0];
rz(2.6047193) q[1];
sx q[1];
rz(-2.8165292) q[1];
sx q[1];
rz(-3.0184025) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7167599) q[0];
sx q[0];
rz(-0.93532978) q[0];
sx q[0];
rz(3.1034971) q[0];
rz(-2.6307879) q[2];
sx q[2];
rz(-1.1996562) q[2];
sx q[2];
rz(0.57194158) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.035288485) q[1];
sx q[1];
rz(-1.0266487) q[1];
sx q[1];
rz(2.418251) q[1];
x q[2];
rz(2.3770272) q[3];
sx q[3];
rz(-0.95029059) q[3];
sx q[3];
rz(2.4303049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3049551) q[2];
sx q[2];
rz(-1.387549) q[2];
sx q[2];
rz(0.48509625) q[2];
rz(-1.1525611) q[3];
sx q[3];
rz(-1.7771114) q[3];
sx q[3];
rz(3.0936354) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3845859) q[0];
sx q[0];
rz(-2.204019) q[0];
sx q[0];
rz(0.49569976) q[0];
rz(-3.0777625) q[1];
sx q[1];
rz(-0.7494691) q[1];
sx q[1];
rz(3.0705423) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249504) q[0];
sx q[0];
rz(-2.0293689) q[0];
sx q[0];
rz(1.8383265) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7456648) q[2];
sx q[2];
rz(-2.1687458) q[2];
sx q[2];
rz(1.1451461) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9378375) q[1];
sx q[1];
rz(-1.3295638) q[1];
sx q[1];
rz(-2.9822442) q[1];
rz(-pi) q[2];
rz(-2.9661353) q[3];
sx q[3];
rz(-0.4532632) q[3];
sx q[3];
rz(-2.3153265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.12585982) q[2];
sx q[2];
rz(-1.6682397) q[2];
sx q[2];
rz(-2.2802172) q[2];
rz(1.8371948) q[3];
sx q[3];
rz(-2.5671037) q[3];
sx q[3];
rz(1.5484352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17643377) q[0];
sx q[0];
rz(-2.6506944) q[0];
sx q[0];
rz(1.4392256) q[0];
rz(3.0610415) q[1];
sx q[1];
rz(-1.4713902) q[1];
sx q[1];
rz(1.7505987) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0223402) q[0];
sx q[0];
rz(-1.0588405) q[0];
sx q[0];
rz(-0.48959022) q[0];
rz(0.020419196) q[2];
sx q[2];
rz(-1.8755185) q[2];
sx q[2];
rz(0.29856506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8409001) q[1];
sx q[1];
rz(-2.024902) q[1];
sx q[1];
rz(0.91030376) q[1];
rz(-pi) q[2];
rz(-0.37064044) q[3];
sx q[3];
rz(-1.9270867) q[3];
sx q[3];
rz(2.4161937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3535658) q[2];
sx q[2];
rz(-1.2949508) q[2];
sx q[2];
rz(3.0143152) q[2];
rz(-2.2150529) q[3];
sx q[3];
rz(-2.8997731) q[3];
sx q[3];
rz(1.658879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.7552898) q[0];
sx q[0];
rz(-2.4918064) q[0];
sx q[0];
rz(1.5007098) q[0];
rz(0.81037784) q[1];
sx q[1];
rz(-0.85226285) q[1];
sx q[1];
rz(-0.77218974) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92896748) q[0];
sx q[0];
rz(-1.0662931) q[0];
sx q[0];
rz(0.8302159) q[0];
rz(-pi) q[1];
rz(-0.29532642) q[2];
sx q[2];
rz(-1.3662896) q[2];
sx q[2];
rz(2.6596065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3153684) q[1];
sx q[1];
rz(-1.5824707) q[1];
sx q[1];
rz(1.0213901) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2507319) q[3];
sx q[3];
rz(-1.5308342) q[3];
sx q[3];
rz(0.077629493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3810252) q[2];
sx q[2];
rz(-2.3367391) q[2];
sx q[2];
rz(0.62925657) q[2];
rz(-0.73090807) q[3];
sx q[3];
rz(-0.84354246) q[3];
sx q[3];
rz(-1.4430911) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44169852) q[0];
sx q[0];
rz(-0.48342539) q[0];
sx q[0];
rz(-0.40406686) q[0];
rz(-0.40456698) q[1];
sx q[1];
rz(-1.7351983) q[1];
sx q[1];
rz(2.4529967) q[1];
rz(0.2307288) q[2];
sx q[2];
rz(-2.8525272) q[2];
sx q[2];
rz(2.0244103) q[2];
rz(2.6647207) q[3];
sx q[3];
rz(-0.58860368) q[3];
sx q[3];
rz(0.5640201) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
