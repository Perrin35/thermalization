OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2332239) q[0];
sx q[0];
rz(-1.5492735) q[0];
sx q[0];
rz(-0.17750658) q[0];
rz(1.5597875) q[1];
sx q[1];
rz(-1.8145365) q[1];
sx q[1];
rz(-1.1787193) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8713407) q[0];
sx q[0];
rz(-0.9862904) q[0];
sx q[0];
rz(2.9181705) q[0];
rz(1.1792405) q[2];
sx q[2];
rz(-2.4610977) q[2];
sx q[2];
rz(3.1093564) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4376154) q[1];
sx q[1];
rz(-1.1363789) q[1];
sx q[1];
rz(2.2056286) q[1];
x q[2];
rz(1.4229544) q[3];
sx q[3];
rz(-0.73972046) q[3];
sx q[3];
rz(1.5802368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.061964758) q[2];
sx q[2];
rz(-2.2140333) q[2];
sx q[2];
rz(2.8947042) q[2];
rz(0.43761349) q[3];
sx q[3];
rz(-1.0381617) q[3];
sx q[3];
rz(-2.8810697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0831809) q[0];
sx q[0];
rz(-2.658168) q[0];
sx q[0];
rz(-1.1046945) q[0];
rz(-0.80766922) q[1];
sx q[1];
rz(-2.3876815) q[1];
sx q[1];
rz(2.5025936) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7154626) q[0];
sx q[0];
rz(-1.3952266) q[0];
sx q[0];
rz(0.28101729) q[0];
rz(0.27490567) q[2];
sx q[2];
rz(-1.6187117) q[2];
sx q[2];
rz(0.78457441) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5377755) q[1];
sx q[1];
rz(-0.81934419) q[1];
sx q[1];
rz(0.50756331) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86231972) q[3];
sx q[3];
rz(-0.42565027) q[3];
sx q[3];
rz(-1.5141534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69940007) q[2];
sx q[2];
rz(-1.7510119) q[2];
sx q[2];
rz(-0.81713444) q[2];
rz(2.6160431) q[3];
sx q[3];
rz(-0.81997973) q[3];
sx q[3];
rz(-1.9513963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84080559) q[0];
sx q[0];
rz(-1.5910933) q[0];
sx q[0];
rz(-0.15980414) q[0];
rz(0.5197168) q[1];
sx q[1];
rz(-2.0964041) q[1];
sx q[1];
rz(-0.91667485) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778939) q[0];
sx q[0];
rz(-0.5710154) q[0];
sx q[0];
rz(-3.0932941) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74059981) q[2];
sx q[2];
rz(-1.3889905) q[2];
sx q[2];
rz(-1.2642853) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3252249) q[1];
sx q[1];
rz(-1.7599517) q[1];
sx q[1];
rz(0.048859096) q[1];
rz(2.5289747) q[3];
sx q[3];
rz(-1.4605986) q[3];
sx q[3];
rz(1.2904905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2436169) q[2];
sx q[2];
rz(-2.6349082) q[2];
sx q[2];
rz(-1.3755196) q[2];
rz(3.1084133) q[3];
sx q[3];
rz(-1.5539919) q[3];
sx q[3];
rz(-2.3110265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0631977) q[0];
sx q[0];
rz(-0.66135186) q[0];
sx q[0];
rz(-0.98264328) q[0];
rz(-0.64535087) q[1];
sx q[1];
rz(-1.9159562) q[1];
sx q[1];
rz(-0.4172079) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66574341) q[0];
sx q[0];
rz(-1.8717093) q[0];
sx q[0];
rz(-2.6074661) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3991314) q[2];
sx q[2];
rz(-2.7560803) q[2];
sx q[2];
rz(-0.76010349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1712492) q[1];
sx q[1];
rz(-1.6303195) q[1];
sx q[1];
rz(-0.67386676) q[1];
rz(0.44261114) q[3];
sx q[3];
rz(-2.5472982) q[3];
sx q[3];
rz(-1.8123008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8365606) q[2];
sx q[2];
rz(-0.3362259) q[2];
sx q[2];
rz(-1.1870556) q[2];
rz(3.0211957) q[3];
sx q[3];
rz(-1.7359366) q[3];
sx q[3];
rz(-0.19336893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.3099986) q[0];
sx q[0];
rz(-3.031142) q[0];
sx q[0];
rz(-0.67685342) q[0];
rz(-1.5325158) q[1];
sx q[1];
rz(-1.1108421) q[1];
sx q[1];
rz(-2.9422876) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53056419) q[0];
sx q[0];
rz(-2.8724551) q[0];
sx q[0];
rz(-1.2236973) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3466269) q[2];
sx q[2];
rz(-1.4557654) q[2];
sx q[2];
rz(0.677106) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3518954) q[1];
sx q[1];
rz(-1.9421862) q[1];
sx q[1];
rz(1.2523734) q[1];
x q[2];
rz(-1.5478743) q[3];
sx q[3];
rz(-0.27697744) q[3];
sx q[3];
rz(0.92844005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4946332) q[2];
sx q[2];
rz(-0.55611742) q[2];
sx q[2];
rz(0.49590084) q[2];
rz(2.2206709) q[3];
sx q[3];
rz(-2.0995188) q[3];
sx q[3];
rz(0.31057772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50716901) q[0];
sx q[0];
rz(-1.2260219) q[0];
sx q[0];
rz(-0.41967151) q[0];
rz(0.93027973) q[1];
sx q[1];
rz(-2.1280961) q[1];
sx q[1];
rz(1.0118265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8735805) q[0];
sx q[0];
rz(-2.4873864) q[0];
sx q[0];
rz(3.0604721) q[0];
rz(-pi) q[1];
x q[1];
rz(2.490245) q[2];
sx q[2];
rz(-0.22964165) q[2];
sx q[2];
rz(-2.1852213) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8817084) q[1];
sx q[1];
rz(-2.1182503) q[1];
sx q[1];
rz(0.24357067) q[1];
rz(-pi) q[2];
rz(-1.3023754) q[3];
sx q[3];
rz(-2.5277349) q[3];
sx q[3];
rz(-0.74210244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73152995) q[2];
sx q[2];
rz(-0.56649929) q[2];
sx q[2];
rz(-1.1140964) q[2];
rz(-1.7690432) q[3];
sx q[3];
rz(-2.1362344) q[3];
sx q[3];
rz(0.70926386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1019185) q[0];
sx q[0];
rz(-1.8721975) q[0];
sx q[0];
rz(1.4920379) q[0];
rz(-0.44293013) q[1];
sx q[1];
rz(-1.8699173) q[1];
sx q[1];
rz(2.2992004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5076601) q[0];
sx q[0];
rz(-2.2325987) q[0];
sx q[0];
rz(0.71983257) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9059431) q[2];
sx q[2];
rz(-0.69219263) q[2];
sx q[2];
rz(-0.38470925) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2017224) q[1];
sx q[1];
rz(-2.2186167) q[1];
sx q[1];
rz(0.52380162) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1792036) q[3];
sx q[3];
rz(-0.11193724) q[3];
sx q[3];
rz(0.051210784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3288021) q[2];
sx q[2];
rz(-0.91965681) q[2];
sx q[2];
rz(0.021050464) q[2];
rz(-0.99714315) q[3];
sx q[3];
rz(-2.5992584) q[3];
sx q[3];
rz(-1.5890582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.3461304) q[0];
sx q[0];
rz(-1.902245) q[0];
sx q[0];
rz(1.1783266) q[0];
rz(-1.0109674) q[1];
sx q[1];
rz(-1.4821472) q[1];
sx q[1];
rz(-2.7706026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3431992) q[0];
sx q[0];
rz(-1.4552552) q[0];
sx q[0];
rz(0.2327524) q[0];
x q[1];
rz(1.7395606) q[2];
sx q[2];
rz(-1.4764364) q[2];
sx q[2];
rz(0.33359087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4996262) q[1];
sx q[1];
rz(-1.8607983) q[1];
sx q[1];
rz(0.15695842) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8274687) q[3];
sx q[3];
rz(-1.1289136) q[3];
sx q[3];
rz(-1.8754547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3498938) q[2];
sx q[2];
rz(-1.0827622) q[2];
sx q[2];
rz(2.2278111) q[2];
rz(-2.6036116) q[3];
sx q[3];
rz(-0.83298433) q[3];
sx q[3];
rz(-0.37916455) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0898042) q[0];
sx q[0];
rz(-1.3642949) q[0];
sx q[0];
rz(0.23189932) q[0];
rz(-1.505835) q[1];
sx q[1];
rz(-1.6927405) q[1];
sx q[1];
rz(2.4969782) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8763257) q[0];
sx q[0];
rz(-2.647222) q[0];
sx q[0];
rz(-0.78440023) q[0];
x q[1];
rz(3.1050131) q[2];
sx q[2];
rz(-2.5365005) q[2];
sx q[2];
rz(0.99672752) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7773103) q[1];
sx q[1];
rz(-0.16723086) q[1];
sx q[1];
rz(-1.4521397) q[1];
rz(-1.1357952) q[3];
sx q[3];
rz(-2.2251872) q[3];
sx q[3];
rz(1.1799348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62442786) q[2];
sx q[2];
rz(-1.6021873) q[2];
sx q[2];
rz(-1.2702764) q[2];
rz(2.7020448) q[3];
sx q[3];
rz(-0.63616532) q[3];
sx q[3];
rz(-0.51072454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9257833) q[0];
sx q[0];
rz(-1.9010811) q[0];
sx q[0];
rz(-1.0260169) q[0];
rz(2.6846474) q[1];
sx q[1];
rz(-2.2508855) q[1];
sx q[1];
rz(2.2073726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6935869) q[0];
sx q[0];
rz(-1.7035055) q[0];
sx q[0];
rz(-0.42610077) q[0];
rz(-pi) q[1];
rz(-1.9719741) q[2];
sx q[2];
rz(-2.8658762) q[2];
sx q[2];
rz(1.3324312) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0353552) q[1];
sx q[1];
rz(-2.3021973) q[1];
sx q[1];
rz(1.3421571) q[1];
rz(-0.85984965) q[3];
sx q[3];
rz(-0.78892869) q[3];
sx q[3];
rz(-3.0956059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.702174) q[2];
sx q[2];
rz(-0.69434035) q[2];
sx q[2];
rz(1.0609421) q[2];
rz(-2.151978) q[3];
sx q[3];
rz(-1.7695534) q[3];
sx q[3];
rz(-0.26452479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2117352) q[0];
sx q[0];
rz(-1.095357) q[0];
sx q[0];
rz(2.3768421) q[0];
rz(1.1872956) q[1];
sx q[1];
rz(-1.7056414) q[1];
sx q[1];
rz(2.5797226) q[1];
rz(-1.5188206) q[2];
sx q[2];
rz(-1.5904273) q[2];
sx q[2];
rz(-1.2640819) q[2];
rz(2.3791947) q[3];
sx q[3];
rz(-2.5350606) q[3];
sx q[3];
rz(1.4192941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
