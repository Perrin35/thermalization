OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90836877) q[0];
sx q[0];
rz(4.6908661) q[0];
sx q[0];
rz(9.6022845) q[0];
rz(1.5597875) q[1];
sx q[1];
rz(4.4686488) q[1];
sx q[1];
rz(8.2460587) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27025199) q[0];
sx q[0];
rz(-2.1553023) q[0];
sx q[0];
rz(2.9181705) q[0];
x q[1];
rz(1.1792405) q[2];
sx q[2];
rz(-0.68049496) q[2];
sx q[2];
rz(0.032236207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70397729) q[1];
sx q[1];
rz(-2.0052138) q[1];
sx q[1];
rz(-2.2056286) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13362515) q[3];
sx q[3];
rz(-2.3006064) q[3];
sx q[3];
rz(-1.3624024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.061964758) q[2];
sx q[2];
rz(-0.92755932) q[2];
sx q[2];
rz(2.8947042) q[2];
rz(0.43761349) q[3];
sx q[3];
rz(-2.103431) q[3];
sx q[3];
rz(2.8810697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-3.0831809) q[0];
sx q[0];
rz(-0.48342466) q[0];
sx q[0];
rz(-1.1046945) q[0];
rz(2.3339234) q[1];
sx q[1];
rz(-2.3876815) q[1];
sx q[1];
rz(2.5025936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.047303) q[0];
sx q[0];
rz(-1.2942137) q[0];
sx q[0];
rz(-1.7533789) q[0];
rz(-pi) q[1];
rz(-2.866687) q[2];
sx q[2];
rz(-1.5228809) q[2];
sx q[2];
rz(2.3570182) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39589992) q[1];
sx q[1];
rz(-1.9338765) q[1];
sx q[1];
rz(-2.3895742) q[1];
rz(-pi) q[2];
rz(1.9023537) q[3];
sx q[3];
rz(-1.8428118) q[3];
sx q[3];
rz(0.60604688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4421926) q[2];
sx q[2];
rz(-1.3905808) q[2];
sx q[2];
rz(-0.81713444) q[2];
rz(0.52554956) q[3];
sx q[3];
rz(-0.81997973) q[3];
sx q[3];
rz(1.9513963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3007871) q[0];
sx q[0];
rz(-1.5910933) q[0];
sx q[0];
rz(-2.9817885) q[0];
rz(0.5197168) q[1];
sx q[1];
rz(-1.0451885) q[1];
sx q[1];
rz(0.91667485) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0063112886) q[0];
sx q[0];
rz(-1.0005299) q[0];
sx q[0];
rz(-1.6018014) q[0];
rz(-pi) q[1];
rz(2.4009928) q[2];
sx q[2];
rz(-1.3889905) q[2];
sx q[2];
rz(-1.2642853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8163678) q[1];
sx q[1];
rz(-1.381641) q[1];
sx q[1];
rz(-0.048859096) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7052207) q[3];
sx q[3];
rz(-2.1791576) q[3];
sx q[3];
rz(-2.9384263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2436169) q[2];
sx q[2];
rz(-0.50668442) q[2];
sx q[2];
rz(-1.766073) q[2];
rz(-3.1084133) q[3];
sx q[3];
rz(-1.5539919) q[3];
sx q[3];
rz(2.3110265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0631977) q[0];
sx q[0];
rz(-2.4802408) q[0];
sx q[0];
rz(-0.98264328) q[0];
rz(0.64535087) q[1];
sx q[1];
rz(-1.2256365) q[1];
sx q[1];
rz(2.7243848) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0785976) q[0];
sx q[0];
rz(-1.0630442) q[0];
sx q[0];
rz(-1.2247472) q[0];
rz(-3.0723802) q[2];
sx q[2];
rz(-1.950351) q[2];
sx q[2];
rz(0.57513851) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5886093) q[1];
sx q[1];
rz(-0.89834301) q[1];
sx q[1];
rz(1.49468) q[1];
x q[2];
rz(2.6989815) q[3];
sx q[3];
rz(-0.59429443) q[3];
sx q[3];
rz(-1.8123008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8365606) q[2];
sx q[2];
rz(-0.3362259) q[2];
sx q[2];
rz(-1.1870556) q[2];
rz(3.0211957) q[3];
sx q[3];
rz(-1.4056561) q[3];
sx q[3];
rz(-2.9482237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3099986) q[0];
sx q[0];
rz(-3.031142) q[0];
sx q[0];
rz(-2.4647392) q[0];
rz(-1.6090769) q[1];
sx q[1];
rz(-1.1108421) q[1];
sx q[1];
rz(2.9422876) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6110285) q[0];
sx q[0];
rz(-0.26913759) q[0];
sx q[0];
rz(-1.9178954) q[0];
rz(-2.8137652) q[2];
sx q[2];
rz(-0.36448267) q[2];
sx q[2];
rz(2.5555697) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0892895) q[1];
sx q[1];
rz(-0.48434231) q[1];
sx q[1];
rz(2.4645096) q[1];
rz(-3.1350769) q[3];
sx q[3];
rz(-1.8476991) q[3];
sx q[3];
rz(-2.2369825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6469595) q[2];
sx q[2];
rz(-0.55611742) q[2];
sx q[2];
rz(-0.49590084) q[2];
rz(-0.9209218) q[3];
sx q[3];
rz(-2.0995188) q[3];
sx q[3];
rz(0.31057772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6344236) q[0];
sx q[0];
rz(-1.9155707) q[0];
sx q[0];
rz(-2.7219211) q[0];
rz(2.2113129) q[1];
sx q[1];
rz(-2.1280961) q[1];
sx q[1];
rz(2.1297661) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8735805) q[0];
sx q[0];
rz(-0.65420628) q[0];
sx q[0];
rz(-3.0604721) q[0];
x q[1];
rz(-2.957785) q[2];
sx q[2];
rz(-1.70924) q[2];
sx q[2];
rz(-3.1174497) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9560555) q[1];
sx q[1];
rz(-2.5474869) q[1];
sx q[1];
rz(-1.1940765) q[1];
rz(-pi) q[2];
rz(-0.18475549) q[3];
sx q[3];
rz(-2.1596381) q[3];
sx q[3];
rz(2.7240804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73152995) q[2];
sx q[2];
rz(-0.56649929) q[2];
sx q[2];
rz(2.0274963) q[2];
rz(1.3725494) q[3];
sx q[3];
rz(-1.0053582) q[3];
sx q[3];
rz(-0.70926386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0396742) q[0];
sx q[0];
rz(-1.8721975) q[0];
sx q[0];
rz(-1.4920379) q[0];
rz(-2.6986625) q[1];
sx q[1];
rz(-1.8699173) q[1];
sx q[1];
rz(-2.2992004) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5076601) q[0];
sx q[0];
rz(-0.90899399) q[0];
sx q[0];
rz(2.4217601) q[0];
rz(-1.2356495) q[2];
sx q[2];
rz(-0.69219263) q[2];
sx q[2];
rz(0.38470925) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7034027) q[1];
sx q[1];
rz(-2.3329607) q[1];
sx q[1];
rz(-0.98677294) q[1];
x q[2];
rz(-1.1792036) q[3];
sx q[3];
rz(-0.11193724) q[3];
sx q[3];
rz(0.051210784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.81279057) q[2];
sx q[2];
rz(-2.2219358) q[2];
sx q[2];
rz(-3.1205422) q[2];
rz(-0.99714315) q[3];
sx q[3];
rz(-0.54233426) q[3];
sx q[3];
rz(-1.5525345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7954623) q[0];
sx q[0];
rz(-1.902245) q[0];
sx q[0];
rz(1.963266) q[0];
rz(2.1306253) q[1];
sx q[1];
rz(-1.4821472) q[1];
sx q[1];
rz(-2.7706026) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31974935) q[0];
sx q[0];
rz(-2.8822063) q[0];
sx q[0];
rz(0.46617561) q[0];
rz(-pi) q[1];
rz(1.402032) q[2];
sx q[2];
rz(-1.4764364) q[2];
sx q[2];
rz(-0.33359087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2579871) q[1];
sx q[1];
rz(-1.7211497) q[1];
sx q[1];
rz(-1.8642051) q[1];
rz(-pi) q[2];
rz(-1.1092126) q[3];
sx q[3];
rz(-1.8538665) q[3];
sx q[3];
rz(0.1666094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79169881) q[2];
sx q[2];
rz(-1.0827622) q[2];
sx q[2];
rz(-2.2278111) q[2];
rz(-0.53798109) q[3];
sx q[3];
rz(-0.83298433) q[3];
sx q[3];
rz(0.37916455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0898042) q[0];
sx q[0];
rz(-1.3642949) q[0];
sx q[0];
rz(2.9096933) q[0];
rz(1.6357577) q[1];
sx q[1];
rz(-1.4488522) q[1];
sx q[1];
rz(0.64461446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.028325) q[0];
sx q[0];
rz(-1.2282983) q[0];
sx q[0];
rz(-1.934608) q[0];
rz(-pi) q[1];
rz(0.03657954) q[2];
sx q[2];
rz(-2.5365005) q[2];
sx q[2];
rz(2.1448651) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.48460173) q[1];
sx q[1];
rz(-1.7368403) q[1];
sx q[1];
rz(3.1216122) q[1];
rz(-1.1357952) q[3];
sx q[3];
rz(-0.91640546) q[3];
sx q[3];
rz(1.9616579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5171648) q[2];
sx q[2];
rz(-1.5394053) q[2];
sx q[2];
rz(1.2702764) q[2];
rz(-2.7020448) q[3];
sx q[3];
rz(-0.63616532) q[3];
sx q[3];
rz(-2.6308681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21580932) q[0];
sx q[0];
rz(-1.2405115) q[0];
sx q[0];
rz(-1.0260169) q[0];
rz(0.45694524) q[1];
sx q[1];
rz(-2.2508855) q[1];
sx q[1];
rz(-2.2073726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1608604) q[0];
sx q[0];
rz(-2.6965158) q[0];
sx q[0];
rz(-2.829192) q[0];
rz(-pi) q[1];
x q[1];
rz(3.031557) q[2];
sx q[2];
rz(-1.8241183) q[2];
sx q[2];
rz(-0.91722721) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9125713) q[1];
sx q[1];
rz(-2.3816436) q[1];
sx q[1];
rz(-0.2473803) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85984965) q[3];
sx q[3];
rz(-2.352664) q[3];
sx q[3];
rz(0.045986764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.702174) q[2];
sx q[2];
rz(-2.4472523) q[2];
sx q[2];
rz(2.0806506) q[2];
rz(0.98961467) q[3];
sx q[3];
rz(-1.7695534) q[3];
sx q[3];
rz(-0.26452479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2117352) q[0];
sx q[0];
rz(-1.095357) q[0];
sx q[0];
rz(2.3768421) q[0];
rz(-1.9542971) q[1];
sx q[1];
rz(-1.7056414) q[1];
sx q[1];
rz(2.5797226) q[1];
rz(-1.9321185) q[2];
sx q[2];
rz(-0.055556297) q[2];
sx q[2];
rz(-0.054097492) q[2];
rz(-0.76239794) q[3];
sx q[3];
rz(-2.5350606) q[3];
sx q[3];
rz(1.4192941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
