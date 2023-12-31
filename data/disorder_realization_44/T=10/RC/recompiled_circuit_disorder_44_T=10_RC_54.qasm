OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(6.7232806) q[0];
sx q[0];
rz(6.4203782) q[0];
rz(-1.7358915) q[1];
sx q[1];
rz(-1.403221) q[1];
sx q[1];
rz(-0.52991968) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0041381) q[0];
sx q[0];
rz(-1.19085) q[0];
sx q[0];
rz(-0.11560346) q[0];
rz(-0.86962236) q[2];
sx q[2];
rz(-2.6343971) q[2];
sx q[2];
rz(1.6091572) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4094028) q[1];
sx q[1];
rz(-0.70530546) q[1];
sx q[1];
rz(-2.3975055) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8205809) q[3];
sx q[3];
rz(-0.82818177) q[3];
sx q[3];
rz(3.0299377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68937504) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-0.33660647) q[2];
rz(-1.5161139) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(-1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7933554) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(-0.021214699) q[0];
rz(-1.1938098) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(0.83591998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7556691) q[0];
sx q[0];
rz(-3.0015411) q[0];
sx q[0];
rz(-1.7007909) q[0];
rz(-pi) q[1];
rz(-2.1199273) q[2];
sx q[2];
rz(-1.9252535) q[2];
sx q[2];
rz(-0.82565386) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4513431) q[1];
sx q[1];
rz(-0.93938821) q[1];
sx q[1];
rz(-2.8020225) q[1];
x q[2];
rz(1.9637945) q[3];
sx q[3];
rz(-0.8702232) q[3];
sx q[3];
rz(-2.8515479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2772284) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(-1.345984) q[2];
rz(2.7820382) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(-2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(1.0536449) q[0];
rz(-1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(2.704481) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29262221) q[0];
sx q[0];
rz(-1.4842352) q[0];
sx q[0];
rz(1.863088) q[0];
x q[1];
rz(-1.1630467) q[2];
sx q[2];
rz(-1.3464658) q[2];
sx q[2];
rz(0.11220223) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8846109) q[1];
sx q[1];
rz(-1.4323438) q[1];
sx q[1];
rz(-2.715766) q[1];
rz(-pi) q[2];
rz(-1.0760355) q[3];
sx q[3];
rz(-1.6238188) q[3];
sx q[3];
rz(0.094129063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.019471021) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(2.1195228) q[2];
rz(-1.2381037) q[3];
sx q[3];
rz(-0.3823897) q[3];
sx q[3];
rz(-2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5220752) q[0];
sx q[0];
rz(-1.2457122) q[0];
sx q[0];
rz(0.98130256) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(-0.19128004) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45368886) q[0];
sx q[0];
rz(-2.9691302) q[0];
sx q[0];
rz(-1.0426636) q[0];
rz(1.26881) q[2];
sx q[2];
rz(-0.75690818) q[2];
sx q[2];
rz(0.41727558) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3595708) q[1];
sx q[1];
rz(-1.6732209) q[1];
sx q[1];
rz(2.3624079) q[1];
rz(-pi) q[2];
rz(2.6077765) q[3];
sx q[3];
rz(-1.2663519) q[3];
sx q[3];
rz(-0.75418762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4613351) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(1.0106687) q[2];
rz(0.7615532) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(-2.9060569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33870944) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(-2.5849735) q[0];
rz(-3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(2.1690878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2302549) q[0];
sx q[0];
rz(-1.5033493) q[0];
sx q[0];
rz(-0.10226843) q[0];
rz(1.2178671) q[2];
sx q[2];
rz(-0.7820411) q[2];
sx q[2];
rz(-2.1176586) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3481808) q[1];
sx q[1];
rz(-0.23985292) q[1];
sx q[1];
rz(-0.94437771) q[1];
x q[2];
rz(-3.1049018) q[3];
sx q[3];
rz(-0.9551691) q[3];
sx q[3];
rz(-0.43945593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82289034) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(1.548432) q[2];
rz(1.7758153) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4218629) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(-1.7156037) q[0];
rz(-1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(2.7672966) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22898856) q[0];
sx q[0];
rz(-1.8870263) q[0];
sx q[0];
rz(-0.4321179) q[0];
rz(-pi) q[1];
rz(0.46005581) q[2];
sx q[2];
rz(-0.54121491) q[2];
sx q[2];
rz(-2.7341026) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7172076) q[1];
sx q[1];
rz(-1.3076412) q[1];
sx q[1];
rz(-0.56916635) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88126392) q[3];
sx q[3];
rz(-2.1173819) q[3];
sx q[3];
rz(-1.5342086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2465683) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(2.1833615) q[2];
rz(-2.9124177) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(0.62098256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(0.36528698) q[0];
sx q[0];
rz(-1.9488652) q[0];
sx q[0];
rz(0.90674415) q[0];
rz(-1.0892185) q[1];
sx q[1];
rz(-1.6420495) q[1];
sx q[1];
rz(1.8315171) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.487264) q[0];
sx q[0];
rz(-1.9582821) q[0];
sx q[0];
rz(1.6261473) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4285165) q[2];
sx q[2];
rz(-1.3165858) q[2];
sx q[2];
rz(1.6377246) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8430427) q[1];
sx q[1];
rz(-1.3805461) q[1];
sx q[1];
rz(-2.802554) q[1];
rz(1.5982315) q[3];
sx q[3];
rz(-1.0155639) q[3];
sx q[3];
rz(2.1813986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2157796) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(-1.6581992) q[2];
rz(0.27967927) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(-3.083995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9782372) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(-0.21959198) q[0];
rz(2.638468) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(0.84987744) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6760611) q[0];
sx q[0];
rz(-1.413835) q[0];
sx q[0];
rz(1.3423052) q[0];
rz(-2.2539027) q[2];
sx q[2];
rz(-1.6437093) q[2];
sx q[2];
rz(-0.86276744) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1666607) q[1];
sx q[1];
rz(-0.76172511) q[1];
sx q[1];
rz(1.5207661) q[1];
rz(-2.0073118) q[3];
sx q[3];
rz(-2.8562162) q[3];
sx q[3];
rz(2.7989822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59051096) q[2];
sx q[2];
rz(-1.3122281) q[2];
sx q[2];
rz(-1.760651) q[2];
rz(0.75602174) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(-2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5091771) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(0.40503043) q[0];
rz(-2.6889154) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(1.2776432) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4687846) q[0];
sx q[0];
rz(-1.0793669) q[0];
sx q[0];
rz(1.3093033) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69127609) q[2];
sx q[2];
rz(-1.1862159) q[2];
sx q[2];
rz(2.861475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8511508) q[1];
sx q[1];
rz(-0.44712198) q[1];
sx q[1];
rz(-1.6563583) q[1];
rz(-0.38655917) q[3];
sx q[3];
rz(-0.87214008) q[3];
sx q[3];
rz(-0.78702918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8032802) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(-0.35783106) q[2];
rz(-1.7221649) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2232067) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(0.11225587) q[0];
rz(0.90011251) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(2.9311438) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8716547) q[0];
sx q[0];
rz(-0.5408113) q[0];
sx q[0];
rz(-0.35738118) q[0];
rz(-1.9824355) q[2];
sx q[2];
rz(-1.4201846) q[2];
sx q[2];
rz(2.6580236) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7604954) q[1];
sx q[1];
rz(-1.7777182) q[1];
sx q[1];
rz(2.8269672) q[1];
rz(-pi) q[2];
rz(0.25037346) q[3];
sx q[3];
rz(-1.3406521) q[3];
sx q[3];
rz(-1.9104513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.18008733) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(1.5853184) q[2];
rz(-1.2735584) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(0.20726985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56959854) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(-2.3251484) q[1];
sx q[1];
rz(-1.888231) q[1];
sx q[1];
rz(2.9838557) q[1];
rz(0.16602892) q[2];
sx q[2];
rz(-2.0353073) q[2];
sx q[2];
rz(-1.4370949) q[2];
rz(-1.9498701) q[3];
sx q[3];
rz(-0.76615292) q[3];
sx q[3];
rz(-2.5828008) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
