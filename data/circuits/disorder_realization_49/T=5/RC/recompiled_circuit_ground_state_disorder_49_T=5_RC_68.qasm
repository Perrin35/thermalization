OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9743118) q[0];
sx q[0];
rz(-0.33742961) q[0];
sx q[0];
rz(0.20198241) q[0];
rz(2.1972411) q[1];
sx q[1];
rz(2.386932) q[1];
sx q[1];
rz(7.9805482) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2371191) q[0];
sx q[0];
rz(-0.54330641) q[0];
sx q[0];
rz(-0.41443698) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0927222) q[2];
sx q[2];
rz(-2.9101924) q[2];
sx q[2];
rz(0.3203985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3814427) q[1];
sx q[1];
rz(-2.0408148) q[1];
sx q[1];
rz(-2.2862741) q[1];
x q[2];
rz(-3.0962288) q[3];
sx q[3];
rz(-1.2063908) q[3];
sx q[3];
rz(1.7045329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28295383) q[2];
sx q[2];
rz(-3.0573461) q[2];
sx q[2];
rz(1.6072744) q[2];
rz(0.18680799) q[3];
sx q[3];
rz(-2.6818633) q[3];
sx q[3];
rz(-1.648858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9635791) q[0];
sx q[0];
rz(-2.3602965) q[0];
sx q[0];
rz(2.436893) q[0];
rz(2.1251382) q[1];
sx q[1];
rz(-1.9489138) q[1];
sx q[1];
rz(-1.1526398) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148852) q[0];
sx q[0];
rz(-0.57327691) q[0];
sx q[0];
rz(1.2529729) q[0];
x q[1];
rz(0.57781685) q[2];
sx q[2];
rz(-2.0426828) q[2];
sx q[2];
rz(-3.0834215) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.071653366) q[1];
sx q[1];
rz(-0.80822435) q[1];
sx q[1];
rz(2.3030242) q[1];
rz(-pi) q[2];
x q[2];
rz(2.393778) q[3];
sx q[3];
rz(-1.6666109) q[3];
sx q[3];
rz(-0.15296061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54027259) q[2];
sx q[2];
rz(-2.2728964) q[2];
sx q[2];
rz(-1.8216088) q[2];
rz(0.071062239) q[3];
sx q[3];
rz(-3.0607405) q[3];
sx q[3];
rz(-2.0474056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0367301) q[0];
sx q[0];
rz(-0.34958378) q[0];
sx q[0];
rz(2.2678251) q[0];
rz(1.8050487) q[1];
sx q[1];
rz(-0.68160325) q[1];
sx q[1];
rz(2.2866586) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9062775) q[0];
sx q[0];
rz(-1.7830669) q[0];
sx q[0];
rz(2.2215823) q[0];
x q[1];
rz(-2.8310815) q[2];
sx q[2];
rz(-1.9726255) q[2];
sx q[2];
rz(-2.7257811) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47006059) q[1];
sx q[1];
rz(-1.6796659) q[1];
sx q[1];
rz(2.837869) q[1];
rz(0.35936648) q[3];
sx q[3];
rz(-2.1443261) q[3];
sx q[3];
rz(-1.9952967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9285589) q[2];
sx q[2];
rz(-1.1308257) q[2];
sx q[2];
rz(2.0862759) q[2];
rz(0.94163752) q[3];
sx q[3];
rz(-1.0712653) q[3];
sx q[3];
rz(0.94282237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4554491) q[0];
sx q[0];
rz(-2.2930155) q[0];
sx q[0];
rz(-2.4461179) q[0];
rz(1.6740359) q[1];
sx q[1];
rz(-1.8753884) q[1];
sx q[1];
rz(-3.0470336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3723264) q[0];
sx q[0];
rz(-2.8594198) q[0];
sx q[0];
rz(-1.7442987) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2945061) q[2];
sx q[2];
rz(-1.0930335) q[2];
sx q[2];
rz(-2.1412882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.27457224) q[1];
sx q[1];
rz(-0.34718514) q[1];
sx q[1];
rz(2.1866261) q[1];
rz(-1.8843205) q[3];
sx q[3];
rz(-2.4289968) q[3];
sx q[3];
rz(-0.37171504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62119421) q[2];
sx q[2];
rz(-2.3475671) q[2];
sx q[2];
rz(-2.3801079) q[2];
rz(2.9705808) q[3];
sx q[3];
rz(-1.2820425) q[3];
sx q[3];
rz(3.0448992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347725) q[0];
sx q[0];
rz(-0.67245317) q[0];
sx q[0];
rz(2.6493454) q[0];
rz(-1.8222088) q[1];
sx q[1];
rz(-1.7183813) q[1];
sx q[1];
rz(-2.6548903) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0224441) q[0];
sx q[0];
rz(-2.7239425) q[0];
sx q[0];
rz(2.8432884) q[0];
x q[1];
rz(-1.9453508) q[2];
sx q[2];
rz(-0.83714467) q[2];
sx q[2];
rz(-2.9840699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8913331) q[1];
sx q[1];
rz(-0.57602611) q[1];
sx q[1];
rz(-2.2800619) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23426849) q[3];
sx q[3];
rz(-2.0741012) q[3];
sx q[3];
rz(0.64555321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95185602) q[2];
sx q[2];
rz(-2.2978013) q[2];
sx q[2];
rz(-2.0437415) q[2];
rz(-2.3073933) q[3];
sx q[3];
rz(-2.4662377) q[3];
sx q[3];
rz(1.437291) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6247691) q[0];
sx q[0];
rz(-2.5747445) q[0];
sx q[0];
rz(-0.045850642) q[0];
rz(0.74371964) q[1];
sx q[1];
rz(-0.36879483) q[1];
sx q[1];
rz(-3.0812982) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1566965) q[0];
sx q[0];
rz(-1.9895574) q[0];
sx q[0];
rz(1.9661222) q[0];
rz(-2.0343389) q[2];
sx q[2];
rz(-2.3899979) q[2];
sx q[2];
rz(-1.7558524) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0963147) q[1];
sx q[1];
rz(-1.8966676) q[1];
sx q[1];
rz(-0.23479825) q[1];
rz(-pi) q[2];
rz(-3.0457959) q[3];
sx q[3];
rz(-2.2471273) q[3];
sx q[3];
rz(-0.47190445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22659773) q[2];
sx q[2];
rz(-1.8151585) q[2];
sx q[2];
rz(-2.162852) q[2];
rz(1.8799051) q[3];
sx q[3];
rz(-2.1224969) q[3];
sx q[3];
rz(-0.88920465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17734811) q[0];
sx q[0];
rz(-1.6151936) q[0];
sx q[0];
rz(2.2571795) q[0];
rz(-2.2824967) q[1];
sx q[1];
rz(-2.0180549) q[1];
sx q[1];
rz(-0.20800796) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73586845) q[0];
sx q[0];
rz(-2.5130773) q[0];
sx q[0];
rz(2.5121538) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.154083) q[2];
sx q[2];
rz(-1.0089968) q[2];
sx q[2];
rz(-0.75129189) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.35321843) q[1];
sx q[1];
rz(-0.93868449) q[1];
sx q[1];
rz(2.8115584) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0020387928) q[3];
sx q[3];
rz(-2.0741794) q[3];
sx q[3];
rz(-2.9640405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38006833) q[2];
sx q[2];
rz(-0.91975776) q[2];
sx q[2];
rz(2.5999787) q[2];
rz(1.2540865) q[3];
sx q[3];
rz(-1.5897635) q[3];
sx q[3];
rz(-2.2327173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8054473) q[0];
sx q[0];
rz(-2.8271524) q[0];
sx q[0];
rz(1.0910777) q[0];
rz(-0.82487851) q[1];
sx q[1];
rz(-0.42366091) q[1];
sx q[1];
rz(-2.0687912) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7110853) q[0];
sx q[0];
rz(-0.46784624) q[0];
sx q[0];
rz(0.43456315) q[0];
rz(-3.1284749) q[2];
sx q[2];
rz(-1.5481536) q[2];
sx q[2];
rz(2.4215339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78847105) q[1];
sx q[1];
rz(-1.7040729) q[1];
sx q[1];
rz(0.2443831) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9735224) q[3];
sx q[3];
rz(-1.3478876) q[3];
sx q[3];
rz(0.18833489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9020646) q[2];
sx q[2];
rz(-0.54486474) q[2];
sx q[2];
rz(2.876335) q[2];
rz(-0.15051633) q[3];
sx q[3];
rz(-1.1082606) q[3];
sx q[3];
rz(2.7786541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4758258) q[0];
sx q[0];
rz(-0.095223991) q[0];
sx q[0];
rz(1.6640523) q[0];
rz(-2.3382969) q[1];
sx q[1];
rz(-1.3731615) q[1];
sx q[1];
rz(1.0172179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68960357) q[0];
sx q[0];
rz(-1.3675121) q[0];
sx q[0];
rz(3.0102121) q[0];
rz(-pi) q[1];
rz(-0.4847277) q[2];
sx q[2];
rz(-0.51186168) q[2];
sx q[2];
rz(-2.8163225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78002659) q[1];
sx q[1];
rz(-1.2793555) q[1];
sx q[1];
rz(2.7855278) q[1];
rz(-0.32633304) q[3];
sx q[3];
rz(-1.2830955) q[3];
sx q[3];
rz(-0.74244754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9251755) q[2];
sx q[2];
rz(-1.567013) q[2];
sx q[2];
rz(0.65640059) q[2];
rz(-2.9260855) q[3];
sx q[3];
rz(-1.7301205) q[3];
sx q[3];
rz(1.6925252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086275252) q[0];
sx q[0];
rz(-1.3441939) q[0];
sx q[0];
rz(-2.2917746) q[0];
rz(-2.1814116) q[1];
sx q[1];
rz(-0.4147059) q[1];
sx q[1];
rz(-0.93920952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4681743) q[0];
sx q[0];
rz(-0.8020173) q[0];
sx q[0];
rz(1.9413663) q[0];
x q[1];
rz(1.5884253) q[2];
sx q[2];
rz(-2.7269502) q[2];
sx q[2];
rz(-1.7830242) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1436436) q[1];
sx q[1];
rz(-0.74926361) q[1];
sx q[1];
rz(0.56203385) q[1];
x q[2];
rz(-0.84374238) q[3];
sx q[3];
rz(-0.66202032) q[3];
sx q[3];
rz(2.6204505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5592929) q[2];
sx q[2];
rz(-1.7377995) q[2];
sx q[2];
rz(-0.55029184) q[2];
rz(-3.0136717) q[3];
sx q[3];
rz(-3.0030799) q[3];
sx q[3];
rz(2.1730455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48548231) q[0];
sx q[0];
rz(-2.4714097) q[0];
sx q[0];
rz(2.4708268) q[0];
rz(-0.35519629) q[1];
sx q[1];
rz(-1.4757481) q[1];
sx q[1];
rz(-0.3955985) q[1];
rz(0.14329362) q[2];
sx q[2];
rz(-1.7096277) q[2];
sx q[2];
rz(-1.1451677) q[2];
rz(-2.0216178) q[3];
sx q[3];
rz(-2.6482441) q[3];
sx q[3];
rz(-1.4239414) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
