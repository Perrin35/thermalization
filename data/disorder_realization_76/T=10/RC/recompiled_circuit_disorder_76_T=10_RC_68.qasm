OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27622142) q[0];
sx q[0];
rz(-0.85715357) q[0];
sx q[0];
rz(0.13248086) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(-0.58499709) q[1];
sx q[1];
rz(2.4490228) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6462631) q[0];
sx q[0];
rz(-1.4709934) q[0];
sx q[0];
rz(-2.0908337) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1039256) q[2];
sx q[2];
rz(-0.97969998) q[2];
sx q[2];
rz(-2.1819654) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9308656) q[1];
sx q[1];
rz(-2.9949246) q[1];
sx q[1];
rz(-1.055483) q[1];
x q[2];
rz(0.2487189) q[3];
sx q[3];
rz(-2.0318444) q[3];
sx q[3];
rz(0.15795262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1074368) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.6050603) q[2];
rz(-1.6202554) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81543106) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(-1.8923627) q[0];
rz(-0.56150395) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(2.5610279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626494) q[0];
sx q[0];
rz(-1.8793086) q[0];
sx q[0];
rz(1.6842151) q[0];
x q[1];
rz(-0.29082362) q[2];
sx q[2];
rz(-0.65867701) q[2];
sx q[2];
rz(2.3386699) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9241582) q[1];
sx q[1];
rz(-2.1719451) q[1];
sx q[1];
rz(2.2648328) q[1];
rz(-pi) q[2];
rz(1.1167691) q[3];
sx q[3];
rz(-1.2516216) q[3];
sx q[3];
rz(-1.9516731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7291752) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(-2.2632329) q[2];
rz(0.39204028) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(2.5382606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47675258) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-0.77366775) q[0];
rz(3.1402918) q[1];
sx q[1];
rz(-1.5258077) q[1];
sx q[1];
rz(3.1087648) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5531909) q[0];
sx q[0];
rz(-1.2598039) q[0];
sx q[0];
rz(-1.3283967) q[0];
x q[1];
rz(0.29067729) q[2];
sx q[2];
rz(-2.103984) q[2];
sx q[2];
rz(2.3454587) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.58363885) q[1];
sx q[1];
rz(-1.7350405) q[1];
sx q[1];
rz(-2.3418531) q[1];
rz(-2.114931) q[3];
sx q[3];
rz(-1.2715724) q[3];
sx q[3];
rz(2.2194089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34439987) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(-0.27734217) q[2];
rz(-0.39595655) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(-2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70401496) q[0];
sx q[0];
rz(-2.8058348) q[0];
sx q[0];
rz(-2.8787956) q[0];
rz(-2.908005) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(2.3707726) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3931657) q[0];
sx q[0];
rz(-1.59032) q[0];
sx q[0];
rz(0.016419134) q[0];
x q[1];
rz(-2.0609444) q[2];
sx q[2];
rz(-1.4785826) q[2];
sx q[2];
rz(2.8871418) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.96344906) q[1];
sx q[1];
rz(-1.932486) q[1];
sx q[1];
rz(0.55808918) q[1];
rz(1.3464438) q[3];
sx q[3];
rz(-1.7566578) q[3];
sx q[3];
rz(2.2546774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9757441) q[2];
sx q[2];
rz(-1.6341012) q[2];
sx q[2];
rz(0.7129933) q[2];
rz(-2.1285848) q[3];
sx q[3];
rz(-0.37390798) q[3];
sx q[3];
rz(-0.95389429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35904303) q[0];
sx q[0];
rz(-1.0721711) q[0];
sx q[0];
rz(1.7011401) q[0];
rz(3.0474995) q[1];
sx q[1];
rz(-2.4021939) q[1];
sx q[1];
rz(-2.9715911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48430303) q[0];
sx q[0];
rz(-0.7938677) q[0];
sx q[0];
rz(-2.8110709) q[0];
rz(2.0001569) q[2];
sx q[2];
rz(-2.6300207) q[2];
sx q[2];
rz(-2.9432952) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.070378455) q[1];
sx q[1];
rz(-2.4134679) q[1];
sx q[1];
rz(0.94707625) q[1];
rz(-pi) q[2];
rz(2.2682297) q[3];
sx q[3];
rz(-2.1057099) q[3];
sx q[3];
rz(-1.3267335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1468982) q[2];
sx q[2];
rz(-2.6165104) q[2];
sx q[2];
rz(-1.4040995) q[2];
rz(1.5935625) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(-1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(2.6830542) q[0];
rz(-2.8857152) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(2.4564254) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93668324) q[0];
sx q[0];
rz(-0.65403599) q[0];
sx q[0];
rz(1.4564287) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14641996) q[2];
sx q[2];
rz(-0.78260566) q[2];
sx q[2];
rz(-2.866982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.203478) q[1];
sx q[1];
rz(-0.96760975) q[1];
sx q[1];
rz(0.60738648) q[1];
rz(0.018304304) q[3];
sx q[3];
rz(-0.98494512) q[3];
sx q[3];
rz(-0.5154807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7489862) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(-1.4292599) q[2];
rz(2.0424992) q[3];
sx q[3];
rz(-0.50656879) q[3];
sx q[3];
rz(-0.18923047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1290865) q[0];
sx q[0];
rz(-1.5456454) q[0];
sx q[0];
rz(-2.4122453) q[0];
rz(2.8485281) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(-1.1475295) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.731819) q[0];
sx q[0];
rz(-1.6454576) q[0];
sx q[0];
rz(0.41418196) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93716623) q[2];
sx q[2];
rz(-1.2680149) q[2];
sx q[2];
rz(-2.5123951) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9784669) q[1];
sx q[1];
rz(-1.6118057) q[1];
sx q[1];
rz(2.6569215) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64671867) q[3];
sx q[3];
rz(-2.5463856) q[3];
sx q[3];
rz(-1.7531542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3532233) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(-2.3925171) q[2];
rz(-2.4979112) q[3];
sx q[3];
rz(-1.0130853) q[3];
sx q[3];
rz(2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1483243) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(-0.83129445) q[0];
rz(1.7656901) q[1];
sx q[1];
rz(-2.3283236) q[1];
sx q[1];
rz(-2.7430699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7300028) q[0];
sx q[0];
rz(-1.411479) q[0];
sx q[0];
rz(-2.0429862) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6767098) q[2];
sx q[2];
rz(-0.97230655) q[2];
sx q[2];
rz(-0.61818365) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9989222) q[1];
sx q[1];
rz(-0.97375662) q[1];
sx q[1];
rz(-3.1139657) q[1];
x q[2];
rz(-2.1274444) q[3];
sx q[3];
rz(-1.1245407) q[3];
sx q[3];
rz(0.92170148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1901671) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(1.297696) q[2];
rz(2.0166345) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(2.4979533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012638906) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(-1.3893611) q[0];
rz(-1.5147491) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(1.0983889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2929045) q[0];
sx q[0];
rz(-1.0875889) q[0];
sx q[0];
rz(2.8781761) q[0];
rz(-pi) q[1];
rz(-2.6838052) q[2];
sx q[2];
rz(-2.2676003) q[2];
sx q[2];
rz(2.2965477) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2201982) q[1];
sx q[1];
rz(-0.99153334) q[1];
sx q[1];
rz(1.5701576) q[1];
x q[2];
rz(1.7421726) q[3];
sx q[3];
rz(-1.3203353) q[3];
sx q[3];
rz(-0.40303883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2417458) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(-1.3174723) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39524233) q[0];
sx q[0];
rz(-2.0210176) q[0];
sx q[0];
rz(0.46646068) q[0];
rz(2.9699504) q[1];
sx q[1];
rz(-1.9263093) q[1];
sx q[1];
rz(2.5126273) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53101978) q[0];
sx q[0];
rz(-1.6197228) q[0];
sx q[0];
rz(-1.1958836) q[0];
x q[1];
rz(1.8877108) q[2];
sx q[2];
rz(-2.4744518) q[2];
sx q[2];
rz(-2.4184879) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.65044636) q[1];
sx q[1];
rz(-1.568734) q[1];
sx q[1];
rz(2.6994929) q[1];
rz(-pi) q[2];
rz(-0.38871308) q[3];
sx q[3];
rz(-2.3832088) q[3];
sx q[3];
rz(-3.0384516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(2.0824599) q[2];
rz(0.6774261) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(-0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.28329904) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(0.25390608) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(2.4181096) q[2];
sx q[2];
rz(-1.6037446) q[2];
sx q[2];
rz(-1.4708191) q[2];
rz(1.4678636) q[3];
sx q[3];
rz(-0.58936215) q[3];
sx q[3];
rz(-2.0274558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
