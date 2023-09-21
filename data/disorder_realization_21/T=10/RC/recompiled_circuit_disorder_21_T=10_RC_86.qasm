OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(-0.50146377) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(-2.8191541) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0599521) q[0];
sx q[0];
rz(-1.8917221) q[0];
sx q[0];
rz(3.0517464) q[0];
rz(-pi) q[1];
rz(-2.2519977) q[2];
sx q[2];
rz(-1.0729562) q[2];
sx q[2];
rz(-1.3537458) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8391708) q[1];
sx q[1];
rz(-1.9170554) q[1];
sx q[1];
rz(-2.195921) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93482165) q[3];
sx q[3];
rz(-1.9536195) q[3];
sx q[3];
rz(-0.32360199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(-0.83267823) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(-0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(-2.9843176) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(-3.0325586) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1125688) q[0];
sx q[0];
rz(-1.7984263) q[0];
sx q[0];
rz(-2.9665222) q[0];
rz(1.3729587) q[2];
sx q[2];
rz(-1.3126144) q[2];
sx q[2];
rz(0.82222647) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2687159) q[1];
sx q[1];
rz(-2.0375588) q[1];
sx q[1];
rz(-0.66653911) q[1];
x q[2];
rz(2.845329) q[3];
sx q[3];
rz(-2.6014572) q[3];
sx q[3];
rz(1.6903282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9033501) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(1.2634574) q[2];
rz(-2.8144828) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(-1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064421244) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(-1.3431312) q[0];
rz(2.893977) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(-0.48167357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1789078) q[0];
sx q[0];
rz(-2.509153) q[0];
sx q[0];
rz(0.90895598) q[0];
rz(0.33294296) q[2];
sx q[2];
rz(-0.98368401) q[2];
sx q[2];
rz(1.3981896) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86320089) q[1];
sx q[1];
rz(-2.26483) q[1];
sx q[1];
rz(-2.4064526) q[1];
rz(-pi) q[2];
rz(-2.5211469) q[3];
sx q[3];
rz(-0.9471604) q[3];
sx q[3];
rz(-2.7321531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3383011) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(-0.49989191) q[2];
rz(2.5806184) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(-1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9445779) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(0.552185) q[0];
rz(-1.5532956) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(-1.8968556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23918505) q[0];
sx q[0];
rz(-0.33948487) q[0];
sx q[0];
rz(0.0048588077) q[0];
x q[1];
rz(-1.3696026) q[2];
sx q[2];
rz(-1.6987213) q[2];
sx q[2];
rz(0.97845562) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1340449) q[1];
sx q[1];
rz(-1.3589077) q[1];
sx q[1];
rz(-2.120554) q[1];
x q[2];
rz(2.7987715) q[3];
sx q[3];
rz(-2.5431513) q[3];
sx q[3];
rz(-0.10859057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2923979) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(-1.1506895) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0634336) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(-3.0601236) q[0];
rz(-0.062462656) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(1.6385471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3969288) q[0];
sx q[0];
rz(-1.6984807) q[0];
sx q[0];
rz(0.98608195) q[0];
rz(-2.224515) q[2];
sx q[2];
rz(-3.0055025) q[2];
sx q[2];
rz(-1.2166785) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7170143) q[1];
sx q[1];
rz(-1.6728405) q[1];
sx q[1];
rz(1.2574408) q[1];
rz(1.8893858) q[3];
sx q[3];
rz(-2.0584403) q[3];
sx q[3];
rz(-2.3209751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.5082671) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(1.2379237) q[2];
rz(2.0189019) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(0.52156633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313107) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(2.8748728) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(0.7985324) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7245367) q[0];
sx q[0];
rz(-2.8310611) q[0];
sx q[0];
rz(2.0470371) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.550225) q[2];
sx q[2];
rz(-0.57833507) q[2];
sx q[2];
rz(0.27507281) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3244464) q[1];
sx q[1];
rz(-1.6231316) q[1];
sx q[1];
rz(0.99919341) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9540326) q[3];
sx q[3];
rz(-2.2700689) q[3];
sx q[3];
rz(-1.9642533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9810527) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(2.7748761) q[2];
rz(-1.8803053) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(-0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40903184) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(2.5352056) q[0];
rz(2.9442893) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(2.6775449) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7921917) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(-2.6130136) q[0];
x q[1];
rz(1.6527962) q[2];
sx q[2];
rz(-0.24157761) q[2];
sx q[2];
rz(-1.8919924) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11394994) q[1];
sx q[1];
rz(-2.0538035) q[1];
sx q[1];
rz(2.7999858) q[1];
rz(-pi) q[2];
rz(0.23630996) q[3];
sx q[3];
rz(-1.650562) q[3];
sx q[3];
rz(-0.41111708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(2.8835473) q[2];
rz(1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(-0.38129693) q[0];
rz(3.0463468) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(1.7000748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7521116) q[0];
sx q[0];
rz(-1.6342667) q[0];
sx q[0];
rz(1.585929) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5289815) q[2];
sx q[2];
rz(-1.5393886) q[2];
sx q[2];
rz(2.999246) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.085233363) q[1];
sx q[1];
rz(-1.7006405) q[1];
sx q[1];
rz(1.6972046) q[1];
rz(-pi) q[2];
rz(0.12318792) q[3];
sx q[3];
rz(-0.33629575) q[3];
sx q[3];
rz(1.5598701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1429446) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(-2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3806234) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.7425849) q[0];
rz(0.31708583) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(0.98659602) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2664468) q[0];
sx q[0];
rz(-0.52030021) q[0];
sx q[0];
rz(-2.8199024) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2676342) q[2];
sx q[2];
rz(-1.1468337) q[2];
sx q[2];
rz(-1.6800113) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0341067) q[1];
sx q[1];
rz(-1.429538) q[1];
sx q[1];
rz(-2.5320401) q[1];
rz(0.17296965) q[3];
sx q[3];
rz(-1.9860455) q[3];
sx q[3];
rz(-3.1239307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9265147) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(-0.40965664) q[2];
rz(-0.26327291) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(2.6221361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(1.7364527) q[0];
rz(-0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.4155037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6575359) q[0];
sx q[0];
rz(-0.41175479) q[0];
sx q[0];
rz(0.62396892) q[0];
rz(2.892898) q[2];
sx q[2];
rz(-1.9242312) q[2];
sx q[2];
rz(2.8949708) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.43511697) q[1];
sx q[1];
rz(-0.26285989) q[1];
sx q[1];
rz(-2.190522) q[1];
rz(-2.9726082) q[3];
sx q[3];
rz(-2.0770693) q[3];
sx q[3];
rz(-0.81912012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4225509) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(-0.12410513) q[2];
rz(-2.1758046) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(-0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.538095) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(2.8812257) q[2];
sx q[2];
rz(-1.7454864) q[2];
sx q[2];
rz(-2.7228552) q[2];
rz(1.7210759) q[3];
sx q[3];
rz(-1.0192623) q[3];
sx q[3];
rz(-2.7912959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
