OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(2.1029544) q[1];
sx q[1];
rz(-1.3488052) q[1];
sx q[1];
rz(-2.2178712) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38102725) q[0];
sx q[0];
rz(-2.0187223) q[0];
sx q[0];
rz(0.3737803) q[0];
rz(-pi) q[1];
x q[1];
rz(2.997424) q[2];
sx q[2];
rz(-1.8509794) q[2];
sx q[2];
rz(-0.35137128) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0337861) q[1];
sx q[1];
rz(-0.68192712) q[1];
sx q[1];
rz(2.4297907) q[1];
rz(-pi) q[2];
rz(-1.9507017) q[3];
sx q[3];
rz(-2.0348747) q[3];
sx q[3];
rz(2.3604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91360056) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733923) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(-1.9447928) q[0];
rz(0.67990047) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(1.4555567) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14917063) q[0];
sx q[0];
rz(-2.143321) q[0];
sx q[0];
rz(-0.19897977) q[0];
x q[1];
rz(-0.20632867) q[2];
sx q[2];
rz(-2.7596139) q[2];
sx q[2];
rz(-2.1260335) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.62555186) q[1];
sx q[1];
rz(-1.5572773) q[1];
sx q[1];
rz(1.0479755) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0807642) q[3];
sx q[3];
rz(-0.41234499) q[3];
sx q[3];
rz(-2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7130647) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(-2.9591566) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(-2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882554) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(1.172539) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81329936) q[0];
sx q[0];
rz(-2.1345703) q[0];
sx q[0];
rz(-1.9378807) q[0];
rz(-1.2433979) q[2];
sx q[2];
rz(-0.84257579) q[2];
sx q[2];
rz(2.3936405) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8923924) q[1];
sx q[1];
rz(-1.3274267) q[1];
sx q[1];
rz(-2.1057486) q[1];
rz(-pi) q[2];
rz(-1.8215239) q[3];
sx q[3];
rz(-1.513895) q[3];
sx q[3];
rz(-1.0292605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0671493) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(0.91119901) q[2];
rz(-2.1905812) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(-2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(-3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(-0.52350837) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9317755) q[0];
sx q[0];
rz(-0.99299252) q[0];
sx q[0];
rz(1.125976) q[0];
rz(-2.8115494) q[2];
sx q[2];
rz(-1.4903307) q[2];
sx q[2];
rz(2.6835494) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7586786) q[1];
sx q[1];
rz(-2.3944693) q[1];
sx q[1];
rz(1.9488571) q[1];
rz(-1.773049) q[3];
sx q[3];
rz(-0.60086717) q[3];
sx q[3];
rz(2.3893389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52544242) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(2.5775487) q[2];
rz(2.8530252) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(-2.585876) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600835) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(-2.2619757) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-0.99194828) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.49682) q[0];
sx q[0];
rz(-0.20475514) q[0];
sx q[0];
rz(0.55069189) q[0];
x q[1];
rz(-1.2679342) q[2];
sx q[2];
rz(-1.5793243) q[2];
sx q[2];
rz(-1.0964583) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76584133) q[1];
sx q[1];
rz(-1.5305133) q[1];
sx q[1];
rz(-0.30893107) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81327849) q[3];
sx q[3];
rz(-1.4762029) q[3];
sx q[3];
rz(2.7300342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1296967) q[2];
sx q[2];
rz(-0.36964881) q[2];
sx q[2];
rz(-2.8707855) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(-2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145684) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(-2.7596966) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.4250925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0483635) q[0];
sx q[0];
rz(-1.3472124) q[0];
sx q[0];
rz(-2.6327052) q[0];
x q[1];
rz(0.09128696) q[2];
sx q[2];
rz(-1.8249776) q[2];
sx q[2];
rz(0.19677481) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6188366) q[1];
sx q[1];
rz(-2.6869876) q[1];
sx q[1];
rz(1.418581) q[1];
rz(-pi) q[2];
rz(-2.1720042) q[3];
sx q[3];
rz(-1.117327) q[3];
sx q[3];
rz(0.78905247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1340593) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-0.40294161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6329704) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-0.41931835) q[0];
rz(1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-2.3197876) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99188995) q[0];
sx q[0];
rz(-2.2737962) q[0];
sx q[0];
rz(-2.4292612) q[0];
rz(-pi) q[1];
rz(0.84341151) q[2];
sx q[2];
rz(-1.1864098) q[2];
sx q[2];
rz(2.9316528) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0710443) q[1];
sx q[1];
rz(-2.3729635) q[1];
sx q[1];
rz(-0.32960906) q[1];
rz(-pi) q[2];
rz(-0.58737289) q[3];
sx q[3];
rz(-1.2578739) q[3];
sx q[3];
rz(1.8250993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(3.0120567) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(-1.6285508) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(-0.82292557) q[0];
rz(-2.8322463) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(-1.3051422) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5236429) q[0];
sx q[0];
rz(-2.4542913) q[0];
sx q[0];
rz(2.6108517) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42758503) q[2];
sx q[2];
rz(-2.3901849) q[2];
sx q[2];
rz(-2.4497355) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0481938) q[1];
sx q[1];
rz(-1.1953925) q[1];
sx q[1];
rz(-0.6742538) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6373789) q[3];
sx q[3];
rz(-0.3399907) q[3];
sx q[3];
rz(0.41223994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(2.819678) q[0];
rz(1.6053258) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(-2.4386491) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84530572) q[0];
sx q[0];
rz(-1.9585113) q[0];
sx q[0];
rz(1.182343) q[0];
rz(-pi) q[1];
x q[1];
rz(0.039673294) q[2];
sx q[2];
rz(-1.0618321) q[2];
sx q[2];
rz(0.85658011) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.34300464) q[1];
sx q[1];
rz(-1.7354021) q[1];
sx q[1];
rz(-2.1144457) q[1];
rz(-pi) q[2];
rz(-0.95548198) q[3];
sx q[3];
rz(-1.2947047) q[3];
sx q[3];
rz(2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2408509) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(1.3396324) q[2];
rz(0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(-1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777578) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(-1.2257858) q[0];
rz(2.2380791) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7077431) q[0];
sx q[0];
rz(-1.9540457) q[0];
sx q[0];
rz(1.7953403) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87601985) q[2];
sx q[2];
rz(-2.6555736) q[2];
sx q[2];
rz(1.8234058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84506449) q[1];
sx q[1];
rz(-1.6390641) q[1];
sx q[1];
rz(-1.8370085) q[1];
rz(-0.43379421) q[3];
sx q[3];
rz(-1.1489831) q[3];
sx q[3];
rz(2.1123561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6580711) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(1.998385) q[2];
rz(0.11463595) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951915) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(0.62190965) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(0.68998228) q[2];
sx q[2];
rz(-0.98914115) q[2];
sx q[2];
rz(3.0422899) q[2];
rz(2.2223496) q[3];
sx q[3];
rz(-1.5083434) q[3];
sx q[3];
rz(-1.2283243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
