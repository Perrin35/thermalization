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
rz(-1.363938) q[0];
sx q[0];
rz(3.5600297) q[0];
sx q[0];
rz(10.726396) q[0];
rz(-2.9786181) q[1];
sx q[1];
rz(-1.6956704) q[1];
sx q[1];
rz(3.1248098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7201286) q[0];
sx q[0];
rz(-1.4961494) q[0];
sx q[0];
rz(-1.9326841) q[0];
rz(1.391259) q[2];
sx q[2];
rz(-1.554536) q[2];
sx q[2];
rz(0.017841466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77474817) q[1];
sx q[1];
rz(-1.5758744) q[1];
sx q[1];
rz(-3.1366411) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34256012) q[3];
sx q[3];
rz(-2.9800219) q[3];
sx q[3];
rz(-1.2662966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5590543) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(-1.5793229) q[2];
rz(-2.9301379) q[3];
sx q[3];
rz(-0.00051694218) q[3];
sx q[3];
rz(-2.9805984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6120537) q[0];
sx q[0];
rz(-2.8654629) q[0];
sx q[0];
rz(1.8147234) q[0];
rz(-2.5621085) q[1];
sx q[1];
rz(-0.0038298413) q[1];
sx q[1];
rz(-2.5025867) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3746637) q[0];
sx q[0];
rz(-2.2013454) q[0];
sx q[0];
rz(-0.85682822) q[0];
x q[1];
rz(1.5875568) q[2];
sx q[2];
rz(-1.4498561) q[2];
sx q[2];
rz(-3.1230833) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.17520576) q[1];
sx q[1];
rz(-1.5537973) q[1];
sx q[1];
rz(1.5592038) q[1];
rz(0.78044807) q[3];
sx q[3];
rz(-1.6326346) q[3];
sx q[3];
rz(-1.058418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4352033) q[2];
sx q[2];
rz(-3.0053164) q[2];
sx q[2];
rz(-1.603568) q[2];
rz(1.5609353) q[3];
sx q[3];
rz(-0.014336421) q[3];
sx q[3];
rz(-3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8921709) q[0];
sx q[0];
rz(-0.51512655) q[0];
sx q[0];
rz(0.38145915) q[0];
rz(-0.70746607) q[1];
sx q[1];
rz(-0.019376945) q[1];
sx q[1];
rz(-2.0170508) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19605787) q[0];
sx q[0];
rz(-1.811341) q[0];
sx q[0];
rz(-1.3105767) q[0];
x q[1];
rz(1.4558305) q[2];
sx q[2];
rz(-1.5966151) q[2];
sx q[2];
rz(1.4173649) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.687011) q[1];
sx q[1];
rz(-0.064116009) q[1];
sx q[1];
rz(-2.9518576) q[1];
rz(-pi) q[2];
rz(-2.5865058) q[3];
sx q[3];
rz(-1.1218117) q[3];
sx q[3];
rz(-0.6787231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6996998) q[2];
sx q[2];
rz(-3.1293588) q[2];
sx q[2];
rz(0.049467889) q[2];
rz(0.60702819) q[3];
sx q[3];
rz(-3.1403465) q[3];
sx q[3];
rz(1.9342669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98168755) q[0];
sx q[0];
rz(-0.1683546) q[0];
sx q[0];
rz(0.017177563) q[0];
rz(2.848564) q[1];
sx q[1];
rz(-2.3511062) q[1];
sx q[1];
rz(1.5471829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80205432) q[0];
sx q[0];
rz(-1.0419163) q[0];
sx q[0];
rz(0.89978973) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54173476) q[2];
sx q[2];
rz(-1.3770665) q[2];
sx q[2];
rz(1.8371179) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77212376) q[1];
sx q[1];
rz(-1.511206) q[1];
sx q[1];
rz(0.12424424) q[1];
rz(-pi) q[2];
rz(0.64597102) q[3];
sx q[3];
rz(-1.5785909) q[3];
sx q[3];
rz(0.92168671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.25345099) q[2];
sx q[2];
rz(-2.6710822) q[2];
sx q[2];
rz(2.4593501) q[2];
rz(3.0864129) q[3];
sx q[3];
rz(-3.1339055) q[3];
sx q[3];
rz(1.8365708) q[3];
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
rz(0.76160112) q[0];
sx q[0];
rz(-0.43552265) q[0];
sx q[0];
rz(0.4250266) q[0];
rz(-1.5399326) q[1];
sx q[1];
rz(-0.48307499) q[1];
sx q[1];
rz(2.3262598) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20629932) q[0];
sx q[0];
rz(-1.5599129) q[0];
sx q[0];
rz(0.061200415) q[0];
rz(-pi) q[1];
rz(1.5940395) q[2];
sx q[2];
rz(-1.574062) q[2];
sx q[2];
rz(-2.4583985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9141691) q[1];
sx q[1];
rz(-2.9939751) q[1];
sx q[1];
rz(-2.5116337) q[1];
x q[2];
rz(-0.32834239) q[3];
sx q[3];
rz(-2.7797065) q[3];
sx q[3];
rz(-0.34235024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0067979) q[2];
sx q[2];
rz(-0.012601348) q[2];
sx q[2];
rz(-1.4726144) q[2];
rz(-0.93004477) q[3];
sx q[3];
rz(-3.127122) q[3];
sx q[3];
rz(-2.3013733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.8188266) q[0];
sx q[0];
rz(-3.1184986) q[0];
sx q[0];
rz(-1.4267138) q[0];
rz(0.73000437) q[1];
sx q[1];
rz(-0.58861029) q[1];
sx q[1];
rz(-2.0402562) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5663213) q[0];
sx q[0];
rz(-2.3744517) q[0];
sx q[0];
rz(-1.1772035) q[0];
rz(-pi) q[1];
rz(2.9726102) q[2];
sx q[2];
rz(-1.3436396) q[2];
sx q[2];
rz(-2.3262466) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37820617) q[1];
sx q[1];
rz(-0.15608938) q[1];
sx q[1];
rz(-2.2264477) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1093842) q[3];
sx q[3];
rz(-1.4463436) q[3];
sx q[3];
rz(-2.2197753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71658984) q[2];
sx q[2];
rz(-3.0807107) q[2];
sx q[2];
rz(1.8343743) q[2];
rz(0.37846765) q[3];
sx q[3];
rz(-0.022947939) q[3];
sx q[3];
rz(-0.65346658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8398447) q[0];
sx q[0];
rz(-1.8740338) q[0];
sx q[0];
rz(-0.92754716) q[0];
rz(-1.357366) q[1];
sx q[1];
rz(-0.83186847) q[1];
sx q[1];
rz(1.5313139) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46857014) q[0];
sx q[0];
rz(-0.044402145) q[0];
sx q[0];
rz(1.0442249) q[0];
x q[1];
rz(2.7448517) q[2];
sx q[2];
rz(-1.7974241) q[2];
sx q[2];
rz(2.6643348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5802637) q[1];
sx q[1];
rz(-3.0124843) q[1];
sx q[1];
rz(-1.5498954) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31503265) q[3];
sx q[3];
rz(-1.0803534) q[3];
sx q[3];
rz(-0.88446188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.36074582) q[2];
sx q[2];
rz(-3.1372034) q[2];
sx q[2];
rz(1.2537664) q[2];
rz(0.69416657) q[3];
sx q[3];
rz(-0.73918754) q[3];
sx q[3];
rz(-2.9085801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6041782) q[0];
sx q[0];
rz(-2.1359213) q[0];
sx q[0];
rz(-2.1042714) q[0];
rz(1.6088156) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(1.467009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64325414) q[0];
sx q[0];
rz(-0.20659978) q[0];
sx q[0];
rz(-1.2363966) q[0];
x q[1];
rz(1.59378) q[2];
sx q[2];
rz(-1.2927552) q[2];
sx q[2];
rz(-0.26593033) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8407708) q[1];
sx q[1];
rz(-1.5696708) q[1];
sx q[1];
rz(-0.00043934396) q[1];
x q[2];
rz(0.45536228) q[3];
sx q[3];
rz(-1.6557367) q[3];
sx q[3];
rz(-2.7698539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1303225) q[2];
sx q[2];
rz(-2.932817) q[2];
sx q[2];
rz(0.043896349) q[2];
rz(-0.52055001) q[3];
sx q[3];
rz(-0.0046516727) q[3];
sx q[3];
rz(1.6827778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540102) q[0];
sx q[0];
rz(-3.1391322) q[0];
sx q[0];
rz(-1.822923) q[0];
rz(-1.7240546) q[1];
sx q[1];
rz(-2.8520165) q[1];
sx q[1];
rz(1.5971378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.330724) q[0];
sx q[0];
rz(-1.6858584) q[0];
sx q[0];
rz(-0.48145357) q[0];
rz(2.804432) q[2];
sx q[2];
rz(-2.4558407) q[2];
sx q[2];
rz(0.27039385) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7629976) q[1];
sx q[1];
rz(-1.7471658) q[1];
sx q[1];
rz(-1.2552983) q[1];
x q[2];
rz(0.68759509) q[3];
sx q[3];
rz(-1.9904226) q[3];
sx q[3];
rz(0.18292038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.80457193) q[2];
sx q[2];
rz(-1.2995517) q[2];
sx q[2];
rz(2.9447832) q[2];
rz(1.9491516) q[3];
sx q[3];
rz(-0.20659031) q[3];
sx q[3];
rz(2.9406252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30109677) q[0];
sx q[0];
rz(-1.7844642) q[0];
sx q[0];
rz(-1.949973) q[0];
rz(-1.5246897) q[1];
sx q[1];
rz(-2.4947417) q[1];
sx q[1];
rz(1.5764538) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7708262) q[0];
sx q[0];
rz(-1.2773371) q[0];
sx q[0];
rz(1.7459041) q[0];
rz(2.9124969) q[2];
sx q[2];
rz(-2.014403) q[2];
sx q[2];
rz(-2.6227621) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8305739) q[1];
sx q[1];
rz(-1.5717447) q[1];
sx q[1];
rz(-1.5702973) q[1];
rz(-3.0281046) q[3];
sx q[3];
rz(-1.7120512) q[3];
sx q[3];
rz(0.065441386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18371753) q[2];
sx q[2];
rz(-0.58791939) q[2];
sx q[2];
rz(-1.6831762) q[2];
rz(3.1111187) q[3];
sx q[3];
rz(-3.1320429) q[3];
sx q[3];
rz(2.9392346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644345) q[0];
sx q[0];
rz(-1.3344593) q[0];
sx q[0];
rz(1.6819171) q[0];
rz(-1.5741813) q[1];
sx q[1];
rz(-1.3290783) q[1];
sx q[1];
rz(-3.0507416) q[1];
rz(1.5051928) q[2];
sx q[2];
rz(-3.0658683) q[2];
sx q[2];
rz(0.22584596) q[2];
rz(0.50441691) q[3];
sx q[3];
rz(-0.88263369) q[3];
sx q[3];
rz(0.42019444) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
