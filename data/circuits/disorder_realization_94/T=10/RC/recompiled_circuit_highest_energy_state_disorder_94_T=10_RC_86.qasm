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
rz(1.2430159) q[0];
sx q[0];
rz(-1.1216811) q[0];
sx q[0];
rz(-2.6774874) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(2.2321489) q[1];
sx q[1];
rz(7.5997054) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038832207) q[0];
sx q[0];
rz(-1.505672) q[0];
sx q[0];
rz(-1.6054356) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0010536) q[2];
sx q[2];
rz(-2.7978483) q[2];
sx q[2];
rz(-0.02862169) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.99440982) q[1];
sx q[1];
rz(-1.8412245) q[1];
sx q[1];
rz(-0.61073357) q[1];
x q[2];
rz(-1.8620231) q[3];
sx q[3];
rz(-1.816245) q[3];
sx q[3];
rz(-3.0994741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2509649) q[2];
sx q[2];
rz(-0.35718063) q[2];
sx q[2];
rz(-0.29699057) q[2];
rz(-2.4885528) q[3];
sx q[3];
rz(-1.5584471) q[3];
sx q[3];
rz(1.0562586) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5121269) q[0];
sx q[0];
rz(-2.1194206) q[0];
sx q[0];
rz(-2.8566991) q[0];
rz(-3.0902872) q[1];
sx q[1];
rz(-2.3785794) q[1];
sx q[1];
rz(0.6368534) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89556614) q[0];
sx q[0];
rz(-2.0478529) q[0];
sx q[0];
rz(-2.0791847) q[0];
rz(-pi) q[1];
rz(-2.7899488) q[2];
sx q[2];
rz(-1.5696916) q[2];
sx q[2];
rz(2.0948834) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6057558) q[1];
sx q[1];
rz(-2.9510289) q[1];
sx q[1];
rz(-0.41175731) q[1];
x q[2];
rz(-1.4376182) q[3];
sx q[3];
rz(-2.2723327) q[3];
sx q[3];
rz(0.12058903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9597943) q[2];
sx q[2];
rz(-2.2836298) q[2];
sx q[2];
rz(-0.96724969) q[2];
rz(-2.8773384) q[3];
sx q[3];
rz(-1.5033009) q[3];
sx q[3];
rz(-1.9056457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018175) q[0];
sx q[0];
rz(-1.5793261) q[0];
sx q[0];
rz(2.0215969) q[0];
rz(-2.6657875) q[1];
sx q[1];
rz(-1.0775403) q[1];
sx q[1];
rz(0.30398223) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1893727) q[0];
sx q[0];
rz(-1.598453) q[0];
sx q[0];
rz(-2.1364487) q[0];
x q[1];
rz(-0.39865785) q[2];
sx q[2];
rz(-0.2572607) q[2];
sx q[2];
rz(-0.21695732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3545565) q[1];
sx q[1];
rz(-0.72595412) q[1];
sx q[1];
rz(1.4158842) q[1];
rz(1.3720361) q[3];
sx q[3];
rz(-0.89722465) q[3];
sx q[3];
rz(-2.4186132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1659282) q[2];
sx q[2];
rz(-2.1037481) q[2];
sx q[2];
rz(0.69592875) q[2];
rz(-1.1337229) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(-0.55293647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7021779) q[0];
sx q[0];
rz(-2.1372097) q[0];
sx q[0];
rz(-2.0942005) q[0];
rz(-1.9678496) q[1];
sx q[1];
rz(-0.67432299) q[1];
sx q[1];
rz(-1.5210927) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0516617) q[0];
sx q[0];
rz(-1.0703846) q[0];
sx q[0];
rz(-1.4335645) q[0];
rz(-pi) q[1];
rz(-1.9416472) q[2];
sx q[2];
rz(-0.74127889) q[2];
sx q[2];
rz(2.0855057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7340659) q[1];
sx q[1];
rz(-2.1868949) q[1];
sx q[1];
rz(2.1156838) q[1];
rz(-pi) q[2];
rz(0.32471809) q[3];
sx q[3];
rz(-1.3007264) q[3];
sx q[3];
rz(-0.76061676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.095470458) q[2];
sx q[2];
rz(-2.0746524) q[2];
sx q[2];
rz(0.53786892) q[2];
rz(-1.6543903) q[3];
sx q[3];
rz(-2.7345149) q[3];
sx q[3];
rz(-0.68638221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122413) q[0];
sx q[0];
rz(-0.23433267) q[0];
sx q[0];
rz(0.60337639) q[0];
rz(0.76250184) q[1];
sx q[1];
rz(-1.9489894) q[1];
sx q[1];
rz(-0.71294436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099261053) q[0];
sx q[0];
rz(-1.6368027) q[0];
sx q[0];
rz(1.1724654) q[0];
rz(1.845417) q[2];
sx q[2];
rz(-2.6634187) q[2];
sx q[2];
rz(-2.630065) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11345574) q[1];
sx q[1];
rz(-1.077768) q[1];
sx q[1];
rz(0.22844577) q[1];
x q[2];
rz(-2.5239046) q[3];
sx q[3];
rz(-2.0632944) q[3];
sx q[3];
rz(0.69128752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9731628) q[2];
sx q[2];
rz(-2.2242686) q[2];
sx q[2];
rz(-1.470835) q[2];
rz(-0.51703185) q[3];
sx q[3];
rz(-2.0972589) q[3];
sx q[3];
rz(1.2928591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9789199) q[0];
sx q[0];
rz(-0.51386583) q[0];
sx q[0];
rz(2.5901219) q[0];
rz(-0.035471352) q[1];
sx q[1];
rz(-2.008581) q[1];
sx q[1];
rz(1.6112526) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12172518) q[0];
sx q[0];
rz(-2.1824153) q[0];
sx q[0];
rz(-2.3082971) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6257203) q[2];
sx q[2];
rz(-0.74797219) q[2];
sx q[2];
rz(2.2750281) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1332273) q[1];
sx q[1];
rz(-0.83725819) q[1];
sx q[1];
rz(2.0627229) q[1];
rz(1.6812165) q[3];
sx q[3];
rz(-2.4856119) q[3];
sx q[3];
rz(0.15925281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0643206) q[2];
sx q[2];
rz(-2.6214226) q[2];
sx q[2];
rz(-1.2426097) q[2];
rz(-0.51314917) q[3];
sx q[3];
rz(-2.7276701) q[3];
sx q[3];
rz(1.7665524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4921017) q[0];
sx q[0];
rz(-0.92925564) q[0];
sx q[0];
rz(-3.0249366) q[0];
rz(-0.77313441) q[1];
sx q[1];
rz(-1.9832858) q[1];
sx q[1];
rz(-1.6366417) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9336108) q[0];
sx q[0];
rz(-0.79219063) q[0];
sx q[0];
rz(1.9879782) q[0];
rz(-pi) q[1];
x q[1];
rz(0.059594056) q[2];
sx q[2];
rz(-1.9589387) q[2];
sx q[2];
rz(1.3898894) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6500068) q[1];
sx q[1];
rz(-1.2131696) q[1];
sx q[1];
rz(-0.8365583) q[1];
rz(1.8582452) q[3];
sx q[3];
rz(-0.93837591) q[3];
sx q[3];
rz(2.8638864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92631212) q[2];
sx q[2];
rz(-2.8947688) q[2];
sx q[2];
rz(0.8872633) q[2];
rz(-2.4617646) q[3];
sx q[3];
rz(-0.66893783) q[3];
sx q[3];
rz(-2.5086856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.55650869) q[0];
sx q[0];
rz(-1.8483138) q[0];
sx q[0];
rz(-0.65628091) q[0];
rz(0.20690021) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(1.2178749) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.098893) q[0];
sx q[0];
rz(-1.6535361) q[0];
sx q[0];
rz(-1.6059884) q[0];
rz(-0.98274173) q[2];
sx q[2];
rz(-1.2763192) q[2];
sx q[2];
rz(2.5786932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6177442) q[1];
sx q[1];
rz(-0.19568014) q[1];
sx q[1];
rz(-2.5168672) q[1];
rz(-pi) q[2];
rz(-0.37533203) q[3];
sx q[3];
rz(-2.4960244) q[3];
sx q[3];
rz(2.4433492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3394341) q[2];
sx q[2];
rz(-1.755244) q[2];
sx q[2];
rz(-0.68186861) q[2];
rz(-1.1784461) q[3];
sx q[3];
rz(-0.75740564) q[3];
sx q[3];
rz(-1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7449529) q[0];
sx q[0];
rz(-1.6093901) q[0];
sx q[0];
rz(2.9314281) q[0];
rz(2.919803) q[1];
sx q[1];
rz(-2.2237325) q[1];
sx q[1];
rz(0.04235696) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51084679) q[0];
sx q[0];
rz(-1.5928942) q[0];
sx q[0];
rz(0.18760292) q[0];
rz(-pi) q[1];
x q[1];
rz(2.508411) q[2];
sx q[2];
rz(-1.0594308) q[2];
sx q[2];
rz(-0.88054576) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90396229) q[1];
sx q[1];
rz(-1.018486) q[1];
sx q[1];
rz(-1.2321074) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1091408) q[3];
sx q[3];
rz(-2.0570786) q[3];
sx q[3];
rz(-1.2580296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8611341) q[2];
sx q[2];
rz(-2.0529842) q[2];
sx q[2];
rz(2.4750278) q[2];
rz(-0.76513964) q[3];
sx q[3];
rz(-0.20604006) q[3];
sx q[3];
rz(1.7407181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6078981) q[0];
sx q[0];
rz(-0.38132897) q[0];
sx q[0];
rz(-0.68699849) q[0];
rz(-1.8580565) q[1];
sx q[1];
rz(-1.9891519) q[1];
sx q[1];
rz(2.5680465) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6476558) q[0];
sx q[0];
rz(-1.0595317) q[0];
sx q[0];
rz(0.6880811) q[0];
x q[1];
rz(-2.0230471) q[2];
sx q[2];
rz(-0.49511038) q[2];
sx q[2];
rz(0.64436382) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0189218) q[1];
sx q[1];
rz(-1.9098567) q[1];
sx q[1];
rz(-1.0215205) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3775842) q[3];
sx q[3];
rz(-0.90029085) q[3];
sx q[3];
rz(1.3585728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0428697) q[2];
sx q[2];
rz(-1.9645773) q[2];
sx q[2];
rz(3.0549808) q[2];
rz(-0.11836554) q[3];
sx q[3];
rz(-1.5990853) q[3];
sx q[3];
rz(3.0983483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8260228) q[0];
sx q[0];
rz(-1.0181027) q[0];
sx q[0];
rz(-2.3636567) q[0];
rz(2.8553873) q[1];
sx q[1];
rz(-2.310391) q[1];
sx q[1];
rz(-1.8117767) q[1];
rz(-0.61997531) q[2];
sx q[2];
rz(-2.740553) q[2];
sx q[2];
rz(0.16410826) q[2];
rz(-0.44888857) q[3];
sx q[3];
rz(-1.3009303) q[3];
sx q[3];
rz(-3.124685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
