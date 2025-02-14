OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.030169686) q[0];
sx q[0];
rz(4.5017894) q[0];
sx q[0];
rz(5.9202249) q[0];
rz(-1.1765923) q[1];
sx q[1];
rz(-1.7791553) q[1];
sx q[1];
rz(-1.4443614) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.022665) q[0];
sx q[0];
rz(-1.3937409) q[0];
sx q[0];
rz(1.4757968) q[0];
rz(2.7411103) q[2];
sx q[2];
rz(-2.8304511) q[2];
sx q[2];
rz(0.36382407) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0054659) q[1];
sx q[1];
rz(-2.695752) q[1];
sx q[1];
rz(1.3492829) q[1];
rz(-pi) q[2];
rz(0.96935848) q[3];
sx q[3];
rz(-1.6043318) q[3];
sx q[3];
rz(2.7905317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.968367) q[2];
sx q[2];
rz(-2.3543365) q[2];
sx q[2];
rz(3.0464029) q[2];
rz(1.7732874) q[3];
sx q[3];
rz(-0.94822001) q[3];
sx q[3];
rz(-0.29909721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.478941) q[0];
sx q[0];
rz(-0.79942411) q[0];
sx q[0];
rz(-0.69315243) q[0];
rz(2.8043546) q[1];
sx q[1];
rz(-1.3062898) q[1];
sx q[1];
rz(0.79280028) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62294856) q[0];
sx q[0];
rz(-2.3272133) q[0];
sx q[0];
rz(1.2252392) q[0];
x q[1];
rz(-1.1801925) q[2];
sx q[2];
rz(-0.55093575) q[2];
sx q[2];
rz(-1.40755) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.63607994) q[1];
sx q[1];
rz(-0.82712251) q[1];
sx q[1];
rz(-1.176314) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4850785) q[3];
sx q[3];
rz(-0.88036075) q[3];
sx q[3];
rz(-0.2218962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1482131) q[2];
sx q[2];
rz(-0.17017636) q[2];
sx q[2];
rz(-1.4519838) q[2];
rz(0.65886894) q[3];
sx q[3];
rz(-1.6784607) q[3];
sx q[3];
rz(-0.39670593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0186998) q[0];
sx q[0];
rz(-2.1227699) q[0];
sx q[0];
rz(-0.10957154) q[0];
rz(-0.84596363) q[1];
sx q[1];
rz(-2.1956317) q[1];
sx q[1];
rz(-0.74839655) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3050675) q[0];
sx q[0];
rz(-2.5782707) q[0];
sx q[0];
rz(1.402424) q[0];
x q[1];
rz(1.0399299) q[2];
sx q[2];
rz(-1.0541087) q[2];
sx q[2];
rz(2.6831107) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9108991) q[1];
sx q[1];
rz(-1.5077796) q[1];
sx q[1];
rz(0.0067733532) q[1];
rz(-2.142228) q[3];
sx q[3];
rz(-1.7212221) q[3];
sx q[3];
rz(2.5932464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0261859) q[2];
sx q[2];
rz(-1.4226961) q[2];
sx q[2];
rz(-2.4288948) q[2];
rz(-1.4767492) q[3];
sx q[3];
rz(-1.289117) q[3];
sx q[3];
rz(-3.0570928) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1348212) q[0];
sx q[0];
rz(-1.1399784) q[0];
sx q[0];
rz(0.086932927) q[0];
rz(2.8839819) q[1];
sx q[1];
rz(-2.0517709) q[1];
sx q[1];
rz(1.2923406) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1202165) q[0];
sx q[0];
rz(-0.82737404) q[0];
sx q[0];
rz(2.8143456) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7316225) q[2];
sx q[2];
rz(-1.8984744) q[2];
sx q[2];
rz(-0.92152061) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4166761) q[1];
sx q[1];
rz(-1.9179763) q[1];
sx q[1];
rz(-2.2031839) q[1];
rz(-0.38141187) q[3];
sx q[3];
rz(-0.3837882) q[3];
sx q[3];
rz(0.82510766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2963691) q[2];
sx q[2];
rz(-1.1901647) q[2];
sx q[2];
rz(-2.2947218) q[2];
rz(1.752468) q[3];
sx q[3];
rz(-1.4877157) q[3];
sx q[3];
rz(-0.59051591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.1883989) q[0];
sx q[0];
rz(-1.1253091) q[0];
sx q[0];
rz(0.97287792) q[0];
rz(2.1994622) q[1];
sx q[1];
rz(-2.3140488) q[1];
sx q[1];
rz(-3.1408302) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98648345) q[0];
sx q[0];
rz(-2.0707978) q[0];
sx q[0];
rz(-0.9802823) q[0];
rz(-pi) q[1];
rz(-2.5070287) q[2];
sx q[2];
rz(-1.8406313) q[2];
sx q[2];
rz(0.25411221) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3734596) q[1];
sx q[1];
rz(-2.4985857) q[1];
sx q[1];
rz(-2.0483584) q[1];
rz(-pi) q[2];
rz(2.5464779) q[3];
sx q[3];
rz(-1.6356118) q[3];
sx q[3];
rz(0.3583554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3804649) q[2];
sx q[2];
rz(-1.7385812) q[2];
sx q[2];
rz(0.38499704) q[2];
rz(0.72832406) q[3];
sx q[3];
rz(-2.5326122) q[3];
sx q[3];
rz(0.42496067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7321135) q[0];
sx q[0];
rz(-0.80687579) q[0];
sx q[0];
rz(2.683486) q[0];
rz(-1.7181646) q[1];
sx q[1];
rz(-0.79045311) q[1];
sx q[1];
rz(3.0480393) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76339611) q[0];
sx q[0];
rz(-2.3468821) q[0];
sx q[0];
rz(2.2346157) q[0];
rz(-pi) q[1];
rz(-1.3504998) q[2];
sx q[2];
rz(-1.4080021) q[2];
sx q[2];
rz(0.3045813) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72274745) q[1];
sx q[1];
rz(-2.0993352) q[1];
sx q[1];
rz(-1.3206119) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5597975) q[3];
sx q[3];
rz(-1.5526532) q[3];
sx q[3];
rz(1.1813576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2576922) q[2];
sx q[2];
rz(-0.42282405) q[2];
sx q[2];
rz(0.60919961) q[2];
rz(-0.60112634) q[3];
sx q[3];
rz(-1.6060035) q[3];
sx q[3];
rz(-2.6570184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5179829) q[0];
sx q[0];
rz(-1.4731151) q[0];
sx q[0];
rz(-1.8286888) q[0];
rz(3.0598705) q[1];
sx q[1];
rz(-1.3879958) q[1];
sx q[1];
rz(1.0645701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.222932) q[0];
sx q[0];
rz(-2.1188508) q[0];
sx q[0];
rz(0.81848829) q[0];
x q[1];
rz(2.1871952) q[2];
sx q[2];
rz(-1.5845044) q[2];
sx q[2];
rz(-0.15132667) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7408167) q[1];
sx q[1];
rz(-0.99641748) q[1];
sx q[1];
rz(3.0568491) q[1];
rz(0.96780583) q[3];
sx q[3];
rz(-1.2295614) q[3];
sx q[3];
rz(-1.5016709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.99819034) q[2];
sx q[2];
rz(-0.4021796) q[2];
sx q[2];
rz(1.9833938) q[2];
rz(-2.4573333) q[3];
sx q[3];
rz(-1.8161769) q[3];
sx q[3];
rz(-1.6392596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.174468) q[0];
sx q[0];
rz(-3.0376349) q[0];
sx q[0];
rz(2.6167468) q[0];
rz(-1.1809433) q[1];
sx q[1];
rz(-0.8371822) q[1];
sx q[1];
rz(-0.7276853) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1215661) q[0];
sx q[0];
rz(-1.354573) q[0];
sx q[0];
rz(-0.89631594) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5843975) q[2];
sx q[2];
rz(-2.3378125) q[2];
sx q[2];
rz(-0.4280099) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7512253) q[1];
sx q[1];
rz(-1.3074682) q[1];
sx q[1];
rz(0.45231426) q[1];
rz(-0.10730174) q[3];
sx q[3];
rz(-0.91917097) q[3];
sx q[3];
rz(-2.1880423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22964792) q[2];
sx q[2];
rz(-1.0934528) q[2];
sx q[2];
rz(-1.3463119) q[2];
rz(-0.66156578) q[3];
sx q[3];
rz(-1.9965568) q[3];
sx q[3];
rz(-0.0865817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10053703) q[0];
sx q[0];
rz(-2.4424398) q[0];
sx q[0];
rz(2.2030785) q[0];
rz(-1.8427294) q[1];
sx q[1];
rz(-0.88911903) q[1];
sx q[1];
rz(3.0079957) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4105281) q[0];
sx q[0];
rz(-2.799187) q[0];
sx q[0];
rz(0.28241856) q[0];
rz(-pi) q[1];
rz(-0.44914896) q[2];
sx q[2];
rz(-2.5173016) q[2];
sx q[2];
rz(2.1438832) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6386968) q[1];
sx q[1];
rz(-1.2112482) q[1];
sx q[1];
rz(1.4877968) q[1];
rz(-pi) q[2];
rz(1.9383921) q[3];
sx q[3];
rz(-2.1798686) q[3];
sx q[3];
rz(-1.6742791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9962697) q[2];
sx q[2];
rz(-2.1349553) q[2];
sx q[2];
rz(-0.88627306) q[2];
rz(-2.8730872) q[3];
sx q[3];
rz(-1.0849181) q[3];
sx q[3];
rz(1.6296384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71235424) q[0];
sx q[0];
rz(-2.7463284) q[0];
sx q[0];
rz(-2.9048753) q[0];
rz(-2.9780544) q[1];
sx q[1];
rz(-1.6103368) q[1];
sx q[1];
rz(2.9530361) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084219639) q[0];
sx q[0];
rz(-1.9002267) q[0];
sx q[0];
rz(-2.2994142) q[0];
x q[1];
rz(-2.677146) q[2];
sx q[2];
rz(-2.6147343) q[2];
sx q[2];
rz(1.7773445) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.982354) q[1];
sx q[1];
rz(-0.51430632) q[1];
sx q[1];
rz(1.6883172) q[1];
rz(-pi) q[2];
rz(-0.27310009) q[3];
sx q[3];
rz(-2.6870603) q[3];
sx q[3];
rz(2.4668193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4208372) q[2];
sx q[2];
rz(-1.5034224) q[2];
sx q[2];
rz(-2.3756964) q[2];
rz(0.01865538) q[3];
sx q[3];
rz(-1.3466287) q[3];
sx q[3];
rz(-2.6115821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50273773) q[0];
sx q[0];
rz(-1.9866332) q[0];
sx q[0];
rz(-1.6992983) q[0];
rz(-0.28932183) q[1];
sx q[1];
rz(-1.0565636) q[1];
sx q[1];
rz(0.4531959) q[1];
rz(-1.4353903) q[2];
sx q[2];
rz(-1.0359905) q[2];
sx q[2];
rz(-0.059992487) q[2];
rz(-2.0404242) q[3];
sx q[3];
rz(-1.6328937) q[3];
sx q[3];
rz(1.9619305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
