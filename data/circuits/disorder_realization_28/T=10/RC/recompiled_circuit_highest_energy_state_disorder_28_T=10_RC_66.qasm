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
rz(0.42143917) q[0];
sx q[0];
rz(5.0285664) q[0];
sx q[0];
rz(9.9528735) q[0];
rz(-0.34673196) q[1];
sx q[1];
rz(4.4693153) q[1];
sx q[1];
rz(9.1413924) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.222966) q[0];
sx q[0];
rz(-0.98969995) q[0];
sx q[0];
rz(1.4544992) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82132116) q[2];
sx q[2];
rz(-1.9275713) q[2];
sx q[2];
rz(-0.58134628) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2436851) q[1];
sx q[1];
rz(-2.7449611) q[1];
sx q[1];
rz(-0.70901386) q[1];
x q[2];
rz(-1.1685405) q[3];
sx q[3];
rz(-1.6682955) q[3];
sx q[3];
rz(1.1643226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9458719) q[2];
sx q[2];
rz(-1.6678145) q[2];
sx q[2];
rz(0.26738581) q[2];
rz(0.17313677) q[3];
sx q[3];
rz(-2.0045529) q[3];
sx q[3];
rz(-0.6828298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46083573) q[0];
sx q[0];
rz(-2.2956678) q[0];
sx q[0];
rz(0.2980921) q[0];
rz(2.8744892) q[1];
sx q[1];
rz(-2.2593081) q[1];
sx q[1];
rz(-0.90046901) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1698477) q[0];
sx q[0];
rz(-1.8119703) q[0];
sx q[0];
rz(1.826188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7212333) q[2];
sx q[2];
rz(-1.7152046) q[2];
sx q[2];
rz(-2.9109746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5067277) q[1];
sx q[1];
rz(-1.7657451) q[1];
sx q[1];
rz(-0.69201236) q[1];
rz(-pi) q[2];
rz(1.9252842) q[3];
sx q[3];
rz(-2.1976314) q[3];
sx q[3];
rz(1.9821061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4233826) q[2];
sx q[2];
rz(-0.98576468) q[2];
sx q[2];
rz(0.59744376) q[2];
rz(-2.1256223) q[3];
sx q[3];
rz(-1.6705325) q[3];
sx q[3];
rz(-1.5731251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.51282561) q[0];
sx q[0];
rz(-0.91575423) q[0];
sx q[0];
rz(0.065091982) q[0];
rz(-1.1029296) q[1];
sx q[1];
rz(-1.8850336) q[1];
sx q[1];
rz(0.38806134) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2723615) q[0];
sx q[0];
rz(-0.80565208) q[0];
sx q[0];
rz(2.8153573) q[0];
x q[1];
rz(1.6460288) q[2];
sx q[2];
rz(-1.2047775) q[2];
sx q[2];
rz(0.96142069) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2414723) q[1];
sx q[1];
rz(-2.6992976) q[1];
sx q[1];
rz(1.9458179) q[1];
rz(-pi) q[2];
rz(-0.89395071) q[3];
sx q[3];
rz(-0.36374796) q[3];
sx q[3];
rz(-1.8574024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1661561) q[2];
sx q[2];
rz(-1.0338975) q[2];
sx q[2];
rz(-0.088509716) q[2];
rz(-2.4837808) q[3];
sx q[3];
rz(-0.43840539) q[3];
sx q[3];
rz(-1.8678166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6863962) q[0];
sx q[0];
rz(-0.96298591) q[0];
sx q[0];
rz(-2.2121867) q[0];
rz(-1.9724253) q[1];
sx q[1];
rz(-2.2678352) q[1];
sx q[1];
rz(-3.015231) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0189782) q[0];
sx q[0];
rz(-0.41176957) q[0];
sx q[0];
rz(-1.3995348) q[0];
x q[1];
rz(-1.3643372) q[2];
sx q[2];
rz(-1.8935618) q[2];
sx q[2];
rz(1.1870804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1388701) q[1];
sx q[1];
rz(-2.5926081) q[1];
sx q[1];
rz(2.5217338) q[1];
rz(1.3170304) q[3];
sx q[3];
rz(-1.9477748) q[3];
sx q[3];
rz(-2.1351817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29752842) q[2];
sx q[2];
rz(-1.0928096) q[2];
sx q[2];
rz(-1.7139942) q[2];
rz(-1.1045688) q[3];
sx q[3];
rz(-1.5966871) q[3];
sx q[3];
rz(0.67232084) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70668689) q[0];
sx q[0];
rz(-0.66449419) q[0];
sx q[0];
rz(-1.8988761) q[0];
rz(-3.0022314) q[1];
sx q[1];
rz(-0.31529537) q[1];
sx q[1];
rz(0.79162663) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46027943) q[0];
sx q[0];
rz(-2.9834207) q[0];
sx q[0];
rz(-1.3307443) q[0];
x q[1];
rz(2.5734719) q[2];
sx q[2];
rz(-1.4213398) q[2];
sx q[2];
rz(-2.7338701) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.086603786) q[1];
sx q[1];
rz(-1.5468925) q[1];
sx q[1];
rz(1.8782645) q[1];
rz(-pi) q[2];
rz(-1.743312) q[3];
sx q[3];
rz(-0.60997395) q[3];
sx q[3];
rz(1.139918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41607729) q[2];
sx q[2];
rz(-1.6107586) q[2];
sx q[2];
rz(-1.3610972) q[2];
rz(2.1040037) q[3];
sx q[3];
rz(-1.6890084) q[3];
sx q[3];
rz(3.1328372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233651) q[0];
sx q[0];
rz(-0.26015493) q[0];
sx q[0];
rz(2.9628229) q[0];
rz(-0.44890064) q[1];
sx q[1];
rz(-1.1652378) q[1];
sx q[1];
rz(-1.9662205) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24488388) q[0];
sx q[0];
rz(-1.5053466) q[0];
sx q[0];
rz(2.3785832) q[0];
rz(-pi) q[1];
x q[1];
rz(1.825338) q[2];
sx q[2];
rz(-2.3582728) q[2];
sx q[2];
rz(-0.50031991) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11131249) q[1];
sx q[1];
rz(-2.2759109) q[1];
sx q[1];
rz(2.6837803) q[1];
rz(-pi) q[2];
rz(-2.6244632) q[3];
sx q[3];
rz(-1.7742566) q[3];
sx q[3];
rz(0.77282016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7281404) q[2];
sx q[2];
rz(-2.7588625) q[2];
sx q[2];
rz(-2.9276796) q[2];
rz(1.4296069) q[3];
sx q[3];
rz(-1.471328) q[3];
sx q[3];
rz(-1.3235929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0837309) q[0];
sx q[0];
rz(-1.2849176) q[0];
sx q[0];
rz(2.4833552) q[0];
rz(-1.4312875) q[1];
sx q[1];
rz(-0.88600102) q[1];
sx q[1];
rz(-1.5586982) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72082544) q[0];
sx q[0];
rz(-2.0090734) q[0];
sx q[0];
rz(1.7780316) q[0];
rz(-pi) q[1];
rz(-1.8020242) q[2];
sx q[2];
rz(-1.8699081) q[2];
sx q[2];
rz(2.0264152) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78608785) q[1];
sx q[1];
rz(-1.2337323) q[1];
sx q[1];
rz(-0.33871594) q[1];
x q[2];
rz(0.57969661) q[3];
sx q[3];
rz(-1.2384185) q[3];
sx q[3];
rz(-1.9975091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7610641) q[2];
sx q[2];
rz(-1.5917799) q[2];
sx q[2];
rz(1.2065678) q[2];
rz(-2.1207464) q[3];
sx q[3];
rz(-2.6267093) q[3];
sx q[3];
rz(0.90768874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1445439) q[0];
sx q[0];
rz(-0.21268614) q[0];
sx q[0];
rz(0.3120684) q[0];
rz(-0.32876217) q[1];
sx q[1];
rz(-1.5212675) q[1];
sx q[1];
rz(-2.388248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8341136) q[0];
sx q[0];
rz(-2.284482) q[0];
sx q[0];
rz(2.3673348) q[0];
rz(-1.5320832) q[2];
sx q[2];
rz(-0.51100327) q[2];
sx q[2];
rz(1.7450361) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6623913) q[1];
sx q[1];
rz(-1.1117142) q[1];
sx q[1];
rz(2.4513112) q[1];
x q[2];
rz(-2.1301792) q[3];
sx q[3];
rz(-2.0205343) q[3];
sx q[3];
rz(1.442329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.83851492) q[2];
sx q[2];
rz(-1.3884156) q[2];
sx q[2];
rz(-2.370749) q[2];
rz(0.54909697) q[3];
sx q[3];
rz(-1.2884459) q[3];
sx q[3];
rz(-0.87636605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2751806) q[0];
sx q[0];
rz(-2.488945) q[0];
sx q[0];
rz(1.8021679) q[0];
rz(-2.0088947) q[1];
sx q[1];
rz(-0.39342543) q[1];
sx q[1];
rz(-0.045230953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6409174) q[0];
sx q[0];
rz(-1.7637327) q[0];
sx q[0];
rz(1.1290324) q[0];
rz(-pi) q[1];
rz(-1.9275749) q[2];
sx q[2];
rz(-1.456481) q[2];
sx q[2];
rz(-3.1014991) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4457827) q[1];
sx q[1];
rz(-0.28430609) q[1];
sx q[1];
rz(0.48281702) q[1];
rz(1.955785) q[3];
sx q[3];
rz(-0.5866881) q[3];
sx q[3];
rz(-1.7314814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.905978) q[2];
sx q[2];
rz(-2.2048042) q[2];
sx q[2];
rz(-2.3255685) q[2];
rz(1.7488545) q[3];
sx q[3];
rz(-2.0413028) q[3];
sx q[3];
rz(-1.6284774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.29869646) q[0];
sx q[0];
rz(-0.5492292) q[0];
sx q[0];
rz(-1.5811051) q[0];
rz(2.9807978) q[1];
sx q[1];
rz(-2.000688) q[1];
sx q[1];
rz(-2.2198832) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.164249) q[0];
sx q[0];
rz(-1.5802339) q[0];
sx q[0];
rz(-1.3816076) q[0];
rz(0.99803136) q[2];
sx q[2];
rz(-1.0596399) q[2];
sx q[2];
rz(-2.9645408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8556577) q[1];
sx q[1];
rz(-1.8536356) q[1];
sx q[1];
rz(-2.9873579) q[1];
rz(-pi) q[2];
rz(-1.7245737) q[3];
sx q[3];
rz(-1.7821081) q[3];
sx q[3];
rz(-0.37195027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.86772743) q[2];
sx q[2];
rz(-0.92550698) q[2];
sx q[2];
rz(-1.072139) q[2];
rz(2.9169361) q[3];
sx q[3];
rz(-1.4916689) q[3];
sx q[3];
rz(2.2255161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12162019) q[0];
sx q[0];
rz(-2.2157123) q[0];
sx q[0];
rz(-2.8257688) q[0];
rz(-3.0674122) q[1];
sx q[1];
rz(-1.0646432) q[1];
sx q[1];
rz(-0.021312996) q[1];
rz(2.3580768) q[2];
sx q[2];
rz(-1.5512237) q[2];
sx q[2];
rz(1.5882391) q[2];
rz(0.40423468) q[3];
sx q[3];
rz(-0.80818292) q[3];
sx q[3];
rz(-0.72738739) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
