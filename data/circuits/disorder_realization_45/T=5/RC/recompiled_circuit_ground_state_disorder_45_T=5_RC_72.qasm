OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0319808) q[0];
sx q[0];
rz(4.1469753) q[0];
sx q[0];
rz(8.4561705) q[0];
rz(-2.3321416) q[1];
sx q[1];
rz(4.6757698) q[1];
sx q[1];
rz(9.4546131) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3650509) q[0];
sx q[0];
rz(-1.8030226) q[0];
sx q[0];
rz(-2.5864629) q[0];
x q[1];
rz(-3.0709463) q[2];
sx q[2];
rz(-2.321817) q[2];
sx q[2];
rz(1.1653314) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1263784) q[1];
sx q[1];
rz(-2.3087037) q[1];
sx q[1];
rz(-0.53453858) q[1];
rz(-pi) q[2];
rz(2.5621595) q[3];
sx q[3];
rz(-0.16210769) q[3];
sx q[3];
rz(2.1539089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6385103) q[2];
sx q[2];
rz(-2.5677887) q[2];
sx q[2];
rz(2.8199675) q[2];
rz(0.62421787) q[3];
sx q[3];
rz(-1.7523242) q[3];
sx q[3];
rz(-1.6577087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6931848) q[0];
sx q[0];
rz(-0.25633651) q[0];
sx q[0];
rz(-2.2053027) q[0];
rz(1.2721277) q[1];
sx q[1];
rz(-1.5446168) q[1];
sx q[1];
rz(1.3880656) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4988588) q[0];
sx q[0];
rz(-1.513645) q[0];
sx q[0];
rz(1.7903891) q[0];
rz(2.3642367) q[2];
sx q[2];
rz(-1.9665951) q[2];
sx q[2];
rz(1.1071902) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1202045) q[1];
sx q[1];
rz(-2.1600381) q[1];
sx q[1];
rz(2.2256149) q[1];
x q[2];
rz(-0.8901435) q[3];
sx q[3];
rz(-2.806744) q[3];
sx q[3];
rz(-1.0687506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0973413) q[2];
sx q[2];
rz(-2.1772431) q[2];
sx q[2];
rz(-0.67330366) q[2];
rz(1.3168969) q[3];
sx q[3];
rz(-1.8484867) q[3];
sx q[3];
rz(0.84002686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3299385) q[0];
sx q[0];
rz(-1.543777) q[0];
sx q[0];
rz(3.0918308) q[0];
rz(0.20637575) q[1];
sx q[1];
rz(-1.3866837) q[1];
sx q[1];
rz(1.4899563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1288717) q[0];
sx q[0];
rz(-1.4507334) q[0];
sx q[0];
rz(-0.61693667) q[0];
x q[1];
rz(-1.1300722) q[2];
sx q[2];
rz(-1.0606546) q[2];
sx q[2];
rz(0.59512344) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1957383) q[1];
sx q[1];
rz(-2.7496413) q[1];
sx q[1];
rz(0.82328226) q[1];
rz(-2.0288896) q[3];
sx q[3];
rz(-1.8761523) q[3];
sx q[3];
rz(-0.28766838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9516307) q[2];
sx q[2];
rz(-0.54533521) q[2];
sx q[2];
rz(0.0541617) q[2];
rz(1.1545898) q[3];
sx q[3];
rz(-1.8535987) q[3];
sx q[3];
rz(2.0258928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7390249) q[0];
sx q[0];
rz(-1.716528) q[0];
sx q[0];
rz(2.2836852) q[0];
rz(0.67445451) q[1];
sx q[1];
rz(-2.1916788) q[1];
sx q[1];
rz(-1.8106921) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.389849) q[0];
sx q[0];
rz(-1.8065592) q[0];
sx q[0];
rz(-1.7111196) q[0];
rz(-pi) q[1];
rz(-1.0687549) q[2];
sx q[2];
rz(-0.34345657) q[2];
sx q[2];
rz(-0.95718996) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.808271) q[1];
sx q[1];
rz(-1.0153607) q[1];
sx q[1];
rz(-2.9395648) q[1];
rz(-pi) q[2];
x q[2];
rz(1.990063) q[3];
sx q[3];
rz(-1.5677088) q[3];
sx q[3];
rz(-3.1379238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3746419) q[2];
sx q[2];
rz(-0.16051126) q[2];
sx q[2];
rz(-3.0373419) q[2];
rz(-2.5637964) q[3];
sx q[3];
rz(-2.0640852) q[3];
sx q[3];
rz(-2.2195063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0084956) q[0];
sx q[0];
rz(-1.6660026) q[0];
sx q[0];
rz(-0.9444899) q[0];
rz(-1.7780444) q[1];
sx q[1];
rz(-1.3282158) q[1];
sx q[1];
rz(-1.69467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773374) q[0];
sx q[0];
rz(-2.2702571) q[0];
sx q[0];
rz(-0.2486458) q[0];
rz(-3.1051272) q[2];
sx q[2];
rz(-1.3205832) q[2];
sx q[2];
rz(1.6584876) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3376278) q[1];
sx q[1];
rz(-1.9586438) q[1];
sx q[1];
rz(-0.66168534) q[1];
rz(-pi) q[2];
rz(-2.1423316) q[3];
sx q[3];
rz(-1.0715228) q[3];
sx q[3];
rz(-1.5443791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.038534433) q[2];
sx q[2];
rz(-2.6699622) q[2];
sx q[2];
rz(0.72059694) q[2];
rz(-0.12901148) q[3];
sx q[3];
rz(-1.7804264) q[3];
sx q[3];
rz(1.8425997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98220372) q[0];
sx q[0];
rz(-0.99697462) q[0];
sx q[0];
rz(0.61656117) q[0];
rz(-2.189134) q[1];
sx q[1];
rz(-1.9472232) q[1];
sx q[1];
rz(-1.2430826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4774507) q[0];
sx q[0];
rz(-1.7074025) q[0];
sx q[0];
rz(1.3090429) q[0];
x q[1];
rz(-2.0427734) q[2];
sx q[2];
rz(-2.3772559) q[2];
sx q[2];
rz(2.7175222) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.19784099) q[1];
sx q[1];
rz(-1.5103522) q[1];
sx q[1];
rz(0.72058679) q[1];
rz(-pi) q[2];
rz(3.0674446) q[3];
sx q[3];
rz(-0.95063248) q[3];
sx q[3];
rz(0.95297232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2359961) q[2];
sx q[2];
rz(-2.8541028) q[2];
sx q[2];
rz(-3.0671885) q[2];
rz(-1.0807886) q[3];
sx q[3];
rz(-1.6483043) q[3];
sx q[3];
rz(-1.0409482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8393132) q[0];
sx q[0];
rz(-1.9957207) q[0];
sx q[0];
rz(-0.062051274) q[0];
rz(1.4150367) q[1];
sx q[1];
rz(-0.77295417) q[1];
sx q[1];
rz(3.0516023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4231071) q[0];
sx q[0];
rz(-0.88750091) q[0];
sx q[0];
rz(-0.32916577) q[0];
rz(-2.2587218) q[2];
sx q[2];
rz(-1.4423048) q[2];
sx q[2];
rz(-1.524432) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94839986) q[1];
sx q[1];
rz(-0.71454778) q[1];
sx q[1];
rz(2.1195861) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2136647) q[3];
sx q[3];
rz(-1.0786329) q[3];
sx q[3];
rz(-1.5318235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0061079582) q[2];
sx q[2];
rz(-1.2697479) q[2];
sx q[2];
rz(0.6832068) q[2];
rz(1.3225383) q[3];
sx q[3];
rz(-0.95219487) q[3];
sx q[3];
rz(-1.3873842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6739864) q[0];
sx q[0];
rz(-1.645393) q[0];
sx q[0];
rz(0.4678539) q[0];
rz(0.6894919) q[1];
sx q[1];
rz(-2.4263224) q[1];
sx q[1];
rz(-2.668344) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5681385) q[0];
sx q[0];
rz(-1.3949419) q[0];
sx q[0];
rz(-0.071608105) q[0];
x q[1];
rz(-0.44395776) q[2];
sx q[2];
rz(-1.7164162) q[2];
sx q[2];
rz(-2.1514926) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2040904) q[1];
sx q[1];
rz(-0.023462208) q[1];
sx q[1];
rz(2.2949982) q[1];
x q[2];
rz(-1.8256447) q[3];
sx q[3];
rz(-1.2989961) q[3];
sx q[3];
rz(2.013828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1271992) q[2];
sx q[2];
rz(-1.8202123) q[2];
sx q[2];
rz(0.42352208) q[2];
rz(1.8262156) q[3];
sx q[3];
rz(-0.49153057) q[3];
sx q[3];
rz(-1.2929085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92726436) q[0];
sx q[0];
rz(-0.94896999) q[0];
sx q[0];
rz(1.8073136) q[0];
rz(-1.3182053) q[1];
sx q[1];
rz(-2.2315836) q[1];
sx q[1];
rz(0.013890161) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2023808) q[0];
sx q[0];
rz(-1.9343543) q[0];
sx q[0];
rz(3.1143051) q[0];
rz(-pi) q[1];
rz(-1.4133164) q[2];
sx q[2];
rz(-1.9183927) q[2];
sx q[2];
rz(1.9633788) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.397032) q[1];
sx q[1];
rz(-2.5066527) q[1];
sx q[1];
rz(1.145361) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.708843) q[3];
sx q[3];
rz(-2.4073753) q[3];
sx q[3];
rz(-2.6003797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58005971) q[2];
sx q[2];
rz(-1.9999028) q[2];
sx q[2];
rz(1.053099) q[2];
rz(0.92755353) q[3];
sx q[3];
rz(-1.2404697) q[3];
sx q[3];
rz(-2.659306) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841566) q[0];
sx q[0];
rz(-2.744839) q[0];
sx q[0];
rz(-0.60047737) q[0];
rz(0.26957574) q[1];
sx q[1];
rz(-0.69200626) q[1];
sx q[1];
rz(-1.762766) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3286014) q[0];
sx q[0];
rz(-2.6909851) q[0];
sx q[0];
rz(-2.2297321) q[0];
x q[1];
rz(1.8767886) q[2];
sx q[2];
rz(-1.2160794) q[2];
sx q[2];
rz(-1.4504832) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8104756) q[1];
sx q[1];
rz(-1.1441191) q[1];
sx q[1];
rz(1.3699378) q[1];
rz(-2.5559101) q[3];
sx q[3];
rz(-0.70428169) q[3];
sx q[3];
rz(1.3552988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1158585) q[2];
sx q[2];
rz(-2.8074042) q[2];
sx q[2];
rz(2.1791229) q[2];
rz(2.1879503) q[3];
sx q[3];
rz(-2.0680659) q[3];
sx q[3];
rz(-1.7474705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3772603) q[0];
sx q[0];
rz(-2.2567516) q[0];
sx q[0];
rz(1.2010038) q[0];
rz(1.3537021) q[1];
sx q[1];
rz(-1.9722912) q[1];
sx q[1];
rz(-0.73696662) q[1];
rz(-0.66966343) q[2];
sx q[2];
rz(-2.1047932) q[2];
sx q[2];
rz(1.3303458) q[2];
rz(0.50696452) q[3];
sx q[3];
rz(-1.6929814) q[3];
sx q[3];
rz(-0.89043989) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
