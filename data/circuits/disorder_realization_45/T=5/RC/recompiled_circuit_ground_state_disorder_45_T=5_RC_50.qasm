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
rz(0.80945102) q[1];
sx q[1];
rz(-1.5341772) q[1];
sx q[1];
rz(3.1117575) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.580509) q[0];
sx q[0];
rz(-2.5445815) q[0];
sx q[0];
rz(2.719814) q[0];
rz(-pi) q[1];
rz(-0.81852976) q[2];
sx q[2];
rz(-1.5191744) q[2];
sx q[2];
rz(2.7843786) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7368841) q[1];
sx q[1];
rz(-0.88062693) q[1];
sx q[1];
rz(1.0600914) q[1];
rz(-1.6601059) q[3];
sx q[3];
rz(-1.7062643) q[3];
sx q[3];
rz(-0.40218807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5030824) q[2];
sx q[2];
rz(-2.5677887) q[2];
sx q[2];
rz(-0.32162515) q[2];
rz(0.62421787) q[3];
sx q[3];
rz(-1.7523242) q[3];
sx q[3];
rz(-1.6577087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44840789) q[0];
sx q[0];
rz(-0.25633651) q[0];
sx q[0];
rz(-0.93628991) q[0];
rz(-1.2721277) q[1];
sx q[1];
rz(-1.5446168) q[1];
sx q[1];
rz(-1.3880656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4988588) q[0];
sx q[0];
rz(-1.513645) q[0];
sx q[0];
rz(-1.3512035) q[0];
x q[1];
rz(-0.53728062) q[2];
sx q[2];
rz(-0.85308077) q[2];
sx q[2];
rz(-3.0514015) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9656877) q[1];
sx q[1];
rz(-2.2908194) q[1];
sx q[1];
rz(0.73890025) q[1];
x q[2];
rz(-1.306698) q[3];
sx q[3];
rz(-1.3624884) q[3];
sx q[3];
rz(-1.1549319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0973413) q[2];
sx q[2];
rz(-0.96434957) q[2];
sx q[2];
rz(0.67330366) q[2];
rz(-1.3168969) q[3];
sx q[3];
rz(-1.8484867) q[3];
sx q[3];
rz(2.3015658) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3299385) q[0];
sx q[0];
rz(-1.5978156) q[0];
sx q[0];
rz(-0.04976186) q[0];
rz(0.20637575) q[1];
sx q[1];
rz(-1.754909) q[1];
sx q[1];
rz(-1.4899563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1288717) q[0];
sx q[0];
rz(-1.6908592) q[0];
sx q[0];
rz(-0.61693667) q[0];
rz(-pi) q[1];
rz(-2.4902053) q[2];
sx q[2];
rz(-2.4804401) q[2];
sx q[2];
rz(-0.1729473) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7327642) q[1];
sx q[1];
rz(-1.2868501) q[1];
sx q[1];
rz(2.8676621) q[1];
rz(-2.8036267) q[3];
sx q[3];
rz(-2.0062048) q[3];
sx q[3];
rz(1.1359648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18996198) q[2];
sx q[2];
rz(-0.54533521) q[2];
sx q[2];
rz(0.0541617) q[2];
rz(1.1545898) q[3];
sx q[3];
rz(-1.8535987) q[3];
sx q[3];
rz(-1.1156999) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
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
rz(-0.94991389) q[1];
sx q[1];
rz(-1.3309006) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9337263) q[0];
sx q[0];
rz(-2.8678992) q[0];
sx q[0];
rz(-0.52723186) q[0];
rz(-1.0687549) q[2];
sx q[2];
rz(-0.34345657) q[2];
sx q[2];
rz(-0.95718996) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96283276) q[1];
sx q[1];
rz(-2.5542027) q[1];
sx q[1];
rz(-1.8835095) q[1];
rz(-pi) q[2];
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
rz(1.3746419) q[2];
sx q[2];
rz(-2.9810814) q[2];
sx q[2];
rz(0.10425076) q[2];
rz(-2.5637964) q[3];
sx q[3];
rz(-2.0640852) q[3];
sx q[3];
rz(-2.2195063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0084956) q[0];
sx q[0];
rz(-1.47559) q[0];
sx q[0];
rz(-2.1971028) q[0];
rz(-1.7780444) q[1];
sx q[1];
rz(-1.3282158) q[1];
sx q[1];
rz(-1.69467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9642553) q[0];
sx q[0];
rz(-0.87133555) q[0];
sx q[0];
rz(-0.2486458) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.036465444) q[2];
sx q[2];
rz(-1.8210095) q[2];
sx q[2];
rz(-1.4831051) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0883184) q[1];
sx q[1];
rz(-0.96573869) q[1];
sx q[1];
rz(-1.0929918) q[1];
rz(-pi) q[2];
rz(0.78131494) q[3];
sx q[3];
rz(-2.4014946) q[3];
sx q[3];
rz(0.61352713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1030582) q[2];
sx q[2];
rz(-2.6699622) q[2];
sx q[2];
rz(-0.72059694) q[2];
rz(-3.0125812) q[3];
sx q[3];
rz(-1.3611662) q[3];
sx q[3];
rz(1.8425997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1593889) q[0];
sx q[0];
rz(-0.99697462) q[0];
sx q[0];
rz(-0.61656117) q[0];
rz(-2.189134) q[1];
sx q[1];
rz(-1.1943694) q[1];
sx q[1];
rz(1.2430826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6641419) q[0];
sx q[0];
rz(-1.7074025) q[0];
sx q[0];
rz(1.3090429) q[0];
x q[1];
rz(-0.86402969) q[2];
sx q[2];
rz(-1.2507157) q[2];
sx q[2];
rz(-2.3478986) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3042708) q[1];
sx q[1];
rz(-2.4189286) q[1];
sx q[1];
rz(-3.0501306) q[1];
rz(2.1922621) q[3];
sx q[3];
rz(-1.5104745) q[3];
sx q[3];
rz(2.5669135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9055966) q[2];
sx q[2];
rz(-0.2874898) q[2];
sx q[2];
rz(3.0671885) q[2];
rz(2.060804) q[3];
sx q[3];
rz(-1.4932884) q[3];
sx q[3];
rz(-2.1006445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3022795) q[0];
sx q[0];
rz(-1.1458719) q[0];
sx q[0];
rz(-0.062051274) q[0];
rz(1.4150367) q[1];
sx q[1];
rz(-2.3686385) q[1];
sx q[1];
rz(0.089990377) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4231071) q[0];
sx q[0];
rz(-0.88750091) q[0];
sx q[0];
rz(-2.8124269) q[0];
rz(-pi) q[1];
rz(-0.16570602) q[2];
sx q[2];
rz(-2.2519654) q[2];
sx q[2];
rz(0.15128862) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.18971009) q[1];
sx q[1];
rz(-1.9196577) q[1];
sx q[1];
rz(2.2079218) q[1];
x q[2];
rz(-0.84109938) q[3];
sx q[3];
rz(-2.3537618) q[3];
sx q[3];
rz(-0.5238409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0061079582) q[2];
sx q[2];
rz(-1.8718448) q[2];
sx q[2];
rz(2.4583859) q[2];
rz(-1.8190544) q[3];
sx q[3];
rz(-0.95219487) q[3];
sx q[3];
rz(-1.3873842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6739864) q[0];
sx q[0];
rz(-1.4961996) q[0];
sx q[0];
rz(0.4678539) q[0];
rz(-0.6894919) q[1];
sx q[1];
rz(-2.4263224) q[1];
sx q[1];
rz(2.668344) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9625585) q[0];
sx q[0];
rz(-0.18973543) q[0];
sx q[0];
rz(1.1880072) q[0];
rz(-0.44395776) q[2];
sx q[2];
rz(-1.7164162) q[2];
sx q[2];
rz(-2.1514926) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6618408) q[1];
sx q[1];
rz(-1.5532232) q[1];
sx q[1];
rz(3.1260465) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28041081) q[3];
sx q[3];
rz(-1.3254998) q[3];
sx q[3];
rz(-2.6287358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1271992) q[2];
sx q[2];
rz(-1.3213804) q[2];
sx q[2];
rz(-0.42352208) q[2];
rz(-1.315377) q[3];
sx q[3];
rz(-0.49153057) q[3];
sx q[3];
rz(-1.2929085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2143283) q[0];
sx q[0];
rz(-2.1926227) q[0];
sx q[0];
rz(1.8073136) q[0];
rz(-1.3182053) q[1];
sx q[1];
rz(-0.91000906) q[1];
sx q[1];
rz(-0.013890161) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5197138) q[0];
sx q[0];
rz(-1.5452928) q[0];
sx q[0];
rz(1.9344781) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7282762) q[2];
sx q[2];
rz(-1.2232) q[2];
sx q[2];
rz(-1.1782139) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1760343) q[1];
sx q[1];
rz(-1.8181043) q[1];
sx q[1];
rz(-2.1618188) q[1];
rz(-pi) q[2];
rz(-2.3002616) q[3];
sx q[3];
rz(-1.6631261) q[3];
sx q[3];
rz(2.009237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5615329) q[2];
sx q[2];
rz(-1.1416898) q[2];
sx q[2];
rz(1.053099) q[2];
rz(-2.2140391) q[3];
sx q[3];
rz(-1.2404697) q[3];
sx q[3];
rz(-2.659306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.7841566) q[0];
sx q[0];
rz(-2.744839) q[0];
sx q[0];
rz(-0.60047737) q[0];
rz(-2.8720169) q[1];
sx q[1];
rz(-2.4495864) q[1];
sx q[1];
rz(1.762766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3286014) q[0];
sx q[0];
rz(-0.4506076) q[0];
sx q[0];
rz(0.91186055) q[0];
rz(-1.264804) q[2];
sx q[2];
rz(-1.2160794) q[2];
sx q[2];
rz(-1.4504832) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8104756) q[1];
sx q[1];
rz(-1.1441191) q[1];
sx q[1];
rz(-1.7716549) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0098712) q[3];
sx q[3];
rz(-1.0008662) q[3];
sx q[3];
rz(0.63907335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1158585) q[2];
sx q[2];
rz(-2.8074042) q[2];
sx q[2];
rz(-2.1791229) q[2];
rz(-0.95364237) q[3];
sx q[3];
rz(-1.0735268) q[3];
sx q[3];
rz(1.7474705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(1.3772603) q[0];
sx q[0];
rz(-2.2567516) q[0];
sx q[0];
rz(1.2010038) q[0];
rz(-1.7878905) q[1];
sx q[1];
rz(-1.9722912) q[1];
sx q[1];
rz(-0.73696662) q[1];
rz(0.66966343) q[2];
sx q[2];
rz(-1.0367994) q[2];
sx q[2];
rz(-1.8112469) q[2];
rz(2.8938724) q[3];
sx q[3];
rz(-0.52023028) q[3];
sx q[3];
rz(0.46432555) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
