OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1096119) q[0];
sx q[0];
rz(-1.0053827) q[0];
sx q[0];
rz(0.96860743) q[0];
rz(0.80945102) q[1];
sx q[1];
rz(-1.5341772) q[1];
sx q[1];
rz(3.1117575) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.077622) q[0];
sx q[0];
rz(-1.0322303) q[0];
sx q[0];
rz(-1.2993815) q[0];
x q[1];
rz(0.81852976) q[2];
sx q[2];
rz(-1.5191744) q[2];
sx q[2];
rz(0.35721401) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.40470859) q[1];
sx q[1];
rz(-2.2609657) q[1];
sx q[1];
rz(-1.0600914) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4814867) q[3];
sx q[3];
rz(-1.7062643) q[3];
sx q[3];
rz(2.7394046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6385103) q[2];
sx q[2];
rz(-2.5677887) q[2];
sx q[2];
rz(-0.32162515) q[2];
rz(-2.5173748) q[3];
sx q[3];
rz(-1.3892684) q[3];
sx q[3];
rz(-1.483884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6931848) q[0];
sx q[0];
rz(-2.8852561) q[0];
sx q[0];
rz(-2.2053027) q[0];
rz(-1.2721277) q[1];
sx q[1];
rz(-1.5446168) q[1];
sx q[1];
rz(-1.3880656) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4988588) q[0];
sx q[0];
rz(-1.513645) q[0];
sx q[0];
rz(-1.7903891) q[0];
rz(0.77735591) q[2];
sx q[2];
rz(-1.1749975) q[2];
sx q[2];
rz(1.1071902) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2876828) q[1];
sx q[1];
rz(-2.1016995) q[1];
sx q[1];
rz(-2.4413051) q[1];
rz(2.2514492) q[3];
sx q[3];
rz(-0.33484866) q[3];
sx q[3];
rz(1.0687506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0442514) q[2];
sx q[2];
rz(-2.1772431) q[2];
sx q[2];
rz(2.468289) q[2];
rz(1.8246957) q[3];
sx q[3];
rz(-1.293106) q[3];
sx q[3];
rz(-2.3015658) q[3];
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
rz(0.3299385) q[0];
sx q[0];
rz(-1.543777) q[0];
sx q[0];
rz(0.04976186) q[0];
rz(0.20637575) q[1];
sx q[1];
rz(-1.754909) q[1];
sx q[1];
rz(-1.4899563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4987652) q[0];
sx q[0];
rz(-0.9589566) q[0];
sx q[0];
rz(1.7176413) q[0];
rz(-0.65138731) q[2];
sx q[2];
rz(-2.4804401) q[2];
sx q[2];
rz(0.1729473) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.94585431) q[1];
sx q[1];
rz(-2.7496413) q[1];
sx q[1];
rz(2.3183104) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0288896) q[3];
sx q[3];
rz(-1.2654403) q[3];
sx q[3];
rz(2.8539243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18996198) q[2];
sx q[2];
rz(-2.5962574) q[2];
sx q[2];
rz(3.087431) q[2];
rz(1.1545898) q[3];
sx q[3];
rz(-1.8535987) q[3];
sx q[3];
rz(-1.1156999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(0.40256777) q[0];
sx q[0];
rz(-1.4250647) q[0];
sx q[0];
rz(2.2836852) q[0];
rz(-0.67445451) q[1];
sx q[1];
rz(-0.94991389) q[1];
sx q[1];
rz(-1.8106921) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9337263) q[0];
sx q[0];
rz(-2.8678992) q[0];
sx q[0];
rz(-2.6143608) q[0];
rz(-pi) q[1];
rz(1.8745918) q[2];
sx q[2];
rz(-1.7335606) q[2];
sx q[2];
rz(3.0050584) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.34506306) q[1];
sx q[1];
rz(-1.3994675) q[1];
sx q[1];
rz(-1.0061128) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5783806) q[3];
sx q[3];
rz(-0.41927734) q[3];
sx q[3];
rz(1.5602001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3746419) q[2];
sx q[2];
rz(-0.16051126) q[2];
sx q[2];
rz(-3.0373419) q[2];
rz(-0.57779622) q[3];
sx q[3];
rz(-1.0775074) q[3];
sx q[3];
rz(0.92208636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0084956) q[0];
sx q[0];
rz(-1.6660026) q[0];
sx q[0];
rz(-0.9444899) q[0];
rz(-1.3635483) q[1];
sx q[1];
rz(-1.3282158) q[1];
sx q[1];
rz(1.69467) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76857012) q[0];
sx q[0];
rz(-1.7602341) q[0];
sx q[0];
rz(-2.2856838) q[0];
rz(-pi) q[1];
rz(-1.4291006) q[2];
sx q[2];
rz(-0.25280127) q[2];
sx q[2];
rz(-1.5122053) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0532743) q[1];
sx q[1];
rz(-0.96573869) q[1];
sx q[1];
rz(-2.0486008) q[1];
rz(-pi) q[2];
rz(2.3602777) q[3];
sx q[3];
rz(-2.4014946) q[3];
sx q[3];
rz(-0.61352713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.038534433) q[2];
sx q[2];
rz(-0.47163042) q[2];
sx q[2];
rz(0.72059694) q[2];
rz(3.0125812) q[3];
sx q[3];
rz(-1.3611662) q[3];
sx q[3];
rz(1.2989929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1593889) q[0];
sx q[0];
rz(-0.99697462) q[0];
sx q[0];
rz(-0.61656117) q[0];
rz(2.189134) q[1];
sx q[1];
rz(-1.1943694) q[1];
sx q[1];
rz(1.8985101) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0117802) q[0];
sx q[0];
rz(-1.3115378) q[0];
sx q[0];
rz(-3.0002322) q[0];
x q[1];
rz(0.41105777) q[2];
sx q[2];
rz(-0.90663821) q[2];
sx q[2];
rz(1.039584) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3042708) q[1];
sx q[1];
rz(-2.4189286) q[1];
sx q[1];
rz(-3.0501306) q[1];
rz(-pi) q[2];
rz(1.6741577) q[3];
sx q[3];
rz(-0.62400104) q[3];
sx q[3];
rz(-2.0614908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.2359961) q[2];
sx q[2];
rz(-2.8541028) q[2];
sx q[2];
rz(-0.074404152) q[2];
rz(-1.0807886) q[3];
sx q[3];
rz(-1.4932884) q[3];
sx q[3];
rz(-2.1006445) q[3];
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
rz(pi/2) q[3];
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
rz(-2.8393132) q[0];
sx q[0];
rz(-1.9957207) q[0];
sx q[0];
rz(-0.062051274) q[0];
rz(1.7265559) q[1];
sx q[1];
rz(-2.3686385) q[1];
sx q[1];
rz(3.0516023) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4231071) q[0];
sx q[0];
rz(-0.88750091) q[0];
sx q[0];
rz(-0.32916577) q[0];
rz(0.88287087) q[2];
sx q[2];
rz(-1.6992879) q[2];
sx q[2];
rz(-1.6171607) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5127232) q[1];
sx q[1];
rz(-0.97755331) q[1];
sx q[1];
rz(-2.7166461) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5513583) q[3];
sx q[3];
rz(-2.1274421) q[3];
sx q[3];
rz(2.7624453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0061079582) q[2];
sx q[2];
rz(-1.2697479) q[2];
sx q[2];
rz(-2.4583859) q[2];
rz(1.8190544) q[3];
sx q[3];
rz(-2.1893978) q[3];
sx q[3];
rz(-1.3873842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4676062) q[0];
sx q[0];
rz(-1.4961996) q[0];
sx q[0];
rz(-2.6737387) q[0];
rz(0.6894919) q[1];
sx q[1];
rz(-2.4263224) q[1];
sx q[1];
rz(-2.668344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1263862) q[0];
sx q[0];
rz(-1.6412982) q[0];
sx q[0];
rz(-1.7470933) q[0];
x q[1];
rz(1.4098005) q[2];
sx q[2];
rz(-1.1318665) q[2];
sx q[2];
rz(2.4919897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0508214) q[1];
sx q[1];
rz(-1.5863401) q[1];
sx q[1];
rz(-1.5883716) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7353294) q[3];
sx q[3];
rz(-2.7711754) q[3];
sx q[3];
rz(-0.35740023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0143934) q[2];
sx q[2];
rz(-1.3213804) q[2];
sx q[2];
rz(2.7180706) q[2];
rz(1.8262156) q[3];
sx q[3];
rz(-2.6500621) q[3];
sx q[3];
rz(-1.8486842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92726436) q[0];
sx q[0];
rz(-2.1926227) q[0];
sx q[0];
rz(-1.8073136) q[0];
rz(1.8233874) q[1];
sx q[1];
rz(-0.91000906) q[1];
sx q[1];
rz(3.1277025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62187884) q[0];
sx q[0];
rz(-1.5962999) q[0];
sx q[0];
rz(-1.2071146) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7330796) q[2];
sx q[2];
rz(-2.7613104) q[2];
sx q[2];
rz(-0.74199669) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74456066) q[1];
sx q[1];
rz(-2.5066527) q[1];
sx q[1];
rz(-1.9962316) q[1];
x q[2];
rz(1.708843) q[3];
sx q[3];
rz(-2.4073753) q[3];
sx q[3];
rz(2.6003797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5615329) q[2];
sx q[2];
rz(-1.9999028) q[2];
sx q[2];
rz(2.0884936) q[2];
rz(2.2140391) q[3];
sx q[3];
rz(-1.901123) q[3];
sx q[3];
rz(-2.659306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841566) q[0];
sx q[0];
rz(-2.744839) q[0];
sx q[0];
rz(0.60047737) q[0];
rz(-2.8720169) q[1];
sx q[1];
rz(-0.69200626) q[1];
sx q[1];
rz(-1.762766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81299128) q[0];
sx q[0];
rz(-0.4506076) q[0];
sx q[0];
rz(-2.2297321) q[0];
rz(-2.7711033) q[2];
sx q[2];
rz(-1.8571808) q[2];
sx q[2];
rz(-3.1305673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32374183) q[1];
sx q[1];
rz(-1.7534337) q[1];
sx q[1];
rz(-2.707213) q[1];
x q[2];
rz(2.5254978) q[3];
sx q[3];
rz(-1.9368251) q[3];
sx q[3];
rz(-0.6835365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.025734162) q[2];
sx q[2];
rz(-2.8074042) q[2];
sx q[2];
rz(-2.1791229) q[2];
rz(2.1879503) q[3];
sx q[3];
rz(-1.0735268) q[3];
sx q[3];
rz(-1.3941221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.7643323) q[0];
sx q[0];
rz(-0.88484103) q[0];
sx q[0];
rz(-1.9405889) q[0];
rz(-1.7878905) q[1];
sx q[1];
rz(-1.9722912) q[1];
sx q[1];
rz(-0.73696662) q[1];
rz(-0.66966343) q[2];
sx q[2];
rz(-2.1047932) q[2];
sx q[2];
rz(1.3303458) q[2];
rz(-0.50696452) q[3];
sx q[3];
rz(-1.4486113) q[3];
sx q[3];
rz(2.2511528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
