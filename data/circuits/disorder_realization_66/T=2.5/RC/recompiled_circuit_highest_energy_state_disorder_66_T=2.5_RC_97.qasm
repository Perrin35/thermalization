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
rz(2.1386327) q[0];
sx q[0];
rz(-0.47337368) q[0];
sx q[0];
rz(2.0546497) q[0];
rz(-0.76397693) q[1];
sx q[1];
rz(-0.44372258) q[1];
sx q[1];
rz(2.3102163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2646541) q[0];
sx q[0];
rz(-1.042141) q[0];
sx q[0];
rz(1.009788) q[0];
rz(-pi) q[1];
rz(1.1016721) q[2];
sx q[2];
rz(-1.7828825) q[2];
sx q[2];
rz(2.9123127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8352141) q[1];
sx q[1];
rz(-1.467538) q[1];
sx q[1];
rz(1.4500965) q[1];
x q[2];
rz(-1.5441203) q[3];
sx q[3];
rz(-1.4988588) q[3];
sx q[3];
rz(-1.5213336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5272556) q[2];
sx q[2];
rz(-1.8680806) q[2];
sx q[2];
rz(0.41198507) q[2];
rz(0.24474239) q[3];
sx q[3];
rz(-0.63260806) q[3];
sx q[3];
rz(2.8542724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43993846) q[0];
sx q[0];
rz(-2.169862) q[0];
sx q[0];
rz(0.42020759) q[0];
rz(0.48928753) q[1];
sx q[1];
rz(-2.2261765) q[1];
sx q[1];
rz(1.9583826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91285465) q[0];
sx q[0];
rz(-1.1381554) q[0];
sx q[0];
rz(2.2590524) q[0];
x q[1];
rz(1.2474044) q[2];
sx q[2];
rz(-1.8696491) q[2];
sx q[2];
rz(-0.0060723532) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.816531) q[1];
sx q[1];
rz(-1.9516212) q[1];
sx q[1];
rz(-0.62398361) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97789852) q[3];
sx q[3];
rz(-0.20836941) q[3];
sx q[3];
rz(0.32754242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6861787) q[2];
sx q[2];
rz(-1.6697465) q[2];
sx q[2];
rz(-0.14032826) q[2];
rz(-0.27648196) q[3];
sx q[3];
rz(-0.77767196) q[3];
sx q[3];
rz(-2.8592143) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8088733) q[0];
sx q[0];
rz(-0.35950279) q[0];
sx q[0];
rz(-1.7669539) q[0];
rz(0.91284347) q[1];
sx q[1];
rz(-2.5990867) q[1];
sx q[1];
rz(-1.0309781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6767117) q[0];
sx q[0];
rz(-1.2455938) q[0];
sx q[0];
rz(-0.19285481) q[0];
rz(-pi) q[1];
rz(-0.39856492) q[2];
sx q[2];
rz(-2.0953278) q[2];
sx q[2];
rz(-0.87510671) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43592855) q[1];
sx q[1];
rz(-1.1586274) q[1];
sx q[1];
rz(3.1306821) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52521962) q[3];
sx q[3];
rz(-0.37074836) q[3];
sx q[3];
rz(1.4300089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1294407) q[2];
sx q[2];
rz(-1.0544216) q[2];
sx q[2];
rz(-0.63808179) q[2];
rz(0.10073999) q[3];
sx q[3];
rz(-2.1295857) q[3];
sx q[3];
rz(2.6428599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3070372) q[0];
sx q[0];
rz(-1.6224253) q[0];
sx q[0];
rz(-1.3822973) q[0];
rz(2.8221829) q[1];
sx q[1];
rz(-1.5089792) q[1];
sx q[1];
rz(2.3841948) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46947458) q[0];
sx q[0];
rz(-1.4961157) q[0];
sx q[0];
rz(-1.5793708) q[0];
x q[1];
rz(2.9992835) q[2];
sx q[2];
rz(-1.6985296) q[2];
sx q[2];
rz(2.3960428) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0075052) q[1];
sx q[1];
rz(-2.0148835) q[1];
sx q[1];
rz(-2.215338) q[1];
x q[2];
rz(-2.5365165) q[3];
sx q[3];
rz(-2.0143442) q[3];
sx q[3];
rz(-0.098066948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77183977) q[2];
sx q[2];
rz(-2.2660793) q[2];
sx q[2];
rz(-0.80779752) q[2];
rz(-1.7317023) q[3];
sx q[3];
rz(-1.8818703) q[3];
sx q[3];
rz(0.65214777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0659502) q[0];
sx q[0];
rz(-2.763651) q[0];
sx q[0];
rz(2.156303) q[0];
rz(-1.6458884) q[1];
sx q[1];
rz(-1.138405) q[1];
sx q[1];
rz(-0.92934242) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1131234) q[0];
sx q[0];
rz(-1.3036393) q[0];
sx q[0];
rz(0.6031674) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1628289) q[2];
sx q[2];
rz(-1.0024602) q[2];
sx q[2];
rz(-1.6108244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95463175) q[1];
sx q[1];
rz(-1.371688) q[1];
sx q[1];
rz(-2.1084) q[1];
rz(-0.014300925) q[3];
sx q[3];
rz(-2.1319444) q[3];
sx q[3];
rz(-2.1545269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0890961) q[2];
sx q[2];
rz(-1.3546983) q[2];
sx q[2];
rz(-2.055114) q[2];
rz(0.9440445) q[3];
sx q[3];
rz(-0.090525301) q[3];
sx q[3];
rz(1.121678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1679967) q[0];
sx q[0];
rz(-2.7041628) q[0];
sx q[0];
rz(-0.22015372) q[0];
rz(1.6379697) q[1];
sx q[1];
rz(-1.817768) q[1];
sx q[1];
rz(-2.5880609) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65992113) q[0];
sx q[0];
rz(-1.685647) q[0];
sx q[0];
rz(0.19821313) q[0];
rz(-pi) q[1];
rz(-0.4398063) q[2];
sx q[2];
rz(-2.7436119) q[2];
sx q[2];
rz(-2.8619638) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.625861) q[1];
sx q[1];
rz(-1.6193832) q[1];
sx q[1];
rz(-2.4481985) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4559654) q[3];
sx q[3];
rz(-1.3925288) q[3];
sx q[3];
rz(-0.36742154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0452051) q[2];
sx q[2];
rz(-1.6738946) q[2];
sx q[2];
rz(0.20827797) q[2];
rz(-2.0221209) q[3];
sx q[3];
rz(-1.8549615) q[3];
sx q[3];
rz(2.313405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49523062) q[0];
sx q[0];
rz(-1.3883075) q[0];
sx q[0];
rz(-1.0736504) q[0];
rz(0.13058361) q[1];
sx q[1];
rz(-1.5675631) q[1];
sx q[1];
rz(-2.1317587) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70390095) q[0];
sx q[0];
rz(-2.7006375) q[0];
sx q[0];
rz(0.44619123) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45744894) q[2];
sx q[2];
rz(-0.1651131) q[2];
sx q[2];
rz(-3.0505863) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4676795) q[1];
sx q[1];
rz(-0.99930489) q[1];
sx q[1];
rz(-1.4614146) q[1];
rz(-pi) q[2];
rz(-1.9319264) q[3];
sx q[3];
rz(-2.2029366) q[3];
sx q[3];
rz(1.2120642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4449571) q[2];
sx q[2];
rz(-0.37142763) q[2];
sx q[2];
rz(2.8182287) q[2];
rz(2.4671386) q[3];
sx q[3];
rz(-0.89265299) q[3];
sx q[3];
rz(3.118012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-3.0167639) q[0];
sx q[0];
rz(-2.6478196) q[0];
sx q[0];
rz(-1.0892518) q[0];
rz(0.59843868) q[1];
sx q[1];
rz(-1.2804223) q[1];
sx q[1];
rz(2.3425897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4161637) q[0];
sx q[0];
rz(-1.874864) q[0];
sx q[0];
rz(-1.8644402) q[0];
rz(-pi) q[1];
rz(-2.2369821) q[2];
sx q[2];
rz(-1.4965726) q[2];
sx q[2];
rz(2.7192309) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.028513718) q[1];
sx q[1];
rz(-1.4220107) q[1];
sx q[1];
rz(2.8362175) q[1];
rz(-pi) q[2];
rz(0.24183065) q[3];
sx q[3];
rz(-1.5655091) q[3];
sx q[3];
rz(-1.9092803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66990718) q[2];
sx q[2];
rz(-2.6354463) q[2];
sx q[2];
rz(1.252582) q[2];
rz(-2.0910828) q[3];
sx q[3];
rz(-0.59058467) q[3];
sx q[3];
rz(-0.93241507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85225409) q[0];
sx q[0];
rz(-1.7409356) q[0];
sx q[0];
rz(2.6259212) q[0];
rz(1.4562666) q[1];
sx q[1];
rz(-1.3412424) q[1];
sx q[1];
rz(0.32404831) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0808602) q[0];
sx q[0];
rz(-1.8243084) q[0];
sx q[0];
rz(0.41187615) q[0];
rz(-pi) q[1];
x q[1];
rz(2.134974) q[2];
sx q[2];
rz(-0.7340275) q[2];
sx q[2];
rz(0.013040868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1157323) q[1];
sx q[1];
rz(-2.6792791) q[1];
sx q[1];
rz(1.2532534) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0131272) q[3];
sx q[3];
rz(-1.5032167) q[3];
sx q[3];
rz(0.24012071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9825762) q[2];
sx q[2];
rz(-2.0197208) q[2];
sx q[2];
rz(2.4328361) q[2];
rz(-2.3580264) q[3];
sx q[3];
rz(-1.2015899) q[3];
sx q[3];
rz(1.5172575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.6056972) q[0];
sx q[0];
rz(-2.504183) q[0];
sx q[0];
rz(-0.15644431) q[0];
rz(-2.5018196) q[1];
sx q[1];
rz(-2.4867058) q[1];
sx q[1];
rz(0.99008647) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7759686) q[0];
sx q[0];
rz(-2.1556615) q[0];
sx q[0];
rz(2.3192899) q[0];
x q[1];
rz(-1.3285816) q[2];
sx q[2];
rz(-2.7822084) q[2];
sx q[2];
rz(-0.49312544) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1875547) q[1];
sx q[1];
rz(-2.1219664) q[1];
sx q[1];
rz(-1.4146845) q[1];
rz(-pi) q[2];
rz(-2.5699432) q[3];
sx q[3];
rz(-2.902973) q[3];
sx q[3];
rz(2.1941136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5618374) q[2];
sx q[2];
rz(-1.729915) q[2];
sx q[2];
rz(0.61335316) q[2];
rz(0.26398811) q[3];
sx q[3];
rz(-1.3365021) q[3];
sx q[3];
rz(2.4700375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092125208) q[0];
sx q[0];
rz(-1.6138074) q[0];
sx q[0];
rz(-2.001979) q[0];
rz(-1.4641948) q[1];
sx q[1];
rz(-0.72036998) q[1];
sx q[1];
rz(1.5710685) q[1];
rz(2.3650305) q[2];
sx q[2];
rz(-0.53755098) q[2];
sx q[2];
rz(2.9992486) q[2];
rz(2.522884) q[3];
sx q[3];
rz(-1.6029463) q[3];
sx q[3];
rz(-1.1986986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
