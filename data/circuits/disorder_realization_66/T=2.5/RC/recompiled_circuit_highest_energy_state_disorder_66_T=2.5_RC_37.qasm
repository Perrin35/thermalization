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
rz(-1.00296) q[0];
sx q[0];
rz(-2.668219) q[0];
sx q[0];
rz(1.0869429) q[0];
rz(2.3776157) q[1];
sx q[1];
rz(-2.6978701) q[1];
sx q[1];
rz(0.8313764) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017480306) q[0];
sx q[0];
rz(-2.3907733) q[0];
sx q[0];
rz(-0.73877891) q[0];
x q[1];
rz(-0.23687266) q[2];
sx q[2];
rz(-2.028596) q[2];
sx q[2];
rz(1.906369) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1727816) q[1];
sx q[1];
rz(-0.15867844) q[1];
sx q[1];
rz(2.2814343) q[1];
rz(-0.3545019) q[3];
sx q[3];
rz(-3.0648764) q[3];
sx q[3];
rz(1.8767954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5272556) q[2];
sx q[2];
rz(-1.8680806) q[2];
sx q[2];
rz(2.7296076) q[2];
rz(-0.24474239) q[3];
sx q[3];
rz(-2.5089846) q[3];
sx q[3];
rz(2.8542724) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43993846) q[0];
sx q[0];
rz(-0.97173062) q[0];
sx q[0];
rz(2.7213851) q[0];
rz(0.48928753) q[1];
sx q[1];
rz(-2.2261765) q[1];
sx q[1];
rz(1.9583826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32589697) q[0];
sx q[0];
rz(-2.1854379) q[0];
sx q[0];
rz(-0.53890284) q[0];
rz(-pi) q[1];
rz(2.8274306) q[2];
sx q[2];
rz(-1.8793686) q[2];
sx q[2];
rz(1.6752288) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3250617) q[1];
sx q[1];
rz(-1.1899714) q[1];
sx q[1];
rz(0.62398361) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.39721) q[3];
sx q[3];
rz(-1.6866444) q[3];
sx q[3];
rz(-0.66044745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6861787) q[2];
sx q[2];
rz(-1.6697465) q[2];
sx q[2];
rz(0.14032826) q[2];
rz(2.8651107) q[3];
sx q[3];
rz(-0.77767196) q[3];
sx q[3];
rz(0.28237835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8088733) q[0];
sx q[0];
rz(-0.35950279) q[0];
sx q[0];
rz(1.3746388) q[0];
rz(-0.91284347) q[1];
sx q[1];
rz(-0.54250598) q[1];
sx q[1];
rz(2.1106145) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6767117) q[0];
sx q[0];
rz(-1.8959989) q[0];
sx q[0];
rz(-2.9487378) q[0];
x q[1];
rz(2.1616251) q[2];
sx q[2];
rz(-0.64729948) q[2];
sx q[2];
rz(1.574263) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.678434) q[1];
sx q[1];
rz(-0.41230508) q[1];
sx q[1];
rz(-1.5458471) q[1];
rz(-1.7632906) q[3];
sx q[3];
rz(-1.88965) q[3];
sx q[3];
rz(-1.1553193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.012152) q[2];
sx q[2];
rz(-2.0871711) q[2];
sx q[2];
rz(0.63808179) q[2];
rz(0.10073999) q[3];
sx q[3];
rz(-1.0120069) q[3];
sx q[3];
rz(-2.6428599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3070372) q[0];
sx q[0];
rz(-1.6224253) q[0];
sx q[0];
rz(-1.3822973) q[0];
rz(0.31940976) q[1];
sx q[1];
rz(-1.5089792) q[1];
sx q[1];
rz(-2.3841948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0396311) q[0];
sx q[0];
rz(-1.5622458) q[0];
sx q[0];
rz(3.0669093) q[0];
rz(2.9992835) q[2];
sx q[2];
rz(-1.443063) q[2];
sx q[2];
rz(-2.3960428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1340874) q[1];
sx q[1];
rz(-1.1267091) q[1];
sx q[1];
rz(0.92625463) q[1];
rz(2.5365165) q[3];
sx q[3];
rz(-2.0143442) q[3];
sx q[3];
rz(0.098066948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77183977) q[2];
sx q[2];
rz(-2.2660793) q[2];
sx q[2];
rz(0.80779752) q[2];
rz(1.7317023) q[3];
sx q[3];
rz(-1.2597224) q[3];
sx q[3];
rz(-2.4894449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0659502) q[0];
sx q[0];
rz(-2.763651) q[0];
sx q[0];
rz(-2.156303) q[0];
rz(1.4957042) q[1];
sx q[1];
rz(-1.138405) q[1];
sx q[1];
rz(2.2122502) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1131234) q[0];
sx q[0];
rz(-1.3036393) q[0];
sx q[0];
rz(-0.6031674) q[0];
rz(-pi) q[1];
rz(-1.9787637) q[2];
sx q[2];
rz(-2.1391324) q[2];
sx q[2];
rz(-1.6108244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8457736) q[1];
sx q[1];
rz(-2.5717178) q[1];
sx q[1];
rz(-1.946158) q[1];
x q[2];
rz(-1.5480488) q[3];
sx q[3];
rz(-2.5802819) q[3];
sx q[3];
rz(-2.1813957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0890961) q[2];
sx q[2];
rz(-1.3546983) q[2];
sx q[2];
rz(2.055114) q[2];
rz(-2.1975482) q[3];
sx q[3];
rz(-3.0510674) q[3];
sx q[3];
rz(2.0199147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.1679967) q[0];
sx q[0];
rz(-0.43742988) q[0];
sx q[0];
rz(-0.22015372) q[0];
rz(-1.6379697) q[1];
sx q[1];
rz(-1.817768) q[1];
sx q[1];
rz(2.5880609) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65992113) q[0];
sx q[0];
rz(-1.685647) q[0];
sx q[0];
rz(-0.19821313) q[0];
x q[1];
rz(0.4398063) q[2];
sx q[2];
rz(-0.39798073) q[2];
sx q[2];
rz(0.27962886) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.625861) q[1];
sx q[1];
rz(-1.5222094) q[1];
sx q[1];
rz(2.4481985) q[1];
x q[2];
rz(-2.8643478) q[3];
sx q[3];
rz(-2.436815) q[3];
sx q[3];
rz(0.98990209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0963875) q[2];
sx q[2];
rz(-1.4676981) q[2];
sx q[2];
rz(-0.20827797) q[2];
rz(1.1194718) q[3];
sx q[3];
rz(-1.2866311) q[3];
sx q[3];
rz(-2.313405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49523062) q[0];
sx q[0];
rz(-1.7532852) q[0];
sx q[0];
rz(-2.0679423) q[0];
rz(-0.13058361) q[1];
sx q[1];
rz(-1.5740296) q[1];
sx q[1];
rz(-2.1317587) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45859858) q[0];
sx q[0];
rz(-1.3855591) q[0];
sx q[0];
rz(-2.7390929) q[0];
rz(-0.14839852) q[2];
sx q[2];
rz(-1.6434533) q[2];
sx q[2];
rz(1.0277445) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.67391318) q[1];
sx q[1];
rz(-0.99930489) q[1];
sx q[1];
rz(-1.4614146) q[1];
rz(0.66422446) q[3];
sx q[3];
rz(-1.8598607) q[3];
sx q[3];
rz(-0.57830304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4449571) q[2];
sx q[2];
rz(-2.770165) q[2];
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
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0167639) q[0];
sx q[0];
rz(-2.6478196) q[0];
sx q[0];
rz(2.0523409) q[0];
rz(-2.543154) q[1];
sx q[1];
rz(-1.2804223) q[1];
sx q[1];
rz(2.3425897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5157817) q[0];
sx q[0];
rz(-0.41951734) q[0];
sx q[0];
rz(-0.74504344) q[0];
x q[1];
rz(-0.90461055) q[2];
sx q[2];
rz(-1.6450201) q[2];
sx q[2];
rz(-0.42236172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0388698) q[1];
sx q[1];
rz(-2.8029222) q[1];
sx q[1];
rz(-0.46249696) q[1];
rz(0.022074583) q[3];
sx q[3];
rz(-0.24188731) q[3];
sx q[3];
rz(-0.31705184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66990718) q[2];
sx q[2];
rz(-0.50614637) q[2];
sx q[2];
rz(1.8890107) q[2];
rz(1.0505098) q[3];
sx q[3];
rz(-0.59058467) q[3];
sx q[3];
rz(-0.93241507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893386) q[0];
sx q[0];
rz(-1.4006571) q[0];
sx q[0];
rz(-0.5156714) q[0];
rz(1.6853261) q[1];
sx q[1];
rz(-1.3412424) q[1];
sx q[1];
rz(-0.32404831) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7406612) q[0];
sx q[0];
rz(-1.9687511) q[0];
sx q[0];
rz(-1.2952572) q[0];
x q[1];
rz(-0.91941763) q[2];
sx q[2];
rz(-1.9371261) q[2];
sx q[2];
rz(-1.1446143) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0258603) q[1];
sx q[1];
rz(-0.46231356) q[1];
sx q[1];
rz(-1.2532534) q[1];
rz(-pi) q[2];
rz(1.4435931) q[3];
sx q[3];
rz(-0.56131828) q[3];
sx q[3];
rz(1.7030681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9825762) q[2];
sx q[2];
rz(-1.1218718) q[2];
sx q[2];
rz(-2.4328361) q[2];
rz(0.7835663) q[3];
sx q[3];
rz(-1.9400027) q[3];
sx q[3];
rz(1.6243352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5358955) q[0];
sx q[0];
rz(-2.504183) q[0];
sx q[0];
rz(0.15644431) q[0];
rz(-0.63977301) q[1];
sx q[1];
rz(-0.6548869) q[1];
sx q[1];
rz(0.99008647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6789843) q[0];
sx q[0];
rz(-0.96742523) q[0];
sx q[0];
rz(2.4067448) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2210571) q[2];
sx q[2];
rz(-1.4863401) q[2];
sx q[2];
rz(-1.8366369) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2459979) q[1];
sx q[1];
rz(-0.57064547) q[1];
sx q[1];
rz(-2.8938608) q[1];
x q[2];
rz(2.5699432) q[3];
sx q[3];
rz(-2.902973) q[3];
sx q[3];
rz(-2.1941136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5618374) q[2];
sx q[2];
rz(-1.729915) q[2];
sx q[2];
rz(0.61335316) q[2];
rz(-0.26398811) q[3];
sx q[3];
rz(-1.3365021) q[3];
sx q[3];
rz(0.67155513) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092125208) q[0];
sx q[0];
rz(-1.6138074) q[0];
sx q[0];
rz(-2.001979) q[0];
rz(-1.6773979) q[1];
sx q[1];
rz(-2.4212227) q[1];
sx q[1];
rz(-1.5705241) q[1];
rz(0.40205545) q[2];
sx q[2];
rz(-1.9378312) q[2];
sx q[2];
rz(2.1294049) q[2];
rz(-1.5313374) q[3];
sx q[3];
rz(-2.1891371) q[3];
sx q[3];
rz(-2.7923765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
