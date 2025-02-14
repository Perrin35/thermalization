OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.225086) q[0];
sx q[0];
rz(-0.14296159) q[0];
sx q[0];
rz(1.6530871) q[0];
rz(-2.0090964) q[1];
sx q[1];
rz(5.8512591) q[1];
sx q[1];
rz(6.9195256) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37381946) q[0];
sx q[0];
rz(-0.99080201) q[0];
sx q[0];
rz(-2.679002) q[0];
x q[1];
rz(2.752334) q[2];
sx q[2];
rz(-1.2137652) q[2];
sx q[2];
rz(0.87212002) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.26633599) q[1];
sx q[1];
rz(-2.7260029) q[1];
sx q[1];
rz(1.8824091) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2872055) q[3];
sx q[3];
rz(-0.99528377) q[3];
sx q[3];
rz(0.42785409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1119614) q[2];
sx q[2];
rz(-2.0165069) q[2];
sx q[2];
rz(-2.4911528) q[2];
rz(1.8247617) q[3];
sx q[3];
rz(-1.7951169) q[3];
sx q[3];
rz(-0.10183798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0541075) q[0];
sx q[0];
rz(-2.5547145) q[0];
sx q[0];
rz(-0.45463872) q[0];
rz(3.1184323) q[1];
sx q[1];
rz(-1.4440447) q[1];
sx q[1];
rz(0.32618162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5989203) q[0];
sx q[0];
rz(-2.4389451) q[0];
sx q[0];
rz(-0.62851466) q[0];
rz(-pi) q[1];
rz(-1.9942547) q[2];
sx q[2];
rz(-1.2988161) q[2];
sx q[2];
rz(0.47951298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7190951) q[1];
sx q[1];
rz(-2.0764253) q[1];
sx q[1];
rz(-1.0984332) q[1];
rz(-pi) q[2];
rz(1.7214016) q[3];
sx q[3];
rz(-2.4753627) q[3];
sx q[3];
rz(-0.008994313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0252016) q[2];
sx q[2];
rz(-1.9390743) q[2];
sx q[2];
rz(-0.95428673) q[2];
rz(-1.8918234) q[3];
sx q[3];
rz(-0.3229177) q[3];
sx q[3];
rz(-0.0013110411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8702451) q[0];
sx q[0];
rz(-1.988669) q[0];
sx q[0];
rz(1.1767607) q[0];
rz(-2.7765043) q[1];
sx q[1];
rz(-1.3026214) q[1];
sx q[1];
rz(1.5286068) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1592401) q[0];
sx q[0];
rz(-1.591658) q[0];
sx q[0];
rz(3.1283035) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0143378) q[2];
sx q[2];
rz(-1.6422286) q[2];
sx q[2];
rz(-2.254666) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0658256) q[1];
sx q[1];
rz(-1.6012553) q[1];
sx q[1];
rz(1.7225313) q[1];
rz(-pi) q[2];
rz(-1.6284451) q[3];
sx q[3];
rz(-1.5022657) q[3];
sx q[3];
rz(1.203361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8927346) q[2];
sx q[2];
rz(-2.5881519) q[2];
sx q[2];
rz(-1.152732) q[2];
rz(-1.24498) q[3];
sx q[3];
rz(-2.1608519) q[3];
sx q[3];
rz(2.8583756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98274851) q[0];
sx q[0];
rz(-11*pi/12) q[0];
sx q[0];
rz(-2.4520279) q[0];
rz(-0.025029643) q[1];
sx q[1];
rz(-1.6233147) q[1];
sx q[1];
rz(1.7207346) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11224225) q[0];
sx q[0];
rz(-2.1922702) q[0];
sx q[0];
rz(1.6981324) q[0];
rz(-pi) q[1];
x q[1];
rz(2.000932) q[2];
sx q[2];
rz(-1.4073512) q[2];
sx q[2];
rz(2.8317833) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5031227) q[1];
sx q[1];
rz(-1.5431738) q[1];
sx q[1];
rz(-2.6383328) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6715325) q[3];
sx q[3];
rz(-0.42123171) q[3];
sx q[3];
rz(-2.4249113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20874061) q[2];
sx q[2];
rz(-0.17540652) q[2];
sx q[2];
rz(-1.9746732) q[2];
rz(-2.6973727) q[3];
sx q[3];
rz(-1.9053713) q[3];
sx q[3];
rz(-2.8319632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8532448) q[0];
sx q[0];
rz(-1.3695559) q[0];
sx q[0];
rz(0.41611588) q[0];
rz(-0.5101282) q[1];
sx q[1];
rz(-0.77113873) q[1];
sx q[1];
rz(2.3663734) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3907523) q[0];
sx q[0];
rz(-1.7699983) q[0];
sx q[0];
rz(-2.5222523) q[0];
x q[1];
rz(0.40480308) q[2];
sx q[2];
rz(-2.7049613) q[2];
sx q[2];
rz(-2.8711808) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8000335) q[1];
sx q[1];
rz(-1.4455035) q[1];
sx q[1];
rz(3.028246) q[1];
rz(-1.1936504) q[3];
sx q[3];
rz(-2.1324649) q[3];
sx q[3];
rz(-2.9557754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2148718) q[2];
sx q[2];
rz(-0.99433172) q[2];
sx q[2];
rz(-2.1057687) q[2];
rz(0.18946798) q[3];
sx q[3];
rz(-2.6697956) q[3];
sx q[3];
rz(-0.50260472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8871317) q[0];
sx q[0];
rz(-0.04763617) q[0];
sx q[0];
rz(2.7857842) q[0];
rz(1.4954781) q[1];
sx q[1];
rz(-0.9551841) q[1];
sx q[1];
rz(-1.1411508) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1352859) q[0];
sx q[0];
rz(-1.7809488) q[0];
sx q[0];
rz(0.89583136) q[0];
rz(-pi) q[1];
rz(-0.95924033) q[2];
sx q[2];
rz(-1.578176) q[2];
sx q[2];
rz(-2.6098698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6427152) q[1];
sx q[1];
rz(-1.4659473) q[1];
sx q[1];
rz(3.0003606) q[1];
rz(-3.1294501) q[3];
sx q[3];
rz(-1.1369161) q[3];
sx q[3];
rz(1.7858684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5100539) q[2];
sx q[2];
rz(-0.71439356) q[2];
sx q[2];
rz(1.7047403) q[2];
rz(-2.0884183) q[3];
sx q[3];
rz(-1.5318233) q[3];
sx q[3];
rz(-2.6385782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.961504) q[0];
sx q[0];
rz(-0.29619521) q[0];
sx q[0];
rz(0.12251138) q[0];
rz(-0.65613121) q[1];
sx q[1];
rz(-1.2686814) q[1];
sx q[1];
rz(1.3414541) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91514523) q[0];
sx q[0];
rz(-1.5350123) q[0];
sx q[0];
rz(-0.045374845) q[0];
rz(-pi) q[1];
rz(0.28320988) q[2];
sx q[2];
rz(-1.4669751) q[2];
sx q[2];
rz(1.6714754) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0602051) q[1];
sx q[1];
rz(-1.3045746) q[1];
sx q[1];
rz(-0.93764825) q[1];
rz(-pi) q[2];
rz(-0.44388598) q[3];
sx q[3];
rz(-1.6735014) q[3];
sx q[3];
rz(3.1004124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5682257) q[2];
sx q[2];
rz(-1.919701) q[2];
sx q[2];
rz(-1.4637671) q[2];
rz(-2.407414) q[3];
sx q[3];
rz(-0.28407431) q[3];
sx q[3];
rz(-1.3045788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-1.3463335) q[0];
sx q[0];
rz(-1.7743552) q[0];
sx q[0];
rz(-0.45860589) q[0];
rz(-2.1268225) q[1];
sx q[1];
rz(-1.9472803) q[1];
sx q[1];
rz(-1.0626622) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39358562) q[0];
sx q[0];
rz(-2.9655122) q[0];
sx q[0];
rz(2.8066471) q[0];
rz(-2.3331679) q[2];
sx q[2];
rz(-1.9701517) q[2];
sx q[2];
rz(-2.6724919) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.96982376) q[1];
sx q[1];
rz(-0.14482982) q[1];
sx q[1];
rz(2.0098196) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0513145) q[3];
sx q[3];
rz(-1.8534582) q[3];
sx q[3];
rz(-1.1804898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8899272) q[2];
sx q[2];
rz(-1.7532316) q[2];
sx q[2];
rz(1.7903222) q[2];
rz(1.9225559) q[3];
sx q[3];
rz(-0.42870298) q[3];
sx q[3];
rz(-0.8333227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625921) q[0];
sx q[0];
rz(-1.3840249) q[0];
sx q[0];
rz(-2.2341527) q[0];
rz(2.9145248) q[1];
sx q[1];
rz(-1.0457958) q[1];
sx q[1];
rz(-2.2423832) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.718841) q[0];
sx q[0];
rz(-2.1678574) q[0];
sx q[0];
rz(2.5121252) q[0];
x q[1];
rz(-0.078446968) q[2];
sx q[2];
rz(-1.6529993) q[2];
sx q[2];
rz(-1.3644827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6688706) q[1];
sx q[1];
rz(-2.5963915) q[1];
sx q[1];
rz(1.9342058) q[1];
rz(-pi) q[2];
rz(1.6257004) q[3];
sx q[3];
rz(-0.74885741) q[3];
sx q[3];
rz(-2.1532173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9565309) q[2];
sx q[2];
rz(-1.9756292) q[2];
sx q[2];
rz(2.5409307) q[2];
rz(-1.0968084) q[3];
sx q[3];
rz(-0.59022248) q[3];
sx q[3];
rz(-2.2685331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7507062) q[0];
sx q[0];
rz(-2.0426671) q[0];
sx q[0];
rz(2.1571958) q[0];
rz(0.4920494) q[1];
sx q[1];
rz(-1.7285408) q[1];
sx q[1];
rz(-1.9459928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22532428) q[0];
sx q[0];
rz(-1.9231503) q[0];
sx q[0];
rz(0.63017643) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1629647) q[2];
sx q[2];
rz(-2.371863) q[2];
sx q[2];
rz(-2.9362188) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3477563) q[1];
sx q[1];
rz(-1.4897921) q[1];
sx q[1];
rz(1.9543314) q[1];
rz(-pi) q[2];
rz(2.0298084) q[3];
sx q[3];
rz(-0.85235607) q[3];
sx q[3];
rz(1.8000923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0857346) q[2];
sx q[2];
rz(-1.5208289) q[2];
sx q[2];
rz(-0.98304191) q[2];
rz(2.8899657) q[3];
sx q[3];
rz(-2.7139137) q[3];
sx q[3];
rz(-1.0035286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.49878237) q[0];
sx q[0];
rz(-0.97220535) q[0];
sx q[0];
rz(-2.4844949) q[0];
rz(1.2596754) q[1];
sx q[1];
rz(-1.4114264) q[1];
sx q[1];
rz(-2.2687601) q[1];
rz(2.2839584) q[2];
sx q[2];
rz(-1.7339755) q[2];
sx q[2];
rz(-0.82790464) q[2];
rz(-0.33369251) q[3];
sx q[3];
rz(-1.8751388) q[3];
sx q[3];
rz(1.4061385) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
