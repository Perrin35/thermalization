OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5767515) q[0];
sx q[0];
rz(-0.54902005) q[0];
sx q[0];
rz(0.18000552) q[0];
rz(4.0408673) q[1];
sx q[1];
rz(4.0087357) q[1];
sx q[1];
rz(8.0198159) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1464935) q[0];
sx q[0];
rz(-1.4425264) q[0];
sx q[0];
rz(1.487182) q[0];
x q[1];
rz(0.55721941) q[2];
sx q[2];
rz(-1.7168593) q[2];
sx q[2];
rz(2.5596325) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9145292) q[1];
sx q[1];
rz(-1.2458015) q[1];
sx q[1];
rz(1.090828) q[1];
rz(-pi) q[2];
rz(-2.6536921) q[3];
sx q[3];
rz(-1.8630233) q[3];
sx q[3];
rz(3.0350181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0286502) q[2];
sx q[2];
rz(-2.3990227) q[2];
sx q[2];
rz(-2.5600625) q[2];
rz(-2.2384426) q[3];
sx q[3];
rz(-1.6553469) q[3];
sx q[3];
rz(-2.7744897) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0105932) q[0];
sx q[0];
rz(-1.7342664) q[0];
sx q[0];
rz(1.1372239) q[0];
rz(-1.3251023) q[1];
sx q[1];
rz(-0.66665998) q[1];
sx q[1];
rz(-2.4773662) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8839845) q[0];
sx q[0];
rz(-0.31766787) q[0];
sx q[0];
rz(2.5619216) q[0];
x q[1];
rz(-0.082809049) q[2];
sx q[2];
rz(-2.2958404) q[2];
sx q[2];
rz(-2.3963181) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8926516) q[1];
sx q[1];
rz(-2.5499857) q[1];
sx q[1];
rz(-0.68343648) q[1];
rz(-1.799753) q[3];
sx q[3];
rz(-1.4959644) q[3];
sx q[3];
rz(-0.018939806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60435158) q[2];
sx q[2];
rz(-1.1553355) q[2];
sx q[2];
rz(-1.7604766) q[2];
rz(0.85754496) q[3];
sx q[3];
rz(-0.92567912) q[3];
sx q[3];
rz(0.05923567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9581167) q[0];
sx q[0];
rz(-0.97267946) q[0];
sx q[0];
rz(-2.3532975) q[0];
rz(-2.6745785) q[1];
sx q[1];
rz(-2.4772418) q[1];
sx q[1];
rz(2.3036387) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88400052) q[0];
sx q[0];
rz(-1.1282451) q[0];
sx q[0];
rz(-1.6167859) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5068717) q[2];
sx q[2];
rz(-0.41213671) q[2];
sx q[2];
rz(-1.2391547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.26688901) q[1];
sx q[1];
rz(-1.500644) q[1];
sx q[1];
rz(-1.8535437) q[1];
rz(-pi) q[2];
rz(-0.19239088) q[3];
sx q[3];
rz(-2.8225401) q[3];
sx q[3];
rz(2.8448679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9517842) q[2];
sx q[2];
rz(-0.34056792) q[2];
sx q[2];
rz(-1.6050485) q[2];
rz(-0.32478452) q[3];
sx q[3];
rz(-1.6580509) q[3];
sx q[3];
rz(-1.270208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.558641) q[0];
sx q[0];
rz(-2.1138209) q[0];
sx q[0];
rz(0.047274832) q[0];
rz(-1.6167697) q[1];
sx q[1];
rz(-0.44175092) q[1];
sx q[1];
rz(1.5025274) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0086347) q[0];
sx q[0];
rz(-2.0300976) q[0];
sx q[0];
rz(1.2643705) q[0];
rz(-0.45452228) q[2];
sx q[2];
rz(-1.3873867) q[2];
sx q[2];
rz(-1.3647788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33519519) q[1];
sx q[1];
rz(-1.3235705) q[1];
sx q[1];
rz(0.25144318) q[1];
x q[2];
rz(1.371869) q[3];
sx q[3];
rz(-1.891948) q[3];
sx q[3];
rz(-2.7407848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9105685) q[2];
sx q[2];
rz(-2.7474521) q[2];
sx q[2];
rz(-2.2310889) q[2];
rz(-2.2855811) q[3];
sx q[3];
rz(-1.4166219) q[3];
sx q[3];
rz(-0.12106171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0360134) q[0];
sx q[0];
rz(-1.3904904) q[0];
sx q[0];
rz(-0.068583071) q[0];
rz(-2.4225281) q[1];
sx q[1];
rz(-1.9870575) q[1];
sx q[1];
rz(-0.4206492) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4040572) q[0];
sx q[0];
rz(-2.1348663) q[0];
sx q[0];
rz(2.2449136) q[0];
x q[1];
rz(0.22847036) q[2];
sx q[2];
rz(-1.4632311) q[2];
sx q[2];
rz(-3.1247849) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9213564) q[1];
sx q[1];
rz(-1.6097798) q[1];
sx q[1];
rz(0.38631744) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9274022) q[3];
sx q[3];
rz(-0.95396368) q[3];
sx q[3];
rz(-2.9717556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61578304) q[2];
sx q[2];
rz(-2.1685648) q[2];
sx q[2];
rz(0.77965492) q[2];
rz(-0.097213216) q[3];
sx q[3];
rz(-1.8153056) q[3];
sx q[3];
rz(-1.9450845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6355316) q[0];
sx q[0];
rz(-2.4748635) q[0];
sx q[0];
rz(0.42688236) q[0];
rz(-1.4510669) q[1];
sx q[1];
rz(-2.0207113) q[1];
sx q[1];
rz(0.60447398) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6647897) q[0];
sx q[0];
rz(-1.2493629) q[0];
sx q[0];
rz(0.037517083) q[0];
rz(2.8858917) q[2];
sx q[2];
rz(-1.8468282) q[2];
sx q[2];
rz(-0.42975858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6410445) q[1];
sx q[1];
rz(-1.5045847) q[1];
sx q[1];
rz(-1.3891298) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5786798) q[3];
sx q[3];
rz(-1.3775702) q[3];
sx q[3];
rz(0.67487992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2957067) q[2];
sx q[2];
rz(-0.87575951) q[2];
sx q[2];
rz(-0.24629822) q[2];
rz(-1.0036428) q[3];
sx q[3];
rz(-1.6946038) q[3];
sx q[3];
rz(3.057737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2698782) q[0];
sx q[0];
rz(-0.074967472) q[0];
sx q[0];
rz(-2.9366034) q[0];
rz(2.9221453) q[1];
sx q[1];
rz(-1.4402025) q[1];
sx q[1];
rz(-0.27935371) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3059275) q[0];
sx q[0];
rz(-0.45991746) q[0];
sx q[0];
rz(-1.826836) q[0];
rz(3.0952697) q[2];
sx q[2];
rz(-1.8221107) q[2];
sx q[2];
rz(0.77071079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1336891) q[1];
sx q[1];
rz(-0.55396307) q[1];
sx q[1];
rz(-3.0336607) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8362813) q[3];
sx q[3];
rz(-1.1141014) q[3];
sx q[3];
rz(-1.8165464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4651935) q[2];
sx q[2];
rz(-2.3602844) q[2];
sx q[2];
rz(0.87731963) q[2];
rz(-1.8026132) q[3];
sx q[3];
rz(-0.78290144) q[3];
sx q[3];
rz(1.7507929) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1124697) q[0];
sx q[0];
rz(-2.5482197) q[0];
sx q[0];
rz(0.33475885) q[0];
rz(0.33084694) q[1];
sx q[1];
rz(-2.0956764) q[1];
sx q[1];
rz(2.5312993) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90972483) q[0];
sx q[0];
rz(-1.8893161) q[0];
sx q[0];
rz(-0.18778778) q[0];
rz(-pi) q[1];
rz(-1.1534821) q[2];
sx q[2];
rz(-1.8886107) q[2];
sx q[2];
rz(1.870188) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5308633) q[1];
sx q[1];
rz(-1.7708774) q[1];
sx q[1];
rz(-0.51053534) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3819954) q[3];
sx q[3];
rz(-2.5196893) q[3];
sx q[3];
rz(1.7523868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.200013) q[2];
sx q[2];
rz(-2.104685) q[2];
sx q[2];
rz(-1.3882136) q[2];
rz(-2.9774104) q[3];
sx q[3];
rz(-2.1037585) q[3];
sx q[3];
rz(-0.29844704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91069094) q[0];
sx q[0];
rz(-2.0257484) q[0];
sx q[0];
rz(0.13634613) q[0];
rz(-0.048010437) q[1];
sx q[1];
rz(-1.0142356) q[1];
sx q[1];
rz(1.8528574) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.890402) q[0];
sx q[0];
rz(-1.8726761) q[0];
sx q[0];
rz(0.22731486) q[0];
rz(-pi) q[1];
rz(2.9057755) q[2];
sx q[2];
rz(-1.1305446) q[2];
sx q[2];
rz(2.4503675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2338841) q[1];
sx q[1];
rz(-1.4310992) q[1];
sx q[1];
rz(2.8833792) q[1];
x q[2];
rz(0.39548042) q[3];
sx q[3];
rz(-0.90944511) q[3];
sx q[3];
rz(-1.9868074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62782225) q[2];
sx q[2];
rz(-3.1352391) q[2];
sx q[2];
rz(-2.9610236) q[2];
rz(1.1394966) q[3];
sx q[3];
rz(-1.5151016) q[3];
sx q[3];
rz(-3.0212121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1861495) q[0];
sx q[0];
rz(-0.97301617) q[0];
sx q[0];
rz(1.0261616) q[0];
rz(-1.8852662) q[1];
sx q[1];
rz(-1.8812814) q[1];
sx q[1];
rz(0.06591448) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8034466) q[0];
sx q[0];
rz(-2.1915276) q[0];
sx q[0];
rz(0.82671637) q[0];
rz(-pi) q[1];
rz(1.3990551) q[2];
sx q[2];
rz(-1.8800003) q[2];
sx q[2];
rz(-1.1698674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.38635284) q[1];
sx q[1];
rz(-0.59659472) q[1];
sx q[1];
rz(0.42930023) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7700636) q[3];
sx q[3];
rz(-2.5209628) q[3];
sx q[3];
rz(0.21658235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5790448) q[2];
sx q[2];
rz(-2.0936091) q[2];
sx q[2];
rz(-0.5272131) q[2];
rz(-0.15050091) q[3];
sx q[3];
rz(-0.42948693) q[3];
sx q[3];
rz(-2.785717) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81446205) q[0];
sx q[0];
rz(-2.3476063) q[0];
sx q[0];
rz(0.052477947) q[0];
rz(1.3606701) q[1];
sx q[1];
rz(-2.4130029) q[1];
sx q[1];
rz(-1.3389814) q[1];
rz(1.1293148) q[2];
sx q[2];
rz(-1.8649351) q[2];
sx q[2];
rz(-1.9492016) q[2];
rz(-0.99450022) q[3];
sx q[3];
rz(-2.4296843) q[3];
sx q[3];
rz(-1.4083023) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
