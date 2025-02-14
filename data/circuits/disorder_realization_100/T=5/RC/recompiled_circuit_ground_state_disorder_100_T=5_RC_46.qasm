OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72803175) q[0];
sx q[0];
rz(-0.95946884) q[0];
sx q[0];
rz(2.5897107) q[0];
rz(-0.38504398) q[1];
sx q[1];
rz(4.9209891) q[1];
sx q[1];
rz(9.8434386) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075567632) q[0];
sx q[0];
rz(-2.2097144) q[0];
sx q[0];
rz(-0.35982168) q[0];
rz(-2.6919882) q[2];
sx q[2];
rz(-1.7309963) q[2];
sx q[2];
rz(-0.9148324) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4837515) q[1];
sx q[1];
rz(-2.9194909) q[1];
sx q[1];
rz(2.7276843) q[1];
x q[2];
rz(3.1244784) q[3];
sx q[3];
rz(-1.601222) q[3];
sx q[3];
rz(-0.78998128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5358413) q[2];
sx q[2];
rz(-1.8391106) q[2];
sx q[2];
rz(0.38899404) q[2];
rz(-0.89748663) q[3];
sx q[3];
rz(-0.56098452) q[3];
sx q[3];
rz(1.8628023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9101343) q[0];
sx q[0];
rz(-2.1899905) q[0];
sx q[0];
rz(2.8125473) q[0];
rz(-1.7973876) q[1];
sx q[1];
rz(-2.3842594) q[1];
sx q[1];
rz(-3.0175041) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24606516) q[0];
sx q[0];
rz(-2.82282) q[0];
sx q[0];
rz(0.84455873) q[0];
rz(2.8555238) q[2];
sx q[2];
rz(-1.7033615) q[2];
sx q[2];
rz(2.2372467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6581003) q[1];
sx q[1];
rz(-2.0112733) q[1];
sx q[1];
rz(-0.3962724) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2749008) q[3];
sx q[3];
rz(-1.4128437) q[3];
sx q[3];
rz(-2.7622836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.2825534) q[2];
sx q[2];
rz(-0.48290792) q[2];
sx q[2];
rz(2.5340951) q[2];
rz(0.88827682) q[3];
sx q[3];
rz(-1.6162623) q[3];
sx q[3];
rz(-1.4265149) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68179503) q[0];
sx q[0];
rz(-1.8495704) q[0];
sx q[0];
rz(-0.23455308) q[0];
rz(-2.081743) q[1];
sx q[1];
rz(-1.1666965) q[1];
sx q[1];
rz(0.92811981) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7674196) q[0];
sx q[0];
rz(-2.1819127) q[0];
sx q[0];
rz(-1.060063) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6675148) q[2];
sx q[2];
rz(-1.9454207) q[2];
sx q[2];
rz(-1.4371536) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27552128) q[1];
sx q[1];
rz(-1.2433941) q[1];
sx q[1];
rz(-1.7129461) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2006525) q[3];
sx q[3];
rz(-1.4273774) q[3];
sx q[3];
rz(-0.36597825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5633391) q[2];
sx q[2];
rz(-1.4207062) q[2];
sx q[2];
rz(-1.2467747) q[2];
rz(-0.16683821) q[3];
sx q[3];
rz(-2.2127547) q[3];
sx q[3];
rz(1.5290574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.700915) q[0];
sx q[0];
rz(-3.0755141) q[0];
sx q[0];
rz(0.75827688) q[0];
rz(1.5757163) q[1];
sx q[1];
rz(-2.5570452) q[1];
sx q[1];
rz(2.27104) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2038531) q[0];
sx q[0];
rz(-2.014262) q[0];
sx q[0];
rz(-0.70758836) q[0];
rz(-pi) q[1];
rz(-2.4192018) q[2];
sx q[2];
rz(-2.4263627) q[2];
sx q[2];
rz(1.720429) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4728001) q[1];
sx q[1];
rz(-1.973098) q[1];
sx q[1];
rz(0.16568664) q[1];
rz(3.0421889) q[3];
sx q[3];
rz(-0.63025852) q[3];
sx q[3];
rz(0.46168345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7190651) q[2];
sx q[2];
rz(-1.0079404) q[2];
sx q[2];
rz(0.75971216) q[2];
rz(2.4921913) q[3];
sx q[3];
rz(-1.5348744) q[3];
sx q[3];
rz(-1.5627129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0287057) q[0];
sx q[0];
rz(-2.4876471) q[0];
sx q[0];
rz(1.1829859) q[0];
rz(0.68823367) q[1];
sx q[1];
rz(-1.8249244) q[1];
sx q[1];
rz(2.7722955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0770586) q[0];
sx q[0];
rz(-2.4774533) q[0];
sx q[0];
rz(-0.044483552) q[0];
rz(-pi) q[1];
rz(-2.4299116) q[2];
sx q[2];
rz(-1.7717012) q[2];
sx q[2];
rz(-0.37362675) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5253882) q[1];
sx q[1];
rz(-1.3862351) q[1];
sx q[1];
rz(2.2730458) q[1];
rz(-pi) q[2];
rz(2.8083388) q[3];
sx q[3];
rz(-1.8849533) q[3];
sx q[3];
rz(-2.9314048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1987622) q[2];
sx q[2];
rz(-1.3254712) q[2];
sx q[2];
rz(-0.99786264) q[2];
rz(-2.5241847) q[3];
sx q[3];
rz(-0.95881763) q[3];
sx q[3];
rz(2.3203885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038641039) q[0];
sx q[0];
rz(-1.2741673) q[0];
sx q[0];
rz(1.8983023) q[0];
rz(-0.78168166) q[1];
sx q[1];
rz(-1.8676753) q[1];
sx q[1];
rz(0.26085687) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7383682) q[0];
sx q[0];
rz(-2.0595831) q[0];
sx q[0];
rz(1.3317778) q[0];
rz(-pi) q[1];
rz(-2.2736808) q[2];
sx q[2];
rz(-2.5943668) q[2];
sx q[2];
rz(0.99582129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1697537) q[1];
sx q[1];
rz(-1.1128281) q[1];
sx q[1];
rz(-2.9724246) q[1];
rz(2.9294037) q[3];
sx q[3];
rz(-2.1916323) q[3];
sx q[3];
rz(2.5356136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6692052) q[2];
sx q[2];
rz(-2.278625) q[2];
sx q[2];
rz(2.6589987) q[2];
rz(-3.0274296) q[3];
sx q[3];
rz(-2.0518905) q[3];
sx q[3];
rz(0.94858661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77940762) q[0];
sx q[0];
rz(-0.055883378) q[0];
sx q[0];
rz(-2.8756397) q[0];
rz(-2.7159122) q[1];
sx q[1];
rz(-1.2454147) q[1];
sx q[1];
rz(2.6735305) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6857306) q[0];
sx q[0];
rz(-1.0294559) q[0];
sx q[0];
rz(3.1336407) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8276026) q[2];
sx q[2];
rz(-2.7903284) q[2];
sx q[2];
rz(-2.6333269) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1359766) q[1];
sx q[1];
rz(-1.3141201) q[1];
sx q[1];
rz(1.0631494) q[1];
rz(-pi) q[2];
rz(0.031074957) q[3];
sx q[3];
rz(-2.090775) q[3];
sx q[3];
rz(1.7620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.257306) q[2];
sx q[2];
rz(-1.2904737) q[2];
sx q[2];
rz(2.5878944) q[2];
rz(2.2181559) q[3];
sx q[3];
rz(-0.90673509) q[3];
sx q[3];
rz(-3.0807909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4629352) q[0];
sx q[0];
rz(-2.904992) q[0];
sx q[0];
rz(2.4635354) q[0];
rz(0.97738114) q[1];
sx q[1];
rz(-1.9934318) q[1];
sx q[1];
rz(-1.4378907) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8554094) q[0];
sx q[0];
rz(-1.1146422) q[0];
sx q[0];
rz(0.21502226) q[0];
x q[1];
rz(3.1372141) q[2];
sx q[2];
rz(-0.56398773) q[2];
sx q[2];
rz(-1.1396688) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9919093) q[1];
sx q[1];
rz(-1.8018197) q[1];
sx q[1];
rz(-0.2607338) q[1];
rz(-2.1071042) q[3];
sx q[3];
rz(-0.89867979) q[3];
sx q[3];
rz(1.514707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4108654) q[2];
sx q[2];
rz(-2.5812456) q[2];
sx q[2];
rz(-0.66696683) q[2];
rz(-0.11735958) q[3];
sx q[3];
rz(-1.8662235) q[3];
sx q[3];
rz(-1.7608775) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9475107) q[0];
sx q[0];
rz(-1.5912594) q[0];
sx q[0];
rz(-2.7981753) q[0];
rz(-1.4944448) q[1];
sx q[1];
rz(-0.98973715) q[1];
sx q[1];
rz(-0.48430482) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5268516) q[0];
sx q[0];
rz(-2.1526045) q[0];
sx q[0];
rz(-1.2528166) q[0];
rz(-pi) q[1];
rz(-0.76700751) q[2];
sx q[2];
rz(-2.5756774) q[2];
sx q[2];
rz(-0.2919251) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2769988) q[1];
sx q[1];
rz(-2.5127257) q[1];
sx q[1];
rz(0.38284812) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8922594) q[3];
sx q[3];
rz(-1.7847381) q[3];
sx q[3];
rz(0.46505022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9453498) q[2];
sx q[2];
rz(-2.0413155) q[2];
sx q[2];
rz(-0.84570447) q[2];
rz(-1.2196994) q[3];
sx q[3];
rz(-1.2518576) q[3];
sx q[3];
rz(0.20488258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4591111) q[0];
sx q[0];
rz(-0.24523188) q[0];
sx q[0];
rz(-0.58897513) q[0];
rz(-2.466195) q[1];
sx q[1];
rz(-0.94869906) q[1];
sx q[1];
rz(1.4640456) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1896473) q[0];
sx q[0];
rz(-2.4485817) q[0];
sx q[0];
rz(-0.89784867) q[0];
rz(-1.6590933) q[2];
sx q[2];
rz(-2.4064734) q[2];
sx q[2];
rz(-2.5500848) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3879536) q[1];
sx q[1];
rz(-1.6170701) q[1];
sx q[1];
rz(-2.5548177) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53873976) q[3];
sx q[3];
rz(-0.77048388) q[3];
sx q[3];
rz(-0.24365261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7318763) q[2];
sx q[2];
rz(-1.3940553) q[2];
sx q[2];
rz(2.0173006) q[2];
rz(2.2935947) q[3];
sx q[3];
rz(-2.336899) q[3];
sx q[3];
rz(3.1112352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59415862) q[0];
sx q[0];
rz(-1.1347329) q[0];
sx q[0];
rz(-1.5190079) q[0];
rz(0.92319725) q[1];
sx q[1];
rz(-1.0293488) q[1];
sx q[1];
rz(1.6346288) q[1];
rz(2.831645) q[2];
sx q[2];
rz(-1.7480231) q[2];
sx q[2];
rz(2.9668948) q[2];
rz(2.2062929) q[3];
sx q[3];
rz(-0.37105303) q[3];
sx q[3];
rz(1.1340352) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
