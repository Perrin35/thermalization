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
rz(-2.5858606) q[0];
sx q[0];
rz(-1.2795804) q[0];
sx q[0];
rz(0.32787856) q[0];
rz(-2.9887587) q[1];
sx q[1];
rz(-2.6522377) q[1];
sx q[1];
rz(2.1305003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0450789) q[0];
sx q[0];
rz(-0.9849087) q[0];
sx q[0];
rz(1.3784598) q[0];
rz(-1.7494124) q[2];
sx q[2];
rz(-1.1682604) q[2];
sx q[2];
rz(2.2312763) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3950306) q[1];
sx q[1];
rz(-0.85828188) q[1];
sx q[1];
rz(-2.0601574) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1527083) q[3];
sx q[3];
rz(-1.3527414) q[3];
sx q[3];
rz(2.5421028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.27935394) q[2];
sx q[2];
rz(-0.78247672) q[2];
sx q[2];
rz(-1.5834825) q[2];
rz(2.8067348) q[3];
sx q[3];
rz(-2.0575276) q[3];
sx q[3];
rz(-2.8607232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8849477) q[0];
sx q[0];
rz(-0.56459752) q[0];
sx q[0];
rz(2.8096492) q[0];
rz(0.360082) q[1];
sx q[1];
rz(-1.8372476) q[1];
sx q[1];
rz(-0.28775451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641787) q[0];
sx q[0];
rz(-2.8197643) q[0];
sx q[0];
rz(-0.68049707) q[0];
x q[1];
rz(0.73189484) q[2];
sx q[2];
rz(-1.8127155) q[2];
sx q[2];
rz(-0.11920028) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8043878) q[1];
sx q[1];
rz(-1.8998977) q[1];
sx q[1];
rz(1.2044524) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23116206) q[3];
sx q[3];
rz(-1.729768) q[3];
sx q[3];
rz(2.8249225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.687279) q[2];
sx q[2];
rz(-1.5328898) q[2];
sx q[2];
rz(1.3909371) q[2];
rz(-2.2823997) q[3];
sx q[3];
rz(-1.8682559) q[3];
sx q[3];
rz(-2.9505762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.745382) q[0];
sx q[0];
rz(-1.6062382) q[0];
sx q[0];
rz(0.23042738) q[0];
rz(2.6241809) q[1];
sx q[1];
rz(-1.0942065) q[1];
sx q[1];
rz(2.8527625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010832615) q[0];
sx q[0];
rz(-1.9790181) q[0];
sx q[0];
rz(-0.064716332) q[0];
rz(-pi) q[1];
rz(0.46858139) q[2];
sx q[2];
rz(-1.6372674) q[2];
sx q[2];
rz(0.85487142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9625712) q[1];
sx q[1];
rz(-1.5177392) q[1];
sx q[1];
rz(1.1607443) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7942811) q[3];
sx q[3];
rz(-1.399125) q[3];
sx q[3];
rz(-0.61644256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1032224) q[2];
sx q[2];
rz(-2.4028845) q[2];
sx q[2];
rz(3.1008516) q[2];
rz(-3.0401958) q[3];
sx q[3];
rz(-1.3359759) q[3];
sx q[3];
rz(-2.0749157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7539702) q[0];
sx q[0];
rz(-0.0059703537) q[0];
sx q[0];
rz(0.23400865) q[0];
rz(0.19800828) q[1];
sx q[1];
rz(-2.104069) q[1];
sx q[1];
rz(0.98349804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9390181) q[0];
sx q[0];
rz(-1.6133623) q[0];
sx q[0];
rz(2.0545511) q[0];
x q[1];
rz(-0.51362546) q[2];
sx q[2];
rz(-1.4911998) q[2];
sx q[2];
rz(0.74515504) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.75960052) q[1];
sx q[1];
rz(-2.1651717) q[1];
sx q[1];
rz(-2.4305953) q[1];
rz(-pi) q[2];
rz(-1.181482) q[3];
sx q[3];
rz(-1.6004171) q[3];
sx q[3];
rz(-2.9148341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6733751) q[2];
sx q[2];
rz(-2.8523291) q[2];
sx q[2];
rz(0.77486983) q[2];
rz(1.3261999) q[3];
sx q[3];
rz(-1.1893136) q[3];
sx q[3];
rz(2.6302122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4020017) q[0];
sx q[0];
rz(-0.87257659) q[0];
sx q[0];
rz(0.032489754) q[0];
rz(1.2292817) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(3.0585739) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5755324) q[0];
sx q[0];
rz(-1.119918) q[0];
sx q[0];
rz(-1.321248) q[0];
rz(0.26814383) q[2];
sx q[2];
rz(-1.4073155) q[2];
sx q[2];
rz(-0.016506052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0159576) q[1];
sx q[1];
rz(-0.86125606) q[1];
sx q[1];
rz(0.63035359) q[1];
rz(-0.015542726) q[3];
sx q[3];
rz(-2.2076026) q[3];
sx q[3];
rz(1.6890845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7897537) q[2];
sx q[2];
rz(-0.80078501) q[2];
sx q[2];
rz(-2.9808673) q[2];
rz(1.1927346) q[3];
sx q[3];
rz(-0.21218097) q[3];
sx q[3];
rz(0.76782697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47650325) q[0];
sx q[0];
rz(-2.022321) q[0];
sx q[0];
rz(0.29092586) q[0];
rz(-1.7581958) q[1];
sx q[1];
rz(-1.29888) q[1];
sx q[1];
rz(-0.49984041) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4668149) q[0];
sx q[0];
rz(-1.5272015) q[0];
sx q[0];
rz(1.6131667) q[0];
x q[1];
rz(0.55576365) q[2];
sx q[2];
rz(-0.92257181) q[2];
sx q[2];
rz(0.52582914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0575057) q[1];
sx q[1];
rz(-1.1228787) q[1];
sx q[1];
rz(-2.491076) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53414102) q[3];
sx q[3];
rz(-2.2075966) q[3];
sx q[3];
rz(2.4780688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2073888) q[2];
sx q[2];
rz(-0.61177212) q[2];
sx q[2];
rz(-2.440051) q[2];
rz(-2.1675341) q[3];
sx q[3];
rz(-2.3249224) q[3];
sx q[3];
rz(0.90132236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3153673) q[0];
sx q[0];
rz(-0.81804818) q[0];
sx q[0];
rz(-2.4819964) q[0];
rz(1.2391799) q[1];
sx q[1];
rz(-2.0678803) q[1];
sx q[1];
rz(-1.0741796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095431134) q[0];
sx q[0];
rz(-3.1218596) q[0];
sx q[0];
rz(0.82206263) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7137382) q[2];
sx q[2];
rz(-1.9249467) q[2];
sx q[2];
rz(-2.4216975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.85107899) q[1];
sx q[1];
rz(-2.004262) q[1];
sx q[1];
rz(2.9032533) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6857008) q[3];
sx q[3];
rz(-2.3237202) q[3];
sx q[3];
rz(2.4867833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2913975) q[2];
sx q[2];
rz(-0.75866282) q[2];
sx q[2];
rz(0.58471739) q[2];
rz(2.0753453) q[3];
sx q[3];
rz(-1.9138347) q[3];
sx q[3];
rz(2.7339981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19121118) q[0];
sx q[0];
rz(-1.9712912) q[0];
sx q[0];
rz(-2.6608652) q[0];
rz(1.9302543) q[1];
sx q[1];
rz(-1.5296661) q[1];
sx q[1];
rz(2.5239677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6473501) q[0];
sx q[0];
rz(-1.2456919) q[0];
sx q[0];
rz(-0.0044429739) q[0];
rz(-2.0938056) q[2];
sx q[2];
rz(-2.1126502) q[2];
sx q[2];
rz(2.7043992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.94719515) q[1];
sx q[1];
rz(-0.36809599) q[1];
sx q[1];
rz(-0.93666623) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9493773) q[3];
sx q[3];
rz(-1.6055067) q[3];
sx q[3];
rz(-2.672003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9628613) q[2];
sx q[2];
rz(-1.3946673) q[2];
sx q[2];
rz(-2.2507131) q[2];
rz(-1.0293055) q[3];
sx q[3];
rz(-1.7512713) q[3];
sx q[3];
rz(1.5911969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7903098) q[0];
sx q[0];
rz(-0.93515486) q[0];
sx q[0];
rz(1.9679605) q[0];
rz(-1.1890746) q[1];
sx q[1];
rz(-0.74780858) q[1];
sx q[1];
rz(-0.48114166) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3279151) q[0];
sx q[0];
rz(-0.60773931) q[0];
sx q[0];
rz(-0.39791664) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16002197) q[2];
sx q[2];
rz(-0.56567398) q[2];
sx q[2];
rz(-0.46505798) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6509241) q[1];
sx q[1];
rz(-1.4238402) q[1];
sx q[1];
rz(0.96687324) q[1];
x q[2];
rz(0.88925006) q[3];
sx q[3];
rz(-0.64683952) q[3];
sx q[3];
rz(-1.3726039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9159307) q[2];
sx q[2];
rz(-1.2292726) q[2];
sx q[2];
rz(-1.6602328) q[2];
rz(0.046317421) q[3];
sx q[3];
rz(-1.9527083) q[3];
sx q[3];
rz(-1.2272629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939821) q[0];
sx q[0];
rz(-1.0337669) q[0];
sx q[0];
rz(0.069742918) q[0];
rz(0.28930411) q[1];
sx q[1];
rz(-1.2846839) q[1];
sx q[1];
rz(1.0522254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3892059) q[0];
sx q[0];
rz(-1.9294039) q[0];
sx q[0];
rz(1.2239271) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2329726) q[2];
sx q[2];
rz(-1.7025456) q[2];
sx q[2];
rz(-0.8140623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3527413) q[1];
sx q[1];
rz(-0.53376275) q[1];
sx q[1];
rz(-2.3351921) q[1];
x q[2];
rz(-2.7131548) q[3];
sx q[3];
rz(-1.1827381) q[3];
sx q[3];
rz(-2.4840499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5436486) q[2];
sx q[2];
rz(-1.9660549) q[2];
sx q[2];
rz(0.0607461) q[2];
rz(2.8525823) q[3];
sx q[3];
rz(-1.776639) q[3];
sx q[3];
rz(-1.8665159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.48988265) q[0];
sx q[0];
rz(-1.4986421) q[0];
sx q[0];
rz(-1.0810252) q[0];
rz(1.0501077) q[1];
sx q[1];
rz(-1.0625912) q[1];
sx q[1];
rz(-1.2407632) q[1];
rz(1.1012668) q[2];
sx q[2];
rz(-1.4360089) q[2];
sx q[2];
rz(0.77965005) q[2];
rz(-1.0084739) q[3];
sx q[3];
rz(-0.49839603) q[3];
sx q[3];
rz(-2.8965542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
