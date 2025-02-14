OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6381792) q[0];
sx q[0];
rz(4.4210202) q[0];
sx q[0];
rz(11.790334) q[0];
rz(-2.4322721) q[1];
sx q[1];
rz(-1.5172989) q[1];
sx q[1];
rz(-2.5425743) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1206643) q[0];
sx q[0];
rz(-2.0428223) q[0];
sx q[0];
rz(3.0042404) q[0];
x q[1];
rz(0.8813192) q[2];
sx q[2];
rz(-1.3506839) q[2];
sx q[2];
rz(2.628919) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4091275) q[1];
sx q[1];
rz(-2.5149087) q[1];
sx q[1];
rz(-2.5123358) q[1];
rz(-pi) q[2];
rz(2.0008068) q[3];
sx q[3];
rz(-2.3035192) q[3];
sx q[3];
rz(0.091146745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0029995) q[2];
sx q[2];
rz(-1.1370167) q[2];
sx q[2];
rz(0.19621672) q[2];
rz(-1.044322) q[3];
sx q[3];
rz(-0.82572562) q[3];
sx q[3];
rz(-2.830982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1021378) q[0];
sx q[0];
rz(-2.4882443) q[0];
sx q[0];
rz(-2.2838604) q[0];
rz(2.6013382) q[1];
sx q[1];
rz(-2.0876355) q[1];
sx q[1];
rz(-0.78261715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5531571) q[0];
sx q[0];
rz(-2.5194476) q[0];
sx q[0];
rz(1.4216656) q[0];
rz(-pi) q[1];
rz(0.23553964) q[2];
sx q[2];
rz(-0.96220926) q[2];
sx q[2];
rz(2.4289301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5021584) q[1];
sx q[1];
rz(-0.63259387) q[1];
sx q[1];
rz(-2.8117287) q[1];
rz(-0.9192809) q[3];
sx q[3];
rz(-0.79943919) q[3];
sx q[3];
rz(-2.8795867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9767849) q[2];
sx q[2];
rz(-2.3471577) q[2];
sx q[2];
rz(2.5505193) q[2];
rz(-1.0278541) q[3];
sx q[3];
rz(-0.63575345) q[3];
sx q[3];
rz(0.13636057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9117821) q[0];
sx q[0];
rz(-2.6214143) q[0];
sx q[0];
rz(2.0853364) q[0];
rz(0.99110574) q[1];
sx q[1];
rz(-1.7678363) q[1];
sx q[1];
rz(-2.3371005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5493889) q[0];
sx q[0];
rz(-2.443586) q[0];
sx q[0];
rz(0.42085676) q[0];
rz(2.1995538) q[2];
sx q[2];
rz(-1.647718) q[2];
sx q[2];
rz(0.49401894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1002297) q[1];
sx q[1];
rz(-1.3303262) q[1];
sx q[1];
rz(2.9036456) q[1];
rz(-2.8896585) q[3];
sx q[3];
rz(-2.478699) q[3];
sx q[3];
rz(-0.54352647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6911917) q[2];
sx q[2];
rz(-0.85019008) q[2];
sx q[2];
rz(2.5737393) q[2];
rz(1.9495226) q[3];
sx q[3];
rz(-2.6046533) q[3];
sx q[3];
rz(-0.15538628) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464226) q[0];
sx q[0];
rz(-1.7429054) q[0];
sx q[0];
rz(1.1871185) q[0];
rz(-2.8861956) q[1];
sx q[1];
rz(-1.4030158) q[1];
sx q[1];
rz(2.752221) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02006836) q[0];
sx q[0];
rz(-0.25857718) q[0];
sx q[0];
rz(2.5917087) q[0];
x q[1];
rz(-0.03002982) q[2];
sx q[2];
rz(-2.2976934) q[2];
sx q[2];
rz(0.50315693) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8592718) q[1];
sx q[1];
rz(-0.92073694) q[1];
sx q[1];
rz(-2.9031624) q[1];
rz(-pi) q[2];
rz(2.161383) q[3];
sx q[3];
rz(-2.906153) q[3];
sx q[3];
rz(2.2058132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7369467) q[2];
sx q[2];
rz(-2.8179822) q[2];
sx q[2];
rz(-1.7741989) q[2];
rz(-0.5395475) q[3];
sx q[3];
rz(-1.5522141) q[3];
sx q[3];
rz(1.3246983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0668199) q[0];
sx q[0];
rz(-2.9434151) q[0];
sx q[0];
rz(-2.5812126) q[0];
rz(2.9226774) q[1];
sx q[1];
rz(-1.0812662) q[1];
sx q[1];
rz(2.970649) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4881011) q[0];
sx q[0];
rz(-0.91918901) q[0];
sx q[0];
rz(2.1733858) q[0];
x q[1];
rz(-0.20410164) q[2];
sx q[2];
rz(-1.4320254) q[2];
sx q[2];
rz(-0.52817527) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.584632) q[1];
sx q[1];
rz(-1.0611532) q[1];
sx q[1];
rz(-2.4181448) q[1];
x q[2];
rz(1.6690977) q[3];
sx q[3];
rz(-2.519033) q[3];
sx q[3];
rz(0.05427256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7818452) q[2];
sx q[2];
rz(-1.1283987) q[2];
sx q[2];
rz(0.95853364) q[2];
rz(0.16252276) q[3];
sx q[3];
rz(-2.2992117) q[3];
sx q[3];
rz(1.5490279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.573134) q[0];
sx q[0];
rz(-3.077226) q[0];
sx q[0];
rz(2.3934613) q[0];
rz(2.023078) q[1];
sx q[1];
rz(-2.7225814) q[1];
sx q[1];
rz(0.31707877) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40561179) q[0];
sx q[0];
rz(-1.5026341) q[0];
sx q[0];
rz(-2.8951485) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5218614) q[2];
sx q[2];
rz(-1.7669919) q[2];
sx q[2];
rz(2.115415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.78817716) q[1];
sx q[1];
rz(-1.4995575) q[1];
sx q[1];
rz(-0.80435462) q[1];
x q[2];
rz(3.1012332) q[3];
sx q[3];
rz(-0.77733126) q[3];
sx q[3];
rz(3.0632927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3743484) q[2];
sx q[2];
rz(-0.92864645) q[2];
sx q[2];
rz(-1.413215) q[2];
rz(3.0610541) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(0.22182375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0291075) q[0];
sx q[0];
rz(-1.2487829) q[0];
sx q[0];
rz(-1.2741733) q[0];
rz(-1.796465) q[1];
sx q[1];
rz(-2.3990264) q[1];
sx q[1];
rz(-1.8003731) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182512) q[0];
sx q[0];
rz(-2.3592383) q[0];
sx q[0];
rz(-0.73114354) q[0];
rz(-pi) q[1];
rz(-0.48168452) q[2];
sx q[2];
rz(-1.8284214) q[2];
sx q[2];
rz(-0.017901808) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6051424) q[1];
sx q[1];
rz(-1.4172232) q[1];
sx q[1];
rz(1.7642412) q[1];
rz(-pi) q[2];
rz(-2.2209211) q[3];
sx q[3];
rz(-1.6643644) q[3];
sx q[3];
rz(-0.9182932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.207927) q[2];
sx q[2];
rz(-1.1873446) q[2];
sx q[2];
rz(-0.46736091) q[2];
rz(0.76006877) q[3];
sx q[3];
rz(-1.4773388) q[3];
sx q[3];
rz(1.2490341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46045983) q[0];
sx q[0];
rz(-2.12119) q[0];
sx q[0];
rz(1.4469752) q[0];
rz(0.32577062) q[1];
sx q[1];
rz(-2.9278946) q[1];
sx q[1];
rz(-1.5209341) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6353519) q[0];
sx q[0];
rz(-1.8367447) q[0];
sx q[0];
rz(2.6145934) q[0];
rz(-pi) q[1];
rz(-1.5704186) q[2];
sx q[2];
rz(-1.6624358) q[2];
sx q[2];
rz(2.3204892) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7829166) q[1];
sx q[1];
rz(-2.7442928) q[1];
sx q[1];
rz(2.0372169) q[1];
x q[2];
rz(-2.2748442) q[3];
sx q[3];
rz(-1.4102077) q[3];
sx q[3];
rz(-0.15204568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1763566) q[2];
sx q[2];
rz(-2.3945645) q[2];
sx q[2];
rz(2.526324) q[2];
rz(-1.8687013) q[3];
sx q[3];
rz(-2.8764909) q[3];
sx q[3];
rz(-2.7191775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.0203005) q[0];
sx q[0];
rz(-1.2739807) q[0];
sx q[0];
rz(0.85025382) q[0];
rz(-2.0203159) q[1];
sx q[1];
rz(-1.3669776) q[1];
sx q[1];
rz(-0.77176315) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5284755) q[0];
sx q[0];
rz(-0.69418797) q[0];
sx q[0];
rz(-1.5334237) q[0];
rz(-pi) q[1];
rz(-2.775101) q[2];
sx q[2];
rz(-2.9482609) q[2];
sx q[2];
rz(-1.9654578) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62006751) q[1];
sx q[1];
rz(-1.3273393) q[1];
sx q[1];
rz(-3.0848178) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2617175) q[3];
sx q[3];
rz(-2.7304683) q[3];
sx q[3];
rz(0.087669186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.11543342) q[2];
sx q[2];
rz(-1.5449646) q[2];
sx q[2];
rz(-0.8405295) q[2];
rz(-3.0606411) q[3];
sx q[3];
rz(-2.5637124) q[3];
sx q[3];
rz(-2.8988083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42846546) q[0];
sx q[0];
rz(-2.9436538) q[0];
sx q[0];
rz(1.1567098) q[0];
rz(-1.0276065) q[1];
sx q[1];
rz(-1.3778069) q[1];
sx q[1];
rz(-2.3311232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76475984) q[0];
sx q[0];
rz(-0.78928052) q[0];
sx q[0];
rz(0.30662243) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9058305) q[2];
sx q[2];
rz(-1.050569) q[2];
sx q[2];
rz(-0.65641415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9864784) q[1];
sx q[1];
rz(-0.41448516) q[1];
sx q[1];
rz(-2.3398967) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4224206) q[3];
sx q[3];
rz(-2.1273489) q[3];
sx q[3];
rz(-2.839193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.345574) q[2];
sx q[2];
rz(-1.2556262) q[2];
sx q[2];
rz(-1.5239117) q[2];
rz(-2.0412622) q[3];
sx q[3];
rz(-1.1805781) q[3];
sx q[3];
rz(2.9617917) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0236459) q[0];
sx q[0];
rz(-1.4143586) q[0];
sx q[0];
rz(2.216862) q[0];
rz(0.042451518) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(1.1645198) q[2];
sx q[2];
rz(-1.5193408) q[2];
sx q[2];
rz(1.2035412) q[2];
rz(-1.6067947) q[3];
sx q[3];
rz(-2.0279584) q[3];
sx q[3];
rz(1.8215836) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
