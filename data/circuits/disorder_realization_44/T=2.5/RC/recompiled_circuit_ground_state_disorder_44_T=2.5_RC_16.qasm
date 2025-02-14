OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39785102) q[0];
sx q[0];
rz(-2.1704817) q[0];
sx q[0];
rz(0.71075332) q[0];
rz(0.96762586) q[1];
sx q[1];
rz(3.1556407) q[1];
sx q[1];
rz(10.217759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4913038) q[0];
sx q[0];
rz(-2.2772335) q[0];
sx q[0];
rz(2.6992348) q[0];
rz(-pi) q[1];
x q[1];
rz(0.014292467) q[2];
sx q[2];
rz(-1.6184752) q[2];
sx q[2];
rz(-0.68916262) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0539581) q[1];
sx q[1];
rz(-2.2184508) q[1];
sx q[1];
rz(2.1235076) q[1];
rz(-0.6019117) q[3];
sx q[3];
rz(-2.0869533) q[3];
sx q[3];
rz(1.8208139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.27796081) q[2];
sx q[2];
rz(-2.7148254) q[2];
sx q[2];
rz(-0.58941907) q[2];
rz(-2.3943118) q[3];
sx q[3];
rz(-3.0487479) q[3];
sx q[3];
rz(0.59982991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0624307) q[0];
sx q[0];
rz(-0.23167647) q[0];
sx q[0];
rz(-0.9011426) q[0];
rz(2.1597553) q[1];
sx q[1];
rz(-0.58126175) q[1];
sx q[1];
rz(0.99220401) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57241523) q[0];
sx q[0];
rz(-2.506034) q[0];
sx q[0];
rz(-1.8028238) q[0];
rz(1.7478701) q[2];
sx q[2];
rz(-2.3858527) q[2];
sx q[2];
rz(0.21142879) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2908615) q[1];
sx q[1];
rz(-2.2358316) q[1];
sx q[1];
rz(2.4348201) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68282376) q[3];
sx q[3];
rz(-1.0299333) q[3];
sx q[3];
rz(1.6554494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.09484) q[2];
sx q[2];
rz(-2.6111626) q[2];
sx q[2];
rz(0.025010427) q[2];
rz(-0.95995861) q[3];
sx q[3];
rz(-0.046006087) q[3];
sx q[3];
rz(0.068232603) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9571824) q[0];
sx q[0];
rz(-3.130641) q[0];
sx q[0];
rz(0.72407323) q[0];
rz(-0.11101668) q[1];
sx q[1];
rz(-0.60581517) q[1];
sx q[1];
rz(-3.1287126) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89319481) q[0];
sx q[0];
rz(-1.7832558) q[0];
sx q[0];
rz(-0.080470632) q[0];
rz(-2.3518042) q[2];
sx q[2];
rz(-0.26726535) q[2];
sx q[2];
rz(0.59556669) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15808039) q[1];
sx q[1];
rz(-1.3338102) q[1];
sx q[1];
rz(-1.7416341) q[1];
rz(-pi) q[2];
rz(2.0184085) q[3];
sx q[3];
rz(-2.4035998) q[3];
sx q[3];
rz(-2.9546705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8956902) q[2];
sx q[2];
rz(-2.4195713) q[2];
sx q[2];
rz(-2.8172909) q[2];
rz(-2.6445828) q[3];
sx q[3];
rz(-2.9091166) q[3];
sx q[3];
rz(-2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4683891) q[0];
sx q[0];
rz(-0.16981801) q[0];
sx q[0];
rz(-0.47624269) q[0];
rz(-0.4854804) q[1];
sx q[1];
rz(-2.6464033) q[1];
sx q[1];
rz(2.9229497) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63594288) q[0];
sx q[0];
rz(-0.46920645) q[0];
sx q[0];
rz(1.0512661) q[0];
x q[1];
rz(0.013979023) q[2];
sx q[2];
rz(-2.4179248) q[2];
sx q[2];
rz(-0.47823634) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2656897) q[1];
sx q[1];
rz(-2.4445718) q[1];
sx q[1];
rz(0.69209309) q[1];
rz(-pi) q[2];
rz(-0.31422887) q[3];
sx q[3];
rz(-1.7933729) q[3];
sx q[3];
rz(-1.6554538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.3525047) q[2];
sx q[2];
rz(-0.75864351) q[2];
sx q[2];
rz(-2.5701806) q[2];
rz(0.93829751) q[3];
sx q[3];
rz(-0.58656991) q[3];
sx q[3];
rz(-0.17294426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57662302) q[0];
sx q[0];
rz(-0.40410703) q[0];
sx q[0];
rz(2.670766) q[0];
rz(0.91104031) q[1];
sx q[1];
rz(-0.43993479) q[1];
sx q[1];
rz(-0.43112531) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.39931) q[0];
sx q[0];
rz(-0.78680187) q[0];
sx q[0];
rz(1.5329225) q[0];
x q[1];
rz(2.2041915) q[2];
sx q[2];
rz(-1.3654815) q[2];
sx q[2];
rz(0.61287731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31897983) q[1];
sx q[1];
rz(-3.0676004) q[1];
sx q[1];
rz(1.1896108) q[1];
rz(-pi) q[2];
rz(-0.2016368) q[3];
sx q[3];
rz(-1.7072225) q[3];
sx q[3];
rz(-2.0755943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8016781) q[2];
sx q[2];
rz(-3.1368003) q[2];
sx q[2];
rz(2.8002296) q[2];
rz(0.3723799) q[3];
sx q[3];
rz(-0.63569331) q[3];
sx q[3];
rz(-0.46975964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.2895198) q[0];
sx q[0];
rz(-2.2884123) q[0];
sx q[0];
rz(3.0186655) q[0];
rz(0.3854824) q[1];
sx q[1];
rz(-2.3434134) q[1];
sx q[1];
rz(-2.4179359) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6521344) q[0];
sx q[0];
rz(-1.7446899) q[0];
sx q[0];
rz(-3.0956718) q[0];
rz(-1.2355818) q[2];
sx q[2];
rz(-1.3928169) q[2];
sx q[2];
rz(-2.2151057) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9229615) q[1];
sx q[1];
rz(-2.1653163) q[1];
sx q[1];
rz(0.84719744) q[1];
rz(-pi) q[2];
rz(-1.9576362) q[3];
sx q[3];
rz(-1.3718492) q[3];
sx q[3];
rz(-1.1038279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.75189292) q[2];
sx q[2];
rz(-2.5395826) q[2];
sx q[2];
rz(-0.67328084) q[2];
rz(0.26345396) q[3];
sx q[3];
rz(-2.6665688) q[3];
sx q[3];
rz(0.30509216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70169705) q[0];
sx q[0];
rz(-2.3450527) q[0];
sx q[0];
rz(-0.51629603) q[0];
rz(-0.75597489) q[1];
sx q[1];
rz(-2.1831903) q[1];
sx q[1];
rz(2.7730673) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8719292) q[0];
sx q[0];
rz(-0.9270398) q[0];
sx q[0];
rz(0.66515775) q[0];
rz(-pi) q[1];
rz(1.6096949) q[2];
sx q[2];
rz(-0.86639154) q[2];
sx q[2];
rz(2.4807841) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8175537) q[1];
sx q[1];
rz(-1.7395258) q[1];
sx q[1];
rz(3.0339452) q[1];
x q[2];
rz(-2.2880408) q[3];
sx q[3];
rz(-1.1677168) q[3];
sx q[3];
rz(2.1385101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64510173) q[2];
sx q[2];
rz(-3.0848905) q[2];
sx q[2];
rz(-2.6594095) q[2];
rz(-2.9218946) q[3];
sx q[3];
rz(-0.69742656) q[3];
sx q[3];
rz(2.4612332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7242303) q[0];
sx q[0];
rz(-2.9927232) q[0];
sx q[0];
rz(0.6361047) q[0];
rz(0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(2.2669534) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88976218) q[0];
sx q[0];
rz(-1.5127851) q[0];
sx q[0];
rz(1.7202352) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6804439) q[2];
sx q[2];
rz(-1.4504408) q[2];
sx q[2];
rz(-2.9936444) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.284901) q[1];
sx q[1];
rz(-1.1686348) q[1];
sx q[1];
rz(-0.42070893) q[1];
rz(-0.6757856) q[3];
sx q[3];
rz(-0.99832557) q[3];
sx q[3];
rz(-1.0990717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26802289) q[2];
sx q[2];
rz(-0.41646725) q[2];
sx q[2];
rz(2.2250309) q[2];
rz(-0.77740866) q[3];
sx q[3];
rz(-2.58367) q[3];
sx q[3];
rz(2.6072445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-3.1189608) q[0];
sx q[0];
rz(-2.4038778) q[0];
sx q[0];
rz(-2.90888) q[0];
rz(-0.64838707) q[1];
sx q[1];
rz(-2.5189923) q[1];
sx q[1];
rz(-0.56786215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0584326) q[0];
sx q[0];
rz(-2.5867903) q[0];
sx q[0];
rz(-2.2834999) q[0];
rz(-2.9265327) q[2];
sx q[2];
rz(-1.8369007) q[2];
sx q[2];
rz(-1.6616481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9514002) q[1];
sx q[1];
rz(-1.7170329) q[1];
sx q[1];
rz(1.7758572) q[1];
rz(-2.4106823) q[3];
sx q[3];
rz(-0.95940351) q[3];
sx q[3];
rz(2.9711593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76967543) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(-0.5522716) q[2];
rz(0.28198379) q[3];
sx q[3];
rz(-0.43006399) q[3];
sx q[3];
rz(-0.10274398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.929739) q[0];
sx q[0];
rz(-3.0775098) q[0];
sx q[0];
rz(3.0098359) q[0];
rz(-0.53648221) q[1];
sx q[1];
rz(-0.15833144) q[1];
sx q[1];
rz(-2.2333455) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2836766) q[0];
sx q[0];
rz(-2.2649797) q[0];
sx q[0];
rz(0.038900872) q[0];
rz(2.0386253) q[2];
sx q[2];
rz(-1.8059633) q[2];
sx q[2];
rz(2.9556592) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7365807) q[1];
sx q[1];
rz(-0.93376505) q[1];
sx q[1];
rz(0.65959658) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.578701) q[3];
sx q[3];
rz(-0.22592446) q[3];
sx q[3];
rz(0.28858063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.50916719) q[2];
sx q[2];
rz(-0.91370344) q[2];
sx q[2];
rz(0.2779648) q[2];
rz(-2.5748504) q[3];
sx q[3];
rz(-0.20213474) q[3];
sx q[3];
rz(2.22866) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3315898) q[0];
sx q[0];
rz(-1.4053874) q[0];
sx q[0];
rz(-0.94595861) q[0];
rz(0.48238659) q[1];
sx q[1];
rz(-1.2336029) q[1];
sx q[1];
rz(-0.69419669) q[1];
rz(-0.24438582) q[2];
sx q[2];
rz(-0.99437154) q[2];
sx q[2];
rz(-2.852358) q[2];
rz(2.0770155) q[3];
sx q[3];
rz(-1.5432705) q[3];
sx q[3];
rz(-1.3936925) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
