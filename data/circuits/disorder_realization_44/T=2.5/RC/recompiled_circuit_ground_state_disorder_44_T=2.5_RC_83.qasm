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
rz(-2.1739668) q[1];
sx q[1];
rz(-0.014048014) q[1];
sx q[1];
rz(-0.79298055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4913038) q[0];
sx q[0];
rz(-0.86435917) q[0];
sx q[0];
rz(-0.44235787) q[0];
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
sx q[0];
rz(-pi/2) q[0];
rz(-2.2574764) q[1];
sx q[1];
rz(-0.82472825) q[1];
sx q[1];
rz(-2.5349239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6019117) q[3];
sx q[3];
rz(-2.0869533) q[3];
sx q[3];
rz(1.8208139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.27796081) q[2];
sx q[2];
rz(-2.7148254) q[2];
sx q[2];
rz(0.58941907) q[2];
rz(-2.3943118) q[3];
sx q[3];
rz(-0.092844754) q[3];
sx q[3];
rz(-0.59982991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079161949) q[0];
sx q[0];
rz(-2.9099162) q[0];
sx q[0];
rz(-0.9011426) q[0];
rz(0.98183739) q[1];
sx q[1];
rz(-2.5603309) q[1];
sx q[1];
rz(0.99220401) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3311148) q[0];
sx q[0];
rz(-1.4338636) q[0];
sx q[0];
rz(-0.94816072) q[0];
rz(-1.7478701) q[2];
sx q[2];
rz(-2.3858527) q[2];
sx q[2];
rz(2.9301639) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8507311) q[1];
sx q[1];
rz(-2.2358316) q[1];
sx q[1];
rz(0.70677251) q[1];
rz(2.2295932) q[3];
sx q[3];
rz(-0.9992632) q[3];
sx q[3];
rz(-2.6603572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.09484) q[2];
sx q[2];
rz(-0.53043008) q[2];
sx q[2];
rz(-0.025010427) q[2];
rz(-0.95995861) q[3];
sx q[3];
rz(-3.0955866) q[3];
sx q[3];
rz(-0.068232603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9571824) q[0];
sx q[0];
rz(-0.010951696) q[0];
sx q[0];
rz(0.72407323) q[0];
rz(-3.030576) q[1];
sx q[1];
rz(-2.5357775) q[1];
sx q[1];
rz(0.012880005) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52791053) q[0];
sx q[0];
rz(-2.9146195) q[0];
sx q[0];
rz(1.9274953) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78978844) q[2];
sx q[2];
rz(-0.26726535) q[2];
sx q[2];
rz(-2.546026) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47564214) q[1];
sx q[1];
rz(-0.29120177) q[1];
sx q[1];
rz(-0.61333527) q[1];
rz(-pi) q[2];
rz(-2.2575081) q[3];
sx q[3];
rz(-1.8662765) q[3];
sx q[3];
rz(-2.0990163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8956902) q[2];
sx q[2];
rz(-2.4195713) q[2];
sx q[2];
rz(2.8172909) q[2];
rz(-0.4970099) q[3];
sx q[3];
rz(-0.23247601) q[3];
sx q[3];
rz(0.23268172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4683891) q[0];
sx q[0];
rz(-0.16981801) q[0];
sx q[0];
rz(0.47624269) q[0];
rz(0.4854804) q[1];
sx q[1];
rz(-0.49518934) q[1];
sx q[1];
rz(-0.21864299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9354585) q[0];
sx q[0];
rz(-1.1674321) q[0];
sx q[0];
rz(2.8950188) q[0];
x q[1];
rz(-0.013979023) q[2];
sx q[2];
rz(-2.4179248) q[2];
sx q[2];
rz(-2.6633563) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2610332) q[1];
sx q[1];
rz(-1.1487242) q[1];
sx q[1];
rz(0.5725533) q[1];
rz(-pi) q[2];
rz(2.8273638) q[3];
sx q[3];
rz(-1.7933729) q[3];
sx q[3];
rz(1.4861388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3525047) q[2];
sx q[2];
rz(-0.75864351) q[2];
sx q[2];
rz(-0.57141203) q[2];
rz(-2.2032951) q[3];
sx q[3];
rz(-2.5550227) q[3];
sx q[3];
rz(0.17294426) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57662302) q[0];
sx q[0];
rz(-0.40410703) q[0];
sx q[0];
rz(-0.47082666) q[0];
rz(-2.2305523) q[1];
sx q[1];
rz(-0.43993479) q[1];
sx q[1];
rz(-0.43112531) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79589383) q[0];
sx q[0];
rz(-2.3568793) q[0];
sx q[0];
rz(3.1036397) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9091102) q[2];
sx q[2];
rz(-0.66146427) q[2];
sx q[2];
rz(2.4544883) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4404802) q[1];
sx q[1];
rz(-1.639469) q[1];
sx q[1];
rz(-3.1140226) q[1];
x q[2];
rz(-2.9399559) q[3];
sx q[3];
rz(-1.4343702) q[3];
sx q[3];
rz(-2.0755943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.33991459) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(-0.3413631) q[2];
rz(-2.7692128) q[3];
sx q[3];
rz(-2.5058993) q[3];
sx q[3];
rz(-2.671833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.85207283) q[0];
sx q[0];
rz(-0.85318035) q[0];
sx q[0];
rz(-0.12292718) q[0];
rz(-0.3854824) q[1];
sx q[1];
rz(-2.3434134) q[1];
sx q[1];
rz(-0.72365671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0892886) q[0];
sx q[0];
rz(-1.6160242) q[0];
sx q[0];
rz(-1.3967229) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9060109) q[2];
sx q[2];
rz(-1.7487757) q[2];
sx q[2];
rz(-0.92648695) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9229615) q[1];
sx q[1];
rz(-2.1653163) q[1];
sx q[1];
rz(2.2943952) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0615929) q[3];
sx q[3];
rz(-0.43269581) q[3];
sx q[3];
rz(3.1263292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3896997) q[2];
sx q[2];
rz(-2.5395826) q[2];
sx q[2];
rz(-2.4683118) q[2];
rz(-0.26345396) q[3];
sx q[3];
rz(-2.6665688) q[3];
sx q[3];
rz(-0.30509216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70169705) q[0];
sx q[0];
rz(-0.79653996) q[0];
sx q[0];
rz(-0.51629603) q[0];
rz(0.75597489) q[1];
sx q[1];
rz(-0.95840234) q[1];
sx q[1];
rz(-0.36852536) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8719292) q[0];
sx q[0];
rz(-2.2145529) q[0];
sx q[0];
rz(2.4764349) q[0];
x q[1];
rz(3.0958648) q[2];
sx q[2];
rz(-2.4362982) q[2];
sx q[2];
rz(0.7208342) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8766899) q[1];
sx q[1];
rz(-1.6769092) q[1];
sx q[1];
rz(1.7404895) q[1];
x q[2];
rz(-2.6266699) q[3];
sx q[3];
rz(-2.2201) q[3];
sx q[3];
rz(-0.89735555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64510173) q[2];
sx q[2];
rz(-3.0848905) q[2];
sx q[2];
rz(0.48218316) q[2];
rz(2.9218946) q[3];
sx q[3];
rz(-0.69742656) q[3];
sx q[3];
rz(-2.4612332) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7242303) q[0];
sx q[0];
rz(-0.14886947) q[0];
sx q[0];
rz(-2.505488) q[0];
rz(0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(2.2669534) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8281404) q[0];
sx q[0];
rz(-0.16022542) q[0];
sx q[0];
rz(-1.1988706) q[0];
rz(-1.6804439) q[2];
sx q[2];
rz(-1.4504408) q[2];
sx q[2];
rz(-0.14794825) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.88747178) q[1];
sx q[1];
rz(-1.1855108) q[1];
sx q[1];
rz(-1.1347358) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88040027) q[3];
sx q[3];
rz(-1.0171618) q[3];
sx q[3];
rz(2.2601489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.26802289) q[2];
sx q[2];
rz(-2.7251254) q[2];
sx q[2];
rz(-0.91656172) q[2];
rz(0.77740866) q[3];
sx q[3];
rz(-0.55792266) q[3];
sx q[3];
rz(2.6072445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1189608) q[0];
sx q[0];
rz(-0.73771483) q[0];
sx q[0];
rz(0.23271261) q[0];
rz(-2.4932056) q[1];
sx q[1];
rz(-2.5189923) q[1];
sx q[1];
rz(-2.5737305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0584326) q[0];
sx q[0];
rz(-0.55480236) q[0];
sx q[0];
rz(0.85809274) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8428734) q[2];
sx q[2];
rz(-1.7781742) q[2];
sx q[2];
rz(-0.033471154) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76946867) q[1];
sx q[1];
rz(-0.25126496) q[1];
sx q[1];
rz(0.94458802) q[1];
x q[2];
rz(0.81553163) q[3];
sx q[3];
rz(-0.99247265) q[3];
sx q[3];
rz(-2.2165143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76967543) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(0.5522716) q[2];
rz(0.28198379) q[3];
sx q[3];
rz(-0.43006399) q[3];
sx q[3];
rz(-0.10274398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.929739) q[0];
sx q[0];
rz(-3.0775098) q[0];
sx q[0];
rz(0.1317568) q[0];
rz(-2.6051104) q[1];
sx q[1];
rz(-0.15833144) q[1];
sx q[1];
rz(-0.90824711) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.344438) q[0];
sx q[0];
rz(-0.69509172) q[0];
sx q[0];
rz(1.6174843) q[0];
x q[1];
rz(-2.0591878) q[2];
sx q[2];
rz(-2.6219212) q[2];
sx q[2];
rz(-0.95275098) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6301292) q[1];
sx q[1];
rz(-0.88246934) q[1];
sx q[1];
rz(-0.87911112) q[1];
x q[2];
rz(-1.796714) q[3];
sx q[3];
rz(-1.5690256) q[3];
sx q[3];
rz(1.8670807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6324255) q[2];
sx q[2];
rz(-0.91370344) q[2];
sx q[2];
rz(2.8636279) q[2];
rz(0.5667423) q[3];
sx q[3];
rz(-2.9394579) q[3];
sx q[3];
rz(0.91293269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81000281) q[0];
sx q[0];
rz(-1.4053874) q[0];
sx q[0];
rz(-0.94595861) q[0];
rz(-2.6592061) q[1];
sx q[1];
rz(-1.2336029) q[1];
sx q[1];
rz(-0.69419669) q[1];
rz(-1.214477) q[2];
sx q[2];
rz(-0.6206442) q[2];
sx q[2];
rz(0.71833687) q[2];
rz(-3.1101221) q[3];
sx q[3];
rz(-2.0768055) q[3];
sx q[3];
rz(0.16184645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
