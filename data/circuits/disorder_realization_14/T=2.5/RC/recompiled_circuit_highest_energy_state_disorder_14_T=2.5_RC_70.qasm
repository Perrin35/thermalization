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
rz(2.4013588) q[0];
sx q[0];
rz(-1.6594247) q[0];
sx q[0];
rz(0.33492127) q[0];
rz(0.51796335) q[1];
sx q[1];
rz(5.2809102) q[1];
sx q[1];
rz(10.032293) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54017055) q[0];
sx q[0];
rz(-2.5937732) q[0];
sx q[0];
rz(1.1146077) q[0];
rz(-pi) q[1];
rz(3.0246441) q[2];
sx q[2];
rz(-1.0643355) q[2];
sx q[2];
rz(-3.1053271) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7360037) q[1];
sx q[1];
rz(-1.0936671) q[1];
sx q[1];
rz(-1.7504343) q[1];
rz(2.7685952) q[3];
sx q[3];
rz(-1.9891051) q[3];
sx q[3];
rz(-1.7278863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5241549) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(-0.66317916) q[2];
rz(-3.0607306) q[3];
sx q[3];
rz(-0.20564779) q[3];
sx q[3];
rz(1.9614356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8993768) q[0];
sx q[0];
rz(-1.7689393) q[0];
sx q[0];
rz(-2.2826165) q[0];
rz(-1.8513177) q[1];
sx q[1];
rz(-1.4651508) q[1];
sx q[1];
rz(1.7346409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55757421) q[0];
sx q[0];
rz(-2.0502591) q[0];
sx q[0];
rz(0.25734253) q[0];
rz(-pi) q[1];
x q[1];
rz(0.036769899) q[2];
sx q[2];
rz(-1.9027793) q[2];
sx q[2];
rz(0.22184243) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7125852) q[1];
sx q[1];
rz(-0.52118329) q[1];
sx q[1];
rz(1.3888478) q[1];
rz(1.9858667) q[3];
sx q[3];
rz(-1.0336116) q[3];
sx q[3];
rz(-1.9946757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0591639) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(1.814369) q[2];
rz(1.0446769) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(2.7248342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5663719) q[0];
sx q[0];
rz(-0.91082585) q[0];
sx q[0];
rz(0.4253934) q[0];
rz(-1.7644024) q[1];
sx q[1];
rz(-1.4981937) q[1];
sx q[1];
rz(1.4345217) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98374983) q[0];
sx q[0];
rz(-0.93261496) q[0];
sx q[0];
rz(-1.4236091) q[0];
x q[1];
rz(-1.6041669) q[2];
sx q[2];
rz(-2.2788308) q[2];
sx q[2];
rz(-1.5607426) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94043865) q[1];
sx q[1];
rz(-2.33719) q[1];
sx q[1];
rz(2.5522347) q[1];
rz(-pi) q[2];
rz(2.8355153) q[3];
sx q[3];
rz(-0.97460213) q[3];
sx q[3];
rz(-1.5709189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44125685) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(-1.2909935) q[2];
rz(-1.357632) q[3];
sx q[3];
rz(-1.5057526) q[3];
sx q[3];
rz(1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.5212379) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(-0.77019101) q[0];
rz(2.1417446) q[1];
sx q[1];
rz(-0.60037535) q[1];
sx q[1];
rz(-1.8005449) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4907341) q[0];
sx q[0];
rz(-0.58877173) q[0];
sx q[0];
rz(-0.38932271) q[0];
rz(0.45939873) q[2];
sx q[2];
rz(-0.93743491) q[2];
sx q[2];
rz(-0.49672302) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51777202) q[1];
sx q[1];
rz(-1.869535) q[1];
sx q[1];
rz(-1.2932528) q[1];
rz(-1.8561268) q[3];
sx q[3];
rz(-1.4754681) q[3];
sx q[3];
rz(-2.8872629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8009214) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(2.3731903) q[2];
rz(-3.0411804) q[3];
sx q[3];
rz(-1.6008335) q[3];
sx q[3];
rz(-1.2787308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95555821) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(2.9146063) q[0];
rz(-1.3817894) q[1];
sx q[1];
rz(-2.558936) q[1];
sx q[1];
rz(-1.7452128) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5433301) q[0];
sx q[0];
rz(-2.6285183) q[0];
sx q[0];
rz(-0.76900478) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3675186) q[2];
sx q[2];
rz(-2.705239) q[2];
sx q[2];
rz(-2.9038716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.13670838) q[1];
sx q[1];
rz(-1.4992104) q[1];
sx q[1];
rz(-0.63386713) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.030524039) q[3];
sx q[3];
rz(-0.82477335) q[3];
sx q[3];
rz(-1.8061226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22623006) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(3.0653595) q[2];
rz(1.3876312) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(1.4001728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48080322) q[0];
sx q[0];
rz(-1.5605518) q[0];
sx q[0];
rz(-1.3355108) q[0];
rz(2.238359) q[1];
sx q[1];
rz(-1.3390373) q[1];
sx q[1];
rz(-2.9551771) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7714027) q[0];
sx q[0];
rz(-1.0433084) q[0];
sx q[0];
rz(0.50748388) q[0];
rz(-pi) q[1];
rz(-1.769563) q[2];
sx q[2];
rz(-1.0093401) q[2];
sx q[2];
rz(1.1793062) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.820436) q[1];
sx q[1];
rz(-1.3296491) q[1];
sx q[1];
rz(-0.92280573) q[1];
rz(-0.78808727) q[3];
sx q[3];
rz(-2.4895146) q[3];
sx q[3];
rz(-1.8502082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4621801) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(1.818044) q[2];
rz(-1.1421674) q[3];
sx q[3];
rz(-2.3565632) q[3];
sx q[3];
rz(-0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46397504) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(2.5073945) q[0];
rz(-1.0448666) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(2.1814836) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2506367) q[0];
sx q[0];
rz(-1.5982369) q[0];
sx q[0];
rz(-1.1720042) q[0];
rz(-pi) q[1];
rz(-3.1158218) q[2];
sx q[2];
rz(-0.30116815) q[2];
sx q[2];
rz(-2.0370551) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83852488) q[1];
sx q[1];
rz(-0.52553287) q[1];
sx q[1];
rz(1.5477033) q[1];
rz(0.76310254) q[3];
sx q[3];
rz(-2.1402485) q[3];
sx q[3];
rz(-0.72009898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1977957) q[2];
sx q[2];
rz(-0.65239492) q[2];
sx q[2];
rz(-1.5956399) q[2];
rz(-1.4173077) q[3];
sx q[3];
rz(-2.3167819) q[3];
sx q[3];
rz(-1.5203169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5328131) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(2.1667495) q[0];
rz(-3.0335562) q[1];
sx q[1];
rz(-1.255475) q[1];
sx q[1];
rz(1.1688165) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0790137) q[0];
sx q[0];
rz(-2.4860408) q[0];
sx q[0];
rz(0.84540875) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5708099) q[2];
sx q[2];
rz(-2.332649) q[2];
sx q[2];
rz(-0.9521614) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4296234) q[1];
sx q[1];
rz(-1.483146) q[1];
sx q[1];
rz(1.3369249) q[1];
x q[2];
rz(2.6498472) q[3];
sx q[3];
rz(-1.5881268) q[3];
sx q[3];
rz(2.8239856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10890266) q[2];
sx q[2];
rz(-2.2638075) q[2];
sx q[2];
rz(2.4857793) q[2];
rz(-0.10410318) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(-0.18812215) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26466894) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(-1.6498097) q[0];
rz(-2.1022294) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(0.46844354) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2830848) q[0];
sx q[0];
rz(-1.383184) q[0];
sx q[0];
rz(-1.5664738) q[0];
x q[1];
rz(2.5889187) q[2];
sx q[2];
rz(-0.50853339) q[2];
sx q[2];
rz(-0.39297152) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5891582) q[1];
sx q[1];
rz(-2.3367408) q[1];
sx q[1];
rz(2.9167152) q[1];
rz(-pi) q[2];
rz(-2.3571068) q[3];
sx q[3];
rz(-2.753559) q[3];
sx q[3];
rz(-1.7540384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.09482065) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(0.92791933) q[2];
rz(1.3377442) q[3];
sx q[3];
rz(-0.51023054) q[3];
sx q[3];
rz(-1.6524338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054166404) q[0];
sx q[0];
rz(-1.0324284) q[0];
sx q[0];
rz(-1.077865) q[0];
rz(-0.36733356) q[1];
sx q[1];
rz(-1.8108188) q[1];
sx q[1];
rz(-1.8574538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8074984) q[0];
sx q[0];
rz(-2.0413412) q[0];
sx q[0];
rz(-0.28326359) q[0];
rz(-0.11576368) q[2];
sx q[2];
rz(-0.9204922) q[2];
sx q[2];
rz(0.34046003) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.974677) q[1];
sx q[1];
rz(-1.484718) q[1];
sx q[1];
rz(2.9831216) q[1];
rz(-pi) q[2];
rz(0.73478847) q[3];
sx q[3];
rz(-1.7289203) q[3];
sx q[3];
rz(-1.2722335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0014570634) q[2];
sx q[2];
rz(-2.6900901) q[2];
sx q[2];
rz(-2.7122811) q[2];
rz(2.9863827) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(-1.3252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24406381) q[0];
sx q[0];
rz(-1.5499935) q[0];
sx q[0];
rz(1.5503379) q[0];
rz(0.86391972) q[1];
sx q[1];
rz(-2.7680631) q[1];
sx q[1];
rz(1.6815129) q[1];
rz(1.3270039) q[2];
sx q[2];
rz(-2.1012123) q[2];
sx q[2];
rz(-0.13441555) q[2];
rz(1.1170837) q[3];
sx q[3];
rz(-2.4632005) q[3];
sx q[3];
rz(1.7076422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
