OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(-0.46407035) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(6.2160677) q[1];
sx q[1];
rz(10.709229) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0425134) q[0];
sx q[0];
rz(-0.49603841) q[0];
sx q[0];
rz(-0.92143671) q[0];
rz(1.60381) q[2];
sx q[2];
rz(-1.0812949) q[2];
sx q[2];
rz(-2.2061493) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1661108) q[1];
sx q[1];
rz(-1.2088641) q[1];
sx q[1];
rz(2.6544177) q[1];
x q[2];
rz(-1.5873317) q[3];
sx q[3];
rz(-1.4201122) q[3];
sx q[3];
rz(2.2076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39711943) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(-2.822067) q[2];
rz(2.5630991) q[3];
sx q[3];
rz(-0.47839034) q[3];
sx q[3];
rz(-2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330866) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(-2.615036) q[0];
rz(2.5780442) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(0.79663509) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731199) q[0];
sx q[0];
rz(-2.1935049) q[0];
sx q[0];
rz(0.40533439) q[0];
x q[1];
rz(1.7700427) q[2];
sx q[2];
rz(-2.1363878) q[2];
sx q[2];
rz(2.9233962) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9990275) q[1];
sx q[1];
rz(-1.5836645) q[1];
sx q[1];
rz(-2.5343115) q[1];
x q[2];
rz(0.49780952) q[3];
sx q[3];
rz(-2.1216765) q[3];
sx q[3];
rz(0.8964955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(0.49318796) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.440381) q[0];
rz(-0.72021833) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(-0.70297855) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5539348) q[0];
sx q[0];
rz(-0.30447391) q[0];
sx q[0];
rz(-1.1712043) q[0];
rz(1.9269283) q[2];
sx q[2];
rz(-2.0424358) q[2];
sx q[2];
rz(-2.4871662) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2762201) q[1];
sx q[1];
rz(-3.0225388) q[1];
sx q[1];
rz(2.7369376) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94743518) q[3];
sx q[3];
rz(-1.4294335) q[3];
sx q[3];
rz(2.3879104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(-2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(0.5660969) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(-2.3655868) q[0];
rz(1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(2.5783096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11382064) q[0];
sx q[0];
rz(-1.8437244) q[0];
sx q[0];
rz(-2.2233637) q[0];
x q[1];
rz(-1.2345384) q[2];
sx q[2];
rz(-2.2328937) q[2];
sx q[2];
rz(-1.7788356) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5717585) q[1];
sx q[1];
rz(-2.8697439) q[1];
sx q[1];
rz(3.1175201) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8654278) q[3];
sx q[3];
rz(-0.8751103) q[3];
sx q[3];
rz(-0.61616117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.48878601) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(2.9768067) q[2];
rz(-2.9131043) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8673458) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(-0.85246032) q[0];
rz(2.7903941) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.16211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64415414) q[0];
sx q[0];
rz(-0.69735202) q[0];
sx q[0];
rz(-2.5028412) q[0];
rz(-pi) q[1];
rz(-0.39761333) q[2];
sx q[2];
rz(-1.0828472) q[2];
sx q[2];
rz(0.60148009) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3106766) q[1];
sx q[1];
rz(-1.984593) q[1];
sx q[1];
rz(-0.78891854) q[1];
x q[2];
rz(-3.1245072) q[3];
sx q[3];
rz(-0.72165976) q[3];
sx q[3];
rz(0.34860308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6953485) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(-1.8943141) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(-0.22918992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86876774) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7346142) q[0];
sx q[0];
rz(-0.66440551) q[0];
sx q[0];
rz(2.2218496) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7734217) q[2];
sx q[2];
rz(-1.1154004) q[2];
sx q[2];
rz(0.95603285) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2893147) q[1];
sx q[1];
rz(-1.5783974) q[1];
sx q[1];
rz(-0.72035933) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2209729) q[3];
sx q[3];
rz(-0.60022012) q[3];
sx q[3];
rz(-0.61629399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7541472) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(2.5409017) q[2];
rz(-2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(0.35282648) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3951185) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(-1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(0.41710645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55506599) q[0];
sx q[0];
rz(-1.546372) q[0];
sx q[0];
rz(-1.7909554) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8344526) q[2];
sx q[2];
rz(-1.8516314) q[2];
sx q[2];
rz(-1.3020696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5187832) q[1];
sx q[1];
rz(-2.2768524) q[1];
sx q[1];
rz(0.6219567) q[1];
rz(0.65317513) q[3];
sx q[3];
rz(-1.8110868) q[3];
sx q[3];
rz(1.0550635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(2.588429) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(-2.6161391) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32507867) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(-1.1897855) q[0];
rz(-1.7143543) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(0.11238012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1988797) q[0];
sx q[0];
rz(-0.50857022) q[0];
sx q[0];
rz(-1.93165) q[0];
rz(-1.4833647) q[2];
sx q[2];
rz(-1.9033252) q[2];
sx q[2];
rz(-2.33193) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1689414) q[1];
sx q[1];
rz(-0.87247889) q[1];
sx q[1];
rz(0.45291839) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12179575) q[3];
sx q[3];
rz(-1.8276916) q[3];
sx q[3];
rz(-1.4956054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(1.3486264) q[2];
rz(-1.2049234) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444721) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(0.034974139) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(2.2299178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1393136) q[0];
sx q[0];
rz(-2.7873435) q[0];
sx q[0];
rz(2.2646963) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4439092) q[2];
sx q[2];
rz(-1.03089) q[2];
sx q[2];
rz(-1.1586231) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94757838) q[1];
sx q[1];
rz(-1.36735) q[1];
sx q[1];
rz(-0.46781637) q[1];
x q[2];
rz(2.0613641) q[3];
sx q[3];
rz(-1.0883696) q[3];
sx q[3];
rz(0.8182984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3578167) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(0.37208474) q[0];
rz(0.58139873) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(-1.4153597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85522643) q[0];
sx q[0];
rz(-0.33272538) q[0];
sx q[0];
rz(2.4871728) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0572817) q[2];
sx q[2];
rz(-2.1272749) q[2];
sx q[2];
rz(1.8884115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49299875) q[1];
sx q[1];
rz(-0.10222888) q[1];
sx q[1];
rz(0.74536721) q[1];
rz(-pi) q[2];
rz(2.1595702) q[3];
sx q[3];
rz(-0.62871274) q[3];
sx q[3];
rz(-2.2002937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(-2.9369205) q[2];
rz(-1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(-2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50080147) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(-1.5851371) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(-1.3265058) q[2];
sx q[2];
rz(-0.46637022) q[2];
sx q[2];
rz(-2.7844219) q[2];
rz(-0.08269357) q[3];
sx q[3];
rz(-1.0001928) q[3];
sx q[3];
rz(2.115888) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];