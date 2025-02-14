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
rz(1.1531416) q[0];
sx q[0];
rz(2.3260131) q[0];
sx q[0];
rz(10.182967) q[0];
rz(-0.012501333) q[1];
sx q[1];
rz(4.9795436) q[1];
sx q[1];
rz(10.995168) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1364098) q[0];
sx q[0];
rz(-1.9290975) q[0];
sx q[0];
rz(-0.65922064) q[0];
rz(-3.0742253) q[2];
sx q[2];
rz(-1.172003) q[2];
sx q[2];
rz(-0.68902868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6627854) q[1];
sx q[1];
rz(-0.36514716) q[1];
sx q[1];
rz(1.7313596) q[1];
rz(0.060293555) q[3];
sx q[3];
rz(-1.2620842) q[3];
sx q[3];
rz(-1.2293881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5300753) q[2];
sx q[2];
rz(-3.1334183) q[2];
sx q[2];
rz(0.44069904) q[2];
rz(-3.0618482) q[3];
sx q[3];
rz(-0.00010448797) q[3];
sx q[3];
rz(-1.124148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9663064) q[0];
sx q[0];
rz(-0.056862406) q[0];
sx q[0];
rz(-2.9735907) q[0];
rz(-0.02027823) q[1];
sx q[1];
rz(-0.30826491) q[1];
sx q[1];
rz(1.5365938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5791721) q[0];
sx q[0];
rz(-1.354166) q[0];
sx q[0];
rz(-0.24026339) q[0];
x q[1];
rz(1.5606784) q[2];
sx q[2];
rz(-1.5732906) q[2];
sx q[2];
rz(-0.025257142) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78950845) q[1];
sx q[1];
rz(-1.5718954) q[1];
sx q[1];
rz(-1.5754682) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9320967) q[3];
sx q[3];
rz(-3.0983546) q[3];
sx q[3];
rz(0.87856495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7961879) q[2];
sx q[2];
rz(-0.91190839) q[2];
sx q[2];
rz(1.4076642) q[2];
rz(1.0450854) q[3];
sx q[3];
rz(-3.0920691) q[3];
sx q[3];
rz(-0.27200562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.3097836) q[0];
sx q[0];
rz(-0.97667664) q[0];
sx q[0];
rz(-0.56104863) q[0];
rz(-0.27684119) q[1];
sx q[1];
rz(-0.012877348) q[1];
sx q[1];
rz(-1.8337839) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7335199) q[0];
sx q[0];
rz(-1.8674984) q[0];
sx q[0];
rz(-0.042859239) q[0];
x q[1];
rz(-1.5707914) q[2];
sx q[2];
rz(-1.5781856) q[2];
sx q[2];
rz(2.0781197) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4958932) q[1];
sx q[1];
rz(-0.99826854) q[1];
sx q[1];
rz(-0.069620274) q[1];
x q[2];
rz(-2.2482292) q[3];
sx q[3];
rz(-1.9382538) q[3];
sx q[3];
rz(1.8830255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3448559) q[2];
sx q[2];
rz(-0.00011809706) q[2];
sx q[2];
rz(-0.5893839) q[2];
rz(-1.2423337) q[3];
sx q[3];
rz(-3.1292249) q[3];
sx q[3];
rz(-1.7813659) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1347443) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(-1.7929329) q[0];
rz(3.1349365) q[1];
sx q[1];
rz(-1.321188) q[1];
sx q[1];
rz(-0.033500813) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.861167) q[0];
sx q[0];
rz(-2.4651732) q[0];
sx q[0];
rz(1.8489714) q[0];
rz(-pi) q[1];
rz(3.1042023) q[2];
sx q[2];
rz(-3.0218736) q[2];
sx q[2];
rz(2.7604288) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12043145) q[1];
sx q[1];
rz(-2.8725) q[1];
sx q[1];
rz(-1.6090367) q[1];
rz(-pi) q[2];
rz(1.4064404) q[3];
sx q[3];
rz(-1.4504315) q[3];
sx q[3];
rz(0.09935483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.032430705) q[2];
sx q[2];
rz(-0.0065294821) q[2];
sx q[2];
rz(-2.8429441) q[2];
rz(-1.8715035) q[3];
sx q[3];
rz(-3.1254369) q[3];
sx q[3];
rz(-0.0531918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.771516) q[0];
sx q[0];
rz(-1.5144441) q[0];
sx q[0];
rz(-0.44684967) q[0];
rz(0.18613786) q[1];
sx q[1];
rz(-3.0797854) q[1];
sx q[1];
rz(1.4131379) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7639358) q[0];
sx q[0];
rz(-1.3400338) q[0];
sx q[0];
rz(3.1070437) q[0];
x q[1];
rz(0.21633001) q[2];
sx q[2];
rz(-1.9842752) q[2];
sx q[2];
rz(-1.0461996) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2286716) q[1];
sx q[1];
rz(-0.067021772) q[1];
sx q[1];
rz(-2.3093501) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2373393) q[3];
sx q[3];
rz(-2.7891141) q[3];
sx q[3];
rz(-0.74732399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83429217) q[2];
sx q[2];
rz(-1.5520381) q[2];
sx q[2];
rz(-2.6429122) q[2];
rz(-0.57791609) q[3];
sx q[3];
rz(-2.6579865) q[3];
sx q[3];
rz(2.5448866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.1805098) q[0];
sx q[0];
rz(-1.120765) q[0];
sx q[0];
rz(0.34641308) q[0];
rz(0.60180426) q[1];
sx q[1];
rz(-1.5608984) q[1];
sx q[1];
rz(-2.3880889) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44418535) q[0];
sx q[0];
rz(-1.7427398) q[0];
sx q[0];
rz(-0.19449046) q[0];
rz(-pi) q[1];
rz(-1.4592811) q[2];
sx q[2];
rz(-1.6454305) q[2];
sx q[2];
rz(2.1837051) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9484472) q[1];
sx q[1];
rz(-2.3420534) q[1];
sx q[1];
rz(-0.97481291) q[1];
rz(1.7338792) q[3];
sx q[3];
rz(-1.661851) q[3];
sx q[3];
rz(2.1084821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5715013) q[2];
sx q[2];
rz(-0.0034905958) q[2];
sx q[2];
rz(1.5241148) q[2];
rz(-3.0213455) q[3];
sx q[3];
rz(-3.1382939) q[3];
sx q[3];
rz(-2.6038468) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26871249) q[0];
sx q[0];
rz(-0.92395067) q[0];
sx q[0];
rz(-0.22126108) q[0];
rz(-1.6809173) q[1];
sx q[1];
rz(-0.9333846) q[1];
sx q[1];
rz(-0.077300765) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50718302) q[0];
sx q[0];
rz(-1.5290676) q[0];
sx q[0];
rz(0.001464837) q[0];
rz(-0.004782025) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(2.9298669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0127924) q[1];
sx q[1];
rz(-0.18928738) q[1];
sx q[1];
rz(0.38893338) q[1];
x q[2];
rz(1.037498) q[3];
sx q[3];
rz(-0.24720705) q[3];
sx q[3];
rz(-2.5105298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7920502) q[2];
sx q[2];
rz(-3.1304066) q[2];
sx q[2];
rz(-2.1816317) q[2];
rz(0.32944426) q[3];
sx q[3];
rz(-3.1335148) q[3];
sx q[3];
rz(0.85028696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9462117) q[0];
sx q[0];
rz(-0.61310261) q[0];
sx q[0];
rz(0.10928133) q[0];
rz(-0.37336135) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(-1.9111309) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5143573) q[0];
sx q[0];
rz(-2.1684596) q[0];
sx q[0];
rz(-2.3742832) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76704278) q[2];
sx q[2];
rz(-2.8708176) q[2];
sx q[2];
rz(2.3363638) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8459863) q[1];
sx q[1];
rz(-1.5822268) q[1];
sx q[1];
rz(1.5039526) q[1];
rz(1.968574) q[3];
sx q[3];
rz(-0.69044411) q[3];
sx q[3];
rz(3.0544314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5750778) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(-1.8129978) q[2];
rz(-1.7447507) q[3];
sx q[3];
rz(-0.003740398) q[3];
sx q[3];
rz(-1.0218792) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0777271) q[0];
sx q[0];
rz(-1.4430178) q[0];
sx q[0];
rz(-0.57300895) q[0];
rz(-2.833448) q[1];
sx q[1];
rz(-0.40987086) q[1];
sx q[1];
rz(-1.0073957) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44237374) q[0];
sx q[0];
rz(-1.6773407) q[0];
sx q[0];
rz(-3.1359948) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2635055) q[2];
sx q[2];
rz(-0.14423926) q[2];
sx q[2];
rz(-0.34473127) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3904265) q[1];
sx q[1];
rz(-1.5301643) q[1];
sx q[1];
rz(1.4875814) q[1];
rz(-pi) q[2];
x q[2];
rz(1.151772) q[3];
sx q[3];
rz(-3.1085827) q[3];
sx q[3];
rz(-2.9275683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8241626) q[2];
sx q[2];
rz(-0.62676668) q[2];
sx q[2];
rz(-2.7516348) q[2];
rz(-0.069084875) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(-2.8100584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020141715) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(0.48594117) q[0];
rz(-2.2700229) q[1];
sx q[1];
rz(-1.3078682) q[1];
sx q[1];
rz(-1.4942687) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93372351) q[0];
sx q[0];
rz(-1.8759125) q[0];
sx q[0];
rz(-2.5795769) q[0];
x q[1];
rz(2.5267692) q[2];
sx q[2];
rz(-1.5346926) q[2];
sx q[2];
rz(1.6143798) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6657662) q[1];
sx q[1];
rz(-1.2428871) q[1];
sx q[1];
rz(-1.8880512) q[1];
x q[2];
rz(2.787519) q[3];
sx q[3];
rz(-0.1205509) q[3];
sx q[3];
rz(0.15124409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5668874) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(3.1097143) q[2];
rz(0.77830642) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(-0.2955029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42302172) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(-3.0146535) q[1];
sx q[1];
rz(-0.23902421) q[1];
sx q[1];
rz(0.21993266) q[1];
rz(1.710464) q[2];
sx q[2];
rz(-1.5763379) q[2];
sx q[2];
rz(-1.3272641) q[2];
rz(-0.056679139) q[3];
sx q[3];
rz(-0.89691848) q[3];
sx q[3];
rz(0.46700029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
