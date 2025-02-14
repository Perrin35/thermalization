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
rz(-1.350116) q[0];
sx q[0];
rz(3.8173563) q[0];
sx q[0];
rz(9.4326333) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(-1.5239198) q[1];
sx q[1];
rz(-0.24267264) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9822916) q[0];
sx q[0];
rz(-0.86573273) q[0];
sx q[0];
rz(0.16601913) q[0];
x q[1];
rz(1.2062293) q[2];
sx q[2];
rz(-2.0040214) q[2];
sx q[2];
rz(-1.792576) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3585904) q[1];
sx q[1];
rz(-1.5471518) q[1];
sx q[1];
rz(0.35332291) q[1];
rz(-0.66565973) q[3];
sx q[3];
rz(-1.375631) q[3];
sx q[3];
rz(1.5100556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0365888) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(0.71199065) q[2];
rz(-2.8365734) q[3];
sx q[3];
rz(-0.22235338) q[3];
sx q[3];
rz(-3.0960848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5104093) q[0];
sx q[0];
rz(-1.864186) q[0];
sx q[0];
rz(-1.3674059) q[0];
rz(-0.039904682) q[1];
sx q[1];
rz(-2.3343562) q[1];
sx q[1];
rz(2.5423999) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78135787) q[0];
sx q[0];
rz(-1.0971945) q[0];
sx q[0];
rz(-1.4417157) q[0];
rz(-pi) q[1];
rz(-2.3223898) q[2];
sx q[2];
rz(-1.1408198) q[2];
sx q[2];
rz(-2.5366304) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.42223052) q[1];
sx q[1];
rz(-1.9047202) q[1];
sx q[1];
rz(-2.4704598) q[1];
x q[2];
rz(-0.33984025) q[3];
sx q[3];
rz(-2.426894) q[3];
sx q[3];
rz(3.0106737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.082531884) q[2];
sx q[2];
rz(-2.579687) q[2];
sx q[2];
rz(-1.833029) q[2];
rz(0.023905309) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(1.5628975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.4557274) q[0];
sx q[0];
rz(-0.66615921) q[0];
sx q[0];
rz(-2.300793) q[0];
rz(-1.0288382) q[1];
sx q[1];
rz(-0.40887555) q[1];
sx q[1];
rz(1.6811446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59305916) q[0];
sx q[0];
rz(-1.4276299) q[0];
sx q[0];
rz(-2.6770704) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47346327) q[2];
sx q[2];
rz(-1.9046648) q[2];
sx q[2];
rz(-0.75227458) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3884133) q[1];
sx q[1];
rz(-2.3824661) q[1];
sx q[1];
rz(2.6469346) q[1];
rz(-pi) q[2];
rz(0.60734235) q[3];
sx q[3];
rz(-1.6173956) q[3];
sx q[3];
rz(-0.49639116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.2568405) q[2];
sx q[2];
rz(-1.2105056) q[2];
sx q[2];
rz(2.9849226) q[2];
rz(1.5001851) q[3];
sx q[3];
rz(-1.2431966) q[3];
sx q[3];
rz(0.18178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0750065) q[0];
sx q[0];
rz(-1.2970507) q[0];
sx q[0];
rz(-0.080667607) q[0];
rz(1.1219885) q[1];
sx q[1];
rz(-0.64067084) q[1];
sx q[1];
rz(1.5544308) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51475731) q[0];
sx q[0];
rz(-2.3957402) q[0];
sx q[0];
rz(-1.6677854) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0253554) q[2];
sx q[2];
rz(-1.6112865) q[2];
sx q[2];
rz(-1.7396648) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5579368) q[1];
sx q[1];
rz(-2.2018345) q[1];
sx q[1];
rz(2.2028752) q[1];
rz(-0.60685632) q[3];
sx q[3];
rz(-1.292406) q[3];
sx q[3];
rz(-0.75137072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2598205) q[2];
sx q[2];
rz(-1.0910923) q[2];
sx q[2];
rz(2.1176977) q[2];
rz(2.3648868) q[3];
sx q[3];
rz(-1.1506162) q[3];
sx q[3];
rz(-1.0030494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40989947) q[0];
sx q[0];
rz(-0.85407805) q[0];
sx q[0];
rz(-0.62752974) q[0];
rz(-0.69951406) q[1];
sx q[1];
rz(-2.1414089) q[1];
sx q[1];
rz(-0.90739179) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9357244) q[0];
sx q[0];
rz(-2.3875934) q[0];
sx q[0];
rz(-0.058490965) q[0];
rz(-pi) q[1];
rz(-1.4110231) q[2];
sx q[2];
rz(-1.3828619) q[2];
sx q[2];
rz(2.2794276) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.34037922) q[1];
sx q[1];
rz(-2.0836909) q[1];
sx q[1];
rz(0.42162924) q[1];
x q[2];
rz(2.1069585) q[3];
sx q[3];
rz(-1.6214633) q[3];
sx q[3];
rz(-1.9127653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13713914) q[2];
sx q[2];
rz(-1.2848102) q[2];
sx q[2];
rz(2.2966906) q[2];
rz(0.6984624) q[3];
sx q[3];
rz(-0.74028492) q[3];
sx q[3];
rz(-0.79160488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7591105) q[0];
sx q[0];
rz(-1.6805205) q[0];
sx q[0];
rz(1.2680898) q[0];
rz(-1.6146487) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(-0.39438927) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0348234) q[0];
sx q[0];
rz(-0.64578694) q[0];
sx q[0];
rz(1.646498) q[0];
x q[1];
rz(-3.0494681) q[2];
sx q[2];
rz(-0.73706223) q[2];
sx q[2];
rz(1.9305522) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.39526597) q[1];
sx q[1];
rz(-1.6155287) q[1];
sx q[1];
rz(-0.6529863) q[1];
x q[2];
rz(-2.9378618) q[3];
sx q[3];
rz(-2.7739848) q[3];
sx q[3];
rz(-2.550761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4322728) q[2];
sx q[2];
rz(-3.0797112) q[2];
sx q[2];
rz(-0.56149948) q[2];
rz(2.0105441) q[3];
sx q[3];
rz(-0.80315042) q[3];
sx q[3];
rz(0.48857442) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5442218) q[0];
sx q[0];
rz(-0.18246305) q[0];
sx q[0];
rz(-0.23319787) q[0];
rz(-0.11416642) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(2.9679969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1602992) q[0];
sx q[0];
rz(-0.7795142) q[0];
sx q[0];
rz(2.1008089) q[0];
rz(-pi) q[1];
rz(-2.0162705) q[2];
sx q[2];
rz(-1.2230754) q[2];
sx q[2];
rz(1.7537774) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50362316) q[1];
sx q[1];
rz(-1.733128) q[1];
sx q[1];
rz(1.3470525) q[1];
rz(1.5830481) q[3];
sx q[3];
rz(-1.7594595) q[3];
sx q[3];
rz(2.0607299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1058098) q[2];
sx q[2];
rz(-2.1828987) q[2];
sx q[2];
rz(-2.2868273) q[2];
rz(-0.0828951) q[3];
sx q[3];
rz(-1.9708743) q[3];
sx q[3];
rz(3.0875201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4891124) q[0];
sx q[0];
rz(-1.2615477) q[0];
sx q[0];
rz(1.1789119) q[0];
rz(1.0014125) q[1];
sx q[1];
rz(-1.8779571) q[1];
sx q[1];
rz(0.15403919) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48705593) q[0];
sx q[0];
rz(-0.77870071) q[0];
sx q[0];
rz(-0.47327431) q[0];
x q[1];
rz(-0.87100864) q[2];
sx q[2];
rz(-0.99620512) q[2];
sx q[2];
rz(-2.1432997) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6045842) q[1];
sx q[1];
rz(-1.8698543) q[1];
sx q[1];
rz(-2.0589776) q[1];
rz(-0.73518153) q[3];
sx q[3];
rz(-2.7594341) q[3];
sx q[3];
rz(-0.89565403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4607294) q[2];
sx q[2];
rz(-1.0793842) q[2];
sx q[2];
rz(-0.66780773) q[2];
rz(-1.4051416) q[3];
sx q[3];
rz(-1.3677771) q[3];
sx q[3];
rz(2.5051129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619096) q[0];
sx q[0];
rz(-1.6661665) q[0];
sx q[0];
rz(-2.5478126) q[0];
rz(1.2807912) q[1];
sx q[1];
rz(-0.96071661) q[1];
sx q[1];
rz(1.2978172) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0827328) q[0];
sx q[0];
rz(-1.642859) q[0];
sx q[0];
rz(-0.58499344) q[0];
rz(1.2601398) q[2];
sx q[2];
rz(-2.1755078) q[2];
sx q[2];
rz(2.8653646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62420708) q[1];
sx q[1];
rz(-2.7605857) q[1];
sx q[1];
rz(2.7966649) q[1];
rz(-0.16420096) q[3];
sx q[3];
rz(-2.3275725) q[3];
sx q[3];
rz(-0.17539737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7456776) q[2];
sx q[2];
rz(-3.119645) q[2];
sx q[2];
rz(1.5094666) q[2];
rz(-0.11219003) q[3];
sx q[3];
rz(-1.0271007) q[3];
sx q[3];
rz(1.7621015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6701732) q[0];
sx q[0];
rz(-2.1355974) q[0];
sx q[0];
rz(-0.26563409) q[0];
rz(2.1521125) q[1];
sx q[1];
rz(-1.2629291) q[1];
sx q[1];
rz(0.33214733) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6270094) q[0];
sx q[0];
rz(-2.8081315) q[0];
sx q[0];
rz(-1.5824806) q[0];
rz(-pi) q[1];
rz(2.0004326) q[2];
sx q[2];
rz(-0.62295914) q[2];
sx q[2];
rz(2.0596383) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89236081) q[1];
sx q[1];
rz(-2.7977825) q[1];
sx q[1];
rz(-1.7147786) q[1];
rz(-pi) q[2];
rz(-3.0933558) q[3];
sx q[3];
rz(-1.7582446) q[3];
sx q[3];
rz(-1.5099883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0997448) q[2];
sx q[2];
rz(-1.0724649) q[2];
sx q[2];
rz(2.9912046) q[2];
rz(2.2930875) q[3];
sx q[3];
rz(-2.7917807) q[3];
sx q[3];
rz(-1.8457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0583508) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(-1.8104443) q[1];
sx q[1];
rz(-0.7970627) q[1];
sx q[1];
rz(2.0373559) q[1];
rz(0.065481117) q[2];
sx q[2];
rz(-0.83233287) q[2];
sx q[2];
rz(-3.0558791) q[2];
rz(2.4348197) q[3];
sx q[3];
rz(-1.5903683) q[3];
sx q[3];
rz(-0.49177468) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
