OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93513918) q[0];
sx q[0];
rz(-2.3606665) q[0];
sx q[0];
rz(0.20679064) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(1.0607399) q[1];
sx q[1];
rz(12.386204) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59158303) q[0];
sx q[0];
rz(-1.9073434) q[0];
sx q[0];
rz(0.49206375) q[0];
rz(1.7018868) q[2];
sx q[2];
rz(-2.0799473) q[2];
sx q[2];
rz(2.5703562) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3611316) q[1];
sx q[1];
rz(-1.6447004) q[1];
sx q[1];
rz(0.88688811) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3869242) q[3];
sx q[3];
rz(-2.2443218) q[3];
sx q[3];
rz(-2.9415188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7303598) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(1.2940548) q[2];
rz(-0.40575746) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(-0.40677235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2186573) q[0];
sx q[0];
rz(-1.0359534) q[0];
sx q[0];
rz(3.0157715) q[0];
rz(0.80548349) q[1];
sx q[1];
rz(-2.3352354) q[1];
sx q[1];
rz(1.3719826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4683511) q[0];
sx q[0];
rz(-1.5090319) q[0];
sx q[0];
rz(-0.67586918) q[0];
x q[1];
rz(-0.54601045) q[2];
sx q[2];
rz(-1.6766607) q[2];
sx q[2];
rz(0.0042303483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0731197) q[1];
sx q[1];
rz(-1.321055) q[1];
sx q[1];
rz(1.5239034) q[1];
x q[2];
rz(-2.6847141) q[3];
sx q[3];
rz(-0.51364726) q[3];
sx q[3];
rz(1.9809686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8877318) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(0.99622336) q[2];
rz(2.0455202) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(2.2912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(-1.0590142) q[0];
rz(1.1478708) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(-2.0770729) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7721467) q[0];
sx q[0];
rz(-1.9659974) q[0];
sx q[0];
rz(-0.69338436) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1883005) q[2];
sx q[2];
rz(-1.7171535) q[2];
sx q[2];
rz(1.0548897) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3435622) q[1];
sx q[1];
rz(-0.31032944) q[1];
sx q[1];
rz(1.6088435) q[1];
x q[2];
rz(-0.84633175) q[3];
sx q[3];
rz(-1.6958106) q[3];
sx q[3];
rz(-0.66305977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.069783) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(-0.26322571) q[2];
rz(-2.0227382) q[3];
sx q[3];
rz(-1.8656732) q[3];
sx q[3];
rz(0.81937218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333106) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(-1.7472965) q[0];
rz(-1.0385723) q[1];
sx q[1];
rz(-1.4515406) q[1];
sx q[1];
rz(1.0669473) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79558668) q[0];
sx q[0];
rz(-2.6104402) q[0];
sx q[0];
rz(1.4941494) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8550451) q[2];
sx q[2];
rz(-0.28038014) q[2];
sx q[2];
rz(-2.7718411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3921515) q[1];
sx q[1];
rz(-1.1679808) q[1];
sx q[1];
rz(-0.72026003) q[1];
rz(-2.9158343) q[3];
sx q[3];
rz(-2.0421713) q[3];
sx q[3];
rz(0.86172047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2944494) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(2.8520544) q[2];
rz(0.52044049) q[3];
sx q[3];
rz(-1.3402904) q[3];
sx q[3];
rz(2.382544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4905869) q[0];
sx q[0];
rz(-3.1327972) q[0];
sx q[0];
rz(2.3068413) q[0];
rz(1.2443776) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(1.7117737) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6288888) q[0];
sx q[0];
rz(-1.5682724) q[0];
sx q[0];
rz(3.0789037) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69447563) q[2];
sx q[2];
rz(-0.36437964) q[2];
sx q[2];
rz(-2.7642872) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5809708) q[1];
sx q[1];
rz(-2.3498658) q[1];
sx q[1];
rz(-3.1091299) q[1];
x q[2];
rz(2.2291282) q[3];
sx q[3];
rz(-1.4294799) q[3];
sx q[3];
rz(-0.055012881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4804068) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(1.8583813) q[2];
rz(-3.0094106) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(2.8018518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(2.6640889) q[0];
rz(1.6409138) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(0.57055155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30734277) q[0];
sx q[0];
rz(-1.2120005) q[0];
sx q[0];
rz(-0.3105727) q[0];
x q[1];
rz(1.19679) q[2];
sx q[2];
rz(-1.1454957) q[2];
sx q[2];
rz(2.5380295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.693336) q[1];
sx q[1];
rz(-2.2168471) q[1];
sx q[1];
rz(-0.64530428) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1122094) q[3];
sx q[3];
rz(-2.4930862) q[3];
sx q[3];
rz(2.9103968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70696124) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(2.3449507) q[2];
rz(-0.32026511) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(-1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758078) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(0.35807034) q[0];
rz(0.2886731) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(-1.8310865) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7983127) q[0];
sx q[0];
rz(-0.57050059) q[0];
sx q[0];
rz(0.35240726) q[0];
x q[1];
rz(0.99408044) q[2];
sx q[2];
rz(-0.95816441) q[2];
sx q[2];
rz(-2.6660369) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9999034) q[1];
sx q[1];
rz(-2.018376) q[1];
sx q[1];
rz(-0.56772851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2745423) q[3];
sx q[3];
rz(-0.82074814) q[3];
sx q[3];
rz(-1.5914608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8075809) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(0.99747783) q[2];
rz(-2.5618662) q[3];
sx q[3];
rz(-0.72967356) q[3];
sx q[3];
rz(1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26577935) q[0];
sx q[0];
rz(-1.4900603) q[0];
sx q[0];
rz(0.34061256) q[0];
rz(1.2365201) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(1.9514726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53100454) q[0];
sx q[0];
rz(-1.6325145) q[0];
sx q[0];
rz(1.8025663) q[0];
rz(-pi) q[1];
rz(2.2731407) q[2];
sx q[2];
rz(-3.0358899) q[2];
sx q[2];
rz(-2.5603574) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5684143) q[1];
sx q[1];
rz(-1.2866932) q[1];
sx q[1];
rz(1.2969639) q[1];
x q[2];
rz(2.1281151) q[3];
sx q[3];
rz(-1.1790457) q[3];
sx q[3];
rz(1.980892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74165806) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(2.7271872) q[2];
rz(-1.7587781) q[3];
sx q[3];
rz(-1.7410991) q[3];
sx q[3];
rz(-2.8167021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25093108) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(2.9898306) q[0];
rz(-1.4258619) q[1];
sx q[1];
rz(-1.3769923) q[1];
sx q[1];
rz(-2.192416) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5715036) q[0];
sx q[0];
rz(-0.46496689) q[0];
sx q[0];
rz(-0.31682195) q[0];
x q[1];
rz(0.29813913) q[2];
sx q[2];
rz(-2.4900644) q[2];
sx q[2];
rz(2.8704314) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6712499) q[1];
sx q[1];
rz(-1.8103231) q[1];
sx q[1];
rz(0.67269477) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2373958) q[3];
sx q[3];
rz(-0.70611806) q[3];
sx q[3];
rz(-2.0437984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(0.43241832) q[2];
rz(1.6010823) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0703053) q[0];
sx q[0];
rz(-0.27305958) q[0];
sx q[0];
rz(0.37316698) q[0];
rz(-0.35692731) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(2.4180791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05941885) q[0];
sx q[0];
rz(-0.84616236) q[0];
sx q[0];
rz(-1.3979777) q[0];
rz(-0.11156545) q[2];
sx q[2];
rz(-2.4590883) q[2];
sx q[2];
rz(-1.1139368) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8404624) q[1];
sx q[1];
rz(-2.3311619) q[1];
sx q[1];
rz(0.37709548) q[1];
rz(-0.27042737) q[3];
sx q[3];
rz(-1.5661097) q[3];
sx q[3];
rz(0.0029431012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2202806) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(0.55345654) q[2];
rz(0.71183318) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.9217459) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(0.28221054) q[1];
sx q[1];
rz(-0.78411513) q[1];
sx q[1];
rz(0.93906739) q[1];
rz(2.0201335) q[2];
sx q[2];
rz(-1.1171787) q[2];
sx q[2];
rz(1.6906307) q[2];
rz(-0.40209963) q[3];
sx q[3];
rz(-1.3105018) q[3];
sx q[3];
rz(-0.49285938) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
