OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9159311) q[0];
sx q[0];
rz(-0.8684648) q[0];
sx q[0];
rz(2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(-2.7116099) q[1];
sx q[1];
rz(-2.4584682) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0177512) q[0];
sx q[0];
rz(-3.0525644) q[0];
sx q[0];
rz(2.428399) q[0];
x q[1];
rz(-2.1630152) q[2];
sx q[2];
rz(-0.41523146) q[2];
sx q[2];
rz(-2.4821971) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5176757) q[1];
sx q[1];
rz(-2.3530934) q[1];
sx q[1];
rz(-2.4967525) q[1];
rz(-1.0115511) q[3];
sx q[3];
rz(-2.4992001) q[3];
sx q[3];
rz(2.0410283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(2.0430298) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(1.1014972) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9803479) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(0.20794491) q[0];
rz(-0.57693276) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.4651441) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4629048) q[0];
sx q[0];
rz(-1.573274) q[0];
sx q[0];
rz(2.4173827) q[0];
rz(-pi) q[1];
rz(1.1128088) q[2];
sx q[2];
rz(-0.53456842) q[2];
sx q[2];
rz(-2.2528258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0456902) q[1];
sx q[1];
rz(-2.0369139) q[1];
sx q[1];
rz(-2.5065266) q[1];
x q[2];
rz(0.37104718) q[3];
sx q[3];
rz(-2.3608532) q[3];
sx q[3];
rz(0.77663976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0186105) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(-0.1097651) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(2.6254568) q[0];
rz(-2.5667045) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(2.3410472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0101937) q[0];
sx q[0];
rz(-0.43361615) q[0];
sx q[0];
rz(-1.2244768) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2764552) q[2];
sx q[2];
rz(-2.1102998) q[2];
sx q[2];
rz(-1.7973961) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86299455) q[1];
sx q[1];
rz(-1.9837712) q[1];
sx q[1];
rz(1.5652656) q[1];
x q[2];
rz(1.21739) q[3];
sx q[3];
rz(-0.73261515) q[3];
sx q[3];
rz(-0.64980799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4642554) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(2.1740186) q[3];
sx q[3];
rz(-1.273497) q[3];
sx q[3];
rz(-1.7165855) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490705) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(-2.676679) q[0];
rz(-0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(2.0565313) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61468716) q[0];
sx q[0];
rz(-2.5589716) q[0];
sx q[0];
rz(0.42838642) q[0];
x q[1];
rz(0.24809804) q[2];
sx q[2];
rz(-1.2144226) q[2];
sx q[2];
rz(2.6756289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9492053) q[1];
sx q[1];
rz(-2.9678223) q[1];
sx q[1];
rz(-2.5237571) q[1];
rz(-pi) q[2];
rz(1.9784847) q[3];
sx q[3];
rz(-0.99916047) q[3];
sx q[3];
rz(2.98416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(0.68112779) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(-0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27914771) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(1.7329247) q[0];
rz(-2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(0.98168215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22589332) q[0];
sx q[0];
rz(-2.3908273) q[0];
sx q[0];
rz(-2.4322926) q[0];
rz(-pi) q[1];
rz(-1.328674) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(-0.15649934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2770734) q[1];
sx q[1];
rz(-2.652087) q[1];
sx q[1];
rz(-2.8237052) q[1];
rz(-pi) q[2];
rz(-0.97335191) q[3];
sx q[3];
rz(-1.6320328) q[3];
sx q[3];
rz(0.51613584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(-1.0148467) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(-1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4846102) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(2.7691675) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(-0.16528027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0342456) q[0];
sx q[0];
rz(-0.099485569) q[0];
sx q[0];
rz(-0.29883595) q[0];
rz(0.084947305) q[2];
sx q[2];
rz(-2.4521378) q[2];
sx q[2];
rz(1.2598318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89989923) q[1];
sx q[1];
rz(-0.69679931) q[1];
sx q[1];
rz(-0.65710575) q[1];
rz(-pi) q[2];
rz(-2.9061756) q[3];
sx q[3];
rz(-1.132292) q[3];
sx q[3];
rz(-2.8188843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(1.4432663) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(-0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845602) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(2.7923287) q[0];
rz(-2.3941984) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-0.73648891) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38024494) q[0];
sx q[0];
rz(-1.5878829) q[0];
sx q[0];
rz(1.5218309) q[0];
rz(-pi) q[1];
rz(-0.94524224) q[2];
sx q[2];
rz(-2.7316299) q[2];
sx q[2];
rz(-1.2298825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0349717) q[1];
sx q[1];
rz(-1.391341) q[1];
sx q[1];
rz(-0.76645318) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6408259) q[3];
sx q[3];
rz(-1.9541249) q[3];
sx q[3];
rz(-0.70481833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(2.6100256) q[2];
rz(-0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.2980365) q[0];
rz(-0.80728665) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(2.2198026) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9471654) q[0];
sx q[0];
rz(-0.88473407) q[0];
sx q[0];
rz(-1.3961193) q[0];
rz(-pi) q[1];
rz(0.28618257) q[2];
sx q[2];
rz(-0.85368644) q[2];
sx q[2];
rz(1.66695) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7505155) q[1];
sx q[1];
rz(-0.32954307) q[1];
sx q[1];
rz(-2.4604172) q[1];
rz(-pi) q[2];
rz(-0.17955762) q[3];
sx q[3];
rz(-0.87053821) q[3];
sx q[3];
rz(1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(-1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(-1.4484891) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885324) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.4755479) q[0];
rz(1.6015923) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(-2.4618861) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98201671) q[0];
sx q[0];
rz(-1.3681108) q[0];
sx q[0];
rz(-2.6152339) q[0];
rz(1.7858511) q[2];
sx q[2];
rz(-0.70677033) q[2];
sx q[2];
rz(3.0422473) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6913773) q[1];
sx q[1];
rz(-2.2274349) q[1];
sx q[1];
rz(-1.1378098) q[1];
rz(-0.9548095) q[3];
sx q[3];
rz(-0.98815742) q[3];
sx q[3];
rz(-0.59347502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0891721) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(-0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.9649327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3424073) q[0];
sx q[0];
rz(-0.17051324) q[0];
sx q[0];
rz(-1.8605581) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5942957) q[2];
sx q[2];
rz(-1.0101057) q[2];
sx q[2];
rz(-1.0123569) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0050051) q[1];
sx q[1];
rz(-2.1274381) q[1];
sx q[1];
rz(0.088076061) q[1];
rz(-pi) q[2];
rz(-1.482974) q[3];
sx q[3];
rz(-1.0613958) q[3];
sx q[3];
rz(-1.0305962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-2.0142377) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(-2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42416278) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(0.39623109) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(0.73666112) q[2];
sx q[2];
rz(-1.7792637) q[2];
sx q[2];
rz(-1.7532495) q[2];
rz(-0.35076326) q[3];
sx q[3];
rz(-1.5784932) q[3];
sx q[3];
rz(-2.2898883) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];