OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6717186) q[0];
sx q[0];
rz(-1.1550386) q[0];
sx q[0];
rz(-1.7116829) q[0];
rz(0.45541304) q[1];
sx q[1];
rz(5.6918511) q[1];
sx q[1];
rz(11.439352) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9874632) q[0];
sx q[0];
rz(-2.9107981) q[0];
sx q[0];
rz(1.5202281) q[0];
rz(-pi) q[1];
rz(1.3052218) q[2];
sx q[2];
rz(-1.2392534) q[2];
sx q[2];
rz(-1.3217373) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6851122) q[1];
sx q[1];
rz(-0.64549533) q[1];
sx q[1];
rz(1.0935893) q[1];
x q[2];
rz(-1.5886878) q[3];
sx q[3];
rz(-1.761529) q[3];
sx q[3];
rz(-2.408556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0226125) q[2];
sx q[2];
rz(-2.4369414) q[2];
sx q[2];
rz(-1.653778) q[2];
rz(2.7704499) q[3];
sx q[3];
rz(-1.1266339) q[3];
sx q[3];
rz(-0.77905542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356075) q[0];
sx q[0];
rz(-1.7636517) q[0];
sx q[0];
rz(-0.66407472) q[0];
rz(-2.5750419) q[1];
sx q[1];
rz(-1.5357176) q[1];
sx q[1];
rz(2.9552592) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.60255) q[0];
sx q[0];
rz(-1.7930174) q[0];
sx q[0];
rz(-2.1150388) q[0];
x q[1];
rz(2.8736224) q[2];
sx q[2];
rz(-1.2023965) q[2];
sx q[2];
rz(0.024193833) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1809606) q[1];
sx q[1];
rz(-0.69968984) q[1];
sx q[1];
rz(1.3153428) q[1];
rz(-1.4545031) q[3];
sx q[3];
rz(-2.6175559) q[3];
sx q[3];
rz(0.91479036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32249195) q[2];
sx q[2];
rz(-1.2995316) q[2];
sx q[2];
rz(1.5619649) q[2];
rz(1.1473848) q[3];
sx q[3];
rz(-2.8473144) q[3];
sx q[3];
rz(-1.3636205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0550302) q[0];
sx q[0];
rz(-1.6494305) q[0];
sx q[0];
rz(-0.30558875) q[0];
rz(-1.2809523) q[1];
sx q[1];
rz(-0.31550229) q[1];
sx q[1];
rz(-1.618128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9394835) q[0];
sx q[0];
rz(-1.8022707) q[0];
sx q[0];
rz(2.7702296) q[0];
rz(-pi) q[1];
rz(-1.5554713) q[2];
sx q[2];
rz(-1.5256679) q[2];
sx q[2];
rz(-0.082782291) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9429607) q[1];
sx q[1];
rz(-0.80201282) q[1];
sx q[1];
rz(2.632891) q[1];
rz(2.9039246) q[3];
sx q[3];
rz(-0.97252995) q[3];
sx q[3];
rz(3.0868343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4898701) q[2];
sx q[2];
rz(-1.2469331) q[2];
sx q[2];
rz(0.22001246) q[2];
rz(-3.0618727) q[3];
sx q[3];
rz(-0.97697512) q[3];
sx q[3];
rz(-0.93130934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43059573) q[0];
sx q[0];
rz(-1.3871223) q[0];
sx q[0];
rz(-0.0084477607) q[0];
rz(-1.9795817) q[1];
sx q[1];
rz(-1.9879257) q[1];
sx q[1];
rz(2.0268424) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8793068) q[0];
sx q[0];
rz(-1.382785) q[0];
sx q[0];
rz(-1.2557058) q[0];
rz(-pi) q[1];
rz(-0.77728072) q[2];
sx q[2];
rz(-2.4156164) q[2];
sx q[2];
rz(-0.41997806) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48753187) q[1];
sx q[1];
rz(-1.2966074) q[1];
sx q[1];
rz(1.2031005) q[1];
rz(-pi) q[2];
rz(-0.42041789) q[3];
sx q[3];
rz(-1.9501996) q[3];
sx q[3];
rz(-2.8914352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.70725012) q[2];
sx q[2];
rz(-1.0184526) q[2];
sx q[2];
rz(-2.2960091) q[2];
rz(-0.25104684) q[3];
sx q[3];
rz(-1.8699402) q[3];
sx q[3];
rz(2.3959851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4276328) q[0];
sx q[0];
rz(-2.7290955) q[0];
sx q[0];
rz(-2.1115671) q[0];
rz(0.59459844) q[1];
sx q[1];
rz(-1.7183036) q[1];
sx q[1];
rz(1.3535708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0145688) q[0];
sx q[0];
rz(-1.5480124) q[0];
sx q[0];
rz(-1.5500853) q[0];
rz(-0.32234742) q[2];
sx q[2];
rz(-1.3953096) q[2];
sx q[2];
rz(-0.86514003) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7723778) q[1];
sx q[1];
rz(-0.90638834) q[1];
sx q[1];
rz(2.496705) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5491675) q[3];
sx q[3];
rz(-1.1325784) q[3];
sx q[3];
rz(-2.2074204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2305962) q[2];
sx q[2];
rz(-1.5065008) q[2];
sx q[2];
rz(0.050749151) q[2];
rz(-2.0283608) q[3];
sx q[3];
rz(-1.9751996) q[3];
sx q[3];
rz(-0.39065233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0473061) q[0];
sx q[0];
rz(-2.081649) q[0];
sx q[0];
rz(-0.37288368) q[0];
rz(-1.8376384) q[1];
sx q[1];
rz(-1.8770437) q[1];
sx q[1];
rz(2.1162927) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.404843) q[0];
sx q[0];
rz(-1.4944634) q[0];
sx q[0];
rz(2.6910604) q[0];
rz(-0.16496678) q[2];
sx q[2];
rz(-1.8968762) q[2];
sx q[2];
rz(3.0086089) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.655788) q[1];
sx q[1];
rz(-0.88690573) q[1];
sx q[1];
rz(-0.25836103) q[1];
x q[2];
rz(2.9168374) q[3];
sx q[3];
rz(-1.9850176) q[3];
sx q[3];
rz(0.97771588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.23124) q[2];
sx q[2];
rz(-0.19008907) q[2];
sx q[2];
rz(-1.4300038) q[2];
rz(1.5405687) q[3];
sx q[3];
rz(-0.68683306) q[3];
sx q[3];
rz(-1.6704667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98825276) q[0];
sx q[0];
rz(-1.3000458) q[0];
sx q[0];
rz(2.7100995) q[0];
rz(-1.2639812) q[1];
sx q[1];
rz(-1.6004205) q[1];
sx q[1];
rz(-2.4914609) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4503195) q[0];
sx q[0];
rz(-2.8421845) q[0];
sx q[0];
rz(-2.8595238) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.442914) q[2];
sx q[2];
rz(-0.68727109) q[2];
sx q[2];
rz(0.56035794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3172315) q[1];
sx q[1];
rz(-1.6859833) q[1];
sx q[1];
rz(-2.2735866) q[1];
x q[2];
rz(2.0286719) q[3];
sx q[3];
rz(-1.5639389) q[3];
sx q[3];
rz(-0.53660652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5905137) q[2];
sx q[2];
rz(-1.7503909) q[2];
sx q[2];
rz(-3.1371269) q[2];
rz(3.1320069) q[3];
sx q[3];
rz(-1.3349345) q[3];
sx q[3];
rz(2.7981304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9850605) q[0];
sx q[0];
rz(-0.28230202) q[0];
sx q[0];
rz(0.16047934) q[0];
rz(1.7566682) q[1];
sx q[1];
rz(-2.3424708) q[1];
sx q[1];
rz(2.8327732) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8978116) q[0];
sx q[0];
rz(-1.7574991) q[0];
sx q[0];
rz(-0.12630442) q[0];
x q[1];
rz(0.43833959) q[2];
sx q[2];
rz(-2.380854) q[2];
sx q[2];
rz(2.1957514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.908573) q[1];
sx q[1];
rz(-1.9601788) q[1];
sx q[1];
rz(-1.977747) q[1];
x q[2];
rz(2.4229523) q[3];
sx q[3];
rz(-2.6860533) q[3];
sx q[3];
rz(-2.9589341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.34362346) q[2];
sx q[2];
rz(-0.60911959) q[2];
sx q[2];
rz(1.5416175) q[2];
rz(0.68495098) q[3];
sx q[3];
rz(-1.5870321) q[3];
sx q[3];
rz(1.5233013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5110382) q[0];
sx q[0];
rz(-1.7031952) q[0];
sx q[0];
rz(2.5901929) q[0];
rz(-2.664387) q[1];
sx q[1];
rz(-2.2255032) q[1];
sx q[1];
rz(-0.83470693) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4906368) q[0];
sx q[0];
rz(-1.0792562) q[0];
sx q[0];
rz(-2.0975344) q[0];
rz(0.70087011) q[2];
sx q[2];
rz(-1.5972023) q[2];
sx q[2];
rz(-0.7757265) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.44025035) q[1];
sx q[1];
rz(-1.6052488) q[1];
sx q[1];
rz(2.1828096) q[1];
x q[2];
rz(-0.61309321) q[3];
sx q[3];
rz(-1.0524155) q[3];
sx q[3];
rz(-1.6676774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7029552) q[2];
sx q[2];
rz(-0.96200395) q[2];
sx q[2];
rz(0.23923242) q[2];
rz(2.8171825) q[3];
sx q[3];
rz(-1.4374461) q[3];
sx q[3];
rz(-1.602406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6259554) q[0];
sx q[0];
rz(-0.13787585) q[0];
sx q[0];
rz(2.4249518) q[0];
rz(0.58249885) q[1];
sx q[1];
rz(-1.3526252) q[1];
sx q[1];
rz(-2.8315721) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9649565) q[0];
sx q[0];
rz(-0.88439098) q[0];
sx q[0];
rz(-0.45748894) q[0];
rz(-pi) q[1];
rz(2.7298195) q[2];
sx q[2];
rz(-1.207282) q[2];
sx q[2];
rz(-0.95136729) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2523732) q[1];
sx q[1];
rz(-2.3205726) q[1];
sx q[1];
rz(-2.0682425) q[1];
rz(-pi) q[2];
rz(-2.4118092) q[3];
sx q[3];
rz(-2.457629) q[3];
sx q[3];
rz(-1.7548949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6325355) q[2];
sx q[2];
rz(-2.5639503) q[2];
sx q[2];
rz(-1.7217815) q[2];
rz(-1.696473) q[3];
sx q[3];
rz(-2.4308379) q[3];
sx q[3];
rz(-2.3783309) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2390908) q[0];
sx q[0];
rz(-0.9466753) q[0];
sx q[0];
rz(-1.7539903) q[0];
rz(-1.6613962) q[1];
sx q[1];
rz(-1.4085242) q[1];
sx q[1];
rz(-2.9396802) q[1];
rz(0.26915941) q[2];
sx q[2];
rz(-2.8648389) q[2];
sx q[2];
rz(1.7493389) q[2];
rz(2.673951) q[3];
sx q[3];
rz(-0.94398879) q[3];
sx q[3];
rz(2.332609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
