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
rz(-1.3136366) q[0];
sx q[0];
rz(-3.0744636) q[0];
sx q[0];
rz(-0.022775291) q[0];
rz(1.0442806) q[1];
sx q[1];
rz(-2.2106946) q[1];
sx q[1];
rz(-3.1276303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3168105) q[0];
sx q[0];
rz(-2.383962) q[0];
sx q[0];
rz(0.05912058) q[0];
rz(-pi) q[1];
rz(-0.44273744) q[2];
sx q[2];
rz(-1.95089) q[2];
sx q[2];
rz(-0.9600823) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8163482) q[1];
sx q[1];
rz(-0.26953408) q[1];
sx q[1];
rz(-2.4052038) q[1];
rz(2.4926643) q[3];
sx q[3];
rz(-2.5960659) q[3];
sx q[3];
rz(0.76577556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17249168) q[2];
sx q[2];
rz(-1.3250985) q[2];
sx q[2];
rz(1.3410404) q[2];
rz(1.7155044) q[3];
sx q[3];
rz(-2.0243578) q[3];
sx q[3];
rz(0.30556998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7476927) q[0];
sx q[0];
rz(-1.9992398) q[0];
sx q[0];
rz(0.72545141) q[0];
rz(0.91616383) q[1];
sx q[1];
rz(-2.3326645) q[1];
sx q[1];
rz(-0.61942548) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325798) q[0];
sx q[0];
rz(-2.469099) q[0];
sx q[0];
rz(2.2605863) q[0];
rz(0.97255637) q[2];
sx q[2];
rz(-2.2248039) q[2];
sx q[2];
rz(2.0218265) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5350772) q[1];
sx q[1];
rz(-1.1454574) q[1];
sx q[1];
rz(0.40014793) q[1];
rz(-pi) q[2];
rz(-2.8579669) q[3];
sx q[3];
rz(-1.1414293) q[3];
sx q[3];
rz(-0.46653433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3118423) q[2];
sx q[2];
rz(-1.0210911) q[2];
sx q[2];
rz(-0.96168438) q[2];
rz(2.035615) q[3];
sx q[3];
rz(-1.032369) q[3];
sx q[3];
rz(2.7963514) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2894534) q[0];
sx q[0];
rz(-1.1281321) q[0];
sx q[0];
rz(-1.3683251) q[0];
rz(-0.013956919) q[1];
sx q[1];
rz(-2.0266666) q[1];
sx q[1];
rz(-1.4143292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.208983) q[0];
sx q[0];
rz(-1.6609207) q[0];
sx q[0];
rz(1.8831253) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3886202) q[2];
sx q[2];
rz(-0.37473703) q[2];
sx q[2];
rz(-0.47606766) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5217375) q[1];
sx q[1];
rz(-2.450469) q[1];
sx q[1];
rz(-0.85994263) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8295123) q[3];
sx q[3];
rz(-1.95041) q[3];
sx q[3];
rz(-2.1624452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76325512) q[2];
sx q[2];
rz(-1.1474643) q[2];
sx q[2];
rz(-0.21667996) q[2];
rz(0.8832776) q[3];
sx q[3];
rz(-2.7933385) q[3];
sx q[3];
rz(2.0047552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388409) q[0];
sx q[0];
rz(-2.1348248) q[0];
sx q[0];
rz(0.83813611) q[0];
rz(-2.5996767) q[1];
sx q[1];
rz(-0.37626615) q[1];
sx q[1];
rz(-0.045104973) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7823001) q[0];
sx q[0];
rz(-1.703039) q[0];
sx q[0];
rz(2.537893) q[0];
x q[1];
rz(-2.0531282) q[2];
sx q[2];
rz(-1.3399897) q[2];
sx q[2];
rz(0.34462356) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1029658) q[1];
sx q[1];
rz(-0.47768394) q[1];
sx q[1];
rz(1.4250907) q[1];
x q[2];
rz(-0.49398936) q[3];
sx q[3];
rz(-1.6574548) q[3];
sx q[3];
rz(-2.983421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0185467) q[2];
sx q[2];
rz(-2.9728153) q[2];
sx q[2];
rz(0.65400845) q[2];
rz(0.6319913) q[3];
sx q[3];
rz(-1.4258823) q[3];
sx q[3];
rz(-0.33647195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80394799) q[0];
sx q[0];
rz(-0.73886442) q[0];
sx q[0];
rz(1.3496572) q[0];
rz(0.82756394) q[1];
sx q[1];
rz(-0.63500985) q[1];
sx q[1];
rz(0.48447022) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1903828) q[0];
sx q[0];
rz(-1.2408549) q[0];
sx q[0];
rz(-1.1879107) q[0];
rz(-pi) q[1];
rz(1.7155649) q[2];
sx q[2];
rz(-0.55895222) q[2];
sx q[2];
rz(-1.1048855) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.81329359) q[1];
sx q[1];
rz(-0.44078953) q[1];
sx q[1];
rz(-0.30118024) q[1];
rz(-pi) q[2];
rz(0.46273939) q[3];
sx q[3];
rz(-0.92904186) q[3];
sx q[3];
rz(2.3433507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95939246) q[2];
sx q[2];
rz(-2.0727938) q[2];
sx q[2];
rz(2.5051795) q[2];
rz(-0.098879769) q[3];
sx q[3];
rz(-1.5560047) q[3];
sx q[3];
rz(-1.4720565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32475489) q[0];
sx q[0];
rz(-0.88961283) q[0];
sx q[0];
rz(-1.1725934) q[0];
rz(1.6088387) q[1];
sx q[1];
rz(-2.1295261) q[1];
sx q[1];
rz(3.1408659) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4945723) q[0];
sx q[0];
rz(-1.9969517) q[0];
sx q[0];
rz(-1.0264574) q[0];
x q[1];
rz(0.41623314) q[2];
sx q[2];
rz(-2.9011167) q[2];
sx q[2];
rz(-1.7822303) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3657942) q[1];
sx q[1];
rz(-1.9025814) q[1];
sx q[1];
rz(1.6901687) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57925333) q[3];
sx q[3];
rz(-1.4046852) q[3];
sx q[3];
rz(-0.36234713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26649228) q[2];
sx q[2];
rz(-2.2537587) q[2];
sx q[2];
rz(-0.20981851) q[2];
rz(-0.80875129) q[3];
sx q[3];
rz(-0.54976141) q[3];
sx q[3];
rz(-0.89890283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631735) q[0];
sx q[0];
rz(-1.2043948) q[0];
sx q[0];
rz(-0.11209442) q[0];
rz(-0.049292715) q[1];
sx q[1];
rz(-0.53105989) q[1];
sx q[1];
rz(1.8045527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1111832) q[0];
sx q[0];
rz(-2.2044535) q[0];
sx q[0];
rz(-0.2588057) q[0];
rz(-pi) q[1];
rz(0.12814204) q[2];
sx q[2];
rz(-2.0754092) q[2];
sx q[2];
rz(-0.59282936) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35838764) q[1];
sx q[1];
rz(-1.9175944) q[1];
sx q[1];
rz(1.5601182) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1483809) q[3];
sx q[3];
rz(-1.3903729) q[3];
sx q[3];
rz(2.5374295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8423395) q[2];
sx q[2];
rz(-0.82483333) q[2];
sx q[2];
rz(-2.4315368) q[2];
rz(-3.1402785) q[3];
sx q[3];
rz(-0.90130663) q[3];
sx q[3];
rz(0.60091758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8173219) q[0];
sx q[0];
rz(-1.7621499) q[0];
sx q[0];
rz(2.5589909) q[0];
rz(2.461589) q[1];
sx q[1];
rz(-0.99907196) q[1];
sx q[1];
rz(1.2635788) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.82915) q[0];
sx q[0];
rz(-1.1151322) q[0];
sx q[0];
rz(-1.2553041) q[0];
rz(-pi) q[1];
rz(-2.3266257) q[2];
sx q[2];
rz(-0.294058) q[2];
sx q[2];
rz(2.5240099) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78681011) q[1];
sx q[1];
rz(-2.218333) q[1];
sx q[1];
rz(-1.3446273) q[1];
rz(2.7177592) q[3];
sx q[3];
rz(-1.6006216) q[3];
sx q[3];
rz(-0.78889293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6363643) q[2];
sx q[2];
rz(-1.6978846) q[2];
sx q[2];
rz(0.12602028) q[2];
rz(1.5382918) q[3];
sx q[3];
rz(-2.3647629) q[3];
sx q[3];
rz(-0.34415027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4890323) q[0];
sx q[0];
rz(-1.5131938) q[0];
sx q[0];
rz(-1.162758) q[0];
rz(1.3105505) q[1];
sx q[1];
rz(-2.4030011) q[1];
sx q[1];
rz(-1.7576677) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2306128) q[0];
sx q[0];
rz(-0.34251102) q[0];
sx q[0];
rz(-1.4369276) q[0];
rz(0.12284367) q[2];
sx q[2];
rz(-1.4990988) q[2];
sx q[2];
rz(-0.724585) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3466665) q[1];
sx q[1];
rz(-0.62267471) q[1];
sx q[1];
rz(-1.7706448) q[1];
rz(-1.7095164) q[3];
sx q[3];
rz(-1.6131716) q[3];
sx q[3];
rz(0.78813533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80013529) q[2];
sx q[2];
rz(-1.2722445) q[2];
sx q[2];
rz(2.6696491) q[2];
rz(2.3927472) q[3];
sx q[3];
rz(-0.92350903) q[3];
sx q[3];
rz(-2.3239465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6227459) q[0];
sx q[0];
rz(-2.2846344) q[0];
sx q[0];
rz(-2.6527606) q[0];
rz(2.7986616) q[1];
sx q[1];
rz(-0.88645005) q[1];
sx q[1];
rz(1.0541281) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3341772) q[0];
sx q[0];
rz(-1.5332959) q[0];
sx q[0];
rz(1.6223909) q[0];
x q[1];
rz(-1.9185877) q[2];
sx q[2];
rz(-2.1556097) q[2];
sx q[2];
rz(1.0725759) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2239687) q[1];
sx q[1];
rz(-1.2792933) q[1];
sx q[1];
rz(-1.0179917) q[1];
rz(-1.6659237) q[3];
sx q[3];
rz(-0.47387487) q[3];
sx q[3];
rz(-2.9364862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9752581) q[2];
sx q[2];
rz(-2.9322093) q[2];
sx q[2];
rz(0.091910467) q[2];
rz(-0.46936834) q[3];
sx q[3];
rz(-0.86463237) q[3];
sx q[3];
rz(-3.0288467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54409201) q[0];
sx q[0];
rz(-1.8202029) q[0];
sx q[0];
rz(1.2910917) q[0];
rz(-0.76580936) q[1];
sx q[1];
rz(-1.7821093) q[1];
sx q[1];
rz(-1.2761188) q[1];
rz(2.1738293) q[2];
sx q[2];
rz(-2.7104678) q[2];
sx q[2];
rz(-0.34061876) q[2];
rz(0.41894434) q[3];
sx q[3];
rz(-0.57804116) q[3];
sx q[3];
rz(-3.0347435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
