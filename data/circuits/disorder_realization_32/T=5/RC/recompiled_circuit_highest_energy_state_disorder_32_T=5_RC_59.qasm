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
rz(-0.028539874) q[0];
sx q[0];
rz(3.9588504) q[0];
sx q[0];
rz(9.1398653) q[0];
rz(0.81975308) q[1];
sx q[1];
rz(-0.88480359) q[1];
sx q[1];
rz(-1.1362145) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4856271) q[0];
sx q[0];
rz(-2.1871217) q[0];
sx q[0];
rz(0.095353145) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66464632) q[2];
sx q[2];
rz(-0.79865361) q[2];
sx q[2];
rz(0.89671521) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2023102) q[1];
sx q[1];
rz(-1.9946672) q[1];
sx q[1];
rz(1.8821299) q[1];
rz(-pi) q[2];
rz(2.6873246) q[3];
sx q[3];
rz(-0.48503808) q[3];
sx q[3];
rz(0.046864101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39584407) q[2];
sx q[2];
rz(-2.6555847) q[2];
sx q[2];
rz(-2.8493472) q[2];
rz(-2.5203943) q[3];
sx q[3];
rz(-2.0665593) q[3];
sx q[3];
rz(0.21875374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.578823) q[0];
sx q[0];
rz(-2.0243473) q[0];
sx q[0];
rz(0.086409464) q[0];
rz(-0.39331618) q[1];
sx q[1];
rz(-1.4540648) q[1];
sx q[1];
rz(-1.4837861) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4230712) q[0];
sx q[0];
rz(-2.4212061) q[0];
sx q[0];
rz(1.2116526) q[0];
rz(1.2582766) q[2];
sx q[2];
rz(-1.7023689) q[2];
sx q[2];
rz(-0.67503345) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.68951974) q[1];
sx q[1];
rz(-1.7240925) q[1];
sx q[1];
rz(0.91233715) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3065763) q[3];
sx q[3];
rz(-1.0162899) q[3];
sx q[3];
rz(-2.4305024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1958486) q[2];
sx q[2];
rz(-1.5986779) q[2];
sx q[2];
rz(-3.1405385) q[2];
rz(2.4347608) q[3];
sx q[3];
rz(-0.89444923) q[3];
sx q[3];
rz(2.8640174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23866776) q[0];
sx q[0];
rz(-2.213573) q[0];
sx q[0];
rz(0.4784041) q[0];
rz(-0.84486419) q[1];
sx q[1];
rz(-1.076315) q[1];
sx q[1];
rz(2.1873651) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85027116) q[0];
sx q[0];
rz(-2.5370165) q[0];
sx q[0];
rz(-0.69285562) q[0];
x q[1];
rz(2.235844) q[2];
sx q[2];
rz(-2.2881788) q[2];
sx q[2];
rz(-0.61418542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0817928) q[1];
sx q[1];
rz(-0.89443086) q[1];
sx q[1];
rz(0.46251071) q[1];
rz(-1.30458) q[3];
sx q[3];
rz(-2.2053218) q[3];
sx q[3];
rz(-2.3887544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2039631) q[2];
sx q[2];
rz(-0.340168) q[2];
sx q[2];
rz(-3.0703239) q[2];
rz(2.1313306) q[3];
sx q[3];
rz(-1.1610169) q[3];
sx q[3];
rz(0.32101139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45692745) q[0];
sx q[0];
rz(-1.2926956) q[0];
sx q[0];
rz(0.36351031) q[0];
rz(2.2495296) q[1];
sx q[1];
rz(-1.0328247) q[1];
sx q[1];
rz(1.8878638) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9867986) q[0];
sx q[0];
rz(-0.61756931) q[0];
sx q[0];
rz(-0.85394959) q[0];
rz(-1.2066098) q[2];
sx q[2];
rz(-1.4763583) q[2];
sx q[2];
rz(-2.0409453) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82361932) q[1];
sx q[1];
rz(-0.25356217) q[1];
sx q[1];
rz(0.16615725) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4917077) q[3];
sx q[3];
rz(-2.657848) q[3];
sx q[3];
rz(-0.080881491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7518647) q[2];
sx q[2];
rz(-2.3675297) q[2];
sx q[2];
rz(-1.6130201) q[2];
rz(-0.45016897) q[3];
sx q[3];
rz(-1.3521786) q[3];
sx q[3];
rz(-1.842513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12440974) q[0];
sx q[0];
rz(-0.15394177) q[0];
sx q[0];
rz(2.2907139) q[0];
rz(1.7393913) q[1];
sx q[1];
rz(-1.5504928) q[1];
sx q[1];
rz(-0.15241399) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7895592) q[0];
sx q[0];
rz(-1.5646345) q[0];
sx q[0];
rz(1.3314817) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28976243) q[2];
sx q[2];
rz(-2.1524977) q[2];
sx q[2];
rz(-1.5982472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85665515) q[1];
sx q[1];
rz(-2.8116715) q[1];
sx q[1];
rz(1.0367212) q[1];
x q[2];
rz(-1.2021303) q[3];
sx q[3];
rz(-1.1970425) q[3];
sx q[3];
rz(2.8522648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.57205814) q[2];
sx q[2];
rz(-0.92430884) q[2];
sx q[2];
rz(2.7657236) q[2];
rz(0.26840261) q[3];
sx q[3];
rz(-2.1208051) q[3];
sx q[3];
rz(1.0387897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1035476) q[0];
sx q[0];
rz(-0.91829848) q[0];
sx q[0];
rz(3.015836) q[0];
rz(-3.0962931) q[1];
sx q[1];
rz(-2.2442975) q[1];
sx q[1];
rz(2.6364141) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4046458) q[0];
sx q[0];
rz(-0.26261273) q[0];
sx q[0];
rz(2.8960431) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0313377) q[2];
sx q[2];
rz(-1.0822191) q[2];
sx q[2];
rz(-1.1854805) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4498429) q[1];
sx q[1];
rz(-2.3556752) q[1];
sx q[1];
rz(-0.20699006) q[1];
rz(-pi) q[2];
rz(-2.1734851) q[3];
sx q[3];
rz(-0.54881964) q[3];
sx q[3];
rz(-0.88879648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3230285) q[2];
sx q[2];
rz(-0.33107859) q[2];
sx q[2];
rz(2.8884812) q[2];
rz(-2.3969635) q[3];
sx q[3];
rz(-1.5537477) q[3];
sx q[3];
rz(-2.871992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.44148663) q[0];
sx q[0];
rz(-1.8550669) q[0];
sx q[0];
rz(-0.043524608) q[0];
rz(1.173191) q[1];
sx q[1];
rz(-1.2930608) q[1];
sx q[1];
rz(2.3399369) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7907958) q[0];
sx q[0];
rz(-2.9614095) q[0];
sx q[0];
rz(0.58191802) q[0];
x q[1];
rz(-2.3771268) q[2];
sx q[2];
rz(-1.0180961) q[2];
sx q[2];
rz(-1.1134195) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5907363) q[1];
sx q[1];
rz(-1.4155861) q[1];
sx q[1];
rz(-0.89977818) q[1];
rz(-1.8160062) q[3];
sx q[3];
rz(-1.2767226) q[3];
sx q[3];
rz(-0.40286703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7511071) q[2];
sx q[2];
rz(-2.0483053) q[2];
sx q[2];
rz(-0.92457479) q[2];
rz(3.0604002) q[3];
sx q[3];
rz(-1.3982747) q[3];
sx q[3];
rz(0.72807062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.28295949) q[0];
sx q[0];
rz(-1.2456243) q[0];
sx q[0];
rz(-1.0754732) q[0];
rz(-1.3772427) q[1];
sx q[1];
rz(-1.3721162) q[1];
sx q[1];
rz(0.65518728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0725043) q[0];
sx q[0];
rz(-3.0720815) q[0];
sx q[0];
rz(0.68060912) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8639815) q[2];
sx q[2];
rz(-0.80481968) q[2];
sx q[2];
rz(2.725986) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97444087) q[1];
sx q[1];
rz(-2.2658684) q[1];
sx q[1];
rz(-1.6593169) q[1];
rz(1.9094134) q[3];
sx q[3];
rz(-1.5430255) q[3];
sx q[3];
rz(0.16781346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7738771) q[2];
sx q[2];
rz(-0.39153063) q[2];
sx q[2];
rz(1.1565304) q[2];
rz(1.1441506) q[3];
sx q[3];
rz(-0.63026989) q[3];
sx q[3];
rz(1.8088079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3211806) q[0];
sx q[0];
rz(-2.2207566) q[0];
sx q[0];
rz(1.17571) q[0];
rz(1.3395478) q[1];
sx q[1];
rz(-2.1422155) q[1];
sx q[1];
rz(-2.7340926) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2605195) q[0];
sx q[0];
rz(-1.1620191) q[0];
sx q[0];
rz(-1.2606032) q[0];
rz(-pi) q[1];
rz(-0.48194285) q[2];
sx q[2];
rz(-2.6066268) q[2];
sx q[2];
rz(0.15967655) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4983771) q[1];
sx q[1];
rz(-2.1581894) q[1];
sx q[1];
rz(-0.024557928) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2506709) q[3];
sx q[3];
rz(-2.8615004) q[3];
sx q[3];
rz(2.1116032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0672062) q[2];
sx q[2];
rz(-1.9830474) q[2];
sx q[2];
rz(-1.0254394) q[2];
rz(-0.49753672) q[3];
sx q[3];
rz(-2.189744) q[3];
sx q[3];
rz(0.40646762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936754) q[0];
sx q[0];
rz(-1.1578639) q[0];
sx q[0];
rz(1.3854618) q[0];
rz(-2.0939743) q[1];
sx q[1];
rz(-1.3448998) q[1];
sx q[1];
rz(0.97283831) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5248643) q[0];
sx q[0];
rz(-1.4848733) q[0];
sx q[0];
rz(1.026154) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49373105) q[2];
sx q[2];
rz(-1.0876473) q[2];
sx q[2];
rz(1.8918623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.23215881) q[1];
sx q[1];
rz(-1.9240018) q[1];
sx q[1];
rz(1.4191737) q[1];
rz(-pi) q[2];
rz(2.2690569) q[3];
sx q[3];
rz(-1.2727238) q[3];
sx q[3];
rz(-2.3301017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73293781) q[2];
sx q[2];
rz(-0.32821566) q[2];
sx q[2];
rz(-2.0396063) q[2];
rz(-1.7105626) q[3];
sx q[3];
rz(-2.3165063) q[3];
sx q[3];
rz(1.2380606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86193209) q[0];
sx q[0];
rz(-1.4623549) q[0];
sx q[0];
rz(-2.097492) q[0];
rz(-0.73061371) q[1];
sx q[1];
rz(-0.63332557) q[1];
sx q[1];
rz(0.81061737) q[1];
rz(-1.2550477) q[2];
sx q[2];
rz(-2.1552255) q[2];
sx q[2];
rz(1.8926839) q[2];
rz(-2.788078) q[3];
sx q[3];
rz(-1.8302038) q[3];
sx q[3];
rz(-1.3448546) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
