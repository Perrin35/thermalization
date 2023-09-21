OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.934259) q[0];
sx q[0];
rz(-0.59036314) q[0];
sx q[0];
rz(-2.7705749) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(2.5420904) q[1];
sx q[1];
rz(11.190344) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3289514) q[0];
sx q[0];
rz(-1.1321804) q[0];
sx q[0];
rz(-0.89226512) q[0];
x q[1];
rz(-2.4260169) q[2];
sx q[2];
rz(-1.1872429) q[2];
sx q[2];
rz(-0.54563145) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8810597) q[1];
sx q[1];
rz(-0.36306371) q[1];
sx q[1];
rz(-1.1313603) q[1];
x q[2];
rz(-1.6269496) q[3];
sx q[3];
rz(-2.3893642) q[3];
sx q[3];
rz(0.55263954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0573037) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(0.98891813) q[2];
rz(0.75254285) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(-2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20724021) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(-1.9616615) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(-0.72431272) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2253101) q[0];
sx q[0];
rz(-0.75964576) q[0];
sx q[0];
rz(2.0347974) q[0];
rz(-2.7472277) q[2];
sx q[2];
rz(-0.21557237) q[2];
sx q[2];
rz(2.8087316) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69246768) q[1];
sx q[1];
rz(-1.7763224) q[1];
sx q[1];
rz(1.9301901) q[1];
x q[2];
rz(-0.23726666) q[3];
sx q[3];
rz(-1.3788584) q[3];
sx q[3];
rz(2.3290079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-2.779707) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996465) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(1.2190855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8273979) q[0];
sx q[0];
rz(-2.6959246) q[0];
sx q[0];
rz(2.218194) q[0];
x q[1];
rz(-0.8692603) q[2];
sx q[2];
rz(-1.3420891) q[2];
sx q[2];
rz(1.2094091) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3238941) q[1];
sx q[1];
rz(-2.5496799) q[1];
sx q[1];
rz(-1.1281668) q[1];
rz(-pi) q[2];
rz(-1.838802) q[3];
sx q[3];
rz(-1.5558814) q[3];
sx q[3];
rz(-2.5915495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7200155) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.4455618) q[2];
rz(-0.56882632) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9179984) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(0.7154243) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(-0.82675654) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2402156) q[0];
sx q[0];
rz(-2.4670521) q[0];
sx q[0];
rz(2.4504689) q[0];
rz(1.4364169) q[2];
sx q[2];
rz(-0.80766404) q[2];
sx q[2];
rz(-0.32854167) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0930867) q[1];
sx q[1];
rz(-2.6566681) q[1];
sx q[1];
rz(0.64694689) q[1];
x q[2];
rz(-2.4171962) q[3];
sx q[3];
rz(-1.8076234) q[3];
sx q[3];
rz(2.608992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15239079) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(-1.3789122) q[2];
rz(3.0692696) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(1.6453843) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82692659) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(2.2391438) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(2.856423) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641159) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(1.478273) q[0];
rz(-pi) q[1];
rz(0.6638078) q[2];
sx q[2];
rz(-1.2592578) q[2];
sx q[2];
rz(0.89154348) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3030745) q[1];
sx q[1];
rz(-2.8285366) q[1];
sx q[1];
rz(-1.5694373) q[1];
rz(0.22484803) q[3];
sx q[3];
rz(-0.41089155) q[3];
sx q[3];
rz(1.9083244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29331648) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(2.6021393) q[2];
rz(-2.8347677) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(-2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25813112) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(0.69333386) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(-1.3670115) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088928662) q[0];
sx q[0];
rz(-3.0658709) q[0];
sx q[0];
rz(-1.1195539) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96732803) q[2];
sx q[2];
rz(-0.8562932) q[2];
sx q[2];
rz(-2.8514903) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4945592) q[1];
sx q[1];
rz(-1.651598) q[1];
sx q[1];
rz(-1.2809491) q[1];
rz(-pi) q[2];
rz(2.7886224) q[3];
sx q[3];
rz(-1.4735231) q[3];
sx q[3];
rz(-0.73247611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29629016) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(0.40346754) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24213174) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(0.8738628) q[0];
rz(2.6938687) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(-1.9708995) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6252977) q[0];
sx q[0];
rz(-1.5769616) q[0];
sx q[0];
rz(0.023029285) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11622073) q[2];
sx q[2];
rz(-1.5201609) q[2];
sx q[2];
rz(-1.9192139) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16329855) q[1];
sx q[1];
rz(-0.98493176) q[1];
sx q[1];
rz(-2.3818124) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4865173) q[3];
sx q[3];
rz(-2.5518637) q[3];
sx q[3];
rz(0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0968904) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(-0.61075413) q[2];
rz(-0.47510535) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(-2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(2.8444667) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(-0.64613211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32983366) q[0];
sx q[0];
rz(-2.9114897) q[0];
sx q[0];
rz(-2.0287201) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7392776) q[2];
sx q[2];
rz(-2.0373166) q[2];
sx q[2];
rz(2.434935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99107332) q[1];
sx q[1];
rz(-2.0690314) q[1];
sx q[1];
rz(1.1789765) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93692245) q[3];
sx q[3];
rz(-0.86808944) q[3];
sx q[3];
rz(-0.015451775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(2.440195) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(2.3440857) q[0];
rz(2.6240255) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(-0.10841766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3603044) q[0];
sx q[0];
rz(-2.4753248) q[0];
sx q[0];
rz(-1.4772619) q[0];
rz(-pi) q[1];
rz(0.046594521) q[2];
sx q[2];
rz(-1.7492883) q[2];
sx q[2];
rz(-2.5537234) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79418102) q[1];
sx q[1];
rz(-1.3550183) q[1];
sx q[1];
rz(-0.17098917) q[1];
rz(-pi) q[2];
rz(2.2750862) q[3];
sx q[3];
rz(-1.4276854) q[3];
sx q[3];
rz(-1.0673616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-2.8016395) q[2];
rz(0.41839504) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6091992) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(-0.6788196) q[0];
rz(2.7774096) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(-3.0864339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.178135) q[0];
sx q[0];
rz(-1.194251) q[0];
sx q[0];
rz(1.6403273) q[0];
x q[1];
rz(2.8137384) q[2];
sx q[2];
rz(-2.556986) q[2];
sx q[2];
rz(-0.84529982) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.060505796) q[1];
sx q[1];
rz(-2.2844237) q[1];
sx q[1];
rz(2.39141) q[1];
rz(0.03667128) q[3];
sx q[3];
rz(-2.026537) q[3];
sx q[3];
rz(2.7367221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(2.6514163) q[2];
rz(3.0040719) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162083) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(-0.31221496) q[2];
sx q[2];
rz(-1.6630465) q[2];
sx q[2];
rz(1.9065471) q[2];
rz(-2.901554) q[3];
sx q[3];
rz(-1.5272899) q[3];
sx q[3];
rz(1.8361113) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
