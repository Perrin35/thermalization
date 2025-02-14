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
rz(0.30492914) q[0];
sx q[0];
rz(-0.066216901) q[0];
sx q[0];
rz(0.15596341) q[0];
rz(1.1671542) q[1];
sx q[1];
rz(-2.4894297) q[1];
sx q[1];
rz(-0.48643938) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64336813) q[0];
sx q[0];
rz(-0.8684477) q[0];
sx q[0];
rz(-0.81612103) q[0];
rz(2.2335792) q[2];
sx q[2];
rz(-2.5129299) q[2];
sx q[2];
rz(-2.6952621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.26775441) q[1];
sx q[1];
rz(-0.43711409) q[1];
sx q[1];
rz(-2.3972191) q[1];
x q[2];
rz(-0.21297314) q[3];
sx q[3];
rz(-2.3031182) q[3];
sx q[3];
rz(-2.4249083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.277694) q[2];
sx q[2];
rz(-0.70964491) q[2];
sx q[2];
rz(-2.7126183) q[2];
rz(-3.0947558) q[3];
sx q[3];
rz(-0.39188477) q[3];
sx q[3];
rz(-0.61280167) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.866339) q[0];
sx q[0];
rz(-2.1493122) q[0];
sx q[0];
rz(-2.476165) q[0];
rz(2.0003419) q[1];
sx q[1];
rz(-1.3803866) q[1];
sx q[1];
rz(0.64249396) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2739853) q[0];
sx q[0];
rz(-1.519975) q[0];
sx q[0];
rz(-0.7342059) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69509721) q[2];
sx q[2];
rz(-1.8029658) q[2];
sx q[2];
rz(-2.1921981) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2765668) q[1];
sx q[1];
rz(-1.5986018) q[1];
sx q[1];
rz(-0.69827484) q[1];
rz(1.8591381) q[3];
sx q[3];
rz(-2.5817462) q[3];
sx q[3];
rz(1.1091055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82681727) q[2];
sx q[2];
rz(-0.78709698) q[2];
sx q[2];
rz(0.55244201) q[2];
rz(-1.4779444) q[3];
sx q[3];
rz(-1.843957) q[3];
sx q[3];
rz(-2.4655931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33573547) q[0];
sx q[0];
rz(-2.4215846) q[0];
sx q[0];
rz(-0.82255256) q[0];
rz(0.81126732) q[1];
sx q[1];
rz(-2.8370116) q[1];
sx q[1];
rz(-1.4623581) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2967591) q[0];
sx q[0];
rz(-0.61723304) q[0];
sx q[0];
rz(0.69302859) q[0];
rz(-2.7285568) q[2];
sx q[2];
rz(-1.6390642) q[2];
sx q[2];
rz(-0.78536805) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2980256) q[1];
sx q[1];
rz(-2.0223534) q[1];
sx q[1];
rz(0.3275015) q[1];
rz(-2.3554166) q[3];
sx q[3];
rz(-0.42359951) q[3];
sx q[3];
rz(2.9860403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2751969) q[2];
sx q[2];
rz(-0.82922816) q[2];
sx q[2];
rz(2.5466476) q[2];
rz(-2.7368937) q[3];
sx q[3];
rz(-2.0345104) q[3];
sx q[3];
rz(-1.8428724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.4054656) q[0];
sx q[0];
rz(-1.2768224) q[0];
sx q[0];
rz(-0.24293105) q[0];
rz(-0.88492197) q[1];
sx q[1];
rz(-0.72598571) q[1];
sx q[1];
rz(-0.0028217908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76317518) q[0];
sx q[0];
rz(-0.88930128) q[0];
sx q[0];
rz(-2.6386518) q[0];
x q[1];
rz(0.52619885) q[2];
sx q[2];
rz(-2.1251273) q[2];
sx q[2];
rz(0.89579158) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62880781) q[1];
sx q[1];
rz(-1.6411726) q[1];
sx q[1];
rz(2.2712399) q[1];
x q[2];
rz(-0.97798621) q[3];
sx q[3];
rz(-2.1007256) q[3];
sx q[3];
rz(-1.1150313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.60488492) q[2];
sx q[2];
rz(-1.9900091) q[2];
sx q[2];
rz(0.73337698) q[2];
rz(0.65507656) q[3];
sx q[3];
rz(-2.8337182) q[3];
sx q[3];
rz(-2.5797599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4578399) q[0];
sx q[0];
rz(-0.34128749) q[0];
sx q[0];
rz(0.14232464) q[0];
rz(-1.3617474) q[1];
sx q[1];
rz(-2.7889377) q[1];
sx q[1];
rz(-0.49837643) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179006) q[0];
sx q[0];
rz(-1.5091584) q[0];
sx q[0];
rz(-0.85718244) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3827078) q[2];
sx q[2];
rz(-1.6079796) q[2];
sx q[2];
rz(2.5868724) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8784762) q[1];
sx q[1];
rz(-0.85678393) q[1];
sx q[1];
rz(-1.962838) q[1];
x q[2];
rz(0.33739319) q[3];
sx q[3];
rz(-0.80885799) q[3];
sx q[3];
rz(-2.0615426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.132823) q[2];
sx q[2];
rz(-1.885773) q[2];
sx q[2];
rz(2.447017) q[2];
rz(-2.3092367) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(2.8913403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4271127) q[0];
sx q[0];
rz(-2.6321593) q[0];
sx q[0];
rz(-2.4795649) q[0];
rz(-1.9561249) q[1];
sx q[1];
rz(-2.1123501) q[1];
sx q[1];
rz(1.1659291) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2070475) q[0];
sx q[0];
rz(-0.70431346) q[0];
sx q[0];
rz(-2.8273316) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2189617) q[2];
sx q[2];
rz(-2.5220519) q[2];
sx q[2];
rz(2.9094686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37167376) q[1];
sx q[1];
rz(-3.0536302) q[1];
sx q[1];
rz(-0.96925737) q[1];
x q[2];
rz(-2.1816129) q[3];
sx q[3];
rz(-2.0997542) q[3];
sx q[3];
rz(2.4774266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.379443) q[2];
sx q[2];
rz(-1.2078441) q[2];
sx q[2];
rz(0.70972788) q[2];
rz(2.659667) q[3];
sx q[3];
rz(-0.47526264) q[3];
sx q[3];
rz(-3.1221534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703281) q[0];
sx q[0];
rz(-2.1605991) q[0];
sx q[0];
rz(1.574466) q[0];
rz(-1.6471242) q[1];
sx q[1];
rz(-2.6946805) q[1];
sx q[1];
rz(2.2692197) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14022889) q[0];
sx q[0];
rz(-1.8793545) q[0];
sx q[0];
rz(-2.9635327) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21193223) q[2];
sx q[2];
rz(-2.4080347) q[2];
sx q[2];
rz(2.9447945) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0371767) q[1];
sx q[1];
rz(-0.1104011) q[1];
sx q[1];
rz(1.8502185) q[1];
rz(-2.2254785) q[3];
sx q[3];
rz(-1.5512084) q[3];
sx q[3];
rz(-0.71924984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.81277728) q[2];
sx q[2];
rz(-2.2810563) q[2];
sx q[2];
rz(-2.013773) q[2];
rz(-0.40337107) q[3];
sx q[3];
rz(-1.4453459) q[3];
sx q[3];
rz(-0.71322125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047852) q[0];
sx q[0];
rz(-1.9984364) q[0];
sx q[0];
rz(-1.8444201) q[0];
rz(0.5212658) q[1];
sx q[1];
rz(-1.2626941) q[1];
sx q[1];
rz(-3.0788132) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2817665) q[0];
sx q[0];
rz(-1.504557) q[0];
sx q[0];
rz(2.0145922) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8028444) q[2];
sx q[2];
rz(-1.9221483) q[2];
sx q[2];
rz(2.6356489) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6160342) q[1];
sx q[1];
rz(-1.8174606) q[1];
sx q[1];
rz(2.9933369) q[1];
x q[2];
rz(2.2388458) q[3];
sx q[3];
rz(-0.76551907) q[3];
sx q[3];
rz(-0.79878858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2433743) q[2];
sx q[2];
rz(-0.90460193) q[2];
sx q[2];
rz(-0.96837366) q[2];
rz(2.6280256) q[3];
sx q[3];
rz(-1.3526724) q[3];
sx q[3];
rz(-2.8106522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53720713) q[0];
sx q[0];
rz(-2.4810915) q[0];
sx q[0];
rz(0.78395098) q[0];
rz(-1.9369269) q[1];
sx q[1];
rz(-2.7684863) q[1];
sx q[1];
rz(-1.7663667) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2581552) q[0];
sx q[0];
rz(-0.34314197) q[0];
sx q[0];
rz(-0.68455066) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.398574) q[2];
sx q[2];
rz(-2.9841514) q[2];
sx q[2];
rz(0.96913183) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8177796) q[1];
sx q[1];
rz(-2.166417) q[1];
sx q[1];
rz(0.88902529) q[1];
x q[2];
rz(-2.9638644) q[3];
sx q[3];
rz(-1.2627708) q[3];
sx q[3];
rz(2.7975688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84460008) q[2];
sx q[2];
rz(-0.30868369) q[2];
sx q[2];
rz(2.6958579) q[2];
rz(-2.5388057) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(-0.84356892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26850253) q[0];
sx q[0];
rz(-0.52274811) q[0];
sx q[0];
rz(-0.52421808) q[0];
rz(1.9408608) q[1];
sx q[1];
rz(-1.2501161) q[1];
sx q[1];
rz(-1.9573617) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6551334) q[0];
sx q[0];
rz(-2.7441027) q[0];
sx q[0];
rz(-2.6349806) q[0];
rz(-pi) q[1];
x q[1];
rz(0.068151926) q[2];
sx q[2];
rz(-1.8859204) q[2];
sx q[2];
rz(-1.5308876) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47329119) q[1];
sx q[1];
rz(-1.4379825) q[1];
sx q[1];
rz(1.5584438) q[1];
rz(0.04751398) q[3];
sx q[3];
rz(-2.5121452) q[3];
sx q[3];
rz(0.87241064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21620096) q[2];
sx q[2];
rz(-2.6915458) q[2];
sx q[2];
rz(-2.0897384) q[2];
rz(2.9205186) q[3];
sx q[3];
rz(-1.8408006) q[3];
sx q[3];
rz(0.93824798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5476407) q[0];
sx q[0];
rz(-1.6407536) q[0];
sx q[0];
rz(-3.092691) q[0];
rz(2.7652057) q[1];
sx q[1];
rz(-2.2793437) q[1];
sx q[1];
rz(1.9752165) q[1];
rz(0.57009956) q[2];
sx q[2];
rz(-1.5440294) q[2];
sx q[2];
rz(3.0837035) q[2];
rz(0.71394271) q[3];
sx q[3];
rz(-1.1957914) q[3];
sx q[3];
rz(1.316432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
