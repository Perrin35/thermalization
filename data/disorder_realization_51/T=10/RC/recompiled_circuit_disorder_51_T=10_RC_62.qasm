OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(0.18145951) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(-2.0884617) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3436514) q[0];
sx q[0];
rz(-2.6290253) q[0];
sx q[0];
rz(-2.0128065) q[0];
rz(1.8671145) q[2];
sx q[2];
rz(-2.1726492) q[2];
sx q[2];
rz(-1.559343) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44838201) q[1];
sx q[1];
rz(-0.795185) q[1];
sx q[1];
rz(2.8661212) q[1];
x q[2];
rz(-1.1374723) q[3];
sx q[3];
rz(-0.37131272) q[3];
sx q[3];
rz(0.565688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.16333214) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(2.1286428) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(-3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0435836) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(-0.59511551) q[0];
rz(2.0960506) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(1.5140623) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1986952) q[0];
sx q[0];
rz(-1.1637582) q[0];
sx q[0];
rz(2.1172949) q[0];
rz(-pi) q[1];
x q[1];
rz(0.02818429) q[2];
sx q[2];
rz(-2.3079254) q[2];
sx q[2];
rz(0.29836269) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4501614) q[1];
sx q[1];
rz(-1.3988004) q[1];
sx q[1];
rz(-1.3732234) q[1];
rz(-2.0239003) q[3];
sx q[3];
rz(-0.11370224) q[3];
sx q[3];
rz(2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4864768) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(0.79616037) q[2];
rz(0.97186175) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(-2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7650836) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(0.31072101) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(1.6832738) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7539833) q[0];
sx q[0];
rz(-1.4235272) q[0];
sx q[0];
rz(3.0940042) q[0];
rz(-pi) q[1];
rz(-2.4086191) q[2];
sx q[2];
rz(-2.3419215) q[2];
sx q[2];
rz(-3.1250931) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2787784) q[1];
sx q[1];
rz(-1.9796951) q[1];
sx q[1];
rz(1.5120974) q[1];
rz(-1.5536669) q[3];
sx q[3];
rz(-1.8855842) q[3];
sx q[3];
rz(2.0103612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1594499) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(0.67908755) q[2];
rz(-1.3726161) q[3];
sx q[3];
rz(-1.8434098) q[3];
sx q[3];
rz(-2.1931271) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452334) q[0];
sx q[0];
rz(-2.5722752) q[0];
sx q[0];
rz(-2.0171719) q[0];
rz(-1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(1.6569998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67498818) q[0];
sx q[0];
rz(-0.059048422) q[0];
sx q[0];
rz(-2.1112291) q[0];
x q[1];
rz(-1.053327) q[2];
sx q[2];
rz(-1.4033068) q[2];
sx q[2];
rz(-2.6094764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.209219) q[1];
sx q[1];
rz(-1.1353496) q[1];
sx q[1];
rz(-2.2030764) q[1];
rz(-pi) q[2];
rz(-0.72317601) q[3];
sx q[3];
rz(-1.3190862) q[3];
sx q[3];
rz(0.46079208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0174039) q[2];
sx q[2];
rz(-1.2839395) q[2];
sx q[2];
rz(-2.1253288) q[2];
rz(1.4034363) q[3];
sx q[3];
rz(-1.4973463) q[3];
sx q[3];
rz(-2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6338585) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(0.39988363) q[0];
rz(-1.1625066) q[1];
sx q[1];
rz(-1.3299273) q[1];
sx q[1];
rz(0.17366017) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7090209) q[0];
sx q[0];
rz(-3.0637494) q[0];
sx q[0];
rz(2.7709333) q[0];
x q[1];
rz(0.47541754) q[2];
sx q[2];
rz(-0.43736514) q[2];
sx q[2];
rz(1.4471444) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34827161) q[1];
sx q[1];
rz(-1.6760567) q[1];
sx q[1];
rz(-0.028491032) q[1];
rz(0.33186121) q[3];
sx q[3];
rz(-1.7680941) q[3];
sx q[3];
rz(-2.8062537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6614723) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(-0.76812569) q[2];
rz(0.85401946) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1059234) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(1.4703898) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(-1.1434198) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2238732) q[0];
sx q[0];
rz(-1.4871162) q[0];
sx q[0];
rz(-2.9907945) q[0];
x q[1];
rz(-1.9009695) q[2];
sx q[2];
rz(-0.99324838) q[2];
sx q[2];
rz(-1.6799048) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1868362) q[1];
sx q[1];
rz(-1.3969159) q[1];
sx q[1];
rz(-1.1980921) q[1];
x q[2];
rz(-2.8361736) q[3];
sx q[3];
rz(-2.0913887) q[3];
sx q[3];
rz(2.6484495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4986971) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(-1.6112304) q[2];
rz(-1.4536084) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(-2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11944184) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(-1.4402333) q[0];
rz(2.4123736) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(-2.008332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5232552) q[0];
sx q[0];
rz(-1.6132857) q[0];
sx q[0];
rz(-2.8190814) q[0];
rz(-1.5845756) q[2];
sx q[2];
rz(-0.85859495) q[2];
sx q[2];
rz(-0.53879246) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1843695) q[1];
sx q[1];
rz(-1.8436699) q[1];
sx q[1];
rz(-0.89783122) q[1];
rz(-pi) q[2];
rz(0.26182884) q[3];
sx q[3];
rz(-1.4187519) q[3];
sx q[3];
rz(2.9933628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7363654) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(2.6531632) q[2];
rz(1.3119665) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(-0.64129889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(-0.07117614) q[0];
rz(0.03216234) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(1.9326928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6541518) q[0];
sx q[0];
rz(-1.7983266) q[0];
sx q[0];
rz(-0.80765101) q[0];
rz(-1.8078631) q[2];
sx q[2];
rz(-1.7500688) q[2];
sx q[2];
rz(2.8693503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68122411) q[1];
sx q[1];
rz(-1.7011233) q[1];
sx q[1];
rz(0.016101109) q[1];
rz(-pi) q[2];
rz(-2.5197221) q[3];
sx q[3];
rz(-1.2688046) q[3];
sx q[3];
rz(-0.57884502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9341087) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(1.5709546) q[2];
rz(-2.2682244) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1677925) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(2.8299676) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(-1.5100381) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.894282) q[0];
sx q[0];
rz(-1.7705435) q[0];
sx q[0];
rz(0.32379177) q[0];
x q[1];
rz(-1.0636343) q[2];
sx q[2];
rz(-2.7535451) q[2];
sx q[2];
rz(-1.0840814) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77126399) q[1];
sx q[1];
rz(-2.1762098) q[1];
sx q[1];
rz(1.3243115) q[1];
rz(-pi) q[2];
rz(-2.908913) q[3];
sx q[3];
rz(-1.7494697) q[3];
sx q[3];
rz(1.0089547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2733549) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(-2.3256425) q[2];
rz(-0.50968918) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2508535) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(2.5157805) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(2.4831916) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31357665) q[0];
sx q[0];
rz(-0.99452924) q[0];
sx q[0];
rz(0.97998951) q[0];
rz(-0.41215956) q[2];
sx q[2];
rz(-0.71752749) q[2];
sx q[2];
rz(1.621643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.876993) q[1];
sx q[1];
rz(-1.0942642) q[1];
sx q[1];
rz(-0.18696733) q[1];
x q[2];
rz(-2.5913521) q[3];
sx q[3];
rz(-0.60866683) q[3];
sx q[3];
rz(-0.88702162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6752424) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(2.8038483) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52453775) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(-1.2383923) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(2.1970489) q[2];
sx q[2];
rz(-1.7164451) q[2];
sx q[2];
rz(-3.0838983) q[2];
rz(0.57237207) q[3];
sx q[3];
rz(-1.5170245) q[3];
sx q[3];
rz(-1.7046884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
