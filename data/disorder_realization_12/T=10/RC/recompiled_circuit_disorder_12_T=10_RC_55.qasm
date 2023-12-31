OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(3.8656524) q[0];
sx q[0];
rz(11.081628) q[0];
rz(4.3447189) q[1];
sx q[1];
rz(0.523518) q[1];
sx q[1];
rz(8.5365774) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9125497) q[0];
sx q[0];
rz(-2.9648844) q[0];
sx q[0];
rz(2.8852709) q[0];
rz(1.8445831) q[2];
sx q[2];
rz(-0.55371504) q[2];
sx q[2];
rz(-0.48534976) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8623212) q[1];
sx q[1];
rz(-2.9746378) q[1];
sx q[1];
rz(0.0032940666) q[1];
rz(-0.69859759) q[3];
sx q[3];
rz(-1.5961831) q[3];
sx q[3];
rz(-2.4217525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.28329864) q[2];
sx q[2];
rz(-2.7259939) q[2];
sx q[2];
rz(-2.0236012) q[2];
rz(2.9962712) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(3.0956691) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5685101) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(2.4867687) q[0];
rz(1.2163935) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(-3.0156946) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1032216) q[0];
sx q[0];
rz(-1.201259) q[0];
sx q[0];
rz(-2.080337) q[0];
rz(-pi) q[1];
rz(0.89181487) q[2];
sx q[2];
rz(-2.0795155) q[2];
sx q[2];
rz(1.5985135) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7014335) q[1];
sx q[1];
rz(-1.6007746) q[1];
sx q[1];
rz(-1.4440126) q[1];
x q[2];
rz(0.84729362) q[3];
sx q[3];
rz(-1.5788659) q[3];
sx q[3];
rz(-1.1524259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(-2.753567) q[2];
rz(-1.7175425) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(0.58732906) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5557142) q[0];
sx q[0];
rz(-0.782574) q[0];
sx q[0];
rz(0.079332381) q[0];
rz(-3.0575867) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(1.9817339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.738049) q[0];
sx q[0];
rz(-1.6532716) q[0];
sx q[0];
rz(1.0595881) q[0];
rz(-0.50767501) q[2];
sx q[2];
rz(-1.0055563) q[2];
sx q[2];
rz(2.450168) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.11324727) q[1];
sx q[1];
rz(-2.5302437) q[1];
sx q[1];
rz(-2.4436185) q[1];
rz(-pi) q[2];
rz(-2.1163164) q[3];
sx q[3];
rz(-1.5226411) q[3];
sx q[3];
rz(0.84850509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7362061) q[2];
sx q[2];
rz(-1.8182886) q[2];
sx q[2];
rz(-3.0333056) q[2];
rz(2.4978499) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9765587) q[0];
sx q[0];
rz(-1.7465916) q[0];
sx q[0];
rz(1.6595586) q[0];
rz(-2.4404793) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(-2.8569417) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4581504) q[0];
sx q[0];
rz(-2.0476641) q[0];
sx q[0];
rz(-0.083806888) q[0];
rz(-pi) q[1];
rz(0.95732032) q[2];
sx q[2];
rz(-0.34733221) q[2];
sx q[2];
rz(2.6383102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9173911) q[1];
sx q[1];
rz(-0.14897878) q[1];
sx q[1];
rz(-1.0148744) q[1];
x q[2];
rz(-2.2399726) q[3];
sx q[3];
rz(-1.1812783) q[3];
sx q[3];
rz(-0.31182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9099137) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(-0.56751928) q[2];
rz(-2.7275758) q[3];
sx q[3];
rz(-2.0040138) q[3];
sx q[3];
rz(2.0295985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.652997) q[0];
sx q[0];
rz(-0.96019205) q[0];
sx q[0];
rz(-1.3265142) q[0];
rz(1.9891706) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(-0.93793905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4612761) q[0];
sx q[0];
rz(-0.094705908) q[0];
sx q[0];
rz(1.3374431) q[0];
rz(-pi) q[1];
rz(-2.3100501) q[2];
sx q[2];
rz(-0.70280308) q[2];
sx q[2];
rz(-2.0895095) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8590947) q[1];
sx q[1];
rz(-1.5364093) q[1];
sx q[1];
rz(0.026245898) q[1];
x q[2];
rz(-2.786817) q[3];
sx q[3];
rz(-1.6924904) q[3];
sx q[3];
rz(1.0697168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9662629) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(2.1002634) q[2];
rz(1.3025618) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(1.3180102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43907169) q[0];
sx q[0];
rz(-1.9938001) q[0];
sx q[0];
rz(0.55737108) q[0];
rz(2.5769261) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(-0.55647892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0628478) q[0];
sx q[0];
rz(-1.9468465) q[0];
sx q[0];
rz(-1.3486805) q[0];
rz(-pi) q[1];
rz(1.7700023) q[2];
sx q[2];
rz(-2.4161985) q[2];
sx q[2];
rz(-0.2142011) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4887052) q[1];
sx q[1];
rz(-0.31213752) q[1];
sx q[1];
rz(0.94241347) q[1];
rz(-pi) q[2];
rz(0.77575404) q[3];
sx q[3];
rz(-1.6072825) q[3];
sx q[3];
rz(-0.072582399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53211987) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(-1.2247941) q[2];
rz(-1.6312284) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(-1.640655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0625793) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(0.042908948) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(-2.6409805) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0450889) q[0];
sx q[0];
rz(-1.5363662) q[0];
sx q[0];
rz(-1.3929429) q[0];
rz(-pi) q[1];
rz(-2.1130354) q[2];
sx q[2];
rz(-1.2985897) q[2];
sx q[2];
rz(0.82424639) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88922933) q[1];
sx q[1];
rz(-1.4503308) q[1];
sx q[1];
rz(-0.40263386) q[1];
rz(-1.9686437) q[3];
sx q[3];
rz(-2.2865191) q[3];
sx q[3];
rz(1.0637103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5381955) q[2];
sx q[2];
rz(-2.0624702) q[2];
sx q[2];
rz(0.67561692) q[2];
rz(-0.44089857) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(2.6202257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(1.0674397) q[0];
rz(0.43287977) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(-1.4656461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059793652) q[0];
sx q[0];
rz(-1.1104715) q[0];
sx q[0];
rz(2.9914809) q[0];
rz(-pi) q[1];
rz(0.18657121) q[2];
sx q[2];
rz(-1.6779043) q[2];
sx q[2];
rz(2.4466116) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2035985) q[1];
sx q[1];
rz(-1.9059056) q[1];
sx q[1];
rz(-2.0957698) q[1];
rz(-pi) q[2];
rz(-2.500324) q[3];
sx q[3];
rz(-2.0740168) q[3];
sx q[3];
rz(-0.36676952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8026768) q[2];
sx q[2];
rz(-2.2854476) q[2];
sx q[2];
rz(-1.6652997) q[2];
rz(0.87351292) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840238) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(2.4043758) q[0];
rz(-0.018741477) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(0.92528701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6773274) q[0];
sx q[0];
rz(-2.263501) q[0];
sx q[0];
rz(0.6662743) q[0];
rz(-1.7032353) q[2];
sx q[2];
rz(-0.91184154) q[2];
sx q[2];
rz(3.1140285) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3146989) q[1];
sx q[1];
rz(-1.3876545) q[1];
sx q[1];
rz(-1.7032743) q[1];
rz(-2.0426644) q[3];
sx q[3];
rz(-2.3657551) q[3];
sx q[3];
rz(2.8692506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3725738) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(-2.1378689) q[2];
rz(-0.090099661) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1019679) q[0];
sx q[0];
rz(-1.573338) q[0];
sx q[0];
rz(1.2596624) q[0];
rz(0.26578495) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(-2.3840747) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19125464) q[0];
sx q[0];
rz(-1.7636824) q[0];
sx q[0];
rz(3.0015776) q[0];
x q[1];
rz(-2.9847758) q[2];
sx q[2];
rz(-2.1553401) q[2];
sx q[2];
rz(0.90442327) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0150891) q[1];
sx q[1];
rz(-1.9941829) q[1];
sx q[1];
rz(-2.8979315) q[1];
rz(-0.032978756) q[3];
sx q[3];
rz(-0.71027256) q[3];
sx q[3];
rz(-2.9041293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1511128) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(-0.63888597) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6486075) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(1.5007301) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(2.7040504) q[2];
sx q[2];
rz(-1.8716639) q[2];
sx q[2];
rz(1.5834091) q[2];
rz(0.99132514) q[3];
sx q[3];
rz(-1.2990868) q[3];
sx q[3];
rz(2.2494863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
