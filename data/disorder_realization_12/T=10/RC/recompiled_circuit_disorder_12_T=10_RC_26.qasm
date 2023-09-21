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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6523525) q[0];
sx q[0];
rz(-1.3999192) q[0];
sx q[0];
rz(-1.6160374) q[0];
rz(-pi) q[1];
rz(-2.1076803) q[2];
sx q[2];
rz(-1.7134588) q[2];
sx q[2];
rz(-2.2906274) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8656617) q[1];
sx q[1];
rz(-1.7377503) q[1];
sx q[1];
rz(1.5702412) q[1];
rz(-pi) q[2];
rz(-2.4429951) q[3];
sx q[3];
rz(-1.5961831) q[3];
sx q[3];
rz(-0.7198402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.858294) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(1.1179914) q[2];
rz(-2.9962712) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(0.045923559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685101) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(2.4867687) q[0];
rz(-1.2163935) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(3.0156946) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1830131) q[0];
sx q[0];
rz(-2.5218681) q[0];
sx q[0];
rz(0.89967863) q[0];
x q[1];
rz(2.2497778) q[2];
sx q[2];
rz(-2.0795155) q[2];
sx q[2];
rz(-1.5985135) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1344578) q[1];
sx q[1];
rz(-1.4440698) q[1];
sx q[1];
rz(3.111372) q[1];
rz(1.5829854) q[3];
sx q[3];
rz(-0.72353957) q[3];
sx q[3];
rz(0.42750588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(0.38802567) q[2];
rz(-1.7175425) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(-0.58732906) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5557142) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(3.0622603) q[0];
rz(0.084005984) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(-1.9817339) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2134316) q[0];
sx q[0];
rz(-2.0800989) q[0];
sx q[0];
rz(0.094497735) q[0];
x q[1];
rz(-0.94295393) q[2];
sx q[2];
rz(-1.1477594) q[2];
sx q[2];
rz(-1.9726276) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11324727) q[1];
sx q[1];
rz(-0.61134895) q[1];
sx q[1];
rz(-0.69797413) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0852774) q[3];
sx q[3];
rz(-2.115613) q[3];
sx q[3];
rz(-0.69308263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40538654) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(-3.0333056) q[2];
rz(-0.64374271) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.165034) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(1.4820341) q[0];
rz(2.4404793) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(-2.8569417) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2157001) q[0];
sx q[0];
rz(-1.6452351) q[0];
sx q[0];
rz(1.0924933) q[0];
rz(0.95732032) q[2];
sx q[2];
rz(-2.7942604) q[2];
sx q[2];
rz(0.50328244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33672562) q[1];
sx q[1];
rz(-1.4443828) q[1];
sx q[1];
rz(0.079041914) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48216362) q[3];
sx q[3];
rz(-0.95957047) q[3];
sx q[3];
rz(1.5507444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.231679) q[2];
sx q[2];
rz(-2.173013) q[2];
sx q[2];
rz(-0.56751928) q[2];
rz(-2.7275758) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(-2.0295985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.652997) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(1.3265142) q[0];
rz(1.9891706) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(2.2036536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4834578) q[0];
sx q[0];
rz(-1.5489274) q[0];
sx q[0];
rz(-1.4786426) q[0];
x q[1];
rz(2.1300975) q[2];
sx q[2];
rz(-2.0213631) q[2];
sx q[2];
rz(0.089103854) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.93463072) q[1];
sx q[1];
rz(-3.0983371) q[1];
sx q[1];
rz(-2.2224777) q[1];
rz(-pi) q[2];
rz(-0.35477562) q[3];
sx q[3];
rz(-1.4491023) q[3];
sx q[3];
rz(-2.0718758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9662629) q[2];
sx q[2];
rz(-2.9523409) q[2];
sx q[2];
rz(2.1002634) q[2];
rz(1.3025618) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.8235824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(2.5842216) q[0];
rz(2.5769261) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(0.55647892) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508931) q[0];
sx q[0];
rz(-1.7771582) q[0];
sx q[0];
rz(2.7569689) q[0];
rz(-pi) q[1];
rz(-1.3715903) q[2];
sx q[2];
rz(-0.7253941) q[2];
sx q[2];
rz(-2.9273916) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.52292697) q[1];
sx q[1];
rz(-1.3892801) q[1];
sx q[1];
rz(-1.8261441) q[1];
x q[2];
rz(0.052080215) q[3];
sx q[3];
rz(-2.3651603) q[3];
sx q[3];
rz(-1.6805502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53211987) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(1.2247941) q[2];
rz(-1.6312284) q[3];
sx q[3];
rz(-1.3204201) q[3];
sx q[3];
rz(1.640655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0790134) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(-3.0986837) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(-2.6409805) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096503784) q[0];
sx q[0];
rz(-1.5363662) q[0];
sx q[0];
rz(1.7486497) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0285573) q[2];
sx q[2];
rz(-1.2985897) q[2];
sx q[2];
rz(2.3173463) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4088879) q[1];
sx q[1];
rz(-1.171247) q[1];
sx q[1];
rz(-1.4399745) q[1];
rz(2.72243) q[3];
sx q[3];
rz(-2.3401642) q[3];
sx q[3];
rz(-0.49405801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5381955) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(2.4659757) q[2];
rz(0.44089857) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.065141) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(2.0741529) q[0];
rz(-2.7087129) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(-1.6759466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6976801) q[0];
sx q[0];
rz(-1.4364103) q[0];
sx q[0];
rz(-1.105955) q[0];
x q[1];
rz(-2.9550214) q[2];
sx q[2];
rz(-1.4636883) q[2];
sx q[2];
rz(0.69498108) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9379942) q[1];
sx q[1];
rz(-1.9059056) q[1];
sx q[1];
rz(2.0957698) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96887178) q[3];
sx q[3];
rz(-1.0191917) q[3];
sx q[3];
rz(-1.5495891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33891588) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(1.476293) q[2];
rz(2.2680797) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(-2.6749558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2575689) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(2.4043758) q[0];
rz(-0.018741477) q[1];
sx q[1];
rz(-0.32967162) q[1];
sx q[1];
rz(-0.92528701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5718482) q[0];
sx q[0];
rz(-1.0751372) q[0];
sx q[0];
rz(0.75832383) q[0];
x q[1];
rz(1.7032353) q[2];
sx q[2];
rz(-0.91184154) q[2];
sx q[2];
rz(0.027564136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.195203) q[1];
sx q[1];
rz(-0.22559799) q[1];
sx q[1];
rz(2.5220847) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0426644) q[3];
sx q[3];
rz(-0.77583757) q[3];
sx q[3];
rz(2.8692506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76901889) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(2.1378689) q[2];
rz(3.051493) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(-2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039624778) q[0];
sx q[0];
rz(-1.573338) q[0];
sx q[0];
rz(-1.2596624) q[0];
rz(-0.26578495) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(2.3840747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.950338) q[0];
sx q[0];
rz(-1.7636824) q[0];
sx q[0];
rz(3.0015776) q[0];
rz(-pi) q[1];
rz(2.1610356) q[2];
sx q[2];
rz(-1.7014116) q[2];
sx q[2];
rz(-2.5622501) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4841008) q[1];
sx q[1];
rz(-1.3490281) q[1];
sx q[1];
rz(1.1360672) q[1];
x q[2];
rz(1.5424472) q[3];
sx q[3];
rz(-2.2806014) q[3];
sx q[3];
rz(-2.860644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.99047986) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(0.79375664) q[2];
rz(2.5027067) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(2.1209774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(1.5007301) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(2.5095148) q[2];
sx q[2];
rz(-2.6161604) q[2];
sx q[2];
rz(-2.5642774) q[2];
rz(0.32140857) q[3];
sx q[3];
rz(-1.0151498) q[3];
sx q[3];
rz(-2.2890454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];