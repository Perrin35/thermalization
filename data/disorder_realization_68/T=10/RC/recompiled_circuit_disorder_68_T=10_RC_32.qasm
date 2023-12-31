OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7473937) q[0];
sx q[0];
rz(-2.6497901) q[0];
sx q[0];
rz(2.9536182) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(4.6586577) q[1];
sx q[1];
rz(12.933856) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4386908) q[0];
sx q[0];
rz(-0.54649788) q[0];
sx q[0];
rz(-1.7706857) q[0];
rz(-pi) q[1];
rz(-0.1349749) q[2];
sx q[2];
rz(-1.0834603) q[2];
sx q[2];
rz(-2.6589573) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.623466) q[1];
sx q[1];
rz(-1.7382442) q[1];
sx q[1];
rz(-1.3256339) q[1];
x q[2];
rz(2.8005881) q[3];
sx q[3];
rz(-1.4031646) q[3];
sx q[3];
rz(1.0238907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.964103) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(-0.5509848) q[2];
rz(1.3059113) q[3];
sx q[3];
rz(-1.4923613) q[3];
sx q[3];
rz(1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630163) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(2.6696894) q[0];
rz(2.7117803) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(-2.205251) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72915709) q[0];
sx q[0];
rz(-1.4834187) q[0];
sx q[0];
rz(2.9108414) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.721644) q[2];
sx q[2];
rz(-0.83050767) q[2];
sx q[2];
rz(-0.74479693) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4221103) q[1];
sx q[1];
rz(-1.4665831) q[1];
sx q[1];
rz(0.69394333) q[1];
rz(-pi) q[2];
rz(1.6132658) q[3];
sx q[3];
rz(-0.50308933) q[3];
sx q[3];
rz(-2.344362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3669746) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(-0.42638391) q[2];
rz(1.9042227) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(-0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8957829) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(-0.93908969) q[0];
rz(-0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(2.5476707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.827841) q[0];
sx q[0];
rz(-1.272164) q[0];
sx q[0];
rz(1.2350425) q[0];
rz(-2.2592696) q[2];
sx q[2];
rz(-1.6154628) q[2];
sx q[2];
rz(-0.67827144) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18332874) q[1];
sx q[1];
rz(-1.3576344) q[1];
sx q[1];
rz(1.1413241) q[1];
rz(-1.2265615) q[3];
sx q[3];
rz(-1.1262745) q[3];
sx q[3];
rz(1.3821186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5014191) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(-1.4397941) q[2];
rz(-2.7539608) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(0.50278062) q[0];
rz(-0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(-2.3847413) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6605646) q[0];
sx q[0];
rz(-1.2111944) q[0];
sx q[0];
rz(1.2147551) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1286131) q[2];
sx q[2];
rz(-2.0364967) q[2];
sx q[2];
rz(-2.4109858) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56843578) q[1];
sx q[1];
rz(-0.70306289) q[1];
sx q[1];
rz(2.1431124) q[1];
x q[2];
rz(-2.0714949) q[3];
sx q[3];
rz(-2.2399733) q[3];
sx q[3];
rz(-2.5239528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7148774) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(1.4871917) q[2];
rz(0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(-0.55707651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-2.3838682) q[0];
rz(1.2879397) q[1];
sx q[1];
rz(-0.92823354) q[1];
sx q[1];
rz(2.0910738) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61705631) q[0];
sx q[0];
rz(-0.40818383) q[0];
sx q[0];
rz(2.257686) q[0];
rz(-0.97425766) q[2];
sx q[2];
rz(-1.9447118) q[2];
sx q[2];
rz(-2.8982382) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1041303) q[1];
sx q[1];
rz(-0.81172746) q[1];
sx q[1];
rz(2.3292259) q[1];
rz(-2.1220783) q[3];
sx q[3];
rz(-0.84421221) q[3];
sx q[3];
rz(0.36908484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22333764) q[2];
sx q[2];
rz(-2.788322) q[2];
sx q[2];
rz(0.63344947) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.6703689) q[3];
sx q[3];
rz(-0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.441992) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(-0.25318405) q[0];
rz(-1.6075915) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(-1.6794499) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7397241) q[0];
sx q[0];
rz(-1.7812294) q[0];
sx q[0];
rz(-1.6519288) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2573104) q[2];
sx q[2];
rz(-2.4682211) q[2];
sx q[2];
rz(1.1416669) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.05789214) q[1];
sx q[1];
rz(-1.9533227) q[1];
sx q[1];
rz(1.3871357) q[1];
rz(-pi) q[2];
rz(-1.2268279) q[3];
sx q[3];
rz(-1.2729537) q[3];
sx q[3];
rz(-2.1732268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9399461) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(-0.80491006) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8686304) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(2.2139363) q[0];
rz(-1.0246798) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(-1.0120846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18400684) q[0];
sx q[0];
rz(-1.1054966) q[0];
sx q[0];
rz(-2.1561949) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8132642) q[2];
sx q[2];
rz(-2.7331181) q[2];
sx q[2];
rz(-0.35818737) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.183179) q[1];
sx q[1];
rz(-1.5378386) q[1];
sx q[1];
rz(-0.54380137) q[1];
rz(-3.0104962) q[3];
sx q[3];
rz(-0.56250611) q[3];
sx q[3];
rz(-1.8545811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6107789) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(1.0144462) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(-0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5291418) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(-0.21417831) q[0];
rz(-2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(-2.8578551) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8003214) q[0];
sx q[0];
rz(-0.30297908) q[0];
sx q[0];
rz(0.11462258) q[0];
x q[1];
rz(2.0869414) q[2];
sx q[2];
rz(-0.37545855) q[2];
sx q[2];
rz(1.6106538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.46854308) q[1];
sx q[1];
rz(-1.0769516) q[1];
sx q[1];
rz(1.3480575) q[1];
rz(-2.489336) q[3];
sx q[3];
rz(-1.0271003) q[3];
sx q[3];
rz(-3.101845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6909137) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(-1.696375) q[2];
rz(1.5971659) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(-0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63672367) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(-3.0723363) q[0];
rz(1.4878558) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(-1.5690631) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78281392) q[0];
sx q[0];
rz(-0.94952119) q[0];
sx q[0];
rz(-0.57399477) q[0];
rz(2.9022129) q[2];
sx q[2];
rz(-2.8068672) q[2];
sx q[2];
rz(-0.3604381) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.20977565) q[1];
sx q[1];
rz(-1.4869542) q[1];
sx q[1];
rz(0.2340338) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2980372) q[3];
sx q[3];
rz(-2.4747304) q[3];
sx q[3];
rz(-0.71344261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9562324) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(0.13988477) q[2];
rz(2.774003) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(-0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(0.64176732) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(0.26783255) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3788911) q[0];
sx q[0];
rz(-2.2305616) q[0];
sx q[0];
rz(0.99878175) q[0];
rz(-pi) q[1];
rz(0.12561663) q[2];
sx q[2];
rz(-1.7526502) q[2];
sx q[2];
rz(0.4609209) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8903058) q[1];
sx q[1];
rz(-0.44154134) q[1];
sx q[1];
rz(0.73015405) q[1];
x q[2];
rz(2.855741) q[3];
sx q[3];
rz(-1.0470069) q[3];
sx q[3];
rz(1.7459735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(1.4769185) q[2];
rz(0.26633513) q[3];
sx q[3];
rz(-2.895152) q[3];
sx q[3];
rz(-0.5464856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1289566) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(0.77990445) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(0.71074625) q[2];
sx q[2];
rz(-1.1087316) q[2];
sx q[2];
rz(1.0089594) q[2];
rz(0.26364636) q[3];
sx q[3];
rz(-0.69477889) q[3];
sx q[3];
rz(-0.033332326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
