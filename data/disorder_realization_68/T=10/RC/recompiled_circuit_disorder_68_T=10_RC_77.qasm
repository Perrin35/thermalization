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
rz(-1.1176874) q[1];
sx q[1];
rz(-1.517065) q[1];
sx q[1];
rz(2.7741073) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6715235) q[0];
sx q[0];
rz(-1.0363665) q[0];
sx q[0];
rz(3.021391) q[0];
x q[1];
rz(0.1349749) q[2];
sx q[2];
rz(-2.0581323) q[2];
sx q[2];
rz(0.48263532) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6403113) q[1];
sx q[1];
rz(-0.2959364) q[1];
sx q[1];
rz(-2.179115) q[1];
rz(-pi) q[2];
rz(0.46842694) q[3];
sx q[3];
rz(-2.763063) q[3];
sx q[3];
rz(-0.98640501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1774896) q[2];
sx q[2];
rz(-2.6314645) q[2];
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
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47857639) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(-0.4719032) q[0];
rz(0.42981237) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(0.93634161) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443003) q[0];
sx q[0];
rz(-2.8951277) q[0];
sx q[0];
rz(0.36578567) q[0];
rz(-pi) q[1];
rz(-1.1510552) q[2];
sx q[2];
rz(-2.310576) q[2];
sx q[2];
rz(1.8120399) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3766107) q[1];
sx q[1];
rz(-0.88135834) q[1];
sx q[1];
rz(-1.4355684) q[1];
rz(2.0735047) q[3];
sx q[3];
rz(-1.5912676) q[3];
sx q[3];
rz(0.81077829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3669746) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(2.7152087) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24580978) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(-0.93908969) q[0];
rz(0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(0.59392196) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31375162) q[0];
sx q[0];
rz(-1.272164) q[0];
sx q[0];
rz(-1.2350425) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8823231) q[2];
sx q[2];
rz(-1.5261298) q[2];
sx q[2];
rz(-2.4633212) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3214026) q[1];
sx q[1];
rz(-2.6650975) q[1];
sx q[1];
rz(-2.0501775) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46843723) q[3];
sx q[3];
rz(-1.2611946) q[3];
sx q[3];
rz(2.7999511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5014191) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(1.4397941) q[2];
rz(-2.7539608) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(0.50278062) q[0];
rz(0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(-0.75685135) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6605646) q[0];
sx q[0];
rz(-1.9303983) q[0];
sx q[0];
rz(-1.9268376) q[0];
rz(-pi) q[1];
rz(-1.5966162) q[2];
sx q[2];
rz(-2.6757247) q[2];
sx q[2];
rz(-0.70170882) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56843578) q[1];
sx q[1];
rz(-2.4385298) q[1];
sx q[1];
rz(2.1431124) q[1];
x q[2];
rz(2.0714949) q[3];
sx q[3];
rz(-2.2399733) q[3];
sx q[3];
rz(-0.6176399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7148774) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(-1.4871917) q[2];
rz(-0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(0.55707651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-2.0563545) q[0];
sx q[0];
rz(-2.3838682) q[0];
rz(1.853653) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(-1.0505189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8331063) q[0];
sx q[0];
rz(-1.8252488) q[0];
sx q[0];
rz(-1.893505) q[0];
rz(-pi) q[1];
x q[1];
rz(2.167335) q[2];
sx q[2];
rz(-1.9447118) q[2];
sx q[2];
rz(-2.8982382) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1615636) q[1];
sx q[1];
rz(-1.0161576) q[1];
sx q[1];
rz(-2.5142923) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1220783) q[3];
sx q[3];
rz(-2.2973804) q[3];
sx q[3];
rz(0.36908484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22333764) q[2];
sx q[2];
rz(-2.788322) q[2];
sx q[2];
rz(-0.63344947) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.6703689) q[3];
sx q[3];
rz(2.4244394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69960064) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(2.8884086) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(1.4621428) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1519449) q[0];
sx q[0];
rz(-1.6501353) q[0];
sx q[0];
rz(-2.9304855) q[0];
rz(-pi) q[1];
rz(1.8842823) q[2];
sx q[2];
rz(-0.67337155) q[2];
sx q[2];
rz(-1.1416669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0837005) q[1];
sx q[1];
rz(-1.9533227) q[1];
sx q[1];
rz(1.3871357) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2268279) q[3];
sx q[3];
rz(-1.2729537) q[3];
sx q[3];
rz(2.1732268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9399461) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(0.80491006) q[2];
rz(-1.9645875) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(-2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8686304) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(-0.92765635) q[0];
rz(1.0246798) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(2.129508) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9818078) q[0];
sx q[0];
rz(-2.4111528) q[0];
sx q[0];
rz(-0.83321379) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32832844) q[2];
sx q[2];
rz(-2.7331181) q[2];
sx q[2];
rz(-2.7834053) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5092897) q[1];
sx q[1];
rz(-1.0273233) q[1];
sx q[1];
rz(-1.6093045) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4885694) q[3];
sx q[3];
rz(-2.1279018) q[3];
sx q[3];
rz(1.132387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53081375) q[2];
sx q[2];
rz(-1.4799708) q[2];
sx q[2];
rz(-2.7116595) q[2];
rz(-2.1271465) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5291418) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(2.9274143) q[0];
rz(2.0902436) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(0.28373757) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2212299) q[0];
sx q[0];
rz(-1.8717248) q[0];
sx q[0];
rz(1.6065341) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2405422) q[2];
sx q[2];
rz(-1.7527765) q[2];
sx q[2];
rz(2.6957126) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
x q[2];
rz(2.489336) q[3];
sx q[3];
rz(-2.1144923) q[3];
sx q[3];
rz(-3.101845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.45067898) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(-1.696375) q[2];
rz(-1.5444267) q[3];
sx q[3];
rz(-1.701136) q[3];
sx q[3];
rz(0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(-2.504869) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(-0.069256393) q[0];
rz(1.4878558) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(-1.5725296) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6206966) q[0];
sx q[0];
rz(-2.3224152) q[0];
sx q[0];
rz(0.92185123) q[0];
rz(-0.23937978) q[2];
sx q[2];
rz(-2.8068672) q[2];
sx q[2];
rz(-0.3604381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.20977565) q[1];
sx q[1];
rz(-1.4869542) q[1];
sx q[1];
rz(2.9075588) q[1];
rz(-pi) q[2];
rz(-0.2089573) q[3];
sx q[3];
rz(-0.93271241) q[3];
sx q[3];
rz(2.0861422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1853603) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(2.774003) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(2.1504413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-2.4998253) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(-0.26783255) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5689241) q[0];
sx q[0];
rz(-0.84416443) q[0];
sx q[0];
rz(0.6092682) q[0];
rz(1.3875302) q[2];
sx q[2];
rz(-1.4472618) q[2];
sx q[2];
rz(1.087041) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8903058) q[1];
sx q[1];
rz(-0.44154134) q[1];
sx q[1];
rz(0.73015405) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28585163) q[3];
sx q[3];
rz(-2.0945858) q[3];
sx q[3];
rz(1.3956192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(-1.6646741) q[2];
rz(-0.26633513) q[3];
sx q[3];
rz(-2.895152) q[3];
sx q[3];
rz(0.5464856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289566) q[0];
sx q[0];
rz(-2.2205882) q[0];
sx q[0];
rz(0.90482774) q[0];
rz(-0.77990445) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(-2.1521679) q[2];
sx q[2];
rz(-0.94716723) q[2];
sx q[2];
rz(-0.92826044) q[2];
rz(-0.6775425) q[3];
sx q[3];
rz(-1.7384221) q[3];
sx q[3];
rz(1.3330028) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];