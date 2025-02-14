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
rz(-3.0946331) q[0];
sx q[0];
rz(-1.5927915) q[0];
sx q[0];
rz(-1.3943075) q[0];
rz(2.6810763) q[1];
sx q[1];
rz(-1.2682275) q[1];
sx q[1];
rz(-2.6654798) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5232485) q[0];
sx q[0];
rz(-1.9501513) q[0];
sx q[0];
rz(0.98486395) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9063375) q[2];
sx q[2];
rz(-1.3089184) q[2];
sx q[2];
rz(2.1507598) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.70041768) q[1];
sx q[1];
rz(-2.9923424) q[1];
sx q[1];
rz(3.0249437) q[1];
rz(1.634811) q[3];
sx q[3];
rz(-1.598017) q[3];
sx q[3];
rz(2.343319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8146879) q[2];
sx q[2];
rz(-1.2219595) q[2];
sx q[2];
rz(-1.2828705) q[2];
rz(-0.13947105) q[3];
sx q[3];
rz(-1.3160416) q[3];
sx q[3];
rz(-0.68525806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71444702) q[0];
sx q[0];
rz(-0.35816631) q[0];
sx q[0];
rz(-2.8156679) q[0];
rz(2.4320995) q[1];
sx q[1];
rz(-2.8804417) q[1];
sx q[1];
rz(2.479898) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0950006) q[0];
sx q[0];
rz(-1.5803156) q[0];
sx q[0];
rz(-0.0074038238) q[0];
rz(-pi) q[1];
rz(-0.93962713) q[2];
sx q[2];
rz(-1.1568312) q[2];
sx q[2];
rz(-0.30218538) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0256309) q[1];
sx q[1];
rz(-0.44646663) q[1];
sx q[1];
rz(-1.0122416) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68485188) q[3];
sx q[3];
rz(-1.0954482) q[3];
sx q[3];
rz(2.1641987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0207396) q[2];
sx q[2];
rz(-1.8031305) q[2];
sx q[2];
rz(2.8815114) q[2];
rz(2.2777879) q[3];
sx q[3];
rz(-1.443202) q[3];
sx q[3];
rz(-1.1687733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.9048555) q[0];
sx q[0];
rz(-2.3598292) q[0];
sx q[0];
rz(2.8535063) q[0];
rz(1.9347363) q[1];
sx q[1];
rz(-2.0286045) q[1];
sx q[1];
rz(-0.77817121) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.208565) q[0];
sx q[0];
rz(-1.7598682) q[0];
sx q[0];
rz(-2.2651432) q[0];
rz(-pi) q[1];
rz(1.6323998) q[2];
sx q[2];
rz(-2.3280316) q[2];
sx q[2];
rz(0.98682994) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6699804) q[1];
sx q[1];
rz(-2.0005352) q[1];
sx q[1];
rz(0.55008908) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0707079) q[3];
sx q[3];
rz(-2.4307561) q[3];
sx q[3];
rz(-0.58037607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8975767) q[2];
sx q[2];
rz(-1.1268758) q[2];
sx q[2];
rz(-2.8934532) q[2];
rz(-2.1085619) q[3];
sx q[3];
rz(-1.4890198) q[3];
sx q[3];
rz(1.2584794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7018062) q[0];
sx q[0];
rz(-2.2070364) q[0];
sx q[0];
rz(-2.5820861) q[0];
rz(-1.4065546) q[1];
sx q[1];
rz(-0.94466698) q[1];
sx q[1];
rz(-1.9962126) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6481768) q[0];
sx q[0];
rz(-0.45873612) q[0];
sx q[0];
rz(3.0381195) q[0];
rz(-pi) q[1];
rz(-0.65088314) q[2];
sx q[2];
rz(-2.0325629) q[2];
sx q[2];
rz(-0.24591638) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2753452) q[1];
sx q[1];
rz(-1.8653578) q[1];
sx q[1];
rz(-0.7742851) q[1];
x q[2];
rz(3.1265852) q[3];
sx q[3];
rz(-1.0010442) q[3];
sx q[3];
rz(0.80591312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8765325) q[2];
sx q[2];
rz(-2.127485) q[2];
sx q[2];
rz(1.8939135) q[2];
rz(2.0606591) q[3];
sx q[3];
rz(-1.7526046) q[3];
sx q[3];
rz(-1.8814253) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0697698) q[0];
sx q[0];
rz(-2.6520196) q[0];
sx q[0];
rz(-2.3967632) q[0];
rz(-1.3818332) q[1];
sx q[1];
rz(-2.0227183) q[1];
sx q[1];
rz(0.098085731) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25135219) q[0];
sx q[0];
rz(-2.6664263) q[0];
sx q[0];
rz(-0.71936468) q[0];
rz(-0.63348153) q[2];
sx q[2];
rz(-2.6085662) q[2];
sx q[2];
rz(0.60163762) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87865767) q[1];
sx q[1];
rz(-2.2644357) q[1];
sx q[1];
rz(1.5469815) q[1];
x q[2];
rz(-1.0897343) q[3];
sx q[3];
rz(-1.4816545) q[3];
sx q[3];
rz(-1.7952023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7006526) q[2];
sx q[2];
rz(-2.7589189) q[2];
sx q[2];
rz(2.0514226) q[2];
rz(2.8160461) q[3];
sx q[3];
rz(-1.6306337) q[3];
sx q[3];
rz(2.8946099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22661041) q[0];
sx q[0];
rz(-2.60422) q[0];
sx q[0];
rz(-1.2363303) q[0];
rz(-0.63713282) q[1];
sx q[1];
rz(-0.56012154) q[1];
sx q[1];
rz(-2.2083652) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0765981) q[0];
sx q[0];
rz(-0.68474283) q[0];
sx q[0];
rz(0.97692722) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9524589) q[2];
sx q[2];
rz(-1.1580843) q[2];
sx q[2];
rz(0.37278644) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5781587) q[1];
sx q[1];
rz(-1.8616042) q[1];
sx q[1];
rz(-1.4842503) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1307484) q[3];
sx q[3];
rz(-2.5820316) q[3];
sx q[3];
rz(-0.51792972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.25202641) q[2];
sx q[2];
rz(-1.3875952) q[2];
sx q[2];
rz(1.8760366) q[2];
rz(-1.4243852) q[3];
sx q[3];
rz(-2.6129183) q[3];
sx q[3];
rz(-2.2831447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.76765656) q[0];
sx q[0];
rz(-2.7601384) q[0];
sx q[0];
rz(0.23442991) q[0];
rz(-0.43342057) q[1];
sx q[1];
rz(-1.0973009) q[1];
sx q[1];
rz(-2.0030599) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075468242) q[0];
sx q[0];
rz(-0.99695092) q[0];
sx q[0];
rz(-2.2908151) q[0];
rz(-pi) q[1];
rz(1.578733) q[2];
sx q[2];
rz(-2.0614377) q[2];
sx q[2];
rz(2.454941) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2039472) q[1];
sx q[1];
rz(-1.0685941) q[1];
sx q[1];
rz(0.03003386) q[1];
rz(-1.7650003) q[3];
sx q[3];
rz(-2.248431) q[3];
sx q[3];
rz(-1.6226131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66702691) q[2];
sx q[2];
rz(-0.67956769) q[2];
sx q[2];
rz(2.1843074) q[2];
rz(-0.75753093) q[3];
sx q[3];
rz(-1.6672971) q[3];
sx q[3];
rz(0.1434513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.8677419) q[0];
sx q[0];
rz(-2.0389281) q[0];
sx q[0];
rz(-0.17909166) q[0];
rz(0.20009072) q[1];
sx q[1];
rz(-0.6691907) q[1];
sx q[1];
rz(2.3650513) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18416737) q[0];
sx q[0];
rz(-1.3901852) q[0];
sx q[0];
rz(-2.1023048) q[0];
x q[1];
rz(-0.9540117) q[2];
sx q[2];
rz(-1.8158834) q[2];
sx q[2];
rz(-0.38604647) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.071826) q[1];
sx q[1];
rz(-0.36301314) q[1];
sx q[1];
rz(-2.9337204) q[1];
rz(-pi) q[2];
rz(-0.051335585) q[3];
sx q[3];
rz(-1.253479) q[3];
sx q[3];
rz(0.3454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.17579235) q[2];
sx q[2];
rz(-2.0696023) q[2];
sx q[2];
rz(-1.3624066) q[2];
rz(-1.8386748) q[3];
sx q[3];
rz(-1.4785654) q[3];
sx q[3];
rz(1.038507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4565444) q[0];
sx q[0];
rz(-1.9467204) q[0];
sx q[0];
rz(-0.1142647) q[0];
rz(-1.9396797) q[1];
sx q[1];
rz(-2.6004531) q[1];
sx q[1];
rz(0.094954403) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9499493) q[0];
sx q[0];
rz(-1.7472634) q[0];
sx q[0];
rz(1.0383738) q[0];
x q[1];
rz(2.6812115) q[2];
sx q[2];
rz(-2.2308439) q[2];
sx q[2];
rz(0.4649064) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11647955) q[1];
sx q[1];
rz(-1.0214865) q[1];
sx q[1];
rz(1.504937) q[1];
rz(0.41422959) q[3];
sx q[3];
rz(-1.0981907) q[3];
sx q[3];
rz(0.35652682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.3007043) q[2];
sx q[2];
rz(-1.3450832) q[2];
sx q[2];
rz(-0.67187205) q[2];
rz(-2.4554409) q[3];
sx q[3];
rz(-2.4785564) q[3];
sx q[3];
rz(1.9239976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9753863) q[0];
sx q[0];
rz(-0.19151846) q[0];
sx q[0];
rz(-1.5331049) q[0];
rz(0.32471049) q[1];
sx q[1];
rz(-1.277781) q[1];
sx q[1];
rz(0.3130354) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9444549) q[0];
sx q[0];
rz(-2.1416683) q[0];
sx q[0];
rz(-3.0464493) q[0];
rz(-pi) q[1];
rz(1.7907421) q[2];
sx q[2];
rz(-1.6306021) q[2];
sx q[2];
rz(-2.4129197) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5165007) q[1];
sx q[1];
rz(-0.91462574) q[1];
sx q[1];
rz(2.28338) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62449091) q[3];
sx q[3];
rz(-1.1905128) q[3];
sx q[3];
rz(0.46562574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2031871) q[2];
sx q[2];
rz(-1.9878191) q[2];
sx q[2];
rz(-1.3333092) q[2];
rz(-2.5238254) q[3];
sx q[3];
rz(-0.94527644) q[3];
sx q[3];
rz(-2.0060523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4278605) q[0];
sx q[0];
rz(-1.328631) q[0];
sx q[0];
rz(-2.791639) q[0];
rz(-0.89823828) q[1];
sx q[1];
rz(-1.0844834) q[1];
sx q[1];
rz(-0.35379298) q[1];
rz(0.66302115) q[2];
sx q[2];
rz(-0.86227476) q[2];
sx q[2];
rz(-3.0515565) q[2];
rz(-3.114936) q[3];
sx q[3];
rz(-2.4116357) q[3];
sx q[3];
rz(0.38105376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
