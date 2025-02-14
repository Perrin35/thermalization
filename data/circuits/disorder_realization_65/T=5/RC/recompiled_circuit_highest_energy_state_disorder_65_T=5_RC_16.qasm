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
rz(1.7841568) q[0];
sx q[0];
rz(-0.62499243) q[0];
sx q[0];
rz(1.5224737) q[0];
rz(1.1433262) q[1];
sx q[1];
rz(3.7781236) q[1];
sx q[1];
rz(14.349024) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56464897) q[0];
sx q[0];
rz(-1.859382) q[0];
sx q[0];
rz(-0.99228575) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8463616) q[2];
sx q[2];
rz(-1.2157939) q[2];
sx q[2];
rz(0.4685185) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7984109) q[1];
sx q[1];
rz(-0.34920563) q[1];
sx q[1];
rz(0.83425547) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3728849) q[3];
sx q[3];
rz(-0.48600125) q[3];
sx q[3];
rz(-0.81102351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6978567) q[2];
sx q[2];
rz(-0.82252684) q[2];
sx q[2];
rz(-2.6653384) q[2];
rz(3.0016628) q[3];
sx q[3];
rz(-1.6401451) q[3];
sx q[3];
rz(-3.068315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3530465) q[0];
sx q[0];
rz(-1.6206425) q[0];
sx q[0];
rz(0.36343685) q[0];
rz(-0.67047554) q[1];
sx q[1];
rz(-1.3251708) q[1];
sx q[1];
rz(-2.061981) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9417861) q[0];
sx q[0];
rz(-2.5528756) q[0];
sx q[0];
rz(-2.0090282) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7893546) q[2];
sx q[2];
rz(-2.5982476) q[2];
sx q[2];
rz(-1.3308412) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0091925) q[1];
sx q[1];
rz(-1.8120488) q[1];
sx q[1];
rz(-3.0962837) q[1];
rz(-pi) q[2];
rz(2.5609301) q[3];
sx q[3];
rz(-1.7911909) q[3];
sx q[3];
rz(3.0739258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67538992) q[2];
sx q[2];
rz(-2.4066996) q[2];
sx q[2];
rz(-0.43851635) q[2];
rz(1.0112666) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(-1.6605759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587875) q[0];
sx q[0];
rz(-0.52849448) q[0];
sx q[0];
rz(0.82861376) q[0];
rz(-1.2476745) q[1];
sx q[1];
rz(-1.1670466) q[1];
sx q[1];
rz(-0.25965986) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8369401) q[0];
sx q[0];
rz(-2.2373709) q[0];
sx q[0];
rz(-2.7949692) q[0];
x q[1];
rz(-2.3724062) q[2];
sx q[2];
rz(-2.63317) q[2];
sx q[2];
rz(1.2505797) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7938215) q[1];
sx q[1];
rz(-1.1060556) q[1];
sx q[1];
rz(-0.79885599) q[1];
rz(-2.2237556) q[3];
sx q[3];
rz(-1.8052832) q[3];
sx q[3];
rz(-0.10374903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9344249) q[2];
sx q[2];
rz(-1.5639037) q[2];
sx q[2];
rz(0.820532) q[2];
rz(-2.3151243) q[3];
sx q[3];
rz(-2.1491437) q[3];
sx q[3];
rz(-1.3124527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8823223) q[0];
sx q[0];
rz(-0.81040183) q[0];
sx q[0];
rz(0.93233863) q[0];
rz(1.8238292) q[1];
sx q[1];
rz(-1.980314) q[1];
sx q[1];
rz(3.0188149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2185635) q[0];
sx q[0];
rz(-1.4866795) q[0];
sx q[0];
rz(0.53972678) q[0];
rz(0.30722801) q[2];
sx q[2];
rz(-2.3005552) q[2];
sx q[2];
rz(2.407069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40298324) q[1];
sx q[1];
rz(-1.5697641) q[1];
sx q[1];
rz(0.82768517) q[1];
rz(-pi) q[2];
rz(1.2776385) q[3];
sx q[3];
rz(-2.1159626) q[3];
sx q[3];
rz(2.4923862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5968898) q[2];
sx q[2];
rz(-0.92146102) q[2];
sx q[2];
rz(-2.5370562) q[2];
rz(-2.4321411) q[3];
sx q[3];
rz(-0.76990288) q[3];
sx q[3];
rz(-2.3179222) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9619047) q[0];
sx q[0];
rz(-3.0178495) q[0];
sx q[0];
rz(-1.3667579) q[0];
rz(2.1429515) q[1];
sx q[1];
rz(-2.3256681) q[1];
sx q[1];
rz(-2.6008115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3872469) q[0];
sx q[0];
rz(-1.7779113) q[0];
sx q[0];
rz(1.7133173) q[0];
rz(0.86414637) q[2];
sx q[2];
rz(-2.1808778) q[2];
sx q[2];
rz(-0.44877258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4879228) q[1];
sx q[1];
rz(-1.2752295) q[1];
sx q[1];
rz(1.771677) q[1];
rz(-pi) q[2];
x q[2];
rz(1.389796) q[3];
sx q[3];
rz(-1.01184) q[3];
sx q[3];
rz(-1.4829092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6775386) q[2];
sx q[2];
rz(-0.81391922) q[2];
sx q[2];
rz(-2.3133254) q[2];
rz(1.6895435) q[3];
sx q[3];
rz(-1.6477081) q[3];
sx q[3];
rz(2.233708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0238817) q[0];
sx q[0];
rz(-1.9458867) q[0];
sx q[0];
rz(1.1718132) q[0];
rz(-0.68012971) q[1];
sx q[1];
rz(-1.9155733) q[1];
sx q[1];
rz(-0.5562869) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.507001) q[0];
sx q[0];
rz(-1.9456353) q[0];
sx q[0];
rz(-0.27347538) q[0];
rz(-1.7173217) q[2];
sx q[2];
rz(-1.5735448) q[2];
sx q[2];
rz(0.78080746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7032675) q[1];
sx q[1];
rz(-0.57619737) q[1];
sx q[1];
rz(0.11771113) q[1];
rz(-pi) q[2];
rz(-0.63594922) q[3];
sx q[3];
rz(-1.1750286) q[3];
sx q[3];
rz(-1.1102939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80643001) q[2];
sx q[2];
rz(-1.282629) q[2];
sx q[2];
rz(0.20509091) q[2];
rz(-2.3388376) q[3];
sx q[3];
rz(-2.0358678) q[3];
sx q[3];
rz(-0.55454379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270808) q[0];
sx q[0];
rz(-2.0682122) q[0];
sx q[0];
rz(-0.22739534) q[0];
rz(-0.79044509) q[1];
sx q[1];
rz(-2.3643654) q[1];
sx q[1];
rz(2.3048185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8604383) q[0];
sx q[0];
rz(-1.3520762) q[0];
sx q[0];
rz(-0.52978306) q[0];
rz(-8/(5*pi)) q[2];
sx q[2];
rz(-0.65213883) q[2];
sx q[2];
rz(0.76667128) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6001125) q[1];
sx q[1];
rz(-1.8842234) q[1];
sx q[1];
rz(-0.42509834) q[1];
x q[2];
rz(-0.38897377) q[3];
sx q[3];
rz(-2.3416061) q[3];
sx q[3];
rz(-0.98158132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5493912) q[2];
sx q[2];
rz(-1.6624007) q[2];
sx q[2];
rz(-2.5679892) q[2];
rz(2.8163689) q[3];
sx q[3];
rz(-2.2461788) q[3];
sx q[3];
rz(2.7371244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0322872) q[0];
sx q[0];
rz(-2.9635297) q[0];
sx q[0];
rz(0.67894116) q[0];
rz(-0.13488723) q[1];
sx q[1];
rz(-2.0811681) q[1];
sx q[1];
rz(1.7178242) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2436062) q[0];
sx q[0];
rz(-2.4084414) q[0];
sx q[0];
rz(2.2871458) q[0];
x q[1];
rz(-2.277454) q[2];
sx q[2];
rz(-1.1306964) q[2];
sx q[2];
rz(0.76264438) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1610802) q[1];
sx q[1];
rz(-0.80787611) q[1];
sx q[1];
rz(3.1001311) q[1];
rz(-1.0342977) q[3];
sx q[3];
rz(-2.6141659) q[3];
sx q[3];
rz(2.3303383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5474825) q[2];
sx q[2];
rz(-2.2308733) q[2];
sx q[2];
rz(-2.9098848) q[2];
rz(-1.0101275) q[3];
sx q[3];
rz(-1.908327) q[3];
sx q[3];
rz(0.89099425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6758839) q[0];
sx q[0];
rz(-2.3864585) q[0];
sx q[0];
rz(0.51374197) q[0];
rz(-1.4328009) q[1];
sx q[1];
rz(-1.2756196) q[1];
sx q[1];
rz(-1.6339711) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20376539) q[0];
sx q[0];
rz(-1.4368334) q[0];
sx q[0];
rz(1.521827) q[0];
rz(-pi) q[1];
rz(-2.5812923) q[2];
sx q[2];
rz(-2.8447084) q[2];
sx q[2];
rz(2.1765944) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3097174) q[1];
sx q[1];
rz(-0.17942218) q[1];
sx q[1];
rz(-2.37876) q[1];
x q[2];
rz(-0.024941878) q[3];
sx q[3];
rz(-1.8224622) q[3];
sx q[3];
rz(-0.97398678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.12019176) q[2];
sx q[2];
rz(-2.2970565) q[2];
sx q[2];
rz(-2.8821442) q[2];
rz(1.1546968) q[3];
sx q[3];
rz(-2.3754933) q[3];
sx q[3];
rz(-1.8416789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65086377) q[0];
sx q[0];
rz(-1.6614953) q[0];
sx q[0];
rz(1.65253) q[0];
rz(-0.64000714) q[1];
sx q[1];
rz(-1.0970486) q[1];
sx q[1];
rz(-2.1515813) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3418808) q[0];
sx q[0];
rz(-0.16946259) q[0];
sx q[0];
rz(-2.5925596) q[0];
x q[1];
rz(2.6578519) q[2];
sx q[2];
rz(-2.3123398) q[2];
sx q[2];
rz(-0.26243789) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4943018) q[1];
sx q[1];
rz(-1.1539946) q[1];
sx q[1];
rz(2.9652262) q[1];
x q[2];
rz(0.02496533) q[3];
sx q[3];
rz(-1.348266) q[3];
sx q[3];
rz(2.3937267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9330357) q[2];
sx q[2];
rz(-1.8427589) q[2];
sx q[2];
rz(2.9162858) q[2];
rz(2.2202668) q[3];
sx q[3];
rz(-0.39042979) q[3];
sx q[3];
rz(-0.050617378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13023547) q[0];
sx q[0];
rz(-2.7738032) q[0];
sx q[0];
rz(-2.2875447) q[0];
rz(1.3182974) q[1];
sx q[1];
rz(-1.5250991) q[1];
sx q[1];
rz(-0.93223882) q[1];
rz(-3.0443424) q[2];
sx q[2];
rz(-1.8293867) q[2];
sx q[2];
rz(0.23501227) q[2];
rz(2.9574838) q[3];
sx q[3];
rz(-1.0297903) q[3];
sx q[3];
rz(-0.31576706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
