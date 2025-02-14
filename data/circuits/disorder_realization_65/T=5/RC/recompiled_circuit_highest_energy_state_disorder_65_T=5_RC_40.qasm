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
rz(-1.3574358) q[0];
sx q[0];
rz(3.7665851) q[0];
sx q[0];
rz(7.9023043) q[0];
rz(1.1433262) q[1];
sx q[1];
rz(3.7781236) q[1];
sx q[1];
rz(14.349024) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5464294) q[0];
sx q[0];
rz(-2.5025285) q[0];
sx q[0];
rz(-2.0682113) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29523103) q[2];
sx q[2];
rz(-1.9257987) q[2];
sx q[2];
rz(0.4685185) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5660043) q[1];
sx q[1];
rz(-1.3145361) q[1];
sx q[1];
rz(-2.9017067) q[1];
rz(-pi) q[2];
rz(-2.0487011) q[3];
sx q[3];
rz(-1.6627668) q[3];
sx q[3];
rz(2.2063279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4437359) q[2];
sx q[2];
rz(-2.3190658) q[2];
sx q[2];
rz(-0.47625429) q[2];
rz(0.13992986) q[3];
sx q[3];
rz(-1.5014476) q[3];
sx q[3];
rz(-3.068315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78854617) q[0];
sx q[0];
rz(-1.6206425) q[0];
sx q[0];
rz(2.7781558) q[0];
rz(2.4711171) q[1];
sx q[1];
rz(-1.8164219) q[1];
sx q[1];
rz(-1.0796116) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1422259) q[0];
sx q[0];
rz(-1.3329263) q[0];
sx q[0];
rz(2.1145942) q[0];
x q[1];
rz(-0.13021664) q[2];
sx q[2];
rz(-1.0417582) q[2];
sx q[2];
rz(1.0769653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9448589) q[1];
sx q[1];
rz(-2.8962039) q[1];
sx q[1];
rz(1.3887482) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3877181) q[3];
sx q[3];
rz(-0.61655515) q[3];
sx q[3];
rz(-1.316837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67538992) q[2];
sx q[2];
rz(-0.73489302) q[2];
sx q[2];
rz(0.43851635) q[2];
rz(1.0112666) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(1.4810168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48280516) q[0];
sx q[0];
rz(-2.6130982) q[0];
sx q[0];
rz(-2.3129789) q[0];
rz(-1.2476745) q[1];
sx q[1];
rz(-1.1670466) q[1];
sx q[1];
rz(-0.25965986) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8369401) q[0];
sx q[0];
rz(-2.2373709) q[0];
sx q[0];
rz(0.34662342) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7607418) q[2];
sx q[2];
rz(-1.9162188) q[2];
sx q[2];
rz(2.7598515) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79154014) q[1];
sx q[1];
rz(-0.87557057) q[1];
sx q[1];
rz(0.94757845) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8494495) q[3];
sx q[3];
rz(-2.202987) q[3];
sx q[3];
rz(1.8504253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20716771) q[2];
sx q[2];
rz(-1.577689) q[2];
sx q[2];
rz(-2.3210607) q[2];
rz(-0.82646838) q[3];
sx q[3];
rz(-2.1491437) q[3];
sx q[3];
rz(-1.8291399) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8823223) q[0];
sx q[0];
rz(-0.81040183) q[0];
sx q[0];
rz(2.209254) q[0];
rz(1.8238292) q[1];
sx q[1];
rz(-1.980314) q[1];
sx q[1];
rz(-0.12277776) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30194382) q[0];
sx q[0];
rz(-1.0331863) q[0];
sx q[0];
rz(-1.4728236) q[0];
rz(2.3244395) q[2];
sx q[2];
rz(-1.3434402) q[2];
sx q[2];
rz(-2.0968693) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9726561) q[1];
sx q[1];
rz(-0.74311174) q[1];
sx q[1];
rz(-1.572322) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5769031) q[3];
sx q[3];
rz(-1.8204692) q[3];
sx q[3];
rz(2.3752729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5968898) q[2];
sx q[2];
rz(-0.92146102) q[2];
sx q[2];
rz(-2.5370562) q[2];
rz(-0.70945159) q[3];
sx q[3];
rz(-0.76990288) q[3];
sx q[3];
rz(-0.82367045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9619047) q[0];
sx q[0];
rz(-3.0178495) q[0];
sx q[0];
rz(1.7748348) q[0];
rz(-2.1429515) q[1];
sx q[1];
rz(-2.3256681) q[1];
sx q[1];
rz(2.6008115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78695145) q[0];
sx q[0];
rz(-1.7102512) q[0];
sx q[0];
rz(-2.9324173) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2774463) q[2];
sx q[2];
rz(-2.1808778) q[2];
sx q[2];
rz(2.6928201) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.65366983) q[1];
sx q[1];
rz(-1.2752295) q[1];
sx q[1];
rz(1.771677) q[1];
rz(-0.56638797) q[3];
sx q[3];
rz(-1.4175804) q[3];
sx q[3];
rz(2.9569616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6775386) q[2];
sx q[2];
rz(-2.3276734) q[2];
sx q[2];
rz(2.3133254) q[2];
rz(-1.6895435) q[3];
sx q[3];
rz(-1.6477081) q[3];
sx q[3];
rz(0.90788466) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0238817) q[0];
sx q[0];
rz(-1.1957059) q[0];
sx q[0];
rz(-1.9697795) q[0];
rz(-0.68012971) q[1];
sx q[1];
rz(-1.2260194) q[1];
sx q[1];
rz(0.5562869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16613516) q[0];
sx q[0];
rz(-1.3167456) q[0];
sx q[0];
rz(-1.182876) q[0];
x q[1];
rz(1.5519737) q[2];
sx q[2];
rz(-0.14655098) q[2];
sx q[2];
rz(0.77136794) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7032675) q[1];
sx q[1];
rz(-2.5653953) q[1];
sx q[1];
rz(0.11771113) q[1];
x q[2];
rz(-2.049796) q[3];
sx q[3];
rz(-0.99078876) q[3];
sx q[3];
rz(-0.73778462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80643001) q[2];
sx q[2];
rz(-1.8589636) q[2];
sx q[2];
rz(0.20509091) q[2];
rz(2.3388376) q[3];
sx q[3];
rz(-2.0358678) q[3];
sx q[3];
rz(-2.5870489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270808) q[0];
sx q[0];
rz(-1.0733805) q[0];
sx q[0];
rz(-2.9141973) q[0];
rz(0.79044509) q[1];
sx q[1];
rz(-2.3643654) q[1];
sx q[1];
rz(-2.3048185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9783426) q[0];
sx q[0];
rz(-1.0549091) q[0];
sx q[0];
rz(1.8229026) q[0];
rz(-pi) q[1];
rz(2.6322962) q[2];
sx q[2];
rz(-2.4894538) q[2];
sx q[2];
rz(-0.76667128) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0322276) q[1];
sx q[1];
rz(-1.9739474) q[1];
sx q[1];
rz(-1.228986) q[1];
x q[2];
rz(1.1985336) q[3];
sx q[3];
rz(-2.296631) q[3];
sx q[3];
rz(-1.6282631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5922015) q[2];
sx q[2];
rz(-1.6624007) q[2];
sx q[2];
rz(-2.5679892) q[2];
rz(0.32522374) q[3];
sx q[3];
rz(-2.2461788) q[3];
sx q[3];
rz(-2.7371244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10930546) q[0];
sx q[0];
rz(-0.17806299) q[0];
sx q[0];
rz(2.4626515) q[0];
rz(-3.0067054) q[1];
sx q[1];
rz(-2.0811681) q[1];
sx q[1];
rz(1.4237684) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37932351) q[0];
sx q[0];
rz(-2.0998619) q[0];
sx q[0];
rz(-0.53405116) q[0];
rz(-pi) q[1];
rz(-2.277454) q[2];
sx q[2];
rz(-2.0108962) q[2];
sx q[2];
rz(2.3789483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0404741) q[1];
sx q[1];
rz(-2.3777739) q[1];
sx q[1];
rz(1.6141255) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1072949) q[3];
sx q[3];
rz(-2.6141659) q[3];
sx q[3];
rz(-0.81125439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5474825) q[2];
sx q[2];
rz(-2.2308733) q[2];
sx q[2];
rz(0.2317079) q[2];
rz(-1.0101275) q[3];
sx q[3];
rz(-1.2332656) q[3];
sx q[3];
rz(2.2505984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4657087) q[0];
sx q[0];
rz(-2.3864585) q[0];
sx q[0];
rz(2.6278507) q[0];
rz(-1.4328009) q[1];
sx q[1];
rz(-1.865973) q[1];
sx q[1];
rz(1.6339711) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3604853) q[0];
sx q[0];
rz(-1.6193266) q[0];
sx q[0];
rz(3.0074708) q[0];
rz(-pi) q[1];
rz(-0.56030031) q[2];
sx q[2];
rz(-0.29688421) q[2];
sx q[2];
rz(2.1765944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1257612) q[1];
sx q[1];
rz(-1.4471701) q[1];
sx q[1];
rz(-0.13036734) q[1];
rz(-pi) q[2];
x q[2];
rz(1.474103) q[3];
sx q[3];
rz(-0.25287273) q[3];
sx q[3];
rz(2.0677572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12019176) q[2];
sx q[2];
rz(-2.2970565) q[2];
sx q[2];
rz(0.25944844) q[2];
rz(1.1546968) q[3];
sx q[3];
rz(-0.76609937) q[3];
sx q[3];
rz(1.8416789) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4907289) q[0];
sx q[0];
rz(-1.6614953) q[0];
sx q[0];
rz(1.65253) q[0];
rz(-0.64000714) q[1];
sx q[1];
rz(-1.0970486) q[1];
sx q[1];
rz(-2.1515813) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8278767) q[0];
sx q[0];
rz(-1.4826688) q[0];
sx q[0];
rz(0.14493305) q[0];
rz(-0.48374071) q[2];
sx q[2];
rz(-0.82925284) q[2];
sx q[2];
rz(-2.8791548) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.23259904) q[1];
sx q[1];
rz(-0.4505583) q[1];
sx q[1];
rz(-1.1934936) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.02496533) q[3];
sx q[3];
rz(-1.348266) q[3];
sx q[3];
rz(-2.3937267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9330357) q[2];
sx q[2];
rz(-1.8427589) q[2];
sx q[2];
rz(2.9162858) q[2];
rz(-2.2202668) q[3];
sx q[3];
rz(-0.39042979) q[3];
sx q[3];
rz(-3.0909753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13023547) q[0];
sx q[0];
rz(-0.36778944) q[0];
sx q[0];
rz(0.854048) q[0];
rz(1.8232952) q[1];
sx q[1];
rz(-1.6164936) q[1];
sx q[1];
rz(2.2093538) q[1];
rz(1.83056) q[2];
sx q[2];
rz(-1.4767892) q[2];
sx q[2];
rz(-1.3108419) q[2];
rz(2.1193567) q[3];
sx q[3];
rz(-1.4132186) q[3];
sx q[3];
rz(1.3506387) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
