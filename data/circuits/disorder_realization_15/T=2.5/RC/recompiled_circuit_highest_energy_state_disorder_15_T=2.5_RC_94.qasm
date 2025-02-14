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
rz(-2.5315142) q[0];
sx q[0];
rz(-0.22992034) q[0];
sx q[0];
rz(-2.4218986) q[0];
rz(-0.62148681) q[1];
sx q[1];
rz(-1.4014333) q[1];
sx q[1];
rz(3.016234) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89594498) q[0];
sx q[0];
rz(-0.60310793) q[0];
sx q[0];
rz(-0.75624098) q[0];
rz(-pi) q[1];
rz(-2.8320931) q[2];
sx q[2];
rz(-0.30213854) q[2];
sx q[2];
rz(3.1378656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6290063) q[1];
sx q[1];
rz(-3.1401445) q[1];
sx q[1];
rz(-1.8667029) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1716648) q[3];
sx q[3];
rz(-2.3775775) q[3];
sx q[3];
rz(1.7872052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3434992) q[2];
sx q[2];
rz(-0.40830475) q[2];
sx q[2];
rz(0.84585345) q[2];
rz(0.79126233) q[3];
sx q[3];
rz(-0.013412272) q[3];
sx q[3];
rz(0.066789269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4132408) q[0];
sx q[0];
rz(-2.6744196) q[0];
sx q[0];
rz(0.043664232) q[0];
rz(1.5751669) q[1];
sx q[1];
rz(-1.3730201) q[1];
sx q[1];
rz(-1.6427737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3960172) q[0];
sx q[0];
rz(-1.5748576) q[0];
sx q[0];
rz(2.1745357) q[0];
rz(-pi) q[1];
rz(2.1421932) q[2];
sx q[2];
rz(-1.5599112) q[2];
sx q[2];
rz(1.5601666) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.0017569204) q[1];
sx q[1];
rz(-0.082232177) q[1];
sx q[1];
rz(0.38767345) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1083634) q[3];
sx q[3];
rz(-1.9150315) q[3];
sx q[3];
rz(2.0118464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.046039) q[2];
sx q[2];
rz(-2.9914896) q[2];
sx q[2];
rz(-0.52774876) q[2];
rz(0.85585099) q[3];
sx q[3];
rz(-3.1400883) q[3];
sx q[3];
rz(1.2214448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-0.39831487) q[0];
sx q[0];
rz(-2.1687431) q[0];
sx q[0];
rz(-1.995218) q[0];
rz(-1.758681) q[1];
sx q[1];
rz(-0.29255602) q[1];
sx q[1];
rz(3.0376099) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87097528) q[0];
sx q[0];
rz(-1.7954602) q[0];
sx q[0];
rz(-1.490834) q[0];
x q[1];
rz(3.1244794) q[2];
sx q[2];
rz(-1.2984973) q[2];
sx q[2];
rz(2.2382617) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5755363) q[1];
sx q[1];
rz(-0.78278274) q[1];
sx q[1];
rz(-2.8498883) q[1];
rz(-pi) q[2];
rz(2.7894734) q[3];
sx q[3];
rz(-2.8330292) q[3];
sx q[3];
rz(1.7087913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18342239) q[2];
sx q[2];
rz(-3.1347745) q[2];
sx q[2];
rz(2.5799694) q[2];
rz(-3.0567452) q[3];
sx q[3];
rz(-0.0054797879) q[3];
sx q[3];
rz(3.1408299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.4168004) q[0];
sx q[0];
rz(-2.8775207) q[0];
sx q[0];
rz(2.6208139) q[0];
rz(0.15781038) q[1];
sx q[1];
rz(-0.66693711) q[1];
sx q[1];
rz(-3.0657943) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9276088) q[0];
sx q[0];
rz(-0.99262041) q[0];
sx q[0];
rz(-1.4347721) q[0];
x q[1];
rz(-0.00064223991) q[2];
sx q[2];
rz(-1.5721605) q[2];
sx q[2];
rz(1.7096007) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0204649) q[1];
sx q[1];
rz(-1.1942099) q[1];
sx q[1];
rz(-0.15213206) q[1];
x q[2];
rz(0.10785477) q[3];
sx q[3];
rz(-2.0216935) q[3];
sx q[3];
rz(-1.0552561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6303404) q[2];
sx q[2];
rz(-3.1255836) q[2];
sx q[2];
rz(1.4207669) q[2];
rz(3.1347647) q[3];
sx q[3];
rz(-0.029191645) q[3];
sx q[3];
rz(-1.487287) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0066978) q[0];
sx q[0];
rz(-2.9338574) q[0];
sx q[0];
rz(-2.8790706) q[0];
rz(-0.94995704) q[1];
sx q[1];
rz(-3.0634395) q[1];
sx q[1];
rz(-0.21771678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9757468) q[0];
sx q[0];
rz(-1.4474156) q[0];
sx q[0];
rz(1.9666155) q[0];
rz(-1.4870013) q[2];
sx q[2];
rz(-1.6577273) q[2];
sx q[2];
rz(-3.0695335) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3847743) q[1];
sx q[1];
rz(-2.9692602) q[1];
sx q[1];
rz(1.5583355) q[1];
rz(-pi) q[2];
rz(1.4520549) q[3];
sx q[3];
rz(-0.53114) q[3];
sx q[3];
rz(-2.0599928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.21506423) q[2];
sx q[2];
rz(-3.1319517) q[2];
sx q[2];
rz(0.24791524) q[2];
rz(-1.7667814) q[3];
sx q[3];
rz(-3.091076) q[3];
sx q[3];
rz(-1.7063399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1011825) q[0];
sx q[0];
rz(-1.8721767) q[0];
sx q[0];
rz(-2.5293479) q[0];
rz(-2.9712037) q[1];
sx q[1];
rz(-3.0590765) q[1];
sx q[1];
rz(-1.4917397) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0675394) q[0];
sx q[0];
rz(-1.051386) q[0];
sx q[0];
rz(-0.75068604) q[0];
x q[1];
rz(1.5749919) q[2];
sx q[2];
rz(-1.5562662) q[2];
sx q[2];
rz(-2.5670403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41497624) q[1];
sx q[1];
rz(-1.6373937) q[1];
sx q[1];
rz(0.29331751) q[1];
rz(1.2475523) q[3];
sx q[3];
rz(-1.7446784) q[3];
sx q[3];
rz(-1.4617048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.49543574) q[2];
sx q[2];
rz(-0.010224552) q[2];
sx q[2];
rz(-0.23908991) q[2];
rz(0.80860364) q[3];
sx q[3];
rz(-3.1294398) q[3];
sx q[3];
rz(-1.3330207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4343774) q[0];
sx q[0];
rz(-0.093952976) q[0];
sx q[0];
rz(2.3648426) q[0];
rz(0.010919318) q[1];
sx q[1];
rz(-0.24591406) q[1];
sx q[1];
rz(-1.4728665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0162139) q[0];
sx q[0];
rz(-2.3998866) q[0];
sx q[0];
rz(1.8893695) q[0];
rz(0.02350925) q[2];
sx q[2];
rz(-1.5644367) q[2];
sx q[2];
rz(1.3964749) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90172988) q[1];
sx q[1];
rz(-2.3119676) q[1];
sx q[1];
rz(1.4338486) q[1];
rz(-pi) q[2];
rz(0.075447791) q[3];
sx q[3];
rz(-1.8561072) q[3];
sx q[3];
rz(2.4782654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9933068) q[2];
sx q[2];
rz(-3.1243656) q[2];
sx q[2];
rz(-0.41603184) q[2];
rz(2.3038583) q[3];
sx q[3];
rz(-3.0031524) q[3];
sx q[3];
rz(1.6637038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96529043) q[0];
sx q[0];
rz(-0.039881341) q[0];
sx q[0];
rz(0.95529977) q[0];
rz(1.924986) q[1];
sx q[1];
rz(-0.33733264) q[1];
sx q[1];
rz(-2.7359656) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55014706) q[0];
sx q[0];
rz(-2.2713775) q[0];
sx q[0];
rz(-1.9931384) q[0];
x q[1];
rz(-1.4816059) q[2];
sx q[2];
rz(-1.4124267) q[2];
sx q[2];
rz(-1.0430481) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.97598398) q[1];
sx q[1];
rz(-2.9635495) q[1];
sx q[1];
rz(-2.6958532) q[1];
rz(-pi) q[2];
rz(1.4838951) q[3];
sx q[3];
rz(-1.358506) q[3];
sx q[3];
rz(-3.1270554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8568628) q[2];
sx q[2];
rz(-0.036981985) q[2];
sx q[2];
rz(-2.3178318) q[2];
rz(3.01037) q[3];
sx q[3];
rz(-0.28990144) q[3];
sx q[3];
rz(0.32740617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4284978) q[0];
sx q[0];
rz(-2.9976124) q[0];
sx q[0];
rz(2.4289828) q[0];
rz(-2.5932942) q[1];
sx q[1];
rz(-0.24990853) q[1];
sx q[1];
rz(2.9329494) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8241103) q[0];
sx q[0];
rz(-1.9696931) q[0];
sx q[0];
rz(2.0092416) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0717355) q[2];
sx q[2];
rz(-3.102157) q[2];
sx q[2];
rz(0.60434228) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5968768) q[1];
sx q[1];
rz(-1.5660389) q[1];
sx q[1];
rz(-2.1265043) q[1];
rz(-2.194368) q[3];
sx q[3];
rz(-1.0290909) q[3];
sx q[3];
rz(0.052963363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.448552) q[2];
sx q[2];
rz(-0.081144944) q[2];
sx q[2];
rz(2.9244259) q[2];
rz(-0.36755696) q[3];
sx q[3];
rz(-0.03438545) q[3];
sx q[3];
rz(2.1448081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9656068) q[0];
sx q[0];
rz(-0.088986926) q[0];
sx q[0];
rz(-2.9678645) q[0];
rz(1.732775) q[1];
sx q[1];
rz(-1.4585739) q[1];
sx q[1];
rz(1.6857612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3575027) q[0];
sx q[0];
rz(-1.0561868) q[0];
sx q[0];
rz(1.7592488) q[0];
x q[1];
rz(1.3633492) q[2];
sx q[2];
rz(-1.4342299) q[2];
sx q[2];
rz(-0.88688021) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5599225) q[1];
sx q[1];
rz(-1.4538723) q[1];
sx q[1];
rz(2.440084) q[1];
rz(-1.3221856) q[3];
sx q[3];
rz(-1.4668307) q[3];
sx q[3];
rz(1.6373529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9659861) q[2];
sx q[2];
rz(-2.6484428) q[2];
sx q[2];
rz(-1.3803587) q[2];
rz(2.4557451) q[3];
sx q[3];
rz(-3.139747) q[3];
sx q[3];
rz(0.67870158) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3995517) q[0];
sx q[0];
rz(-0.69235943) q[0];
sx q[0];
rz(1.6211744) q[0];
rz(-1.5581268) q[1];
sx q[1];
rz(-1.635066) q[1];
sx q[1];
rz(0.20546694) q[1];
rz(1.6995399) q[2];
sx q[2];
rz(-0.12643554) q[2];
sx q[2];
rz(0.34496107) q[2];
rz(-1.2711502) q[3];
sx q[3];
rz(-1.4574403) q[3];
sx q[3];
rz(1.311655) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
