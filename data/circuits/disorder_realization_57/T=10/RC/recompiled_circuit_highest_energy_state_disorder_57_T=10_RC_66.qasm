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
rz(1.1433831) q[0];
sx q[0];
rz(1.5513865) q[0];
sx q[0];
rz(9.0976465) q[0];
rz(1.7786572) q[1];
sx q[1];
rz(-2.5442446) q[1];
sx q[1];
rz(2.8391431) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0344447) q[0];
sx q[0];
rz(-1.4054789) q[0];
sx q[0];
rz(1.3470696) q[0];
rz(-2.2313868) q[2];
sx q[2];
rz(-1.1288092) q[2];
sx q[2];
rz(3.0750753) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0634583) q[1];
sx q[1];
rz(-1.7094104) q[1];
sx q[1];
rz(1.8891508) q[1];
x q[2];
rz(2.6685295) q[3];
sx q[3];
rz(-1.6636208) q[3];
sx q[3];
rz(-1.0582093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3737619) q[2];
sx q[2];
rz(-0.094002873) q[2];
sx q[2];
rz(2.8685699) q[2];
rz(-2.7772969) q[3];
sx q[3];
rz(-1.7184869) q[3];
sx q[3];
rz(3.1209768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.453422) q[0];
sx q[0];
rz(-2.0437129) q[0];
sx q[0];
rz(1.706644) q[0];
rz(-0.20225987) q[1];
sx q[1];
rz(-0.63019284) q[1];
sx q[1];
rz(0.30466255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019783171) q[0];
sx q[0];
rz(-2.2643548) q[0];
sx q[0];
rz(0.32525639) q[0];
rz(2.8424047) q[2];
sx q[2];
rz(-1.169489) q[2];
sx q[2];
rz(-0.68986675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2167336) q[1];
sx q[1];
rz(-1.2302515) q[1];
sx q[1];
rz(1.3751283) q[1];
rz(-pi) q[2];
rz(-0.10778569) q[3];
sx q[3];
rz(-1.6906212) q[3];
sx q[3];
rz(0.15080632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.27663407) q[2];
sx q[2];
rz(-0.43143299) q[2];
sx q[2];
rz(1.5383833) q[2];
rz(2.3254584) q[3];
sx q[3];
rz(-2.3147801) q[3];
sx q[3];
rz(-0.89474595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.5012539) q[0];
sx q[0];
rz(-0.29359874) q[0];
sx q[0];
rz(-1.5119934) q[0];
rz(-1.3410131) q[1];
sx q[1];
rz(-2.2684542) q[1];
sx q[1];
rz(0.27892932) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68088433) q[0];
sx q[0];
rz(-1.6795923) q[0];
sx q[0];
rz(0.1291424) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2401593) q[2];
sx q[2];
rz(-2.3877904) q[2];
sx q[2];
rz(1.9447226) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8388357) q[1];
sx q[1];
rz(-2.347337) q[1];
sx q[1];
rz(-1.1357186) q[1];
rz(-3.1380021) q[3];
sx q[3];
rz(-0.85455214) q[3];
sx q[3];
rz(0.2936756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5035847) q[2];
sx q[2];
rz(-1.6562485) q[2];
sx q[2];
rz(-0.33835641) q[2];
rz(1.407297) q[3];
sx q[3];
rz(-1.1275848) q[3];
sx q[3];
rz(-0.97051632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9426743) q[0];
sx q[0];
rz(-0.022138683) q[0];
sx q[0];
rz(-0.66377798) q[0];
rz(-0.55555073) q[1];
sx q[1];
rz(-1.5063565) q[1];
sx q[1];
rz(-1.3551855) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91265011) q[0];
sx q[0];
rz(-1.602409) q[0];
sx q[0];
rz(-2.4949882) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6333601) q[2];
sx q[2];
rz(-2.782397) q[2];
sx q[2];
rz(-2.0421093) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0993886) q[1];
sx q[1];
rz(-2.4335833) q[1];
sx q[1];
rz(0.96570496) q[1];
rz(-pi) q[2];
rz(-1.7029477) q[3];
sx q[3];
rz(-2.6717917) q[3];
sx q[3];
rz(-2.2061004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8252455) q[2];
sx q[2];
rz(-2.1630042) q[2];
sx q[2];
rz(2.6157731) q[2];
rz(2.53287) q[3];
sx q[3];
rz(-2.5059097) q[3];
sx q[3];
rz(2.4354602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71259251) q[0];
sx q[0];
rz(-1.1845931) q[0];
sx q[0];
rz(-2.1767966) q[0];
rz(-2.7305799) q[1];
sx q[1];
rz(-0.95714584) q[1];
sx q[1];
rz(1.1119276) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2913366) q[0];
sx q[0];
rz(-1.5550858) q[0];
sx q[0];
rz(-1.5474623) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0093727) q[2];
sx q[2];
rz(-2.2242332) q[2];
sx q[2];
rz(-2.6613622) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7605684) q[1];
sx q[1];
rz(-2.3443017) q[1];
sx q[1];
rz(-2.3098195) q[1];
x q[2];
rz(-0.26740243) q[3];
sx q[3];
rz(-1.4286213) q[3];
sx q[3];
rz(-0.09718516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8370221) q[2];
sx q[2];
rz(-1.8111753) q[2];
sx q[2];
rz(2.0988317) q[2];
rz(1.9581095) q[3];
sx q[3];
rz(-2.2730998) q[3];
sx q[3];
rz(-2.1709501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(2.9633164) q[0];
sx q[0];
rz(-1.2238418) q[0];
sx q[0];
rz(0.29603145) q[0];
rz(0.0010679642) q[1];
sx q[1];
rz(-1.8031392) q[1];
sx q[1];
rz(1.0329049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6955687) q[0];
sx q[0];
rz(-3.1349413) q[0];
sx q[0];
rz(0.59799303) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0216713) q[2];
sx q[2];
rz(-1.9397221) q[2];
sx q[2];
rz(-3.1110087) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5653648) q[1];
sx q[1];
rz(-1.8180483) q[1];
sx q[1];
rz(1.5134211) q[1];
x q[2];
rz(2.37866) q[3];
sx q[3];
rz(-1.2733394) q[3];
sx q[3];
rz(0.92253387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7615776) q[2];
sx q[2];
rz(-1.6753316) q[2];
sx q[2];
rz(-1.638691) q[2];
rz(-2.5401529) q[3];
sx q[3];
rz(-2.1419339) q[3];
sx q[3];
rz(0.39826605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10694207) q[0];
sx q[0];
rz(-1.5926188) q[0];
sx q[0];
rz(-0.73367992) q[0];
rz(-0.94379464) q[1];
sx q[1];
rz(-1.5993886) q[1];
sx q[1];
rz(-0.61818799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9551463) q[0];
sx q[0];
rz(-1.6206564) q[0];
sx q[0];
rz(-1.7503034) q[0];
x q[1];
rz(2.3011977) q[2];
sx q[2];
rz(-1.5638509) q[2];
sx q[2];
rz(0.23887979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.88981121) q[1];
sx q[1];
rz(-1.5446583) q[1];
sx q[1];
rz(-2.2356116) q[1];
rz(-pi) q[2];
rz(2.5060081) q[3];
sx q[3];
rz(-1.8653324) q[3];
sx q[3];
rz(2.9377191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9435297) q[2];
sx q[2];
rz(-1.9518096) q[2];
sx q[2];
rz(2.5243536) q[2];
rz(-1.9700358) q[3];
sx q[3];
rz(-2.7660683) q[3];
sx q[3];
rz(2.530781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1307369) q[0];
sx q[0];
rz(-2.3541088) q[0];
sx q[0];
rz(-2.6686344) q[0];
rz(-1.4594151) q[1];
sx q[1];
rz(-1.727203) q[1];
sx q[1];
rz(0.032729538) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6115295) q[0];
sx q[0];
rz(-0.034688799) q[0];
sx q[0];
rz(-2.7682224) q[0];
rz(-0.76489438) q[2];
sx q[2];
rz(-1.9143606) q[2];
sx q[2];
rz(-2.9091121) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53712691) q[1];
sx q[1];
rz(-0.49931881) q[1];
sx q[1];
rz(-2.0313655) q[1];
rz(1.8440014) q[3];
sx q[3];
rz(-3.0298067) q[3];
sx q[3];
rz(-0.32958406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98256436) q[2];
sx q[2];
rz(-0.50243598) q[2];
sx q[2];
rz(-1.095613) q[2];
rz(-0.76121965) q[3];
sx q[3];
rz(-1.8571564) q[3];
sx q[3];
rz(0.42263862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6765167) q[0];
sx q[0];
rz(-0.07971555) q[0];
sx q[0];
rz(-2.4327143) q[0];
rz(-2.830016) q[1];
sx q[1];
rz(-1.0919002) q[1];
sx q[1];
rz(0.30837217) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4080183) q[0];
sx q[0];
rz(-2.7112506) q[0];
sx q[0];
rz(-0.4056613) q[0];
rz(-pi) q[1];
rz(-2.4089085) q[2];
sx q[2];
rz(-0.53712979) q[2];
sx q[2];
rz(2.2718475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7180135) q[1];
sx q[1];
rz(-2.1732959) q[1];
sx q[1];
rz(1.275255) q[1];
rz(-pi) q[2];
rz(-2.8838168) q[3];
sx q[3];
rz(-0.80398241) q[3];
sx q[3];
rz(2.3641158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0333905) q[2];
sx q[2];
rz(-1.4806925) q[2];
sx q[2];
rz(0.22076503) q[2];
rz(2.0681785) q[3];
sx q[3];
rz(-2.150841) q[3];
sx q[3];
rz(-2.3555135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6440426) q[0];
sx q[0];
rz(-2.2798517) q[0];
sx q[0];
rz(-1.0368689) q[0];
rz(1.4700302) q[1];
sx q[1];
rz(-2.122888) q[1];
sx q[1];
rz(1.9482313) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84813335) q[0];
sx q[0];
rz(-2.6984697) q[0];
sx q[0];
rz(-0.070290914) q[0];
rz(-0.5035053) q[2];
sx q[2];
rz(-0.76447884) q[2];
sx q[2];
rz(-1.6703005) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2630849) q[1];
sx q[1];
rz(-0.89205884) q[1];
sx q[1];
rz(-1.3287067) q[1];
rz(-1.8032298) q[3];
sx q[3];
rz(-1.6940191) q[3];
sx q[3];
rz(-0.49613813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4736453) q[2];
sx q[2];
rz(-1.2845984) q[2];
sx q[2];
rz(0.13772193) q[2];
rz(1.600945) q[3];
sx q[3];
rz(-0.23156229) q[3];
sx q[3];
rz(-0.93129492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4731428) q[0];
sx q[0];
rz(-1.4746329) q[0];
sx q[0];
rz(2.0869577) q[0];
rz(-2.6558381) q[1];
sx q[1];
rz(-1.9172485) q[1];
sx q[1];
rz(1.6419372) q[1];
rz(-0.68436868) q[2];
sx q[2];
rz(-1.3791313) q[2];
sx q[2];
rz(-2.287938) q[2];
rz(-1.0578591) q[3];
sx q[3];
rz(-2.2649962) q[3];
sx q[3];
rz(-1.0640127) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
