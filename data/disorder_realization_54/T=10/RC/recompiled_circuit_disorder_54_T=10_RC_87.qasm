OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78046507) q[0];
sx q[0];
rz(-2.4523003) q[0];
sx q[0];
rz(0.33049345) q[0];
rz(-2.8117872) q[1];
sx q[1];
rz(-2.2916315) q[1];
sx q[1];
rz(-0.70911521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1256975) q[0];
sx q[0];
rz(-0.96244922) q[0];
sx q[0];
rz(-0.35981052) q[0];
rz(-pi) q[1];
rz(1.4719226) q[2];
sx q[2];
rz(-0.31497248) q[2];
sx q[2];
rz(-2.0801983) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4830556) q[1];
sx q[1];
rz(-0.72209789) q[1];
sx q[1];
rz(2.2621821) q[1];
rz(-pi) q[2];
rz(0.96592824) q[3];
sx q[3];
rz(-2.5336207) q[3];
sx q[3];
rz(-2.0917497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4101397) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(-1.5343792) q[2];
rz(-0.93506995) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(-2.3527761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7222897) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(2.9653446) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(2.2252749) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57740649) q[0];
sx q[0];
rz(-2.2803377) q[0];
sx q[0];
rz(-2.2400411) q[0];
rz(-pi) q[1];
rz(1.4257405) q[2];
sx q[2];
rz(-2.3999891) q[2];
sx q[2];
rz(1.8650101) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20827046) q[1];
sx q[1];
rz(-1.9083438) q[1];
sx q[1];
rz(-0.34668215) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.022577062) q[3];
sx q[3];
rz(-1.9687679) q[3];
sx q[3];
rz(0.28385362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.24094412) q[2];
sx q[2];
rz(-2.1449461) q[2];
sx q[2];
rz(-2.3201578) q[2];
rz(-0.017283043) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(-0.78330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9538552) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(1.0082555) q[0];
rz(-0.035765212) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(-0.52454138) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7099972) q[0];
sx q[0];
rz(-1.4072197) q[0];
sx q[0];
rz(3.1113935) q[0];
rz(-pi) q[1];
rz(-2.8183297) q[2];
sx q[2];
rz(-2.4369536) q[2];
sx q[2];
rz(-1.7636253) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9279889) q[1];
sx q[1];
rz(-0.065918006) q[1];
sx q[1];
rz(-2.7830178) q[1];
rz(1.2479765) q[3];
sx q[3];
rz(-2.7016692) q[3];
sx q[3];
rz(2.3390714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3699469) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.6463722) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(-0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111073) q[0];
sx q[0];
rz(-0.91902584) q[0];
sx q[0];
rz(-0.088949732) q[0];
rz(-2.6308909) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(-0.68960062) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8423858) q[0];
sx q[0];
rz(-0.33873522) q[0];
sx q[0];
rz(-2.6329106) q[0];
rz(-pi) q[1];
rz(1.7961411) q[2];
sx q[2];
rz(-2.1225404) q[2];
sx q[2];
rz(-1.4622886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5213485) q[1];
sx q[1];
rz(-0.84940956) q[1];
sx q[1];
rz(1.1299302) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8702909) q[3];
sx q[3];
rz(-0.50008431) q[3];
sx q[3];
rz(-3.1262731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(0.62292567) q[2];
rz(2.0056491) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(0.46666551) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85161197) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(2.5581397) q[0];
rz(1.1460229) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(1.6437644) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4125815) q[0];
sx q[0];
rz(-2.3238365) q[0];
sx q[0];
rz(-1.2110932) q[0];
rz(2.3439212) q[2];
sx q[2];
rz(-1.5016218) q[2];
sx q[2];
rz(2.5728512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6795265) q[1];
sx q[1];
rz(-0.70409173) q[1];
sx q[1];
rz(0.042298869) q[1];
rz(-pi) q[2];
rz(-1.7329526) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(-1.4354524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(1.5779457) q[2];
rz(-2.2359713) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-2.5031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6435796) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(3.0946099) q[0];
rz(-2.9934096) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(1.4354338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5934138) q[0];
sx q[0];
rz(-2.3313064) q[0];
sx q[0];
rz(-0.20472783) q[0];
x q[1];
rz(1.6666744) q[2];
sx q[2];
rz(-0.55170689) q[2];
sx q[2];
rz(0.72223896) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.012949) q[1];
sx q[1];
rz(-1.1638906) q[1];
sx q[1];
rz(0.13601555) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8085603) q[3];
sx q[3];
rz(-1.5441582) q[3];
sx q[3];
rz(0.8386855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0753714) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(-3.0701239) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(2.215109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28618318) q[0];
sx q[0];
rz(-1.3278642) q[0];
sx q[0];
rz(2.563971) q[0];
rz(1.8619934) q[1];
sx q[1];
rz(-1.3459233) q[1];
sx q[1];
rz(-2.1320027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9855232) q[0];
sx q[0];
rz(-2.0445604) q[0];
sx q[0];
rz(-1.9457293) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6206714) q[2];
sx q[2];
rz(-2.705057) q[2];
sx q[2];
rz(1.5182564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13751444) q[1];
sx q[1];
rz(-0.61811781) q[1];
sx q[1];
rz(-2.532258) q[1];
rz(-pi) q[2];
rz(2.3791802) q[3];
sx q[3];
rz(-2.0250642) q[3];
sx q[3];
rz(1.7319958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0499095) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(1.1996777) q[2];
rz(2.6692634) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6427479) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(-2.902466) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(-2.696864) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037576588) q[0];
sx q[0];
rz(-2.4009279) q[0];
sx q[0];
rz(0.91233493) q[0];
rz(-0.66419454) q[2];
sx q[2];
rz(-2.1014629) q[2];
sx q[2];
rz(1.2150089) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7263033) q[1];
sx q[1];
rz(-0.88978926) q[1];
sx q[1];
rz(-2.9773832) q[1];
rz(-2.8430311) q[3];
sx q[3];
rz(-0.42793722) q[3];
sx q[3];
rz(0.096979389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42177054) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-0.81531173) q[2];
rz(2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(1.3172654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8274882) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(-2.8544193) q[0];
rz(2.9526967) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(-2.8093991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6824326) q[0];
sx q[0];
rz(-0.81572616) q[0];
sx q[0];
rz(-2.0072323) q[0];
rz(-pi) q[1];
rz(3.095093) q[2];
sx q[2];
rz(-0.87247712) q[2];
sx q[2];
rz(2.7547714) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13739535) q[1];
sx q[1];
rz(-1.157837) q[1];
sx q[1];
rz(-1.2910299) q[1];
rz(-1.166415) q[3];
sx q[3];
rz(-1.0746733) q[3];
sx q[3];
rz(-2.8465084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2087848) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(1.8509289) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(-0.65657842) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0861417) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(-1.0593876) q[0];
rz(-2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(2.8870781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33728889) q[0];
sx q[0];
rz(-0.24233195) q[0];
sx q[0];
rz(-1.2605577) q[0];
rz(-pi) q[1];
rz(-0.43912402) q[2];
sx q[2];
rz(-0.52999485) q[2];
sx q[2];
rz(-2.8119171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3855615) q[1];
sx q[1];
rz(-2.8061211) q[1];
sx q[1];
rz(2.8940593) q[1];
rz(-1.6700527) q[3];
sx q[3];
rz(-2.0072862) q[3];
sx q[3];
rz(1.0222767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.09482) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(-2.0937031) q[2];
rz(-1.3607599) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(-0.60539436) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142674) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(-0.36021532) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(-2.0965626) q[2];
sx q[2];
rz(-1.7797995) q[2];
sx q[2];
rz(-1.9893653) q[2];
rz(1.8923106) q[3];
sx q[3];
rz(-2.5026863) q[3];
sx q[3];
rz(1.7659059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
