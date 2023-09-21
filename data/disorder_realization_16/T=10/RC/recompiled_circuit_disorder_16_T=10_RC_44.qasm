OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(-2.3119976) q[0];
sx q[0];
rz(-0.15396804) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(-2.1492465) q[1];
sx q[1];
rz(-0.33831236) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4727288) q[0];
sx q[0];
rz(-0.40574408) q[0];
sx q[0];
rz(2.2292482) q[0];
rz(0.89262427) q[2];
sx q[2];
rz(-1.2783588) q[2];
sx q[2];
rz(-0.48722789) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0677883) q[1];
sx q[1];
rz(-1.5177625) q[1];
sx q[1];
rz(-2.8552613) q[1];
rz(-pi) q[2];
rz(0.015720856) q[3];
sx q[3];
rz(-1.058488) q[3];
sx q[3];
rz(-2.6920464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14264318) q[2];
sx q[2];
rz(-0.34037408) q[2];
sx q[2];
rz(-1.9677229) q[2];
rz(-3.0657892) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(-3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3409815) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(-0.064963438) q[0];
rz(2.5669572) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(1.2423135) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500538) q[0];
sx q[0];
rz(-2.8613052) q[0];
sx q[0];
rz(-2.6602402) q[0];
rz(-pi) q[1];
x q[1];
rz(1.785948) q[2];
sx q[2];
rz(-2.4977376) q[2];
sx q[2];
rz(1.3607963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2989267) q[1];
sx q[1];
rz(-1.0416404) q[1];
sx q[1];
rz(1.7556721) q[1];
rz(-pi) q[2];
rz(2.4317125) q[3];
sx q[3];
rz(-1.9544365) q[3];
sx q[3];
rz(-2.3185454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3339281) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(0.58369613) q[2];
rz(2.5675473) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(0.42049256) q[0];
sx q[0];
rz(-0.89389602) q[0];
sx q[0];
rz(-2.4131391) q[0];
rz(1.4942253) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(-1.0167936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3166312) q[0];
sx q[0];
rz(-1.5350071) q[0];
sx q[0];
rz(0.081598452) q[0];
rz(-pi) q[1];
rz(2.617308) q[2];
sx q[2];
rz(-1.1384083) q[2];
sx q[2];
rz(-2.1081032) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2336725) q[1];
sx q[1];
rz(-0.83041149) q[1];
sx q[1];
rz(-2.9552712) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.73909) q[3];
sx q[3];
rz(-2.8561391) q[3];
sx q[3];
rz(1.1807549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3399405) q[2];
sx q[2];
rz(-1.5321956) q[2];
sx q[2];
rz(2.0920848) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(-1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291572) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(1.4720434) q[0];
rz(-2.4064348) q[1];
sx q[1];
rz(-2.3627294) q[1];
sx q[1];
rz(0.24681117) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26537672) q[0];
sx q[0];
rz(-2.2582158) q[0];
sx q[0];
rz(-0.19296293) q[0];
rz(0.39101379) q[2];
sx q[2];
rz(-0.51119971) q[2];
sx q[2];
rz(-0.15399394) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88627316) q[1];
sx q[1];
rz(-1.6233994) q[1];
sx q[1];
rz(-2.7707151) q[1];
x q[2];
rz(-1.9645343) q[3];
sx q[3];
rz(-2.4871832) q[3];
sx q[3];
rz(-0.27458336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4776769) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(1.5412615) q[2];
rz(0.70704308) q[3];
sx q[3];
rz(-2.0440846) q[3];
sx q[3];
rz(0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33070579) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(1.3274308) q[0];
rz(1.5785626) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(-0.24838233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7029593) q[0];
sx q[0];
rz(-0.54134936) q[0];
sx q[0];
rz(-0.066141733) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18647285) q[2];
sx q[2];
rz(-1.4154134) q[2];
sx q[2];
rz(1.965167) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0866962) q[1];
sx q[1];
rz(-0.61453648) q[1];
sx q[1];
rz(2.8413248) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31110839) q[3];
sx q[3];
rz(-2.4690383) q[3];
sx q[3];
rz(-2.6274519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7053232) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-0.79745897) q[2];
rz(-2.752839) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(-2.6388772) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0080863) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(1.3457993) q[0];
rz(1.0812409) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(-0.1246917) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7170982) q[0];
sx q[0];
rz(-1.4846804) q[0];
sx q[0];
rz(2.9647102) q[0];
x q[1];
rz(-2.5325534) q[2];
sx q[2];
rz(-1.8134724) q[2];
sx q[2];
rz(0.86953029) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9166959) q[1];
sx q[1];
rz(-1.5142913) q[1];
sx q[1];
rz(2.4042261) q[1];
x q[2];
rz(-0.74495875) q[3];
sx q[3];
rz(-0.31674851) q[3];
sx q[3];
rz(1.6572286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6283915) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(0.61895269) q[2];
rz(2.0882873) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(-0.9296023) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597647) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(1.0429617) q[0];
rz(-0.46328059) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(-1.0707062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8585513) q[0];
sx q[0];
rz(-1.0867449) q[0];
sx q[0];
rz(-1.5714684) q[0];
rz(-2.7878739) q[2];
sx q[2];
rz(-0.36834799) q[2];
sx q[2];
rz(-3.0596717) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0740944) q[1];
sx q[1];
rz(-2.3405582) q[1];
sx q[1];
rz(1.4317516) q[1];
rz(-pi) q[2];
rz(0.33220746) q[3];
sx q[3];
rz(-1.5515965) q[3];
sx q[3];
rz(1.7393877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7730007) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(0.32361844) q[2];
rz(2.1598024) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(-1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48802808) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(-1.9352242) q[0];
rz(1.9288829) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(-0.94747296) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3377209) q[0];
sx q[0];
rz(-2.5944355) q[0];
sx q[0];
rz(-1.7659811) q[0];
x q[1];
rz(-1.4119554) q[2];
sx q[2];
rz(-1.3961332) q[2];
sx q[2];
rz(2.5607002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.84711134) q[1];
sx q[1];
rz(-1.0080907) q[1];
sx q[1];
rz(0.28114762) q[1];
x q[2];
rz(-1.2097589) q[3];
sx q[3];
rz(-2.3449538) q[3];
sx q[3];
rz(2.3494997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56132135) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(-2.6718111) q[2];
rz(1.8404768) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(-2.8619213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.8895421) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(-1.3289733) q[0];
rz(2.3503616) q[1];
sx q[1];
rz(-2.8104517) q[1];
sx q[1];
rz(2.9387617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348939) q[0];
sx q[0];
rz(-1.5968423) q[0];
sx q[0];
rz(0.039020122) q[0];
x q[1];
rz(0.22325309) q[2];
sx q[2];
rz(-1.3666144) q[2];
sx q[2];
rz(2.5399361) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.054246) q[1];
sx q[1];
rz(-0.11212238) q[1];
sx q[1];
rz(-2.0263158) q[1];
rz(-2.8392302) q[3];
sx q[3];
rz(-0.48067579) q[3];
sx q[3];
rz(-2.9722948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.41708502) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(2.1450796) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0614232) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(-0.18173519) q[0];
rz(-0.043047992) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(-2.8607686) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94175324) q[0];
sx q[0];
rz(-1.117525) q[0];
sx q[0];
rz(1.0928632) q[0];
rz(1.3781204) q[2];
sx q[2];
rz(-1.780605) q[2];
sx q[2];
rz(-3.0378621) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.909543) q[1];
sx q[1];
rz(-0.068040158) q[1];
sx q[1];
rz(-2.1310991) q[1];
x q[2];
rz(1.0481846) q[3];
sx q[3];
rz(-2.7624353) q[3];
sx q[3];
rz(-0.89695938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8250371) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(-2.5184856) q[2];
rz(-2.1394219) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(-0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6476718) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(2.2676246) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(0.59861029) q[2];
sx q[2];
rz(-2.241588) q[2];
sx q[2];
rz(-0.66551756) q[2];
rz(2.279083) q[3];
sx q[3];
rz(-0.52048341) q[3];
sx q[3];
rz(-2.4181548) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
