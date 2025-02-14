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
rz(0.8402549) q[0];
sx q[0];
rz(-2.1036966) q[0];
sx q[0];
rz(0.84501803) q[0];
rz(0.71212274) q[1];
sx q[1];
rz(4.1432015) q[1];
sx q[1];
rz(10.92034) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63461958) q[0];
sx q[0];
rz(-1.15043) q[0];
sx q[0];
rz(0.93376055) q[0];
x q[1];
rz(-2.7833013) q[2];
sx q[2];
rz(-2.7346277) q[2];
sx q[2];
rz(0.55487061) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.481166) q[1];
sx q[1];
rz(-1.1200953) q[1];
sx q[1];
rz(-0.53498241) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.846662) q[3];
sx q[3];
rz(-1.1893236) q[3];
sx q[3];
rz(-1.671553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40319765) q[2];
sx q[2];
rz(-1.1824111) q[2];
sx q[2];
rz(2.5851868) q[2];
rz(-0.49499908) q[3];
sx q[3];
rz(-0.35111108) q[3];
sx q[3];
rz(-2.2569412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.86968386) q[0];
sx q[0];
rz(-0.10232919) q[0];
sx q[0];
rz(-0.62865692) q[0];
rz(2.2643845) q[1];
sx q[1];
rz(-2.6429206) q[1];
sx q[1];
rz(2.8884851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.077674) q[0];
sx q[0];
rz(-1.5627994) q[0];
sx q[0];
rz(-1.5863938) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5649389) q[2];
sx q[2];
rz(-0.51311427) q[2];
sx q[2];
rz(2.6000044) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3971218) q[1];
sx q[1];
rz(-1.6872477) q[1];
sx q[1];
rz(-0.20040234) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.289648) q[3];
sx q[3];
rz(-2.3901403) q[3];
sx q[3];
rz(0.40159097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46193281) q[2];
sx q[2];
rz(-1.2075281) q[2];
sx q[2];
rz(-2.0750462) q[2];
rz(1.8343532) q[3];
sx q[3];
rz(-2.6756838) q[3];
sx q[3];
rz(2.2333142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2147373) q[0];
sx q[0];
rz(-0.94811386) q[0];
sx q[0];
rz(2.1562449) q[0];
rz(-2.7477879) q[1];
sx q[1];
rz(-0.99291283) q[1];
sx q[1];
rz(2.6148112) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.55888) q[0];
sx q[0];
rz(-2.0979019) q[0];
sx q[0];
rz(3.0452791) q[0];
x q[1];
rz(2.5622732) q[2];
sx q[2];
rz(-2.0181106) q[2];
sx q[2];
rz(1.8471225) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64766769) q[1];
sx q[1];
rz(-0.45803775) q[1];
sx q[1];
rz(-1.1801932) q[1];
rz(-pi) q[2];
rz(2.6963077) q[3];
sx q[3];
rz(-1.0202583) q[3];
sx q[3];
rz(-0.30552826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3452611) q[2];
sx q[2];
rz(-2.3806206) q[2];
sx q[2];
rz(2.995028) q[2];
rz(2.2740299) q[3];
sx q[3];
rz(-0.87351322) q[3];
sx q[3];
rz(0.7578907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4124311) q[0];
sx q[0];
rz(-2.6456092) q[0];
sx q[0];
rz(-2.5878986) q[0];
rz(1.6665392) q[1];
sx q[1];
rz(-0.29860425) q[1];
sx q[1];
rz(-0.24436229) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8575154) q[0];
sx q[0];
rz(-1.2173442) q[0];
sx q[0];
rz(-1.9259324) q[0];
rz(-pi) q[1];
rz(0.8342488) q[2];
sx q[2];
rz(-1.8015727) q[2];
sx q[2];
rz(-2.1342056) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1235413) q[1];
sx q[1];
rz(-1.9885716) q[1];
sx q[1];
rz(-1.17735) q[1];
rz(-3.1126106) q[3];
sx q[3];
rz(-2.4292415) q[3];
sx q[3];
rz(-2.1952352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.060140572) q[2];
sx q[2];
rz(-1.1271366) q[2];
sx q[2];
rz(0.78090182) q[2];
rz(0.39685708) q[3];
sx q[3];
rz(-3.1271827) q[3];
sx q[3];
rz(-2.0675596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2395372) q[0];
sx q[0];
rz(-0.1665512) q[0];
sx q[0];
rz(-2.8848414) q[0];
rz(0.17770879) q[1];
sx q[1];
rz(-2.5959028) q[1];
sx q[1];
rz(-2.1977052) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0832053) q[0];
sx q[0];
rz(-2.5542066) q[0];
sx q[0];
rz(0.87736954) q[0];
rz(-0.98942049) q[2];
sx q[2];
rz(-2.7983449) q[2];
sx q[2];
rz(-2.1892669) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5748716) q[1];
sx q[1];
rz(-1.5751667) q[1];
sx q[1];
rz(-0.53968756) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6750608) q[3];
sx q[3];
rz(-1.0634907) q[3];
sx q[3];
rz(3.0477085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2029767) q[2];
sx q[2];
rz(-2.1111033) q[2];
sx q[2];
rz(2.9317648) q[2];
rz(-3.0948011) q[3];
sx q[3];
rz(-1.6822466) q[3];
sx q[3];
rz(0.34745026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(3.0824579) q[0];
sx q[0];
rz(-2.3248398) q[0];
sx q[0];
rz(1.9385852) q[0];
rz(-1.5164392) q[1];
sx q[1];
rz(-2.7639183) q[1];
sx q[1];
rz(-0.032940544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0432949) q[0];
sx q[0];
rz(-1.8445208) q[0];
sx q[0];
rz(2.0541463) q[0];
rz(0.52741427) q[2];
sx q[2];
rz(-1.048812) q[2];
sx q[2];
rz(-2.5185713) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6194544) q[1];
sx q[1];
rz(-0.65548766) q[1];
sx q[1];
rz(-2.5966899) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3108211) q[3];
sx q[3];
rz(-0.83163315) q[3];
sx q[3];
rz(-2.3297854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.55543) q[2];
sx q[2];
rz(-1.0032434) q[2];
sx q[2];
rz(-0.45247751) q[2];
rz(0.14719506) q[3];
sx q[3];
rz(-0.23284027) q[3];
sx q[3];
rz(-1.2991306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.074987) q[0];
sx q[0];
rz(-2.7050278) q[0];
sx q[0];
rz(0.35705745) q[0];
rz(-2.1287411) q[1];
sx q[1];
rz(-1.9722936) q[1];
sx q[1];
rz(0.036651932) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9215611) q[0];
sx q[0];
rz(-0.71418327) q[0];
sx q[0];
rz(2.4490687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42958625) q[2];
sx q[2];
rz(-1.4768147) q[2];
sx q[2];
rz(0.89236255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0554296) q[1];
sx q[1];
rz(-2.2698672) q[1];
sx q[1];
rz(-0.12075874) q[1];
rz(2.996309) q[3];
sx q[3];
rz(-0.031153208) q[3];
sx q[3];
rz(0.6980074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4770294) q[2];
sx q[2];
rz(-0.9534812) q[2];
sx q[2];
rz(-0.089740962) q[2];
rz(1.2612032) q[3];
sx q[3];
rz(-1.829105) q[3];
sx q[3];
rz(1.4242273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.66202128) q[0];
sx q[0];
rz(-0.57858545) q[0];
sx q[0];
rz(-1.4805967) q[0];
rz(-1.4501694) q[1];
sx q[1];
rz(-2.4822576) q[1];
sx q[1];
rz(0.75993842) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93913883) q[0];
sx q[0];
rz(-1.6869776) q[0];
sx q[0];
rz(-0.010384801) q[0];
rz(-pi) q[1];
rz(-0.6282232) q[2];
sx q[2];
rz(-1.6464273) q[2];
sx q[2];
rz(0.13720195) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7260598) q[1];
sx q[1];
rz(-0.59846717) q[1];
sx q[1];
rz(0.80962555) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8283056) q[3];
sx q[3];
rz(-1.0777359) q[3];
sx q[3];
rz(-2.8185533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8370168) q[2];
sx q[2];
rz(-1.2695856) q[2];
sx q[2];
rz(-0.059583511) q[2];
rz(2.7130821) q[3];
sx q[3];
rz(-0.84463745) q[3];
sx q[3];
rz(-2.5394411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2483599) q[0];
sx q[0];
rz(-2.8587274) q[0];
sx q[0];
rz(0.34348139) q[0];
rz(0.19613014) q[1];
sx q[1];
rz(-2.2562512) q[1];
sx q[1];
rz(2.3416065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752026) q[0];
sx q[0];
rz(-2.3227442) q[0];
sx q[0];
rz(-2.1455392) q[0];
rz(0.56234151) q[2];
sx q[2];
rz(-1.3571171) q[2];
sx q[2];
rz(0.7299698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5243036) q[1];
sx q[1];
rz(-0.56131786) q[1];
sx q[1];
rz(-0.08331068) q[1];
rz(-pi) q[2];
rz(-1.6558582) q[3];
sx q[3];
rz(-0.68455212) q[3];
sx q[3];
rz(0.61699293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34667748) q[2];
sx q[2];
rz(-0.76378834) q[2];
sx q[2];
rz(2.9128892) q[2];
rz(-1.4165357) q[3];
sx q[3];
rz(-1.718947) q[3];
sx q[3];
rz(2.4712839) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5423841) q[0];
sx q[0];
rz(-0.51775652) q[0];
sx q[0];
rz(1.9589348) q[0];
rz(0.54284894) q[1];
sx q[1];
rz(-0.12195568) q[1];
sx q[1];
rz(-2.6861232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96129721) q[0];
sx q[0];
rz(-1.1396798) q[0];
sx q[0];
rz(2.6102553) q[0];
rz(-pi) q[1];
rz(1.1588016) q[2];
sx q[2];
rz(-0.48589009) q[2];
sx q[2];
rz(-1.1166355) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7389981) q[1];
sx q[1];
rz(-1.7070424) q[1];
sx q[1];
rz(-0.46936492) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9712168) q[3];
sx q[3];
rz(-1.1569303) q[3];
sx q[3];
rz(-3.1088157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9246284) q[2];
sx q[2];
rz(-0.89322007) q[2];
sx q[2];
rz(-0.42579892) q[2];
rz(2.9567772) q[3];
sx q[3];
rz(-0.63822377) q[3];
sx q[3];
rz(-3.1137915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4375147) q[0];
sx q[0];
rz(-1.5298433) q[0];
sx q[0];
rz(-0.035506305) q[0];
rz(-0.47954814) q[1];
sx q[1];
rz(-1.5621114) q[1];
sx q[1];
rz(1.5893804) q[1];
rz(-1.3479955) q[2];
sx q[2];
rz(-0.98735129) q[2];
sx q[2];
rz(-2.9410887) q[2];
rz(-1.3666338) q[3];
sx q[3];
rz(-2.5208409) q[3];
sx q[3];
rz(2.3552166) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
