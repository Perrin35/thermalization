OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73206168) q[0];
sx q[0];
rz(-1.7763897) q[0];
sx q[0];
rz(2.1172297) q[0];
rz(-2.536474) q[1];
sx q[1];
rz(-2.6095698) q[1];
sx q[1];
rz(-1.1693118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071872358) q[0];
sx q[0];
rz(-0.9099996) q[0];
sx q[0];
rz(-1.1138492) q[0];
x q[1];
rz(-2.3697882) q[2];
sx q[2];
rz(-2.3183841) q[2];
sx q[2];
rz(-0.31847218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4989657) q[1];
sx q[1];
rz(-1.1945915) q[1];
sx q[1];
rz(2.7514003) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8291679) q[3];
sx q[3];
rz(-1.4989304) q[3];
sx q[3];
rz(-1.3952599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26596507) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(1.3226091) q[2];
rz(2.8406075) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(0.47505501) q[0];
rz(-1.3985727) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(1.0377201) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1514725) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(-2.5734076) q[0];
rz(-pi) q[1];
rz(0.77387626) q[2];
sx q[2];
rz(-2.4485588) q[2];
sx q[2];
rz(-2.3351923) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2270826) q[1];
sx q[1];
rz(-2.7343379) q[1];
sx q[1];
rz(1.6362908) q[1];
rz(0.59229895) q[3];
sx q[3];
rz(-2.2364738) q[3];
sx q[3];
rz(2.388282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(0.084687106) q[2];
rz(2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5304853) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(-2.5684165) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5983551) q[0];
sx q[0];
rz(-2.7063473) q[0];
sx q[0];
rz(1.7217365) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8163221) q[2];
sx q[2];
rz(-2.3420482) q[2];
sx q[2];
rz(2.6550967) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0108311) q[1];
sx q[1];
rz(-1.2847932) q[1];
sx q[1];
rz(-0.39949135) q[1];
rz(-pi) q[2];
rz(-1.1099986) q[3];
sx q[3];
rz(-1.2013544) q[3];
sx q[3];
rz(0.24584578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(2.9476681) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8594584) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(-0.55066806) q[0];
rz(-2.0129054) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(-0.36270025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4578611) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(1.9283717) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.018410725) q[2];
sx q[2];
rz(-1.5393179) q[2];
sx q[2];
rz(1.6056431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1945222) q[1];
sx q[1];
rz(-2.2010989) q[1];
sx q[1];
rz(1.6987726) q[1];
x q[2];
rz(1.6710715) q[3];
sx q[3];
rz(-0.62697151) q[3];
sx q[3];
rz(-0.37563045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(-2.4827042) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(3.127393) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(1.4594706) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0844903) q[0];
sx q[0];
rz(-1.4243037) q[0];
sx q[0];
rz(-1.3296207) q[0];
rz(-pi) q[1];
x q[1];
rz(2.622501) q[2];
sx q[2];
rz(-2.1862098) q[2];
sx q[2];
rz(1.5965243) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7534605) q[1];
sx q[1];
rz(-1.3452483) q[1];
sx q[1];
rz(-1.3688341) q[1];
x q[2];
rz(-1.0829809) q[3];
sx q[3];
rz(-1.8564463) q[3];
sx q[3];
rz(-2.3398427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(2.1250336) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(3.1275948) q[0];
sx q[0];
rz(-0.70677775) q[0];
sx q[0];
rz(2.5573964) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(-3.0029283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4045227) q[0];
sx q[0];
rz(-1.5911907) q[0];
sx q[0];
rz(1.5741332) q[0];
x q[1];
rz(0.3615985) q[2];
sx q[2];
rz(-0.21441678) q[2];
sx q[2];
rz(2.8117361) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0181959) q[1];
sx q[1];
rz(-1.0180078) q[1];
sx q[1];
rz(-2.4559896) q[1];
rz(-pi) q[2];
rz(0.74216446) q[3];
sx q[3];
rz(-1.8932749) q[3];
sx q[3];
rz(3.1165251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77506322) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(3.0464723) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(2.6050674) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(-2.5172863) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37492232) q[0];
sx q[0];
rz(-2.2793988) q[0];
sx q[0];
rz(1.2929582) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1820071) q[2];
sx q[2];
rz(-2.4215536) q[2];
sx q[2];
rz(0.50146539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99164167) q[1];
sx q[1];
rz(-2.2551564) q[1];
sx q[1];
rz(2.409694) q[1];
rz(0.47847139) q[3];
sx q[3];
rz(-1.0020743) q[3];
sx q[3];
rz(-2.3346321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6162993) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-2.7590511) q[2];
rz(-0.031490695) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.265825) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-3.0016622) q[0];
rz(-1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(0.12891842) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873868) q[0];
sx q[0];
rz(-1.6763858) q[0];
sx q[0];
rz(1.9681853) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56477408) q[2];
sx q[2];
rz(-2.161705) q[2];
sx q[2];
rz(1.5608982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1210632) q[1];
sx q[1];
rz(-1.9951092) q[1];
sx q[1];
rz(1.8706277) q[1];
rz(2.4017879) q[3];
sx q[3];
rz(-1.5590258) q[3];
sx q[3];
rz(1.7108325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-0.62409419) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36214608) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.8918442) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(-1.3508266) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6791145) q[0];
sx q[0];
rz(-0.20972855) q[0];
sx q[0];
rz(-1.0366584) q[0];
rz(-1.9753014) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(1.5664958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2495888) q[1];
sx q[1];
rz(-2.2715886) q[1];
sx q[1];
rz(2.4904817) q[1];
rz(0.25165598) q[3];
sx q[3];
rz(-0.76905426) q[3];
sx q[3];
rz(0.8151527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.089036971) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(-3.0272711) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(-0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(2.5949123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28323805) q[0];
sx q[0];
rz(-2.3843117) q[0];
sx q[0];
rz(-1.2454633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7231862) q[2];
sx q[2];
rz(-2.0414464) q[2];
sx q[2];
rz(-1.6413123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9676799) q[1];
sx q[1];
rz(-0.2914857) q[1];
sx q[1];
rz(3.013054) q[1];
rz(-pi) q[2];
rz(1.7095079) q[3];
sx q[3];
rz(-1.7943534) q[3];
sx q[3];
rz(0.15299882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13835779) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(-2.1440078) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(0.51013851) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99988408) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(-1.8757204) q[2];
sx q[2];
rz(-0.049449895) q[2];
sx q[2];
rz(0.54686875) q[2];
rz(-2.0426345) q[3];
sx q[3];
rz(-1.2755339) q[3];
sx q[3];
rz(-0.76225029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];